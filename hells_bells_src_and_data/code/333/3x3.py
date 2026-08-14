import torch
from torch import optim
from model import train_step, eval_loss, BilinearMatMul
from verification import snap, verify
import argparse
from hyperparams import TrainingArgs
from hyperparams import _print_hparams

import sys
import io
sys.stdin.reconfigure(encoding='utf-8')
sys.stdout.reconfigure(encoding='utf-8')

_parser = argparse.ArgumentParser(add_help=False)
_parser.add_argument("--seed", type=int)
_parser.add_argument("--lr", type=float, dest="lr_init")
_parser.add_argument("--lambda_ternary_end", type=float)
_parser.add_argument("--anneal_steps", type=int)
_parser.add_argument("--enable_ternary_at", type=float)
_parser.add_argument("--epochs", type=int)
_parser.add_argument("--lr_reduction_threshold", type=float)
_parser.add_argument("--converge_loss_threshold", type=float)
_parser.add_argument("--converge_ternary_percent", type=float)
_parser.add_argument("--lr_reduction_value", type=float)
_CLI_ARGS, _ = _parser.parse_known_args()

args = TrainingArgs(_CLI_ARGS)
_print_hparams(args)
torch.manual_seed(args.seed)
net = BilinearMatMul(args.n, args.r, dtype=args.dtype).to(args.device)

if args.optimizer == "Adam":
    opt = optim.Adam(net.parameters(), lr=args.lr_init, betas=(args.beta_1, args.beta_2), weight_decay=args.adam_weight_decay)
elif args.optimizer == "AdamW":
    opt = optim.AdamW(net.parameters(), lr=args.lr_init, betas=(args.beta_1, args.beta_2), weight_decay=args.adam_weight_decay)

ema = float("inf")
best = float("inf")
ternary_start_step = None

for step in range(1, args.epochs):
    # Generate a fresh random batch each step - idea is that this encourages generalization but no idea if this actually works
    A = torch.randn(args.batch, args.n, args.n, device=args.device, dtype=args.dtype)
    B = torch.randn_like(A)

    # Normalise to unit Frobenius norm – keeps targets O(1) and stabilises MSE scale
    A = A / torch.linalg.norm(A, ord='fro', dim=(1, 2), keepdim=True)
    B = B / torch.linalg.norm(B, ord='fro', dim=(1, 2), keepdim=True)

    # Decide if ternary should be enabled
    if ternary_start_step is not None:
        ternary_progress = step - ternary_start_step
        λ_T, temp = args.coeffs(step, ternary_progress)
    else:
       λ_T, temp = 0.0, args.temperature_start
    sparse_coeff = args.lambda_sparse if ternary_start_step is not None else 0.0

    loss = train_step(net, opt, A, B, lambda_sparse=sparse_coeff, lambda_ternary=λ_T, temperature=temp)

    if step % 100 == 0:
        with torch.no_grad():
            val_loss = eval_loss(net, A, B)

            if (val_loss < args.lr_reduction_threshold) and (ternary_start_step is None):
                for g in opt.param_groups:
                    g["lr"] = args.lr_reduction_value

            w = torch.cat([p.flatten() for p in net.parameters()])
            weight_norm = torch.linalg.norm(w, ord=2)
            ternary = 100 * (1 - torch.min(torch.stack([w.abs(), (w - 1).abs(), (w + 1).abs()]), 0)[0].mean())
            current_lr = opt.param_groups[0]["lr"]
            print(f"[Epoch {step:,}, Val Loss: {val_loss:.2e}, Weight norm: {weight_norm:.2e}, Ternary: {ternary:.2f}, LR: {current_lr:.2e}")

        # Switch on ternary regularisation once we cross the accuracy bar
        if ternary_start_step is None and val_loss < args.enable_ternary_at:
            ternary_start_step = step
            for g in opt.param_groups:
                g["lr"] = args.lr_init
            print(f"\n→ Enabled ternary regularisation at step {step:,} (val_loss={val_loss:.2e})")

        best = min(best, val_loss)

        # Additional simple convergence check (use validation loss)
        if (val_loss < args.converge_loss_threshold) and (ternary > args.converge_ternary_percent): # This can definitely be optimised
            break

literals = snap(net.WA.weight.detach().numpy(), 
                net.WB.weight.detach().numpy(), 
                net.WC.weight.detach().numpy())

try:
    verify(literals[0], literals[1], literals[2])
except AssertionError:
    raise