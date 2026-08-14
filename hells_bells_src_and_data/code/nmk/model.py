import torch, torch.nn as nn
import torch.nn.functional as F
import math

DTYPE = torch.float64

def ternary_ready_init(w: nn.Parameter, p_non_zero: float):
    """
    Xavier-uniform core + sparse ternary impulse.
      • p_non_zero   - fraction of weights we pre-place at ±1
    """
    fan_in, _ = nn.init._calculate_fan_in_and_fan_out(w)
    bound = math.sqrt(6.0 / fan_in)          # Xavier uniform ±sqrt(6/fan_in)
    nn.init.uniform_(w, -bound, bound)       # keeps variance ~ 1/fan_in

    # Sprinkle a few exact ±1 values to speed up ternarisation - don't know if this does anything
    mask = torch.rand_like(w) < p_non_zero
    signs = torch.where(torch.rand_like(w) < 0.5, -1.0, 1.0)
    w.data = torch.where(mask, signs, w)     # replace subset with ±1


class BilinearMatMul(nn.Module):
    def __init__(self, n, m, k, r, dtype=torch.float32):
        super().__init__()
        self.n = n
        self.m = m
        self.k = k
        self.r = r
        self.dtype = dtype
        
        # Input: A (batch, n, m) -> (batch, n*m)
        #        B (batch, m, k) -> (batch, m*k)
        # Weights shape: (r, n*m) and (r, m*k)
        self.WA = nn.Linear(n*m, r, bias=False, dtype=dtype)
        self.WB = nn.Linear(m*k, r, bias=False, dtype=dtype)
        self.WC = nn.Linear(r, n*k, bias=False, dtype=dtype)
        
        # Initialize with small random values
        for w in (self.WA.weight, self.WB.weight, self.WC.weight):
            ternary_ready_init(w, p_non_zero=0.1)   #   <-- key line, imparts selection bias
        
    def forward(self, x):
        A, B = x
        bsz = A.shape[0]
        
        # Flatten matrices
        A_flat = A.reshape(bsz, -1)
        B_flat = B.reshape(bsz, -1)
        
        # Apply linear transformations
        a_features = self.WA(A_flat)
        b_features = self.WB(B_flat)
        
        # Element-wise multiplication: (batch, r)
        products = a_features * b_features
        
        # Final transformation: (batch, r) -> (batch, n*k)
        C_flat = self.WC(products)
        
        # Reshape back: (batch, n*k) -> (batch, n, k)
        C = C_flat.reshape(bsz, self.n, self.k)
        
        return C


def ternary_regularization(model, *, lambda_sparse: float, lambda_ternary: float, temperature: float):
    total_reg = 0.0
    
    for name, param in model.named_parameters():
        if 'weight' not in name:
            # Biases are left unconstrained (there are none in BilinearMatMul, but good practice)
            continue

        sparsity_term = lambda_sparse * param.abs().sum()

        # Euclidean distances to every ternary anchor
        # dist_to_zero = param.pow(2)
        # dist_to_one = (param - 1).pow(2)
        # dist_to_neg_one = (param + 1).pow(2)

        # can also do absolute:
        dist_to_zero = param.abs()
        dist_to_one = (param - 1).abs()
        dist_to_neg_one = (param + 1).abs()

        # Stack & compute a smooth minimum of the three distances
        distances = torch.stack([dist_to_zero, dist_to_one, dist_to_neg_one], dim=0)
        smooth_min = -temperature * torch.logsumexp(-distances / temperature, dim=0)
        ternary_term = lambda_ternary * smooth_min.sum()
        total_reg = total_reg + ternary_term + sparsity_term
    
    return total_reg


def eval_loss(net, A, B):  # evaluation only (no grads to ensure max efficiency)
    was_training = net.training
    net.eval()
    with torch.no_grad():
        logits = net((A, B))
        truth = torch.matmul(A, B)
        loss = torch.mean((logits - truth) ** 2)
    # Restore original training state
    net.train(was_training)
    return loss

def train_step(net, opt, A, B, *, lambda_sparse: float, lambda_ternary: float, temperature: float):
    opt.zero_grad()
    
    # Compute the forward pass
    C_pred = net((A, B))
    C_true = torch.matmul(A, B)
    
    # Mean squared error loss
    mse_loss = F.mse_loss(C_pred, C_true)
    # struct_loss = structured_mat_loss(C_pred, C_true)
    
    # Add (scheduled) ternary regularization
    reg_loss = ternary_regularization( # could we add the determinant in the loss? Or some other informant of matrix product?
        net,
        lambda_sparse=lambda_sparse,
        lambda_ternary=lambda_ternary,
        temperature=temperature,
    )
    
    # Total loss
    loss = mse_loss + reg_loss
    loss.backward()
    opt.step()
    return loss


