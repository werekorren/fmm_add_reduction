import torch
import os

# Helpers to normalise container override values
def _parse_dtype(name: str):
    s = (name or "").strip().lower()
    if s in {"float32", "fp32", "torch.float32", "f32"}:
        return torch.float32
    if s in {"bfloat16", "bf16", "torch.bfloat16"}:
        return torch.bfloat16
    if s in {"float16", "half", "fp16", "torch.float16", "f16"}:
        return torch.float16
    return torch.float32

def _normalize_optimizer(name: str) -> str:
    s = (name or "").strip().lower()
    if s in {"adamw", "adam_w", "adam-w"}:
        return "AdamW"
    if s == "adam":
        return "Adam"
    return name or "AdamW"

def _stringify_value(v):
    if isinstance(v, torch.dtype):
        return str(v).replace("torch.", "")
    return v

def _print_hparams(args_obj) -> None:
    print("=== Run configuration (hyperparameters) ===")
    for key in sorted(vars(args_obj).keys()):
        val = getattr(args_obj, key)
        print(f"- {key}: {_stringify_value(val)}")

# Adaptive schedule - bump up weight decay if weight norm is exloding, bump up lr if ternarity is plateuing etc.
# Is it possible to give Adam memory such that it knows where to return to?

class TrainingArgs:  
    def __init__(self, _CLI_ARGS):
        self.n, self.m, self.k, self.r = 2, 3, 3, 15
        self.epochs = 25_000
        self.batch = 256
        self.optimizer = "AdamW"
        # self.optimizer = "Muon"

        # beta_1 is momentum
        self.beta_1 = 0.9
        self.beta_2 = 0.97
        self.adam_weight_decay = 0.01

        self.anneal_steps = 150_000
        self.enable_ternary_at = 3e-6

        # Optimiser / LR schedule
        self.lr_init = 8e-3

        # Thresholds
        self.lr_reduction_threshold = 1e-5
        self.lr_reduction_value = 3e-3
        self.converge_loss_threshold = 3e-6
        self.converge_ternary_percent = 99.0

        # Sparse penalty value once ternary regularisation starts (0.0 before)
        self.lambda_sparse = 3e-6 
        self.lambda_ternary_end = 1e-4
        self.temperature_start = 1e-3
        self.temperature_end = 1e-3
        
        # Misc
        self.device = "cuda" if torch.cuda.is_available() else "cpu"
        self.dtype = torch.float32
        self.seed = 42

        # Container-level overrides via environment variables (AWS Batch overrides)
        # Accept uppercase names like EPOCHS, BATCH, LR_INIT, etc. CLI flags still take precedence later.
        env_map_int = {
            "N": "n",
            "M": "m",
            "K": "k",
            "R": "r",
            "EPOCHS": "epochs",
            "BATCH": "batch",
            "ANNEAL_STEPS": "anneal_steps",
            "SEED": "seed",
        }
        env_map_float = {
            "BETA_1": "beta_1",
            "BETA_2": "beta_2",
            "ADAM_WEIGHT_DECAY": "adam_weight_decay",
            "LR_REDUCTION_THRESHOLD": "lr_reduction_threshold",
            "LR_REDUCTION_VALUE": "lr_reduction_value",
            "CONVERGE_LOSS_THRESHOLD": "converge_loss_threshold",
            "CONVERGE_TERNARY_PERCENT": "converge_ternary_percent",
            "ENABLE_TERNARY_AT": "enable_ternary_at",
            "LAMBDA_SPARSE": "lambda_sparse",
            "LAMBDA_TERNARY_END": "lambda_ternary_end",
            "TEMPERATURE_START": "temperature_start",
            "TEMPERATURE_END": "temperature_end",
        }

        # Simple typed overrides
        for env_name, attr in env_map_int.items():
            val = os.getenv(env_name)
            if val is not None:
                try:
                    setattr(self, attr, int(val))
                except ValueError:
                    pass

        for env_name, attr in env_map_float.items():
            val = os.getenv(env_name)
            if val is not None:
                try:
                    setattr(self, attr, float(val))
                except ValueError:
                    pass

        # Learning rate supports LR_INIT or LR
        if os.getenv("LR_INIT") is not None:
            try:
                self.lr_init = float(os.getenv("LR_INIT"))
            except ValueError:
                pass
        elif os.getenv("LR") is not None:
            try:
                self.lr_init = float(os.getenv("LR"))
            except ValueError:
                pass

        # String/enum-like overrides
        opt_env = os.getenv("OPTIMIZER")
        if opt_env is not None:
            self.optimizer = _normalize_optimizer(opt_env)

        device_env = os.getenv("DEVICE")
        if device_env is not None:
            # No validation here; torch will error if invalid when used.
            self.device = device_env

        dtype_env = os.getenv("DTYPE")
        if dtype_env is not None:
            self.dtype = _parse_dtype(dtype_env)

        # Apply CLI overrides if provided
        if _CLI_ARGS.seed is not None:
            self.seed = _CLI_ARGS.seed
        if _CLI_ARGS.lr_init is not None:
            self.lr_init = _CLI_ARGS.lr_init
        if _CLI_ARGS.lambda_ternary_end is not None:
            self.lambda_ternary_end = _CLI_ARGS.lambda_ternary_end
        if _CLI_ARGS.anneal_steps is not None:
            self.anneal_steps = _CLI_ARGS.anneal_steps
        if _CLI_ARGS.enable_ternary_at is not None:
            self.enable_ternary_at = _CLI_ARGS.enable_ternary_at
        if _CLI_ARGS.epochs is not None:
            self.epochs = _CLI_ARGS.epochs
        if _CLI_ARGS.lr_reduction_threshold is not None:
            self.lr_reduction_threshold = _CLI_ARGS.lr_reduction_threshold
        if _CLI_ARGS.lr_reduction_value is not None:
            self.lr_reduction_value = _CLI_ARGS.lr_reduction_value
        if _CLI_ARGS.converge_loss_threshold is not None:
            self.converge_loss_threshold = _CLI_ARGS.converge_loss_threshold
        if _CLI_ARGS.converge_ternary_percent is not None:
            self.converge_ternary_percent = _CLI_ARGS.converge_ternary_percent

    def coeffs(self, step: int, ternary_progress: float):
        """Linearly interpolate λ_T and T once ternary kicks in."""
        progress_ratio = ternary_progress / self.anneal_steps
        λ_T = progress_ratio * self.lambda_ternary_end

        t = min(1.0, progress_ratio)
        T = self.temperature_start * (1 - t) + self.temperature_end * t
        return λ_T, T


