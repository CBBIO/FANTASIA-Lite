"""From the repo root: apptainer exec --nv IMAGE python containers/smoke_gpu.py."""
import json
import torch
from transformers import EsmModel, T5EncoderModel, T5Tokenizer
from esm.models.esmc import ESMC

assert torch.cuda.is_available(), "No CUDA GPU is visible inside the container"
torch.manual_seed(123)
left = torch.randn(256, 256)
right = torch.randn(256, 256)
actual = (left.cuda() @ right.cuda()).cpu()
torch.testing.assert_close(actual, left @ right, atol=1e-4, rtol=1e-4)
torch.cuda.synchronize()
print(json.dumps({
    "torch": torch.__version__,
    "cuda": torch.version.cuda,
    "gpu": torch.cuda.get_device_name(0),
    "capability": torch.cuda.get_device_capability(0),
    "compiled_architectures": torch.cuda.get_arch_list(),
    "matrix_multiply": "passed",
    "model_runtime_imports": "passed",
}, indent=2))
