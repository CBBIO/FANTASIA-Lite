"""Check dependencies, allowing only ESM's unused vision/text requirements."""
import subprocess
import sys

result = subprocess.run([sys.executable, "-m", "pip", "check"], capture_output=True, text=True)
allowed = {
    "esm 3.2.3 requires torchvision, which is not installed.",
    "esm 3.2.3 requires torchtext, which is not installed.",
}
issues = set(result.stdout.strip().splitlines())
if result.returncode and (not issues or issues - allowed):
    sys.exit(result.stdout + result.stderr)
print(result.stdout)
