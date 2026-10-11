#!/opt/fantasia/venv/bin/python
"""Validate the legacy installer's requests against the immutable runtime."""
import sys
from importlib.metadata import PackageNotFoundError, version
from packaging.requirements import Requirement


def main():
    args = sys.argv[1:]
    if not args or args.pop(0) != "install":
        raise RuntimeError("The container pip adapter only supports install validation.")
    while args:
        item = args.pop(0)
        if item in {"--upgrade", "--no-deps"}:
            continue
        if item == "--index-url":
            args.pop(0)
            continue
        requirement = Requirement(item)
        installed = version(requirement.name)
        if requirement.specifier and installed not in requirement.specifier:
            raise RuntimeError(f"Preinstalled {requirement.name}=={installed} does not satisfy {item}")
    print("Container dependencies verified; no packages installed.")


if __name__ == "__main__":
    try:
        main()
    except (RuntimeError, PackageNotFoundError, ValueError, IndexError) as exc:
        sys.exit(str(exc))
