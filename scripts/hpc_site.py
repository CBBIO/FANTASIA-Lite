"""Validation and path mapping for a FANTASIA Lite HPC site configuration."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path, PurePosixPath
from typing import Any, Dict, Iterable


class SiteConfigError(ValueError):
    """Raised when an HPC site YAML file is incomplete or unsafe."""


@dataclass(frozen=True)
class Mount:
    name: str
    host: Path
    container: PurePosixPath
    mode: str

    def bind_argument(self, *, mode: str | None = None) -> str:
        return f"{self.host}:{self.container}:{mode or self.mode}"


@dataclass(frozen=True)
class SiteConfig:
    image: Path
    apptainer: str
    apptainer_cache: Path
    use_nv: bool
    mounts: Dict[str, Mount]
    lookup_url: str | None
    lookup_sha256: str | None
    lookup_archive_name: str

    def mount(self, name: str) -> Mount:
        try:
            return self.mounts[name]
        except KeyError as exc:
            raise SiteConfigError(f"Missing storage.{name} in site configuration.") from exc


REQUIRED_MOUNTS = ("lookup", "models", "inputs", "outputs", "scratch", "staging")


def _as_mapping(value: Any, location: str) -> Dict[str, Any]:
    if not isinstance(value, dict):
        raise SiteConfigError(f"{location} must be a YAML mapping.")
    return value


def _absolute_path(value: Any, location: str) -> Path:
    if not isinstance(value, str) or not value.strip():
        raise SiteConfigError(f"{location} must be a non-empty absolute path.")
    path = Path(value).expanduser()
    if not path.is_absolute():
        raise SiteConfigError(f"{location} must be absolute: {value!r}")
    return path


def _container_path(value: Any, location: str) -> PurePosixPath:
    if not isinstance(value, str) or not value.startswith("/"):
        raise SiteConfigError(f"{location} must be an absolute container path.")
    path = PurePosixPath(value)
    if ".." in path.parts:
        raise SiteConfigError(f"{location} cannot contain '..'.")
    return path


def load_site_config(path: Path) -> SiteConfig:
    """Read and validate the portable YAML contract used by the HPC launcher."""
    try:
        import yaml
    except ImportError as exc:  # pragma: no cover - host installation issue
        raise SiteConfigError(
            "PyYAML is required by scripts/fantasia-hpc. Install it in the host "
            "Python used to launch Apptainer."
        ) from exc

    try:
        raw = yaml.safe_load(path.read_text(encoding="utf-8"))
    except OSError as exc:
        raise SiteConfigError(f"Cannot read site configuration {path}: {exc}") from exc
    except yaml.YAMLError as exc:
        raise SiteConfigError(f"Invalid YAML in {path}: {exc}") from exc
    raw = _as_mapping(raw, "site configuration")
    runtime = _as_mapping(raw.get("runtime"), "runtime")
    storage = _as_mapping(raw.get("storage"), "storage")
    assets = _as_mapping(raw.get("assets", {}), "assets")
    lookup = _as_mapping(assets.get("lookup", {}), "assets.lookup")

    mounts: Dict[str, Mount] = {}
    container_paths: set[PurePosixPath] = set()
    for name in REQUIRED_MOUNTS:
        item = _as_mapping(storage.get(name), f"storage.{name}")
        host = _absolute_path(item.get("host"), f"storage.{name}.host")
        container = _container_path(item.get("container"), f"storage.{name}.container")
        mode = item.get("mode")
        if mode not in {"ro", "rw"}:
            raise SiteConfigError(f"storage.{name}.mode must be 'ro' or 'rw'.")
        if container in container_paths:
            raise SiteConfigError(f"Container mount path {container} is repeated.")
        container_paths.add(container)
        mounts[name] = Mount(name, host, container, mode)

    for name in ("lookup", "models", "inputs"):
        if mounts[name].mode != "ro":
            raise SiteConfigError(f"storage.{name}.mode must be 'ro' for compute jobs.")
    for name in ("outputs", "scratch"):
        if mounts[name].mode != "rw":
            raise SiteConfigError(f"storage.{name}.mode must be 'rw' for compute jobs.")

    image = _absolute_path(runtime.get("image"), "runtime.image")
    cache = _absolute_path(runtime.get("apptainer_cache"), "runtime.apptainer_cache")
    command = runtime.get("apptainer", "apptainer")
    if not isinstance(command, str) or not command.strip():
        raise SiteConfigError("runtime.apptainer must be a command name or path.")
    use_nv = runtime.get("use_nv", True)
    if not isinstance(use_nv, bool):
        raise SiteConfigError("runtime.use_nv must be true or false.")

    lookup_url = lookup.get("url")
    lookup_sha256 = lookup.get("sha256")
    archive_name = lookup.get("archive_name", "fantasia_lite_data_folder.tgz")
    if lookup_url is not None and not isinstance(lookup_url, str):
        raise SiteConfigError("assets.lookup.url must be a string.")
    if lookup_sha256 is not None and (
        not isinstance(lookup_sha256, str) or len(lookup_sha256) != 64
    ):
        raise SiteConfigError("assets.lookup.sha256 must be a 64-character SHA-256.")
    if not isinstance(archive_name, str) or not archive_name:
        raise SiteConfigError("assets.lookup.archive_name must be a non-empty filename.")

    return SiteConfig(
        image=image,
        apptainer=command,
        apptainer_cache=cache,
        use_nv=use_nv,
        mounts=mounts,
        lookup_url=lookup_url,
        lookup_sha256=lookup_sha256,
        lookup_archive_name=archive_name,
    )


def mapped_input(site: SiteConfig, input_path: Path) -> PurePosixPath:
    """Map an input below storage.inputs.host to its container destination."""
    mount = site.mount("inputs")
    candidate = input_path.expanduser().resolve()
    try:
        relative = candidate.relative_to(mount.host.resolve())
    except ValueError as exc:
        raise SiteConfigError(
            f"Input {candidate} is outside configured storage.inputs.host ({mount.host})."
        ) from exc
    return mount.container / PurePosixPath(relative.as_posix())


def run_mounts(site: SiteConfig, run_name: str) -> Iterable[Mount]:
    """Yield compute mounts, isolating writable files below a run name."""
    if not run_name or Path(run_name).name != run_name:
        raise SiteConfigError("run name must be a single relative path component.")
    for name in ("lookup", "models", "inputs"):
        yield site.mount(name)
    for name in ("outputs", "scratch"):
        mount = site.mount(name)
        yield Mount(name, mount.host / run_name, mount.container, "rw")
