#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Install a conda-compatible package manager and create the SRARec environment.

Supported package managers (auto-detected, or forced with --manager):
  micromamba, mamba, conda
If none is found, Miniconda is installed (or micromamba with --manager micromamba).

Supported targets:
  --system wsl        -> Linux x86_64, recommended for WSL/Ubuntu
  --system linux      -> Linux x86_64
  --system mac-intel  -> macOS x86_64 (also the recommended choice for Apple Silicon, via Rosetta 2)
  --system mac-arm    -> macOS Apple Silicon arm64 (native; ngsngs and smalt/shiver have no osx-arm64 build)

Examples:
  python3 install_srarec_env.py --system wsl
  python3 install_srarec_env.py --system mac-intel          # Intel Macs AND Apple Silicon (Rosetta 2)
  python3 install_srarec_env.py --system linux --manager micromamba

The environment contains:
  sra-tools, entrez-direct, cutadapt, bwa, minimap2, samtools, bcftools,
  picard, spades, bbmap, shiver, ngsngs and the Python libraries used by
  SRARec1.3.py / generate_graphs.py (pysam, tqdm, biopython, pandas, numpy,
  matplotlib, plotly, python-kaleido).
"""

import argparse
import os
import platform
import shutil
import subprocess
import sys
import tarfile
import urllib.request
from pathlib import Path


MINICONDA_URLS = {
    "wsl": "https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh",
    "linux": "https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh",
    "mac-intel": "https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh",
    "mac-arm": "https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh",
}

# Platform string used by https://micro.mamba.pm/api/micromamba/<platform>/latest
MICROMAMBA_PLATFORMS = {
    "wsl": "linux-64",
    "linux": "linux-64",
    "mac-intel": "osx-64",
    "mac-arm": "osx-arm64",
}

DEFAULT_INSTALL_DIRS = {
    "miniconda": Path.home() / "miniconda3",
    "micromamba": Path.home() / "micromamba",
}

# Executable names probed for each manager, in order of preference when --manager auto.
MANAGERS = ("micromamba", "mamba", "conda")

# Channels in priority order (conda-forge first). 'defaults' is intentionally left out:
# mixing it with conda-forge/bioconda is a common source of solver conflicts.
CHANNELS = ["conda-forge", "bioconda"]

CONTACT_MESSAGE = (
    "If you have problems installing the dependencies, contact us and we will help you:\n"
    "  - Luis Daniel González-Vázquez: luisdaniel.gonzalez@uvigo.es\n"
    "  - Paula Iglesias-Rivas:         paula.iglesias.rivas@uvigo.gal\n"
    "Please include your operating system, the package manager used and the full error message."
)


def print_contact_message():
    print(f"\n[HELP] {CONTACT_MESSAGE}")


SRAREC_PACKAGES = [
    "python=3.10",
    # External tools
    "sra-tools",
    "entrez-direct",
    "cutadapt",
    "bwa",
    "minimap2",
    "samtools",
    "bcftools",
    "picard",
    "spades",
    "bbmap",
    "shiver",
    "ngsngs",
    # Python libraries used by SRARec1.3.py and generate_graphs.py
    "pysam",
    "tqdm",
    "biopython",
    "pandas",
    "numpy",
    "matplotlib",
    "plotly",
    "python-kaleido",
]


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #

def run(cmd, *, env=None, check=True):
    print("\n[RUN]", " ".join(str(x) for x in cmd))
    result = subprocess.run([str(x) for x in cmd], env=env, text=True)
    if check and result.returncode != 0:
        raise RuntimeError(f"Command failed with exit code {result.returncode}: {' '.join(map(str, cmd))}")
    return result


def manager_kind(manager_bin: Path) -> str:
    """Return 'micromamba', 'mamba' or 'conda' from the executable name."""
    name = Path(manager_bin).name.lower()
    for kind in MANAGERS:
        if name.startswith(kind):
            return kind
    return "conda"


def is_micromamba(manager_bin: Path) -> bool:
    return manager_kind(manager_bin) == "micromamba"


def manager_env(manager_bin: Path, subdir=None):
    """Environment for subprocesses. micromamba needs a root prefix; Rosetta needs CONDA_SUBDIR."""
    env = os.environ.copy()
    default_root = DEFAULT_INSTALL_DIRS["micromamba"]
    if is_micromamba(manager_bin) and "MAMBA_ROOT_PREFIX" not in env and default_root in Path(manager_bin).parents:
        env["MAMBA_ROOT_PREFIX"] = str(default_root)  # micromamba installed by this script
    if subdir:
        env["CONDA_SUBDIR"] = subdir
    return env


def download_file(url: str, output_path: Path):
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if output_path.exists() and output_path.stat().st_size > 0:
        print(f"[INFO] File already exists: {output_path}")
        return

    print(f"[INFO] Downloading:\n  {url}\n  -> {output_path}")
    urllib.request.urlretrieve(url, output_path)

    if not output_path.exists() or output_path.stat().st_size == 0:
        raise RuntimeError(f"Download failed or produced empty file: {output_path}")


# --------------------------------------------------------------------------- #
# Manager detection / installation
# --------------------------------------------------------------------------- #

def detect_existing_manager(install_dir: Path, preferred: str = "auto"):
    """Find micromamba / mamba / conda. Looks in install_dir first, then in PATH."""
    candidates = MANAGERS if preferred == "auto" else (preferred,)

    for name in candidates:
        in_prefix = install_dir / "bin" / name
        if in_prefix.exists():
            return in_prefix
        from_path = shutil.which(name)
        if from_path:
            return Path(from_path)

    return None


def install_miniconda(system_name: str, install_dir: Path, force: bool = False):
    if install_dir.exists() and force:
        print(f"[WARN] Removing existing Miniconda directory: {install_dir}")
        shutil.rmtree(install_dir)

    url = MINICONDA_URLS[system_name]
    installer = Path.cwd() / f"Miniconda3-latest-{system_name}.sh"
    download_file(url, installer)

    install_dir.parent.mkdir(parents=True, exist_ok=True)
    run(["bash", installer, "-b", "-p", install_dir])

    conda_bin = install_dir / "bin" / "conda"
    if not conda_bin.exists():
        raise RuntimeError(f"Conda was not found after installation: {conda_bin}")

    return conda_bin


def install_micromamba(system_name: str, install_dir: Path, force: bool = False):
    if install_dir.exists() and force:
        print(f"[WARN] Removing existing micromamba directory: {install_dir}")
        shutil.rmtree(install_dir)

    plat = MICROMAMBA_PLATFORMS[system_name]
    url = f"https://micro.mamba.pm/api/micromamba/{plat}/latest"
    archive = Path.cwd() / f"micromamba-{plat}.tar.bz2"
    download_file(url, archive)

    install_dir.mkdir(parents=True, exist_ok=True)
    with tarfile.open(archive, "r:bz2") as tar:
        member = tar.getmember("bin/micromamba")
        member.name = "micromamba"  # extract only the binary, no path traversal
        (install_dir / "bin").mkdir(parents=True, exist_ok=True)
        tar.extract(member, path=install_dir / "bin")

    micromamba_bin = install_dir / "bin" / "micromamba"
    micromamba_bin.chmod(0o755)
    if not micromamba_bin.exists():
        raise RuntimeError(f"micromamba was not found after installation: {micromamba_bin}")

    return micromamba_bin


def ensure_manager(system_name, install_dir, preferred, force, skip_install):
    if not force:
        existing = detect_existing_manager(install_dir, preferred)
        if existing:
            print(f"[INFO] Using existing package manager: {existing}")
            return existing

    if skip_install:
        raise RuntimeError(
            "No existing micromamba/mamba/conda found. Remove --skip-conda-install or provide --install-dir."
        )

    if preferred in ("auto", "conda"):
        return install_miniconda(system_name, install_dir, force=force)
    if preferred == "micromamba":
        return install_micromamba(system_name, install_dir, force=force)

    raise RuntimeError(
        f"'{preferred}' was not found and cannot be installed automatically. "
        "Install it first, or use --manager auto/conda/micromamba."
    )


# --------------------------------------------------------------------------- #
# Configuration and environment handling
# --------------------------------------------------------------------------- #

def configure_conda(conda_bin: Path, env=None):
    """Configure channels and a flexible channel priority.

    micromamba and mamba 2.x: `config append channels X` (they do not accept --add).
    conda and mamba 1.x:      `config --add channels X`, which puts X at the TOP of the list, so the
                              channels are added in reverse order to end with conda-forge first.
    For mamba the syntax depends on the version, so it tries one style and falls back to the other.
    Flexible priority avoids solver deadlocks when old Bioconda builds meet recent conda-forge libs.
    """
    def append_style():
        for channel in CHANNELS:
            run([conda_bin, "config", "append", "channels", channel], env=env, check=False)
        run([conda_bin, "config", "set", "channel_priority", "flexible"], env=env)

    def add_style(strict=True):
        for i, channel in enumerate(reversed(CHANNELS)):
            result = run([conda_bin, "config", "--add", "channels", channel], env=env, check=strict)
            if result.returncode != 0:
                return False
        run([conda_bin, "config", "--set", "channel_priority", "flexible"], env=env, check=strict)
        return True

    kind = manager_kind(conda_bin)
    if kind == "micromamba":
        append_style()
    elif kind == "mamba":
        if not add_style(strict=False):
            print("[INFO] This mamba does not accept 'config --add'; using 'config append' (mamba 2.x).")
            append_style()
    else:
        add_style()


def env_exists(conda_bin: Path, env_name: str, env=None):
    result = subprocess.run(
        [str(conda_bin), "env", "list"],
        text=True,
        capture_output=True,
        env=env,
    )
    if result.returncode != 0:
        return False

    for line in result.stdout.splitlines():
        fields = line.split()
        # conda and micromamba both print the environment name first ("name  [*]  /path")
        if fields and fields[0] == env_name:
            return True
    return False


def remove_env(conda_bin: Path, env_name: str, env=None):
    if env_exists(conda_bin, env_name, env):
        print(f"[WARN] Removing existing environment: {env_name}")
        run([conda_bin, "env", "remove", "-n", env_name, "-y"], env=env)


def create_srarec_env(conda_bin: Path, env_name: str, force_env: bool = False, env=None):
    if force_env:
        remove_env(conda_bin, env_name, env)
    elif env_exists(conda_bin, env_name, env):
        print(f"[INFO] Environment already exists: {env_name}")
        print("[INFO] Use --force-env to recreate it.")
        return False

    cmd = [conda_bin, "create", "-n", env_name, "-y", "--override-channels"]
    for channel in CHANNELS:
        cmd += ["-c", channel]
    cmd += SRAREC_PACKAGES

    run(cmd, env=env)
    return True


def pin_env_subdir(conda_bin: Path, env_name: str, subdir: str, env=None):
    """Persist the platform (e.g. osx-64 under Rosetta) inside the environment's own .condarc,
    so later `conda install` calls in that environment keep resolving for the same architecture."""
    result = subprocess.run(
        [str(conda_bin), "run", "-n", env_name, "python", "-c", "import sys; print(sys.prefix)"],
        text=True, capture_output=True, env=env,
    )
    prefix = result.stdout.strip().splitlines()[-1] if result.returncode == 0 and result.stdout.strip() else None
    if not prefix:
        print(f"[WARN] Could not locate the environment prefix; run this manually after activating it:\n"
              f"  conda config --env --set subdir {subdir}")
        return
    condarc = Path(prefix) / ".condarc"
    condarc.write_text(f"subdir: {subdir}\n")
    print(f"[INFO] Environment pinned to {subdir} ({condarc})")


def verify_tools(conda_bin: Path, env_name: str, env=None):
    checks = [
        ["python", "--version"],
        ["python", "-c", "import pysam, tqdm, Bio, pandas, numpy, matplotlib, plotly"],
        ["python", "-c", "import kaleido"],
        ["prefetch", "--version"],
        ["fasterq-dump", "--version"],
        ["esearch", "-version"],
        ["cutadapt", "--version"],
        ["which", "bwa"],  # bwa/picard/ngsngs return exit code 1 on --help: just check they exist
        ["minimap2", "--version"],
        ["samtools", "--version"],
        ["bcftools", "--version"],
        ["which", "picard"],
        ["spades.py", "--version"],
        ["repair.sh", "--help"],
        ["shiver_init.sh", "--help"],
        ["which", "ngsngs"],
    ]

    print("\n[INFO] Verifying installed tools...")
    failed = []

    def ok(check_cmd):
        result = subprocess.run(
            [str(x) for x in [conda_bin, "run", "-n", env_name] + check_cmd],
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, text=True, env=env,
        )
        return result.returncode == 0

    for check_cmd in checks:
        if check_cmd[0] == "python" and len(check_cmd) >= 3:
            tool = f"python ({check_cmd[2]})"
        elif check_cmd[0] == "which":
            tool = check_cmd[1]
        else:
            tool = check_cmd[0]

        if ok(check_cmd):
            print(f"[OK] {tool}")
        elif tool == "ngsngs" and ok(["which", "NGSNGS"]):  # sometimes uppercase binary
            print("[OK] NGSNGS")
        else:
            print(f"[WARN] Could not verify: {tool}")
            failed.append(tool)

    if failed:
        print("\n[WARN] Some tools could not be verified automatically:")
        for tool in failed:
            print(f"  - {tool}")
        print("This does not always mean installation failed; some tools return non-zero for --help.")
        print_contact_message()
    else:
        print("\n[INFO] All tools were verified successfully.")


def print_activation_instructions(manager_bin: Path, env_name: str, subdir=None):
    kind = manager_kind(manager_bin)

    print("\n[DONE] Installation finished.")

    if kind == "micromamba":
        root = os.environ.get("MAMBA_ROOT_PREFIX") or str(
            DEFAULT_INSTALL_DIRS["micromamba"]
            if DEFAULT_INSTALL_DIRS["micromamba"] in Path(manager_bin).parents
            else Path.home() / ".local" / "share" / "mamba"
        )
        print("\nTo activate the environment in this terminal, run:")
        print(f'  eval "$({manager_bin} shell hook --shell bash)"   # zsh: --shell zsh')
        print(f"  micromamba activate {env_name}")
        print("\nTo initialize micromamba permanently, run:")
        print(f"  {manager_bin} shell init --shell bash --root-prefix {root}   # zsh: --shell zsh")
    else:
        # <prefix>/bin/conda -> <prefix>; for conda in PATH this is still the right base.
        base = Path(manager_bin).resolve().parent.parent
        conda_sh = base / "etc" / "profile.d" / "conda.sh"
        print("\nTo activate the environment in this terminal, run:")
        if conda_sh.exists():
            print(f"  source {conda_sh}")
        print(f"  {'mamba' if kind == 'mamba' else 'conda'} activate {env_name}")
        print("\nTo initialize conda permanently for bash, run:")
        print(f"  {manager_bin} init bash")
        print("  source ~/.bashrc")
        if platform.system() == "Darwin":
            print("\nOn macOS with zsh, use:")
            print(f"  {manager_bin} init zsh")
            print("  source ~/.zshrc")

    if subdir:
        print(f"\nThis environment targets {subdir} (Rosetta 2) and is pinned to it, so later installs")
        print("inside it keep resolving for the same architecture.")


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #

def parse_args():
    parser = argparse.ArgumentParser(
        description="Install a conda-compatible manager and create the SRARec environment."
    )

    parser.add_argument(
        "--system",
        required=True,
        choices=["wsl", "linux", "mac-intel", "mac-arm"],
        help="Target system. 'wsl' for Ubuntu/WSL, 'mac-intel' for Intel Macs and (recommended) "
             "Apple Silicon via Rosetta 2, 'mac-arm' for native Apple Silicon (some packages unavailable)."
    )

    parser.add_argument(
        "--manager",
        choices=["auto", "conda", "mamba", "micromamba"],
        default="auto",
        help="Package manager to use. Default: auto (micromamba > mamba > conda; installs Miniconda if none is found)."
    )

    parser.add_argument(
        "--env-name", "--env_name",
        dest="env_name",
        default="srarec",
        help="Name of the environment to create. Default: srarec"
    )

    parser.add_argument(
        "--install-dir",
        default=None,
        help="Installation directory for a newly installed manager. Default: ~/miniconda3 (~/micromamba for micromamba)"
    )

    parser.add_argument(
        "--force-conda",
        action="store_true",
        help="Remove and reinstall the manager if the installation directory already exists."
    )

    parser.add_argument(
        "--force-env",
        action="store_true",
        help="Remove and recreate the environment if it already exists."
    )

    parser.add_argument(
        "--skip-conda-install",
        action="store_true",
        help="Do not install anything; use an existing micromamba/mamba/conda from PATH or --install-dir."
    )

    parser.add_argument(
        "--skip-verify",
        action="store_true",
        help="Skip tool verification at the end."
    )

    return parser.parse_args()


def main():
    args = parse_args()

    system_name = args.system
    default_dir = DEFAULT_INSTALL_DIRS["micromamba" if args.manager == "micromamba" else "miniconda"]
    install_dir = Path(args.install_dir).expanduser().resolve() if args.install_dir else default_dir

    # Apple Silicon host running the Intel (osx-64) target -> resolve packages for osx-64 under Rosetta 2.
    subdir = None
    if system_name == "mac-intel" and platform.system() == "Darwin" and platform.machine() == "arm64":
        subdir = "osx-64"
        print(
            "[INFO] Apple Silicon detected with --system mac-intel: using osx-64 (Rosetta 2). "
            "If Rosetta is not installed, run: softwareupdate --install-rosetta --agree-to-license"
        )

    if system_name == "mac-arm":
        print(
            "[WARN] You selected mac-arm (native osx-arm64). ngsngs and smalt (a SHIVER dependency) have no "
            "osx-arm64 build, so the environment will probably fail to solve. Use '--system mac-intel' "
            "instead: it works on Apple Silicon through Rosetta 2."
        )

    manager_bin = ensure_manager(
        system_name=system_name,
        install_dir=install_dir,
        preferred=args.manager,
        force=args.force_conda,
        skip_install=args.skip_conda_install,
    )

    env = manager_env(manager_bin, subdir)

    configure_conda(manager_bin, env)
    created = create_srarec_env(manager_bin, args.env_name, force_env=args.force_env, env=env)

    if created and subdir:
        pin_env_subdir(manager_bin, args.env_name, subdir, env)

    if not args.skip_verify:
        verify_tools(manager_bin, args.env_name, env)

    print_activation_instructions(manager_bin, args.env_name, subdir)


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\n[ERROR] Interrupted by user.")
        sys.exit(130)
    except Exception as exc:
        print(f"\n[ERROR] {exc}")
        print_contact_message()
        sys.exit(1)
