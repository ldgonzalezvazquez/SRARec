#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Install Miniconda and create a SRARec conda environment.

Supported targets:
  --system wsl        -> Linux x86_64 installer, recommended for WSL/Ubuntu
  --system linux      -> Linux x86_64 installer
  --system mac-intel  -> macOS x86_64 installer
  --system mac-arm    -> macOS Apple Silicon arm64 installer

Example:
  python3 install_srarec_env.py --system wsl
  python3 install_srarec_env.py --system mac-intel
  python3 install_srarec_env.py --system mac-arm

The script installs:
  sra-tools, entrez-direct, cutadapt, bwa, minimap2, samtools,
  bcftools, picard, spades, bbmap, shiver, ngsngs
"""

import argparse
import os
import platform
import shutil
import subprocess
import sys
import urllib.request
from pathlib import Path


MINICONDA_URLS = {
    "wsl": "https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh",
    "linux": "https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh",
    "mac-intel": "https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh",
    "mac-arm": "https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh",
}

DEFAULT_INSTALL_DIRS = {
    "wsl": Path.home() / "miniconda3",
    "linux": Path.home() / "miniconda3",
    "mac-intel": Path.home() / "miniconda3",
    "mac-arm": Path.home() / "miniconda3",
}

SRAREC_PACKAGES = [
    "python=3.10",
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
    # Python libraries used by SRARec
    "pysam",
    "tqdm",
    "biopython",
    "pandas",
    "numpy",
]


def run(cmd, *, env=None, check=True):
    print("\n[RUN]", " ".join(str(x) for x in cmd))
    result = subprocess.run(
        [str(x) for x in cmd],
        env=env,
        text=True,
    )
    if check and result.returncode != 0:
        raise RuntimeError(f"Command failed with exit code {result.returncode}: {' '.join(map(str, cmd))}")
    return result


def download_file(url: str, output_path: Path):
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if output_path.exists() and output_path.stat().st_size > 0:
        print(f"[INFO] Installer already exists: {output_path}")
        return

    print(f"[INFO] Downloading Miniconda installer:\n  {url}\n  -> {output_path}")
    urllib.request.urlretrieve(url, output_path)

    if not output_path.exists() or output_path.stat().st_size == 0:
        raise RuntimeError(f"Download failed or produced empty file: {output_path}")


def detect_existing_conda(install_dir: Path):
    conda_in_prefix = install_dir / "bin" / "conda"
    if conda_in_prefix.exists():
        return conda_in_prefix

    conda_from_path = shutil.which("conda")
    if conda_from_path:
        return Path(conda_from_path)

    return None


def install_miniconda(system_name: str, install_dir: Path, force: bool = False):
    if install_dir.exists() and force:
        print(f"[WARN] Removing existing Miniconda directory: {install_dir}")
        shutil.rmtree(install_dir)

    existing_conda = detect_existing_conda(install_dir)
    if existing_conda and not force:
        print(f"[INFO] Conda already found: {existing_conda}")
        return existing_conda

    url = MINICONDA_URLS[system_name]
    installer = Path.cwd() / f"Miniconda3-latest-{system_name}.sh"

    download_file(url, installer)

    install_dir.parent.mkdir(parents=True, exist_ok=True)

    run([
        "bash",
        installer,
        "-b",
        "-p",
        install_dir,
    ])

    conda_bin = install_dir / "bin" / "conda"
    if not conda_bin.exists():
        raise RuntimeError(f"Conda was not found after installation: {conda_bin}")

    return conda_bin


def configure_conda(conda_bin: Path):
    # Do not depend on shell activation. Use conda binary directly.
    run([conda_bin, "config", "--add", "channels", "conda-forge"])
    run([conda_bin, "config", "--add", "channels", "bioconda"])
    run([conda_bin, "config", "--add", "channels", "defaults"])
    run([conda_bin, "config", "--set", "channel_priority", "strict"])


def env_exists(conda_bin: Path, env_name: str):
    result = subprocess.run(
        [str(conda_bin), "env", "list"],
        text=True,
        capture_output=True,
    )
    if result.returncode != 0:
        return False

    for line in result.stdout.splitlines():
        fields = line.split()
        if fields and fields[0] == env_name:
            return True
    return False


def remove_env(conda_bin: Path, env_name: str):
    if env_exists(conda_bin, env_name):
        print(f"[WARN] Removing existing environment: {env_name}")
        run([conda_bin, "env", "remove", "-n", env_name, "-y"])


def create_srarec_env(conda_bin: Path, env_name: str, force_env: bool = False):
    if force_env:
        remove_env(conda_bin, env_name)
    elif env_exists(conda_bin, env_name):
        print(f"[INFO] Environment already exists: {env_name}")
        print("[INFO] Use --force-env to recreate it.")
        return

    cmd = [
        conda_bin,
        "create",
        "-n",
        env_name,
        "-y",
        "-c",
        "conda-forge",
        "-c",
        "bioconda",
        "--strict-channel-priority",
    ] + SRAREC_PACKAGES

    run(cmd)


def verify_tools(conda_bin: Path, env_name: str):
    checks = [
        ["python", "--version"],
        ["prefetch", "--version"],
        ["fasterq-dump", "--version"],
        ["esearch", "-version"],
        ["cutadapt", "--version"],
        ["bwa"],
        ["minimap2", "--version"],
        ["samtools", "--version"],
        ["bcftools", "--version"],
        ["picard", "--version"],
        ["spades.py", "--version"],
        ["repair.sh", "--help"],
        ["shiver_init.sh", "--help"],
        ["ngsngs", "-h"],
    ]

    print("\n[INFO] Verifying installed tools...")
    failed = []

    for check_cmd in checks:
        cmd = [conda_bin, "run", "-n", env_name] + check_cmd
        result = subprocess.run(
            [str(x) for x in cmd],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            text=True,
        )

        tool = check_cmd[0]
        if result.returncode == 0:
            print(f"[OK] {tool}")
        else:
            # NGSNGS sometimes appears with uppercase binary name
            if tool == "ngsngs":
                alt = [conda_bin, "run", "-n", env_name, "NGSNGS", "-h"]
                alt_result = subprocess.run(
                    [str(x) for x in alt],
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL,
                    text=True,
                )
                if alt_result.returncode == 0:
                    print("[OK] NGSNGS")
                    continue

            print(f"[WARN] Could not verify: {tool}")
            failed.append(tool)

    if failed:
        print("\n[WARN] Some tools could not be verified automatically:")
        for tool in failed:
            print(f"  - {tool}")
        print("This does not always mean installation failed; some tools return non-zero for --help.")
    else:
        print("\n[INFO] All tools were verified successfully.")


def print_activation_instructions(install_dir: Path, env_name: str):
    conda_sh = install_dir / "etc" / "profile.d" / "conda.sh"

    print("\n[DONE] Installation finished.")
    print("\nTo activate the environment in this terminal, run:")
    print(f"  source {conda_sh}")
    print(f"  conda activate {env_name}")

    print("\nTo initialize conda permanently for bash, run:")
    print(f"  {install_dir}/bin/conda init bash")
    print("  source ~/.bashrc")

    if platform.system() == "Darwin":
        print("\nOn macOS with zsh, use:")
        print(f"  {install_dir}/bin/conda init zsh")
        print("  source ~/.zshrc")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Install Miniconda and create the SRARec conda environment."
    )

    parser.add_argument(
        "--system",
        required=True,
        choices=["wsl", "linux", "mac-intel", "mac-arm"],
        help="Target system. Use 'wsl' for Ubuntu/WSL, 'mac-intel' for Intel Mac, 'mac-arm' for Apple Silicon."
    )

    parser.add_argument(
        "--env-name",
        default="srarec",
        help="Name of the conda environment to create. Default: srarec"
    )

    parser.add_argument(
        "--install-dir",
        default=None,
        help="Miniconda installation directory. Default: ~/miniconda3"
    )

    parser.add_argument(
        "--force-conda",
        action="store_true",
        help="Remove and reinstall Miniconda if the installation directory already exists."
    )

    parser.add_argument(
        "--force-env",
        action="store_true",
        help="Remove and recreate the conda environment if it already exists."
    )

    parser.add_argument(
        "--skip-conda-install",
        action="store_true",
        help="Do not install Miniconda; use existing conda from PATH or --install-dir."
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
    install_dir = Path(args.install_dir).expanduser().resolve() if args.install_dir else DEFAULT_INSTALL_DIRS[system_name]

    if system_name == "mac-arm":
        print(
            "[WARN] You selected mac-arm. Most packages should resolve through Bioconda/conda-forge, "
            "but some older bioinformatics tools may be easier to install on WSL/Linux or macOS Intel/Rosetta."
        )

    if args.skip_conda_install:
        conda_bin = detect_existing_conda(install_dir)
        if not conda_bin:
            raise RuntimeError(
                "No existing conda found. Remove --skip-conda-install or provide --install-dir."
            )
    else:
        conda_bin = install_miniconda(
            system_name=system_name,
            install_dir=install_dir,
            force=args.force_conda,
        )

    configure_conda(conda_bin)
    create_srarec_env(conda_bin, args.env_name, force_env=args.force_env)

    if not args.skip_verify:
        verify_tools(conda_bin, args.env_name)

    print_activation_instructions(install_dir, args.env_name)


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\n[ERROR] Interrupted by user.")
        sys.exit(130)
    except Exception as exc:
        print(f"\n[ERROR] {exc}")
        sys.exit(1)
