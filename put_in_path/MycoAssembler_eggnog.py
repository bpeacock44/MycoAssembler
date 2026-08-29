#!/usr/bin/env python3
import os
import shutil
import subprocess
import argparse
from packaging import version

# --- Helpers to replicate funannotate ---
def get_emapper_version():
    """Return Eggnog-mapper version as string."""
    try:
        out = subprocess.check_output(["emapper.py", "--version"], universal_newlines=True)
        return out.strip().split()[0].replace("emapper-", "")
    except Exception:
        return "0.0.0"

def get_diamond_version():
    """Return Diamond version as string."""
    try:
        out = subprocess.check_output(["diamond", "--version"], universal_newlines=True)
        return out.strip().split()[-1]  # usually last token
    except Exception:
        return "0.0.0"

def memory_gb():
    """Return total system memory in GB (rough approximation)."""
    import psutil
    return psutil.virtual_memory().total / 1e9

# --- Parse arguments ---
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", required=True, help="Protein FASTA file")
parser.add_argument("-o", "--outdir", required=True, help="Output directory")
parser.add_argument("--cpus", type=int, default=8)
parser.add_argument("--tmpdir", default="tmp_eggnog")
parser.add_argument("--scratch_dir", default="scratch_eggnog")
args = parser.parse_args()

os.makedirs(args.outdir, exist_ok=True)

# --- Build command ---
cmd = [
    "emapper.py",
    "-m", "diamond",
    "-i", args.input,
    "-o", args.outdir,
    "--cpu", str(args.cpus),
    "--override"
]

emap_ver = version.parse(get_emapper_version())
diamond_ver = version.parse(get_diamond_version())

if emap_ver >= version.parse("2.1.0"):
    os.makedirs(args.tmpdir, exist_ok=True)
    os.makedirs(args.scratch_dir, exist_ok=True)
    cmd += ["--scratch_dir", args.scratch_dir, "--temp_dir", args.tmpdir]
    if memory_gb() >= 48:
        cmd.append("--dbmem")

if emap_ver >= version.parse("2.1.4"):
    if diamond_ver < version.parse("2.0.11"):
        cmd += ["--dmnd_iterate", "no"]

# --- Run Eggnog-mapper ---
print("Running command:", " ".join(cmd))
subprocess.run(cmd, check=True)

# --- Cleanup ---
if os.path.isdir(args.scratch_dir):
    shutil.rmtree(args.scratch_dir)
