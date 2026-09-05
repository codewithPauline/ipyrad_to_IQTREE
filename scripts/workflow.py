#!/usr/bin/env python3
"""Portable entry points for pre-demultiplexed single-end RADseq analysis."""
import argparse
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path


def positive(value):
    number = int(value)
    if number < 1:
        raise argparse.ArgumentTypeError("must be a positive integer")
    return number


def existing(value):
    path = Path(value).expanduser().resolve()
    if not path.is_file() or path.stat().st_size == 0:
        raise argparse.ArgumentTypeError(f"missing or empty file: {path}")
    return str(path)


def phylip(path):
    """Validate relaxed sequential DNA PHYLIP (one complete sequence per line)."""
    with open(path) as handle:
        header = handle.readline().split()
        if len(header) != 2 or not all(x.isdigit() for x in header):
            raise ValueError("PHYLIP header must contain sample and site counts")
        samples, sites = map(int, header)
        if samples < 4 or sites < 1:
            raise ValueError("at least four samples and one site are required")
        names = set()
        for line in handle:
            if not line.strip():
                continue
            fields = line.split()
            if len(fields) < 2:
                raise ValueError("expected a sample ID and complete sequence on each line")
            name, sequence = fields[0], ''.join(fields[1:]).upper()
            if name in names:
                raise ValueError(f"duplicate sample ID: {name}")
            if len(sequence) != sites or re.search(r"[^ACGTRYSWKMBDHVN?\-]", sequence):
                raise ValueError(f"invalid length or DNA symbols for {name}")
            if not any(base in "ACGT" for base in sequence):
                raise ValueError(f"no resolved nucleotide bases for {name}")
            names.add(name)
        if len(names) != samples:
            raise ValueError(f"header declares {samples} samples; found {len(names)}")
    return {"samples": samples, "sites": sites}


def run(command, dry_run):
    print(json.dumps({"command": command}), flush=True)
    if dry_run:
        return
    if not shutil.which(command[0]):
        raise ValueError(f"{command[0]} not found; activate the analysis environment")
    subprocess.run(command, check=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="stage", required=True)
    check = sub.add_parser("check", help="validate relaxed sequential PHYLIP")
    check.add_argument("alignment", type=existing)
    qc = sub.add_parser("qc", help="FastQC on demultiplexed FASTQ files")
    qc.add_argument("reads", nargs="+", type=existing)
    qc.add_argument("--out", required=True)
    trim = sub.add_parser("trim", help="fastp on one demultiplexed single-end sample")
    trim.add_argument("reads", type=existing)
    trim.add_argument("--out", required=True, help="new output directory for this sample")
    trim.add_argument("--adapter", required=True, help="known adapter sequence, from library protocol")
    trim.add_argument("--min-length", type=positive, default=35)
    assembly = sub.add_parser("assemble", help="run ipyrad with a user-reviewed parameter file")
    assembly.add_argument("params", type=existing)
    assembly.add_argument("--steps", default="1234567")
    infer = sub.add_parser("infer", help="IQ-TREE model selection and branch support")
    infer.add_argument("alignment", type=existing)
    infer.add_argument("--out", required=True, help="IQ-TREE output prefix")
    infer.add_argument("--model", default="MFP")
    infer.add_argument("--seed", type=positive, default=2026)
    infer.add_argument("--binary", default="iqtree3", help="IQ-TREE 3 executable or path")
    for entry in (qc, trim, assembly, infer):
        entry.add_argument("--threads", type=positive, default=1)
        entry.add_argument("--dry-run", action="store_true", help="validate inputs and print command without running")
    args = parser.parse_args()
    if args.stage == "check":
        print(json.dumps(phylip(args.alignment)))
        return
    if args.stage == "qc":
        out = Path(args.out).expanduser().resolve()
        if out.exists() and any(out.iterdir()):
            raise ValueError("QC output directory must be empty or new")
        basenames = [Path(x).name for x in args.reads]
        if len(set(basenames)) != len(basenames):
            raise ValueError("FASTQ basenames must be unique to avoid report collisions")
        command = ["fastqc", "--threads", str(args.threads), "--outdir", str(out), *args.reads]
        if not args.dry_run:
            out.mkdir(parents=True, exist_ok=True)
    elif args.stage == "trim":
        if re.fullmatch(r"[ACGTacgt]+", args.adapter) is None:
            raise ValueError("adapter must be a DNA sequence containing only A, C, G, T")
        if args.threads > 16:
            raise ValueError("fastp supports at most 16 threads")
        out = Path(args.out).expanduser().resolve()
        if out.exists() and any(out.iterdir()):
            raise ValueError("trim output directory must be empty or new")
        command = ["fastp", "--in1", args.reads, "--out1", str(out / "reads.fastq.gz"),
                   "--html", str(out / "fastp.html"), "--json", str(out / "fastp.json"),
                   "--thread", str(args.threads), "--adapter_sequence", args.adapter.upper(),
                   "--length_required", str(args.min_length)]
        if not args.dry_run:
            out.mkdir(parents=True, exist_ok=True)
    elif args.stage == "assemble":
        if re.fullmatch(r"1?2?3?4?5?6?7?", args.steps) is None or not args.steps:
            raise ValueError("steps must be an ordered, nonrepeating subset of 1234567")
        text = Path(args.params).read_text()
        if "REPLACE_" in text or "/shared/jezkovt_shared/" in text:
            raise ValueError("replace historical paths and review parameters before running")
        command = ["ipyrad", "-p", args.params, "-s", args.steps, "-c", str(args.threads)]
    else:
        phylip(args.alignment)
        prefix = Path(args.out).expanduser().resolve()
        if prefix == Path(args.alignment):
            raise ValueError("use an output prefix separate from the input alignment")
        # IQ-TREE handles checkpoint resumption and refuses completed runs; never force redo.
        command = [args.binary, "-s", args.alignment, "-st", "DNA", "-m", args.model,
                   "-B", "1000", "--alrt", "1000", "-T", str(args.threads),
                   "--prefix", str(prefix), "-seed", str(args.seed)]
        if not args.dry_run:
            prefix.parent.mkdir(parents=True, exist_ok=True)
    run(command, args.dry_run)


if __name__ == "__main__":
    try:
        main()
    except (ValueError, OSError, subprocess.CalledProcessError) as error:
        print(f"Error: {error}", file=sys.stderr)
        sys.exit(1)
