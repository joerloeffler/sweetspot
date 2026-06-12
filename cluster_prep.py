#!/usr/bin/env python3
import argparse
from collections import OrderedDict
from pathlib import Path
import yaml


SNAKEFILE_TEXT = r"""
import os

JOBDIR = os.path.dirname(workflow.snakefile)
configfile: os.path.join(JOBDIR, "config.yaml")

FAA = os.path.join(JOBDIR, "msa_queries.fasta")
YAML = os.path.join(JOBDIR, "input.yaml")
MSA_DIR = os.path.join(JOBDIR, "boltz2output/msas")
OUT_DIR = os.path.join(JOBDIR, "boltz2output/predictions")

COLABFOLD_BIN = config["colabfold_bin"]
MMSEQS_BIN = config["mmseqs_bin"]
DB1 = config["db1"]
DB_DIR = config["db_dir"]
THREADS = int(config["msa_threads"])
USE_ENV = str(config["use_env"])
TMP_ROOT = config["tmp_root"]
DIFFUSION_SAMPLES = int(config["diffusion_samples"])

MSA_OUTPUTS = [os.path.join(JOBDIR, path) for path in config["msa_outputs"]]
MSA_MOVES = "\n        ".join(
    f'test -s "$RUN_DIR/{i}.a3m"\n        mv -f "$RUN_DIR/{i}.a3m" "{path}"'
    for i, path in enumerate(MSA_OUTPUTS)
)

rule all:
    input:
        MSA_OUTPUTS,
        os.path.join(OUT_DIR, "DONE.txt")

rule msa:
    output:
        MSA_OUTPUTS
    shell:
        r'''
        set -euo pipefail
        mkdir -p "{MSA_DIR}" "{TMP_ROOT}"

        RUN_DIR="$(mktemp -d "{TMP_ROOT}/run_XXXXXX")"
        trap 'rm -rf "$RUN_DIR"' EXIT

        "{COLABFOLD_BIN}" \
          --threads {THREADS} \
          --mmseqs "{MMSEQS_BIN}" \
          --use-env {USE_ENV} \
          "{FAA}" \
          --db1 "{DB1}" \
          "{DB_DIR}" \
          "$RUN_DIR"

        {MSA_MOVES}
        '''

rule predict:
    input:
        MSA_OUTPUTS,
        YAML
    output:
        os.path.join(OUT_DIR, "DONE.txt")
    shell:
        r'''
        set -euo pipefail
        mkdir -p "{OUT_DIR}"

        "{config[boltz_bin]}" predict "{YAML}" \
          --out_dir "{OUT_DIR}" \
          --diffusion_samples {DIFFUSION_SAMPLES}

        copied_any=0
        while IFS= read -r -d '' cif; do
          bn="$(basename "$cif")"
          cp -f "$cif" "{JOBDIR}/$bn"
          copied_any=1
        done < <(find "{OUT_DIR}" -type f \( -name '*.cif' -o -name '*.mmcif' \) -print0)

        if [ "$copied_any" -eq 0 ]; then
          echo "ERROR: no CIF/mmCIF files found under {OUT_DIR}" >&2
          exit 1
        fi

        echo "ok" > "{output}"
        '''
""".lstrip()


def chain_ids_from_protein_id(raw_id):
    if isinstance(raw_id, list):
        return [str(chain_id) for chain_id in raw_id]
    return [str(raw_id if raw_id is not None else "A")]


def classify_system(proteins):
    unique_sequences = {sequence for _, sequence in proteins}
    if len(proteins) == 1:
        return "monomer"
    if len(unique_sequences) == 1:
        return "homomer"
    return "heteromer"


def extract_sequences_and_patch_yaml(yaml_file: Path, mode: str):
    data = yaml.safe_load(yaml_file.read_text())

    proteins = []

    for entry in data.get("sequences", []):
        if "protein" not in entry:
            continue

        protein = entry["protein"]
        chain_ids = chain_ids_from_protein_id(protein.get("id", "A"))
        sequence = protein.get("sequence")

        if not sequence:
            raise SystemExit(f"Protein chain(s) {', '.join(chain_ids)} have no sequence field.")

        for chain_id in chain_ids:
            proteins.append((chain_id, sequence))

    if not proteins:
        raise SystemExit("No protein entries found in YAML.")

    system_type = classify_system(proteins)

    if mode == "monomer" and system_type != "monomer":
        raise SystemExit(
            f"--mode monomer expects exactly 1 protein chain, found {len(proteins)}."
        )
    if mode == "complex" and len(proteins) == 1:
        raise SystemExit("--mode complex expects more than 1 protein chain.")

    representatives_by_sequence = OrderedDict()
    for chain_id, sequence in proteins:
        representatives_by_sequence.setdefault(sequence, chain_id)

    msa_by_sequence = {
        sequence: f"./boltz2output/msas/{representative}.a3m"
        for sequence, representative in representatives_by_sequence.items()
    }

    for entry in data.get("sequences", []):
        if "protein" not in entry:
            continue
        protein = entry["protein"]
        protein["msa"] = msa_by_sequence[protein["sequence"]]

    fasta_lines = []
    msa_outputs = []
    msa_plan = []
    for sequence, representative in representatives_by_sequence.items():
        fasta_lines.append(f">{representative}")
        fasta_lines.append(sequence)
        msa_outputs.append(f"boltz2output/msas/{representative}.a3m")

        chains = [chain_id for chain_id, chain_sequence in proteins if chain_sequence == sequence]
        msa_plan.append({
            "msa": f"boltz2output/msas/{representative}.a3m",
            "representative": representative,
            "chains": chains,
        })

    fasta_text = "\n".join(fasta_lines) + "\n"

    patched_yaml_text = yaml.safe_dump(
        data,
        sort_keys=False,
        default_flow_style=False,
    )

    metadata = {
        "system_type": system_type,
        "chain_count": len(proteins),
        "unique_sequence_count": len(representatives_by_sequence),
        "msa_outputs": msa_outputs,
        "msa_plan": msa_plan,
    }

    return fasta_text, patched_yaml_text, metadata


def make_sbatch(job_name, cores, mem, partition, gres, nodelist, conda_env):
    gres_line = f"#SBATCH --gres={gres}" if gres else ""
    nodelist_line = f"#SBATCH --nodelist={nodelist}" if nodelist else ""

    return f"""#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --cpus-per-task={cores}
#SBATCH --nodes=1
#SBATCH --mem={mem}
#SBATCH --partition={partition}
{gres_line}
{nodelist_line}

set -euo pipefail

JOBDIR="$SLURM_SUBMIT_DIR"
cd "$JOBDIR"

eval "$(/z/linux/bin/micromamba shell hook -s posix)"
micromamba activate "{conda_env}"

echo "Starting {job_name} on $(hostname) at $(date)"
snakemake --snakefile "$JOBDIR/Snakefile" --directory "$JOBDIR" --cores {cores}
echo "Finished {job_name} at $(date)"
"""


def main():
    p = argparse.ArgumentParser(
        description="Set up a Snakemake + Slurm Boltz2 job from an existing Boltz YAML."
    )

    p.add_argument("-y", "--yaml", required=True, help="Existing Boltz YAML file")
    p.add_argument("-o", "--outdir", required=True, help="Job directory to create")
    p.add_argument(
        "--mode",
        choices=["auto", "monomer", "complex"],
        default="auto",
        help="Validate as monomer/complex, or infer automatically. Default: auto.",
    )
    p.add_argument("--job-name", default=None)

    p.add_argument("--boltz-bin", default="boltz")
    p.add_argument("--diffusion-samples", type=int, default=1)

    p.add_argument("--tmp-root", default="/var/tmp/boltz_tmp")
    p.add_argument(
        "--colabfold-bin",
        default="/z/bio/biotools/miniconda/envs/colabfold_beta2/bin/colabfold_search",
    )
    p.add_argument("--mmseqs-bin", default="/z/bio/biotools/mmseqs")
    p.add_argument("--db1", default="/var/tmp/databases/colabfold/uniref30_2103")
    p.add_argument("--db-dir", default="/var/tmp/databases/colabfold/")
    p.add_argument("--msa-threads", type=int, default=64)
    p.add_argument("--use-env", type=int, default=0)

    p.add_argument("--sbatch-cores", type=int, default=16)
    p.add_argument("--sbatch-mem", default="16gb")
    p.add_argument("--sbatch-partition", default="low")
    p.add_argument("--sbatch-gres", default="shard:30")
    p.add_argument("--sbatch-nodelist", default="ai,pika")
    p.add_argument("--conda-env", default="/z/bio/biotools/miniconda/envs/boltz2/")

    args = p.parse_args()

    yaml_in = Path(args.yaml).resolve()
    jobdir = Path(args.outdir).resolve()
    jobdir.mkdir(parents=True, exist_ok=True)

    if not yaml_in.exists():
        raise SystemExit(f"YAML file not found: {yaml_in}")

    job_name = args.job_name or jobdir.name

    fasta_text, patched_yaml_text, metadata = extract_sequences_and_patch_yaml(
        yaml_in,
        mode=args.mode,
    )

    (jobdir / "input.yaml").write_text(patched_yaml_text)
    (jobdir / "msa_queries.fasta").write_text(fasta_text)
    (jobdir / "Snakefile").write_text(SNAKEFILE_TEXT)

    config = {
        "mode": args.mode,
        "system_type": metadata["system_type"],
        "chain_count": metadata["chain_count"],
        "unique_sequence_count": metadata["unique_sequence_count"],
        "msa_outputs": metadata["msa_outputs"],
        "job_name": job_name,
        "boltz_bin": args.boltz_bin,
        "diffusion_samples": args.diffusion_samples,
        "tmp_root": args.tmp_root,
        "colabfold_bin": args.colabfold_bin,
        "mmseqs_bin": args.mmseqs_bin,
        "db1": args.db1,
        "db_dir": args.db_dir,
        "msa_threads": args.msa_threads,
        "use_env": args.use_env,
    }
    (jobdir / "config.yaml").write_text(yaml.safe_dump(config, sort_keys=False))

    msa_plan = ["MSA plan", "========", ""]
    msa_plan.append(
        f"Detected {metadata['system_type']} with {metadata['chain_count']} chain(s) "
        f"and {metadata['unique_sequence_count']} unique protein sequence(s)."
    )
    msa_plan.append("")
    for item in metadata["msa_plan"]:
        msa_plan.append(f"{item['msa']}: chains {', '.join(item['chains'])}")
    msa_plan.append("")
    (jobdir / "MSA_PLAN.txt").write_text("\n".join(msa_plan))

    (jobdir / "run.sbatch").write_text(
        make_sbatch(
            job_name=job_name,
            cores=args.sbatch_cores,
            mem=args.sbatch_mem,
            partition=args.sbatch_partition,
            gres=args.sbatch_gres,
            nodelist=args.sbatch_nodelist,
            conda_env=args.conda_env,
        )
    )

    print(f"Created Boltz2 Snakemake job: {jobdir}")
    print(f"Patched YAML written to: {jobdir / 'input.yaml'}")
    print(f"MSA query FASTA written to: {jobdir / 'msa_queries.fasta'}")
    print(f"MSA plan written to: {jobdir / 'MSA_PLAN.txt'}")
    print(
        f"Detected {metadata['system_type']} with {metadata['chain_count']} chain(s) "
        f"and {metadata['unique_sequence_count']} unique protein sequence(s)."
    )
    print()
    print("Submit with:")
    print(f"  cd {jobdir}")
    print("  sbatch run.sbatch")


if __name__ == "__main__":
    main()
