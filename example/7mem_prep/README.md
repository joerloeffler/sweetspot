# 7MEM Multichain Example

This example uses `example/7MEM.pdb`, a 2-chain system, to exercise Sweetspot's
multichain glycoprotein setup.

## Generate Folding Inputs

From the repository root:

```bash
python sweetspot.py -f example/7MEM.pdb -m all -o example/7mem_prep
```

Sweetspot detected 7 N-glycosylation sequons across chains A and C and assigned
glycan chain IDs B and D through I.

Generated setup folders:

- `sites/`: detected glycosylation sites
- `af3/`: AlphaFold3 JSON plus linked-glycan sanitizer files
- `boltz2/`: Boltz2 YAML plus linked-glycan sanitizer files
- `rfaa/`: RFAA/RF3-style YAML, per-chain FASTA files, and run helper
- `esm/`: protein-only FASTA

## Multichain MSA Setup

The Boltz2 cluster scaffold was generated with:

```bash
python cluster_prep.py \
  --yaml example/7mem_prep/boltz2/7MEM_all_nglyc_sites_boltz2.yaml \
  --outdir example/7mem_prep/boltz2/cluster_job
```

`cluster_prep.py` classified the system as a heteromer with 2 protein chains
and 2 unique protein sequences:

```text
boltz2output/msas/A.a3m: chains A
boltz2output/msas/C.a3m: chains C
```

The generated `boltz2/cluster_job/input.yaml` patches each protein entry to the
appropriate shared MSA path, and `boltz2/cluster_job/msa_queries.fasta` contains
one query per unique sequence.

## Run Boltz2 On The Cluster

From the generated cluster job directory:

```bash
cd example/7mem_prep/boltz2/cluster_job
sbatch run.sbatch
```

The Slurm job runs MSA search, then `boltz predict` against `input.yaml`.

## Postprocess Folded PDBs

After folding, sanitize the predictor-emitted PDB before GLYCAM/OpenMM handoff.
Run this from the `boltz2/` output directory, using the generated manifest:

```bash
cd example/7mem_prep/boltz2
python sanitize_linked_glycans.py <boltz_output.pdb> \
  -m linked_glycans.json \
  -o 7MEM_all_nglyc_sites.cleaned.pdb
```

For AF3, use the same pattern from `example/7mem_prep/af3/` after exporting the
folded PDB.

## Notes

- This repository contains generated setup files, not completed AF3/Boltz2/RFAA
  folding results for 7MEM.
- RFAA uses one merged glycan SDF per glycan chain; place the SDF expected by
  `rfaa/README.txt` before running RFAA.
- ESM output is protein-only and does not explicitly model glycans.
