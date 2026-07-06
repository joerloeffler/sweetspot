<p align="center">
  <img src="sweetspot.png" alt="Sweetspot logo" width="420">
</p>

# Sweetspot

Sweetspot prepares glycoprotein prediction inputs from an existing PDB, mmCIF,
or FASTA input. It scans protein chains for canonical N-glycosylation sequons
(`N-X-S/T`, where `X` is not proline), assigns a glycan chain to each detected
site, and writes backend-specific input files for structure prediction workflows.

Supported output targets:

### Disclaimer: currently only evaluated for boltz2!

- RoseTTAFold-All-Atom / RF3-style covalent input (`rfaa`)
- AlphaFold3 JSON (`af3`)
- Boltz2 YAML (`boltz2`)
- ESMFold / ESM2 protein-only FASTA (`esm`)

The default glycan is `NAG-NAG-MAN`. High-mannose presets are available with
`-g man5` and `-g man9`.

## Repository Layout

```text
sweetspot.py                 Main CLI for detecting sites and writing inputs
sanitize_linked_glycans.py   Standalone helper for predictor-emitted PDB files
cluster_prep.py              Boltz2 Snakemake + Slurm job scaffold helper
example/5m8n.pdb             Example input structure
example/sweet_prep/          Example generated outputs
sweetspot.png                README logo
```

## Install

Sweetspot is a small Python CLI. A virtual environment is recommended.

```bash
python -m venv .venv
source .venv/bin/activate
pip install biopython pyyaml
```

`PyYAML` is required for `rfaa` and `boltz2` output. `biopython` is required for
parsing PDB/mmCIF structures.

## Quick Start

Generate the default RFAA input from the included example:

```bash
python sweetspot.py -f example/5m8n.pdb
```

Generate every supported backend input:

```bash
python sweetspot.py -f example/5m8n.pdb -m all
```

Write outputs to a custom directory:

```bash
python sweetspot.py -f example/5m8n.pdb -m boltz2 -o example/sweet_prep
```

Run from a single-chain FASTA instead of a structure:

```bash
python sweetspot.py -f protein.fasta -m boltz2
```

Run a multi-record FASTA as a complex:

```bash
python sweetspot.py -f complex.fasta -m boltz2 --complex
```

Run a homomultimer from a single-record FASTA, for example a homodimer:

```bash
python sweetspot.py -f protein.fasta -m boltz2 --multimer 2
```

Use a custom glycan:

```bash
python sweetspot.py -f protein.pdb -m all -g NAG-NAG-MAN-MAN
```

Use a high-mannose preset:

```bash
python sweetspot.py -f protein.pdb -m boltz2 -g man5
python sweetspot.py -f protein.pdb -m boltz2 -g man9
```

Use a branched glycan layout:

```bash
python sweetspot.py -f protein.pdb -m boltz2 -t branched -g NAG-NAG-MAN-MAN
```

List supported glycan residue names:

```bash
python sweetspot.py --list-glycans
```

## Output

By default, Sweetspot writes to:

```text
<input_stem>_glyco_inputs/
```

Common output:

- `sites/nglyc_sites.tsv`: detected sequons and assigned glycan chain IDs

Backend-specific output:

- `af3/<name>_all_nglyc_sites_af3.json`: AlphaFold3 input JSON
- `boltz2/<name>_all_nglyc_sites_boltz2.yaml`: Boltz2 input YAML
- `rfaa/rfaa.yaml`: RFAA configuration
- `rfaa/run_rfaa.sh`: convenience launcher for an RFAA environment
- `esm/protein.fasta`: protein-only FASTA for ESM workflows

For AF3 and Boltz2, Sweetspot also writes:

- `linked_glycans.json`: intended protein-glycan and glycan-glycan bonds
- `sanitize_linked_glycans.py`: helper to clean predictor-emitted PDB files
- `README_glycan_handoff.txt`: notes for converting predictor output toward
  GLYCAM/OpenMM-ready linked glycans

## Glycan Handoff

AF3 and Boltz2 inputs represent glycans as CCD monosaccharides plus covalent
bond constraints. Some predictor outputs may retain standalone CCD leaving-group
atoms, such as anomeric `O1`, even when the sugar is covalently linked.

Treat sanitization as the postprocessing step after folding. Once AF3 or Boltz2
has produced a PDB, run the generated sanitizer from the same backend output
directory:

```bash
python sanitize_linked_glycans.py <predictor_output.pdb> \
  -m linked_glycans.json \
  -o <cleaned.pdb>
```

The sanitizer removes linked-sugar leaving-group atoms described in the manifest
and rewrites `CONECT` records from Sweetspot's intended covalent bonds.

A typical Boltz2 flow looks like:

```bash
python sweetspot.py -f protein.pdb -m boltz2
cd protein_glyco_inputs/boltz2
boltz predict protein_all_nglyc_sites_boltz2.yaml --out_dir boltz2output
python sanitize_linked_glycans.py boltz2output/<folded_model>.pdb \
  -m linked_glycans.json \
  -o protein_all_nglyc_sites.cleaned.pdb
```

For AF3, run the same sanitizer after downloading or exporting the folded PDB
from the AF3 job output.

## RFAA Notes

RFAA covalent glycan input uses one merged glycan SDF per site rather than the
CCD glycan residue list directly. If a file named `glycan.sdf` exists when you
run Sweetspot, it is copied into the generated `rfaa/` directory. Otherwise,
`rfaa/README.txt` explains where to place the SDF before running RFAA.

Useful RFAA options:

```bash
python sweetspot.py -f protein.pdb -m rfaa \
  --rfaa-glycan-sdf glycan.sdf \
  --rfaa-link-atom-index 1 \
  --rfaa-link-chirality CW
```

## Boltz2 Cluster Helper

`cluster_prep.py` can turn an existing Boltz2 YAML into a Snakemake + Slurm job
directory. It detects monomers, homomers, and heteromers, writes an MSA query
FASTA for each unique protein sequence, reuses the same MSA for identical
chains, patches MSA paths into the YAML, writes a `Snakefile`, and creates
`run.sbatch`.

```bash
python cluster_prep.py \
  --yaml example/sweet_prep/boltz2/5m8n_all_nglyc_sites_boltz2.yaml \
  --outdir boltz_job
```

Submit from the generated job directory:

```bash
cd boltz_job
sbatch run.sbatch
```

## License

No license file is currently included in this repository.
