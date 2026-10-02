<p align="center">
  <img src="sweetspot.png" alt="Sweetspot logo" width="420">
</p>

# Sweetspot

Sweetspot prepares and runs glycoprotein and molecular-complex structure
prediction jobs. It supports proteins, multimers, glycans, ligands, structural
ions, and modified or noncanonical amino acids.

The recommended interface is the JSON workflow in `sweetspot_workflow.py`. It
creates inspectable, portable job directories for Boltz2, ESMFold2, and
AlphaFold 3. The original `sweetspot.py` CLI remains available for generating
standalone predictor inputs from PDB, mmCIF, or FASTA.

## Install Sweetspot

Preparation requires Python 3.10 or newer, Biopython, and PyYAML:

```bash
python -m venv .venv
source .venv/bin/activate
pip install biopython pyyaml
```

Prediction environments and model data are backend-specific and are not
installed by Sweetspot. A reusable cluster profile can activate those
environments in the generated Slurm job.

## Install prediction backends with Mamba or Conda

Use separate environments for Boltz2 and ESMFold2. The backend projects
currently publish their Python packages through PyPI, so Mamba/Conda manages
the isolated environment and pip installs the official package.

### Boltz2

The [official Boltz repository](https://github.com/jwohlwend/boltz) recommends a
fresh environment and `boltz[cuda]` for NVIDIA GPU systems:

```bash
micromamba create -n sweetspot-boltz2 -c conda-forge python=3.11 pip -y
micromamba activate sweetspot-boltz2
python -m pip install --upgrade pip
python -m pip install --upgrade "boltz[cuda]"

boltz --help
python -c "import boltz; print('Boltz import OK')"
```

With Conda, replace `micromamba` with `conda`. For a CPU-only installation,
install `boltz` without the `[cuda]` extra, although prediction will be much
slower. The model and CCD cache are downloaded on first use unless `BOLTZ_CACHE`
points to a populated shared cache.

### ESMFold2

The [official Biohub ESM repository](https://github.com/Biohub/esm) publishes
ESMFold2 through the `esm` package:

```bash
micromamba create -n sweetspot-esmfold2 -c conda-forge python=3.12 pip -y
micromamba activate sweetspot-esmfold2
python -m pip install --upgrade pip
python -m pip install --upgrade esm

python -c "from esm.models.esmfold2 import ESMFold2InputBuilder; print('ESMFold2 import OK')"
```

ESMFold2 downloads its model weights from Hugging Face on first use. On a
cluster, set `HF_HOME` to a shared cache, populate it on a networked node, and
set `HF_HUB_OFFLINE=1` only after the complete ESMFold2 and ESMC snapshots are
present. GPU/PyTorch compatibility is cluster-specific; if your site provides
a tested CUDA-enabled PyTorch module or wheel, install that before `esm`.

Finally, edit the repository-level `submission.json` so its environment
activation commands, Python launcher, caches, Slurm resources, and local MSA
paths match your cluster. The committed profile is the in-house deployment
example, not a portable cluster autodetection mechanism.

## Recommended JSON workflow

Scientific choices belong in a job JSON. Machine-specific paths, Slurm
resources, environment activation, and local-MSA database settings belong in
the repository-level `submission.json`.

A concise Boltz2 job can look like this:

```json
{
  "name": "target_binder_man9",
  "backend": "boltz2",
  "proteins": [
    {
      "id": "target",
      "fasta": "inputs/target.fasta",
      "chains": ["A"],
      "glycosylation": {
        "enabled": true,
        "sites": "auto",
        "glycan": "man9",
        "tree": "branched"
      }
    },
    {
      "id": "binder",
      "fasta": "inputs/binder.fasta",
      "chains": ["B"]
    }
  ],
  "align": "local",
  "prediction": {"diffusion_samples": 20}
}
```

Prepare the portable job directory:

```bash
python sweetspot_workflow.py prepare \
  --settings path/to/job.json \
  --outdir path/to/prepared_job
```

The repository-level `submission.json` is loaded automatically. A
`submission.json` beside a job takes precedence, and
`--submission path/to/submission.json` explicitly selects another profile.

Inspect `resolved_settings.json` and the generated predictor input before
running. Then either run directly:

```bash
python path/to/prepared_job/workflow.py run \
  --job-dir path/to/prepared_job
```

or submit through Slurm:

```bash
cd path/to/prepared_job
sbatch run.sbatch
```

### Prepared job files

| File | Purpose |
| --- | --- |
| `settings.json` | Original concise job definition |
| `submission.json` | Cluster and environment profile |
| `effective_settings.json` | Job settings merged with profile defaults |
| `resolved_settings.json` | Expanded chains, glycans, bonds, MSAs, and defaults |
| `input.yaml` | Boltz2 input |
| `alphafold3_input.json` | AlphaFold 3 input |
| `esmfold2_input.json` | ESMFold2 input |
| `msa_queries.fasta` | Local-MSA queries; omitted for server/no-MSA jobs |
| `workflow.py` | Portable runner copied into the job |
| `run.sbatch` | Generated Slurm wrapper |

Only the input matching the selected backend is generated.

## Alignment policy

The job-facing control is `align`:

| Value | Behavior |
| --- | --- |
| `false` or `"none"` | Sequence-only inference |
| `true` | Use configured local settings; Boltz2 falls back to its server when none are configured |
| `"local"` | Run the configured local ColabFold/MMseqs search |
| `"server"` | Let Boltz2 use its MSA server |

Boltz2 defaults to server MSA when no alignment setting is present. ESMFold2
defaults to sequence-only inference. The older `"msa": "local|server|none"`
syntax remains supported.

Local Boltz2 heteromers use one A3M per unique protein sequence, matching the
established cluster pipeline. Identical copies reuse the same alignment. These
independent MSAs are unpaired; use Boltz2's MSA server when paired interface
information is important. Explicit paired local mode requires an MMseqs
database with the necessary taxonomy data.

## Glycosylation

Automatic N-glycosylation detects canonical `N-X-S/T` sequons where `X` is not
proline:

```json
"glycosylation": {
  "enabled": true,
  "sites": "auto",
  "glycan": "man9",
  "tree": "branched"
}
```

The default glycan is `NAG-NAG-MAN`. `man5` and `man9` are built-in branched
high-mannose presets. Sites can also be selected explicitly or excluded; see
the JSON examples under `examples/workflow/`.

## Samples and memory

`prediction.diffusion_samples` is the total number of structures requested.
`prediction.max_parallel_samples` controls how many samples are generated at
once and defaults to `1` to limit peak GPU memory use:

```json
"prediction": {
  "diffusion_samples": 20,
  "max_parallel_samples": 1,
  "seed": 0
}
```

For ESMFold2, keep the full model in `float32` unless the installed release has
been validated for whole-model BF16.

## Noncanonical amino acids

Use `X` as the FASTA placeholder and map each occurrence to a wwPDB Chemical
Component Dictionary code. The canonical `base_residue` is used for the model
sequence and MSA; the CCD modification supplies the intended chemistry.

```json
"noncanonical_amino_acids": {
  "mode": "auto",
  "ccds": ["AIB"],
  "base_residue": "A"
}
```

For explicit positions:

```json
"noncanonical_amino_acids": {
  "mode": "explicit",
  "base_residue": "A",
  "items": [
    {"chain": "B", "position": 16, "ccd": "AIB"}
  ]
}
```

Positions are one-based. CCD availability depends on the exact backend release
and local CCD cache; recognition of a code does not guarantee that its geometry
has been benchmarked by a model.

## Clean examples and deployment test

`examples/workflow/` contains schema-valid examples for monomers, protein
complexes, ligands, ions, Man5/Man9 glycans, and NCAAs. Prepare the minimal
Boltz2 example with:

```bash
python sweetspot_workflow.py prepare \
  --settings examples/workflow/minimal_boltz2.json \
  --outdir minimal_job
```

`dev/workflow_prep_test/` is the untracked EK0080 plus binder deployment matrix:

- standard binder and AIB16 binder;
- Boltz2 with local MSA;
- ESMFold2 with local MSA and without MSA;
- branched Man9 at all detected target sites;
- 20 diffusion samples per job.

Prepare all six jobs:

```bash
bash dev/workflow_prep_test/prepare_all.sh
```

See `dev/workflow_prep_test/README.md` for exact files and submission commands.

Only source structures, FASTA/SMILES inputs, and JSON configurations are kept
in the tracked example directories. Prepared jobs, predictions, MSAs, Slurm
logs, molecular-dynamics files, and visualization sessions are intentionally
not committed.

## Legacy input generator

Use `sweetspot.py` when only backend input generation is needed:

```bash
python sweetspot.py -f protein.pdb -m boltz2 -g man9 -t branched
python sweetspot.py -f complex.fasta -m boltz2 --complex
python sweetspot.py -f protein.fasta -m boltz2 --multimer 3
python sweetspot.py --list-glycans
```

Use `-m all` to generate every legacy target. By default, output is written to
`<input_stem>_glyco_inputs/`. AlphaFold 3 and Boltz2 outputs also include a
`linked_glycans.json` manifest and `sanitize_linked_glycans.py` helper for
post-prediction glycan cleanup.

`cluster_prep.py` is the older Snakemake/Slurm wrapper for an existing Boltz2
YAML. New workflows should normally use `sweetspot_workflow.py`.

## Repository layout

```text
sweetspot_workflow.py       Recommended JSON preparation and execution workflow
sweetspot.py                Legacy input generator
cluster_prep.py             Legacy Boltz2 Snakemake/Slurm helper
sanitize_linked_glycans.py  Predictor-output glycan cleanup helper
schemas/                    JSON schemas for jobs and submission profiles
submission.json             Default cluster, environment, and local-MSA profile
examples/workflow/          Reusable schema-valid job examples
dev/workflow_prep_test/     Untracked six-job EK0080 deployment test matrix
example/                    Raw legacy PDB example inputs only
test_workflow.py            Workflow unit tests
```

## Verification

Run the dependency-light test suite with:

```bash
python -m unittest -v
```

These tests validate configuration resolution and generated inputs but do not
replace one real smoke prediction for each installed backend/environment.

## Current limitations

- Local paired heteromer MSAs require taxonomy-enabled MMseqs databases.
- Structural ions are explicit entities, not a bulk solvent concentration.
- ESMFold2 and AlphaFold 3 support depends on the installed backend APIs and
  model assets.
- Glycan atom names and custom CCDs should be checked against the CCD bundle
  used by the prediction environment.

## License

No license file is currently included in this repository.
