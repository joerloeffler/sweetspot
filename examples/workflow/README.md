# Sweetspot workflow examples

Each numbered JSON file uses the shared repository-level `submission.json`:

```bash
python ../../sweetspot_workflow.py prepare \
  --settings 01_monomer.json \
  --outdir jobs/01_monomer
```

The examples cover:

0. `minimal_boltz2.json`: the smallest typical Boltz2 job; its name and local
   MSA configuration are inferred from the filename and shared profile.
1. `01_monomer.json`: ESMFold2-Fast monomer in single-sequence mode.
2. `02_protein_protein.json`: Boltz2 protein heterocomplex from one multirecord FASTA.
3. `03_protein_ligand.json`: Boltz2 protein with ligands read from a SMILES file.
4. `04_protein_protein_ligand.json`: AF3 two-protein complex with a CCD ligand.
5. `05_protein_ions.json`: ESMFold2 protein with explicit catalytic/structural ions.
6. `06_complex_man9.json`: two glycoproteins, with Man9 added at every detected sequon.
7. `07_trimer_man5.json`: homotrimer with Man5 and a single reused local MSA.
8. `08_ncaa_auto_x.json`: two `X` residues mapped, in order, to AIB and NLE.
9. `09_ncaa_explicit_complex.json`: four explicitly positioned NCAAs across two chains.

These sequences are deliberately short configuration examples, not biological
benchmarks. Inspect each generated `resolved_settings.json` and predictor input
before production use.

Cluster-specific MSA executables, databases, caches, and temporary directories
belong in the root `submission.json` under `workflow_defaults`. It is found
automatically. A normal settings file can therefore use `"align": true`
or `false`; use `"local"` or `"server"` only when a job needs to override the
profile/default source. Likewise, `"copies": 3` creates an automatically
named homotrimer, while `"chains": ["A", "B", "C"]` remains available when
explicit IDs matter.

Prepared jobs contain:

- `settings.json`: the concise user input, unchanged.
- `effective_settings.json`: user input merged with submission defaults.
- `resolved_settings.json`: fully expanded chains, glycans, bonds, and MSA plan.
- `input.yaml`: native, readable Boltz YAML (for Boltz2 jobs).
