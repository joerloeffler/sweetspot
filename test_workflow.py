import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import yaml

import sweetspot_workflow as workflow


class WorkflowTests(unittest.TestCase):
    def test_repository_submission_is_default_for_external_settings(self):
        with tempfile.TemporaryDirectory() as tmp:
            settings = Path(tmp) / "job.json"
            settings.write_text("{}")
            path = workflow.default_submission_path(settings)
            self.assertIsNotNone(path)
            self.assertEqual(path.name, "submission.json")
            self.assertEqual(path.parent, Path(workflow.__file__).resolve().parent)

    def test_homotrimer_reuses_msa_and_glycosylates_every_copy(self):
        settings = {
            "schema_version": 1,
            "name": "trimer",
            "backend": {"name": "boltz2"},
            "proteins": [{
                "id": "p", "sequence": "ANVTAA", "chains": ["A", "B", "C"],
                "glycosylation": {"enabled": True, "sites": "auto", "glycan": "NAG-NAG-MAN"},
            }],
            "msa": {"source": "local"},
        }
        resolved = workflow.resolve_settings(settings, Path.cwd())
        resolved["msa_plan"] = workflow.build_msa_plan(resolved)
        self.assertEqual(resolved["system_type"], "homomer")
        self.assertEqual(len(resolved["glycans"]), 3)
        self.assertEqual(len(resolved["msa_plan"]), 1)
        self.assertEqual(resolved["msa_plan"][0]["chains"], ["A", "B", "C"])
        with tempfile.TemporaryDirectory() as tmp:
            name = workflow.write_backend_input(Path(tmp), resolved)
            payload = yaml.safe_load((Path(tmp) / name).read_text())
        proteins = [item["protein"] for item in payload["sequences"] if "protein" in item]
        self.assertEqual(proteins[0]["id"], ["A", "B", "C"])

    def test_heteromer_can_glycosylate_both_proteins(self):
        settings = {
            "schema_version": 1, "name": "heteromer", "backend": "alphafold3",
            "proteins": [
                {"sequence": "ANVT", "chains": ["A"], "glycosylation": {"enabled": True}},
                {"sequence": "MNSTA", "chains": ["B"], "glycosylation": {"enabled": True}},
            ],
            "msa": {"source": "none"},
        }
        resolved = workflow.resolve_settings(settings, Path.cwd())
        self.assertEqual(resolved["system_type"], "heteromer")
        self.assertEqual({x["chain"] for x in resolved["glycosylation_sites"]}, {"A", "B"})

    def test_backend_aware_alignment_defaults_and_simple_switch(self):
        base = {"schema_version": 1, "proteins": [{"sequence": "AAAA"}]}
        boltz = workflow.resolve_settings({**base, "backend": "boltz2"}, Path.cwd())
        esm = workflow.resolve_settings({**base, "backend": "esmfold2"}, Path.cwd())
        local = workflow.resolve_settings({
            **base, "backend": "boltz2", "align": True,
            "msa": {"database": "/db/u30", "database_directory": "/db"},
        }, Path.cwd())
        disabled = workflow.resolve_settings({
            **base, "backend": "boltz2", "align": False,
        }, Path.cwd())
        self.assertEqual(boltz["msa"]["source"], "server")
        self.assertEqual(esm["msa"]["source"], "none")
        self.assertEqual(local["msa"]["source"], "local")
        self.assertEqual(disabled["msa"]["source"], "none")

    def test_af3_modifications_use_af3_field_names(self):
        settings = {
            "schema_version": 1, "backend": "alphafold3",
            "proteins": [{
                "sequence": "AAAA", "modifications": [{"position": 2, "ccd": "MSE"}],
            }],
            "align": False,
        }
        resolved = workflow.resolve_settings(settings, Path.cwd())
        resolved["msa_plan"] = workflow.build_msa_plan(resolved)
        entities = workflow.backend_entities(resolved, af3=True)
        self.assertEqual(entities[0]["protein"]["modifications"], [
            {"ptmType": "MSE", "ptmPosition": 2},
        ])

    def test_submission_backend_override(self):
        submission = {
            "slurm": {"job_name": "base", "gres": "gpu:1"},
            "environment": {"setup": ["base"]},
            "backends": {"esmfold2": {
                "slurm": {"job_name": "esm", "gres": "shard:60"},
                "environment": {"setup": ["esm-env"], "python": "/tools/run_python.sh"},
            }},
        }
        text = workflow.make_sbatch(submission, "esmfold2")
        self.assertIn("#SBATCH --job-name=esm", text)
        self.assertIn("#SBATCH --gres=shard:60", text)
        self.assertIn("esm-env", text)
        self.assertIn('/tools/run_python.sh "$JOBDIR/workflow.py"', text)
        self.assertNotIn("\nbase\n", text)

    def test_multirecord_fasta_becomes_complex(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "complex.fasta").write_text(">receptor\nAAAA\n>binder\nCCCC\n")
            resolved = workflow.resolve_settings({
                "schema_version": 1, "name": "complex", "backend": "boltz2",
                "proteins": [{"fasta": "complex.fasta", "chains": ["A", "B"]}],
                "msa": {"source": "none"},
            }, root)
            self.assertEqual(resolved["system_type"], "heteromer")
            self.assertEqual([p["chain"] for p in resolved["proteins"]], ["A", "B"])

    def test_prepare_writes_inspectable_job(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            settings = root / "settings.json"
            submission = root / "submission.json"
            settings.write_text(json.dumps({
                "schema_version": 1, "name": "one", "backend": "boltz2",
                "proteins": [{"sequence": "AAAA", "chains": ["A"]}],
                "msa": {"source": "none"},
            }))
            submission.write_text(json.dumps({"slurm": {"cpus_per_task": 2}}))
            workflow.prepare(settings, submission, root / "job")
            self.assertTrue((root / "job/resolved_settings.json").exists())
            self.assertTrue((root / "job/input.yaml").exists())
            self.assertTrue((root / "job/run.sbatch").exists())
            self.assertTrue((root / "job/effective_settings.json").exists())
            resolved = json.loads((root / "job/resolved_settings.json").read_text())
            self.assertEqual(resolved["prediction"]["max_parallel_samples"], 1)

    def test_submission_workflow_defaults_keep_settings_concise(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            settings = root / "settings.json"
            submission = root / "submission.json"
            settings.write_text(json.dumps({
                "schema_version": 1, "name": "one", "backend": "boltz2",
                "proteins": [{"sequence": "AAAA"}], "msa": "local",
            }))
            submission.write_text(json.dumps({
                "workflow_defaults": {"msa": {"database": "/db/u30", "database_directory": "/db"}}
            }))
            workflow.prepare(settings, submission, root / "job")
            effective = json.loads((root / "job/effective_settings.json").read_text())
            self.assertEqual(effective["msa"]["source"], "local")
            self.assertEqual(effective["msa"]["database"], "/db/u30")

    def test_boltz_command_uses_quality_and_fresh_run_flags(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            resolved = {
                "backend_input": "input.yaml",
                "backend_settings": {"executable": "boltz", "use_potentials": True},
                "prediction": {
                    "diffusion_samples": 20, "max_parallel_samples": 1,
                    "sampling_steps": 200, "recycling_steps": 10,
                },
                "msa": {"source": "server"},
                "output": {"directory": "results", "format": "mmcif"},
            }

            def produce_output(*args, **kwargs):
                output = root / "results/predictions/job"
                output.mkdir(parents=True)
                (output / "job_model_0.cif").write_text("data_job\n")

            with mock.patch("sweetspot_workflow.subprocess.run", side_effect=produce_output) as run:
                workflow.run_boltz(root, resolved)
            command = run.call_args.args[0]
            self.assertIn("--use_msa_server", command)
            self.assertIn("--use_potentials", command)
            self.assertIn("--override", command)
            self.assertEqual(command[command.index("--recycling_steps") + 1], "10")

    def test_paired_a3m_is_split_into_keyed_boltz_csvs(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "complex.a3m"
            source.write_text(">query\nABCDEF\n>hit\nABcC-EF\n")
            plan = [
                {"sequence": "ABC", "path": "A.csv"},
                {"sequence": "DEF", "path": "B.csv"},
            ]
            workflow.split_paired_a3m(source, plan, root)
            a = (root / "A.csv").read_text().splitlines()
            b = (root / "B.csv").read_text().splitlines()
            self.assertEqual(a, ["sequence,key", "ABC,0", "ABC,1"])
            self.assertEqual(b, ["sequence,key", "DEF,0", "-EF,1"])

    def test_paired_local_heteromer_uses_csv_paths(self):
        settings = {
            "schema_version": 1, "name": "pair", "backend": "boltz2",
            "proteins": [
                {"id": "one", "sequence": "AAAA", "chains": ["A"]},
                {"id": "two", "sequence": "CCCC", "chains": ["B"]},
            ],
            "msa": {"source": "local", "pairing": "paired"},
        }
        resolved = workflow.resolve_settings(settings, Path.cwd())
        resolved["msa_plan"] = workflow.build_msa_plan(resolved)
        self.assertEqual([item["path"] for item in resolved["msa_plan"]], ["msas/A.csv", "msas/B.csv"])
        self.assertEqual(resolved["warnings"], [])

    def test_unpaired_local_boltz_heteromer_uses_original_a3m_behavior(self):
        settings = {
            "backend": "boltz2",
            "proteins": [{"sequence": "AAAA"}, {"sequence": "CCCC"}],
            "msa": {"source": "local"},
        }
        resolved = workflow.resolve_settings(settings, Path.cwd())
        resolved["msa_plan"] = workflow.build_msa_plan(resolved)
        self.assertEqual([item["path"] for item in resolved["msa_plan"]], ["msas/A.a3m", "msas/B.a3m"])

    def test_esmfold2_input_keeps_named_glycan_bonds(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            settings = {
                "schema_version": 1, "name": "esm", "backend": {"name": "esmfold2", "model": "quality"},
                "proteins": [{"sequence": "ANVT", "chains": ["A"], "glycosylation": {"enabled": True}}],
                "msa": {"source": "none"},
            }
            resolved = workflow.resolve_settings(settings, root)
            resolved["msa_plan"] = workflow.build_msa_plan(resolved)
            filename = workflow.write_backend_input(root, resolved)
            payload = json.loads((root / filename).read_text())
            self.assertEqual(payload["sequences"][1]["ccd"], ["NAG", "NAG", "MAN"])
            self.assertEqual(payload["covalent_bonds_by_name"][0]["atom1"], ["A", 2, "ND2"])

    def test_automatic_ncaa_mapping_follows_x_order(self):
        sequence, mods = workflow.resolve_ncaas(
            {"mode": "auto", "ccds": ["AIB", "NLE"], "base_residue": "G"},
            "MAXKX", "A",
        )
        self.assertEqual(sequence, "MAGKG")
        self.assertEqual(mods, [{"position": 3, "ccd": "AIB"}, {"position": 5, "ccd": "NLE"}])

    def test_explicit_ncaas_are_filtered_per_chain(self):
        config = {"mode": "explicit", "items": [
            {"chain": "A", "position": 2, "ccd": "AIB"},
            {"chain": "B", "position": 3, "ccd": "ORN"},
            {"chain": "B", "position": 5, "ccd": "SAR"},
        ]}
        seq_a, mods_a = workflow.resolve_ncaas(config, "AAAAA", "A")
        seq_b, mods_b = workflow.resolve_ncaas(config, "AAAAA", "B")
        self.assertEqual(seq_a, "AAAAA")
        self.assertEqual([x["position"] for x in mods_a], [2])
        self.assertEqual([x["position"] for x in mods_b], [3, 5])

    def test_esm_cache_resolution_prefers_complete_snapshot(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            snapshot = root / "models--Biohub--ESMFold2/snapshots/abc"
            snapshot.mkdir(parents=True)
            (snapshot / "model.safetensors").write_text("")
            resolved = Path(workflow.resolve_cached_esm_model("Biohub/ESMFold2", [str(root)]))
            self.assertTrue(resolved.samefile(snapshot))


if __name__ == "__main__":
    unittest.main()
