#!/usr/bin/env python3
"""JSON-driven molecular-complex preparation and execution for Sweetspot.

The preparation stage is dependency-light and produces resolved, inspectable
inputs.  The run stage may be executed directly or through generated Slurm.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
from collections import OrderedDict
from pathlib import Path
from typing import Any

try:
    import yaml
except ImportError:
    yaml = None  # type: ignore

try:
    import sweetspot
except ImportError:
    # A prepared job carries resolved inputs and can run without the source tree.
    sweetspot = None  # type: ignore


BACKENDS = {"boltz2", "alphafold3", "esmfold2"}
ION_CODES = {"MG", "ZN", "CL", "CA", "NA", "MN", "K", "FE", "CU", "CO"}
CHAIN_IDS = list("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789")


class ConfigError(ValueError):
    pass


def read_json(path: Path) -> dict[str, Any]:
    try:
        data = json.loads(path.read_text())
    except FileNotFoundError as exc:
        raise ConfigError(f"File not found: {path}") from exc
    except json.JSONDecodeError as exc:
        raise ConfigError(f"Invalid JSON in {path}: {exc}") from exc
    if not isinstance(data, dict):
        raise ConfigError(f"Top level of {path} must be a JSON object")
    return data


def write_json(path: Path, data: Any) -> None:
    path.write_text(json.dumps(data, indent=2) + "\n")


def merge_settings(defaults: dict[str, Any], settings: dict[str, Any]) -> dict[str, Any]:
    """Recursively merge portable job defaults with user-facing settings."""
    merged = dict(defaults)
    for key, value in settings.items():
        if isinstance(value, dict) and isinstance(merged.get(key), dict):
            merged[key] = merge_settings(merged[key], value)
        else:
            merged[key] = value
    return merged


def resolve_path(base: Path, value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else base / path


def require_keys(obj: dict[str, Any], keys: set[str], where: str) -> None:
    missing = sorted(keys - obj.keys())
    if missing:
        raise ConfigError(f"{where} is missing required field(s): {', '.join(missing)}")


def load_sequences(protein: dict[str, Any], base: Path) -> list[tuple[str, str]]:
    has_sequence = bool(protein.get("sequence"))
    has_fasta = bool(protein.get("fasta"))
    if has_sequence == has_fasta:
        raise ConfigError("Each protein needs exactly one of 'sequence' or 'fasta'")
    if has_sequence:
        return [(str(protein.get("id", "protein")), re.sub(r"\s+", "", protein["sequence"]).upper())]
    path = resolve_path(base, protein["fasta"])
    try:
        return sweetspot.parse_fasta(path)
    except FileNotFoundError as exc:
        raise ConfigError(f"Protein FASTA not found: {path}") from exc


def allocate_chain_ids(requested: Any, count: Any, used: set[str]) -> list[str]:
    if isinstance(count, bool):
        raise ConfigError("Protein 'copies' must be a positive integer")
    if requested is not None:
        ids = [requested] if isinstance(requested, str) else list(requested)
        if count is not None and len(ids) != int(count):
            raise ConfigError("'chains' and 'copies' describe different numbers of protein chains")
    else:
        try:
            copies = int(count or 1)
        except (TypeError, ValueError) as exc:
            raise ConfigError("Protein 'copies' must be a positive integer") from exc
        if copies < 1 or isinstance(count, bool):
            raise ConfigError("Protein 'copies' must be a positive integer")
        ids = []
        for _ in range(copies):
            cid = next((x for x in CHAIN_IDS if x not in used and x not in ids), "")
            ids.append(cid)
    if not ids or any(not isinstance(x, str) or not x for x in ids):
        raise ConfigError("Could not allocate protein chain IDs")
    if len(set(ids)) != len(ids) or any(x in used for x in ids):
        raise ConfigError(f"Duplicate chain ID in {ids}")
    used.update(ids)
    return ids


def normalize_modifications(raw: Any, sequence: str, chain: str) -> list[dict[str, Any]]:
    result = []
    for mod in raw or []:
        require_keys(mod, {"position", "ccd"}, f"modification on chain {chain}")
        pos = int(mod["position"])
        if not 1 <= pos <= len(sequence):
            raise ConfigError(f"Modification {chain}:{pos} is outside sequence length {len(sequence)}")
        result.append({"position": pos, "ccd": str(mod["ccd"]).upper()})
    return result


def resolve_ncaas(
    config: dict[str, Any] | None, sequence: str, chain: str
) -> tuple[str, list[dict[str, Any]]]:
    """Replace NCAA positions with a canonical MSA placeholder and annotate CCDs."""
    if not config:
        return sequence, []
    mode = config.get("mode", "explicit")
    base = str(config.get("base_residue", "A")).upper()
    if len(base) != 1 or base not in set("ACDEFGHIKLMNPQRSTVWY"):
        raise ConfigError("noncanonical_amino_acids.base_residue must be one canonical amino-acid letter")
    modifications: list[dict[str, Any]] = []
    if mode == "auto":
        ccds = [str(x).upper() for x in config.get("ccds", [])]
        positions = [i + 1 for i, residue in enumerate(sequence) if residue == "X"]
        if len(ccds) != len(positions):
            raise ConfigError(
                f"Automatic NCAA mapping on chain {chain} found {len(positions)} X residue(s) "
                f"but received {len(ccds)} CCD code(s)"
            )
        modifications = [{"position": pos, "ccd": ccd} for pos, ccd in zip(positions, ccds)]
    elif mode == "explicit":
        for item in config.get("items", []):
            if item.get("chain") not in (None, chain):
                continue
            require_keys(item, {"position", "ccd"}, f"explicit NCAA on chain {chain}")
            modifications.append({"position": int(item["position"]), "ccd": str(item["ccd"]).upper()})
    else:
        raise ConfigError("noncanonical_amino_acids.mode must be 'auto' or 'explicit'")
    chars = list(sequence)
    seen = set()
    for mod in modifications:
        pos = mod["position"]
        if not 1 <= pos <= len(chars):
            raise ConfigError(f"NCAA {chain}:{pos} is outside sequence length {len(chars)}")
        if pos in seen:
            raise ConfigError(f"More than one NCAA was assigned to {chain}:{pos}")
        seen.add(pos)
        chars[pos - 1] = base
    if mode == "explicit":
        unresolved = [i + 1 for i, residue in enumerate(chars) if residue == "X"]
        if unresolved:
            raise ConfigError(f"Explicit NCAA mapping left unresolved X residue(s) on chain {chain}: {unresolved}")
    return "".join(chars), modifications


def choose_glycosites(chain: str, sequence: str, config: dict[str, Any]) -> list[dict[str, Any]]:
    if not config.get("enabled", False):
        return []
    residues = [(str(i), aa) for i, aa in enumerate(sequence, 1)]
    detected = sweetspot.find_nglyc_sites({chain: residues})
    requested = config.get("sites", "auto")
    excluded = {int(x) for x in config.get("exclude_positions", [])}
    if requested == "auto":
        return [site for site in detected if int(site["asn_chain_index"]) not in excluded]
    if not isinstance(requested, list):
        raise ConfigError("glycosylation.sites must be 'auto' or a list")
    by_pos = {int(site["asn_chain_index"]): site for site in detected}
    selected = []
    for entry in requested:
        entry = {"position": entry} if isinstance(entry, int) else entry
        pos = int(entry["position"])
        if not 1 <= pos <= len(sequence):
            raise ConfigError(f"Glycosylation site {chain}:{pos} is outside the sequence")
        site = dict(by_pos.get(pos, {
            "site_id": "explicit", "chain": chain, "asn_resid": str(pos),
            "asn_chain_index": pos, "sequon": sequence[pos - 1:pos + 2],
            "x_resid": str(pos + 1), "st_resid": str(pos + 2),
        }))
        site["explicit"] = True
        if isinstance(entry, dict) and entry.get("glycan"):
            site["glycan_override"] = entry["glycan"]
        selected.append(site)
    return selected


def parse_smiles_file(path: Path) -> list[dict[str, Any]]:
    ligands = []
    try:
        lines = path.read_text().splitlines()
    except FileNotFoundError as exc:
        raise ConfigError(f"Ligand SMILES file not found: {path}") from exc
    for line_no, raw in enumerate(lines, 1):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        fields = line.split(maxsplit=1)
        ligands.append({"smiles": fields[0], "name": fields[1] if len(fields) > 1 else f"ligand_{line_no}"})
    return ligands


def allocate_nonprotein_id(used: set[str], prefix: str = "L") -> str:
    for cid in CHAIN_IDS:
        if cid not in used:
            used.add(cid)
            return cid
    index = 1
    while f"{prefix}{index}" in used:
        index += 1
    cid = f"{prefix}{index}"
    used.add(cid)
    return cid


def normalize_msa_settings(settings: dict[str, Any], backend: str) -> dict[str, Any]:
    """Resolve the concise ``align`` switch and the advanced MSA settings.

    ``align`` is the job-facing control. The ``msa`` object can remain in a
    reusable profile and hold local paths, commands, or a preferred source.
    Existing ``msa`` strings/objects remain supported.
    """
    raw = settings.get("msa", {})
    if isinstance(raw, str):
        msa: dict[str, Any] = {"source": raw}
    elif isinstance(raw, dict):
        msa = dict(raw)
    else:
        raise ConfigError("msa must be a source string or an object")

    align = settings.get("align")
    if align is not None:
        if isinstance(align, bool):
            if not align:
                msa["source"] = "none"
            else:
                # A configured source wins. With no profile, Boltz's server is
                # the zero-setup choice; other adapters use local alignment.
                has_local_profile = any(
                    key in msa for key in ("database", "database_directory", "command", "mmseqs")
                )
                msa.setdefault(
                    "source", "local" if has_local_profile else "server" if backend == "boltz2" else "local"
                )
        elif isinstance(align, str) and align in {"local", "server", "none"}:
            msa["source"] = align
        else:
            raise ConfigError("align must be true, false, 'local', 'server', or 'none'")
    elif "source" not in msa:
        # Safe, useful defaults with no per-job MSA input at all.
        msa["source"] = "server" if backend == "boltz2" else "none" if backend == "esmfold2" else "local"

    msa.setdefault("reuse_identical_sequences", True)
    if msa["source"] not in {"local", "server", "none"}:
        raise ConfigError("MSA source must be local, server, or none")
    if msa["source"] == "server" and backend != "boltz2":
        raise ConfigError("Server MSA is currently supported only by the Boltz2 adapter")
    return msa


def resolve_settings(settings: dict[str, Any], base: Path) -> dict[str, Any]:
    if sweetspot is None:
        raise ConfigError("Preparation requires sweetspot.py on PYTHONPATH")
    backend_cfg = settings.get("backend", {})
    backend = backend_cfg if isinstance(backend_cfg, str) else backend_cfg.get("name")
    if backend not in BACKENDS:
        raise ConfigError(f"backend.name must be one of: {', '.join(sorted(BACKENDS))}")
    proteins_in = settings.get("proteins")
    if not isinstance(proteins_in, list) or not proteins_in:
        raise ConfigError("settings.json needs a non-empty 'proteins' list")

    used: set[str] = set()
    proteins = []
    glycosites = []
    for number, raw in enumerate(proteins_in, 1):
        records = load_sequences(raw, base)
        if len(records) > 1 and raw.get("copies"):
            raise ConfigError("'copies' cannot be combined with a multirecord FASTA")
        requested_chains = raw.get("chains")
        if len(records) > 1 and requested_chains is not None and len(requested_chains) != len(records):
            raise ConfigError("For a multirecord FASTA, 'chains' must match the number of records")
        record_groups = []
        if len(records) == 1:
            record_groups.append((records[0], allocate_chain_ids(requested_chains, raw.get("copies"), used)))
        else:
            allocated = allocate_chain_ids(requested_chains, len(records), used)
            record_groups.extend((record, [chain]) for record, chain in zip(records, allocated))
        for record_index, ((header, sequence), chains) in enumerate(record_groups, 1):
            invalid = sorted(set(sequence) - sweetspot.VALID_FASTA_AA)
            if not sequence or invalid:
                raise ConfigError(f"Invalid protein sequence for {raw.get('id', header)}: {invalid}")
            entity_id = str(raw.get("id") or header or f"protein_{number}_{record_index}")
            if len(records) > 1 and raw.get("id"):
                entity_id = f"{entity_id}:{header}"
            glyco_cfg = raw.get("glycosylation", {"enabled": False})
            skip_chains = set(glyco_cfg.get("exclude_chains", []))
            for chain in chains:
                chain_sequence, ncaa_mods = resolve_ncaas(raw.get("noncanonical_amino_acids"), sequence, chain)
                legacy_mods = normalize_modifications(raw.get("modifications"), chain_sequence, chain)
                occupied = {x["position"] for x in ncaa_mods}
                duplicate = occupied.intersection(x["position"] for x in legacy_mods)
                if duplicate:
                    raise ConfigError(f"Duplicate modification/NCAA positions on chain {chain}: {sorted(duplicate)}")
                proteins.append({
                    "entity_id": entity_id, "chain": chain, "sequence": chain_sequence,
                    "modifications": ncaa_mods + legacy_mods,
                })
                if chain not in skip_chains:
                    for site in choose_glycosites(chain, chain_sequence, glyco_cfg):
                        glycan_name = site.pop("glycan_override", glyco_cfg.get("glycan", sweetspot.DEFAULT_GLYCAN))
                        residues = sweetspot.parse_glycan(glycan_name)
                        sweetspot.validate_glycan(residues)
                        site.update({
                            "glycan": glycan_name, "residues": residues,
                            "tree": glyco_cfg.get("tree", "linear"),
                            "preset": sweetspot.glycan_preset_name(glycan_name),
                        })
                        glycosites.append(site)

    ligands = []
    raw_ligands = settings.get("ligands", [])
    if isinstance(raw_ligands, dict):
        ligand_list = list(raw_ligands.get("items", []))
        if raw_ligands.get("smiles_file"):
            ligand_list.extend(parse_smiles_file(resolve_path(base, raw_ligands["smiles_file"])))
    else:
        ligand_list = list(raw_ligands or [])
    for ligand in ligand_list:
        if bool(ligand.get("smiles")) == bool(ligand.get("ccd")):
            raise ConfigError("Each ligand needs exactly one of 'smiles' or 'ccd'")
        raw_count = ligand.get("count", 1)
        if isinstance(raw_count, bool):
            raise ConfigError("Ligand 'count' must be a positive integer")
        try:
            copies = int(raw_count)
        except (TypeError, ValueError) as exc:
            raise ConfigError("Ligand 'count' must be a positive integer") from exc
        if copies < 1:
            raise ConfigError("Ligand 'count' must be a positive integer")
        requested_ids = ligand.get("chains")
        ids = allocate_chain_ids(requested_ids, copies, used) if requested_ids else [allocate_nonprotein_id(used) for _ in range(copies)]
        for cid in ids:
            item = {"chain": cid, "name": ligand.get("name", cid)}
            item["smiles" if ligand.get("smiles") else "ccd"] = ligand.get("smiles") or ligand.get("ccd")
            ligands.append(item)

    ions = []
    ion_cfg = settings.get("ions", {"mode": "none"})
    if ion_cfg.get("mode", "none") not in {"none", "explicit"}:
        raise ConfigError("ions.mode must be 'none' or 'explicit'")
    if ion_cfg.get("mode") == "explicit":
        for species in ion_cfg.get("species", []):
            code = str(species.get("ccd", "")).upper()
            if not code:
                raise ConfigError("Every explicit ion requires a CCD code")
            for _ in range(int(species.get("count", 1))):
                ions.append({"chain": allocate_nonprotein_id(used, "I"), "ccd": code,
                              "role": species.get("role", "structural"),
                              "recognized_common_ion": code in ION_CODES})

    glycans = []
    bonds = []
    for index, site in enumerate(glycosites, 1):
        gid = allocate_nonprotein_id(used, "G")
        site["site_id"] = f"site{index}"
        site["glycan_chain"] = gid
        glycans.append({"chain": gid, "ccd": site["residues"], "site_id": site["site_id"]})
        protein_atom, glycan_atom = sweetspot.protein_glycan_link_atoms(site["residues"])
        bonds.append({"atom1": [site["chain"], int(site["asn_chain_index"]), protein_atom],
                      "atom2": [gid, 1, glycan_atom], "kind": "protein-glycan"})
        for parent, child, atom1, atom2 in sweetspot.glycan_bonds(
            site["residues"], site["tree"], site["preset"]
        ):
            bonds.append({"atom1": [gid, parent + 1, atom1],
                          "atom2": [gid, child + 1, atom2], "kind": "glycan-glycan"})

    msa = normalize_msa_settings(settings, backend)
    warnings = []
    distinct_proteins = len({protein["sequence"] for protein in proteins})
    if backend == "boltz2" and msa["source"] == "local" and distinct_proteins > 1 and msa.get("pairing") != "paired":
        warnings.append(
            "Boltz2 received separate local A3Ms for a heteromer. These are unpaired; "
            "for interface-sensitive prediction prefer msa.source='server' or provide a "
            "Boltz paired CSV MSA."
        )
    prediction = dict(settings.get("prediction", {}))
    prediction.setdefault("diffusion_samples", 1)
    prediction.setdefault("max_parallel_samples", 1)
    prediction.setdefault("seed", 0)
    prediction.setdefault("num_loops", 20)
    prediction.setdefault("recycling_steps", 3)
    prediction.setdefault("sampling_steps", 100)
    if int(prediction["diffusion_samples"]) < 1:
        raise ConfigError("prediction.diffusion_samples must be at least 1")
    if int(prediction["max_parallel_samples"]) < 1:
        raise ConfigError("prediction.max_parallel_samples must be at least 1")
    return {
        "schema_version": 1, "name": settings.get("name", "sweetspot_job"),
        "backend": backend, "backend_settings": backend_cfg if isinstance(backend_cfg, dict) else {},
        "system_type": "monomer" if len(proteins) == 1 else ("homomer" if len({p['sequence'] for p in proteins}) == 1 else "heteromer"),
        "proteins": proteins, "glycosylation_sites": glycosites, "glycans": glycans,
        "ligands": ligands, "ions": ions, "bonds": bonds, "msa": msa,
        "prediction": prediction, "output": settings.get("output", {"directory": "results", "format": "mmcif"}),
        "warnings": warnings,
    }


def build_msa_plan(resolved: dict[str, Any]) -> list[dict[str, Any]]:
    representatives: OrderedDict[str, list[str]] = OrderedDict()
    for protein in resolved["proteins"]:
        representatives.setdefault(protein["sequence"], []).append(protein["chain"])
    plan = []
    paired_csv = (
        resolved["backend"] == "boltz2"
        and resolved["msa"]["source"] == "local"
        and resolved["msa"].get("pairing") == "paired"
        and len(representatives) > 1
    )
    for seq, chains in representatives.items():
        representative = chains[0]
        plan.append({"representative": representative, "chains": chains, "sequence": seq,
                     "path": f"msas/{representative}.{'csv' if paired_csv else 'a3m'}"})
    return plan


def backend_entities(
    resolved: dict[str, Any], *, af3: bool = False, group_boltz_copies: bool = False
) -> list[dict[str, Any]]:
    entities = []
    protein_groups: list[list[dict[str, Any]]] = []
    for protein in resolved["proteins"]:
        if not group_boltz_copies:
            protein_groups.append([protein])
            continue
        matching = next((group for group in protein_groups if (
            group[0]["entity_id"] == protein["entity_id"]
            and group[0]["sequence"] == protein["sequence"]
            and group[0]["modifications"] == protein["modifications"]
        )), None)
        if matching is None:
            protein_groups.append([protein])
        else:
            matching.append(protein)
    for group in protein_groups:
        protein = group[0]
        chain_ids = [item["chain"] for item in group]
        body: dict[str, Any] = {
            "id": chain_ids[0] if len(chain_ids) == 1 else chain_ids,
            "sequence": protein["sequence"],
        }
        if protein["modifications"]:
            if af3:
                body["modifications"] = [
                    {"ptmType": item["ccd"], "ptmPosition": item["position"]}
                    for item in protein["modifications"]
                ]
            else:
                body["modifications"] = protein["modifications"]
        if resolved["msa"]["source"] == "none":
            if af3:
                body["unpairedMsa"] = ""
                body["pairedMsa"] = ""
            else:
                body["msa"] = "empty"
        elif resolved["msa"]["source"] == "server":
            pass
        else:
            plan = next(x for x in resolved["msa_plan"] if protein["chain"] in x["chains"])
            if af3:
                body["unpairedMsaPath"] = plan["path"]
                body["pairedMsa"] = ""
            else:
                body["msa"] = plan["path"]
        entities.append({"protein": body})
    for item in resolved["glycans"] + resolved["ligands"] + resolved["ions"]:
        body = {"id": item["chain"]}
        if "smiles" in item:
            body["smiles"] = item["smiles"]
        else:
            code = item["ccd"]
            body["ccdCodes" if af3 else "ccd"] = code if isinstance(code, list) else ([code] if af3 else code)
        entities.append({"ligand": body})
    return entities


def write_backend_input(jobdir: Path, resolved: dict[str, Any]) -> str:
    backend = resolved["backend"]
    if backend == "boltz2":
        payload = {"version": 1, "sequences": backend_entities(resolved, group_boltz_copies=True),
                   "constraints": [{"bond": bond} for bond in (
                       {"atom1": x["atom1"], "atom2": x["atom2"]} for x in resolved["bonds"]
                   )]}
        if yaml is None:
            raise ConfigError("Boltz2 preparation requires PyYAML: pip install pyyaml")
        path = jobdir / "input.yaml"
    elif backend == "alphafold3":
        payload = {"name": resolved["name"], "modelSeeds": [resolved["prediction"]["seed"]],
                   "sequences": backend_entities(resolved, af3=True),
                   "bondedAtomPairs": [[x["atom1"], x["atom2"]] for x in resolved["bonds"]],
                   "dialect": "alphafold3", "version": 3}
        path = jobdir / "alphafold3_input.json"
    else:
        payload = {"sequences": [], "covalent_bonds_by_name": resolved["bonds"]}
        for protein in resolved["proteins"]:
            payload["sequences"].append({"type": "protein", "id": protein["chain"],
                                          "sequence": protein["sequence"],
                                          "modifications": protein["modifications"]})
        for item in resolved["glycans"] + resolved["ligands"] + resolved["ions"]:
            entry = {"type": "ligand", "id": item["chain"]}
            entry["smiles" if "smiles" in item else "ccd"] = item.get("smiles", item.get("ccd"))
            payload["sequences"].append(entry)
        path = jobdir / "esmfold2_input.json"
    if backend == "boltz2":
        path.write_text(yaml.safe_dump(payload, sort_keys=False))
    else:
        write_json(path, payload)
    return path.name


def shell_quote_command(command: list[str]) -> str:
    return " ".join(shlex.quote(str(x)) for x in command)


def merged_submission(submission: dict[str, Any], backend: str) -> dict[str, Any]:
    """Merge common submission values with optional per-backend overrides."""
    result = {key: value for key, value in submission.items() if key != "backends"}
    override = submission.get("backends", {}).get(backend, {})
    for section in ("slurm", "environment"):
        merged = dict(result.get(section, {}))
        merged.update(override.get(section, {}))
        result[section] = merged
    return result


def make_sbatch(submission: dict[str, Any], backend: str) -> str:
    submission = merged_submission(submission, backend)
    slurm = submission.get("slurm", submission)
    mapping = [
        ("job_name", "job-name"), ("cpus_per_task", "cpus-per-task"), ("nodes", "nodes"),
        ("memory", "mem"), ("partition", "partition"), ("gres", "gres"),
        ("nodelist", "nodelist"), ("exclude", "exclude"), ("time", "time"),
        ("output", "output"), ("account", "account"), ("qos", "qos"),
    ]
    lines = ["#!/bin/bash"]
    for key, directive in mapping:
        value = slurm.get(key)
        if value is not None and value != "":
            lines.append(f"#SBATCH --{directive}={value}")
    lines.extend(["", "set -euo pipefail", "", 'JOBDIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")" && pwd)}"', 'cd "$JOBDIR"'])
    env = submission.get("environment", {})
    if env.get("shell"):
        lines.append(f"export SHELL={shlex.quote(str(env['shell']))}")
    for key, value in env.get("variables", {}).items():
        if not re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*", key):
            raise ConfigError(f"Invalid environment variable name: {key}")
        lines.append(f"export {key}={shlex.quote(str(value))}")
    lines.extend(str(x) for x in env.get("setup", []))
    python_runner = env.get("python", "python")
    lines.extend(["", 'echo "Starting ${SLURM_JOB_NAME:-sweetspot} on $(hostname) at $(date)"',
                  f'{shlex.quote(str(python_runner))} "$JOBDIR/workflow.py" run --job-dir "$JOBDIR"',
                  'echo "Finished ${SLURM_JOB_NAME:-sweetspot} at $(date)"', ""])
    return "\n".join(lines)


def prepare(settings_path: Path, submission_path: Path | None, jobdir: Path) -> None:
    settings_path = settings_path.resolve()
    user_settings = read_json(settings_path)
    submission = read_json(submission_path.resolve()) if submission_path else {}
    backend_cfg = user_settings.get("backend", {})
    backend = backend_cfg if isinstance(backend_cfg, str) else backend_cfg.get("name")
    defaults = submission.get("workflow_defaults", {})
    backend_defaults = submission.get("backends", {}).get(backend, {}).get("workflow_defaults", {})
    settings = merge_settings(merge_settings({}, defaults), backend_defaults)
    mergeable_user_settings = dict(user_settings)
    mergeable_user_settings.setdefault("name", settings_path.stem)
    if isinstance(mergeable_user_settings.get("msa"), str):
        mergeable_user_settings["msa"] = {"source": mergeable_user_settings["msa"]}
    settings = merge_settings(settings, mergeable_user_settings)
    resolved = resolve_settings(settings, settings_path.parent)
    resolved["msa_plan"] = build_msa_plan(resolved)
    jobdir.mkdir(parents=True, exist_ok=True)
    (jobdir / "msas").mkdir(exist_ok=True)
    backend_input = write_backend_input(jobdir, resolved)
    resolved["backend_input"] = backend_input
    write_json(jobdir / "settings.json", user_settings)
    write_json(jobdir / "effective_settings.json", settings)
    write_json(jobdir / "submission.json", submission)
    write_json(jobdir / "resolved_settings.json", resolved)
    if resolved["msa"]["source"] == "local":
        with (jobdir / "msa_queries.fasta").open("w") as handle:
            paired = resolved["msa"].get("pairing") == "paired" and len(resolved["msa_plan"]) > 1
            if paired:
                handle.write(">sweetspot_complex\n")
                handle.write(":".join(item["sequence"] for item in resolved["msa_plan"]) + "\n")
            else:
                for item in resolved["msa_plan"]:
                    handle.write(f">{item['representative']}\n{item['sequence']}\n")
    shutil.copy2(Path(__file__).resolve(), jobdir / "workflow.py")
    sbatch = jobdir / "run.sbatch"
    sbatch.write_text(make_sbatch(submission, resolved["backend"]))
    sbatch.chmod(0o755)
    print(f"Prepared {resolved['system_type']} {resolved['backend']} job: {jobdir}")
    print(f"Proteins: {len(resolved['proteins'])}; glycans: {len(resolved['glycans'])}; ligands: {len(resolved['ligands'])}; ions: {len(resolved['ions'])}")
    for warning in resolved.get("warnings", []):
        print(f"WARNING: {warning}", file=sys.stderr)


def split_paired_a3m(source: Path, plan: list[dict[str, Any]], jobdir: Path) -> None:
    """Convert ColabFold's concatenated complex A3M into keyed Boltz CSVs."""
    records: list[tuple[str, str]] = []
    header = None
    chunks: list[str] = []
    for raw in source.read_text().splitlines():
        if raw.startswith(">"):
            if header is not None:
                records.append((header, "".join(chunks)))
            header, chunks = raw[1:].strip(), []
        elif raw.strip():
            chunks.append(raw.strip())
    if header is not None:
        records.append((header, "".join(chunks)))
    if not records:
        raise RuntimeError(f"Paired MSA is empty: {source}")
    lengths = [len(item["sequence"]) for item in plan]
    total = sum(lengths)
    rows: list[list[str]] = [[] for _ in plan]
    keys: list[int] = []
    for row_number, (record_header, sequence) in enumerate(records):
        aligned = "".join(char for char in sequence if not char.islower() and char != ":")
        if len(aligned) != total:
            raise RuntimeError(
                f"Paired MSA row {record_header!r} has {len(aligned)} aligned columns; expected {total}"
            )
        # Boltz stores CSV keys as integers; strings fail in parse_csv when
        # converted to its structured MSASequence dtype.
        key = row_number
        keys.append(key)
        start = 0
        for index, length in enumerate(lengths):
            rows[index].append(aligned[start:start + length])
            start += length
    for item, sequences in zip(plan, rows):
        with (jobdir / item["path"]).open("w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(["sequence", "key"])
            writer.writerows(zip(sequences, keys))


def run_local_msa(jobdir: Path, resolved: dict[str, Any]) -> None:
    msa = resolved["msa"]
    if msa["source"] == "none":
        return
    if msa["source"] == "server":
        if resolved["backend"] == "boltz2":
            return
        raise ConfigError("msa.source='server' is currently delegated only to Boltz2; use local or none for this backend")
    expected = [jobdir / item["path"] for item in resolved["msa_plan"]]
    if expected and all(path.is_file() and path.stat().st_size > 0 for path in expected):
        print("Reusing existing local MSA files: " + ", ".join(str(path) for path in expected))
        return
    command = msa.get("command", "colabfold_search")
    mmseqs = msa.get("mmseqs", "mmseqs")
    database = msa.get("database")
    database_dir = msa.get("database_directory")
    if not database or not database_dir:
        raise ConfigError("Local MSA generation requires msa.database and msa.database_directory")
    tmp_root = Path(msa.get("tmp_directory", os.environ.get("TMPDIR", "/tmp")))
    tmp_root.mkdir(parents=True, exist_ok=True)
    run_dir = Path(tempfile.mkdtemp(prefix="sweetspot_msa_", dir=tmp_root))
    try:
        cmd = [command, "--threads", str(msa.get("threads", 16)), "--mmseqs", mmseqs,
               "--use-env", "1" if msa.get("use_environmental_database", False) else "0",
               str(jobdir / "msa_queries.fasta"), "--db1", database, database_dir, str(run_dir)]
        subprocess.run(cmd, check=True)
        paired = msa.get("pairing") == "paired" and len(resolved["msa_plan"]) > 1
        if paired:
            source = run_dir / "0.a3m"
            if not source.exists():
                raise RuntimeError(f"colabfold_search did not create paired complex MSA {source}")
            split_paired_a3m(source, resolved["msa_plan"], jobdir)
        else:
            for index, item in enumerate(resolved["msa_plan"]):
                source = run_dir / f"{index}.a3m"
                if not source.exists():
                    raise RuntimeError(f"colabfold_search did not create {source}")
                shutil.move(str(source), jobdir / item["path"])
    finally:
        shutil.rmtree(run_dir, ignore_errors=True)


def run_boltz(jobdir: Path, resolved: dict[str, Any]) -> None:
    cfg = resolved["backend_settings"]
    command = [cfg.get("executable", "boltz"), "predict", resolved["backend_input"],
               "--out_dir", resolved["output"].get("directory", "results"),
               "--diffusion_samples", str(resolved["prediction"]["diffusion_samples"]),
               "--max_parallel_samples", str(resolved["prediction"]["max_parallel_samples"]),
               "--sampling_steps", str(resolved["prediction"]["sampling_steps"]),
               "--recycling_steps", str(resolved["prediction"]["recycling_steps"]),
               "--output_format", str(resolved["output"].get("format", "mmcif"))]
    if cfg.get("use_potentials", False):
        command.append("--use_potentials")
    if cfg.get("override", True):
        command.append("--override")
    if cfg.get("step_scale") is not None:
        command.extend(["--step_scale", str(cfg["step_scale"])])
    if resolved["msa"]["source"] == "server":
        command.append("--use_msa_server")
        if resolved["msa"].get("server_url"):
            command.extend(["--msa_server_url", resolved["msa"]["server_url"]])
    subprocess.run(command, cwd=jobdir, check=True)
    outdir = jobdir / resolved["output"].get("directory", "results")
    if not any(outdir.rglob("*.cif")) and not any(outdir.rglob("*.mmcif")):
        raise RuntimeError(f"Boltz2 completed without producing a CIF/mmCIF under {outdir}")


def run_af3(jobdir: Path, resolved: dict[str, Any]) -> None:
    cfg = resolved["backend_settings"]
    executable = cfg.get("executable", "run_alphafold.py")
    command = [executable, f"--json_path={resolved['backend_input']}",
               f"--output_dir={resolved['output'].get('directory', 'results')}"]
    for key in ("model_dir", "db_dir"):
        if cfg.get(key):
            command.append(f"--{key}={cfg[key]}")
    subprocess.run(command, cwd=jobdir, check=True)


def resolve_cached_esm_model(repo: str, extra_roots: list[str] | None = None) -> str:
    """Prefer a complete local Hugging Face snapshot, including the site cache."""
    roots = [Path(os.environ.get("HF_HOME", "/var/tmp/hf_cache")), Path("/z/bio/caches/esmfold2/hf_cache")]
    roots.extend(Path(x) for x in (extra_roots or []))
    variants = {repo}
    if "/" in repo:
        owner, name = repo.split("/", 1)
        variants.update({f"{owner.lower()}/{name}", f"{owner.capitalize()}/{name}"})
    for root in roots:
        for variant in variants:
            cache = root / f"models--{variant.replace('/', '--')}"
            snapshots = cache / "snapshots"
            candidates = []
            ref = cache / "refs/main"
            if ref.is_file():
                candidates.append(snapshots / ref.read_text().strip())
            if snapshots.is_dir():
                candidates.extend(sorted(snapshots.iterdir(), reverse=True))
            for candidate in candidates:
                if candidate.is_dir() and (
                    (candidate / "model.safetensors").is_file()
                    or (candidate / "model.safetensors.index.json").is_file()
                ):
                    return str(candidate)
    return repo


def run_esmfold2(jobdir: Path, resolved: dict[str, Any]) -> None:
    # Imports stay runtime-only so preparation does not require GPU dependencies.
    import torch  # type: ignore
    from esm.models.esmfold2 import CovalentBond, ESMFold2InputBuilder  # type: ignore
    from esm.utils.msa import MSA  # type: ignore
    from esm.utils.structure.input_builder import (  # type: ignore
        LigandInput, Modification, ProteinInput, StructurePredictionInput,
    )
    from transformers.models.esmfold2.modeling_esmfold2 import ESMFold2Model  # type: ignore
    from esm.models.esmfold2.prepare_input import build_chains_from_input  # type: ignore

    data = read_json(jobdir / resolved["backend_input"])
    sequences = []
    for item in data["sequences"]:
        if item["type"] == "protein":
            mods = [Modification(position=int(x["position"]) - 1, ccd=x["ccd"]) for x in item.get("modifications", [])]
            msa_obj = None
            if resolved["msa"]["source"] != "none":
                plan = next(x for x in resolved["msa_plan"] if item["id"] in x["chains"])
                msa_obj = MSA.from_a3m(
                    str(jobdir / plan["path"]), remove_insertions=True,
                    max_sequences=int(resolved["msa"].get("max_sequences", 256)),
                )
            sequences.append(ProteinInput(id=item["id"], sequence=item["sequence"], modifications=mods or None, msa=msa_obj))
        else:
            ccd = item.get("ccd")
            sequences.append(LigandInput(id=item["id"], ccd=ccd if isinstance(ccd, list) else ([ccd] if ccd else None), smiles=item.get("smiles")))
    named = data.get("covalent_bonds_by_name", [])
    provisional = [CovalentBond(chain_id1=x["atom1"][0], res_idx1=int(x["atom1"][1]) - 1, atom_idx1=0,
                                chain_id2=x["atom2"][0], res_idx2=int(x["atom2"][1]) - 1, atom_idx2=0) for x in named]
    spi = StructurePredictionInput(sequences=sequences, covalent_bonds=provisional or None)
    chains, tokens, atoms = build_chains_from_input(spi, seed=int(resolved["prediction"]["seed"]))
    chain_by_asym = {c.asym_id: c.chain_id for c in chains}
    residue_atoms: dict[tuple[str, int], list[Any]] = {}
    for atom in atoms:
        token = tokens[atom.token_index]
        residue_atoms.setdefault((chain_by_asym[token.asym_id], token.residue_index), []).append(atom)
    bonds = []
    for item in named:
        refs = []
        for chain, residue, atom_name in (item["atom1"], item["atom2"]):
            atom_list = residue_atoms.get((chain, int(residue) - 1), [])
            index = next((i for i, atom in enumerate(atom_list) if atom.name == atom_name), None)
            if index is None:
                raise ConfigError(f"ESMFold2 could not resolve atom {chain}:{residue}:{atom_name}")
            refs.append((chain, int(residue) - 1, index))
        bonds.append(CovalentBond(chain_id1=refs[0][0], res_idx1=refs[0][1], atom_idx1=refs[0][2],
                                  chain_id2=refs[1][0], res_idx2=refs[1][1], atom_idx2=refs[1][2]))
    spi.covalent_bonds = bonds or None
    cfg = resolved["backend_settings"]
    choice = cfg.get("model", "quality").lower()
    model_names = {"fast": "Biohub/ESMFold2-Fast", "quality": "Biohub/ESMFold2", "hq-msa": "Biohub/ESMFold2"}
    model_name = model_names.get(choice, cfg.get("model"))
    if not model_name:
        raise ConfigError("No ESMFold2 model was selected")
    dtype_name = cfg.get("dtype", "float32")
    dtype = torch.bfloat16 if dtype_name == "bfloat16" else torch.float32
    device = cfg.get("device", "cuda" if torch.cuda.is_available() else "cpu")
    cache_roots = cfg.get("model_cache_roots", [])
    model_source = resolve_cached_esm_model(model_name, cache_roots)
    print(f"Loading ESMFold2 {model_name} from {model_source}")
    model = ESMFold2Model.from_pretrained(model_source, dtype=dtype, load_esmc=False).to(device=device, dtype=dtype).eval()
    if os.environ.get("ESMFOLD2_SKIP_ESMC", "0") != "1":
        precision = cfg.get("esmc_precision", "bf16" if dtype == torch.bfloat16 else "fp32")
        esmc_source = resolve_cached_esm_model(model.config.esmc_id, cache_roots)
        print(f"Loading ESMC from {esmc_source}")
        model.load_esmc(esmc_source, precision=precision)
    chunk_size = cfg.get("chunk_size", 64)
    if hasattr(model, "set_chunk_size"):
        model.set_chunk_size(None if chunk_size is None or str(chunk_size).lower() == "none" else int(chunk_size))
    outdir = jobdir / resolved["output"].get("directory", "results")
    outdir.mkdir(parents=True, exist_ok=True)
    total_samples = int(resolved["prediction"]["diffusion_samples"])
    parallel_samples = min(int(resolved["prediction"]["max_parallel_samples"]), total_samples)
    base_seed = int(resolved["prediction"]["seed"])
    metrics_all = []
    completed = 0
    while completed < total_samples:
        batch_size = min(parallel_samples, total_samples - completed)
        batch_seed = base_seed + completed
        print(f"Generating samples {completed + 1}-{completed + batch_size} of {total_samples} (seed {batch_seed})")
        with torch.inference_mode():
            batch = ESMFold2InputBuilder().fold(
                model, spi, num_loops=int(resolved["prediction"]["num_loops"]),
                num_sampling_steps=int(resolved["prediction"]["sampling_steps"]),
                num_diffusion_samples=batch_size,
                seed=batch_seed,
            )
        results = batch if isinstance(batch, list) else [batch]
        if len(results) != batch_size:
            raise RuntimeError(f"ESMFold2 returned {len(results)} result(s) for a requested batch of {batch_size}")
        for offset, result in enumerate(results):
            sample_number = completed + offset + 1
            output_name = f"{resolved['name']}_model_{sample_number:03d}.cif"
            (outdir / output_name).write_text(result.complex.to_mmcif())
            metrics = {
                "sample": sample_number,
                "seed": batch_seed + offset,
                "file": output_name,
            }
            metrics.update({key: float(getattr(result, key)) for key in ("ptm", "iptm") if getattr(result, key, None) is not None})
            if getattr(result, "plddt", None) is not None:
                metrics["mean_plddt"] = float(result.plddt.mean())
            metrics_all.append(metrics)
        completed += batch_size
        del batch, results
        if device.startswith("cuda"):
            torch.cuda.empty_cache()
    write_json(outdir / "confidence.json", metrics_all)


def run_job(jobdir: Path) -> None:
    jobdir = jobdir.resolve()
    resolved = read_json(jobdir / "resolved_settings.json")
    run_local_msa(jobdir, resolved)
    if resolved["backend"] == "boltz2":
        run_boltz(jobdir, resolved)
    elif resolved["backend"] == "alphafold3":
        run_af3(jobdir, resolved)
    else:
        run_esmfold2(jobdir, resolved)


def default_submission_path(settings_path: Path) -> Path | None:
    """Find a job-local profile first, then the repository-wide profile."""
    candidates = [
        settings_path.resolve().parent / "submission.json",
        Path(__file__).resolve().parent / "submission.json",
    ]
    return next((path for path in candidates if path.is_file()), None)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    prep = sub.add_parser("prepare", help="Resolve JSON and create a portable job directory")
    prep.add_argument("--settings", required=True, type=Path)
    prep.add_argument(
        "--submission", type=Path,
        help=("Reusable Slurm/profile JSON (default: submission.json beside "
              "--settings, then repository root)"),
    )
    prep.add_argument("--outdir", required=True, type=Path)
    run = sub.add_parser("run", help="Run a prepared job locally or inside Slurm")
    run.add_argument("--job-dir", required=True, type=Path)
    args = parser.parse_args()
    try:
        if args.command == "prepare":
            submission = args.submission
            if submission is None:
                submission = default_submission_path(args.settings)
            prepare(args.settings, submission, args.outdir.resolve())
        else:
            run_job(args.job_dir)
    except (ConfigError, ValueError) as exc:
        raise SystemExit(f"error: {exc}") from exc


if __name__ == "__main__":
    main()
