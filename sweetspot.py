#!/usr/bin/env python3
"""
sweetspot.py

Detect N-glycosylation sequons in PDB/mmCIF and write inputs for:
  - RFAA / RF3-style RoseTTAFold-All-Atom
  - AlphaFold3
  - Boltz2
  - ESMFold/ESM2 protein-only FASTA

AF3/Boltz2/RFAA outputs are aggregate jobs: one input file contains all detected
N-glycosylation sites, each with its own glycan ligand.

Usage:
  python sweetspot.py -f protein.cif
  python sweetspot.py -f protein.pdb -m all
  python sweetspot.py -f protein.pdb -m all -t branched -g NAG-NAG-MAN-MAN
  python sweetspot.py -f protein.pdb -m boltz2 -g man9
  python sweetspot.py -f protein.pdb -m rfaa --rfaa-glycan-sdf glycan.sdf
"""

from __future__ import annotations

import argparse
import json
import re
import shutil
from pathlib import Path
from typing import Dict, List, Tuple

from Bio.PDB import PDBParser, MMCIFParser
from Bio.PDB.Polypeptide import protein_letters_3to1

try:
    import yaml
except ImportError:
    yaml = None


DEFAULT_GLYCAN = "NAG-NAG-MAN"

# Man9GlcNAc2 high-mannose preset.
# Residue order / IDs inside each ligand chain:
#   1 NAG  ASN-linked reducing GlcNAc
#   2 NAG  second GlcNAc
#   3 MAN  core beta-Man
#   4 MAN  alpha1-3 arm Man
#   5 MAN  alpha1-6 arm Man
#   6 MAN  D1 arm Man
#   7 MAN  terminal D1 Man
#   8 MAN  D2 arm Man
#   9 MAN  terminal D2 Man
#  10 MAN  D3 arm Man
#  11 MAN  terminal D3 Man
GLYCAN_PRESETS = {
    "man9": [
        "NAG", "NAG",
        "MAN",
        "MAN", "MAN",
        "MAN", "MAN",
        "MAN", "MAN",
        "MAN", "MAN",
    ],
    "man9glcnac2": [
        "NAG", "NAG",
        "MAN",
        "MAN", "MAN",
        "MAN", "MAN",
        "MAN", "MAN",
        "MAN", "MAN",
    ],
}

# Explicit covalent topology for Man9GlcNAc2. Atom naming assumes CCD-style
# sugar atoms where the child anomeric carbon is C1.
# Linkage labels are comments only; atom pairs define the generated constraints.
GLYCAN_PRESET_BONDS = {
    "man9": [
        (0, 1, "O4", "C1"),   # GlcNAc beta1-4 GlcNAc
        (1, 2, "O4", "C1"),   # GlcNAc beta1-4 Man core
        (2, 3, "O3", "C1"),   # Man alpha1-3 arm
        (2, 4, "O6", "C1"),   # Man alpha1-6 arm
        (3, 5, "O2", "C1"),   # D1 extension
        (5, 6, "O2", "C1"),   # D1 terminal
        (4, 7, "O3", "C1"),   # D2 extension
        (7, 8, "O2", "C1"),   # D2 terminal
        (4, 9, "O6", "C1"),   # D3 extension
        (9, 10, "O2", "C1"),  # D3 terminal
    ],
    "man9glcnac2": [
        (0, 1, "O4", "C1"),
        (1, 2, "O4", "C1"),
        (2, 3, "O3", "C1"),
        (2, 4, "O6", "C1"),
        (3, 5, "O2", "C1"),
        (5, 6, "O2", "C1"),
        (4, 7, "O3", "C1"),
        (7, 8, "O2", "C1"),
        (4, 9, "O6", "C1"),
        (9, 10, "O2", "C1"),
    ],
}

# Synced from /Users/joe/Software/mojoeMD/mojoeMD.py glycoprotein helpers.
PDB_TO_GLYCAM_SUFFIXES = {
    "ABE": ("AE", "AF"),
    "ADA": ("OA", "OB"),
    "AHR": ("AA", "AB", "AD", "AU"),
    "ALL": ("NA", "NB"),
    "ALT": ("EA", "EB"),
    "ARA": ("AA", "AB", "AD", "AU"),
    "BAC": ("BC",),
    "BGC": ("GB",),
    "BMA": ("MB",),
    "FRU": ("CA", "CB", "CD", "CU"),
    "FUC": ("FA",),
    "FUL": ("FA",),
    "GAL": ("LA", "LB"),
    "GCS": ("YN", "YS"),
    "GCU": ("ZA", "ZB"),
    "GLA": ("LA", "LB"),
    "GLC": ("GA", "GB"),
    "GL0": ("KA", "KB"),
    "GUL": ("KA", "KB"),
    "GUP": ("KA", "KB"),
    "IDS": ("UA", "UB"),
    "LYX": ("DA", "DB", "DD", "DU"),
    "LXC": ("XA", "XB", "XD", "XU"),
    "MAL": ("GA",),
    "MAN": ("MA",),
    "NAG": ("YB",),
    "NDG": ("YB",),
    "NGA": ("VA", "VB"),
    "QUI": ("QA", "QB"),
    "PSI": ("PA", "PB", "PD", "PU"),
    "RAM": ("HA", "HB"),
    "RIB": ("RA", "RB", "RD", "RU"),
    "SGN": ("YB",),
    "SIA": ("SA", "SB"),
    "SOR": ("BA", "BB", "BD", "BU"),
    "TAG": ("JA", "JB", "JD", "JU"),
    "TAL": ("TA", "TB"),
    "TYV": ("TV", "Tv"),
    "XYP": ("XA", "XB", "XD", "XU"),
    "XYS": ("XA", "XB", "XD", "XU"),
    "XYL": ("XA", "XB", "XD", "XU"),
}

GLYCAN_RESNAMES = set(PDB_TO_GLYCAM_SUFFIXES) | {
    "API", "BEM", "COL", "DHA", "DIG", "GMH", "IDO", "KDN", "KDO", "LGU", "MAV", "MUR", "NEU",
    "OLI", "PAR",
}

PROTEIN_GLYCAN_LINKS = {
    ("ASN", "ND2", "NAG", "C1"),
    ("ASN", "ND2", "NDG", "C1"),
    ("ASN", "ND2", "BMA", "C1"),
}

PROTEIN_GLYCAN_LINK_ATOMS = {
    glycan_res: (protein_atom, glycan_atom)
    for _, protein_atom, glycan_res, glycan_atom in PROTEIN_GLYCAN_LINKS
}

GLYCAN_RING_ATOM = "C1"

GLYCAN_ANOMERIC_ATOMS = {
    "SIA": "C2",
    "NEU": "C2",
    "KDN": "C2",
    "KDO": "C2",
}

GLYCAN_LINK_ATOMS = {
    ("NAG", "NAG"): ("O4", "C1"),
    ("NAG", "MAN"): ("O4", "C1"),
    ("NAG", "FUC"): ("O6", "C1"),
    ("MAN", "MAN"): ("O3", "C1"),
}

CHAIN_ID_POOL = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"


GLYCAN_LINK_MANIFEST_NAME = "linked_glycans.json"
GLYCAN_SANITIZER_NAME = "sanitize_linked_glycans.py"


def parse_args():
    p = argparse.ArgumentParser(
        description="Detect N-glycosylation sites and generate AF3/Boltz2/RFAA/ESM inputs."
    )
    p.add_argument("-f", "--file", help="Input PDB/mmCIF file.")
    p.add_argument(
        "-m",
        "--model",
        default="rfaa",
        choices=["rfaa", "rf3", "rosettafold_all_atom", "af3", "boltz2", "esm", "esmfold", "esm2", "all"],
        help="Backend input to generate. Default: rfaa.",
    )
    p.add_argument(
        "-g",
        "--glycans",
        default=DEFAULT_GLYCAN,
        help="Glycan string, e.g. NAG-NAG-MAN or NAG-NAG-MAN-MAN.",
    )
    p.add_argument(
        "-t",
        "--tree",
        default="linear",
        choices=["linear", "branched"],
        help="linear = NAG-NAG-MAN; branched = branch after second sugar.",
    )
    p.add_argument(
        "--rfaa-glycan-sdf",
        "--glycan-sdf",
        dest="rfaa_glycan_sdf",
        default="glycan.sdf",
        help=(
            "Merged glycan SDF for RFAA/RF3-style covalent input. "
            "If it exists, it is copied next to each rfaa.yaml."
        ),
    )
    p.add_argument(
        "--rfaa-link-atom-index",
        type=int,
        default=1,
        help="1-indexed atom in the merged RFAA glycan SDF that bonds to ASN ND2. Default: 1.",
    )
    p.add_argument(
        "--rfaa-link-chirality",
        default="CW",
        choices=["CW", "CCW", "null"],
        help="RFAA chirality field for the glycan link atom. Default follows the RFAA glycan example: CW.",
    )
    p.add_argument(
        "-o",
        "--outdir",
        default=None,
        help="Output dir. Default: <input_stem>_glyco_inputs.",
    )
    p.add_argument(
        "--list-glycans",
        action="store_true",
        help="List supported glycan CCD/PDB residue names imported from mojoeMD and exit.",
    )
    return p.parse_args()


def normalize_model(model: str) -> str:
    if model in {"rf3", "rosettafold_all_atom"}:
        return "rfaa"
    if model in {"esmfold", "esm2"}:
        return "esm"
    return model


def parse_structure(path: Path):
    if path.suffix.lower() in {".cif", ".mmcif"}:
        return MMCIFParser(QUIET=True).get_structure(path.stem, str(path))
    return PDBParser(QUIET=True).get_structure(path.stem, str(path))


def aa3_to_aa1(resname: str) -> str:
    return protein_letters_3to1.get(resname.upper().strip(), "X")


def residue_label(residue) -> str:
    _, resseq, icode = residue.id
    return f"{resseq}{icode.strip()}" if icode.strip() else str(resseq)


def residue_int(label: str) -> int:
    m = re.match(r"\d+", str(label))
    if not m:
        raise ValueError(f"Cannot extract residue number from {label}")
    return int(m.group())


def extract_chains(structure) -> Dict[str, List[Tuple[str, str]]]:
    chains = {}
    model = next(structure.get_models())

    for chain in model:
        cid = chain.id.strip() or "A"
        residues = []

        for res in chain:
            hetflag, _, _ = res.id
            if hetflag.strip():
                continue

            aa = aa3_to_aa1(res.resname)
            if aa != "X":
                residues.append((residue_label(res), aa))

        if residues:
            chains[cid] = residues

    return chains


def find_nglyc_sites(chains):
    sites = []

    for cid, residues in chains.items():
        seq = "".join(aa for _, aa in residues)

        for i in range(len(seq) - 2):
            tri = seq[i:i + 3]
            if tri[0] == "N" and tri[1] != "P" and tri[2] in {"S", "T"}:
                sites.append({
                    "site_id": f"site{len(sites) + 1}",
                    "chain": cid,
                    "asn_resid": residues[i][0],
                    "asn_chain_index": i + 1,
                    "sequon": tri,
                    "x_resid": residues[i + 1][0],
                    "st_resid": residues[i + 2][0],
                })

    return sites


def glycan_preset_name(glycan: str) -> str | None:
    key = glycan.strip().lower().replace("-", "").replace("_", "")
    aliases = {
        "man9": "man9",
        "man9glcnac2": "man9",
        "highmannose9": "man9",
        "highman9": "man9",
    }
    return aliases.get(key)


def parse_glycan(glycan: str) -> List[str]:
    preset = glycan_preset_name(glycan)
    if preset:
        return list(GLYCAN_PRESETS[preset])
    return [x.strip().upper() for x in re.split(r"[-,>]+", glycan) if x.strip()]


def require_glycan(glycan: List[str]):
    if not glycan:
        raise ValueError("Glycan cannot be empty. Example: -g NAG-NAG-MAN")


def validate_glycan(glycan: List[str]):
    require_glycan(glycan)

    unsupported = [res for res in glycan if res not in GLYCAN_RESNAMES]
    if unsupported:
        supported = ", ".join(sorted(GLYCAN_RESNAMES))
        raise ValueError(
            f"Unsupported glycan residue code(s): {', '.join(unsupported)}. "
            f"Use --list-glycans to inspect supported names. Supported: {supported}"
        )

    first = glycan[0]
    if first not in PROTEIN_GLYCAN_LINK_ATOMS:
        allowed = ", ".join(sorted(PROTEIN_GLYCAN_LINK_ATOMS))
        raise ValueError(
            f"N-glycan attachment to ASN is only configured for first residue {allowed}; got {first}."
        )


def glycan_anomeric_atom(resname: str) -> str:
    return GLYCAN_ANOMERIC_ATOMS.get(resname, GLYCAN_RING_ATOM)


def protein_glycan_link_atoms(glycan: List[str]) -> Tuple[str, str]:
    return PROTEIN_GLYCAN_LINK_ATOMS[glycan[0]]


def glycan_link_atoms(parent_res: str, child_res: str) -> Tuple[str, str]:
    return GLYCAN_LINK_ATOMS.get((parent_res, child_res), ("O4", glycan_anomeric_atom(child_res)))


def glycan_edges(glycan: List[str], tree: str):
    if len(glycan) <= 1:
        return []

    if tree == "linear":
        return [(i, i + 1) for i in range(len(glycan) - 1)]

    edges = [(0, 1)]
    edges.extend((1, i) for i in range(2, len(glycan)))
    return edges


def glycan_bonds(glycan: List[str], tree: str, preset_name: str | None = None):
    if preset_name and preset_name in GLYCAN_PRESET_BONDS:
        return list(GLYCAN_PRESET_BONDS[preset_name])

    bonds = []
    for parent, child in glycan_edges(glycan, tree):
        parent_res = glycan[parent]
        child_res = glycan[child]
        atom1, atom2 = glycan_link_atoms(parent_res, child_res)
        bonds.append((parent, child, atom1, atom2))
    return bonds


def atom_ref(chain: str, residue: int, atom: str, role: str, resname: str | None = None) -> dict:
    ref = {
        "chain": chain,
        "residue": residue,
        "atom": atom,
        "role": role,
    }
    if resname:
        ref["resname"] = resname
    return ref


def build_linked_glycan_manifest(sites, glycan, tree, name: str, backend: str, preset_name: str | None = None) -> dict:
    bonds = []
    c1_linked = set()
    protein_atom, glycan_atom = protein_glycan_link_atoms(glycan)

    for site in sites:
        gly_id = site["glycan_id"]
        bonds.append({
            "atom1": atom_ref(site["chain"], site["asn_chain_index"], protein_atom, "protein", "ASN"),
            "atom2": atom_ref(gly_id, 1, glycan_atom, "glycan", glycan[0]),
            "kind": "protein-glycan",
        })
        if glycan_atom == glycan_anomeric_atom(glycan[0]):
            c1_linked.add((gly_id, 1, glycan[0]))

        for parent, child, atom1, atom2 in glycan_bonds(glycan, tree, preset_name):
            parent_res = glycan[parent]
            child_res = glycan[child]
            bonds.append({
                "atom1": atom_ref(gly_id, parent + 1, atom1, "glycan", parent_res),
                "atom2": atom_ref(gly_id, child + 1, atom2, "glycan", child_res),
                "kind": "glycan-glycan",
            })
            if atom1 == glycan_anomeric_atom(parent_res):
                c1_linked.add((gly_id, parent + 1, parent_res))
            if atom2 == glycan_anomeric_atom(child_res):
                c1_linked.add((gly_id, child + 1, child_res))

    remove_atoms = [
        {
            "chain": chain,
            "residue": residue,
            "resname": resname,
            "atom": "O1",
            "reason": "Remove anomeric O1 from residues whose anomeric carbon is externally bonded.",
        }
        for chain, residue, resname in sorted(c1_linked)
    ]

    return {
        "version": 1,
        "name": f"{name}_all_nglyc_sites",
        "backend": backend,
        "glycan": "-".join(glycan),
        "tree": tree,
        "preset": preset_name,
        "notes": [
            "Prediction inputs use CCD monomers plus covalent constraints.",
            "Many predictors emit standalone CCD sugar atom sets; this manifest condenses linked sugars for GLYCAM/OpenMM handoff.",
            "Residue numbers are prediction-input sequence indices, not necessarily source PDB residue labels.",
        ],
        "bonds": bonds,
        "remove_atoms": remove_atoms,
    }


def write_linked_glycan_manifest(sites, glycan, tree, outdir: Path, name: str, backend: str, preset_name: str | None = None) -> Path:
    manifest = build_linked_glycan_manifest(sites, glycan, tree, name, backend, preset_name)
    path = outdir / GLYCAN_LINK_MANIFEST_NAME
    path.write_text(json.dumps(manifest, indent=2) + "\n")
    return path


def write_glycan_sanitizer(outdir: Path) -> Path:
    script = outdir / GLYCAN_SANITIZER_NAME
    script.write_text("""#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path


def atom_key_from_line(line: str):
    return (line[21].strip() or " ", int(line[22:26]), line[12:16].strip())


def atom_serial_from_line(line: str) -> int:
    return int(line[6:11])


def renumber_atom_line(line: str, serial: int) -> str:
    return f"{line[:6]}{serial:5d}{line[11:]}"


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Condense predictor-emitted CCD glycans into GLYCAM/OpenMM-ready linked glycans."
    )
    parser.add_argument("input_pdb", help="Boltz/AF3/RFAA PDB to sanitize.")
    parser.add_argument("-m", "--manifest", default="linked_glycans.json", help="Linked glycan manifest from sweetspot.py.")
    parser.add_argument("-o", "--output", required=True, help="Sanitized PDB path.")
    args = parser.parse_args()

    manifest = json.loads(Path(args.manifest).read_text())
    remove = {
        (entry["chain"], int(entry["residue"]), entry["atom"])
        for entry in manifest.get("remove_atoms", [])
    }

    atom_lines = []
    other_lines = []
    old_conect_pairs = set()
    serial_to_key = {}
    removed_serials = set()

    for raw in Path(args.input_pdb).read_text().splitlines():
        line = raw.rstrip("\\n")
        record = line[:6].strip()
        if record in {"ATOM", "HETATM"}:
            serial = atom_serial_from_line(line)
            key = atom_key_from_line(line)
            serial_to_key[serial] = key
            if key in remove:
                removed_serials.add(serial)
                continue
            atom_lines.append((serial, line, key))
        elif record == "CONECT":
            fields = line.split()
            if len(fields) >= 3:
                src = int(fields[1])
                for dst_text in fields[2:]:
                    dst = int(dst_text)
                    if src != dst:
                        old_conect_pairs.add(tuple(sorted((src, dst))))
        elif record not in {"END", "TER"}:
            other_lines.append(line)

    old_to_new = {}
    key_to_new = {}
    renumbered_atoms = []
    for new_serial, (old_serial, line, key) in enumerate(atom_lines, 1):
        old_to_new[old_serial] = new_serial
        key_to_new[key] = new_serial
        renumbered_atoms.append(renumber_atom_line(line, new_serial))

    conect_pairs = set()
    for src, dst in old_conect_pairs:
        if src in removed_serials or dst in removed_serials:
            continue
        if src in old_to_new and dst in old_to_new:
            conect_pairs.add(tuple(sorted((old_to_new[src], old_to_new[dst]))))

    missing_manifest_bonds = []
    for bond in manifest.get("bonds", []):
        k1 = (bond["atom1"]["chain"], int(bond["atom1"]["residue"]), bond["atom1"]["atom"])
        k2 = (bond["atom2"]["chain"], int(bond["atom2"]["residue"]), bond["atom2"]["atom"])
        s1 = key_to_new.get(k1)
        s2 = key_to_new.get(k2)
        if s1 is None or s2 is None:
            missing_manifest_bonds.append((k1, k2))
            continue
        conect_pairs.add(tuple(sorted((s1, s2))))

    grouped = defaultdict(list)
    for src, dst in sorted(conect_pairs):
        grouped[src].append(dst)
        grouped[dst].append(src)

    output = []
    output.extend(other_lines)
    output.extend(renumbered_atoms)
    for src in sorted(grouped):
        partners = sorted(set(grouped[src]))
        for i in range(0, len(partners), 4):
            chunk = partners[i:i + 4]
            output.append("CONECT" + f"{src:5d}" + "".join(f"{dst:5d}" for dst in chunk))
    output.append("END")
    Path(args.output).write_text("\\n".join(output) + "\\n")

    print(f"Wrote {args.output}")
    if remove:
        print(f"Removed {len(remove)} linked-glycan leaving-group atom definition(s): " + ", ".join(f"{c}:{r}:{a}" for c, r, a in sorted(remove)))
    if missing_manifest_bonds:
        print("WARNING: Some manifest bonds could not be mapped onto the PDB:")
        for k1, k2 in missing_manifest_bonds:
            print(f"  {k1} -- {k2}")


if __name__ == "__main__":
    main()
""")
    script.chmod(0o755)
    return script


def assign_glycan_chain_ids(sites, chains):
    used_chain_ids = set(chains)
    available = [cid for cid in CHAIN_ID_POOL if cid not in used_chain_ids]

    if len(sites) > len(available):
        raise ValueError(
            f"Detected {len(sites)} glycosites but only {len(available)} PDB-safe "
            "single-character glycan chain IDs are available."
        )

    for site, glycan_id in zip(sites, available):
        site["glycan_id"] = glycan_id


def rfaa_sdf_for_dir(glycan_sdf: Path, outdir: Path) -> Tuple[str, str]:
    if glycan_sdf.exists():
        target = outdir / glycan_sdf.name
        source = glycan_sdf.resolve()
        if target.resolve() != source:
            shutil.copy2(source, target)
        return target.name, f"Copied RFAA glycan SDF from {source}."

    return (
        glycan_sdf.name,
        (
            f"RFAA glycan SDF was not found at {glycan_sdf}. "
            f"Place a merged glycan SDF named {glycan_sdf.name} in this directory "
            "or edit sm_inputs.B.input before running."
        ),
    )


def print_supported_glycans():
    print("Supported glycan residue names from mojoeMD:")
    for i, resname in enumerate(sorted(GLYCAN_RESNAMES), 1):
        suffixes = PDB_TO_GLYCAM_SUFFIXES.get(resname)
        detail = f" GLYCAM families: {', '.join(suffixes)}" if suffixes else " detected only"
        print(f"  {resname}{detail}")
        if i % 999 == 0:
            print()

    anchored = ", ".join(sorted(PROTEIN_GLYCAN_LINK_ATOMS))
    print(f"\nConfigured ASN-linked first residues: {anchored}")
    print("Default glycan-glycan linkage fallback: parent O4 to child anomeric atom.")


def write_fasta(chains, path: Path):
    with path.open("w") as f:
        for cid, residues in chains.items():
            seq = "".join(aa for _, aa in residues)
            f.write(f">{cid}\n")
            for i in range(0, len(seq), 80):
                f.write(seq[i:i + 80] + "\n")


def write_sites(sites, outdir: Path):
    outdir.mkdir(parents=True, exist_ok=True)
    path = outdir / "nglyc_sites.tsv"

    with path.open("w") as f:
        f.write("site_id\tglycan_id\tchain\tasn_resid\tsequon\tx_resid\tst_resid\tasn_chain_index\n")
        for s in sites:
            f.write(
                f"{s['site_id']}\t{s['glycan_id']}\t{s['chain']}\t{s['asn_resid']}\t{s['sequon']}\t"
                f"{s['x_resid']}\t{s['st_resid']}\t{s['asn_chain_index']}\n"
            )

    return path


def write_esm(chains, sites, outdir: Path):
    outdir.mkdir(parents=True, exist_ok=True)
    fasta = outdir / "protein.fasta"
    write_fasta(chains, fasta)

    note = outdir / "README.txt"
    note.write_text(
        "ESMFold/ESM2 is protein-sequence-only and cannot explicitly model glycans.\n"
        "Use this FASTA for a protein-only backbone/reference model.\n"
        "Detected glycosylation sites are in ../sites/nglyc_sites.tsv.\n"
    )
    return [fasta, note]


def write_af3(chains, sites, glycan, tree, outdir: Path, name: str, preset_name: str | None = None):
    outdir.mkdir(parents=True, exist_ok=True)
    sequences = []
    bonds = []

    for cid, residues in chains.items():
        sequences.append({
            "protein": {
                "id": cid,
                "sequence": "".join(aa for _, aa in residues),
            }
        })

    protein_atom, glycan_atom = protein_glycan_link_atoms(glycan)

    for site in sites:
        gly_id = site["glycan_id"]
        sequences.append({
            "ligand": {
                "id": gly_id,
                "ccdCodes": list(glycan),
            }
        })
        bonds.append([
            [site["chain"], site["asn_chain_index"], protein_atom],
            [gly_id, 1, glycan_atom],
        ])

        for parent, child, atom1, atom2 in glycan_bonds(glycan, tree, preset_name):
            a = glycan[parent]
            b = glycan[child]

            bonds.append([
                [gly_id, parent + 1, atom1],
                [gly_id, child + 1, atom2],
            ])

    payload = {
        "name": f"{name}_all_nglyc_sites",
        "modelSeeds": [1],
        "sequences": sequences,
        "bondedAtomPairs": bonds,
        "dialect": "alphafold3",
        "version": 2,
    }

    path = outdir / f"{name}_all_nglyc_sites_af3.json"
    path.write_text(json.dumps(payload, indent=2))

    manifest = write_linked_glycan_manifest(sites, glycan, tree, outdir, name, "af3", preset_name)
    sanitizer = write_glycan_sanitizer(outdir)
    readme = outdir / "README_glycan_handoff.txt"
    readme.write_text(
        "AF3 glycan handoff\n"
        "==================\n\n"
        "The AF3 JSON uses CCD monosaccharides plus bondedAtomPairs because that is the "
        "native predictor representation. Some predictor outputs retain standalone "
        "CCD leaving-group atoms such as anomeric O1 even when C1 is covalently linked.\n\n"
        "Before using an AF3 PDB with GLYCAM/OpenMM, run:\n\n"
        f"  python {GLYCAN_SANITIZER_NAME} <af3_output.pdb> -m {GLYCAN_LINK_MANIFEST_NAME} -o <cleaned.pdb>\n\n"
        "The sanitizer removes O1 from C1-linked sugars and rewrites CONECT records "
        "from the intended covalent bonds.\n"
    )

    return [path, manifest, sanitizer, readme]


def write_boltz2(chains, sites, glycan, tree, outdir: Path, name: str, preset_name: str | None = None):
    if yaml is None:
        raise RuntimeError("PyYAML missing. Install with: pip install pyyaml")

    outdir.mkdir(parents=True, exist_ok=True)
    sequences = []
    constraints = []
    site_metadata = []

    for cid, residues in chains.items():
        sequences.append({
            "protein": {
                "id": cid,
                "sequence": "".join(aa for _, aa in residues),
            }
        })

    protein_atom, glycan_atom = protein_glycan_link_atoms(glycan)

    for site in sites:
        gly_id = site["glycan_id"]
        sequences.append({
            "ligand": {
                "id": gly_id,
                "ccd": list(glycan),
            }
        })
        constraints.append({
            "bond": {
                "atom1": [site["chain"], site["asn_chain_index"], protein_atom],
                "atom2": [gly_id, 1, glycan_atom],
            }
        })
        site_metadata.append({
            "site_id": site["site_id"],
            "chain": site["chain"],
            "asn_resid": site["asn_resid"],
            "sequon": site["sequon"],
            "glycan_id": gly_id,
        })

        for parent, child, atom1, atom2 in glycan_bonds(glycan, tree, preset_name):
            a = glycan[parent]
            b = glycan[child]

            constraints.append({
                "bond": {
                    "atom1": [gly_id, parent + 1, atom1],
                    "atom2": [gly_id, child + 1, atom2],
                }
            })

    payload = {
        "version": 1,
        "sequences": sequences,
        "constraints": constraints,
        "metadata": {
            "name": f"{name}_all_nglyc_sites",
            "glycan": "-".join(glycan),
            "tree": tree,
            "preset": preset_name,
            "sites": site_metadata,
            "note": "Check CCD atom names/linkages before production.",
        },
    }

    path = outdir / f"{name}_all_nglyc_sites_boltz2.yaml"
    with path.open("w") as f:
        yaml.safe_dump(payload, f, sort_keys=False)

    manifest = write_linked_glycan_manifest(sites, glycan, tree, outdir, name, "boltz2", preset_name)
    sanitizer = write_glycan_sanitizer(outdir)
    readme = outdir / "README_glycan_handoff.txt"
    readme.write_text(
        "Boltz2 glycan handoff\n"
        "=====================\n\n"
        "The Boltz2 YAML uses CCD monosaccharides plus covalent bond constraints because "
        "that is the native predictor representation. Boltz may emit standalone CCD "
        "sugar atom sets in PDB output, including anomeric O1 atoms that should be "
        "absent once C1 is covalently linked.\n\n"
        "Before using a Boltz PDB with GLYCAM/OpenMM, run:\n\n"
        f"  python {GLYCAN_SANITIZER_NAME} <boltz_output.pdb> -m {GLYCAN_LINK_MANIFEST_NAME} -o <cleaned.pdb>\n\n"
        "The sanitizer removes O1 from C1-linked sugars and rewrites CONECT records "
        "from the intended covalent bonds.\n"
    )

    return [path, manifest, sanitizer, readme]


def write_rfaa(
    chains,
    sites,
    glycan,
    tree,
    outdir: Path,
    name: str,
    glycan_sdf: Path,
    link_atom_index: int,
    link_chirality: str,
    preset_name: str | None = None,
):
    if yaml is None:
        raise RuntimeError("PyYAML missing. Install with: pip install pyyaml")
    if link_atom_index < 1:
        raise ValueError("--rfaa-link-atom-index must be 1 or greater.")

    outdir.mkdir(parents=True, exist_ok=True)
    protein_inputs = {}

    for cid, residues in chains.items():
        fasta = outdir / f"chain_{cid}.fasta"
        with fasta.open("w") as f:
            seq = "".join(aa for _, aa in residues)
            f.write(f">{cid}\n{seq}\n")

        protein_inputs[cid] = {"fasta_file": str(fasta.name)}

    sdf_input, sdf_note = rfaa_sdf_for_dir(glycan_sdf, outdir)
    sm_inputs = {}
    covalent_bonds = []
    site_metadata = []

    for site in sites:
        gly_id = site["glycan_id"]
        sm_inputs[gly_id] = {
            "input": sdf_input,
            "input_type": "sdf",
        }
        covalent_bonds.append(
            f"{site['chain']},{residue_int(site['asn_resid'])},ND2:"
            f"{gly_id},{link_atom_index}:{link_chirality},null"
        )
        site_metadata.append({
            "site_id": site["site_id"],
            "chain": site["chain"],
            "asn_resid": site["asn_resid"],
            "sequon": site["sequon"],
            "glycan_id": gly_id,
        })

    covale_inputs = ";".join(covalent_bonds)

    config = {
        "defaults": ["base"],
        "job_name": f"{name}_all_nglyc_sites",
        "protein_inputs": protein_inputs,
        "sm_inputs": sm_inputs,
        "covale_inputs": covale_inputs,
        "loader_params": {
            "MAXCYCLE": 10,
        },
        "metadata": {
            "glycan": "-".join(glycan),
            "tree": tree,
            "preset": preset_name,
            "sites": site_metadata,
            "rfaa_link_atom_index": link_atom_index,
            "rfaa_link_chirality": link_chirality,
            "note": (
                "RFAA covalent glycan input must use one merged SDF per glycan. "
                "The same configured SDF is reused for every detected site."
            ),
        },
    }

    cfg = outdir / "rfaa.yaml"
    with cfg.open("w") as f:
        yaml.safe_dump(config, f, sort_keys=False)

    run = outdir / "run_rfaa.sh"
    run.write_text(
        "#!/bin/bash\n"
        "set -euo pipefail\n\n"
        "cd \"$(dirname \"$0\")\"\n\n"
        "# Run from your RoseTTAFold-All-Atom environment.\n"
        "python -m rf2aa.run_inference -cd . --config-name rfaa\n"
    )
    run.chmod(0o755)

    readme = outdir / "README.txt"
    readme.write_text(
        "RFAA/RF3-style covalent glycan input\n"
        "====================================\n\n"
        "This directory contains one aggregate RFAA job for every detected N-glycosylation sequon.\n"
        "Each site gets its own small-molecule chain, all reusing the configured merged glycan SDF.\n\n"
        "RFAA cannot use the CCD glycan residue list directly for covalent glycans.\n"
        "Use a single merged SDF for the whole glycan, with leaving groups removed.\n\n"
        f"Intended glycan per site: {'-'.join(glycan)} ({tree})\n"
        f"Configured SDF: {sdf_input}\n"
        f"Configured covale_inputs: {covale_inputs if covale_inputs else '-'}\n\n"
        f"{sdf_note}\n\n"
        "Check the SDF atom order before production runs. The small-molecule atom "
        "index in covale_inputs is 1-indexed and must point to the ASN-linking atom.\n"
        "RFAA documents semicolon-separated multiple covalent bonds, but notes that "
        "multiple-bond input has not been heavily tested upstream.\n"
    )

    return [cfg, run, readme]


def main():
    args = parse_args()

    if args.list_glycans:
        print_supported_glycans()
        return

    if not args.file:
        raise SystemExit("error: -f/--file is required unless --list-glycans is used")

    infile = Path(args.file)
    model = normalize_model(args.model)
    preset_name = glycan_preset_name(args.glycans)
    glycan = parse_glycan(args.glycans)
    try:
        validate_glycan(glycan)
        if args.rfaa_link_atom_index < 1:
            raise ValueError("--rfaa-link-atom-index must be 1 or greater.")
    except ValueError as exc:
        raise SystemExit(f"error: {exc}") from exc

    root = Path(args.outdir) if args.outdir else Path(f"{infile.stem}_glyco_inputs")
    root.mkdir(parents=True, exist_ok=True)

    structure = parse_structure(infile)
    chains = extract_chains(structure)
    sites = find_nglyc_sites(chains)

    if not chains:
        raise RuntimeError("No protein chains found.")

    try:
        assign_glycan_chain_ids(sites, chains)
    except ValueError as exc:
        raise SystemExit(f"error: {exc}") from exc

    written = []
    written.append(write_sites(sites, root / "sites"))

    if model in {"esm", "all"}:
        written.extend(write_esm(chains, sites, root / "esm"))

    if model in {"af3", "all"}:
        written.extend(write_af3(chains, sites, glycan, args.tree, root / "af3", infile.stem, preset_name))

    if model in {"boltz2", "all"}:
        written.extend(write_boltz2(chains, sites, glycan, args.tree, root / "boltz2", infile.stem, preset_name))

    if model in {"rfaa", "all"}:
        written.extend(
            write_rfaa(
                chains,
                sites,
                glycan,
                args.tree,
                root / "rfaa",
                infile.stem,
                Path(args.rfaa_glycan_sdf),
                args.rfaa_link_atom_index,
                args.rfaa_link_chirality,
                preset_name,
            )
        )

    print(f"\nOutput folder: {root}")
    print(f"Detected {len(sites)} N-glycosylation sequon(s):")
    for s in sites:
        print(f"  {s['site_id']}: chain {s['chain']} ASN {s['asn_resid']} sequon {s['sequon']}")

    print("\nWrote:")
    for p in written:
        print(f"  {p}")

    print("\nNotes:")
    print("  ESM is protein-only.")
    if preset_name:
        print(f"  Glycan preset used: {preset_name} ({'-'.join(glycan)}).")
    print("  AF3/Boltz2 glycan atom names should be checked against CCD definitions.")
    print("  RFAA covalent glycans need one merged SDF; the CCD glycan string is only metadata there.")


if __name__ == "__main__":
    main()