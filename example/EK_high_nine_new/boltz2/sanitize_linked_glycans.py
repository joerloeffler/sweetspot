#!/usr/bin/env python3
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
        line = raw.rstrip("\n")
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
    Path(args.output).write_text("\n".join(output) + "\n")

    print(f"Wrote {args.output}")
    if remove:
        print(f"Removed {len(remove)} linked-glycan leaving-group atom definition(s): " + ", ".join(f"{c}:{r}:{a}" for c, r, a in sorted(remove)))
    if missing_manifest_bonds:
        print("WARNING: Some manifest bonds could not be mapped onto the PDB:")
        for k1, k2 in missing_manifest_bonds:
            print(f"  {k1} -- {k2}")


if __name__ == "__main__":
    main()
