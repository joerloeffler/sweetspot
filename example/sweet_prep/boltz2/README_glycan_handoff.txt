Boltz2 glycan handoff
=====================

The Boltz2 YAML uses CCD monosaccharides plus covalent bond constraints because that is the native predictor representation. Boltz may emit standalone CCD sugar atom sets in PDB output, including anomeric O1 atoms that should be absent once C1 is covalently linked.

Before using a Boltz PDB with GLYCAM/OpenMM, run:

  python sanitize_linked_glycans.py <boltz_output.pdb> -m linked_glycans.json -o <cleaned.pdb>

The sanitizer removes O1 from C1-linked sugars and rewrites CONECT records from the intended covalent bonds.
