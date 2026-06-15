AF3 glycan handoff
==================

The AF3 JSON uses CCD monosaccharides plus bondedAtomPairs because that is the native predictor representation. Some predictor outputs retain standalone CCD leaving-group atoms such as anomeric O1 even when C1 is covalently linked.

Before using an AF3 PDB with GLYCAM/OpenMM, run:

  python sanitize_linked_glycans.py <af3_output.pdb> -m linked_glycans.json -o <cleaned.pdb>

The sanitizer removes O1 from C1-linked sugars and rewrites CONECT records from the intended covalent bonds.
