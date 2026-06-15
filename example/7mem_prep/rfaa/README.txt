RFAA/RF3-style covalent glycan input
====================================

This directory contains one aggregate RFAA job for every detected N-glycosylation sequon.
Each site gets its own small-molecule chain, all reusing the configured merged glycan SDF.

RFAA cannot use the CCD glycan residue list directly for covalent glycans.
Use a single merged SDF for the whole glycan, with leaving groups removed.

Intended glycan per site: NAG-NAG-MAN (linear)
Configured SDF: glycan.sdf
Configured covale_inputs: A,154,ND2:B,1:CW,null;C,20,ND2:D,1:CW,null;C,21,ND2:E,1:CW,null;C,33,ND2:F,1:CW,null;C,94,ND2:G,1:CW,null;C,278,ND2:H,1:CW,null;C,289,ND2:I,1:CW,null

RFAA glycan SDF was not found at glycan.sdf. Place a merged glycan SDF named glycan.sdf in this directory or edit sm_inputs.B.input before running.

Check the SDF atom order before production runs. The small-molecule atom index in covale_inputs is 1-indexed and must point to the ASN-linking atom.
RFAA documents semicolon-separated multiple covalent bonds, but notes that multiple-bond input has not been heavily tested upstream.
