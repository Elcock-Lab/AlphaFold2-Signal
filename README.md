# AlphaFold2-Signal
Simple Fortran code that analyzes atomic contacts between a predicted signal and the rest of the protein in an AlphaFold2-predicted protein structure

# Compilation
The code is very simple so should be easy to compile with any Fortran compiler such as gfortran. e.g.:

`gfortran analyze_alphafold2_make_contact_map_signal_for_release.f -o analyze_alphafold2_make_contact_map_signal_for_release.exe`

# Using the code

The code has five input arguments and produces one output file. Argument are:

(1) name of the AlphaFold2 .pdb file - must contain full-length protein, i.e. including signal
(2) distance cutoff for defining atomic contacts (we used 4 Angstroms in our work)
(3) number of residues in the signal - the code assumes the signal is at the N-terminus
(4) pLDDT threshold for including residues in the atomic contact calculation
(5) number of residues either side of the cleavage site to exclude from atomic contact calculations (we used 4 in our work)

With all the above in mind, the code could be run thus:

`./analyze_alphafold2_make_contact_map_signal_for_release.exe your-alphafold2-pdb-here 4.0 25 90.0 4`

The code generates a single output file that is self-explanatory - this is:

`signal_contacts.txt`
