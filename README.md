# interdomain_hamiltonian
# Hamiltonian calculation for inter-domain pairs 
## Fastahamiltonian.m
Calculates the Hamiltonian using a subset of the couplings.  That is, it only uses the couplings across two segments of the sequence and ignores the couplings within those segments. 
## Generalhamiltonian_corrected.m
Calculates the Hamiltonian for any number of sequences.  Input is an alignment of sequences in number format, the couplings and local fields. 
## completeDCA.m
Uses mean field DCA to estimate couplings and local field parameters for an input MSA in fasta format. Also outputs the alignment in number format. 
## concat_scramble_unconcat_pfamdatabase.m
Mixes and matches user-specified segments of an alignment and outputs a new fasta file. 
