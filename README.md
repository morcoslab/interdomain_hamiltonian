# interdomain_hamiltonian
# Hamiltonian calculation for inter-domain pairs 
## completeDCA.m
MATLAB code, used to extract coevolutionary parameters, eij and hi, from homologous multiple sequence alignments (MSAs) of interested protein family. Also outputs the DI values and the alignment in number format. 
## Generalhamiltonian_corrected.m
MATLAB code. Calculates the Hamiltonian for any number of sequences.  Input is an alignment of sequences in number format, the couplings and local fields. 
## Fastahamiltonian.m
MATLAB code. Calculates the Hamiltonian using a subset of the couplings.  That is, it only uses the couplings across two segments of the sequence and ignores the couplings within those segments. 
## concat_scramble_unconcat_pfamdatabase.m
MATLAB code. Mixes and matches user-specified segments of an alignment and outputs a new fasta file. 
