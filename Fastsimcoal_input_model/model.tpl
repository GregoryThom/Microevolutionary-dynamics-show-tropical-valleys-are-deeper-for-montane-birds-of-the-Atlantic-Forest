//Parameters for the coalescence simulation program : fastsimcoal.exe
3 samples to simulate :
//Population effective sizes (number of genes) pop0 = M, pop1 = S, pop2 = A
NPOP0
NPOP1
//Samples sizes and samples age
10
14
//Growth rates	: negative growth implies population expansion
0
0
//Number of migration matrices : 0 implies no migration between demes
2
//Migration matrix 0
0	M10
M01	0
//Migration matrix 1
0	0
0	0
//historical event: time, source, sink, migrants, new deme size, new growth rate, migration matrix index
3 historical event
tch1 0 0 1 SC1 0 1
tch2 1 1 1 SC2 0 1
TDIV1 0 1 1 RESIZE1 0 1
//Number of independent loci [chromosome] 
1 0
//Per chromosome: Number of contiguous linkage Block: a block is a set of contiguous loci
1
//per Block:data type, number of loci, per generation recombination and mutation rates and optional parameters
FREQ 1 0.0000 2.5e-9 OUTEXP



























