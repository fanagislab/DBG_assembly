1. algorithm: similar to ruanjue's APE in many aspects:

(1) All the processes are performed on the kmer de bruijn graphs, this is 
similar to APE; but different from velvet and grape, which firstly simplified
the kmer graph into edge graph, and then mainly processed on the edge graph. 

(2) One remained question in short-read assembly method, is how to use all the available
information. Previous program such as velvet and grape mainly focused on the 
topology structure, such as tip and bubble, but do not fully use the coverage depth
information. This program considers both graph topology and coverage depth (stored
from 0-127), which will be benifit in distinguishing repeatitive and heterozygous sequences. 

(3) You are recommended to use the "-d" option, to remove the extrem low-frequency kmers
on the graph, which are mostly caused by sequencing errors. This action can also break down
many bubble structures into tips, which make it easy for error cleaning in later steps.

(4) By a combination of removing tips, bubbles, and low-coverage linear edges, most of the 
sequencing errors will be removed, and then the linear parts of the graph grows significantly. 
Note that all these graph-simplifying steps have coverage depth and length cutoff,
in order to avoid removing the genine sequences. 

(5) The bubble processing functions can be used to remove sequencing errors, they can also 
be used to merge heterozygous SNPs and Indels. The repeatitive and heterozgyous sequences will
both form bubble structures, however, repeatitive sequences will have higher coverage depth
than heterozygous sequences. You should be cautious to merge the bubbles, the spirit is that,
you can merge heterozygous sequences, but you can't merge repeatitive sequences.

(6) The final output is contig sequnces, which is defined as combination of linear kmer nodes
on the graph, i.e, continous unique Kmer regions in the genome. Contig will break down on 
the branching kmer nodes, which are formed by repeat structures.

(7) Along with the contig sequences, the coverage depth for each base is also output, the coverage
depth can be viewed as the qualitly evaluation of a base, given that high coverage depth reflects
high accuracy. The leftward and rightward ending kmers were also output as the contig sequences,
all the contig squences that starting or ending on a branching kmer node, will share the ending
kmers. You can use this information to know which contigs are neighbors.


2. Parameters and usage:

   build_contig  <reads_file_list>
   
   version: 0.3

   One input file: the address list of reads file£¬each reads files take a line, do not consider
   pair-end relation in this program. The reads file should be in fasta format, after sequencing
   error correction.
   
   Parmaters to deal with kmer size and frequency cutoff:
   -k <int>   set kmer size, default=31
   -d <int>   delete kmers with frequency no larger than, default=1
   
   Parameter to deal with removing tips:
   -t <int>   set the max allowed tip length, default=100
   -i <float>  set the max allowed tip depth, default=3
   
   Parameters to deal with removing bubbles:
   -b <int>   set the max allowed bubble length, default=100
   -u <int>   set the max allowed difference length between the two bubbles, default=10
   -l <float>  set the max allowed base difference rate between the two bubbles, default=0.1
   -e <float>  set the max allowed bubble depth, default=3
   
   Parameters to deal with removing low coverage edges
   -c <int>    set the max allowed length for low coverage linear fragments, default=100
   -v <float>  set the max allowed depth for low coverage linear fragments, default=3
   
   Other parametes:
   -o <str>   set the output prefix, default = output
   -h         get the help information


