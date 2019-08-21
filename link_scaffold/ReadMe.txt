In this package, there are four programs: map_pair and link_scaffold are co-used for pair-end and mate-pair mapping;
while map_reads and link_contig are co-used for single reads mapping.


1. map_pair, this function maps pairs of reads onto the contig sequences.

   This program use a seed-and-extension globle alignment method, and caluclate the identity between reads and contigs

2. link_scaffold, this function convert read PEs into contig relations, links the contigs into scaffolds.

  (1) If two contigs are neighbor, and they are in the same DNA strand, a pair of PE read will mapped in 
   the F and R way, that is, one read mapped forwardly on one contig, and the other read mapped reversely
   on the other contig. We use this character to decide the strand relations between two neighbor contigs.

  (2) Each contig has two strand forms, we stored both of them in the memory, and each contig node only
   has one link, the 3'-direction. This can make the strand problem quite easy, and do not increase the
   memory usage too much.

  (3) This program require more than a specified number of read pairs to support two neighbor contigs,
   and build a graph using neighboring-link-technology, take contig as nodes, and the 3'-link as arcs.

  (4) The scaffold sequences were read out from a linear contig, that is, with only one ingoing arc and one
   outgoing arc, and scaffolding stops at breaking or branching nodes.


3. map_reads, this function maps single reads onto the contig sequences, used similar alignment method with map_pair


4. link_contig, this function links the contigs by read-ends that mapped to two different contigs, and fill the gaps with 
   the conensus sequence inferred from all the reads mapped to the gap. The link strategy is similar to link_scaffold,
   the great advantage here is that all the inside gaps can be accurately filled.

