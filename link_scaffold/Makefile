all: map_reads  map_pair link_scaffold link_contig  link_supertig

map_reads: kmerSet.cpp  map_func.cpp  seqKmer.cpp  gzstream.cpp map_reads.cpp
	g++ -O2  -o $@  $^  -lz -lpthread


map_pair: kmerSet.cpp  map_func.cpp  seqKmer.cpp  gzstream.cpp map_pair.cpp
	g++ -O2  -o $@  $^  -lz -lpthread


link_scaffold: link_func.cpp  link_scaffold.cpp  seqKmer.cpp gzstream.cpp
	g++ -O2  -o $@  $^ -lz

link_contig: link_func.cpp  link_contig.cpp  seqKmer.cpp gzstream.cpp
	g++ -O2  -o $@  $^ -lz

link_supertig: link_func.cpp  link_supertig.cpp  seqKmer.cpp gzstream.cpp
	g++ -O2  -o $@  $^ -lz


clean:
	rm map_reads  map_pair link_scaffold link_scafftig link_supertig

