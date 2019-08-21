all: debruijn_contig

debruijn_contig: main.cpp kmerSet.cpp seqKmer.cpp  DBGgraph.cpp gzstream.cpp contig.cpp global_aligning.cpp
	g++ -O2  -o $@  $^ -lz -lpthread

clean:
	rm debruijn_contig

