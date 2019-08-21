all: correct_error correct_error_parallel correct_error_reads simulate_lowfreq_kmer 

correct_error: seqKmer.cpp main.cpp correct.cpp gzstream.cpp
	g++  -O3  -o $@  $^ -lz

correct_error_parallel: seqKmer.cpp main_pthread_read.cpp correct.cpp gzstream.cpp
	g++ -O3  -o $@  $^ -lz -lpthread 

correct_error_reads: seqKmer.cpp main_parallel_senior.cpp correct.cpp gzstream.cpp
	g++ -O3  -o $@  $^ -lz -lpthread 


simulate_lowfreq_kmer: simulate_lowfreq_kmer.cpp seqKmer.cpp gzstream.cpp
	g++  -O3  -o $@  $^ -lz
	
clean:
	rm  correct_error correct_error_parallel simulate_lowfreq_kmer  correct_error_reads

