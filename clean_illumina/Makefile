all:  clean_adapter  clean_lowqual

clean_adapter: clean_adapter.cpp  gzstream.cpp
	g++ -O3 -o $@ $^ -lz -lpthread

clean_lowqual: clean_lowqual.cpp  gzstream.cpp
	g++ -O3 -o $@ $^ -lz -lpthread

clean: 
	rm clean_adapter  clean_lowqual

