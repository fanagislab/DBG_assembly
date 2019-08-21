#ifndef __MAP_FUNC_H_
#define __MAP_FUNC_H_

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <set>
#include <cmath>
#include <ctime>
#include "seqKmer.h"
#include "kmerSet.h"
#include <inttypes.h>
#include "gzstream.h"




using namespace std;

extern int KmerSize;               //kmer size
extern double MinMapIdentity;      //minimum mapping identity allowed to determine a true mapping
extern int SeedKmerNum;            //number of kmers in the seed-alignment
extern int MinReadLen;             //minimum read length, discard shorter reads
extern int MinCtgLen;               //minimum contig length, ignore short contigs
extern int Input_file_format;      //input file format: fq, 1; fa,2;
extern string Output_prefix;       //input file format: 1,fq; 2,fa

extern uint64_t KmerHeadMaskVal;    //used to get the kmer bit value on the forward strand from the previous one 
extern uint64_t KmerRCOrVal[4];     //used to get the kmer bit value on the reverse strand from the previous one 


//split a string into a vector
void split(string &strLine, vector<string>& tokens, const char* delim);

//read file_list into a vector
void reading_lib_file(string &para_file, vector<string> &read_files);

//read multiple-fasta format contig file
void read_contig_file(string &contig_seq_file, vector<string> &contig_ids, vector <string> &contig_seqs);

//chop contig sequences into kmers, and add the kmers into the hash
void chop_contig_to_kmerset(KmerSet *kset, vector <string> &contig_seqs);


//find the exact-match alignment of a specified region between read and contig,  and determine the align direction
//search the seed alignment in a specified read region with given start and end positions
//use two kmers to locate the seed-alignment, use the first kmer to decide location and direction, use the second kmer to confirm
void get_align_seed(KmerSet *kset, string &read, int search_start, int search_end, int &contig_id_index, int &seed_contig_start, int &seed_contig_end, int &seed_read_start, int &seed_read_end, char &seed_direct);


//based on the seed alignment, extend the alignment regions until the read ends or contig ends, and calculate the base identity value
void extend_align_region(string &read, string &contig, float &identity, int &align_contig_start, int &align_contig_end, int &align_read_start, int &align_read_end, char &direct);

//split scaffold into contigs and stored in vector
void scaffold_to_contig(string &scaffold, vector <string> &contigs, vector <int> &starts);


#endif
