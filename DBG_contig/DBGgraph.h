//This file include variables and routines to build the kmer graph which is stored in the KmerSet
//Author: Wei Fan, fanweiagis@126.com
//Date: 2015-11-4


#ifndef __DBGRAPH_H_
#define __DBGRAPH_H_


#include <iostream>
#include <cmath>
#include <ctime>
#include <inttypes.h>
#include <string.h>
#include <stdlib.h>
#include <inttypes.h>
#include <pthread.h>
#include "gzstream.h"
#include "seqKmer.h"
#include "kmerSet.h"

using namespace std;


extern int KmerSize;             //size of kmer
extern int maxReadLen;           //max read length allowed, accosiate with memory allocation, reads lengther longer than this value will be trimmed off
extern int KmerNumInRead;        //kmer number in each reads with max allowed read length, equal to: maxReadLen - KmerSize + 1;  
extern int Input_file_format;    //input file format: fq, 1; fa,2;
extern string Output_prefix;      //prefix of all the output files
extern int threadNum;            //thread number used for parallel computation
extern KmerSet *kset;            //the major data structure that store the kmer deburijn graph
extern double initHashSize;      //1 G array size in the hash, as the the unit size is 16 byte, so the hash takes 16G memory
extern uint64_t maxDoubleHashTimes;  //2^10 times larger, 1T array size, 16T memory
extern uint64_t doubleHashTimes;     //current doubling times of the kmerset
extern float hashLoadFactor;         //loading factor of the hash
extern int BufferNum;            //number of reads stored each time in the memeory buffer 
extern string *RawReads;         //buffer memory to store the reads sequence
extern uint64_t *StoreKmer;      //buffer memory to store the chopped kmers
extern uint8_t *StoreLeftBase;   //buffer memory to store the left base of the kmer
extern uint8_t *StoreRightBase;  //buffer memory to store the right base of the kmer
extern uint8_t *Signal;          //buffer memory to store the processing status: 0, empty; 1, filled; 2, end of file;
extern uint64_t Kmer_total_num;     //the total number of processed kmers from input files
extern uint64_t Total_reads_num;    //the total number of processed reads from input files
extern uint64_t KmerHeadMaskVal;    //used to get the kmer bit value on the forward strand from the previous one 
extern uint64_t KmerRCOrVal[4];     //used to get the kmer bit value on the reverse strand from the previous one 
extern KmerNode *PolyA;             //store the information for polyA and polyT kmers

extern clock_t time_start;
extern clock_t time_end;

//this is the thread routine to convert reads to kmers and get the leftbase and rightbase parallely
//each reads is assigned to a specific thread, i.e. different threads parse different reads
void *thread_parseBlock(void* threadId_p);


//this is the thread routine to convert reads to kmers parallely
//each specific kmer species is assigned to a specific thread, i.e. different threads parse different kmer species
void *thread_updatekmers(void* threadId_p);


//parse one reads file in parallel, controlling the children threads
void parse_one_reads_file(string &reads_file);


//initialization, load each input files and construct the kmer de bruijn graph
void build_debruijn_graph( vector<string> &reads_files);



#endif



