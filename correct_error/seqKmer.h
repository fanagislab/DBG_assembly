#ifndef __SEQ_KMER_H_
#define __SEQ_KMER_H_

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <set>
#include <cmath>
#include <zlib.h>
#include <inttypes.h>
#include <memory.h>
using namespace std;

extern char alphabet[128];
extern char bases[5];
extern char c_bases[5];
extern uint8_t bitAll[8];
extern uint64_t bitLeft[4];

//convert kmer-seq to kmer-bit��64bit�����װ32bp
//��Ҫ����ȷ��������ֻ����ACGT(��acgt)4�����
uint64_t seq2bit(string &kseq);

//convert kmer-bit to kmer-seq
//�˴��������kmer�ĳ���ֵ
string bit2seq(uint64_t kbit, int kmerSize);

//check whether a sequence contain non base characters, such as "N"
int check_seq (string &seq);

//get the reverse and complement sequence
void reverse_complement (string &in_str, string &out_str);

//used in routine: make_kmerFreq_2bit_table
//get the frequence value for a given index
int get_freq(uint8_t *freq, uint64_t idx);

//read file_list into a vector
void reading_file_list(string &file_list, vector<string> &files);

//get the reverse and complement kbit, independent of the sequence 
uint64_t get_rev_com_kbit(uint64_t kbit, uint8_t ksize);


#endif
