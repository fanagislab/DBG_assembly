//author: Fan Wei, email: fanweiagis@126.com, date: 2015-11-4
//collect useful C/C++ subroutines related with Kmers
//注意：本文件内的函数和全局变量均可以供外部程序调用

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
#include <inttypes.h>

using namespace std;

//由ＡＣＧＴ到ASCII码到０１２３，能自动处理大小写
extern char alphabet[128];

//由０１２３ 4到ＡＣＧＴ N
extern char bases[5];

//由０１２３ 4到ＡＣＧＴ N的互补碱基
extern char c_bases[5];

//convert kmer-seq to kmer-bit，64bit最多能装32bp
//需要事先确定序列中只含有ACGT(或acgt)4个碱基
uint64_t seq2bit(string &kseq);

//convert kmer-bit to kmer-seq
//此处必须给定kmer的长度值
string bit2seq(uint64_t kbit, int kmerSize);

//check whether a sequence contain non base characters, such as "N"
int check_seq (string &seq);

//get the complement base
inline char complement_base(char base)
{	return c_bases[alphabet[base]];
}

//get the reverse and complement sequence
void reverse_complement (string &in_str, string &out_str);

//read file_list into a vector
void reading_file_list(string &file_list, vector<string> &files);

//get the reverse and complement kbit, independent of the sequence 
uint64_t get_rev_com_kbit(uint64_t kbit, uint8_t ksize);

//display a number in bit format
void display_num_in_bits(uint64_t num, int len);

//chengfang calculation for small integers
uint64_t pow_integer(int base, int exponent);

//get the complement sequence
void complement_sequence (string &str);


#endif
