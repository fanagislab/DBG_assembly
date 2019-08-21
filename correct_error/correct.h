#ifndef __CORRECT_H_
#define __CORRECT_H_

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <set>
#include <cmath>
#include <inttypes.h>
#include "seqKmer.h"
#include "global.h"

using namespace std;


//get continuous low or high frequence kmers on reads, this is a cluster function
void get_cont_kmerfreq_region(string &read, uint8_t *freq, vector <FreqContReg> &continuous_regions);

//当平均测序错误率为1%时，这类错误造成的低频区域占总共的80%，而这一步的纠错效率大约有90%，因此总体纠错效率(两者乘积)约为70%，
//由于这类错误中仅含有1一个错误base，相对平均来说算最少的，紧靠这一步大搞能纠正40-50%的错误bases
int correct_one_base(string &read, uint8_t *freq,  int error_pos, int check_start, int check_end);

//连接相邻的高频区域，同时去除小的低频区域(length < HighFreqRegLenCutoff)
void get_high_freq_region(vector <FreqContReg> &freqContRegs, vector <FreqContReg> &highFreqContRegs, int &num_in_highFreqRegs, int &num_out_highFreqRegs);

//To correct the error base, modify on the read directly
void correct_one_read(string &read_head, string &read, uint8_t *freq, int &correct_one_base_score, int &correct_multi_base_score, int &is_deleted, int &trim_left_end_len, int &trim_right_end_len);

//检查方向为：高频区kmer尾(左)->低频区kmer尾(右)
int correct_multi_bases_rightward(string &read, uint8_t *freq, int check_start, int check_end, int &len_need_trim, int is_modify_trimmed_reads, int max_allowed_change, int &last_change_pos);

//get the kmer bit by trace back to the right forward
void get_kmer_rightward(uint64_t &kbit, uint8_t cur_base, vector <TreeNode> &node_vec, uint32_t pos, uint64_t start_point_bit);

//检查方向为：高频区kmer头(右)->低频区kmer头(左)
int correct_multi_bases_leftward(string &read, uint8_t *freq, int check_start, int check_end, int &len_need_trim, int is_modify_trimmed_reads, int max_allowed_change, int &last_change_pos);

//get the kmer bit by trace back to the left forward
void get_kmer_leftward (uint64_t &kbit, uint8_t cur_base, vector <TreeNode> &node_vec, uint32_t pos, uint64_t start_point_bit);

//get the max highFreq region
void get_max_highFreq_region(vector <FreqContReg> &highFreqContRegs, FreqContReg &maxHighFreqReg, vector<int> &failCorrectReg_ids);

//parse one fasta format file
void parse_one_reads_fa_file(string &raw_reads_file);

//parse one fastaq format file
void parse_one_reads_fq_file(string &raw_reads_file);

//merge the read1 and read2 files into pair and single file
void merge_two_corr_files(string &read1_file, string &read2_file);

#endif
