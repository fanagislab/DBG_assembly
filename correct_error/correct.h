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

//��ƽ�����������Ϊ1%ʱ�����������ɵĵ�Ƶ����ռ�ܹ���80%������һ���ľ���Ч�ʴ�Լ��90%������������Ч��(���߳˻�)ԼΪ70%��
//������������н�����1һ������base�����ƽ����˵�����ٵģ�������һ������ܾ���40-50%�Ĵ���bases
int correct_one_base(string &read, uint8_t *freq,  int error_pos, int check_start, int check_end);

//�������ڵĸ�Ƶ����ͬʱȥ��С�ĵ�Ƶ����(length < HighFreqRegLenCutoff)
void get_high_freq_region(vector <FreqContReg> &freqContRegs, vector <FreqContReg> &highFreqContRegs, int &num_in_highFreqRegs, int &num_out_highFreqRegs);

//To correct the error base, modify on the read directly
void correct_one_read(string &read_head, string &read, uint8_t *freq, int &correct_one_base_score, int &correct_multi_base_score, int &is_deleted, int &trim_left_end_len, int &trim_right_end_len);

//��鷽��Ϊ����Ƶ��kmerβ(��)->��Ƶ��kmerβ(��)
int correct_multi_bases_rightward(string &read, uint8_t *freq, int check_start, int check_end, int &len_need_trim, int is_modify_trimmed_reads, int max_allowed_change, int &last_change_pos);

//get the kmer bit by trace back to the right forward
void get_kmer_rightward(uint64_t &kbit, uint8_t cur_base, vector <TreeNode> &node_vec, uint32_t pos, uint64_t start_point_bit);

//��鷽��Ϊ����Ƶ��kmerͷ(��)->��Ƶ��kmerͷ(��)
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
