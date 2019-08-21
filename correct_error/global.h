#ifndef __GLOBAL_H_
#define __GLOBAL_H_

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <set>
#include <cmath>
#include<pthread.h>
#include <inttypes.h>

using namespace std;

extern int KmerSize;
extern int Low_freq_cutoff;
extern int HighFreqRegLenCutoff;
extern int Max_change_in_one_read;
extern int Max_node_in_BB_tree;
extern int Min_trimmed_read_len;
extern uint8_t *KmerFreq;
extern int Further_trim_len;
extern int Join_read1_read2_file;


typedef struct {
	int start;
	int end;
	int status;
} FreqContReg;


typedef struct {
        uint32_t pointer:26,base:2,change:3,same:1;
} TreeNode;


#endif
