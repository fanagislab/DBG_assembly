#ifndef __GLOBALALIGN_H_
#define __GLOBALALIGN_H_

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

//当分数相等的时候，优先选择match/mismatch而不是gap, 这样得出来的gap会变得很碎，但保证比对上的碱基最多
void get_max_score (int sub_score, int gapi_score, int gapj_score, int &max_score, int &direction);

//dealing global varibles,递归函数， direction stands for path 0:subs, 1:gapi; 2:gapj;
void trace_back (int *DPdirect, string &seq_i, string &seq_j, string &align_i, string &align_j, int &pos_i, int &pos_j, int &j_size);

//用来输出DP分值矩阵，以查错
void output_DPmatrix (string &seq_i, string &seq_j, int *DPscore, int *DPdirect);

//外部调用函数
void global_aligning (string &seq_i, string &seq_j, string &align_i, string &align_j, int &final_score);


#endif
