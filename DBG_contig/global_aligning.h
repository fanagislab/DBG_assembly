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

//��������ȵ�ʱ������ѡ��match/mismatch������gap, �����ó�����gap���ú��飬����֤�ȶ��ϵļ�����
void get_max_score (int sub_score, int gapi_score, int gapj_score, int &max_score, int &direction);

//dealing global varibles,�ݹ麯���� direction stands for path 0:subs, 1:gapi; 2:gapj;
void trace_back (int *DPdirect, string &seq_i, string &seq_j, string &align_i, string &align_j, int &pos_i, int &pos_j, int &j_size);

//�������DP��ֵ�����Բ��
void output_DPmatrix (string &seq_i, string &seq_j, int *DPscore, int *DPdirect);

//�ⲿ���ú���
void global_aligning (string &seq_i, string &seq_j, string &align_i, string &align_j, int &final_score);


#endif
