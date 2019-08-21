//
//this is a simple sequence alignment program, using dynamic programming algorith.
//only record one best path, and can only generate one best alignment result
//version 2

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <set>
#include "global_aligning.h"

using namespace std;



//当分数相等的时候，优先选择match/mismatch而不是gap, 这样得出来的gap会变得很碎，但保证比对上的碱基最多
void get_max_score (int sub_score, int gapi_score, int gapj_score, int &max_score, int &direction)
{	
	if (sub_score >= gapi_score && sub_score >= gapj_score)
	{	max_score = sub_score;
		direction = 0;
	}
	else if (gapi_score > sub_score && gapi_score >= gapj_score)
	{	max_score = gapi_score;
		direction = 1;
	}
	else
	{	max_score = gapj_score;
		direction = 2;
	}
	
}


//dealing global varibles,递归函数， direction stands for path 0:subs, 1:gapi; 2:gapj;
void trace_back (int *DPdirect, string &seq_i, string &seq_j, string &align_i, string &align_j, int &pos_i, int &pos_j, int &j_size)
{	
	int direction = DPdirect[pos_i * j_size + pos_j];

	if (direction == 0)
	{ 
		align_i.push_back( seq_i[pos_i - 1] );
		align_j.push_back( seq_j[pos_j - 1] );
		pos_i --;
		pos_j --;
	}
	else if (direction == 1)
	{ 
		align_i.push_back( '-' );
		align_j.push_back( seq_j[pos_j - 1] );
		pos_j --;
	}
	else if (direction == 2)
	{ 
		align_i.push_back( seq_i[pos_i - 1] );
		align_j.push_back( '-' ); 
		
		pos_i --;
	}

	if (pos_i > 0 || pos_j > 0)
	{	trace_back(DPdirect, seq_i, seq_j, align_i, align_j, pos_i, pos_j, j_size);
	}

}

//用来输出DP分值矩阵，以查错
void output_DPmatrix (string &seq_i, string &seq_j, int *DPscore, int *DPdirect)
{	
	int j_size = seq_j.size() + 1;

	cout << " \t" << " \t";
	for (int j=0; j<seq_j.size(); j++)
	{	cout << seq_j[j] << "\t";
	}
	cout << "\n";

	for (int i=0; i<=seq_i.size(); i++)
	{	if ( i == 0)
		{	cout << " \t";
		}
		else
		{	cout << seq_i[i-1] << "\t";
		}
		
		for (int j=0; j<=seq_j.size(); j++)
		{
			cout <<  DPscore[i * j_size + j] << "\t";			
		}
		cout << "\n";
	}
}

//外部调用函数
void global_aligning (string &seq_i, string &seq_j, string &align_i, string &align_j, int &final_score)
{	
	int gapPenalty = -5;
	
	char alphabet[128] =
	{
	 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	 4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 
	 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	 4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
	 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4
	};

	 //A(0), C(1), G(2), T(3), N(4)
	int scoreMatrix[5][5] = 
	{	{3, -5, -5, -5, -5},
		{-5, 3, -5, -5, -5},
		{-5, -5, 3, -5, -5},
		{-5, -5, -5, 3, -5},
		{-5, -5, -5, -5, 3}
	}; //the matrix score only used for assembly alignment


//	//A(0), C(1), G(2), T(3), N(4)
//	int scoreMatrix[5][5] = 
//	{	{2, -7, -5, -7, -5},
//		{-7, 2, -7, -5, -5},
//		{-5, -7, 2, -7, -5},
//		{-7, -5, -7, 2, -5},
//		{-5, -5, -5, -5, -5}
//	}; //this is the blast score matrix

	int *DPscore, *DPdirect;
	int j_size = seq_j.size() + 1;
	int array_size = (seq_i.size() + 1) * (seq_j.size() + 1);
	DPscore = new int[array_size];
	DPdirect = new int[array_size];
	
	//construct DP matrix
	DPscore[0] = 0; //stands for score
	DPdirect[0] = 0; //stands for path 0:subs, 1:gapi; 2:gapj;
	for (int j=1; j<=seq_j.size(); j++)
	{	
		DPscore[0 * j_size + j] = gapPenalty * j;  
		DPdirect[0 * j_size + j] = 1;
	}

	for (int i=1; i<=seq_i.size(); i++)
	{	
		DPscore[i * j_size + 0] = gapPenalty * i;
		DPdirect[i * j_size + 0] = 2;
	}
	
	//i: y-axis; j: x-axis;
	for (int i=1; i<=seq_i.size(); i++)
	{	for (int j=1; j<=seq_j.size(); j++)
		{	int subsScore = DPscore[(i-1) * j_size + (j-1)] + scoreMatrix[ alphabet[ seq_i[i-1] ] ][ alphabet[ seq_j[j-1] ] ];
			int gapiScore = DPscore[i * j_size + (j-1)] + gapPenalty;
			int gapjScore = DPscore[(i-1) * j_size + j] + gapPenalty;
			get_max_score(subsScore,gapiScore,gapjScore,DPscore[i * j_size + j],DPdirect[i * j_size + j]);
		}
	}
	
	//cerr << "after assign all DP values\n";
	
	int pos_i = seq_i.size(); //初值为seqi序列长度
	int pos_j = seq_j.size(); //初值为seqj序列长度
	
	final_score = DPscore[pos_i * j_size + pos_j];

	trace_back(DPdirect, seq_i, seq_j, align_i, align_j, pos_i, pos_j, j_size);
	reverse( align_i.begin(), align_i.end() );
	reverse( align_j.begin(), align_j.end() );
	//cerr << "after trace back\n";

	//output the DP score and direction matrix
	//output_DPmatrix (seq_i, seq_j, DPscore, DPdirect);

	delete []DPscore;
	delete []DPdirect;

}


