//author: Fan Wei, email: fanw@genomics.org.cn, date: 2010-2-1
//collect useful C/C++ subroutines related with Kmers
//注意：本文件内的函数和全局变量均可以供外部程序调用
//Date: 2011-7-18
//version: 2.0

#include "seqKmer.h"

//由ＡＣＧＴ到ASCII码到０１２３，能自动处理大小写
//256个字母表alphabet数组,用8bit的char型存储,A=a=N=n=0,C=c=1,G=g=2,T=t=3,其他的字母都为4
char alphabet[128] =
{
 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
 4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 0, 4, 
 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
 4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 0, 4,
 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4
};

//由０１２３ 4到ＡＣＧＴ N
char bases[5] ={
		'A', 'C', 'G', 'T', 'N'
};

//由０１２３ 4到ＡＣＧＴ N
char c_bases[5] ={
		'T', 'G', 'C', 'A', 'N'
};

//used for bit-operation
uint8_t bitAll[8] = {128,64,32,16,8,4,2,1};

uint64_t bitLeft[4] = {
	0x0000000000000000,
	0x4000000000000000,
	0x8000000000000000,
	0xc000000000000000
};


//convert kmer-seq to kmer-bit，64bit最多能装32bp
//需要事先确定序列中只含有ACGT(或acgt)4个碱基
uint64_t seq2bit(string &kseq)
{
	uint64_t kbit=0;
	for(int i=0; i<kseq.size(); i++) {
			kbit=(kbit<<2)|alphabet[kseq[i]];
	}
	return kbit;
}

//convert kmer-bit to kmer-seq
//此处必须给定kmer的长度值
string bit2seq(uint64_t kbit, int kmerSize)
{
	string kseq;
	for(int i=0; i<kmerSize; i++) {
			kseq.push_back(bases[(kbit>>(kmerSize-1-i)*2)&0x3]);
	}
	return kseq;
}

//check whether a sequence contain non base characters, such as "N"
int check_seq (string &seq)
{       
	int is_good = 1;
	for (int i = 0; i < seq.size(); i++)
	{       if (alphabet[seq[i]] == 4)
			{   is_good = 0;
				break;
			}
	}
	return is_good;
}

//get the reverse and complement sequence
void reverse_complement (string &in_str, string &out_str)
{	
	for (int64_t i=in_str.size()-1; i>=0; i--)
	{	
		out_str.push_back(c_bases[alphabet[in_str[i]]]);
	}
}

//get the reverse and complement kbit, independent of the sequence 
uint64_t get_rev_com_kbit(uint64_t kbit, uint8_t ksize){
	kbit = ~kbit;
	kbit = ((kbit & 0x3333333333333333LLU)<<  2) | ((kbit & 0xCCCCCCCCCCCCCCCCLLU)>>  2);
	kbit = ((kbit & 0x0F0F0F0F0F0F0F0FLLU)<<  4) | ((kbit & 0xF0F0F0F0F0F0F0F0LLU)>>  4);
	kbit = ((kbit & 0x00FF00FF00FF00FFLLU)<<  8) | ((kbit & 0xFF00FF00FF00FF00LLU)>>  8);
	kbit = ((kbit & 0x0000FFFF0000FFFFLLU)<< 16) | ((kbit & 0xFFFF0000FFFF0000LLU)>> 16);
	kbit = ((kbit & 0x00000000FFFFFFFFLLU)<< 32) | ((kbit & 0xFFFFFFFF00000000LLU)>> 32);
	return kbit >> (64 - (ksize<<1));
}


//used in routine: make_kmerFreq_2bit_table
//get the frequence value for a given index
int get_freq(uint8_t *freq, uint64_t idx) 
{
	int value = ( freq[idx/8] >> (7-idx%8) ) & 0x1u;
	return value;
}


//read file_list into a vector
void reading_file_list(string &file_list, vector<string> &files)
{	
	ifstream infile ( file_list.c_str() );
	if ( ! infile )
	{	cerr << "fail to open input file" << file_list << endl;
	}

	string line_str;
	while ( getline( infile, line_str, '\n' ) )
	{	string file_name;
		for(int i=0; i<line_str.size();i++)
		{	if(line_str[i] != ' ' && line_str[i] != '\t' && line_str[i] != '\n')
			{	file_name.push_back(line_str[i]);
			}
		}
		if (file_name.size())
		{	files.push_back(file_name);
		}
	}
}


			
			
