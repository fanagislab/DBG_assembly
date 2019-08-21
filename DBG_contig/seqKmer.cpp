//author: Fan Wei, email: fanweiagis@126.com, date: 2015-11-4
//collect useful C/C++ subroutines related with Kmers
//ע�⣺���ļ��ڵĺ�����ȫ�ֱ��������Թ��ⲿ�������

#include "seqKmer.h"

//�ɣ��ãǣԵ�ASCII�뵽�������������Զ������Сд
//256����ĸ��alphabet����,��8bit��char�ʹ洢,A=a=N=n=0,C=c=1,G=g=2,T=t=3,��������ĸ��Ϊ4
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

//�ɣ������� 4�����ãǣ� N
char bases[5] ={
		'A', 'C', 'G', 'T', 'N'
};

//�ɣ������� 4�����ãǣ� N�Ļ������
char c_bases[5] ={
		'T', 'G', 'C', 'A', 'N'
};


//convert kmer-seq to kmer-bit��64bit�����װ32bp
//��Ҫ����ȷ��������ֻ����ACGT(��acgt)4�����
uint64_t seq2bit(string &kseq)
{
	uint64_t kbit=0;
	for(int i=0; i<kseq.size(); i++) {
			kbit=(kbit<<2)|alphabet[kseq[i]];
	}
	return kbit;
}

//convert kmer-bit to kmer-seq
//�˴��������kmer�ĳ���ֵ
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

//get the complement sequence
void complement_sequence (string &str)
{	
	for (int64_t i=0; i<str.size(); i++)
	{	
		str[i] = c_bases[alphabet[str[i]]];
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


//read file_list into a vector
void reading_file_list(string &file_list, vector<string> &files)
{	
	ifstream infile ( file_list.c_str() );
	if ( ! infile )
	{	cerr << "fail to open input file" << file_list << endl;
	}

	string file_name;
	while ( getline( infile, file_name, '\n' ) )
	{	if (file_name.size())
		{	files.push_back(file_name);
		}
	}
}

//display a number in bit format
void display_num_in_bits(uint64_t num, int len)
{	
	cout << num << "\t" << len << "\t";
	len -= 1;
	for (int i=0; i<=len; i++)
	{	int bit = (num >> (len-i)) & 1;
		cerr << bit;
	}
	cerr << endl;
}


//chengfang calculation for small integers
uint64_t pow_integer(int base, int exponent)
{	uint64_t result = 1;
	for (int i = 1; i<=exponent; i++)
	{	result *= base;
	}
	return result;
}
