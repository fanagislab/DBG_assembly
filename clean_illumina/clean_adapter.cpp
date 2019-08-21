//This program is designed to clean the contamination in illumina sequencing data
//Trimming the adapters induced in illumina library construction, default is genomic DNA adapter 33-bp
//Find adapter contamination in reads by ungapped local dynamic programming alignment, only report the best hit with max score
//The output read length after trimmings should be >= Minimum_trimmed_read_len bp, else the whole read replaced by a empty string
//Process one end file at a time, do not consider pair-end relations
//This program can run in parallel using multiple threads
//Both the forward and reverse complementary strands are considered when doing alignment

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <set>
#include <cmath>
#include <inttypes.h>
#include <zlib.h>
#include <pthread.h>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include "gzstream.h"

using namespace std;

int Alignment_score_cutoff = 12;  //minimum alignment score 

int Minimum_trimmed_read_len = 75; //the required shortest length of resulting trimmed reads

int Also_use_adapterRC = 0;

vector<string> AdapterVec;
vector<string> AdapterVecId;
string illumina_adapter_file = "Both-adapter";

//global variables used for the threads
int threadNum = 3;
uint64_t bufferNum = 1000000;  //read bufferNum reads as one block into memory each time
int *gotReads;
uint64_t readNum;
string *StoreReads;
string *StoreHeads;
string *StoreQuals;
void* thread_trimReads(void* threadId_p);

uint64_t *num_clean_reads ;
uint64_t *num_clean_bases ;
uint64_t *num_adapter_trimmed_reads ;
uint64_t *num_adapter_trimmed_bases ;
uint64_t *num_Nmasked_short_reads ;
uint64_t *num_Nmasked_short_bases ;


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

//用2-bit表示一个base：A(0), C(1), G(2), T(3), N(4)
int scoreMatrix[5][5] = 
{	{1, -2, -2, -2, -2},
	{-2, 1, -2, -2, -2},
	{-2, -2, 1, -2, -2},
	{-2, -2, -2, 1, -2},
	{-2, -2, -2, -2, -2}
}; //matrix for DNA alignment with 95% identity


void usage() 
{	cout << "Description:\nclean_adapter identifies and trims adapter sequence in raw reads by ungapped local dynamic programming alignment. The program loads adapter sequences as aligning target from a default self-taken or specially user-defined multiple-fasta format file by parameter (-a), in which more records can be given at the same time.  In theory, this program can be used to filter any contaminant sequences besides adapters, however, this version is specifically written for trimming adapter, which considers adapters locating on the tail part of reads, moreover, the alignment search will stop when find the first qualified hit (>=minimum alignment score) for each reads. The alignment algorithm used ungapped dynamic programming local alignment, with a score matrix [match:1; mismatch:-2], and reports only the best hit, and the aligning score cutoff can be set by a parameter (-s). The program runs in a multiple thread mode (-t). The input file should be fastq or gzip-fastq format, and there are two resulting files: one is the clean reads file, and the other one is a statistics file.\n";
	cout << "\nUsage:\n  clean_adapter  <Input.fq.gz> <Output.clean.gz> <Output.clean.stat> " << endl;
	cout << "   Author: Fanwei, fanweiagis@126.com" << endl;
	cout << "   Version 1.1;" << endl;
	cout << "   -a <str>   contaminant sequence file: Both-adapter, R1-adapter, R2-adapter for default adapter files, otherwise for user-defined cotaminant file, default=" << illumina_adapter_file << endl;
	cout << "   -b <int>   use both strands of sequence for alignment, 0: no; 1:yes; default=" << Also_use_adapterRC << endl;
	cout << "   -s <int>   minimum alignment score, score matrix [match:1; mismatch:-2],default=" << Alignment_score_cutoff << endl;
	cout << "   -r <int>   minimum read length after trimming, default=" << Minimum_trimmed_read_len << endl;
	cout << "   -t <int>   number of threads to run, default=" << threadNum << endl;
	cout << "   -h         get help information" << endl;
	cout << "\nExample:\n  clean_adapter  -a Both-adapter -r 75 -s 12 sequencing_data_R1.fq.gz sequencing_data_R1.fq.nonAdapter.gz sequencing_data_R1.fq.nonAdapter.stat\n" << endl;

	exit(0);
}


//不容gap的动态规划局部比对程序，仅找到分值最大的一个比对结果作为输出
void local_ungapped_aligning (string &seq_i, string &seq_j, int &max_score, int &align_i_start, int &align_j_start, int &align_i_end, int &align_j_end)
{	
	int Given_read_length = seq_i.size();
	//cout << "Given_read_length: " << Given_read_length << endl;
	int adapter_size = seq_j.size();
	//cout << "adapter_size: " << adapter_size << endl;

	int array_size = (Given_read_length + 1) * (adapter_size + 1);
	int *DPscore = new int[array_size];
	//construct DP matrix
	DPscore[0] = 0; //stands for score
	for (int j=1; j<=adapter_size; j++)
	{	DPscore[j] = 0;  
	}
	
	for (int i=1; i<=Given_read_length; i++)
	{	DPscore[i * (adapter_size+1)] = 0;
	}
	
	int i_size = seq_i.size() + 1;
	int j_size = seq_j.size() + 1;
	max_score = 0;
	
	//i: y-axis; j: x-axis;
	int pre_i;
	int pre_j;
	for (int i=1; i<i_size; i++)
	{	for (int j=1; j<j_size; j++)
		{	pre_i = i - 1;
			pre_j = j - 1;
			int this_score = DPscore[pre_i * j_size + pre_j] + scoreMatrix[ alphabet[ seq_i[pre_i] ] ][ alphabet[ seq_j[pre_j] ] ];
			if (this_score < 0)
			{	this_score = 0;
			}
			DPscore[i * j_size + j] = this_score; //(i+1), (j+1)代表1-based order
			if (this_score > max_score)
			{	max_score = this_score;
				align_i_end = i;
				align_j_end = j;
			}
		}
	}
	
	//trace back by looping instead of recursion
	int pos_i = align_i_end;
	int pos_j = align_j_end;
	int score = max_score;
	while (score > 0)
	{	
		//align_i.push_back( seq_i[pos_i - 1] );
		//align_j.push_back( seq_j[pos_j - 1] );
		pos_i --;
		pos_j --;
		score = DPscore[pos_i * j_size + pos_j];
	}
	align_i_start = pos_i + 1;
	align_j_start = pos_j + 1;
	//reverse( align_i.begin(), align_i.end() );
	//reverse( align_j.begin(), align_j.end() );
	//cerr << "after trace back\n";
	
	delete []DPscore;

}

//get the reverse and complement sequence
void reverse_complement (string &in_str, string &out_str)
{	//由０１２３ 4到ＡＣＧＴ N
	char c_bases[5] ={
			'T', 'G', 'C', 'A', 'N'
	};

	for (int64_t i=in_str.size()-1; i>=0; i--)
	{	
		out_str.push_back(c_bases[alphabet[in_str[i]]]);
	}
}

//this is the thread routine, parse each read parallely
//多线程控制函数，任务分配的机制是各线程均分所有的任务并独立完成
void* thread_trimReads(void* threadId_p)
{
	int threadId = *((int*)threadId_p);

	while (1)
	{
		usleep(1);

		if (gotReads[threadId] == 1)
		{
			for (uint64_t i=0; i<readNum; i++)
			{
				if (i%threadNum == threadId)
				{
					//trim the illumina adapter in the tail end of reads
					for (int j=0; j<AdapterVec.size(); j++)
					{	
						int align_score, seq_start, seq_end, contaminant_start, contaminant_end;
						local_ungapped_aligning (StoreReads[i], AdapterVec[j], align_score, seq_start, contaminant_start, seq_end, contaminant_end);
						if (align_score >= Alignment_score_cutoff)
						{	
							int read_len = StoreReads[i].size();
							int trimmed_read_len = seq_start - 1;
							StoreReads[i] = StoreReads[i].substr(0,trimmed_read_len);
							StoreQuals[i] = StoreQuals[i].substr(0,trimmed_read_len);
							StoreHeads[i] += "   Aligned to adapter " + AdapterVecId[j] + ", ";
							StoreHeads[i] += " reads_pos: " + boost::lexical_cast<std::string>(seq_start) + "-" + boost::lexical_cast<std::string>(seq_end) + ", ";
							StoreHeads[i] += "adapter_pos: " + boost::lexical_cast<std::string>(contaminant_start) + "-" + boost::lexical_cast<std::string>(contaminant_end) + ", ";
							StoreHeads[i] += "  score: " + boost::lexical_cast<std::string>(align_score);
							num_adapter_trimmed_reads[threadId] ++;
							num_adapter_trimmed_bases[threadId] += read_len - seq_start + 1;
							break;
						}
					}
					

					//trim the extreme short reads to empty string
					if (StoreReads[i].size() < Minimum_trimmed_read_len)
					{	num_Nmasked_short_reads[threadId] ++;
						num_Nmasked_short_bases[threadId] += StoreReads[i].size();
						StoreReads[i] = "";
						StoreQuals[i] = "";
						StoreHeads[i] += "   RemoveShort";
					}else{
						num_clean_reads[threadId] ++;
						num_clean_bases[threadId] += StoreReads[i].size();
					}
				}
			}

			gotReads[threadId] = 0;
		}
		else if (gotReads[threadId] == 2)
		{
			return NULL;
		}
	}
}

//reading the sequences from fasta-format file
void read_fasta ( string &file_name, vector<string> &seqs, vector<string> &ids )
{	
	ifstream infile( file_name.c_str(), ios::in );
	if ( ! infile )
	{	cerr << "fail to open input file: " << file_name << endl;
		exit( -1 );
	}

	string textline;
	getline( infile, textline, '>' );
	while ( getline( infile, textline, '\n' ) )
	{	vector<string> vechead;
		boost::split(vechead, textline, boost::is_any_of(" \t\n"), boost::token_compress_on);
		
		getline( infile, textline, '>' );
		string purestr;
		for (int i=0; i<textline.size(); i++)
		{	if (textline[i] != '\n' && textline[i] != ' ' && textline[i] != '\t')
			{	
				purestr.push_back( textline[i] );
			} 
		}
		seqs.push_back(purestr);
		ids.push_back(vechead[0]);
		
		if (Also_use_adapterRC == 1)
		{
			string rc_str;
			reverse_complement(purestr,rc_str);
			seqs.push_back(rc_str);
			ids.push_back(vechead[0]+" minus-strand");
		}
	}
	
}



int main(int argc, char *argv[])
{	
	//get options from command line
	int c;
	while((c=getopt(argc, argv, "s:r:a:b:t:h")) !=-1) {
		switch(c) {
			case 's': Alignment_score_cutoff = atoi(optarg); break;
			case 'r': Minimum_trimmed_read_len = atoi(optarg); break;
			case 'a': illumina_adapter_file = optarg; break;
			case 'b': Also_use_adapterRC = atoi(optarg); break;
			case 't': threadNum=atoi(optarg); break;
			case 'h': usage(); break;
			default: usage();
		}
	}
	if (argc < 4) usage();
	
	string in_reads1_file = argv[optind++]; //optind, argv[optind++]顺序指向非option的参数
	string out_reads1_file = argv[optind++]; 
	string out_stat_file = argv[optind++]; 	

	cerr << "\nAlignment score cutoff: " << Alignment_score_cutoff << endl;
	cerr << "Minimum trimmed read length: >=" << Minimum_trimmed_read_len << endl;
	
	clock_t time_start, time_end;
	time_start = clock();
	
/////////////////////////////////////////////////////////////////////////////////////////////////
	
	//if type eqauls Both-adapter, R1-adapter, R2-adapter, use the default adapter sequence file
	if (illumina_adapter_file == "Both-adapter")
	{	illumina_adapter_file = "/qdata1/public/software/install/clean_illumina/illumina_NEB_adapter.fa";
	}else if(illumina_adapter_file == "R1-adapter")
	{	illumina_adapter_file = "/qdata1/public/software/install/clean_illumina/illumina_NEB_adapter_R1.fa";
	}else if(illumina_adapter_file == "R2-adapter")
	{	illumina_adapter_file = "/qdata1/public/software/install/clean_illumina/illumina_NEB_adapter_R2.fa";
	}
	cerr << "\nLoad adapter sequences from " << illumina_adapter_file << endl;

	//load the illumina Adapter sequences into memory
	read_fasta(illumina_adapter_file, AdapterVec, AdapterVecId);

	for (int i=0; i<AdapterVec.size(); i++)
	{	cerr << "Used illumina adapter: " << AdapterVecId[i] << " :   " << AdapterVec[i] << endl;
	}
	if (! AdapterVec.size())
	{	cerr << "Died because illumina_adapter_file does not exist or no contamination sequences are provided\n";
		exit(0);
	}

	cerr << "\nInput Reads  file:  " << in_reads1_file << endl << endl;


/////////////////////////////////////////////////////////////////////////////////////////////////

	igzstream infile1 (in_reads1_file.c_str());
	ogzstream cleanfile1 (out_reads1_file.c_str());	
	ofstream statfile (out_stat_file.c_str());
	
	uint64_t total_raw_reads = 0;
	uint64_t total_raw_bases = 0;
	uint64_t total_clean_reads = 0;
	uint64_t total_clean_bases = 0;
	uint64_t total_adapter_trimmed_reads = 0;
	uint64_t total_adapter_trimmed_bases = 0;
	uint64_t total_Nmasked_short_reads = 0;
	uint64_t total_Nmasked_short_bases = 0;

	StoreReads = new string[bufferNum];
	StoreHeads  = new string[bufferNum];
	StoreQuals  = new string[bufferNum];
	num_clean_reads  = new uint64_t[threadNum];
	num_clean_bases  = new uint64_t[threadNum];
	num_adapter_trimmed_reads  = new uint64_t[threadNum];
	num_adapter_trimmed_bases  = new uint64_t[threadNum];
	num_Nmasked_short_reads  = new uint64_t[threadNum];
	num_Nmasked_short_bases  = new uint64_t[threadNum];

	pthread_t *pthread = new pthread_t[threadNum];
	int *pthreadId = new int[threadNum];
	gotReads = new int[threadNum];
	for (int i=0; i<threadNum; i++)
	{	
		num_clean_reads[i] = 0;
		num_clean_bases[i] = 0;
		num_adapter_trimmed_reads[i] = 0;
		num_adapter_trimmed_bases[i] = 0;
		num_Nmasked_short_reads[i] = 0;
		num_Nmasked_short_bases[i] = 0;
		gotReads[i] = 0;
		pthreadId[i] = i;
		pthread_create((pthread+i), NULL, thread_trimReads, (void*)(pthreadId+i));
	}
	cerr << threadNum << " threads created!\n" << endl;


	//读取每一块数据，在内存中进行并行处理
	while (1)
	{	
		readNum = 0;
		
		//读取bufferNum个reads
		//读取fq格式文件
		string empty_line;
		while ( readNum < bufferNum && getline( infile1, StoreHeads[readNum], '\n' ) )
		{	
			if (StoreHeads[readNum][0] == '@') 
			{	
				getline( infile1, StoreReads[readNum], '\n');
				getline( infile1, empty_line, '\n');
				getline( infile1, StoreQuals[readNum], '\n');
				total_raw_reads ++;
				total_raw_bases += StoreReads[readNum].size();
				readNum ++;
			}	
		}
		
		cerr << "reading num " << readNum << endl;

		//线程运行，直到Loop_id等于readNum, 然后进入子线程等待状态
		//放行主线程，进入下一个循环，读数据
		for (int i=0; i<threadNum; i++)
		{	gotReads[i] = 1;
		}
		
		//等子线程，直到gotReads都是0，说明任务完成，可以继续读下一批数据
		while (1)
		{
			usleep(1);
			int i=0;
			for (; i<threadNum; i++)
			{	if (gotReads[i] == 1)
				{	break;
				}
			}
			if (i == threadNum)
			{	break;
			}
		}
		
		//输出比对结果
		for (int i=0; i<readNum; i++)
		{	cleanfile1 << StoreHeads[i] << "\n" << StoreReads[i] << "\n" << "+" << "\n" << StoreQuals[i] << endl;			
		}
		
		//如果已经到了文件尾，则终止全部子线程，并且退出读文件的循环
		if (readNum < bufferNum)
		{	for (int i=0; i<threadNum; i++)
			{	gotReads[i] = 2;
			}
			break; //读到文件尾退出
		}
	}
	
	//等待全部子线程结束
	for (int i=0; i<threadNum; i++)
	{
		pthread_join(pthread[i], NULL);
	}
	
	//综合统计结果
	for (int i=0; i<threadNum; i++)
	{	total_clean_reads += num_clean_reads[i];
		total_clean_bases += num_clean_bases[i];
		total_adapter_trimmed_reads += num_adapter_trimmed_reads[i];
		total_adapter_trimmed_bases += num_adapter_trimmed_bases[i];
		total_Nmasked_short_reads += num_Nmasked_short_reads[i];
		total_Nmasked_short_bases += num_Nmasked_short_bases[i];
	}

	delete[] num_clean_reads;
	delete[] num_clean_bases;
	delete[] num_adapter_trimmed_reads;
	delete[] num_adapter_trimmed_bases;
	delete[] num_Nmasked_short_reads;
	delete[] num_Nmasked_short_bases;
	delete[] pthread;
	delete[] pthreadId;
	delete[] gotReads;
	delete[] StoreReads;
	delete[] StoreHeads;
	delete[] StoreQuals;


/////////////////////////////////////////////////////////////////////////////////////////////////

	statfile << "total_raw_reads:  " << total_raw_reads << endl;
	statfile << "total_raw_bases:  " << total_raw_bases << endl;
	
	double adapter_ratio = (double)total_adapter_trimmed_bases / total_raw_bases;
	statfile << "total_adapter_trimmed_reads:  " << total_adapter_trimmed_reads << endl;
	statfile << "total_adapter_trimmed_bases:  " << total_adapter_trimmed_bases << "\t" << adapter_ratio << endl;
		
	double Nmasked_ratio = (double)total_Nmasked_short_bases / total_raw_bases;
	statfile << "total_short_trimmed_reads:  " << total_Nmasked_short_reads << endl;
	statfile << "total_short_trimmed_bases:  " << total_Nmasked_short_bases << "\t" << Nmasked_ratio  << endl;
	
	double clean_ratio = (double)total_clean_bases / total_raw_bases;
	statfile << "total_clean_reads:  " << total_clean_reads << endl;
	statfile << "total_clean_bases:  " << total_clean_bases  << "\t" << clean_ratio  << endl;

	time_end = clock();
	cerr << "\nAll jobs finished\n";
	cerr << "Run time: " << double(time_end - time_start) / CLOCKS_PER_SEC << endl;

}

