//纠错程序：根据kmer频率的高低，低频kmer绝大多数由测序错误产生，可以被高频kmer(正确的)进行纠正。
//纠错程序成功的关键，在于kmerfreq表的低频高频区分分明，低频代表错误，高频代表正确
//纠错的关键在于控制trim_len和error_num的平衡，在一条read上累积计算error数，在头或者尾巴部位trim

//KmerSize可用范围预计为15-19，KmerSize越大，适用于的基因组越大，所占内存也越大
//Alert: Treat N as A (sequencing error) automatically, and will correct most of the N bases with normal bases

//Alert: By a combination of correction and trimming, the result of this program will be the sequences which do not contain low-frequency Kmers (default: frequency < 10)
//i.e, there will be no low-frequency Kmers in the result sequences, all the low-frequency-Kmer sequencess have been either corrected or trimmed/deleted.

//This program works well on the repeat-less genome (fruitfly), but it works bad on the repeat-rich genome (maize).
//The trimming and deleting function will remove some sequences, if the genome is repeat-rich, then the repeat sequences may be removed more than unique sequences, due to the difficulty for correction when alternative minimum_change_path may be appeared.

//author: Fan Wei, email: fanw@genomics.org.cn, date: 2011-1-20
//author: Li zhenyu, email: lizhenyu@genomics.org.cn, date: 2011-1-23
//author: Yuan Jianying, email: yuanjianying@genomics.org.cn, date: 2011-1-20
//version: 1.5


#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <set>
#include <cmath>
#include <ctime>
#include <inttypes.h>
#include "seqKmer.h"
#include "correct.h"
#include "global.h"
#include "gzstream.h"

using namespace std;

//set global variables
int KmerSize = 17;            //当KmerSize=17时，4^17 = 2^34 = 16G, store 1bit for a kmer, kmer频率表将占内存2G
int Low_freq_cutoff = 10; 
int HighFreqRegLenCutoff = 0; //在main主函数数里赋值=KmerSize
int Max_change_in_one_read = 2;
int Further_trim_len = KmerSize / 2;  //trim the low-quality end instead of correction
int Min_trimmed_read_len = 50;
int Max_node_in_BB_tree = 15000000;  //一个branch and bound tree将占15M x 4 x 2 = 120 M bytes memeory, 8个线程占1G
int Input_file_format = 1;
int Join_read1_read2_file = 1;
uint8_t *KmerFreq;
uint64_t SrcBlockSize = 8*1024*1024;

//read the kmers from 8bit gz file into memory, store frequency value in 1 bit, range 0-1
uint8_t* make_kmerFreq_1bit_table_from_8BitGz(string &kmer_freq_file, int Ksize, uint64_t &total, uint64_t &num_total_kmers, uint64_t &num_effect_kmers, int low_freq_cutoff, double &low_freq_ratio);


void usage() 
{	cout << "\ncorrect_error (v2.2) <kmer_freq.cz> <reads_file_list>\n" 
		 << "   -k <int>   set kmer size (be same in kmer_freq and correct_error), default=" << KmerSize << endl
		 << "   -l <int>   set the low frequency cutoff(sites need correct), default=" << Low_freq_cutoff << endl 
		 << "   -m <int>   set the minimum length of a continuous high-freq-kmer region, default=" << HighFreqRegLenCutoff << endl
		 << "   -c <int>   set the maximum change allowed in one read, default=" << Max_change_in_one_read << endl
		 << "   -x <int>   set the trimmed length at low-quality ends instead of correct them, default=" <<  Further_trim_len << endl
		 << "   -n <int>   set the maximum node number allowed in the branch and bound tree, default="  << Max_node_in_BB_tree << endl
	     << "   -r <int>   set the minimum length of trimmed read, default=" << Min_trimmed_read_len << endl
		 << "   -p <str>   set the prefix of the output files, default=output" << endl
		 << "   -f <int>   set the input file format: 1: fq, 2: fa, default=" << Input_file_format << endl
		 << "   -j <int>   set whether join read1 and read2 result file: 1, join, 0, no; default=" << Join_read1_read2_file << endl
		 << "   -h         get the help information\n" << endl;
	exit(0);
}


int main(int argc, char *argv[])
{	
	string prefix = "output";

	//get options from command line
	int c;
	while((c=getopt(argc, argv, "k:l:m:c:x:n:r:p:f:j:h")) !=-1) {
		switch(c) {
			case 'k': KmerSize=atoi(optarg); break;
			case 'l': Low_freq_cutoff=atoi(optarg); break;
			case 'm': HighFreqRegLenCutoff=atoi(optarg); break;
			case 'c': Max_change_in_one_read=atoi(optarg); break;
			case 'x': Further_trim_len=atoi(optarg); break;
			case 'n': Max_node_in_BB_tree=atoi(optarg); break;
			case 'r': Min_trimmed_read_len=atoi(optarg); break;
			case 'p': prefix=optarg; break;
			case 'f': Input_file_format=atoi(optarg); break;
			case 'j': Join_read1_read2_file=atoi(optarg); break;
			case 'h': usage(); break;
			default: usage();
		}
	}
	if (HighFreqRegLenCutoff == 0)
	{	HighFreqRegLenCutoff = KmerSize;
	}

	if (argc < 2) usage();

	clock_t time_start, time_end;
	time_start = clock();
	
	string kmer_freq_file = argv[optind++]; //optind, argv[optind++]顺序指向非option的参数
	string reads_file_list = argv[optind++];
	uint64_t Kmer_theory_total = 0;
	uint64_t Kmer_total_num = 0;
	uint64_t Kmer_effect_num = 0;
	double low_freq_kmer_ratio = 0.0;

	KmerFreq = make_kmerFreq_1bit_table_from_8BitGz(kmer_freq_file, KmerSize, Kmer_theory_total, Kmer_total_num, Kmer_effect_num, Low_freq_cutoff, low_freq_kmer_ratio);
	
	cerr << "Load the kmer frequency table completed" << endl;
	cerr << "KmerSize " << KmerSize << endl;
	cerr << "Kmer_effect_num " << Kmer_effect_num << endl;
	cerr << "low_freq_kmer_ratio " << low_freq_kmer_ratio << endl;
	time_end = clock();
	cerr << "Run time: " << double(time_end - time_start) / CLOCKS_PER_SEC << endl;

	//Correct errors and do statistics file by file
	vector <string> reads_files;
	reading_file_list(reads_file_list, reads_files);
	
	for (int i = 0; i < reads_files.size(); i++)
	{	
		if (Input_file_format == 1)
		{	
			parse_one_reads_fq_file(reads_files[i]);
		}
		else
		{
			parse_one_reads_fa_file(reads_files[i]);
		}
	}
	
	delete[] KmerFreq;
	
	cerr << "Paresed all the reads files completed" << endl;
	time_end = clock();
	cerr << "Run time: " << double(time_end - time_start) / CLOCKS_PER_SEC << endl;
	
	//merge the read1 and read2 files into pair and single file
	if (Join_read1_read2_file != 1)
	{	return 0;
	}

	for (int i = 0; i < reads_files.size(); i+=2)
	{	string read1_file = reads_files[i] + ".cor";
		string read2_file = reads_files[i+1] + ".cor";
		merge_two_corr_files(read1_file, read2_file);
		string rm_cmd = "rm " + read1_file + " " + read2_file;
		system(rm_cmd.c_str());
	}
	
	cerr << "Merge all the paired reads files completed" << endl;
	time_end = clock();
	cerr << "Run time: " << double(time_end - time_start) / CLOCKS_PER_SEC << endl;

}


//read the kmers from the compressed kmer frequency talbe file into memory, store frequency value in 1 bit, range 0-1
uint8_t* make_kmerFreq_1bit_table_from_8BitGz(string &kmer_freq_file, int Ksize, uint64_t &total, uint64_t &num_total_kmers, uint64_t &num_effect_kmers, int low_freq_cutoff, double &low_freq_ratio)
{	
	//得到理论上最大的kmer个数
	total=0;  
	for(int i=0; i<Ksize; i++) {
		total=(total<<2)|0x3;
	}
	total++;
	
	//动态分配的数组必须赋初值
	uint64_t array_size = total/8;
	uint8_t *freq=new uint8_t[array_size]; 
	memset(freq, 0, array_size);
	
	//读取kmer-frequency-table压缩文件的长度信息文件
	string kmer_freq_len_file = kmer_freq_file + ".len";
	vector<uint64_t> BlockSizeVec;
	uint8_t *buffer_compress = new uint8_t[SrcBlockSize];
	uint8_t *buffer_uncompress = new uint8_t[SrcBlockSize];
	
	ifstream LenInfoFile (kmer_freq_len_file.c_str());
	string Line_str;
	while (getline( LenInfoFile, Line_str, '\n' ))
	{	BlockSizeVec.push_back(atoi(Line_str.c_str()));
	}
	
	//读取kmer-frequency-table压缩文件, 并完成高频选择，和由byte到bit的转换
	uint64_t low_freq_num = 0;
	ifstream CompDataFile (kmer_freq_file.c_str(), ifstream::binary);
	for (uint64_t i = 0; i < BlockSizeVec.size(); i ++)
	{	uint64_t BlockStartPos = i * SrcBlockSize;
		CompDataFile.read((char *)buffer_compress, BlockSizeVec[i]);
		uLongf uncompLen = SrcBlockSize;
		uncompress((Bytef *)buffer_uncompress, &uncompLen, (const Bytef*)buffer_compress, BlockSizeVec[i]);
		
		for (uint64_t j = 0; j < uncompLen; j ++)
		{	uint64_t idx = BlockStartPos + j;
			if (buffer_uncompress[j] > 0)
			{	num_total_kmers += buffer_uncompress[j]; 
				num_effect_kmers ++;
				
				if(buffer_uncompress[j] > low_freq_cutoff)
				{	
					freq[idx/8] |= bitAll[idx%8];
					
					//只分成两个档,低频和高频, 在此处将高频转换成正负链都有的形式
					//能自动处理KmerSize为偶数的情况，正负链形式存在一个地方
					uint64_t rc_idx = get_rev_com_kbit(idx, Ksize);
					freq[rc_idx/8] |= bitAll[rc_idx%8];
				}else{
					low_freq_num ++;
				}
			}
		}
	}

	low_freq_ratio = (double)low_freq_num / (double)num_effect_kmers;

	return freq;
}

