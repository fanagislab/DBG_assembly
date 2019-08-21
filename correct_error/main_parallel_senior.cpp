//纠错程序：根据kmer频率的高低，低频kmer绝大多数由测序错误产生，可以被高频kmer(正确的)进行纠正。
//纠错程序成功的关键，在于kmerfreq表的低频高频区分分明，低频代表错误，高频代表正确
//纠错的关键在于控制trim_len和error_num的平衡，在一条read上累积计算error数，在头或者尾巴部位trim

//KmerSize可用范围预计为15-19，KmerSize越大，适用于的基因组越大，所占内存也越大
//Alert: Treat N as A (sequencing error) automatically, and will correct most of the N bases with normal bases

//Alert: By a combination of correction and trimming, the result of this program will be the sequences which do not contain low-frequency Kmers 
//i.e, there will be no low-frequency Kmers in the result sequences, all the low-frequency-Kmer sequencess have been either corrected or trimmed/deleted.

//This program works well on the repeat-less genome (fruitfly), but it works bad on the repeat-rich genome (maize).
//The trimming and deleting function will remove some sequences, if the genome is repeat-rich, then the repeat sequences may be removed more than unique sequences, due to the difficulty for correction when alternative minimum_change_path may be appeared.

//author: Fan Wei, email: fanw@genomics.org.cn, date: 2015-11-13
//version: 2.0

//This is the multiple-threads version of correct_error program.
//Read bufferNum reads into memory each time, and parsed in parallel with threadNum threads.

//Work process: 
//(1) divide all the kmers into low frequency and high frequency regions;
//(2) correct simple errors (length of low frequency kmers equal to KmerSize) by the simple correcting algorithm
//    after correction, the corresponding low frequency region becomes high frequency regions
//(3) get high frequnecy regions with length larger than a cutoff, smaller high frequency regions are considered not reliable which may contain sequencing error
//(4) correct the regions between two large high frequency regions by the branch and bond tree method,  if fail from one direction, then try from the other direction
//(5) correct the read ends (read heads and read tails), if they are low frequency regions, by the branch and bond tree method, the uncorrected regions are trimmed off. To increase the result accuray, when trimming at the read heads and tails happens, the program trim a futher length at the read ends in addtion.


//New Functions:
//This program version read the kmer frequency table in one-bit format with low and high frequency differetiated, which is generated from kmerfreq_twobyte.


#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <set>
#include <cmath>
#include <ctime>
#include <pthread.h>
#include <inttypes.h>
#include "seqKmer.h"
#include "correct.h"
#include "global.h"
#include "gzstream.h"

using namespace std;

//set global variables, used by other cpp files
int KmerSize = 17;       //当KmerSize=17时，4^17 = 2^34 = 16G, store 1bit for a kmer, kmer频率表将占内存2G
int HighFreqRegLenCutoff = KmerSize; //在main主函数数里赋值=KmerSize
int Max_change_in_one_read = 2;  //每条reads最多允许修改的碱基数目，当设置为0时，则程序对于低频区域不做任何纠错，而是将低频区域全部删除掉
int Further_trim_len = KmerSize;  //After correction by the branch and bond method, Further trim the low-quality read heads and ends to increase result accuary
int Min_trimmed_read_len = 75;
int Max_node_in_BB_tree = 5000000;  //一个branch and bound tree将占5M x 4 x 2 = 40 M bytes memeory, 25个线程同时发生最大情况将共占1G，此情况属于极限,几乎不可能发生
int Input_file_format = 1;
int Join_read1_read2_file = 0;
int threadNum = 10;
uint8_t *KmerFreq;

uint64_t Kmer_theory_total = 0;
uint64_t Kmer_hifreq_num = 0;


//thread variables and routines to load the kmers from compressed file into memory
uint8_t* make_kmerFreq_1bit_table_from_1BitGz_pthread(string &kmer_freq_file, int Ksize, uint64_t &total);
uint64_t bufBlockId = 0;
uint64_t bufBlockNum = 0; //calculate in the program
uint64_t SrcBlockSize = 8*1024*1024;   //8M kmers for each compressing block, must be set the same size as variable bufBlockSize in program kmerfreq_twobyte;
vector<uint64_t> BlockSizeVec;
uint8_t *buffer_compress;
uLongf *uncompLen;
ifstream Compr_KmerFreq_File;
pthread_mutex_t Block_mutex = PTHREAD_RECURSIVE_MUTEX_INITIALIZER_NP;
void *thread_uncompressKmerFreq(void* threadId_p);
void *thread_setrevcompkmer(void* threadId_p);


//thread variables used to parse reads blocks parallely
uint64_t bufferNum = 1000000;  //read bufferNum reads as one block into memory each time
int *gotReads;
uint64_t readNum;
string *RawReads;
string *HeadStr;
int *Score;
int *LeftTrim;
int *RightTrim;
int *IsDeleted;
uint64_t *numResReads;
uint64_t *numResBases;
uint64_t *numTrimmedReads;
uint64_t *numTrimmedBases;
uint64_t *numDeletedReads;
uint64_t *OneBaseCorrectScore;
uint64_t *MultiBaseCorrectScore;
void* thread_correctReads(void* threadId_p);
void parse_one_reads_file( string &reads_file);


//thread variables used to merge corrected paired files parallely
pthread_mutex_t File_mutex = PTHREAD_RECURSIVE_MUTEX_INITIALIZER_NP;
int File_num = 0;
void *thread_mergeTwoFile(void *data);


void usage() 
{	
	cout << "\n1. Function introduction:\
	\nThis program corrects the sequencing errors by using the kmer frequency information. It assume that most low-frequency Kmers were generated by sequencing errors, so the key of error correction is that the distinguish rate of the low-frequency and high-frequency Kmers, the larger Kmer size, the better of this effect. In order to get an extreme high accuracy result, we balanced the trimmed length and delete ratio with the accuracy level. The practical Kmer size is from 15 to 19, which should be chosen based on the genome size. Alerts: (1) This program treats base N as A (sequencing error) automatically, and will correct most of the N bases with normal bases. (2) By a combination of correction and trimming, the result of this program will be the sequences which do not contain low-frequency Kmers (default: frequency < 10), i.e, there will be no low-frequency Kmers in the result sequences, all the low-frequency-Kmer sequencess have been either corrected or trimmed/deleted, this will filter some real reads from low-coverage genomic regions, and decrease contig length in assembly. (3)This program works well on the repeat-less genome (fruitfly), but it works badly on the repeat-rich genome (maize). The trimming and deleting function will remove some sequences, if the genome is repeat-rich, then the repeat sequences may be removed more than unique sequences, due to the difficulty for correction when alternative minimum_change_path may be appeared. \
\n\n2.Input and output\
	\nThe first input is the kmer frequency table storing high and low frequency kmers, which was generated by kmerfreq_16bit. The second input is a libary file that contains the path of all the sequencing reads files in fastq format or one-line fasta format. The program parse each sequencing reads file by the given order, and make a whole statistics. Output are corrected reads files in fasta format, as well as statistics files. The program run in a multiple-threads way.\n\n";
	
	
	cout << "\ncorrect_error_reads <kmer_freq.cz> <reads_file.lib>\n" 
		 << "   -k <int>   kmer size, must be set same with that in in kmerfreq_16bit, default=" << KmerSize << endl
		 << "   -m <int>   the minimum length of a continuous high-freq-kmer region, default=KmerSize: " << HighFreqRegLenCutoff << endl
		 << "   -c <int>   the maximum change allowed in one read, default=" << Max_change_in_one_read << endl
		 << "   -x <int>   the trimmed length at low-quality ends instead of correcting them, default=KmerSize: " <<  Further_trim_len << endl
		 << "   -n <int>   the maximum node number allowed in the branch and bound tree, default="  << Max_node_in_BB_tree << endl
	     << "   -r <int>   the minimum length of trimmed read, default=" << Min_trimmed_read_len << endl
		 << "   -t <int>   thread number to use, default=" << threadNum << endl
		 //<< "   -p <str>   set the prefix of the output files, default=output" << endl
		 << "   -f <int>   the input file format: 1: fq, 2: fa, default=" << Input_file_format << endl
		 << "   -j <int>   whether join read1 and read2 result file: 1, join, 0, no; default=" << Join_read1_read2_file << endl
		 << "   -h         get the help information\n" << endl
		 << "Examples: " << endl
		 << "1. Remove all the reads regions with low kmer frequency \n"
		 << "   correct_error_reads -c 0  random1M.freq.cz  random1M.lib\n"
		 << "2. Allow correct no more than 2 bases in one reads, and remove the other read regions with low kmer frequency\n"
		 << "   correct_error_reads -c 2  random1M.freq.cz  random1M.lib\n"
		 << "\nNote: " << endl
		 << "   THe file kmer_freq.cz used here is generated from program kmerfreq_16bit,  the .cz.len file should be put in the same directory\n"
		 << endl;


	exit(0);
}


int main(int argc, char *argv[])
{	
	string prefix = "output";

	//get options from command line
	int c;
	while((c=getopt(argc, argv, "k:m:c:x:n:r:t:p:f:j:h")) !=-1) {
		switch(c) {
			case 'k': KmerSize=atoi(optarg); break;
			case 'm': HighFreqRegLenCutoff=atoi(optarg); break;
			case 'c': Max_change_in_one_read=atoi(optarg); break;
			case 'x': Further_trim_len=atoi(optarg); break;
			case 'n': Max_node_in_BB_tree=atoi(optarg); break;
			case 'r': Min_trimmed_read_len=atoi(optarg); break;
			case 't': threadNum=atoi(optarg); break;
			case 'p': prefix=optarg; break;
			case 'f': Input_file_format=atoi(optarg); break;
			case 'j': Join_read1_read2_file=atoi(optarg); break;
			case 'h': usage(); break;
			default: usage();
		}
	}
	if (argc < 2) usage();
	
	clock_t time_start, time_end;
	time_start = clock();
	
	string kmer_freq_file = argv[optind++]; //optind, argv[optind++]顺序指向非option的参数
	string reads_file_list = argv[optind++];
	
	//load compressed kmer-frequency-table into memory parallely
	make_kmerFreq_1bit_table_from_1BitGz_pthread(kmer_freq_file, KmerSize, Kmer_theory_total);
	
	cerr << "Load the kmer frequency table completed" << endl;
	cerr << "KmerSize " << KmerSize << endl;
	cerr << "Kmer_theory_total " << Kmer_theory_total << endl;
	cerr << "Kmer_hifreq_num   " << Kmer_hifreq_num << endl;
	time_end = clock();
	cerr << "Run time: " << double(time_end - time_start) / CLOCKS_PER_SEC << endl;

///////////////////////////////////////////////////////////////////////////////////

	//开指定个子线程,并分配相关资源
	RawReads = new string[bufferNum];
	HeadStr	= new string[bufferNum];
	Score = new int[bufferNum];
	LeftTrim = new int[bufferNum];
	RightTrim = new int[bufferNum];
	IsDeleted = new int[bufferNum];

	numResReads = new uint64_t[threadNum];
	numResBases = new uint64_t[threadNum];
	numTrimmedReads = new uint64_t[threadNum];
	numTrimmedBases  = new uint64_t[threadNum];
	numDeletedReads  = new uint64_t[threadNum];
	OneBaseCorrectScore = new uint64_t[threadNum];
	MultiBaseCorrectScore  = new uint64_t[threadNum];
	gotReads = new int[threadNum];	
	
	//Correct errors and do statistics file by file
	vector <string> reads_files;
	reading_file_list(reads_file_list, reads_files);
	
	for (int file_i = 0; file_i < reads_files.size(); file_i++)
	{	
		parse_one_reads_file(reads_files[file_i]);
	}
	
	//释放动态分配的内存
	delete[] KmerFreq;
	delete[] RawReads;
	delete[] HeadStr;
	delete[] Score;
	delete[] LeftTrim;
	delete[] RightTrim;
	delete[] IsDeleted;
	delete[] numResReads;
	delete[] numResBases;
	delete[] numTrimmedReads;
	delete[] numTrimmedBases;
	delete[] numDeletedReads;
	delete[] OneBaseCorrectScore;
	delete[] MultiBaseCorrectScore;
	delete[] gotReads;

	cerr << "Parsed all the reads files completed" << endl;
	time_end = clock();
	cerr << "Run time: " << double(time_end - time_start) / CLOCKS_PER_SEC << endl;

///////////////////////////////////////////////////////////////////////////////////

	//merge the read1 and read2 files into pair and single file
	if (Join_read1_read2_file != 1)
	{	return 0;
	}

	File_num = 0;
	pthread_t *pthread = new pthread_t[threadNum];
	
	if (pthread == NULL)
	{	cerr << "Out of memory!" << endl;
		exit(0);
	}
	
	for (int i=0; i<threadNum; i++)
	{ /* create children threads */
		pthread_create(&pthread[i], NULL, thread_mergeTwoFile, (void*)(&reads_files));
	}

	/* wait threads to exit and free resource after threads exit */
	for (int i=0; i<threadNum; i++)
	{ 
		pthread_join(pthread[i], NULL);
	}
	
	//delete the old files of read1 and read2
	for (int i = 0; i < reads_files.size(); i+=2)
	{	string read1_file = reads_files[i] + ".correct.fa.gz";
		string read2_file = reads_files[i+1] + ".correct.fa.gz";
		string rm_cmd = "rm " + read1_file + " " + read2_file;
		system(rm_cmd.c_str());
	}

	cerr << "Merge all the paired reads files completed" << endl;
	time_end = clock();
	cerr << "Run time: " << double(time_end - time_start) / CLOCKS_PER_SEC << endl;

}



void *thread_uncompressKmerFreq(void* threadId_p)
{
	int threadId = *((int*)threadId_p);
	
	while (1)
	{	
		usleep(1);
		pthread_mutex_lock(&Block_mutex);
		
		if (bufBlockId < bufBlockNum)
		{	
			uint64_t this_block_id = bufBlockId ++;  //because of multi-threads, one thread firstly get the original value, and then bufBlockId add one.
			Compr_KmerFreq_File.read((char *)(buffer_compress+threadId*SrcBlockSize/8), BlockSizeVec[this_block_id]);

			pthread_mutex_unlock(&Block_mutex);

			///////////////////////////////////////////////////////////////////////////////////////
			
			//对读入缓存的内容进行解压缩
			uncompLen[threadId] = SrcBlockSize/8;  //Must do not forget to assign initital value to this variable
			
			
			uncompress((Bytef *)(KmerFreq+this_block_id*SrcBlockSize/8), &uncompLen[threadId], (const Bytef*)(buffer_compress+threadId*SrcBlockSize/8), BlockSizeVec[this_block_id]);
			
			
			///////////////////////////////////////////////////////////////////////////////////////


		}
		else
		{	pthread_mutex_unlock(&Block_mutex);
			break;
		}
	}
}


void *thread_setrevcompkmer(void* threadId_p)
{
	int threadId = *((int*)threadId_p);

	for (uint64_t i = threadId; i< Kmer_theory_total; i += threadNum )
	{	
		if ( get_freq( KmerFreq, i ) )
		{	
			uint64_t rc_i = get_rev_com_kbit(i, KmerSize);
			
			if (i <= rc_i)
			{
				__sync_or_and_fetch( &KmerFreq[rc_i/8], bitAll[rc_i%8] );  //use atomic-process
				__sync_add_and_fetch( &Kmer_hifreq_num, 1 );
			}

		}
	}

}



//read the kmers from the compressed kmer frequency talbe file into memory, store frequency value in 1 bit, range 0-1
uint8_t* make_kmerFreq_1bit_table_from_1BitGz_pthread(string &kmer_freq_file, int Ksize, uint64_t &total)
{	
	//得到理论上最大的kmer个数
	total=0;  
	for(int i=0; i<Ksize; i++) {
		total=(total<<2)|0x3;
	}
	total++;

	
	//动态分配的数组必须赋初值
	uint64_t array_size = total/8;  //一定能够整除的
	KmerFreq = new uint8_t[array_size]; 
	memset(KmerFreq, 0, array_size);
	
	//读取kmer-frequency-table压缩文件的长度信息文件
	string kmer_freq_len_file = kmer_freq_file + ".len";
	ifstream LenInfoFile (kmer_freq_len_file.c_str());
	if (! LenInfoFile)
	{	cerr << "Fail to open " << kmer_freq_len_file << endl;
		exit(0);
	}
	string Line_str;
	while (getline( LenInfoFile, Line_str, '\n' ))
	{	BlockSizeVec.push_back(atoi(Line_str.c_str()));
	}
	bufBlockNum = BlockSizeVec.size();

	//读取kmer-frequency-table压缩文件，执行并行解压缩过程
	buffer_compress = new uint8_t[SrcBlockSize/8*threadNum];
	uncompLen = new uLongf[threadNum];

	Compr_KmerFreq_File.open(kmer_freq_file.c_str(), ifstream::binary);

	pthread_t *pthread = new pthread_t[threadNum];
	int *pthreadId = new int[threadNum];
	
	//create children threads
	for (int i=0; i<threadNum; i++)
	{	
		pthreadId[i] = i;
		uncompLen[i] = SrcBlockSize/8;  //Must do not forget to assign initital value to this variable
		pthread_create(pthread+i, NULL, thread_uncompressKmerFreq, (void*)(pthreadId+i));
	}

	//等待全部子线程结束
	for (int i=0; i<threadNum; i++)
	{	pthread_join(pthread[i], NULL);
	}
	
			
	delete[] buffer_compress;
	delete[] uncompLen;
	
	Compr_KmerFreq_File.close();


	//set the reverse and complement kmer (rc_kbit)
	pthread = new pthread_t[threadNum];
	pthreadId = new int[threadNum];
	
	//create children threads
	for (int i=0; i<threadNum; i++)
	{	
		pthreadId[i] = i;
		pthread_create(pthread+i, NULL,  thread_setrevcompkmer, (void*)(pthreadId+i));
	}

	//等待全部子线程结束
	for (int i=0; i<threadNum; i++)
	{	pthread_join(pthread[i], NULL);
	}
	

}



//this is the thread routine, parse each read parallely
void* thread_correctReads(void* threadId_p)
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
					int one_base_correct_score = 0;
					int multi_base_correct_score = 0;
					int all_base_correct_score = 0;
					int is_deleted = 0;
					int trim_left_end_len = 0;
					int trim_right_end_len = 0;
					int final_read_len = 0;
					
					correct_one_read(HeadStr[i], RawReads[i], KmerFreq, one_base_correct_score, multi_base_correct_score, is_deleted, trim_left_end_len, trim_right_end_len);

					all_base_correct_score = one_base_correct_score + multi_base_correct_score;
					final_read_len = RawReads[i].size() - trim_left_end_len - trim_right_end_len;
					Score[i] = all_base_correct_score;
					LeftTrim[i] = trim_left_end_len;
					RightTrim[i] = trim_right_end_len;
					IsDeleted[i] = is_deleted;

					if ( ! is_deleted )
					{	
						OneBaseCorrectScore[threadId] += one_base_correct_score;
						MultiBaseCorrectScore[threadId] += multi_base_correct_score;

						if (trim_left_end_len > 0 || trim_right_end_len > 0)
						{	
							RawReads[i] = RawReads[i].substr(trim_left_end_len, final_read_len);
							numTrimmedReads[threadId] ++;
							numTrimmedBases[threadId] += trim_left_end_len + trim_right_end_len;
						}
						numResReads[threadId] ++;
						numResBases[threadId] += final_read_len;
					}
					else
					{	numDeletedReads[threadId] ++;
						RawReads[i] = "";
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


//this is the thread routine to merge the read1 and read2 corrected files parallely
void *thread_mergeTwoFile(void *data)
{

	vector <string> *file_vec_p = (vector <string>*)data;
	int vec_size = (*file_vec_p).size();

	while (1)
	{
		pthread_mutex_lock(&File_mutex);
	
		if (File_num < vec_size)
		{ 
			/* get an un-processed file*/
			string read1_file = (*file_vec_p)[File_num++] + ".correct.fa.gz";
			string read2_file = (*file_vec_p)[File_num++] + ".correct.fa.gz";
			
			pthread_mutex_unlock(&File_mutex);
			
			merge_two_corr_files(read1_file, read2_file);
		}
		else
		{
			pthread_mutex_unlock(&File_mutex);
			break;
		}
	}

	return NULL;
}


//parse one input reads file
void parse_one_reads_file( string &reads_file)
{	
	
	igzstream infile ( reads_file.c_str() );
	if ( ! infile )
	{	cerr << "fail to open input file " << reads_file << endl;
	}
	
	string reads_file_cor = reads_file + ".correct.fa.gz";
	ogzstream corFile (reads_file_cor.c_str());
	if ( ! corFile )
	{	cerr << "fail to open output file " << reads_file_cor << endl;
	}
	
	//Statistic numbers for a input reads file
	uint64_t num_raw_reads = 0;
	uint64_t num_raw_bases = 0;
	uint64_t num_res_reads = 0;
	uint64_t num_res_bases = 0;
	uint64_t num_trimmed_reads = 0;
	uint64_t num_trimmed_bases = 0;
	uint64_t num_deleted_reads = 0;
	uint64_t total_one_base_correct_score = 0;
	uint64_t total_multi_base_correct_score = 0;
	uint64_t total_all_base_correct_score = 0;

	pthread_t *pthread = new pthread_t[threadNum];
	int *pthreadId = new int[threadNum];
	for (int i=0; i<threadNum; i++)
	{	
		numResReads[i] = 0;
		numResBases[i] = 0;
		numTrimmedReads[i] = 0;
		numTrimmedBases[i] = 0;
		numDeletedReads[i] = 0;
		OneBaseCorrectScore[i] = 0;
		MultiBaseCorrectScore[i] = 0;
		gotReads[i] = 0;
		pthreadId[i] = i;
		pthread_create((pthread+i), NULL, thread_correctReads, (void*)(pthreadId+i));
	}
	cerr << threadNum << " threads creation done!\n" << endl;

	//读取每一块数据，在内存中进行并行处理
	while (1)
	{	
		readNum = 0;
		
		//读取bufferNum个reads
		if (Input_file_format == 1)
		{	//读取fq格式文件
			string empty_line;
			while ( readNum < bufferNum && getline( infile, HeadStr[readNum], '\n' ) )
			{	
				if (HeadStr[readNum][0] == '@') 
				{	
					HeadStr[readNum][0] = '>'; //the output is fasta format
					getline( infile, RawReads[readNum], '\n');
					getline( infile, empty_line, '\n');
					getline( infile, empty_line, '\n');
					Score[readNum] = 0;
					LeftTrim[readNum] = 0;
					RightTrim[readNum] = 0;
					IsDeleted[readNum] = 0;

					num_raw_reads ++;
					num_raw_bases += RawReads[readNum].size();
					readNum ++;
				}	
			}
		}
		else
		{	//读取fa格式文件
			while ( readNum < bufferNum && getline( infile, HeadStr[readNum], '\n' ) )
			{	
				if (HeadStr[readNum][0] == '>') 
				{	
					getline( infile, RawReads[readNum], '\n');
					Score[readNum] = 0;
					LeftTrim[readNum] = 0;
					RightTrim[readNum] = 0;
					IsDeleted[readNum] = 0;

					num_raw_reads ++;
					num_raw_bases += RawReads[readNum].size();
					readNum ++;
				}	
			}
		}
		
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
		
		//输出纠错结果
		for (int i=0; i<readNum; i++)
		{	
			corFile << HeadStr[i]  << "\tModifiedBaseNum: " << Score[i] << "\tFinalReadLength: " << RawReads[i].size() << "\tLeftEndTrim: " << LeftTrim[i] << "\tRightEndTrim: " << RightTrim[i] << "\tIsDeleted: " << IsDeleted[i]  << "\n" << RawReads[i] << endl;
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

	//综合各线程的统计结果
	for (int i=0; i<threadNum; i++)
	{	num_res_reads += numResReads[i];
		num_res_bases  += numResBases[i];            
		num_trimmed_reads  += numTrimmedReads[i];             
		num_trimmed_bases  += numTrimmedBases[i];          
		num_deleted_reads  += numDeletedReads[i];          
		total_one_base_correct_score += OneBaseCorrectScore[i];
		total_multi_base_correct_score  += MultiBaseCorrectScore[i];
	}
	total_all_base_correct_score = total_one_base_correct_score + total_multi_base_correct_score;

	//output statistic result
	string reads_file_cor_stat = reads_file + ".correct.stat";
	ofstream staFile (reads_file_cor_stat.c_str());
	if ( ! staFile )
	{	cerr << "fail to open output file " << reads_file_cor_stat << endl;
	}

	double data_filter_ratio = (double)(num_raw_bases - num_res_bases) / (double)num_raw_bases;
	double corrected_error_ratio = (double)total_all_base_correct_score / (double)num_res_bases;

	staFile << "num_raw_reads " << num_raw_reads << endl;
	staFile << "num_raw_bases " << num_raw_bases << endl;
	staFile << "num_result_reads " << num_res_reads << endl;
	staFile << "num_result_bases " << num_res_bases << endl;
	
	staFile << "\nnum_trimmed_reads " << num_trimmed_reads << endl;
	staFile << "num_trimmed_bases " << num_trimmed_bases << endl;
	staFile << "num_deleted_reads " << num_deleted_reads << endl;
	
	staFile << "\nnum_corrected_bases_by_Fast_method " << total_one_base_correct_score << endl;
	staFile << "num_corrected_bases_by_BBtree_method " << total_multi_base_correct_score << endl;
	staFile << "num_corrected_bases_by_two_methods " << total_all_base_correct_score << endl;
	
	staFile << "\nfilter_ratio: (num_raw_bases - num_res_bases) / num_raw_bases " << data_filter_ratio << endl;
	staFile << "correct_ratio: total_all_base_correct_score / num_res_bases " << corrected_error_ratio << endl;
	
	cerr << "Finished to parse reads file: " << reads_file << endl;

}


