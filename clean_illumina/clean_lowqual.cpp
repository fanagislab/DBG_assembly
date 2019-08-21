// get high quality reads from raw sequencing data by trimming the low-quality parts.
// output the largest read block in which the average error rate is lower than the Error_rate_cutoff.
// The final left read length should be larger than Min_read_len, else it will be filtered out.
// The average read error rate (RQ) will be added in the read heads. 
// Run in multi-threads mode: read bufferNum reads as one block into memory each time, followed by parallell processing, and single-channel output

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <set>
#include <cmath>
#include <inttypes.h>
#include <pthread.h>
#include<zlib.h>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include "gzstream.h"

using namespace std;

double Error_rate_cutoff = 0.001;
int Min_read_len = 75;
int Quality_shift = 33;
double Qual2Err[128];
int threadNum = 3;
uint64_t bufferNum = 1000000;  //read bufferNum reads as one block into memory each time
uint64_t readNum;
string *StoreHeads;
string *StoreReads;
string *StoreQuals;
void *thread_cleanlowqual(void* threadId_p);

//在这里初始值赋值为0不管用，必须在主程序中赋初始值
uint64_t total_raw_reads;
uint64_t total_raw_bases;
uint64_t total_filtered_lowqual_reads;
uint64_t total_filtered_lowqual_bases;
uint64_t total_filtered_short_reads;
uint64_t total_filtered_short_bases;
uint64_t total_clean_reads;
uint64_t total_clean_bases;


void usage() 
{	
	cout << "Description:\nclean_lowqual detects and trims the low quality parts in a reads, and output a largest block in which the average error rate is lower than a given cutoff(-e). The programs gets all the blocks in a reads with average error rate lower than the cutoff (-e), and choose the longest block as the final result (trimmed reads).   The illumina sequencing machine produces reads with average error rate of 1%, however, the error bases are not distributed evenly. By set a cutoff (-e) of 0.1% for the average error rate, we can exclude most of(90%) of the sequencing errors by filtering out only a small ratio read sequences. The program also has a function to filter the extreme short reads by using a cutoff (-r); The program runs in a multiple thread mode (-t). The input file should be fastq or gzipped-fasts format, and there are two result files: one is the high-quality reads file, and the other one is a statistics file.\n";
	cout << "\nUsage:\n  clean_lowqual <input.fq.gz>  <output.fq.gz>  <output.stat>" << endl;
	cout << "   Author: Fanwei, fanweiagis@126.com" << endl;
	cout << "   Version 1.0;" << endl;
	cout << "   -e <float>  average error rate cutoff for a read, default=" << Error_rate_cutoff << endl;
	cout << "   -q <int>    base quality shift value, default=" << Quality_shift << endl;
	cout << "   -r <int>    minimum read length for output, default=" << Min_read_len << endl;
	cout << "   -t <int>    thread number to run in parallel, default=" << threadNum << endl;
	cout << "   -h          get help information" << endl << endl;
	cout << "Example:\n  ../clean_lowqual -e 0.001 -r 75 sequencing_data_R1.fq.gz sequencing_data_R1.fq.gz.nonLowQual.gz sequencing_data_R1.fq.gz.nonLowQual.stat\n" << endl;
	exit(0);
}


//this is the thread routine, parse each read parallely
//多线程控制函数，任务分配的机制是各线程均分所有的任务并独立完成
void* thread_cleanlowqual(void* threadId_p)
{
	int threadId = *((int*)threadId_p);
	
	for (uint64_t i=0; i<readNum; i++)
	{	if (i%threadNum == threadId)
		{	
			//cerr << StoreReads[i] << endl;
			//cerr << StoreQuals[i] << endl;
			if(StoreReads[i].size() != StoreQuals[i].size()){
				StoreReads[i] = "";
				StoreQuals[i] = "";
			}

			int seq_len = StoreReads[i].size();
			
			//cerr << "seq_len:" << seq_len << endl;
			//cerr << StoreReads[i] << endl;
			//cerr << StoreQuals[i] << endl;
			//calculate average error rate for the read
			double error_rate = 0;
			int N_count = 0;
			vector <int> breakpoints;
			
			for (int j=0; j<seq_len; j++)
			{	if (StoreReads[i][j] == 'N')
				{	N_count ++;
					StoreQuals[i][j] = Quality_shift;   
				}
				error_rate += Qual2Err[StoreQuals[i][j]];
			}
			
			StoreHeads[i] += "    RQ: " + boost::lexical_cast<std::string>(error_rate/seq_len*100) + "%";
			
			//cerr << "Error rate: " << error_rate << endl;
			
			//trim the low quality reads
			if (error_rate > Error_rate_cutoff * seq_len)
			{	
				double accum_error = 0;
				int accum_length = 0;
				int breakpos = 0;
				int breakpos_last = 0;
				int start_in_maxblock = 0;
				int end_in_maxblock = 0;
				int length_maxblock = 0;
				int start_in_block = 0;
				int end_in_block = 0;
				int length_block = 0;
				
				//get all blocks in a reads with average error rate lower than the cutoff, and find out the longest block as the final output result
				for (int j=0; j<seq_len; j++)
				{	
					accum_error += Qual2Err[StoreQuals[i][j]];
					accum_length ++;
					if (accum_error > Error_rate_cutoff * accum_length) //find a breakpoint
					{	
						breakpos = j + 1;
						start_in_block = breakpos_last + 1;
						end_in_block = breakpos - 1;
						length_block = end_in_block - start_in_block + 1;
						//cerr << start_in_block << "\t" << end_in_block << "\t" << length_block << endl;
						if (length_block > length_maxblock)  //replace the max record
						{	length_maxblock = length_block;
							start_in_maxblock = start_in_block;
							end_in_maxblock = end_in_block;
						}
						accum_error = 0;
						accum_length = 0;
						breakpos_last = j + 1;
					}
				}
				//cerr << "before" << endl;
				//get the last block
				breakpos = seq_len + 1;
				start_in_block = breakpos_last + 1;
				end_in_block = breakpos - 1;
				length_block = end_in_block - start_in_block + 1;
				//cerr << start_in_block << "\t" << end_in_block << "\t" << length_block << endl;
				if (length_block > length_maxblock)
				{	length_maxblock = length_block;
					start_in_maxblock = start_in_block;
					end_in_maxblock = end_in_block;
				}
				//cerr << "maxblock: " << start_in_maxblock << "\t" << end_in_maxblock << "\t" << length_maxblock << endl;
				
				StoreHeads[i] += "  TrimLowQual";
				if(start_in_maxblock >= 1 && start_in_maxblock <= seq_len)
				{
					StoreReads[i] = StoreReads[i].substr(start_in_maxblock-1, length_maxblock);
					StoreQuals[i] = StoreQuals[i].substr(start_in_maxblock-1, length_maxblock);
				}else{
					StoreReads[i] = "";
					StoreQuals[i] = "";
					length_maxblock = 0;
				}		
				__sync_add_and_fetch(&total_filtered_lowqual_reads,1);  //原子计算，多线程可同时操作内存变量
				__sync_add_and_fetch(&total_filtered_lowqual_bases, seq_len - length_maxblock);  //原子计算，多线程可同时操作内存变量
			}
			
			//cerr << "trimmed low qulaity" << endl;

			//filter short reads after trimming
			if (StoreReads[i].size() < Min_read_len)
			{	
				__sync_add_and_fetch(&total_filtered_short_reads,1);  //原子计算，多线程可同时操作内存变量
				__sync_add_and_fetch(&total_filtered_short_bases,StoreReads[i].size());  //原子计算，多线程可同时操作内存变量
				StoreHeads[i] += "  FilterShort";
				StoreReads[i] = "";
				StoreQuals[i] = "";
				seq_len = 0;
			}
			//cerr << "filter short" << endl;			
			//make statistics of resulting clean reads
			if (StoreReads[i].size())
			{
				__sync_add_and_fetch(&total_clean_reads,1);   //原子计算，多线程可同时操作内存变量
				__sync_add_and_fetch(&total_clean_bases,StoreReads[i].size());  //原子计算，多线程可同时操作内存变量
			}
			//cerr << "stat done\n" << endl;
			
		}
	}
}


int main(int argc, char *argv[])
{	
	//get options from command line
	int c;
	while((c=getopt(argc, argv, "e:r:q:t:h")) !=-1) {
		switch(c) {
			case 'e': Error_rate_cutoff=atof(optarg); break;
			case 'r': Min_read_len=atoi(optarg); break;
			case 'q': Quality_shift=atoi(optarg); break;
			case 't': threadNum=atoi(optarg); break;
			case 'h': usage(); break;
			default: usage();
		}
	}
	if (argc < 4) usage();
	
	string in_reads1_file = argv[optind++]; //optind, argv[optind++]顺序指向非option的参数
	string out_reads1_file = argv[optind++];
	string out_stat_file = argv[optind++];

	clock_t time_start, time_end;
	time_start = clock();
	
	time_end = clock();
	cerr << "\nProgram starting\n";
	cerr << "Run time: " << double(time_end - time_start) / CLOCKS_PER_SEC << endl;
	

	//assign initial values to Qual2Err
	for(int i=0; i<100; i++)
	{	Qual2Err[i+Quality_shift] = pow(10.0,-i/10.0);
	}
	
	igzstream infile1 (in_reads1_file.c_str());
	ogzstream cleanfile1 (out_reads1_file.c_str());
	ofstream statfile (out_stat_file.c_str());
	
	StoreReads = new string[bufferNum];
	StoreHeads  = new string[bufferNum];
	StoreQuals  = new string[bufferNum];
		
	total_raw_reads = 0;
	total_raw_bases = 0;
	total_filtered_lowqual_reads = 0;
	total_filtered_lowqual_bases = 0;
	total_filtered_short_reads = 0;
	total_filtered_short_bases = 0;
	total_clean_reads = 0;
	total_clean_bases = 0;

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

		//定义children thread数组
		pthread_t *pthread = new pthread_t[threadNum];
		int *pthreadId = new int[threadNum];
		
		if (pthread == NULL || pthreadId == NULL)
		{	cerr << "Out of memory!" << endl;
			exit(0);
		}
		
		/* create children threads */
		for (int i=0; i<threadNum; i++)
		{ 
			pthreadId[i] = i;
			pthread_create(pthread+i, NULL, thread_cleanlowqual, (void*)(pthreadId+i));
		}
		cerr << threadNum << " threads were created!" << endl;

		

		/* wait threads to exit and free resource after threads exit */
		for (int i=0; i<threadNum; i++)
		{ 
			pthread_join(pthread[i], NULL);
		}
		
		cerr << "All threads of this batch finished" << endl;
		
		
		//动态分配的数组需要用delete主动删除，否则内存可能不会释放
		delete[] pthread;
		delete[] pthreadId;

		//output the results
		for (int i = 0; i < readNum; i++)
		{	cleanfile1 << StoreHeads[i]   << "\n" << StoreReads[i] << "\n+\n" <<  StoreQuals[i] << endl;
		}
		
		//如果已经到了文件尾，退出读文件的循环
		if (readNum < bufferNum)
		{	break; //读到文件尾退出
		}
	}

	//动态分配的数组需要用delete主动删除，否则内存可能不会释放
	delete[] StoreReads;
	delete[] StoreHeads;
	delete[] StoreQuals;

	//output the statistic result

	statfile << "#total_raw_reads:   " << total_raw_reads << endl;
	statfile << "#total_raw_bases:   " << total_raw_bases << endl;


	statfile << "#filtered_lowqual_reads: " << total_filtered_lowqual_reads << "\t" << (double)total_filtered_lowqual_reads/total_raw_reads*100 << "%"  << endl;
	statfile << "#filtered_lowqual_bases: " << total_filtered_lowqual_bases << "\t" << (double)total_filtered_lowqual_bases/total_raw_bases*100 << "%"  << endl;
	
	statfile << "#filtered_short_reads: " << total_filtered_short_reads << "\t" << (double)total_filtered_short_reads/total_raw_reads*100 << "%"  << endl;
	statfile << "#filtered_short_bases: " << total_filtered_short_bases << "\t" << (double)total_filtered_short_bases/total_raw_bases*100 << "%"  << endl;
	
	statfile << "#total_clean_reads: " << total_clean_reads << "\t" << total_clean_reads/(double)total_raw_reads*100 << "%" << endl;
	statfile << "#total_clean_bases: " << total_clean_bases << "\t" << total_clean_bases/(double)total_raw_bases*100 << "%" << endl;
	

	time_end = clock();
	cerr << "\nAll jobs finishd done\n";
	cerr << "Run time: " << double(time_end - time_start) / CLOCKS_PER_SEC << endl;

}

