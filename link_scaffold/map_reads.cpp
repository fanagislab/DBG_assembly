/*
This function maps single reads onto the contig sequences.

Author: Fan Wei, fanweiagis@126.com
Date: 2015-12-12
*/


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
#include "gzstream.h"
#include "seqKmer.h"
#include "kmerSet.h"
#include "map_func.h"


//global variables
KmerSet *kset;
vector <string> contig_ids;
vector <string> contig_seqs;

//variables and routines in parse the read files parallely
int threadNum = 10;
int BufferNum = 10000;  //number of reads stored each time in the memeory buffer 
string *ReadsID;  //store reads id
string *RawReads;  //store reads sequence
uint8_t *Signal; //store the status: 0, empty; 1, filled; 2, end of file;
int ReadsNum = 0;   //real number of reads stored in the buffer block
void *thread_parseBlock(void* threadId_p);
void parse_one_reads_file(string &reads_file);


//statistics for the alignment result, set as global variables because there are too many variables difficult to pass to a routine with local variables
int *contig_id_index ;
int *seed_contig_start ;
int *seed_contig_end ;
int *seed_read_start ;
int *seed_read_end ;
char *align_direct ;
		
int *align_contig_start ;
int *align_contig_end ;
int *align_read_start ;
int *align_read_end ;
float *align_identity ;

int *contig2_id_index ;
int *seed2_contig_start ;
int *seed2_contig_end ;
int *seed2_read_start ;
int *seed2_read_end ;
char *align2_direct ;

int *align2_contig_start ;
int *align2_contig_end ;
int *align2_read_start ;
int *align2_read_end ;
float *align2_identity ;


void usage()
{
	cout << "\nFunction instruction:\
\n\nmap_reads, maps single illumina reads onto the contig sequences, using similar alignment method with map_pair. \
One read could at most be mapped to two contigs, and only those reads mapped to two different contigs could be used \
to link contigs and fill the gaps between contigs. Tt is recommended to set parameter MinCtgLen (-l) to be 1/2 * MinReadLen (-r). \ 
The output of map_reads is the input of link_scafftig. \n";

	cout << "\nmap_reads  <contig_file.fa>  <reads_files.lib>\n"
		 << "   Function:  map single reads onto contigs"  << endl 
		 << "   Version: 1.0"  << endl
		 << "   -k <int>     kmer size (construct hash), default=" << KmerSize << endl
	 	 << "   -s <int>     seed size (number of contained kmers in a seed), default=" << SeedKmerNum << endl
		 << "   -l <int>     contigs not shorter than this cutoff are used for mapping [and scaffolding], default=" << MinCtgLen << endl
		 << "   -r <int>     reads not shorter than this cutoff are used for mapping [and scaffolding], default=" << MinReadLen << endl
		 << "   -i <float>   minimum mapping identity, default=" << MinMapIdentity << endl
	     << "   -f <int>     input file format: 1: fq|gz(one-line), 2: fa|gz(one-line), default=" << Input_file_format << endl
		 << "   -o <str>     output direcotry, default = " << Output_prefix << endl
		 << "   -t <int>     number of threads to run in parallel, default=" << threadNum << endl
		 << "   -h           get the help information\n" << endl
		 << "Example: map_reads  -l 125 -r 250 -t 10 -o ./maping_results/  Ecoli.contig.fa illumina_reads.lib" << endl
		 << endl;
	exit(0);
}


int main(int argc, char *argv[])
{	
	//get options from command line
	int c;
	while((c=getopt(argc, argv, "k:s:l:r:i:f:o:t:h")) !=-1) {
		switch(c) {
			case 'k': KmerSize=atoi(optarg); break;
			case 's': SeedKmerNum=atoi(optarg); break;
			case 'l': MinCtgLen=atoi(optarg); break;
			case 'r': MinReadLen=atoi(optarg); break;
			case 'i': MinMapIdentity=atof(optarg); break;
			case 'f': Input_file_format=atoi(optarg); break;
			case 'o': Output_prefix=optarg; break;
			case 't': threadNum=atoi(optarg); break;
			case 'h': usage(); break;
			default: usage();
		}
	}
	
	if (argc < 3) usage();
	
	string contig_seq_file = argv[optind++];
	string reads_lib_file = argv[optind++];
	
	clock_t time_start, time_end;
	time_start = clock();

	cerr << "\nProgram start ............" << endl;
	time_end = clock();
	cerr << "Run time: " << double(time_end - time_start) / CLOCKS_PER_SEC << endl;
	
	
	vector <string> reads_files;
	reading_lib_file(reads_lib_file, reads_files);
	
	cerr << "\nInput reads file number: " << reads_files.size() << endl;
	for ( int i = 0; i < reads_files.size(); i++)
	{	cerr << reads_files[i] << endl;
	}
	
	cerr << "\nLoad the reads lib file finished" << endl;
	time_end = clock();
	cerr << "Run time: " << double(time_end - time_start) / CLOCKS_PER_SEC << endl;

	cerr << "\nCreat the mapped 2ctg libary file" << endl;
	string mapped_2ctg_file = reads_lib_file + ".map_reads.2ctg.lib";
	ofstream TwoCtgFile (mapped_2ctg_file.c_str());
	for ( int i = 0; i < reads_files.size(); i++)
	{	TwoCtgFile << Output_prefix << "/" << reads_files[i] << ".map_reads.2ctg.gz" << endl;
	}

	read_contig_file(contig_seq_file, contig_ids, contig_seqs);
	uint64_t total_contig_num = 0;
	uint64_t total_contig_len = 0;
	for ( int i = 0; i < contig_seqs.size(); i ++)
	{	if (contig_seqs[i].size() >= MinCtgLen)
		{
			total_contig_num ++;
			total_contig_len += contig_seqs[i].size();
		}else
		{	contig_seqs[i] = "";    //remove the contig sequence if its length is shorter than the MinCtgLen cutoff
		}

	}
	cerr << "\nInput contig sequence number: " << total_contig_num << endl;
	cerr << "Total contig sequence length: " << total_contig_len << endl;
	cerr << "Load contigs into memory finished" << endl;
	time_end = clock();
	cerr << "Run time: " << double(time_end - time_start) / CLOCKS_PER_SEC << endl;
	
	string mkdir_cmd = "mkdir " + Output_prefix;
	system(mkdir_cmd.c_str());

	uint64_t hash_size = total_contig_len * 3;
	double load_factor = 0.5;
	kset = init_kmerset(hash_size, load_factor);
	chop_contig_to_kmerset(kset, contig_seqs);
	
	cerr << "\nThe hash initialization space(array) size: " << hash_size << endl;
	cerr << "The hash loading factor :   " << load_factor << endl;
	cerr << "Build the kmer hash finished" << endl;
	time_end = clock();
	cerr << "Run time: " << double(time_end - time_start) / CLOCKS_PER_SEC << endl;

	print_kmerset_parameter(kset);

	
	cerr << "\nAlign input reads to the kmer-hash: " << endl; 
	for (int i = 0; i < reads_files.size(); i++)
	{	
		cerr << "\n\tParse " << reads_files[i] << " ............." << endl;
		parse_one_reads_file(reads_files[i]);
	}

	cerr << "\nProgram finished !" << endl;
	time_end = clock();
	cerr << "Run time: " << double(time_end - time_start) / CLOCKS_PER_SEC << endl;
	
}



//parse one reads file, invoked by the main function
void parse_one_reads_file(string &reads_file)
{
		
	int pos = reads_file.find_last_of('/');
	string reads_file_name = reads_file.substr(pos+1, reads_file.size()-pos-1);

	//creat result files
	string reads_map_ctg_diff_file = Output_prefix + "/" + reads_file_name + ".map_reads.2ctg.gz";
	ogzstream MapCtgDiff ( reads_map_ctg_diff_file.c_str() );
	if ( ! MapCtgDiff )
	{       cerr << "fail to open output file " << reads_map_ctg_diff_file << endl;
	}

	string reads_map_ctg_diff_seq_file = Output_prefix + "/" + reads_file_name + ".map_reads.2ctg.gz.reads.fa.gz";
	ogzstream MapCtgDiffSeq ( reads_map_ctg_diff_seq_file.c_str() );
	if ( ! MapCtgDiffSeq )
	{       cerr << "fail to open output file " << reads_map_ctg_diff_seq_file << endl;
	}
	
	string reads_map_ctg_same_file = Output_prefix + "/" + reads_file_name + ".map_reads.1ctg.gz";
	ogzstream MapCtgSame ( reads_map_ctg_same_file.c_str() );
	if ( ! MapCtgSame )
	{       cerr << "fail to open output file " << reads_map_ctg_same_file << endl;
	}
	

	string reads_map_ctg_stat_file = Output_prefix + "/" + reads_file_name + ".map_reads.stat";
	ofstream MapCtgStat ( reads_map_ctg_stat_file.c_str() );
	if ( ! MapCtgStat )
	{       cerr << "fail to open output file " << reads_map_ctg_stat_file << endl;
	}
	
	MapCtgDiff << "#read_id\tread_length\talign_read_start\talign_read_end\tcontig_id\tcontig_length\talign_contig_start\talign_contig_end\talign_direct\talign_identity%\tread_id\tread_length\talign2_read_start\talign2_read_end\tcontig2_id\tcontig2_length\talign2_contig_start\talign2_contig_end\talign2_direct\talign2_identity%" << endl;

	MapCtgSame << "#read_id\tread_length\talign_read_start\talign_read_end\tcontig_id\tcontig_length\talign_contig_start\talign_contig_end\talign_direct\talign_identity%" << endl;

	uint64_t total_read_num = 0;
	uint64_t map_ctg_diff_num = 0;
	uint64_t map_ctg_same_num = 0;
	uint64_t map_no_no_num = 0;
	uint64_t error_map_num = 0;
	
	string LineStr;
	igzstream currentFile; //Caution: gzstream.cpp has a bug, which can be overcomed by gcc/g++ optimation (-O1, -O2, -O3).
	currentFile.open(reads_file.c_str());
	

	//allocate memory
	ReadsID = new string[BufferNum];
	RawReads = new string[BufferNum];
	Signal = new uint8_t[BufferNum];

	contig_id_index = new int[BufferNum];
	seed_contig_start  = new int[BufferNum];
	seed_contig_end  = new int[BufferNum];;
	seed_read_start  = new int[BufferNum];
	seed_read_end  = new int[BufferNum];
	align_direct = new char[BufferNum];
		
	align_contig_start  = new int[BufferNum];
	align_contig_end  = new int[BufferNum];
	align_read_start  = new int[BufferNum];
	align_read_end  = new int[BufferNum];
	align_identity = new float[BufferNum];

	contig2_id_index  = new int[BufferNum];
	seed2_contig_start  = new int[BufferNum];
	seed2_contig_end  = new int[BufferNum];
	seed2_read_start  = new int[BufferNum];
	seed2_read_end  = new int[BufferNum];
	align2_direct  = new char[BufferNum];

	align2_contig_start  = new int[BufferNum];
	align2_contig_end  = new int[BufferNum];
	align2_read_start  = new int[BufferNum];
	align2_read_end  = new int[BufferNum];
	align2_identity = new float[BufferNum];
	
	//load a block of reads into memory and then parse them block by block
	while (1)
	{	
		//asign initial values
		ReadsNum = 0;
		memset(Signal, 0, BufferNum);

		//create and run children threads to parse the reads that have loaded into the memory buffer
		pthread_t *pthread = new pthread_t[threadNum];
		int *pthreadId = new int[threadNum];
		for (int i=0; i<threadNum; i++)
		{
			pthreadId[i] = i;
			pthread_create((pthread+i), NULL, thread_parseBlock, (void*)(pthreadId+i));

		}
		//cerr << threadNum << " threads creation done!\n" << endl;
		
		//load the reads data into memory buffer by the father thread
		if (Input_file_format == 1)
		{	//读取fq格式文件, support one-line fastq format
			while ( ReadsNum < BufferNum && getline( currentFile, ReadsID[ReadsNum], '\n' ) )
			{	
				if (ReadsID[ReadsNum][0] == '@') 
				{	
					getline( currentFile, RawReads[ReadsNum], '\n');			
					getline( currentFile, LineStr, '\n');
					getline( currentFile, LineStr, '\n');
					Signal[ReadsNum] = 1;
					ReadsNum ++;
				}	
			}
		}
		else
		{	//读取fa格式文件, support one-line fasta format
			while ( ReadsNum < BufferNum && getline( currentFile, ReadsID[ReadsNum], '\n' )  )
			{	
				if (ReadsID[ReadsNum][0] == '>') 
				{
					getline( currentFile, RawReads[ReadsNum], '\n');	
					Signal[ReadsNum] = 1;	
					ReadsNum ++;
				}	
			}
		}
		
		//cerr << "load " << ReadsNum << " reads into memory done\n" << endl;
		
		//judge the end of file, and make signal for the children threads
		if ( ReadsNum < BufferNum )
		{
			for (int i = ReadsNum; i < BufferNum; i ++) 
			{	
				Signal[i] = 2;
			}
		}

		//当父线程读取数据结束后，等待全部子线程分析数据结束
		for (int i=0; i<threadNum; i++)
		{
			pthread_join(pthread[i], NULL);
		}
		
		
		//cerr << "align these reads to contigs done\n" << endl;

		//output the result
		for (int i = 0; i < ReadsNum; i ++ )
		{	
		//////////////////////////////////////////compare the align result of left-end and right-end and make decision
			if (RawReads[i].size() < MinReadLen)
			{	continue;
			}
			
			total_read_num ++;

			if (contig_id_index[i] != -1 )
			{	
				if (contig2_id_index[i] != -1 )
				{	
					if (contig_id_index[i] != contig2_id_index[i])
					{
						map_ctg_diff_num ++;
						MapCtgDiff << ReadsID[i] << "\t" << RawReads[i].size() << "\t" << align_read_start[i] << "\t" << align_read_end[i] << "\t" 
						<< contig_ids[contig_id_index[i]] << "\t" << contig_seqs[contig_id_index[i]].size() << "\t" << align_contig_start[i] << "\t" << align_contig_end[i] << "\t" << align_direct[i] << "\t" << align_identity[i]*100 << "%\t" << ReadsID[i] << "\t" << RawReads[i].size() 
						<< "\t" << align2_read_start[i] << "\t" << align2_read_end[i] << "\t" 
						<< contig_ids[contig2_id_index[i]] << "\t"<< contig_seqs[contig2_id_index[i]].size() << "\t" << align2_contig_start[i] << "\t" << align2_contig_end[i] << "\t" << align2_direct[i] << "\t" << align2_identity[i]*100 << "%" << endl;
					
						MapCtgDiffSeq << ">" << ReadsID[i]  << "\n" << RawReads[i] << "\n";

					}else
					{	error_map_num ++;
					}

				}else
				{	
					map_ctg_same_num ++;
					MapCtgSame << ReadsID[i] << "\t" << RawReads[i].size() << "\t" << align_read_start[i] << "\t" << align_read_end[i] << "\t" 
					<< contig_ids[contig_id_index[i]] << "\t" << contig_seqs[contig_id_index[i]].size() << "\t"<< align_contig_start[i] << "\t" << align_contig_end[i] << "\t" << align_direct[i] << "\t" << align_identity[i]*100 << "%" << endl;

				}
			}else
			{	map_no_no_num ++;
			}

		//////////////////////////////////////////
			
		}
		
		//cerr << "output the aligning results to file done\n" << endl;
		

		//when reach the end of file, and all children threads finished
		if ( ReadsNum < BufferNum )
		{	break;
		}
		
	}

	currentFile.close();
	
	MapCtgStat << "\ttotal_read_num: " << total_read_num << endl;
	MapCtgStat << "\tmap_ctg_diff_num: " << map_ctg_diff_num << "  " <<  (double)map_ctg_diff_num/total_read_num*100 << "%" << endl;
	MapCtgStat << "\tmap_ctg_same_num: " << map_ctg_same_num << "  " << (double)map_ctg_same_num/total_read_num*100 << "%" << endl;
	MapCtgStat << "\tmap_no_no_num: " << map_no_no_num << "  " << (double)map_no_no_num/total_read_num*100 << "%" << endl;
	MapCtgStat << "\terror_map_num: " << error_map_num << "  " << (double)error_map_num/total_read_num*100 << "%" << endl;

}



//this is the thread routine to parse files parallely
void *thread_parseBlock(void* threadId_p)
{
	int threadId = *((int*)threadId_p);

	for (uint64_t i = threadId; i < BufferNum; i += threadNum)
	{		
		while (1)
		{
			if (Signal[i] == 1)
			{  
				//get the kmers from a read sequence
				string &read = RawReads[i];
				
				//make validate read_id
				vector<string> vec_head;
				split(ReadsID[i], vec_head, ">@ \t\n");
				ReadsID[i] = vec_head[0];
				if (vec_head.size() > 1)
				{	ReadsID[i] += "-" + vec_head[1];
				}
				
				//asign initial value
				contig_id_index[i] = -1;
				seed_contig_start[i] = -1;
				seed_contig_end[i] = -1;
				seed_read_start[i] = -1;
				seed_read_end[i] = -1;
				align_direct[i] = 'N';
			
				align_contig_start[i] = -1;
				align_contig_end[i] = -1;
				align_read_start[i] = -1;
				align_read_end[i] = -1;
				align_identity[i] = 0.0;

				contig2_id_index[i] = -1;
				seed2_contig_start[i] = -1;
				seed2_contig_end[i] = -1;
				seed2_read_start[i] = -1;
				seed2_read_end[i] = -1;
				align2_direct[i] = 'N';

				align2_contig_start[i] = -1;
				align2_contig_end[i] = -1;
				align2_read_start[i] = -1;
				align2_read_end[i] = -1;
				align2_identity[i] = 0.0;
				
				if (read.size() < MinReadLen)
				{	break;
				}

				/////////////////////////////////////////////////////////////////////////////map one reads

				//align the left part of a reads, seed and extension and calculate identity
				if(read.size() >= KmerSize+SeedKmerNum)
				{	get_align_seed(kset, read, 1, read.size(), contig_id_index[i], seed_contig_start[i], seed_contig_end[i], seed_read_start[i], seed_read_end[i], align_direct[i]);
					if (contig_id_index[i] != -1)
					{	align_contig_start[i] = seed_contig_start[i];
						align_contig_end[i] = seed_contig_end[i];
						align_read_start[i] = seed_read_start[i];
						align_read_end[i] = seed_read_end[i];

						extend_align_region(read, contig_seqs[contig_id_index[i]], align_identity[i], align_contig_start[i], align_contig_end[i], align_read_start[i], align_read_end[i], align_direct[i]);
						if (align_identity[i] < MinMapIdentity)
						{	contig_id_index[i] = -1;
						}
					}
				}


				//align the right part of a reads, seed and extension and calculate identity, when the left part of the reads has been positioned to a contig while its right part is not aligned
				if (contig_id_index[i] != -1 && align_read_end[i] < read.size())
				{
		
					if( (read.size() - align_read_end[i]) >= KmerSize+SeedKmerNum)
					{	get_align_seed(kset, read, align_read_end[i]+1, read.size(), contig2_id_index[i], seed2_contig_start[i], seed2_contig_end[i], seed2_read_start[i], seed2_read_end[i], align2_direct[i]);
						if (contig2_id_index[i] != -1)
						{	align2_contig_start[i] = seed2_contig_start[i];
							align2_contig_end[i] = seed2_contig_end[i];
							align2_read_start[i] = seed2_read_start[i];
							align2_read_end[i] = seed2_read_end[i];

							extend_align_region(read, contig_seqs[contig2_id_index[i]], align2_identity[i], align2_contig_start[i], align2_contig_end[i], align2_read_start[i], align2_read_end[i], align2_direct[i]);
							if (align2_identity[i] < MinMapIdentity)
							{	contig2_id_index[i] = -1;
							}
						}
					}

				}

				/////////////////////////////////////////////////////////////////////////////map one reads
				
				
				break;  //finished parsing kmers from one reads, and break loop to parse the next reads
		
			}else
			{
				if (Signal[i] == 0)
				{	usleep(1); // 相当于0.01s, 休息0.01秒，不占计算资源
				}else  
				{	break;   //when signal[i] equals 2, reach the end of file
				}
			}
		}


	}		
			
	return NULL;
}



