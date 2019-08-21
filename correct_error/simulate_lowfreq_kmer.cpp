//研究kmer大小，对于鉴别error Kmer的影响，多大比例为低频Kmer(非real Kmer),来帮助决定对于特定基因组选取合适的Kmer大小
//将参考基因组中的所有kmer装入内存中，然后再基因组上每隔一定距离(如100bp)模拟一个测序错误(碱基突变)，看
//其所影响的周边Kmer中低频Kmer(在参考Kmer集之外)的比例，以及统计低频比例的分布。

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
#include "gzstream.h"

using namespace std;

uint8_t *construct_ref_kmer_table(string &genome_seq_file, int KmerSize, uint64_t &total);
void stat_mutated_lowfreq_kmers(string &genome_seq_file, uint8_t *RefKmer, int KmerSize, int Skip_dist);
void mutate_one_base(string &frag_seq, int mutate_site);

//set global variables
int KmerSize = 17;            //当KmerSize=17时，4^17 = 2^34 = 16G, store 1bit for a kmer, kmer频率表将占内存2G
int Skip_dist = 100;

void usage() 
{	cout << "\nsimulate_lowfreq_kmer <genome_seq*.fa | *.fa.gz>\n" 
		 << "   -k <int>   set kmer size, default=" << KmerSize << endl
		 << "   -s <int>   set skip distance of muation on reference, default=" << Skip_dist << endl
		 << "   -h         get the help information\n" << endl
		 << "\n\nInstructions: We simulate each mutation (as sequencing error) along the reference genome with distance 100 bp, one mutation will affect a set of Kmers crossed that site, the number of these Kmers is equal to the Kmer size. We defined low-frequency as a Kmer do not exist in the reference genome, then we calculate the ratio of low-frequency Kmers in each set, and make distribution statistics, to illustrate the difficulty of error correction for each type of genomes.\n" << endl;
	exit(0);
}


int main(int argc, char *argv[])
{	
	//get options from command line
	int c;
	while((c=getopt(argc, argv, "k:s:h")) !=-1) {
		switch(c) {
			case 'k': KmerSize=atoi(optarg); break;
			case 's': Skip_dist=atoi(optarg); break;
			case 'h': usage(); break;
			default: usage();
		}
	}
	if (argc < 2) usage();
	
	clock_t time_start, time_end;
	time_start = clock();
	
	string genome_seq_file = argv[optind++]; //optind, argv[optind++]顺序指向非option的参数	
	cerr << "\nInput file is: " << genome_seq_file << "\n\n" << endl;
	
	cerr << "Begin to construct the reference kmer table:" << endl;
	uint64_t total=0;
	uint8_t *RefKmer = construct_ref_kmer_table(genome_seq_file, KmerSize, total);
	time_end = clock();
	cerr << "\nFinished time: " << double(time_end - time_start) / CLOCKS_PER_SEC << endl;

	cerr << "\nBegin to analyze the mutated Kmers:" << endl;
	stat_mutated_lowfreq_kmers(genome_seq_file, RefKmer, KmerSize, Skip_dist);
	time_end = clock();
	cerr << "\nFinished time: " << double(time_end - time_start) / CLOCKS_PER_SEC << endl;

}

void stat_mutated_lowfreq_kmers(string &genome_seq_file, uint8_t *RefKmer, int KmerSize, int Skip_dist)
{
	//set variables and do initialization
	int Frag_len = KmerSize + KmerSize - 1;
	int mutate_site = KmerSize - 1;
	uint64_t *freq_nums = new uint64_t[KmerSize+1];
	for (int i=0; i<=KmerSize; i++)
	{	freq_nums[i] = 0;
	}
	
	igzstream infile;
	infile.open(genome_seq_file.c_str());
	if ( ! infile )
	{	cerr << "fail to open input file" << genome_seq_file << endl;
	}
	
	//construct reference kmer table, parse one chromosome at a time
	string textline;
	getline( infile, textline, '>' );
	while ( getline( infile, textline, '\n' ) )
	{	
		//cerr << "parsing  " << textline << endl;

		//读取一条染色体的序列
		getline( infile, textline, '>' );
		string seq;
		for (int i=0; i<textline.size(); i++)
		{	if (textline[i] != '\n' && textline[i] != ' ')
			{	
				seq.push_back( textline[i] );
			} 
		}
		
		//generate muated kmers and analyze frequency in reference kmer table
		for (int i = 0; i <= seq.size()- Frag_len; i += Skip_dist)
		{	string frag_seq = seq.substr(i, Frag_len);
			mutate_one_base(frag_seq, mutate_site);
			uint8_t num_lowfreq_kmers = 0;
			for (int j = 0; j < KmerSize; j ++)
			{	string kmer_seq = frag_seq.substr(j, KmerSize);
				uint64_t kmer_bit = seq2bit(kmer_seq);
				uint8_t kmer_freq = get_freq(RefKmer, kmer_bit);
				if (kmer_freq == 0)
				{	num_lowfreq_kmers ++;
				}
			}
			freq_nums[num_lowfreq_kmers] ++;
		}
	}

	infile.close();

	//make final statisitcs
	uint64_t total_low_kmers = 0;
	uint64_t total_group_kmers = 0;
	uint64_t total_group100_kmers = 0;
	uint64_t total_group80_kmers = 0;
	uint64_t total_group50_kmers = 0;
	uint64_t total_group20_kmers = 0;
	uint64_t total_group1_kmers = 0;
	double low_stat_ratio = 0.0;
	double group100_stat_ratio = 0.0;
	double group80_stat_ratio = 0.0;
	double group50_stat_ratio = 0.0;
	double group20_stat_ratio = 0.0;
	double group1_stat_ratio = 0.0;

	for (int i = 0; i <= KmerSize; i ++)
	{	
		total_group_kmers += freq_nums[i];
		total_low_kmers += i * freq_nums[i];
		
		double group_low_ratio = (double)i / KmerSize;
		//cerr << "## " << group_low_ratio << endl;
		if ( i == KmerSize)
		{	total_group100_kmers += freq_nums[i];
		}
		if ( group_low_ratio >= 0.8  )
		{	total_group80_kmers += freq_nums[i];
		}
		if ( group_low_ratio >= 0.5  )
		{	total_group50_kmers += freq_nums[i];
		}
		if ( group_low_ratio >= 0.2  )
		{	total_group20_kmers += freq_nums[i];
		}
		if ( i >= 1)
		{	total_group1_kmers += freq_nums[i];
		}
	}
	group100_stat_ratio = (double)total_group100_kmers / total_group_kmers;
	group1_stat_ratio = (double)total_group1_kmers / total_group_kmers;
	group80_stat_ratio = (double)total_group80_kmers / total_group_kmers;
	group50_stat_ratio = (double)total_group50_kmers / total_group_kmers;
	group20_stat_ratio = (double)total_group20_kmers / total_group_kmers;
	low_stat_ratio = (double)total_low_kmers / (total_group_kmers*KmerSize);

	cout << "\nKmer size: " << KmerSize << endl;
	cout << "\nRatio of low-freq kmers in all kmers by muation : " << low_stat_ratio << endl;
	cout << "\nRatio of mutations with 100% low-freq kmers:  " << group100_stat_ratio << endl;
	cout << "\nRatio of mutations with >=80% low-freq kmers: " << group80_stat_ratio << endl;
	cout << "\nRatio of mutations with >=50% low-freq kmers: " << group50_stat_ratio << endl;
	cout << "\nRatio of mutations with >=20% low-freq kmers: " << group20_stat_ratio << endl;
	cout << "\nRatio of mutations with >= 1 low-freq kmers:  " << group1_stat_ratio << endl;
}


void mutate_one_base(string &frag_seq, int mutate_site)
{
	char raw_base = frag_seq[mutate_site];
	uint8_t new_base_num = alphabet[raw_base] + 1;
	if (new_base_num == 4)
	{	new_base_num = 0;
	}
	char new_base = bases[new_base_num];
	frag_seq[mutate_site] = new_base;
}

uint8_t *construct_ref_kmer_table(string &genome_seq_file, int KmerSize, uint64_t &total)
{
	//set variables and initialize values
	total=0;  
	for(int i=0; i<KmerSize; i++) {
		total=(total<<2)|0x3;
	}
	
	//动态分配的数组必须赋初值
	uint64_t array_size = total/8+1;
	uint8_t *freq=new uint8_t[array_size]; 
	for(uint64_t i=0; i<array_size; i++) {
		freq[i]=0; 
	}
	
	igzstream infile;
	infile.open(genome_seq_file.c_str());
	if ( ! infile )
	{	cerr << "fail to open input file" << genome_seq_file << endl;
	}
	
	//construct reference kmer table, parse one chromosome at a time
	uint64_t genome_size_bp = 0;
	uint64_t kmer_total_number = 0;
	uint64_t kmer_species_number = 0;

	string textline;
	getline( infile, textline, '>' );
	while ( getline( infile, textline, '\n' ) )
	{	
		//cerr << "parsing  " << textline << endl;

		//读取一条染色体的序列
		getline( infile, textline, '>' );
		string seq;
		for (int i=0; i<textline.size(); i++)
		{	
			if (textline[i] != '\n' && textline[i] != ' ')
			{	
				seq.push_back( textline[i] );
			} 
		}
		
		genome_size_bp += seq.size();
		kmer_total_number += seq.size()-KmerSize+1;

		//读取一条染色体上的kmer, 1-bit stored in freq table
		for (int i = 0; i<=seq.size()-KmerSize; i++)
		{	string kseq = seq.substr(i, KmerSize);
			uint64_t idx = seq2bit(kseq);
			freq[idx/8] |= bitAll[idx%8];

			uint64_t rc_idx = get_rev_com_kbit(idx, KmerSize);
			freq[rc_idx/8] |= bitAll[rc_idx%8];
		}
		
	}
	infile.close();

	for (uint64_t idx=0; idx<total; idx++)
	{	if (get_freq(freq, idx))
		{	kmer_species_number ++;
		}
	}
	
	//output the statistics result
	cout << "The Genome size is:  " << genome_size_bp << endl;
	cout << "Kmer total number:   " << kmer_total_number << endl;
	cout << "Kmer species number: " << kmer_species_number << endl << endl;

	return freq;
}

