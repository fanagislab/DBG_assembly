/*
This function maps pairs of reads onto the contig sequences.

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
#include "seqKmer.h"
#include "kmerSet.h"
#include "map_func.h"
#include "gzstream.h"

using namespace std;

void parse_two_paired_reads_file(string &input_reads_file, string &input2_reads_file,  KmerSet *kset, vector <string> &contig_ids, vector <string> &contig_seqs);

void usage()
{   

	cout << "\nFunction instruction:\n\nmap_pair aligns pair of reads to the contigs or scafftigs, using a seed-and-extension globle alignment method. Note that the cutoff for contig size (MinCtgLen -l) is not critical, because link_scaffold will automatically determine the filtering cutoff of small contigs by the insert size of pair-end or mate-paired reads. It is OK to use all or the major part of the contigs for mapping. The recommeded settings is: MinCtgLen (-l) = 1/2 * MinReadLen (-r)\n";
		
	cout << "\nmap_pair  <contig_file.fa>  <reads_files.lib>\n" 
		<< "   Function:  map pair-end or mate-pair reads onto contigs"  << endl 
		 << "   Version: 1.0"  << endl
 
		 << "   -k <int>     kmer size (construct hash), default=" << KmerSize << endl
	 	 << "   -s <int>     seed size (number of contained kmers in a seed), default=" << SeedKmerNum << endl
		 << "   -l <int>     contigs not shorter than this cutoff are used for mapping [and scaffolding], default=" << MinCtgLen << endl
		 << "   -r <int>     reads not shorter than this cutoff are used for mapping [and scaffolding], default=" << MinReadLen << endl
		 << "   -i <float>   minimum mapping identity, default=" << MinMapIdentity << endl
	     << "   -f <int>     input file format: 1: fq|gz(one-line), 2: fa|gz(one-line), default=" << Input_file_format << endl
		 << "   -o <str>     output directory, default = " << Output_prefix << endl
		 << "   -h           get the help information\n" << endl
		 << "Example: map_pair  -l 125 -r 250 -o ./maping_results/  Ecoli.contig.fa illumina_reads.lib" << endl
		 << endl;

	exit(0);
}


int main(int argc, char *argv[])
{	
	//get options from command line
	int c;
	while((c=getopt(argc, argv, "k:s:l:r:i:f:o:h")) !=-1) {
		switch(c) {
			case 'k': KmerSize=atoi(optarg); break;
			case 's': SeedKmerNum=atoi(optarg); break;
			case 'l': MinCtgLen=atoi(optarg); break;
			case 'r': MinReadLen=atoi(optarg); break;
			case 'i': MinMapIdentity=atof(optarg); break;
			case 'f': Input_file_format=atoi(optarg); break;
			case 'o': Output_prefix=optarg; break;
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

	cerr << "\nCreat the mapped 2ctg libary file" << endl;
	string mapped_2ctg_file = reads_lib_file + ".map_pair.2ctg.lib";
	ofstream TwoCtgFile (mapped_2ctg_file.c_str());
	for ( int i = 0; i < reads_files.size(); i+=2)
	{	TwoCtgFile << Output_prefix << "/" << reads_files[i] << ".map_pair.2ctg.gz" << endl;
	}


	vector <string> contig_ids, contig_seqs;
	read_contig_file(contig_seq_file, contig_ids, contig_seqs);
	
	uint64_t total_contig_num = 0;
	uint64_t total_contig_len = 0;
	for ( int i = 0; i < contig_seqs.size(); i ++)
	{	
		if (contig_seqs[i].size() >= MinCtgLen)
		{	total_contig_num ++;
			total_contig_len += contig_seqs[i].size();
		}else
		{	contig_seqs[i] = "";    //remove the contig sequence if its length is shorter than the MinCtgLen cutoff
		}
	}
	cerr << "\nInput contig sequence number: " << total_contig_num << endl;
	cerr << "Total contig sequence length: " << total_contig_len << endl;

	cerr << "\nRead contigs into memory finished" << endl;
	time_end = clock();
	cerr << "Run time: " << double(time_end - time_start) / CLOCKS_PER_SEC << endl;
	
	string mkdir_cmd = "mkdir " + Output_prefix;
	system(mkdir_cmd.c_str());

	
	uint64_t hash_size = total_contig_len * 3;
	double load_factor = 0.5;
	KmerSet *kset = init_kmerset(hash_size, load_factor);	
	chop_contig_to_kmerset(kset, contig_seqs);
	
	cerr << "\nThe hash initialization space(array) size: " << hash_size << endl;
	cerr << "The hash loading factor :   " << load_factor << endl;
	cerr << "\nBuild the kmer hash finished" << endl;
	time_end = clock();
	cerr << "Run time: " << double(time_end - time_start) / CLOCKS_PER_SEC << endl;
	
	//read data into memory, construct kmer de bruijn graph
	cerr << "\nparse input reads files: " << endl; 
	for (int i = 0; i < reads_files.size(); i+=2)
	{	
		cerr << "\n\t" << reads_files[i] << " ............." << endl;
		cerr << "\n\t" << reads_files[i+1] << " ............." << endl;

		parse_two_paired_reads_file(reads_files[i], reads_files[i+1], kset, contig_ids, contig_seqs);
	}

	cerr << "\nProgram finished !" << endl;
	time_end = clock();
	cerr << "Run time: " << double(time_end - time_start) / CLOCKS_PER_SEC << endl;
	
}




void parse_two_paired_reads_file(string &input_reads_file, string &input2_reads_file,  KmerSet *kset, vector <string> &contig_ids, vector <string> &contig_seqs)
{
		
	igzstream infile ( input_reads_file.c_str() );
	if ( ! infile )
	{       cerr << "fail to open input file " << input_reads_file << endl;
	}

	igzstream infile2 ( input2_reads_file.c_str() );
	if ( ! infile2 )
	{       cerr << "fail to open input file " << input2_reads_file << endl;
	}
	


	int pos = input_reads_file.find_last_of('/');
	string reads_file_name = input_reads_file.substr(pos+1, input_reads_file.size()-pos-1);


	string reads_map_ctg_diff_file =  Output_prefix + "/" + reads_file_name + ".map_pair.2ctg.gz";
	ogzstream MapCtgDiff ( reads_map_ctg_diff_file.c_str() );
	if ( ! MapCtgDiff )
	{       cerr << "fail to open output file " << reads_map_ctg_diff_file << endl;
	}
	
	string reads_map_ctg_same_file =  Output_prefix + "/" + reads_file_name + ".map_pair.1ctg.gz";
	ogzstream MapCtgSame ( reads_map_ctg_same_file.c_str() );
	if ( ! MapCtgSame )
	{       cerr << "fail to open output file " << reads_map_ctg_same_file << endl;
	}
	
	string reads_map_ctg_gap_file =  Output_prefix + "/" + reads_file_name + ".map_pair.gap.gz";
	ogzstream MapCtgGap ( reads_map_ctg_gap_file.c_str() );
	if ( ! MapCtgGap )
	{       cerr << "fail to open output file " << reads_map_ctg_gap_file << endl;
	}
	
	string reads_map_ctg_stat_file =  Output_prefix + "/" + reads_file_name + ".map_pair.stat";
	ofstream MapCtgStat ( reads_map_ctg_stat_file.c_str() );
	if ( ! MapCtgStat )
	{       cerr << "fail to open output file " << reads_map_ctg_stat_file << endl;
	}

	uint64_t total_read_pair_num = 0;
	uint64_t map_ctg_diff_num = 0;
	uint64_t map_ctg_same_num = 0;
	uint64_t map_ctg_gap_num = 0;
	uint64_t map_no_no_num = 0;


	MapCtgDiff << "#read_id\tread_length\talign_read_start\talign_read_end\tcontig_id\tcontig_length\talign_contig_start\talign_contig_end\talign_direct\talign_identity%\tread_id\tread_length\talign2_read_start\talign2_read_end\tcontig2_id\tcontig2_length\talign2_contig_start\talign2_contig_end\talign2_direct\talign2_identity%" << endl;

	MapCtgSame << "#read_id\tread_length\talign_read_start\talign_read_end\tcontig_id\tcontig_length\talign_contig_start\talign_contig_end\talign_direct\talign_identity%" << endl;

	MapCtgGap << "#read_id\tread_length\talign_read_start\talign_read_end\tcontig_id\tcontig_length\talign_contig_start\talign_contig_end\talign_direct\talign_identity%" << endl;


	//parse each read sequence, fa or fq format
	string read, read_head, read_id;
	string read2, read2_head, read2_id;
	string line_no_use;
	while ( getline( infile, read_head, '\n' ) )
	{

		if (Input_file_format == 2 && read_head[0] == '>')   //fasta format 
		{  
			vector<string> vec_head;
			split(read_head, vec_head, "> \t");
			read_id = vec_head[0];
			if (vec_head.size() > 1)
			{	read_id	+= "-" + vec_head[1];
			}

			getline( infile, read, '\n' );
			
			getline( infile2, read2_head, '\n' );
			vector<string> vec2_head;
			split(read2_head, vec2_head, "> \t");
			read2_id = vec2_head[0];
			if (vec2_head.size() > 1)
			{	read2_id	+= "-" + vec2_head[1];
			}

			getline( infile2, read2, '\n' );

		}
		
		if (Input_file_format == 1 && read_head[0] == '@')   //fastq format 
		{  
			vector<string> vec_head;
			split(read_head, vec_head, "@ \t");
			read_id = vec_head[0];
			if (vec_head.size() > 1)
			{	read_id	+= "-" + vec_head[1];
			}

			getline( infile, read, '\n' );
			getline( infile, line_no_use, '\n' );
			getline( infile, line_no_use, '\n' );


			getline( infile2, read2_head, '\n' );
			vector<string> vec2_head;
			split(read2_head, vec2_head, "@ \t");
			read2_id = vec2_head[0];
			if (vec2_head.size() > 1)
			{	read2_id	+= "-" + vec2_head[1];
			}

			getline( infile2, read2, '\n' );
			getline( infile2, line_no_use, '\n' );
			getline( infile2, line_no_use, '\n' );

			
		}
		
		//cerr << read.size()  << "\t" << read2.size() << "\t" <<  MinReadLen << endl;
		//require both reads length pass a given cutoff
		if (read.size() < MinReadLen || read2.size() < MinReadLen)
		{	continue;
		}
		
		
		int contig_id_index = -1;
		int seed_contig_start = -1;
		int seed_contig_end = -1;
		int seed_read_start = -1;
		int seed_read_end = -1;
		char seed_direct = 'N';
		float align_identity = 0.0;
		
		//seed and extension and calculate identity
		if(read.size() >= KmerSize+SeedKmerNum)
		{	get_align_seed(kset, read, 1, read.size(), contig_id_index, seed_contig_start, seed_contig_end, seed_read_start, seed_read_end, seed_direct);
			if (contig_id_index != -1)
			{	extend_align_region(read, contig_seqs[contig_id_index], align_identity, seed_contig_start, seed_contig_end, seed_read_start, seed_read_end, seed_direct);
				if (align_identity < MinMapIdentity)
				{	contig_id_index = -1;
				}
			}
		}
		
		int contig2_id_index = -1;
		int seed2_contig_start = -1;
		int seed2_contig_end = -1;
		int seed2_read_start = -1;
		int seed2_read_end = -1;
		char seed2_direct = 'N';
		float align2_identity = 0.0;
		
		//seed and extension and calculate identity
		if(read2.size() >= KmerSize+SeedKmerNum)
		{	get_align_seed(kset, read2, 1, read2.size(), contig2_id_index, seed2_contig_start, seed2_contig_end, seed2_read_start, seed2_read_end, seed2_direct);
			if (contig2_id_index != -1)
			{	extend_align_region(read2, contig_seqs[contig2_id_index], align2_identity, seed2_contig_start, seed2_contig_end, seed2_read_start, seed2_read_end, seed2_direct);
				if (align2_identity < MinMapIdentity)
				{	contig2_id_index = -1;
				}
			}
		}
		
		//compare the two alignment of read1 and read2, to see whether they mapped on to different contigs and output
		total_read_pair_num ++;
		if (contig_id_index != -1 && contig2_id_index != -1 )
		{	
			if (contig_id_index != contig2_id_index)
			{	map_ctg_diff_num ++;
				MapCtgDiff << read_id << "\t" << read.size() << "\t" << seed_read_start << "\t" << seed_read_end << "\t" 
				<< contig_ids[contig_id_index] << "\t" << contig_seqs[contig_id_index].size() << "\t" << seed_contig_start << "\t" << seed_contig_end << "\t" << seed_direct << "\t"  << align_identity*100 << "%\t"
				<< read2_id << "\t" << read2.size() << "\t" << seed2_read_start << "\t" << seed2_read_end << "\t"
				<< contig_ids[contig2_id_index] << "\t" << contig_seqs[contig2_id_index].size() << "\t" << seed2_contig_start << "\t" << seed2_contig_end << "\t" << seed2_direct << "\t" 
 << align2_identity*100 << "%" << endl;
			}else
			{	map_ctg_same_num ++;
				MapCtgSame << read_id << "\t" << read.size() << "\t" << seed_read_start << "\t" << seed_read_end << "\t" 
				<< contig_ids[contig_id_index] << "\t" << contig_seqs[contig_id_index].size() << "\t" << seed_contig_start << "\t" << seed_contig_end << "\t" << seed_direct << "\t" << align_identity*100 << "%\t"
				<< read2_id << "\t" << read2.size() << "\t" << seed2_read_start << "\t" << seed2_read_end << "\t" 
				<< contig_ids[contig2_id_index] << "\t" << contig_seqs[contig2_id_index].size() << "\t" << seed2_contig_start << "\t" << seed2_contig_end << "\t" << seed2_direct << "\t" << align2_identity*100 << "%" << endl;
		}
	}else if(contig_id_index != -1 || contig2_id_index != -1)
	{	map_ctg_gap_num ++;
			if (contig_id_index != -1)
			{	MapCtgGap << read_id << "\t" << read.size() << "\t" << seed_read_start << "\t" << seed_read_end << "\t"
				<< contig_ids[contig_id_index] << "\t" << contig_seqs[contig_id_index].size() << "\t" << seed_contig_start << "\t" << seed_contig_end << "\t" << seed_direct << "\t" 
 << align_identity*100 << "%" << endl;
			}
			if (contig2_id_index != -1)
			{	MapCtgGap << read2_id << "\t" << read2.size() << "\t" << seed2_read_start << "\t" << seed2_read_end << "\t"
				<< contig_ids[contig2_id_index] << "\t" << contig_seqs[contig2_id_index].size() << "\t" << seed2_contig_start << "\t" << seed2_contig_end << "\t" << seed2_direct << "\t" 
 << align2_identity*100 << "%" << endl;
			}
		}else
		{	map_no_no_num ++;
		}
	}

	MapCtgStat << "\ttotal_read_pair_num: " << total_read_pair_num << endl;
	MapCtgStat << "\tmap_ctg_diff_num: " << map_ctg_diff_num << "  " <<  (double)map_ctg_diff_num/total_read_pair_num*100 << "%" << endl;
	MapCtgStat << "\tmap_ctg_same_num: " << map_ctg_same_num << "  " << (double)map_ctg_same_num/total_read_pair_num*100 << "%" << endl;
	MapCtgStat << "\tmap_ctg_gap_num: " << map_ctg_gap_num << "  " <<  (double)map_ctg_gap_num/total_read_pair_num*100 << "%" << endl;
	MapCtgStat << "\tmap_no_no_num: " << map_no_no_num << "  " << (double)map_no_no_num/total_read_pair_num*100 << "%" << endl;

}




