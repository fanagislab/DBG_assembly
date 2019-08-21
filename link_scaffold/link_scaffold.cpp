/*
This program links the contigs into scaffolds, can use 2 types of data as input:
a. pair-end reads (the two reads from a pair-end mapped to different contigs)
b. mated-pair reads (the two reads from a mated-pair mapped to different contigs)

Author: Fan Wei, fanweiagis@126.com
Date: 2015-12-6
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
#include "gzstream.h"
#include "link_func.h"
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
using namespace std;


//get the scaffold sequences and related contig and gap information
void read_out_scaffold( vector <string> &contig_seqs, vector<string> &contig_ids, CtgLink **ctgLink,  LinkStat *linkStat,  int contig_num );

//generate one scaffold sequence and its related contig and gap information
void generate_scaffold( vector <string> &contig_seqs, vector<int> &vec_combined, vector<int> &scaff_ids, string &scaff_seq, string &scaff_directs);


int IsMatePair = 0;                     //mapping data type: 0,pair-end; 1, mate-pair

uint64_t total_scaffold_num = 0;       //output scaffold number
uint64_t total_scaffold_len = 0;       //output scaffold length
uint64_t total_scaffold_lenwogap = 0;  //output scaffold length without gap "N"

uint64_t contig_total_num = 0;         //total input contig number
uint64_t contig_total_len = 0;         //total input contig length
uint64_t contig_included_num = 0;      //number of contigs included into the scaffolds
uint64_t contig_included_len = 0;      //length of contigs included into the scaffolds
uint64_t contig_excluded_num = 0;      //excluded repeat and small contig number
uint64_t contig_excluded_len = 0;      //excluded repeat and small contig length

int scaffold_id = -1;


void usage()
{
	cout << "\nFunction instruction:\nlink_scaffold, converts read PEs into contig relations, links the contigs into scaffolds.\
\n(1) If two contigs are neighbor, and they are in the same DNA strand, a pair of PE read will mapped in\
   the F and R way, that is, one read mapped forwardly on one contig, and the other read mapped reversely\
   on the other contig. We use this character to decide the strand relations between two neighbor contigs.\
\n(2) Each contig has two strand forms, we stored both of them in the memory, and each contig node only\
   has one link, the 3'-direction. This can make the strand problem quite easy, and do not increase the\
   memory usage too much.\
\n(3) This program require more than a specified number of read pairs to support two neighbor contigs,\
   and build a graph using neighboring-link-technology, take contig as nodes, and the 3'-link as arcs.\
\n(4) The scaffold sequences were read out from a linear contig, that is, with only one ingoing arc and one\
   outgoing arc, and scaffolding stops at breaking or branching nodes.\n";

	
	cout << "\nlink_scaffold  <contig_file.fa>  <mapping_twoctg_files.lib>\n" 
		<< "   Function: link contigs into scaffolds by pair-end or mate-pair reads, inside gap are not filled"  << endl 
		<< "   Version: 1.0"  << endl
		 << "   -m <int>   input mapping data type: 0, pair-ends; 1. mated-pair,  default=" << IsMatePair << endl		 		 
		 << "   -n <int>   the minimum number of read-pairs required to support a link between two contigs, default=" << PairNumCut << endl	 
	 	 << "   -i <int>   mean insert size for pair-ends or mated-pair reads, default=" << InsertSize << endl
		 << "   -o <str>   the output file prefix, set in commond-line, default = " << Output_prefix << endl
		 << "   -h         get the help information\n" << endl;
	exit(0);
}


void output_parameters()
{
	cerr << "link_scaffold  [version 1.0]"  << endl
		 << "   -m <int>   input mapping data type: 0, pair-ends; 1. mated-pair,  default=" << IsMatePair << endl		 		 
		 << "   -n <int>   the minimum number of read pairs required to support a link between two contigs, default=" << PairNumCut << endl	 
	 	 << "   -i <int>   mean insert size for pair-ends or mated-pair reads, default=" << InsertSize << endl
		 << "   -o <str>   the output prefix, set in commond-line, default = " << Output_prefix << endl
		 << "   -h         get the help information\n" << endl;
			
}


int main(int argc, char *argv[])
{	
	//get options from command line
	int c;
	while((c=getopt(argc, argv, "m:n:i:o:h")) !=-1) {
		switch(c) {
			case 'm': IsMatePair=atoi(optarg); break;
			case 'n': PairNumCut=atoi(optarg); break;
			case 'i': InsertSize=atoi(optarg); break;
			case 'o': Output_prefix=optarg; break;
			case 'h': usage(); break;
			default: usage();
		}
	}
	
	if (argc < 3) usage();

	output_parameters();
	
	string contig_seq_file = argv[optind++];
	string para_map_file = argv[optind++];
	
	clock_t time_start, time_end;
	time_start = clock();

	cerr << "\nProgram start ............" << endl;
	time_end = clock();
	cerr << "Run time: " << double(time_end - time_start) / CLOCKS_PER_SEC << endl;
	
	vector <string> contig_ids;
	vector <string> contig_seqs;
	read_contig_file(contig_seq_file, contig_seqs, contig_ids);
	int contig_num = contig_seqs.size();
	
	for ( int i = 0; i < contig_num; i ++ )
	{	if ( i % 2 == 1)
		{	contig_total_num ++;
			contig_total_len += contig_seqs[i].size();
		}
	}


	cerr << "\nInput contig number: " << contig_total_num << endl;
	cerr << "Input contig length: " << contig_total_len << endl;
	cerr << "Read contigs into memory finished !" << endl;
	
	vector<string> Paired_map_files;
	reading_para_file(para_map_file, Paired_map_files);
	cerr << "\nInput reads mapping files number: " << Paired_map_files.size() << endl;
	
	time_end = clock();
	cerr << "Run time: " << double(time_end - time_start) / CLOCKS_PER_SEC << endl;
	
	CtgLink **ctgLink = new CtgLink* [contig_num]; //二级指针，需要用两个*来定义
	
	//initialization
	for (int i=0; i<contig_num; i++)
	{	ctgLink[i] = NULL;
	}
	
	//load PE map relations into memory
	for (int i=0; i<Paired_map_files.size(); i++)
	{	cerr << "\nparse map file: " << Paired_map_files[i] << endl;
		
		if (IsMatePair == 0)
		{	parse_pair_ends_map_file(Paired_map_files[i], ctgLink, contig_seqs);
		}
		else if(IsMatePair == 1)
		{	parse_mate_pairs_map_file(Paired_map_files[i], ctgLink, contig_seqs);
		}
	}

	cerr << "\nParsed the map files done !" << endl;
	time_end = clock();
	cerr << "Run time: " << double(time_end - time_start) / CLOCKS_PER_SEC << endl;
	
	int Effect_link_num = FR_link_num + RF_link_num + FF_link_num + RR_link_num;
	cerr << "\nFR_link_num: " << FR_link_num << endl;
	cerr << "RF_link_num: " << RF_link_num << endl;
	cerr << "FF_link_num: " << FF_link_num << endl;
	cerr << "RR_link_num: " << RR_link_num << endl;
	cerr << "Effect_link_num: " << Effect_link_num << endl;
	cerr << "Wrong_link_num: " << Wrong_link_num << endl;
	

	//remove low frequency links and stat the link number
	LinkStat *linkStat = new LinkStat[contig_num];
	for ( int i=0; i<contig_num; i++)
	{	linkStat[i].link = 0;  
		linkStat[i].inlink = 0;
		linkStat[i].del = 0;   
	}
	remove_lowfreq_link_and_stat(ctgLink, contig_num, linkStat);
	
	cerr << "\nRemoved LowFreq link num: " << LowFreq_link_num << endl;

	//statistics of links
	uint64_t total_link_num = 0;
	uint64_t uniq_link_num = 0;
	uint64_t multiple_link_num = 0;  //repeat or interleaving caused
	uint64_t empty_link_num = 0;
	for ( int i = 0; i < contig_num; i ++ )
	{	if ( i % 2 == 1)
		{	int links_num = linkStat[i].link;
			if (links_num == 0)
			{	empty_link_num ++;
			}else if (links_num == 1)
			{	uniq_link_num ++;
			}else if (links_num > 1)
			{	multiple_link_num ++;
			}
			total_link_num ++;
		}
	}
	cerr << "Number and ratio of contigs having a unique 3'-link: " << uniq_link_num << "  " << (float)uniq_link_num / total_link_num << endl;
	cerr << "Number and ratio of contigs having multiple 3'-link: " << multiple_link_num << "  " << (float)multiple_link_num / total_link_num << endl;
	cerr << "Number and ratio of contigs having zero 3'-link:     " << empty_link_num << "  " << (float)empty_link_num / total_link_num << endl;

	
	//output the links data
	string link_out_file = Output_prefix + ".insert" + boost::lexical_cast<string>(InsertSize) + ".scaffold.links.all";
	display_data_in_link(ctgLink, linkStat, contig_num, link_out_file);

	//solve the interleaving problem for pair-end and mated-pair mapping data
	//must remove interleaving links before remove repeat nodes
	remove_interleaving_links(ctgLink, contig_num, linkStat, contig_seqs);
	cerr << "\nRemoved interleave links num: " << Interleave_link_num << endl;

	//remove repeat nodes
	vector<int> repeat_nodes_vec;
	remove_repeat_nodes(contig_num, linkStat, repeat_nodes_vec, contig_seqs);

	cerr << "\nRemoved repeat nodes num: " << repeat_nodes_vec.size() / 2 << endl;

	
	//to resolve interleaving problem, remove small nodes whose length are less than half of the Insert_size
	//vector<int> small_nodes_vec;
	//remove_small_nodes(contig_num, linkStat, small_nodes_vec, contig_seqs);
	//cerr << "\nRemoved small nodes num: " << small_nodes_vec.size() / 2 << endl;


	//remove all the links belong to deleted nodes
	remove_links_from_deleted_nodes(ctgLink, contig_num, linkStat);
	cerr << "\nRemoved links [related with repeat or small nodes] num: " << Deleted_link_num << endl;

	
	//output the links data
	link_out_file = Output_prefix + ".insert" + boost::lexical_cast<string>(InsertSize) + ".scaffold.links.uniq";
	display_data_in_link(ctgLink, linkStat, contig_num, link_out_file);
	
	read_out_scaffold( contig_seqs, contig_ids, ctgLink,  linkStat, contig_num );
	
	//output singlet sequences not been incorpated into scaffold sequences
	string singlet_seq_file = Output_prefix + ".insert" + boost::lexical_cast<string>(InsertSize) + ".scaffold_repeat.seq.fa";
	ofstream SingletFile ( singlet_seq_file.c_str() );
	if ( ! SingletFile )
	{   cerr << "fail to open file" << singlet_seq_file << endl;
	}
	string singlet_pos_file = Output_prefix + ".insert" + boost::lexical_cast<string>(InsertSize) + ".scaffold_repeat.pos.tab";
	ofstream SingletPosFile ( singlet_pos_file.c_str() );
	if ( ! SingletPosFile )
	{   cerr << "fail to open file" << singlet_pos_file << endl;
	}

	vector<LenAndSeq> LenSeqVec;
	for (int i = 0; i < repeat_nodes_vec.size(); i ++)
	{	int ctg_id = repeat_nodes_vec[i];
		if ( ctg_id % 2 == 1 ) 
		{	
			LenAndSeq this_node;
			this_node.len = contig_seqs[ctg_id].size();
			this_node.seq = contig_seqs[ctg_id];
			this_node.pos = "\t" + contig_ids[ctg_id] +  "\t" + boost::lexical_cast<string>(1)  +  "\t" + boost::lexical_cast<string>(contig_seqs[ctg_id].size())
 +  "\t" + boost::lexical_cast<string>(contig_seqs[ctg_id].size())  +  "\t" + "F" + "\n";
			LenSeqVec.push_back(this_node);

			contig_excluded_num ++;
			contig_excluded_len += contig_seqs[ctg_id].size();
		}
	}
	sort(LenSeqVec.begin(), LenSeqVec.end(), cmpSeqByLen);

	for(int i = 0; i < LenSeqVec.size(); i ++)
	{	scaffold_id += 2;
		SingletFile << ">scf_" + boost::lexical_cast<string>(scaffold_id) + "   fragment_num:" + boost::lexical_cast<string>(1) + "   length:" + boost::lexical_cast<string>(LenSeqVec[i].len) + "   lenwogap:" + boost::lexical_cast<string>(LenSeqVec[i].len) + "   RepeatNode\n" << LenSeqVec[i].seq << "\n";
		SingletPosFile << ">scf_" + boost::lexical_cast<string>(scaffold_id) + "\n" << LenSeqVec[i].pos;
	}

	

		
	cerr << "\nRead out scaffold sequence done" << endl;
	cerr << "\nTotal scaffold number:          " << total_scaffold_num << endl;
	cerr << "Total scaffold length[WithGap]: " << total_scaffold_len << endl;
	cerr << "Total scaffold length[NoGap]:   " << total_scaffold_lenwogap << endl;
	
	cerr << "\nIncluded contig number: " << contig_included_num << "  "<< (float)contig_included_num / contig_total_num  << endl;
	cerr << "Included contig length: " << contig_included_len << "  "<< (float)contig_included_len / contig_total_len  << endl;
	cerr << "Excluded repeat contig number: " << contig_excluded_num << "  "<< (float)contig_excluded_num / contig_total_num  << endl;
	cerr << "Excluded repeat contig length: " << contig_excluded_len << "  "<< (float)contig_excluded_len / contig_total_len  << endl;


	cerr << "\nProgram finished !" << endl;
	time_end = clock();
	cerr << "Run time: " << double(time_end - time_start) / CLOCKS_PER_SEC << endl;
	
}



//get the scaffold sequences and related contig and gap information
void read_out_scaffold( vector <string> &contig_seqs, vector<string> &contig_ids, CtgLink **ctgLink,  LinkStat *linkStat,  int contig_num )
{	
	string scaff_seq_file = Output_prefix + ".insert" + boost::lexical_cast<string>(InsertSize) + ".scaffold.seq.fa";
	string scaff_pos_file = Output_prefix + ".insert" + boost::lexical_cast<string>(InsertSize) + ".scaffold.pos.tab";

	ofstream ScafPosFile ( scaff_pos_file.c_str() );
	if ( ! ScafPosFile )
	{   cerr << "fail to open file" << scaff_pos_file << endl;
	}

	ofstream ScafSeqFile (scaff_seq_file.c_str());
	if (! ScafSeqFile)
	{   cerr << "fail to open file" << scaff_seq_file << endl;
	}
	
	vector<LenAndSeq> LenSeqVec;
	
	for ( int i=1; i<contig_num; i++)
	{	
		if ( linkStat[i].del == 1 || i % 2 == 0  )
		{	continue;
		}
			
		vector<int> vec_left_id_and_gap, vec_right_id_and_gap;
		vector<int> combined_ids_and_gaps;
		
		//先删除当前节点，以防出现环路时，两端为同一个contig
		linkStat[i].del = 1;
		linkStat[get_pair_id(i)].del = 1;

		//get rightward set of non-branching contigs
		if ( linkStat[i].link == 1 )
		{	get_linear_seq( ctgLink, linkStat, i, vec_right_id_and_gap);
		}
		
		//get leftward set of non-branching contigs
		int paired_id = get_pair_id(i);
		if (linkStat[paired_id].link == 1 )
		{	
			get_linear_seq( ctgLink, linkStat, paired_id, vec_left_id_and_gap);
			
			//get reverse and complement strand
			reverse( vec_left_id_and_gap.begin(), vec_left_id_and_gap.end() );
			for (int k=0; k<vec_left_id_and_gap.size(); k++)
			{	if(k % 2 == 0)
				{	vec_left_id_and_gap[k] = get_pair_id(vec_left_id_and_gap[k]);
				}
			}
		}
		
		//combine the leftward, this, and rightward contig ids and inside gaps into one vector
		for (int z =0; z < vec_left_id_and_gap.size(); z++)
		{	combined_ids_and_gaps.push_back( vec_left_id_and_gap[z] );
		}
		combined_ids_and_gaps.push_back(i);
		for (int z =0; z < vec_right_id_and_gap.size(); z++)
		{	combined_ids_and_gaps.push_back( vec_right_id_and_gap[z] );
		}


		//generate and output the scaffold sequence
		string scaff_seq;  //store the final result scaffold sequence, fill gaps with Ns.
		vector<int> scaff_ids; //use minus value to store gap size, use plus value to store contig_id
		string scaff_directs;  //store the directions of each contigs
		generate_scaffold(contig_seqs, combined_ids_and_gaps, scaff_ids, scaff_seq, scaff_directs);
		
		int scaf_ctg_num = 0;
		int scaf_len = 0;
		int scaf_lenwogap = 0;
		int block_start = 0;
		int block_size = 0;
		
		string output_scafpos_res;
		string output_scafseq_res;

		for (int i=0; i<scaff_ids.size(); i++)
		{	
			if ( i % 2 == 0 )   //store contig ids
			{	scaf_ctg_num ++;
				block_start = scaf_len + 1;
				block_size = contig_seqs[scaff_ids[i]].size();
				scaf_len += block_size;
				scaf_lenwogap += block_size;
				
			    output_scafpos_res += "\t" + contig_ids[scaff_ids[i]] + "\t" + boost::lexical_cast<string>(block_start) + "\t" + boost::lexical_cast<string>(scaf_len) + "\t" + boost::lexical_cast<string>(block_size) + "\t" + scaff_directs[i] + "\n";
				contig_included_num ++;
				contig_included_len += block_size;
			}else               //store gap sizes
			{	block_size = scaff_ids[i];
				block_start = scaf_len + 1;
				scaf_len += block_size;
				output_scafpos_res += "\tgap\t"  + boost::lexical_cast<string>(block_start) + "\t" + boost::lexical_cast<string>(scaf_len) + "\t" + boost::lexical_cast<string>(block_size)  + "\t" + scaff_directs[i] + "\n";
			
			}
		}
		
		output_scafseq_res += "   fragment_num:" + boost::lexical_cast<string>(scaf_ctg_num) + "   length:" + boost::lexical_cast<string>(scaf_len) + "   lenwogap:" + boost::lexical_cast<string>(scaf_lenwogap) + "\n" + scaff_seq + "\n";
		
		LenAndSeq this_node;
		this_node.len = scaf_len;
		this_node.seq = output_scafseq_res;
		this_node.pos = output_scafpos_res;
		LenSeqVec.push_back(this_node);

		total_scaffold_num ++;
		total_scaffold_len += scaf_len;
		total_scaffold_lenwogap += scaf_lenwogap;
	}
	
	//////////////////////////
	//sort the scafftig results by length
	sort(LenSeqVec.begin(), LenSeqVec.end(), cmpSeqByLen);

	for(int i = 0; i < LenSeqVec.size(); i ++)
	{	scaffold_id += 2;
		ScafSeqFile << ">scf_" + boost::lexical_cast<string>(scaffold_id) << LenSeqVec[i].seq;
		ScafPosFile << ">scf_" + boost::lexical_cast<string>(scaffold_id) + "\n" << LenSeqVec[i].pos;
	}
	
	ScafPosFile.close();
	ScafSeqFile.close();


}


//generate one scaffold sequence and its related contig and gap information
void generate_scaffold( vector <string> &contig_seqs, vector<int> &vec_combined, vector<int> &scaff_ids, string &scaff_seq, string &scaff_directs)
{	
	//parse the leftward contig ids
	for (int i=0; i<vec_combined.size(); i++)
	{	
		if ( i % 2 == 0)   // vec_combined[i] store contig ids
		{	
			int ctg_id = vec_combined[i];
			if ( ctg_id % 2 == 1 )
			{	scaff_ids.push_back(ctg_id);
				scaff_seq += contig_seqs[ctg_id];
				scaff_directs.push_back('F');
			} 
			else
			{	int paired_id = get_pair_id(ctg_id);
				scaff_ids.push_back( paired_id );
				string contig_rc;
				reverse_complement(contig_seqs[paired_id] , contig_rc);
				scaff_seq += contig_rc;
				scaff_directs.push_back('R');

			}
		
		}else              // vec_combined[i] store gap sizes
		{	
			int gap_size = (vec_combined[i] > 1) ? vec_combined[i] : 1;  //For visual purpose, the minimum output gap size is 1
			scaff_ids.push_back(gap_size);
			string Nstr;
			generate_Nstr(gap_size , Nstr );
			
			scaff_seq += Nstr;
			scaff_directs.push_back('N');
		}

	}
	
}


