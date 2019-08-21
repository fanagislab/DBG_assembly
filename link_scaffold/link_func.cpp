/*
This function convert read-end and read-pair relations into contig relations, links the contigs into scaffolds.
Deal with 3 types of input data: a. single reads data; b. illumina pair-end data; c. illumina mate-pair data.


(1) If two contigs are neighbor, and they are in the same DNA strand, a pair of PE read will mapped in 
   the F and R way, that is, one read mapped forwardly on one contig, and the other read mapped reversely
   on the other contig. However, a pair of mated-PE will be mapped in the R and F way. Moreover, using a single
   long reads, the two ends from a single reads mapped to two separate contigs will be F and F (or R and R).
   We use this character to decide the strand relations between two neighbor contigs.

(2) Each contig has two strand forms, we stored both of them in the memory, and each contig node only
   has one link, the 3'-direction. This can make the strand problem quite easy, and do not increase the
   memory usage too much.

(3) This program require more than a specified number of read pairs to support two neighbor contigs,
   and build a graph using neighboring-link-technology, take contig as nodes, and the 3'-link as arcs.

(4) The scaffold sequences were read out from a linear contig, that is, with only one ingoing arc and one
   outgoing arc, and scaffolding stops at breaking or branching nodes.

(5) The key point for scaffolding, is to deal with the repeat nodes, and interleaving links caused by
   small fragments, which have length less than the library insert sizes. We removed the repeat nodes
   with >=2 incoming and >=2 outgoing links, and small nodes with length <= InsertSize/2, as well as
   the related links.

(6) Another thing to remind, is the required pair number for each insert sizes, this paramter affects
   the assembled scaffold length seriously. For small insert 170bp, 500bp, 800 bp, it is suitable to
   set 3 pairs. For large insert 2kb, 5kb, and 10bp, this parameter should set larger, in order to 
   distinguish real links and random links. In a testing, we get good result, with the following
   parameter file.


Author: Fan Wei, fanw@genomics.org.cn
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

using namespace std;

string Output_prefix = "Output";          //output file prefix
int InsertSize = 400;                     //insert size of pair-end and mate-pair data
int PairNumCut = 3;                       //cutoff for number of read-pairs to support a contig link

int FR_link_num = 0;                // read1 forward mapping, read2 reverse mapping
int RF_link_num = 0;                // read1 reverse mappin,g read2 forward mapping
int FF_link_num = 0;                // both read1 and read2 forward mapping
int RR_link_num = 0;                // both read1 and read2 reverse mapping
int Wrong_link_num = 0;             // wrong links which do not belong to the above four mapping types
int LowFreq_link_num = 0;           // deleted links number which have supporting read pair number less than PairNumCut
int Deleted_link_num = 0;           // deleted links number derived from deleted repeat and small contig nodes
int Interleave_link_num = 0;        // deleted interleaving links number



//used for sorting the gap-seqs by its length
bool cmpSeqByLen (LenAndSeq a, LenAndSeq b)
{       return b.len < a.len;
}


//read file_list into a vector
void reading_para_file(string &para_file, vector<string> &Paired_map_files)
{	
	ifstream infile ( para_file.c_str() );
	if ( ! infile )
	{	cerr << "fail to open input file" << para_file << endl;
	}

	string line;
	while ( getline( infile, line, '\n' ) )
	{	
		if (line[0] == '#')
		{	continue;
		}
		vector<string> vec_line;
		split(line, vec_line, " \t\n");
		if (vec_line.size() == 0)
		{	continue;
		}
		Paired_map_files.push_back(vec_line[0]);
	}
}


//read multiple-fasta format contig file
void read_contig_file(string &contig_seq_file, vector <string> &contig_seqs, vector <string> &contig_ids)
{
	ifstream infile ( contig_seq_file.c_str() );
	if ( ! infile )
	{       cerr << "fail to open input file " << contig_seq_file << endl;
	}
	
	contig_seqs.push_back("");
	contig_ids.push_back("");

	string contig_str, id_str, line;
	while ( getline( infile, line, '\n' ) )
	{       
		if (line[0] == '>') 
		{	
			vector<string> vec_head;
			split(line, vec_head, "> \t");
			id_str = vec_head[0];
			contig_ids.push_back(id_str);
			contig_ids.push_back("");

			if(contig_str.size() > 0){
				contig_seqs.push_back(contig_str);
				contig_seqs.push_back("");	
			}
			contig_str.clear();
		}
		else
		{	contig_str += line;     
		}
	}

	if(contig_str.size() >0){
		contig_seqs.push_back(contig_str);
		contig_seqs.push_back("");
	}

}

//parse the read-ends mapping file, add the links information into the memory
//gap size is estimated by read aligning coordinates
//read_400_3/1    250     1       138     ctg_763 706     569     706     F       100%    read_400_3/1    250     109     250     ctg_563 7685    1       142     F       10
void parse_read_ends_map_file(string &PE_map_file, CtgLink **ctgLink, vector <string> &contig_seqs)
{
	igzstream infile ( PE_map_file.c_str() );
	if ( ! infile )
	{       cerr << "fail to open input file " << PE_map_file << endl;
	}
	
	string line;
	while ( getline( infile, line, '\n' ) )
	{	
		if ( line[0] == '#' )
		{	continue;
		}

		string line_store = line;
		vector<string> vec_line;
		split(line, vec_line, " \t\n");
		
		int ctg1_id = 0;
		int ctg2_id = 0;
		int ctg3_id = 0;
		int ctg4_id = 0;
		int gap_size = 0;  //estimated gap size

		
		//根据字段意思，提取构建scaffold需要的信息
		string direct1 = vec_line[8];
		string direct2 = vec_line[18];
		string contig1 = vec_line[4];
		string contig2 = vec_line[14];

		string read_align1_end = vec_line[3];
		string read_align2_start = vec_line[12];

		if (direct1 == "F" && direct2 == "F")
		{	
			ctg1_id = ctgStr2Id(contig1);
			ctg2_id = ctg1_id + 1;
			ctg3_id = ctgStr2Id(contig2);
			ctg4_id = ctg3_id + 1;
			gap_size = atoi(read_align2_start.c_str()) - atoi(read_align1_end.c_str()) - 1;
			FF_link_num ++;
		}
		else if (direct1 == "R" && direct2 == "R")
		{	
			ctg2_id = ctgStr2Id(contig1);
			ctg1_id = ctg2_id + 1;
			ctg4_id = ctgStr2Id(contig2);
			ctg3_id = ctg4_id + 1;
			gap_size = atoi(read_align2_start.c_str()) - atoi(read_align1_end.c_str()) - 1;
			RR_link_num ++;
		}
		else if (direct1 == "F" && direct2 == "R")
		{	
			ctg1_id = ctgStr2Id(contig1);
			ctg2_id = ctg1_id + 1;
			ctg4_id = ctgStr2Id(contig2);
			ctg3_id = ctg4_id + 1;
			gap_size = atoi(read_align2_start.c_str()) - atoi(read_align1_end.c_str()) - 1;
			FR_link_num ++;
		}
		else if (direct1 == "R" && direct2 == "F")
		{
			ctg2_id = ctgStr2Id(contig1);
			ctg1_id = ctg2_id + 1;
			ctg3_id = ctgStr2Id(contig2);
			ctg4_id = ctg3_id + 1;
			gap_size = atoi(read_align2_start.c_str()) - atoi(read_align1_end.c_str()) - 1;
			RF_link_num ++;
		}
		else
		{	Wrong_link_num ++;
			continue;
		}
		
		add_data_into_link(ctgLink, ctg1_id, ctg3_id, gap_size);
		add_data_into_link(ctgLink, ctg4_id, ctg2_id, gap_size);
	}

}


//parse the pair-end mapping file, add the links information into the memory
//gap size is estimated by contig aligning coordinates//
//read_400_3/1    250     1       138     ctg_763 706     569     706     F       100%    read_400_3/1    250     109     250     ctg_563 7685    1       142     F       10
void parse_pair_ends_map_file(string &PE_map_file, CtgLink **ctgLink, vector <string> &contig_seqs)
{
	igzstream infile ( PE_map_file.c_str() );
	if ( ! infile )
	{       cerr << "fail to open input file " << PE_map_file << endl;
	}
	
	string line;
	while ( getline( infile, line, '\n' ) )
	{	
		if ( line[0] == '#' )
		{	continue;
		}

		string line_store = line;
		vector<string> vec_line;
		split(line, vec_line, " \t\n");
		
		int ctg1_id = 0;
		int ctg2_id = 0;
		int ctg3_id = 0;
		int ctg4_id = 0;
		int gap_size = 0;  //estimated gap size
		int ctg1_start = 0;
		int ctg3_end = 0;
		
		//根据字段意思，提取构建scaffold需要的信息
		string direct1 = vec_line[8];
		string direct2 = vec_line[18];
		string contig1 = vec_line[4];
		string contig2 = vec_line[14];
		string contig1_start = vec_line[6];
		string contig1_end = vec_line[7];
		string contig2_start = vec_line[16];
		string contig2_end = vec_line[17];

		if (direct1 == "F" && direct2 == "R")
		{	
			ctg1_id = ctgStr2Id(contig1);
			ctg2_id = ctg1_id + 1;
			ctg3_id = ctgStr2Id(contig2);
			ctg4_id = ctg3_id + 1;
			ctg1_start = atoi(contig1_start.c_str());
			ctg3_end = atoi(contig2_end.c_str());
			gap_size = InsertSize - (contig_seqs[ctg1_id].size() -ctg1_start) - ctg3_end; 
			//cerr << "##gap_size: " << gap_size << endl;
			FR_link_num ++;
		}
		else if (direct1 == "R" && direct2 == "F")
		{	
			ctg1_id = ctgStr2Id(contig2);
			ctg2_id = ctg1_id + 1;
			ctg3_id = ctgStr2Id(contig1);
			ctg4_id = ctg3_id + 1;
			ctg1_start = atoi(contig2_start.c_str());
			ctg3_end = atoi(contig1_end.c_str());
			gap_size = InsertSize - (contig_seqs[ctg1_id].size() - ctg1_start) - ctg3_end;
			//cerr << "##gap_size: " << gap_size << endl;
			RF_link_num ++;
		}
		else if (direct1 == "F" && direct2 == "F")
		{	
			ctg1_id = ctgStr2Id(contig1);
			ctg2_id = ctg1_id + 1;
			ctg4_id = ctgStr2Id(contig2);
			ctg3_id = ctg4_id + 1;
			ctg1_start = atoi(contig1_start.c_str());
			ctg3_end = contig_seqs[ctg4_id].size() - atoi(contig2_start.c_str());
			gap_size = InsertSize - (contig_seqs[ctg1_id].size() - ctg1_start) - ctg3_end;
			//cerr << "##gap_size: " << gap_size << endl;
			FF_link_num ++;
		}
		else if (direct1 == "R" && direct2 == "R")
		{
			ctg2_id = ctgStr2Id(contig1);
			ctg1_id = ctg2_id + 1;
			ctg3_id = ctgStr2Id(contig2);
			ctg4_id = ctg3_id + 1;
			ctg1_start = contig_seqs[ctg2_id].size() -  atoi(contig1_end.c_str());
			ctg3_end = atoi(contig2_end.c_str());
			gap_size = InsertSize - (contig_seqs[ctg2_id].size() - ctg1_start) - ctg3_end;
			//cerr << "##gap_size: " << gap_size << endl;
			RR_link_num ++;
		}
		else
		{	Wrong_link_num ++;
			continue;
		}
		
		//the maximum calculated gap_size value is equal to InsertSize
		//when the neighboring contigs have overlap, the calculated gap_size is minus value
		if(gap_size > -InsertSize/2 && gap_size <= InsertSize)
		{
			add_data_into_link(ctgLink, ctg1_id, ctg3_id, gap_size);
			add_data_into_link(ctgLink, ctg4_id, ctg2_id, gap_size);
		}
	}

}

//parse the mate-pair mapping file, add the links information into the memory
//gap size is estimated by contig aligning coordinates
//read_400_3/1    250     1       138     ctg_763 706     569     706     F       100%    read_400_3/1    250     109     250     ctg_563 7685    1       142     F       10
void parse_mate_pairs_map_file(string &PE_map_file, CtgLink **ctgLink, vector <string> &contig_seqs)
{
	igzstream infile ( PE_map_file.c_str() );
	if ( ! infile )
	{       cerr << "fail to open input file " << PE_map_file << endl;
	}
	
	string line;
	while ( getline( infile, line, '\n' ) )
	{	
		if ( line[0] == '#' )
		{	continue;
		}

		string line_store = line;
		vector<string> vec_line;
		split(line, vec_line, " \t\n");
		
		int ctg1_id = 0;
		int ctg2_id = 0;
		int ctg3_id = 0;
		int ctg4_id = 0;
		int gap_size = 0;  //estimated gap size
		int ctg1_start = 0;
		int ctg3_end = 0;
		
		//根据字段意思，提取构建scaffold需要的信息
		string direct1 = vec_line[8];
		string direct2 = vec_line[18];
		string contig1 = vec_line[4];
		string contig2 = vec_line[14];
		string contig1_start = vec_line[6];
		string contig1_end = vec_line[7];
		string contig2_start = vec_line[16];
		string contig2_end = vec_line[17];


		//mate pair各种情况与non-mate pair正好完全相反
		if (direct1 == "F" && direct2 == "R")
		{	
			ctg1_id = ctgStr2Id(contig2);
			ctg2_id = ctg1_id + 1;
			ctg3_id = ctgStr2Id(contig1);
			ctg4_id = ctg3_id + 1;
			ctg1_start = atoi(contig2_start.c_str());
			ctg3_end = atoi(contig1_end.c_str());
			gap_size = InsertSize - (contig_seqs[ctg1_id].size() - ctg1_start) - ctg3_end;
			FR_link_num ++;
		}
		else if (direct1 == "R" && direct2 == "F")
		{	
			ctg1_id = ctgStr2Id(contig1);
			ctg2_id = ctg1_id + 1;
			ctg3_id = ctgStr2Id(contig2);
			ctg4_id = ctg3_id + 1;
			ctg1_start = atoi(contig1_start.c_str());
			ctg3_end = atoi(contig2_end.c_str());
			gap_size = InsertSize - (contig_seqs[ctg1_id].size() - ctg1_start) - ctg3_end; 
			RF_link_num ++;
		}
		else if (direct1 == "F" && direct2 == "F")
		{	
			ctg2_id = ctgStr2Id(contig1);
			ctg1_id = ctg2_id + 1;
			ctg3_id = ctgStr2Id(contig2);
			ctg4_id = ctg3_id + 1;
			ctg1_start = contig_seqs[ctg2_id].size() - atoi(contig1_end.c_str());
			ctg3_end = atoi(contig2_end.c_str());
			gap_size = InsertSize -  (contig_seqs[ctg2_id].size() - ctg1_start) - ctg3_end;

			FF_link_num ++;
		}
		else if (direct1 == "R" && direct2 == "R")
		{	
			ctg1_id = ctgStr2Id(contig1);
			ctg2_id = ctg1_id + 1;
			ctg4_id = ctgStr2Id(contig2);
			ctg3_id = ctg4_id + 1;
			ctg1_start = atoi(contig1_start.c_str());
			ctg3_end = contig_seqs[ctg4_id].size() - atoi(contig2_start.c_str());
			gap_size = InsertSize - (contig_seqs[ctg1_id].size() - ctg1_start) - ctg3_end;
			RR_link_num ++;
		}
		else
		{	Wrong_link_num ++;
			continue;
		}
		
		//the maximum calculated gap_size value is equal to InsertSize
		//when the neighboring contigs have overlap, the calculated gap_size is minus value
		if(gap_size > -InsertSize/2 && gap_size <= InsertSize)
		{
			add_data_into_link(ctgLink, ctg1_id, ctg3_id, gap_size);
			add_data_into_link(ctgLink, ctg4_id, ctg2_id, gap_size);
		}
	}

}


//add one link information into the link data structure
void add_data_into_link(CtgLink **ctgLink, int source_id, int target_id, int gap_size)
{	
	int target_exist = 0;
	
	//add the first element
	if (ctgLink[source_id] == NULL)
	{	CtgLink *e = new CtgLink;
		e->id = target_id;
		e->freq = 1;
		e->size = gap_size;
		e->pointer = NULL;
		ctgLink[source_id] = e;
		return;
	}
	
	CtgLink *p = ctgLink[source_id];

	while (1)
	{	if (p->id == target_id)
		{	target_exist = 1;
			break;
		}
		if (p->pointer == NULL)
		{	break;
		}
		p = p->pointer;
	}

	if (target_exist == 1)
	{	if (p->freq < 1023)
		{	
			p->freq ++;
			p->size += gap_size;
		}
		
	}else
	{	CtgLink *e = new CtgLink;
		e->id = target_id;
		e->freq = 1;
		e->size = gap_size;
		e->pointer = NULL;
		p->pointer = e;
	}
}


//remove the links with low frequency < PairNumCut
void remove_lowfreq_link_and_stat(CtgLink **ctgLink, int contig_num, LinkStat *linkStat)
{	
	for ( int i=0; i<contig_num; i++)
	{	
		if (ctgLink[i] == NULL)
		{	continue;	
		}

		int link_num = 0;
		CtgLink *p = ctgLink[i];

		while (1)
		{	
			if (p->freq < PairNumCut)
			{	
				p->id =0;
				p->freq = 0;
				p->size = 0;
				LowFreq_link_num ++;
			}else{
				link_num ++;
				if (linkStat[p->id].inlink < 255)
				{	linkStat[p->id].inlink ++ ;
				}
			}

			if (p->pointer == NULL)
			{	break;
			}
			p = p->pointer;
		}
		
		linkStat[i].link = (link_num < 255) ? link_num : 255;
	}
}


//output data in the link data-structure
void display_data_in_link(CtgLink **ctgLink, LinkStat *linkStat, int contig_num, string &link_out_file)
{	
	ofstream outFile ( link_out_file.c_str() );
	if ( ! outFile )
	{	cerr << "fail to open file" << link_out_file << endl;
	}
	
	outFile << "ctg_id\tincoming_link_num\toutgoing_link_num\tlinked_id,pair_num,sum_size,avg_size;\n";
	for (int i=1; i<contig_num; i++)
	{
		outFile << i << "\t" << int(linkStat[i].inlink) << "\t" << int(linkStat[i].link);

		CtgLink *p = ctgLink[i];
		while (p != NULL)
		{	if (p->freq > 0)
			{	outFile << "\t" << p->id << "," << p->freq << "," << p->size << "," << p->size / p->freq;
			}

			p = p->pointer;
		}
		outFile << endl;
	}
}


//mainly used for solving pair-end or mate-pair data, because single-read mapping can avoid interleaving problem
//solve the interleaving problem, by removing the interleaving links
//the middle contig in a interleaving case, has only one inlink and only one outlink
void remove_interleaving_links(CtgLink **ctgLink, int contig_num, LinkStat *linkStat, vector<string> &contig_seqs)
{	
	for (int i=1; i<contig_num; i++)
	{	
		int start_node = i;
		if (linkStat[start_node].del == 0 && linkStat[start_node].link == 2)
		{	vector<int> linked_ids; 
			vector<int> gap_sizes;
			get_all_linked_ids(ctgLink, start_node, linked_ids, gap_sizes);
			
			if (linkStat[linked_ids[0]].link == 1 && linkStat[linked_ids[0]].inlink == 1)
			{	int middle_ctg_id = (linked_ids[0] % 2 == 1) ? linked_ids[0] : (linked_ids[0] - 1);
				int judge_len = gap_sizes[1] * 2;
				int end_insert = 0;
				int end_node = get_next_linked_id( ctgLink, linked_ids[0], end_insert );

				if (end_node == linked_ids[1] && gap_sizes[0] < judge_len && end_insert < judge_len && contig_seqs[middle_ctg_id].size() < judge_len )
				{	//cerr << "##**## " << start_node << "\t" << end_node << endl;
					delete_linked_id( ctgLink, linkStat, start_node,  end_node);
					Interleave_link_num ++;
				}
			}

			if (linkStat[linked_ids[1]].link == 1 && linkStat[linked_ids[1]].inlink == 1)
			{	int middle_ctg_id = (linked_ids[1] % 2 == 1) ? linked_ids[1] : (linked_ids[1] - 1);
				int judge_len = gap_sizes[0] * 2;
				int end_insert = 0;
				int end_node = get_next_linked_id( ctgLink, linked_ids[1], end_insert );
				
				if (end_node == linked_ids[0]  && gap_sizes[1] < judge_len && end_insert < judge_len && contig_seqs[middle_ctg_id].size() < judge_len )
				{	//cerr << "#**#**## " << start_node << "\t" << end_node << endl;
					delete_linked_id( ctgLink, linkStat, start_node,  end_node);
					Interleave_link_num ++;
				}
			}

		}
	}
}


//this routine is specially designed for pacbio scaffoling, 
//because of unspecific mapping, pacbio reads mapping to contigs will generate many missing hits that cause interleaving problem
//search all the possible path to find interleaving potential
void strong_remove_interleaving_links(CtgLink **ctgLink, int contig_num, LinkStat *linkStat, vector<string> &contig_seqs)
{	
	int Rank_Num = 2;
	
	for (int i=1; i<contig_num; i++)
	{	
		
		int start_node = i;


		if (linkStat[start_node].del == 0 && linkStat[start_node].link >= 2 && linkStat[start_node].link <= 3)
		{	
			
			vector<int> linked_ids; 
			vector<int> gap_sizes;
			get_all_linked_ids(ctgLink, start_node, linked_ids, gap_sizes);
			

			vector<vector<int> > children_ids;

			for (int j = 0; j < linked_ids.size(); j ++)
			{	
				vector<int> children_vec;

				//get all the children ids for a starting father node
				int ori_id = linked_ids[j];
				vector<int> cur_vec;
				cur_vec.push_back(ori_id);
				int rank_num = 0;
				while (rank_num < Rank_Num)
				{	
					vector<int> tmp_vec; //store all the ids for this rank

					for (int k = 0; k < cur_vec.size(); k ++ )
					{
						vector<int> linkedids;
						vector<int> gapsizes;
						get_all_linked_ids(ctgLink, cur_vec[k], linkedids, gapsizes);
						for (int m = 0; m < linkedids.size(); m ++)
						{	tmp_vec.push_back(linkedids[m]);
							children_vec.push_back(linkedids[m]);
						}
					}
					cur_vec = tmp_vec;
					rank_num ++;
				}
				/////////////////////////////////////
				
				children_ids.push_back(children_vec);
			
				
			}

			//decide which child-node should be deleted
			map<int,int> deleted_ids;
			for (int j = 0; j < linked_ids.size(); j ++)
			{	
				vector<int> &childrenids = children_ids[j];
				for (int k = 0; k < linked_ids.size(); k ++)
				{	
					vector<int>::iterator s = find( childrenids.begin(), childrenids.end(), linked_ids[k] );
					 if( s != childrenids.end() ) //找到
					 {	deleted_ids[linked_ids[k]] = 1;
					 }
				}
			}
			
			//delete the related links and reduce related stats
			map<int,int>::iterator myiter;
			for (myiter = deleted_ids.begin(); myiter != deleted_ids.end(); myiter ++ )
			{	int source_id = start_node;
				int target_id = myiter->first;
				delete_linked_id( ctgLink, linkStat, source_id,  target_id);
				Interleave_link_num ++;

			}
			
		}
	}
}



//delete a specified link between source node and target node
void delete_linked_id( CtgLink **ctgLink, LinkStat *linkStat, int source_id, int target_id )
{
	CtgLink *p = ctgLink[source_id];
	while (p != NULL) 
	{	if (p->freq > 0) 
		{	
			if (p->id == target_id)
			{	p->id = 0;
				p->freq = 0;
				p->size = 0;

				if (linkStat[source_id].link > 0)
				{	linkStat[source_id].link --;
				}
				if (linkStat[target_id].inlink > 0)
				{	linkStat[target_id].inlink --;
				}
				break;
			}
		}
		p = p->pointer;
	}
	
}


//get all the linked contig ids for a given contig node
void get_all_linked_ids( CtgLink **ctgLink, int source_id, vector<int> &linked_ids, vector<int> &gap_sizes )
{
	CtgLink *p = ctgLink[source_id];
	while (p != NULL) 
	{	if (p->freq > 0) 
		{	linked_ids.push_back(p->id);
			int gap_size = p->size / p->freq;
			gap_sizes.push_back(gap_size);
		}
		p = p->pointer;
	}
	
}

//remove repeat nodes which have multiple incoming and outgoing links
void remove_repeat_nodes(int contig_num, LinkStat *linkStat, vector<int> &repeat_nodes, vector<string> &contig_seqs)
{
	for (int i=1; i<contig_num; i++)
	{	
		if (linkStat[i].del == 0 &&  (linkStat[i].inlink >= 2 || linkStat[i].link >= 2) )
		{	repeat_nodes.push_back(i);
			linkStat[i].del = 1;
			int paired_id = get_pair_id(i);
			linkStat[paired_id].del = 1;
			repeat_nodes.push_back(paired_id);
		}
		
	}
}

//remove small nodes whose length are less than half of the Insert_size
void remove_small_nodes(int contig_num, LinkStat *linkStat, vector<int> &small_nodes, vector<string> &contig_seqs)
{
	for (int i=1; i<contig_num; i++)
	{	
		//also remove the small fragments, which may cause interleaving in later steps
		if ( (linkStat[i].del == 0) && (i % 2 == 1) && (contig_seqs[i].size() < InsertSize / 2)  )
		{	
			small_nodes.push_back(i);
			linkStat[i].del = 1;
			int paired_id = get_pair_id(i);
			linkStat[paired_id].del = 1;
			small_nodes.push_back(paired_id);
		}
	}
}


//remove all the links that related with repeat nodes
void remove_links_from_deleted_nodes(CtgLink **ctgLink, int contig_num, LinkStat *linkStat)
{	
	for ( int i=0; i<contig_num; i++)
	{	
		if (ctgLink[i] == NULL)
		{	continue;	
		}

		CtgLink *p = ctgLink[i];

		while (1)
		{	
			int source_node = i;
			int target_node = p->id;

			if (linkStat[source_node].del == 1 || linkStat[target_node].del == 1)
			{	
				
				p->id =0;
				p->freq = 0;
				p->size = 0;
				Deleted_link_num ++;

				if (linkStat[source_node].link > 0)
				{	linkStat[source_node].link -- ;
				}
				if (linkStat[target_node].inlink > 0)
				{	linkStat[target_node].inlink -- ;
				}
			}

			if (p->pointer == NULL)
			{	break;
			}
			p = p->pointer;
		}
		
	}
}



//make a string with specified number of Ns
void generate_Nstr(int num, string &Nstr)
{	
	for (int i=0; i<num; i++)
	{	Nstr.push_back('N');
	}
}


//get the set of linear nodes from a given start node
void get_linear_seq( CtgLink **ctgLink,  LinkStat *linkStat, int start_id, vector<int> &vec_id_and_gap )
{
	int next_id = start_id;
	int gap_size = 0;
	while(1)
	{
		next_id = get_next_linked_id( ctgLink, next_id, gap_size );
		
		if(linkStat[next_id].del != 1)
		{	vec_id_and_gap.push_back( gap_size );
			vec_id_and_gap.push_back(next_id);
		}else
		{	break;
		}
		
		linkStat[next_id].del = 1;
		linkStat[get_pair_id(next_id)].del = 1;
		
		if(linkStat[next_id].link != 1)
		{	break;
		}
		
	}
}


//get the next linked contig id for current contig node [.link == 1]
int get_next_linked_id( CtgLink **ctgLink, int source_id, int &gap_size )
{
	int next_id = 0;
	CtgLink *p = ctgLink[source_id];
	while (p != NULL) 
	{	if (p->freq > 0) 
		{	next_id = p->id;
			gap_size = p->size / p->freq;
			break;
		}
		p = p->pointer;
	}
	
	return next_id;
}




//split scaffold into contigs and stored in vector
void scaffold_to_contig(string &scaffold, vector <string> &contigs, vector <int> &starts)
{
	int i = 0;
	while (i < scaffold.size())
	{
		while ( i < scaffold.size() && scaffold[i] == 'N' )
		{	i++;
		} // pass the N regions
		
		string contig_str;
		while ( i < scaffold.size() && scaffold[i] != 'N' )
		{	
			contig_str.push_back(scaffold[i]);			
			i++;
		} // get one continuous base region

		if (contig_str.size() > 0)
		{	contigs.push_back(contig_str);
			starts.push_back( i - contig_str.size() );
		}
	}
}


//release the memory of ctgLink
void free_ctgLink_memory(CtgLink **ctgLink, int contig_num)
{	
	for ( int i=0; i<contig_num; i++)
	{	
		if (ctgLink[i] == NULL)
		{	continue;	
		}
		CtgLink *p = ctgLink[i];
		while (1)
		{	CtgLink *p_next = p->pointer;
			delete p;
			if (p_next == NULL)
			{	break;
			}
			p = p_next;
		}
	}
	delete [] ctgLink;

}


//function similar to perl's routine split()
void split(string &strLine, vector<string>& tokens, const char* delim) 
{
	int count = 0;
	for(;;) {
		//erase delimiter
	   int i = strLine.find_first_not_of(delim);
	   if(i == -1)
	    	break;
	   strLine.erase(0, i);
	   
	   i = strLine.find_first_of(delim);
	   if(i == -1) {
		    tokens.push_back(strLine);
		    break;
	   } else {
		    string token = strLine.substr(0, i);
		    strLine.erase(0, i);
		    tokens.push_back(token);
	   }
	}
}




