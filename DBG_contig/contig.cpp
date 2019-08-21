/* 
This file include variables and routines to prune the kmer debruijn graph, and read out the contig sequence and depth information
Author: Wei Fan, fanweiagis@126.com
Date: 2015-11-4

Note: either tandem repeat and dispersed repeat will resuilt in branching nodes, so contig can be read out by extending linear nodes and broken at the branching nodes

Altert: If there are two neighbor kmers, both of them are branching, the current algorithm will miss processing them 

*/


#include "contig.h"

using namespace std;


KmerLink *klink;                        //store the link related information for each kmer node
int KmerFreqCutoff = 2;                 //kmer link with frequnency less than this cutoff will be deleted
vector <uint64_t> tip_nodes;            //store the tip nodes, used to accelarte speed 
vector <uint64_t> branch_nodes;         //store the branching nodes, used to accelarte speed

int Tip_len_cutoff = 100;               //the max allowed length for a tip
double Tip_depth_cutoff = 3.0;          //the max allowed depth for a tip
int Bubble_len_cutoff = 100;            //the max allowed length for a bubble
double Bubble_len_diff_rate_cutoff  = 0.1;    //the max allowed length difference rate between the two edges in a bubble
double Bubble_base_diff_rate_cutoff = 0.1;   //the max allowed base difference rate between the two edges in a bubble
int LowCovEdge_len_cutoff = 100;        //the max allowed length for small and low coverage edges between two branching kmer-nodes
double LowCovEdge_depth_cutoff = 3.0;   //the max allowed dpeth for small and low coverage edges between two branching kmer-nodes

uint32_t BitMaskVal[4] = {0xFFFFFFu, 0xFF00FFFFu, 0xFFFF00FFu, 0xFFFFFF00u};  //used to clear the link-depth of specific base, by bit & process

int is_remove_tip = 1;              // 1: do it;  0: not do it;
int is_remove_bubble = 1;           // 1: do it;  0: not do it;
int is_remove_lowedge = 1;          // 1: do it;  0: not do it;

int Contig_len_cutoff = 125;


//used for sorting the gap-seqs by its length
typedef struct{
        uint64_t len;
        string seq;
		string depth;
}LenAndSeq;

//used for sorting the gap-seqs by its length
bool cmpSeqByLen (LenAndSeq a, LenAndSeq b)
{       return b.len < a.len;
}


//prune the kmer debruijn graph, and read out the contig sequence
void build_contig_sequence()
{

	//calculate the link status for each kmer node
	cerr << "\nStart to calulate kmer links information!" << endl;
	klink = (KmerLink *)malloc(kset->size * sizeof(KmerLink));
	memset_parallel( klink, 0, kset->size * sizeof(KmerLink), threadNum );

	calculate_kmer_links(kset, klink, tip_nodes, branch_nodes);
	time_end = clock();
	cerr << "Finished! Run time: " << double(time_end - time_start) / CLOCKS_PER_SEC << endl;


	//remove error tips
	if (is_remove_tip)
	{	cerr << "\nStart to remove tips caused by sequencing error!" << endl;
		remove_error_tips(kset, klink, tip_nodes, Tip_len_cutoff, Tip_depth_cutoff);
		tip_nodes.clear();
		time_end = clock();
		cerr << "Finished! Run time: " << double(time_end - time_start) / CLOCKS_PER_SEC << endl;
	}

	
	//remove small lowCoverage edges between two branching nodes
	if (is_remove_lowedge)
	{	cerr << "\nStart to remove small low coverage edges between two branching nodes!" << endl;
		remove_lowCov_edges(kset, klink,  branch_nodes, LowCovEdge_len_cutoff, LowCovEdge_depth_cutoff);
		time_end = clock();
		cerr << "Finshed! Run time: " << double(time_end - time_start) / CLOCKS_PER_SEC << endl;
	}
 
	
	//remove bubbles mainly caused by large repeats and heterozygous SNPs and indels
	if (is_remove_bubble)
	{	cerr << "\nStart to remove bubbles caused by repeats and heterozygotes!" << endl;
		remove_hetero_bubbles(kset, klink,  branch_nodes, Bubble_len_cutoff, Bubble_len_diff_rate_cutoff, Bubble_base_diff_rate_cutoff);
		time_end = clock();
		branch_nodes.clear();
		cerr << "Finished! Run time: " << double(time_end - time_start) / CLOCKS_PER_SEC << endl;
	}

	
	//read out the contig sequence
	cerr << "\nStart to read out contig sequence and the depth information!" << endl;
	read_out_contig(kset, klink);
	time_end = clock();
	cerr << "Finished! Run time: " << double(time_end - time_start) / CLOCKS_PER_SEC << endl;
	
}


//calculate the leftward and rightward link number for each kmer, store the data in klink for later usage
//also stat the depth of kmer links, and get the linear, tip and bubble nodes
void calculate_kmer_links(KmerSet *kset, KmerLink *klink, vector <uint64_t> &tip_nodes, vector <uint64_t> &branch_nodes)
{	
	int64_t total_kmer_speceis_num = 0;
	int64_t deleted_lowFreq_kmer_num = 0;
	int64_t linear_kmer_node_num = 0;
	int64_t DepthStat[256];
	for (int i=0; i<256; i++)
	{	DepthStat[i] = 0;
	}

	KmerNode *array = kset->array;

	for (uint64_t i=0; i<kset->size; i++)
	{	
		if (! is_entity_null(kset->nul_flag, i))
		{	
			//klink[i].l_link_num = 0;
			//klink[i].l_link_base = 0;  //store the base with max depth
			//klink[i].linear = 0;
			int max_kmer_depth = 0;
				
			for (uint8_t j=0; j<4; j++)
			{	int link_base_depth = get_next_kmer_depth(array[i].l_link, j);
				DepthStat[link_base_depth] ++;

				if ( link_base_depth > KmerFreqCutoff )
				{	if (klink[i].l_link_num < 3)
					{	klink[i].l_link_num ++;
					}
					if(max_kmer_depth < link_base_depth) 
					{	max_kmer_depth = link_base_depth; 
						klink[i].l_link_base = j; 
					}
				}
			}
			
			//klink[i].r_link_num = 0;
			//klink[i].r_link_base = 0;  //store the base with max depth
			max_kmer_depth = 0;
				
			for (uint8_t j=0; j<4; j++)
			{	int link_base_depth = get_next_kmer_depth(array[i].r_link, j);
				DepthStat[link_base_depth] ++;

				if ( link_base_depth > KmerFreqCutoff )
				{	if (klink[i].r_link_num < 3)
					{	klink[i].r_link_num ++;
					}
					if(max_kmer_depth < link_base_depth) 
					{	max_kmer_depth = link_base_depth; 
						klink[i].r_link_base = j; 
					}
				}
			}
			
			total_kmer_speceis_num ++;

			if (klink[i].l_link_num == 0 && klink[i].r_link_num == 0)
			{	set_entity_delete(kset->del_flag, i);
				deleted_lowFreq_kmer_num ++;
			}

			if (klink[i].l_link_num == 1 && klink[i].r_link_num == 1)
			{	klink[i].linear = 1;
				linear_kmer_node_num ++;
			}
			if (klink[i].l_link_num + klink[i].r_link_num == 1)
			{	tip_nodes.push_back(i);
			}
			if (klink[i].l_link_num > 1 || klink[i].r_link_num > 1)
			{	branch_nodes.push_back(i);
			}

		}
	}



	//output the kmer links data
	string kmer_depth_file = Output_prefix + ".contig.kmer.freq";
	ofstream depthFile (kmer_depth_file.c_str());
	if ( ! depthFile )
	{	cerr << "fail to open file " << kmer_depth_file << endl;
	}
	
	cerr << "\nTotal kmer nodes number:    " << total_kmer_speceis_num << endl;
	cerr << "Deleted lowfreq kmer nodes: " << deleted_lowFreq_kmer_num << "\t"<< (double)deleted_lowFreq_kmer_num/total_kmer_speceis_num << endl;
	cerr << "Used linear kmer nodes:     " << linear_kmer_node_num << "\t"<< (double)linear_kmer_node_num/total_kmer_speceis_num << endl;
	cerr << "Used tip kmer nodes:        " << tip_nodes.size() << "\t"<< (double)tip_nodes.size()/total_kmer_speceis_num << endl;
	cerr << "Used branching kmer nodes:  " << branch_nodes.size() << "\t"<< (double)branch_nodes.size()/total_kmer_speceis_num << endl;


	depthFile << "Kmer_depth\tAppear_times\n";
	for(int i=1; i<=255; i++)
	{	depthFile << i << "\t" << DepthStat[i] << endl;
	}
	depthFile.close();

}



//re-calculate klinks information and re-set the l_link and r_link for a given kmer node, used after the tip/bubble/edge deleting process
void recalculate_kmer_links(KmerSet *kset, KmerLink *klink, uint64_t idx)
{
	if (idx == kset->size)
	{	return;
	}

	KmerNode *array = kset->array;
	
	klink[idx].l_link_num = 0;
	klink[idx].l_link_base = 0;  //store the base with max depth
	klink[idx].linear = 0;
	int max_kmer_depth = 0;

	for (uint8_t j=0; j<4; j++)
	{	int link_base_depth = get_next_kmer_depth(array[idx].l_link, j);
		if ( link_base_depth > KmerFreqCutoff )
		{	
			uint64_t next_kmer, next_kmer_rc, next_kmer_final;
			next_kmer = get_next_kmer_leftward(array[idx].kmer, j);
			next_kmer_rc = get_rev_com_kbit(next_kmer, KmerSize);
			next_kmer_final = (next_kmer < next_kmer_rc) ? next_kmer : next_kmer_rc;
			if( exist_kmerset(kset,next_kmer_final) != kset->size )
			{	if (klink[idx].l_link_num < 3)
				{	klink[idx].l_link_num ++;
				}
				if(max_kmer_depth < link_base_depth) 
				{	max_kmer_depth = link_base_depth; 
					klink[idx].l_link_base = j; 
				}
			}else{
				array[idx].l_link &= BitMaskVal[j];  //clear the link-depth of specific base 
			}
		}
	}
	
	klink[idx].r_link_num = 0;
	klink[idx].r_link_base = 0;  //store the base with max depth
	max_kmer_depth = 0;
		
	for (uint8_t j=0; j<4; j++)
	{	int link_base_depth = get_next_kmer_depth(array[idx].r_link, j);
		if ( link_base_depth > KmerFreqCutoff )
		{	
			uint64_t next_kmer, next_kmer_rc, next_kmer_final;
			next_kmer = get_next_kmer_rightward(array[idx].kmer, j);
			next_kmer_rc = get_rev_com_kbit(next_kmer, KmerSize);
			next_kmer_final = (next_kmer < next_kmer_rc) ? next_kmer : next_kmer_rc;
			
			if( exist_kmerset(kset,next_kmer_final) != kset->size )
			{
				if (klink[idx].r_link_num < 3)
				{	klink[idx].r_link_num ++;
				}
				if(max_kmer_depth < link_base_depth) 
				{	max_kmer_depth = link_base_depth; 
					klink[idx].r_link_base = j; 
				}
			}else{
				array[idx].r_link &= BitMaskVal[j];  //clear the link-depth of specific base 
			}
		}
	}
	
	if (klink[idx].l_link_num == 1 && klink[idx].r_link_num == 1)
	{	klink[idx].linear = 1;
	}
	
}

//remove dead-end tips on the kmer de bruijn graph, tips are caused by sequencing errors 
//search tip from the break-end kmer node to the branching kmer
void remove_error_tips(KmerSet *kset, KmerLink *klink,  vector <uint64_t> &tip_nodes, int tip_len_cutoff, double tip_depth_cutoff)
{	
	uint64_t total_tip_num = 0;
	uint64_t total_tip_len = 0;

	//output the kmer tips
	string tip_file = Output_prefix + ".contig.tip.fa";
	ofstream tipFile (tip_file.c_str());
	if ( ! tipFile )
	{	cerr << "fail to open file " << tip_file << endl;
	}

	KmerNode *array = kset->array;
	for (uint64_t i=0; i<tip_nodes.size(); i++)
	{	uint64_t idx = tip_nodes[i];
		int walk_direct = (klink[idx].l_link_num == 1) ? -1 : 1;  //direction -1: left; 1: right; because only one side of the node has a branch
		int tip_len = 0;
		int tip_depth = 0;
		string tip_str;
		vector <uint64_t> vec_idx;
		uint64_t last_idx; //the break point idx
		string break_or_branch;
		//the tip kmer nodes has already been checked (one end is break-point, the other end is linear)
		get_linear_path(idx, walk_direct, tip_len_cutoff, tip_len, tip_depth, vec_idx, tip_str, last_idx, break_or_branch);
		double tip_depth_avg = (double)tip_depth / tip_len;

		//judge a tip, delete tip nodes, and recover the links of last kmer
		if (tip_depth_avg <= tip_depth_cutoff && tip_len <= tip_len_cutoff)
		{	
			total_tip_num ++;
			total_tip_len += tip_len;

			for (int i = 0; i < vec_idx.size(); i++)
			{	set_entity_delete(kset->del_flag, vec_idx[i]);
			}

			recalculate_kmer_links(kset, klink, last_idx);
			klink[last_idx].in_tip = 1;
			
			uint64_t leftendkmer, rightendkmer;
			string leftendmark, rightendmark;
			if (walk_direct == 1)
			{	leftendkmer = array[idx].kmer;
				leftendmark = "break";
				rightendkmer = array[last_idx].kmer;
				rightendmark = break_or_branch;
			}else
			{	rightendkmer = array[idx].kmer;
				rightendmark = "break";
				leftendkmer = array[last_idx].kmer;
				leftendmark = break_or_branch;
			}


			string kmer_str = bit2seq(array[idx].kmer, KmerSize);
			string out_str;
			if (walk_direct == 1)
			{	out_str = kmer_str + tip_str;
			}else
			{	reverse( tip_str.begin(), tip_str.end() );
				out_str = tip_str + kmer_str;
			}

			tipFile << ">tip_" <<  total_tip_num << "\tlength: " << tip_len + KmerSize << "\tavgDepth: " << tip_depth_avg << "\tLeftEndKmer: " << leftendkmer << " " << leftendmark << "\tRightEndKmer: " << rightendkmer << " " << rightendmark << "\n"  << out_str << "\n";

		}


	}

	tipFile.close();
	
	cerr << "\nremove total tip number:  " << total_tip_num << endl;
	cerr << "remove total tip length:  " << total_tip_len << endl;
}




//get the bases and depths from a given link data-structure
void get_branch_bases(uint32_t link, vector <uint8_t> &vec_bases, vector <uint8_t> &vec_depths)
{	
	for (uint8_t j=0; j<4; j++)
	{	int link_base_depth = get_next_kmer_depth(link, j);
		if ( link_base_depth > KmerFreqCutoff )
		{	vec_bases.push_back(j);
			vec_depths.push_back(link_base_depth);
		}
	}
}

//used after removing tips and lowedges, the purpose is to remove bi-bubbles caused by heterozygous SNPs and indels as well as repeats
//this routine can also remove some low-depth bubble edges caused by sequencing errors, but these errors have been processed by the remove_lowedge function
//only process simple bubbles, which do not have other branches on each bubble edges, after removing tip and lowedges, the graph will get much clearer
void remove_hetero_bubbles(KmerSet *kset, KmerLink *klink,  vector <uint64_t> &branch_nodes, int bubble_len_cutoff, double bubble_len_diff_rate_cutoff, double bubble_base_diff_rate_cutoff)
{
	
	//output the kmer bubbles
	string bubble_file = Output_prefix + ".contig.bubble.fa";
	ofstream bubbleFile (bubble_file.c_str());
	if ( ! bubbleFile )
	{	cerr << "fail to open file " << bubble_file << endl;
	}

	uint64_t total_bubble_num = 0;
	uint64_t total_bubble_len = 0;

	KmerNode *array = kset->array;
	for (uint64_t i=0; i<branch_nodes.size(); i++)
	{	
		//set information of current branching kmer node
		uint64_t idx = branch_nodes[i];
		

		int walk_direct = 0;
		vector <uint8_t> vec_bases;
		vector <uint8_t> vec_depths;

		//the conditions for both bubble and tiny repeats (uint length less than Kmer length)
		if (klink[idx].l_link_num == 2 && klink[idx].r_link_num == 1)
		{	walk_direct = -1;
			get_branch_bases(array[idx].l_link, vec_bases, vec_depths);
		}
		else if (klink[idx].l_link_num == 1 && klink[idx].r_link_num == 2)
		{	walk_direct = 1;
			get_branch_bases(array[idx].r_link, vec_bases, vec_depths);
		}
		else
		{	continue;
		}
		

		//get the next 2 branched kmer nodes
		uint64_t next_kmer1, next_kmer1_rc, next_kmer1_final;
		uint64_t next_kmer2, next_kmer2_rc, next_kmer2_final;
		int walk_direct1, walk_direct2;
		if (walk_direct == 1)
		{	next_kmer1 = get_next_kmer_rightward(array[idx].kmer, vec_bases[0]);
			next_kmer2 = get_next_kmer_rightward(array[idx].kmer, vec_bases[1]);
		}else{
			next_kmer1 = get_next_kmer_leftward(array[idx].kmer, vec_bases[0]);
			next_kmer2 = get_next_kmer_leftward(array[idx].kmer, vec_bases[1]);
		}
		next_kmer1_rc = get_rev_com_kbit(next_kmer1, KmerSize);
		next_kmer2_rc = get_rev_com_kbit(next_kmer2, KmerSize);
		
		if (next_kmer1 < next_kmer1_rc)
		{	next_kmer1_final = next_kmer1;
			walk_direct1 = walk_direct;
		}else
		{	next_kmer1_final = next_kmer1_rc;
			walk_direct1 = -walk_direct;
		}

		if (next_kmer2 < next_kmer2_rc)
		{	next_kmer2_final = next_kmer2;
			walk_direct2 = walk_direct;
		}else
		{	next_kmer2_final = next_kmer2_rc;
			walk_direct2 = -walk_direct;
		}

		uint64_t idx1 = exist_kmerset(kset, next_kmer1_final);
		uint64_t idx2 = exist_kmerset(kset, next_kmer2_final);

		//get the two extending paths
		int bubble_len1, bubble_len2;
		int bubble_depth1, bubble_depth2;
		vector <uint64_t> vec_idx1, vec_idx2;
		string bubble_str1, bubble_str2;
		string align_bubble_str1, align_bubble_str2;
		uint64_t last_idx1, last_idx2;
		int diff_base_num = 0;
		double diff_base_rate = 0;
		string break_or_branch1, break_or_branch2;
		string bubble_type;

		if(klink[idx1].linear != 1 || klink[idx2].linear != 1)
		{	continue;
		}
		

		get_linear_path(idx1, walk_direct1, bubble_len_cutoff, bubble_len1, bubble_depth1, vec_idx1, bubble_str1, last_idx1, break_or_branch1);
		get_linear_path(idx2, walk_direct2, bubble_len_cutoff, bubble_len2, bubble_depth2, vec_idx2, bubble_str2, last_idx2, break_or_branch2);

		
		double bubble_depth1_avg = (double)bubble_depth1 / bubble_len1;
		double bubble_depth2_avg = (double)bubble_depth2 / bubble_len2; 

		if(last_idx1 != last_idx2)
		{	if (  bubble_depth1_avg >  LowCovEdge_depth_cutoff  &&  bubble_depth2_avg > LowCovEdge_depth_cutoff  )
			{	klink[idx].in_repeat = 1;   //not bubble structure, but is tiny-repeat branching structure
			}
			continue;
		}

		
		string kmer_str1 = bit2seq(array[idx1].kmer, KmerSize);
		if (walk_direct1 == 1)
		{	bubble_str1 = kmer_str1 + bubble_str1;
		}else
		{	reverse( bubble_str1.begin(), bubble_str1.end() );
			bubble_str1 = bubble_str1 + kmer_str1;
		}
		
		string kmer_str2 = bit2seq(array[idx2].kmer, KmerSize);
		if (walk_direct2 == 1)
		{	bubble_str2 = kmer_str2 + bubble_str2;
		}else
		{	reverse( bubble_str2.begin(), bubble_str2.end() );
			bubble_str2 = bubble_str2 + kmer_str2;
		}
		
		if ( walk_direct1 != walk_direct2 )
		{	reverse( bubble_str1.begin(), bubble_str1.end() );
			complement_sequence(bubble_str1);
		}

		//Add the first character 
		bubble_len1 ++;
		bubble_len2 ++;
		bubble_depth1 += vec_depths[0];
		bubble_depth2 += vec_depths[1];
		
		//compare the two bubble paths
		if (bubble_len1 == bubble_len2) //consider only seq-error and hetero-SNP, but not consider hetero-indel at this moment
		{	diff_base_num = compare_two_seq_simple(bubble_str1, bubble_str2);
			
			diff_base_rate = (double)diff_base_num / bubble_len1;
			bubble_type = "SNP";
		}
		if(bubble_len1 != bubble_len2 || diff_base_rate > bubble_base_diff_rate_cutoff)  //whne no indel or two nearby indel happens
		{	
			int final_dp_score = 0;
			global_aligning(bubble_str1, bubble_str2, align_bubble_str1, align_bubble_str2, final_dp_score); //for different length aligning, use the dynamic programming
			bubble_str1 = align_bubble_str1;
			bubble_str2 = align_bubble_str2;

			diff_base_num = compare_two_seq_simple(bubble_str1, bubble_str2);
			 
			diff_base_rate = (double)diff_base_num / bubble_len1;

			bubble_type = "INDEL";
		}


		//remove the lower depth bubble edge
		if( diff_base_rate  <  bubble_base_diff_rate_cutoff && abs(bubble_len1 - bubble_len2) < bubble_len_cutoff * bubble_len_diff_rate_cutoff  && (bubble_len1 <= bubble_len_cutoff && bubble_len2 <= bubble_len_cutoff) )
		{	
			int removed_bubble_id = -1;
			if (bubble_depth1_avg < bubble_depth2_avg)
			{	for (int i = 0; i < vec_idx1.size(); i++)
			    {   set_entity_delete(kset->del_flag, vec_idx1[i]);
				}
				recalculate_kmer_links(kset, klink, last_idx1);
				recalculate_kmer_links(kset, klink, idx);
				total_bubble_num ++;
				total_bubble_len += bubble_len1;
				removed_bubble_id = 1;
				
			}else
			{	for (int i = 0; i < vec_idx2.size(); i++) 
				{	set_entity_delete(kset->del_flag, vec_idx2[i]);
				}
				recalculate_kmer_links(kset, klink, last_idx2);
				recalculate_kmer_links(kset, klink, idx);
				total_bubble_num ++;
				total_bubble_len += bubble_len2;
				removed_bubble_id = 2;

			}


			uint64_t leftendkmer, rightendkmer;
			string leftendmark, rightendmark;
			if (walk_direct == 1)
			{	leftendkmer = array[idx].kmer;
				leftendmark = "branch";
				rightendkmer = array[last_idx1].kmer;
				rightendmark = break_or_branch1;
			}else
			{	rightendkmer = array[idx].kmer;
				rightendmark = "branch";
				leftendkmer = array[last_idx1].kmer;
				leftendmark = break_or_branch1;
			}

			bubbleFile << ">bubble_" << total_bubble_num << "\ttype: " << bubble_type  << "\tlength1: " << bubble_len1  + KmerSize << "\tavgDepth1: " << bubble_depth1_avg << "\tlength2: " << bubble_len2 + KmerSize  << "\tavgDepth2: " << bubble_depth2_avg << "\tremoved: " << removed_bubble_id  <<     "\tLeftEndKmer: " << leftendkmer << " " << leftendmark << "\tRightEndKmer: " << rightendkmer << " " << rightendmark  << "\n" << bubble_str1 << "\n" << bubble_str2 << "\n";


			klink[idx].in_bubble = 1;         //start point of the bubble
			klink[last_idx1].in_bubble = 1;   //end point of the bubble, last_idx1 is equal to last_idx2, so only used one is enough

		}	
	}

	cerr << "\nremove total bubble number: " << total_bubble_num << endl;
	cerr << "remove total bubble length: " << total_bubble_len << endl;

	bubbleFile.close();

}



//compare two bubble sequences with equal length
int compare_two_seq_simple(string &str1, string &str2)
{	int diff_base_num = 0;
	for (int i=0; i<str1.size(); i++)
	{	if (str1[i] != str2[i] && str1[i] != '-' && str2[i] != '-')  //indel position does not count difference, because indel length is much bigger than SNP
		{	diff_base_num ++;
		}
	}
	return diff_base_num;
}



//Further remove small and low-depth edges, caused by sequencing-errors, including bi-bubbles and mis-link structures
//In fact, this function can remove all the sequencing-error caused tips, bi-bubbles and mis-link structures  
void remove_lowCov_edges(KmerSet *kset, KmerLink *klink,  vector <uint64_t> &branch_nodes, int edge_len_cutoff, double edge_depth_cutoff)
{	
	int lowCovEdge_num = 0;
	int lowCovEdge_len =0;
	
	//output the small low depth edges
	string lowedge_file = Output_prefix + ".contig.lowedge.fa";
	ofstream lowedgeFile (lowedge_file.c_str());
	if ( ! lowedgeFile )
	{	cerr << "fail to open file " << lowedge_file << endl;
	}

	KmerNode *array = kset->array;
	for (uint64_t i=0; i<branch_nodes.size(); i++)
	{	uint64_t idx = branch_nodes[i];
		

		if(klink[idx].r_link_num >= 2)
		{
			//search rightward, get the first branching bases
			int walk_direct = 1;
			vector <uint8_t> vec_bases;
			vector <uint8_t> vec_depths;
			get_branch_bases(array[idx].r_link, vec_bases, vec_depths);
			
			
			for(int j=0; j<vec_bases.size(); j++)
			{	uint64_t next_kmer, next_kmer_rc, next_kmer_final;
				next_kmer = get_next_kmer_rightward(array[idx].kmer, vec_bases[j]);
				next_kmer_rc = get_rev_com_kbit(next_kmer, KmerSize);
				int walk_direct1 = 0;
				if (next_kmer < next_kmer_rc)
				{	next_kmer_final = next_kmer;
					walk_direct1 = walk_direct;
				}else{
					next_kmer_final = next_kmer_rc;
					walk_direct1 = -walk_direct;
				}
				uint64_t idx1 = exist_kmerset(kset, next_kmer_final);
				int edge_len1 =  0;
				int edge_depth1 = 0;
				vector <uint64_t> vec_idx1;
				string edge_str1;
				uint64_t last_idx1;
				string break_or_branch1;
				

				if(klink[idx1].linear != 1)
				{	continue;
				}				
				get_linear_path(idx1, walk_direct1, edge_len_cutoff, edge_len1, edge_depth1, vec_idx1, edge_str1, last_idx1, break_or_branch1);
				

				edge_len1 ++;
				edge_depth1 += vec_depths[j];
				double edge_depth1_avg = (double)edge_depth1 / edge_len1;

				if( edge_len1 <= edge_len_cutoff && edge_depth1_avg <= edge_depth_cutoff && klink[last_idx1].linear != 1)
				{	
					lowCovEdge_num ++;
					lowCovEdge_len += edge_len1;
					for (int i = 0; i < vec_idx1.size(); i++)
					{	set_entity_delete(kset->del_flag, vec_idx1[i]);
					}
					recalculate_kmer_links(kset, klink, last_idx1);
					recalculate_kmer_links(kset, klink, idx);
					klink[idx].in_lowedge = 1;
					klink[last_idx1].in_lowedge = 1;


					uint64_t leftendkmer, rightendkmer;
					string leftendmark, rightendmark;
					leftendkmer = array[idx].kmer;
					leftendmark = "branch";
					rightendkmer = array[last_idx1].kmer;
					rightendmark = break_or_branch1;	

					string kmer_str1 = bit2seq(array[idx1].kmer, KmerSize);
					string out_str1;
					if (walk_direct1 == 1)
					{	out_str1 = kmer_str1 + edge_str1;
					}else
					{	reverse( edge_str1.begin(), edge_str1.end() );
						out_str1 = edge_str1 + kmer_str1;
					}				

					lowedgeFile << ">lowedge_" << lowCovEdge_num << "\tlength: " << edge_len1  + KmerSize << "\tavgDepth: " << edge_depth1_avg  << "\tLeftEndKmer: " << leftendkmer << " " << leftendmark << "\tRightEndKmer: " << rightendkmer << " " << rightendmark   << "\n" << out_str1 << "\n";


				}
			}
		}
       
		 if(klink[idx].l_link_num >= 2)
         {
             //search leftward
             int walk_direct = -1;
             vector <uint8_t> vec_bases;
             vector <uint8_t> vec_depths;
             get_branch_bases(array[idx].l_link, vec_bases, vec_depths);

 
             for(int j=0; j<vec_bases.size(); j++)
             {   uint64_t next_kmer, next_kmer_rc, next_kmer_final;
                 next_kmer = get_next_kmer_leftward(array[idx].kmer, vec_bases[j]);
                 next_kmer_rc = get_rev_com_kbit(next_kmer, KmerSize);
                 int walk_direct1 = 0;
                 if (next_kmer < next_kmer_rc)
                 {   next_kmer_final = next_kmer;
                     walk_direct1 = walk_direct;
                 }else{
                     next_kmer_final = next_kmer_rc;
                     walk_direct1 = -walk_direct;
                 }
                 uint64_t idx1 = exist_kmerset(kset, next_kmer_final);
                 int edge_len1 =  0;
                 int edge_depth1 = 0;
                 vector <uint64_t> vec_idx1;
                 string edge_str1;
                 uint64_t last_idx1;
				 string break_or_branch1;

				
				 if(klink[idx1].linear != 1)
				 {	continue;
				 }	
                 get_linear_path(idx1, walk_direct1, edge_len_cutoff, edge_len1, edge_depth1, vec_idx1, edge_str1, last_idx1, break_or_branch1);

				 
                 edge_len1 ++;
                 edge_depth1 += vec_depths[j];
				 double edge_depth1_avg = (double)edge_depth1 / edge_len1;

                 if( edge_len1 <= edge_len_cutoff && edge_depth1_avg <= edge_depth_cutoff && klink[last_idx1].linear != 1)
                 {   
					 lowCovEdge_num ++;
				 	 lowCovEdge_len += edge_len1;
					 for (int i = 0; i < vec_idx1.size(); i++)
                     {   set_entity_delete(kset->del_flag, vec_idx1[i]);
                     }
                     recalculate_kmer_links(kset, klink, last_idx1);
                     recalculate_kmer_links(kset, klink, idx);
					 klink[idx].in_lowedge = 1;
					 klink[last_idx1].in_lowedge = 1;


					uint64_t leftendkmer, rightendkmer;
					string leftendmark, rightendmark;
					rightendkmer = array[idx].kmer;
					rightendmark = "branch";
					leftendkmer = array[last_idx1].kmer;
					leftendmark = break_or_branch1;	

					string kmer_str1 = bit2seq(array[idx1].kmer, KmerSize);
					string out_str1;
					 if (walk_direct1 == 1)
					 {	out_str1 = kmer_str1 + edge_str1;
					 }else
					 {	reverse( edge_str1.begin(), edge_str1.end() );
					 	out_str1 = edge_str1 + kmer_str1;
					 }

					 lowedgeFile << ">lowedge_" << lowCovEdge_num << "    length:" << edge_len1  + KmerSize << "    avgDepth:" << edge_depth1_avg  << "\tLeftEndKmer: " << leftendkmer << " " << leftendmark << "\tRightEndKmer: " << rightendkmer << " " << rightendmark  << "\n" << out_str1 << "\n";
 
                 }
             }
         }

	}

	cerr << "\nremove total lowCovEdge number: " << lowCovEdge_num  << endl;
	cerr << "remove total lowCovEdge length: " << lowCovEdge_len  << endl;
	
	lowedgeFile.close();

}

//get all available information on a linear path, which starts from a linear or break-end kmer node, and ending at a branching kmer node
void get_linear_path(uint64_t idx, int walk_direct, int len_cutoff, int &path_len, int &path_depth, vector <uint64_t> &vec_idx, string &path_str, uint64_t &last_idx, string &break_or_branch)
{	
	KmerNode *array = kset->array;
	int original_walk_direct = walk_direct;
	break_or_branch = "linear";   //1: linear;  -1 for break;  2 for branch;
	
	path_len = 0;
	path_depth = 0;
	while (1)
	{	path_len ++;
		vec_idx.push_back(idx);

		uint64_t next_kmer, next_kmer_rc, next_kmer_final;
		if (walk_direct == 1)
		{	next_kmer = get_next_kmer_rightward(array[idx].kmer, klink[idx].r_link_base);
			path_depth += get_next_kmer_depth(array[idx].r_link, klink[idx].r_link_base);
			path_str.push_back( (original_walk_direct == 1) ? bases[klink[idx].r_link_base] : c_bases[klink[idx].r_link_base] );
		}else{
			next_kmer = get_next_kmer_leftward(array[idx].kmer, klink[idx].l_link_base);
			path_depth += get_next_kmer_depth(array[idx].l_link, klink[idx].l_link_base);
			path_str.push_back( (original_walk_direct == 1) ? c_bases[klink[idx].l_link_base] : bases[klink[idx].l_link_base]  );
		}
		next_kmer_rc = get_rev_com_kbit(next_kmer, KmerSize);
		if (next_kmer < next_kmer_rc)
		{	next_kmer_final = next_kmer;
		}else
		{	next_kmer_final = next_kmer_rc;
			walk_direct = -walk_direct;
		}
		idx = exist_kmerset(kset, next_kmer_final); //exist require not empty and not deleted

		if (klink[idx].linear != 1 || idx == kset->size || path_len >= len_cutoff)
		{	
			last_idx = idx;
			
			if(idx == kset->size)
			{	break_or_branch = "break";
			}else
			{	if(klink[idx].l_link_num == 0 || klink[idx].r_link_num == 0)
				{	break_or_branch = "break";
				}else
				{	break_or_branch = "branch";
				}
			}

			break;
		}
	}
}

//used to get the contig sequence, start from a linear kmer node, and extending leftward and rightward seperately
//the delete mechanism can deal with circular path automatically
//use the depth information to qualify the accuracy of each base
void get_linear_seq(uint64_t idx, int walk_direct, int &seq_len, int &seq_depth,  string &seq_str, uint64_t &last_idx, string &break_or_branch, string &seq_depths, string &IsRepeat)
{	
	KmerNode *array = kset->array;
	int original_walk_direct = walk_direct;
	seq_len = 0;
	seq_depth = 0;
	last_idx = 0;
	break_or_branch = "linear"; //1: linear;  -1 for break;  2 for branch;

	while (1)
	{	seq_len ++;

		uint64_t next_kmer, next_kmer_rc, next_kmer_final;
		if (walk_direct == 1)
		{	next_kmer = get_next_kmer_rightward(array[idx].kmer, klink[idx].r_link_base);
			int this_depth = get_next_kmer_depth(array[idx].r_link, klink[idx].r_link_base);
			seq_depth += this_depth;
			if (this_depth == 10 || this_depth == 62) //the new line "\n" and fasta key character ">" should not be used
			{	this_depth --;
			}
			seq_depths.push_back((char)this_depth);
			seq_str.push_back( (original_walk_direct == 1) ? bases[klink[idx].r_link_base] : c_bases[klink[idx].r_link_base] );
		}else{
			next_kmer = get_next_kmer_leftward(array[idx].kmer, klink[idx].l_link_base);
			int this_depth  = get_next_kmer_depth(array[idx].l_link, klink[idx].l_link_base);
			seq_depth += this_depth;
			if (this_depth == 10 || this_depth == 62) //the new line "\n" and fasta key character ">" should not be used
			{	this_depth --;
			}
			seq_depths.push_back((char)this_depth);
			seq_str.push_back( (original_walk_direct == 1) ? c_bases[klink[idx].l_link_base] : bases[klink[idx].l_link_base] );
		}
		next_kmer_rc = get_rev_com_kbit(next_kmer, KmerSize);
		if (next_kmer < next_kmer_rc)
		{	next_kmer_final = next_kmer;
		}else
		{	next_kmer_final = next_kmer_rc;
			walk_direct = -walk_direct;
		}
		idx = exist_kmerset(kset, next_kmer_final);

		if (klink[idx].linear != 1 || idx == kset->size)
		{	last_idx = idx;
			if(idx == kset->size)
			{	break_or_branch = "break";
			}else
			{	if(klink[idx].l_link_num == 0 || klink[idx].r_link_num == 0)
				{	break_or_branch = "break";
				}else
				{	break_or_branch = "branch";
					if ( (walk_direct == 1 && klink[idx].r_link_num > 1) || (walk_direct == -1 && klink[idx].l_link_num > 1)  )
					{	IsRepeat = "Repeat";
					}else
					{	IsRepeat = "Unique";
					}

					//cerr << "#######\t" << walk_direct << "\t" << (int)klink[idx].l_link_num << "\t" << (int)klink[idx].r_link_num  << "\t" << IsRepeat << endl;
				}
			}
			break;
		}else{
			set_entity_delete(kset->del_flag, idx);
		}
	}
}

//read out all the contig sequences, which are defined as continuous non-branching kmer regions.
//start from any linear kmer nodes, extending to both sides, merge sequences from two sides to form the final contig
void read_out_contig(KmerSet *kset, KmerLink *klink)
{	
	KmerNode *array = kset->array;

	uint64_t break_point_num = 0;   //start from linear node, end in break-point or branching node
	uint64_t branch_point_num = 0;  //start from linear node, end in break-point or branching node

	uint64_t total_contig_num = 0;
	uint64_t total_contig_len = 0;
	uint64_t total_small_num = 0;
	uint64_t total_small_len = 0;

	string contig_str_file = Output_prefix + ".contig.seq.fa";
	string contig_depth_file = Output_prefix + ".contig.seq.depth";
	ofstream CtgFile (contig_str_file.c_str());
	ofstream CtgDepthFile (contig_depth_file.c_str());
	if(!CtgFile || !CtgDepthFile)
	{	cerr << "fail to open contig file " << contig_str_file << "\t" <<  contig_depth_file << endl;
	}

	string small_str_file = Output_prefix + ".contig.small.fa";
	string small_depth_file = Output_prefix + ".contig.small.depth";
	ofstream SmallFile (small_str_file.c_str());
	ofstream SmallDepthFile (small_depth_file.c_str());
	if(!SmallFile || !SmallDepthFile)
	{	cerr << "fail to open small file " << small_str_file << "\t" <<  small_depth_file << endl;
	}

	vector<LenAndSeq> CtgLenAndSeq;
	
	for (uint64_t i=0; i<kset->size; i++)
	{	if (! is_entity_null(kset->nul_flag, i) && ! is_entity_delete(kset->del_flag,i) && klink[i].linear == 1)
		{	
			//get one contig sequence from both sides, leftward and rightward;
			string contig_str;
			int contig_len = 0;
			double contig_depth = 0.0;
			string kmer_str = bit2seq(array[i].kmer, KmerSize);

			//以当前kmer为基点，向右遍历linear kmer nodes，每次循环得出一个idx和direction
			uint64_t idx = i;
			int walk_direct = 1;
			string contig_rightward;
			string contig_rightward_depths;
			int contig_rightward_len = 0;
			int contig_rightward_depth = 0;
			uint64_t rightward_last_idx = 0;
			string rightward_end_mark;   //break-point or branching node: -1 for break; 2 for branch;
			string rightward_IsRepeat = "Unknown";

			get_linear_seq(idx, walk_direct, contig_rightward_len, contig_rightward_depth,  contig_rightward, rightward_last_idx, rightward_end_mark, contig_rightward_depths, rightward_IsRepeat);
			
			//以当前kmer为基点，向左遍历linear kmer nodes，每次循环得出一个idx和direction
			idx = i;
		    walk_direct = -1;
			string contig_leftward;
			string contig_leftward_depths;
			int contig_leftward_len = 0;
			int contig_leftward_depth = 0;
			uint64_t leftward_last_idx = 0;
			string leftward_end_mark;      //break-point or branching node: -1 for break; 2 for branch;
			string leftward_IsRepeat = "Unknown";


			get_linear_seq(idx, walk_direct, contig_leftward_len, contig_leftward_depth,  contig_leftward, leftward_last_idx, leftward_end_mark, contig_leftward_depths, leftward_IsRepeat);
			
			//default is unique node, if both the leftward and rightward is repeat-like branching, then call it a repeat node
			string contig_type;
			if ( leftward_IsRepeat == "Repeat" && rightward_IsRepeat == "Repeat" )
			{	contig_type = "RepeatNode";
			}

			//merge the leftward and rightward contigs
			set_entity_delete(kset->del_flag, i);
			reverse( contig_leftward.begin(), contig_leftward.end() );
			reverse( contig_leftward_depths.begin(), contig_leftward_depths.end() );
			contig_str = contig_leftward + kmer_str + contig_rightward;
			contig_len = contig_leftward_len + KmerSize + contig_rightward_len;
			contig_depth = contig_leftward_depth + contig_rightward_depth;
			int contig_len_noKmer = contig_leftward_len + contig_rightward_len;
		   	contig_depth /= (double)contig_len_noKmer;
			
			string contig_middlekmer_depths;  //the bases of start kmer use average depth from the contig
			for(int j=0; j<KmerSize; j++)
			{	char this_depth = (char)contig_depth;
				if (this_depth == 10 || this_depth == 62) //the new line "\n" and fasta key character ">" should not be used
				{	this_depth --;
				}
				contig_middlekmer_depths.push_back(this_depth);
			}

			if(rightward_end_mark != "branch")
			{	break_point_num ++;
			}else
			{	branch_point_num ++;
			}
			if(leftward_end_mark!= "branch")
			{   break_point_num ++;
			}else
			{   branch_point_num ++;
			}
			
			///////////////////////////
			LenAndSeq one_ctg;
			
			one_ctg.len = contig_str.size();
			one_ctg.seq =  "\tlength: " +  boost::lexical_cast<string>(contig_len) + "\tavgDepth: " + boost::lexical_cast<string>(contig_depth) + "\tLeftEndKmer: " + boost::lexical_cast<string>(array[leftward_last_idx].kmer) + " " + leftward_end_mark + "-" + leftward_IsRepeat + "\tRightEndKmer: " + boost::lexical_cast<string>(array[rightward_last_idx].kmer) + " " + rightward_end_mark + "-" + rightward_IsRepeat + "\t" + contig_type + "\n"  + contig_str + "\n";
			one_ctg.depth = contig_leftward_depths + contig_middlekmer_depths + contig_rightward_depths;
			
			CtgLenAndSeq.push_back(one_ctg);
			
		}
	}
	
	sort(CtgLenAndSeq.begin(), CtgLenAndSeq.end(), cmpSeqByLen);
	
	uint64_t contig_id = 1;
	for (uint64_t i = 0; i < CtgLenAndSeq.size(); i ++ )
	{	if (CtgLenAndSeq[i].len >= Contig_len_cutoff)
		{	

			CtgFile << ">ctg_" + boost::lexical_cast<string>(contig_id) + CtgLenAndSeq[i].seq;				
			CtgDepthFile << ">ctg_" + boost::lexical_cast<string>(contig_id) + "\n" + CtgLenAndSeq[i].depth + "\n";
			total_contig_num ++;
			total_contig_len += CtgLenAndSeq[i].len ;

		}else
		{	
			SmallFile << ">ctg_" + boost::lexical_cast<string>(contig_id) + CtgLenAndSeq[i].seq;				
			SmallDepthFile << ">ctg_" + boost::lexical_cast<string>(contig_id) + "\n" + CtgLenAndSeq[i].depth + "\n";
			total_small_num ++;
			total_small_len += CtgLenAndSeq[i].len ;
		}
		contig_id += 2;
	}
	
	cerr << "\ncontig break-point number:     " << break_point_num << endl;
	cerr << "contig branch-point number:    " << branch_point_num << endl;

	cerr << "\nTotal contig number:   " << total_contig_num << endl;
	cerr << "Total contig length:   " << total_contig_len << endl;

	cerr << "\nTotal small edge number:   " << total_small_num << endl;
	cerr << "Total small edge length:   " << total_small_len << endl;


}



