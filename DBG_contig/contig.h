/* 
This file include variables and routines to prune the kmer debruijn graph, and read out the contig sequence and depth information
Author: Wei Fan, fanweiagis@126.com
Date: 2015-11-4

Note: either tandem repeat and dispersed repeat will resuilt in branching nodes, so contig can be read out by extending linear nodes and broken at the branching nodes
*/

#ifndef __CONTIG_H_
#define __CONTIG_H_

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <set>
#include <cmath>
#include "seqKmer.h"
#include "kmerSet.h"
#include "DBGgraph.h"
#include "global_aligning.h"
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>


using namespace std;


typedef struct {
    uint8_t l_link_num:2;           //number of branches in the leftward of kmer
	uint8_t l_link_base:2;          //the base with max frequency in the leftward of kmer
	uint8_t r_link_num:2;           //number of branches in the rightward of kmer
	uint8_t r_link_base:2;          //the base with max frequency in the rightward of kmer
	uint8_t linear:1;               //whether this kmer is linear, only one branch at each side
	uint8_t in_tip:1;
	uint8_t in_bubble:1;
	uint8_t in_lowedge:1;
	uint8_t in_repeat:1;
	//still left 3 bits not used
} KmerLink;


extern KmerLink *klink;             //store the link related information for each kmer node
extern int KmerFreqCutoff;          //kmer link with frequnency less than this cutoff will be deleted
extern vector <uint64_t> tip_nodes;        //store the tip nodes, used to accelarte speed 
extern vector <uint64_t> branch_nodes;     //store the branching nodes, used to accelarte speed

extern int Tip_len_cutoff;                 //the max allowed length for a tip
extern double Tip_depth_cutoff;            //the max allowed depth for a tip
extern int Bubble_len_cutoff;              //the max allowed length for a bubble
extern double Bubble_len_diff_rate_cutoff;      //the max allowed length difference between the two edges in a bubble
extern double Bubble_base_diff_rate_cutoff;     //the max allowed base difference rate between the two edges in a bubble
extern int LowCovEdge_len_cutoff;          //the max allowed length for small and low coverage edges between two branching kmer-nodes
extern double LowCovEdge_depth_cutoff;     //the max allowed dpeth for small and low coverage edges between two branching kmer-nodes
extern uint32_t BitMaskVal[4];                    //used to clear the link-depth of specific base, by bit & process

extern int is_remove_tip;                // 1: do it;  0: not do it;
extern int is_remove_bubble;             // 1: do it;  0: not do it;
extern int is_remove_lowedge;            // 1: do it;  0: not do it;

extern int Contig_len_cutoff;


//prune the kmer debruijn graph, and read out the contig sequence
void build_contig_sequence();


//calculate the leftward and rightward link number for each kmer, store the data in klink for later usage
//also stat the depth of kmer links, and get the linear, tip and bubble nodes
void calculate_kmer_links(KmerSet *kset, KmerLink *klink, vector <uint64_t> &tip_nodes, vector <uint64_t> &branch_nodes);


//re-calculate klinks information and re-set the l_link and r_link for a given kmer node, used after the tip/bubble/edge deleting process
void recalculate_kmer_links(KmerSet *kset, KmerLink *klink, uint64_t idx);


//remove dead-end tips on the kmer de bruijn graph, tips are caused by sequencing errors 
//search tip from the break-end kmer node to the branching kmer
void remove_error_tips(KmerSet *kset, KmerLink *klink,  vector <uint64_t> &tip_nodes, int tip_len_cutoff, double tip_depth_cutoff);


//the major purpose is to remove bubbles caused by large repeats, heterozygous SNPs and indels;
//only process simple bubbles, which do not have other branches on each bubble edges
//this routine can also remove some low-depth bubble edges caused by sequencing errors 
void remove_hetero_bubbles(KmerSet *kset, KmerLink *klink,  vector <uint64_t> &branch_nodes, int bubble_len_cutoff, double bubble_len_diff_rate_cutoff, double bubble_base_diff_rate_cutoff);

//Further remove small and low-depth edges, mainly caused by sequencing-errors
//these edges do not have tip or bubbble structures, they are small and low-depth between two branching kmer nodes
void remove_lowCov_edges(KmerSet *kset, KmerLink *klink,  vector <uint64_t> &branch_nodes, int edge_len_cutoff, double edge_depth_cutoff);


//get the bases and depths from a given link data-structure
void get_branch_bases(uint32_t link, vector <uint8_t> &vec_bases, vector <uint8_t> &vec_depths);


//read out all the contig sequences, which are defined as continuous non-branching kmer regions.
//start from any linear kmer nodes, extending to both sides, merge sequences from two sides to form the final contig
void read_out_contig(KmerSet *kset, KmerLink *klink);



//get all available information on a linear path, which starts from a linear or break-end kmer node, and ending at a branching kmer node
void get_linear_path(uint64_t idx, int walk_direct, int len_cutoff, int &path_len, int &path_depth, vector <uint64_t> &vec_idx, string &path_str, uint64_t &last_idx, string &break_or_branch);


//used to get the contig sequence, start from a linear kmer node, and extending leftward and rightward seperately
//the delete mechanism can deal with circular path automatically
void get_linear_seq(uint64_t idx, int walk_direct, int &seq_len, int &seq_depth,  string &seq_str, uint64_t &last_idx, string &break_or_branch, string &seq_depths, string &IsRepeat);


//compare two bubble sequences with equal length
int compare_two_seq_simple(string &str1, string &str2);


//get the leftward kmer by removing one base at right side and add one base at left side
//the inline routine code must be written in the header file, so other source file can invoke it
inline uint64_t get_next_kmer_leftward( uint64_t current_kmer , uint8_t next_base)
{	
	uint64_t assist_kmer = next_base;
	return (current_kmer >> 2) + (assist_kmer << (KmerSize - 1) * 2);
}

//get the rightward kmer by removing one base at left side and add one base at right side
//the inline routine code must be written in the header file, so other source file can invoke it
inline uint64_t get_next_kmer_rightward(uint64_t current_kmer , uint8_t next_base)
{	
	return ((current_kmer << 2) | next_base) & KmerHeadMaskVal;
}



#endif
