#ifndef __LINKFUNC_H_
#define __LINKFUNC_H_

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <set>
#include <cmath>
#include "seqKmer.h"
#include <inttypes.h>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
using namespace std;

extern string Output_prefix;     //output file prefix
extern int InsertSize;           //insert size of pair-end and mate-pair data
extern int PairNumCut;           //cutoff for number of read-pairs to support a contig link

extern int FR_link_num ;         // read1 forward mapping, read2 reverse mapping
extern int RF_link_num ;         // read1 reverse mappin,g read2 forward mapping
extern int FF_link_num ;         // both read1 and read2 forward mapping
extern int RR_link_num ;         // both read1 and read2 reverse mapping
extern int Wrong_link_num ;      // wrong links which do not belong to the above four mapping types
extern int LowFreq_link_num ;    // deleted links number which have supporting read pair number less than PairNumCut
extern int Deleted_link_num ;    // deleted links number derived from deleted repeat and small contig nodes
extern int Interleave_link_num ; // deleted interleaving links number

//define the CtgLink data structure
typedef struct CtgLink {
	uint64_t id:54;        //the 3¡¯-linked contig id, range unlimited
	uint64_t freq:10;      //the total number of read-ends or read-pairs supported this link, range 1023
	int64_t size;          //the total accumulated sum value of all estimated gap sizes, range unlimited
	CtgLink *pointer; //the pointer to next linked node
} CtgLink;

typedef struct {
	uint8_t link;      // number of 3'-links for a contig node, range 255
	uint8_t inlink;    // number of 5'-links for a contig node, range 255
	uint8_t del:1;     // delete flag: 0, not deleted; 1, been deleted
} LinkStat;


//used for sorting the gap-seqs by its length
typedef struct{
        uint64_t len;
        string seq;
		string pos;
}LenAndSeq;



bool cmpSeqByLen (LenAndSeq a, LenAndSeq b);


//read file_list into a vector
void reading_para_file(string &para_file,  vector<string> &Paired_map_files);

//read multiple-fasta format contig file
void read_contig_file(string &contig_seq_file, vector <string> &contig_seqs, vector <string> &contig_ids);


//parse the read-ends mapping file, add the links information into the memory
void parse_read_ends_map_file(string &PE_map_file, CtgLink **ctgLink, vector <string> &contig_seqs);

//parse the pair-end mapping file, add the links information into the memory
void parse_pair_ends_map_file(string &PE_map_file, CtgLink **ctgLink, vector <string> &contig_seqs);

//parse the mate-pair mapping file, add the links information into the memory
void parse_mate_pairs_map_file(string &PE_map_file, CtgLink **ctgLink, vector <string> &contig_seqs);

//add one link information into the link data structure
void add_data_into_link(CtgLink **ctgLink, int source_id, int target_id, int insert_size);

//output data in the link data-structure
void display_data_in_link(CtgLink **ctgLink, LinkStat *linkStat, int contig_num, string &link_out_file);


//remove the links with low frequency < PairNumCut
void remove_lowfreq_link_and_stat(CtgLink **ctgLink, int contig_num, LinkStat *linkStat);

//mainly used for solving pair-end or mate-pair data, because single-read mapping can avoid interleaving problem
//solve the interleaving problem, by removing the interleaving links
//the middle contig in a interleaving case, has only one inlink and only one outlink
void remove_interleaving_links(CtgLink **ctgLink, int contig_num, LinkStat *linkStat, vector<string> &contig_seqs);
void strong_remove_interleaving_links(CtgLink **ctgLink, int contig_num, LinkStat *linkStat, vector<string> &contig_seqs);

//get all the linked contig ids for a given contig node
void get_all_linked_ids( CtgLink **ctgLink, int source_id, vector<int> &linked_ids, vector<int> &insert_sizes );

//delete a specified link between source node and target node
void delete_linked_id( CtgLink **ctgLink, LinkStat *linkStat, int source_id, int target_id );


//remove repeat nodes which have multiple incoming and outgoing links
void remove_repeat_nodes(int contig_num, LinkStat *linkStat, vector<int> &repeat_nodes, vector<string> &contig_seqs);

//remove small nodes whose length are less than half of the Insert_size
void remove_small_nodes(int contig_num, LinkStat *linkStat, vector<int> &small_nodes, vector<string> &contig_seqs);


//remove all the links that related with repeat nodes
void remove_links_from_deleted_nodes(CtgLink **ctgLink, int contig_num, LinkStat *linkStat);


//get the set of linear nodes from a given start node
void get_linear_seq( CtgLink **ctgLink,  LinkStat *linkStat, int start_id, vector<int> &vec_id );

//get the next linked contig id for current contig node [.link == 1]
int get_next_linked_id( CtgLink **ctgLink, int source_id, int &insert_size);

//make a string with specified number of Ns
void generate_Nstr(int num, string &Nstr);


//split scaffold into contigs and stored in vector
void scaffold_to_contig(string &scaffold, vector <string> &contigs, vector <int> &starts);

//release the memory of ctgLink
void free_ctgLink_memory(CtgLink **ctgLink, int contig_num);


//split a string into a vector<string>
void split(string &strLine, vector<string>& tokens, const char* delim);


//convert "ctg_NNN" to "NNN"
inline int ctgStr2Id(string ctg_id_str)
{	return atoi( (ctg_id_str.substr(4,ctg_id_str.size()-4)).c_str() );
}

//get the paired contig id
inline int get_pair_id( int id)
{   return ( id % 2 ==0 ) ? (id - 1) : (id + 1);
}


#endif
