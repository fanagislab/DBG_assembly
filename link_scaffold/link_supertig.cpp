/*
This program links the contigs into scafftigs, can use one type of data as input:
a. single long reads mapping (two ends of a single reads mapped to different contigs)

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



//get the scafftig related contig id and gap information
void read_out_scaffinfo( vector <string> &contig_seqs, CtgLink **ctgLink,  LinkStat *linkStat,  int contig_num, vector < vector<int> > &ScafInfo);


//load the reads mapped to two different contigs, ready for gap filling
void load_reads_fa_file(string &mapped_twoctg_reads_file, map <string,string> &ReadsInfo);


//load the read-end mapping file, store the mapped reads information for each gap between two contigs
void load_map_twoctg_file(string &map_file, map <string, vector<vector<string> > >   &MapInfo);


//decide the gap size information [mode_size, mode_freq, total_freq, average_variance] from all the mapped reads for each gap
void decide_gap_size(map <string, vector<vector<string> > > &MapInfo, map <string,vector<int> > &GapSize);


//fill the gaps by two methods
//(1) if two neighboring contigs are overlapped, i.e. gap_size <= 0, cut out the overlapped region from the upstream contig, and linked with the downstream contig
//(2) if the neighboring contigs are not overlappe, i.e. gap_size > 0, calculate the consensus bases for the gap from all the reads mapped to this gap
//the results are output into 2 files "*scafftig*" and "*smalltig*" based on the length cutoff
void fill_gaps_inside_scaffold(vector <string> &contig_seqs, vector<string> &contig_ids, vector < vector<int> > &ScafInfo, map <string,string> &ReadsInfo, map <string, vector<vector<string> > > &MapInfo, map <string,vector<int> > &GapSize);



uint64_t total_supertig_num = 0;        //output scafftig  number
uint64_t total_supertig_len = 0;        //output scafftig   length
uint64_t total_supertig_lenwogap = 0;   //output scafftig  length without gap "N"


uint64_t contig_total_num = 0;           //total input contig number
uint64_t contig_total_len = 0;           //total input contig length
uint64_t contig_included_num = 0;        //number of contigs included into the scafftigs
uint64_t contig_included_len = 0;        //length of contigs included into the scafftigs
uint64_t contig_excluded_num = 0;        //excluded repeat contig number
uint64_t contig_excluded_len = 0;        //excluded repeat contig length

int supertig_id = -1;


void usage()
{
	cerr << "\nlink_supertig  <contig|scafftig_file.fa>  <mapping_twoctg_files.lib>\n" 
		 << "   Function: link illumina-derived scafftigs into super-contigs by pacbio reads, inside gap are filled"  << endl
		 << "   Version: 1.0"  << endl
		 << "   -n <int>   the minimum number of read-ends required to support a link, default=" << PairNumCut << endl
		 << "   -o <str>   the output prefix, set in commond-line, default = " << Output_prefix << endl
		 << "   -h         get the help information\n" << endl
		 << "Example:    link_supertig Ecoli.scafftig.seq.fa  pacbio_mapping.lib\n" << endl;
	exit(0);
}


void output_parameters()
{
	cerr  << "link_supertig   [version 1.0]"  << endl
		 << "   -n <int>   the minimum number of read-ends required to support a link, default=" << PairNumCut << endl
		 << "   -o <str>   the output prefix, set in commond-line, default = " << Output_prefix << endl
		 << "   -h         get the help information\n" << endl;
	
}


int main(int argc, char *argv[])
{	
	//get options from command line
	int c;
	while((c=getopt(argc, argv, "n:o:h")) !=-1) {
		switch(c) {
			case 'n': PairNumCut=atoi(optarg); break;
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
		//cout << i << "\t" <<contig_ids[i] << "\t" << contig_seqs[i].size() << endl;
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
		parse_read_ends_map_file(Paired_map_files[i], ctgLink, contig_seqs);
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
	string link_out_file = Output_prefix + ".supertig.links.all";
	display_data_in_link(ctgLink, linkStat, contig_num, link_out_file);
	

	remove_interleaving_links(ctgLink, contig_num, linkStat, contig_seqs);
	cerr << "\nRemoved interleave links num: " << Interleave_link_num << endl;
	
	//strong_remove_interleaving_links(ctgLink, contig_num, linkStat, contig_seqs);
	//cerr << "\nRemoved interleave links num: " << Interleave_link_num << endl;

	//remove repeat nodes
	vector<int> repeat_nodes_vec;
	remove_repeat_nodes(contig_num, linkStat, repeat_nodes_vec, contig_seqs);


	cerr << "\nRemoved repeat nodes num: " << repeat_nodes_vec.size() / 2 << endl;


	//remove all the links belong to deleted nodes
	remove_links_from_deleted_nodes(ctgLink, contig_num, linkStat);
	cerr << "\nRemoved links [related with repeat or small nodes] num: " << Deleted_link_num << endl;


	//output the links data
	link_out_file = Output_prefix + ".supertig.links.uniq";
	display_data_in_link(ctgLink, linkStat, contig_num, link_out_file);
	
	//get the scaffold infofomration, including a set of contig_ids, directions, and gap sizes
	vector < vector<int> > ScafInfo;
	read_out_scaffinfo( contig_seqs, ctgLink,  linkStat, contig_num, ScafInfo);


	//free the dynamic assigned memory in the Linked-list
	free_ctgLink_memory(ctgLink, contig_num);
	delete [] linkStat;
	


	//Start to fill the gap between neiboring contigs 
	
	//Load reads sequence which mapped to two contigs
	map <string,string> ReadsInfo;
	for (int i=0; i<Paired_map_files.size(); i++)
	{	string reads_file = Paired_map_files[i] + ".reads.fa.gz";
		cerr << "\nparse reads file: " << reads_file << endl;
		load_reads_fa_file(reads_file, ReadsInfo);
	}
	
	cerr << "load reads used to fill gaps done\n" << endl;


	//load the mapping information into memory
	map <string, vector<vector<string> > > MapInfo;
	for (int i=0; i<Paired_map_files.size(); i++)
	{	string map_file = Paired_map_files[i];
		cerr << "\nparse map file: " << map_file << endl;
		load_map_twoctg_file(map_file, MapInfo);
	}

	cerr << "load reads mapping results done\n" << endl;
	
	//decide the gap size information
	map <string,vector<int> > GapSize; 
    decide_gap_size(MapInfo, GapSize);

	cerr << "Decide the gap sizes done\n" << endl;

	
	//fill the gaps and generate scafftig sequences	
	fill_gaps_inside_scaffold(contig_seqs, contig_ids, ScafInfo, ReadsInfo, MapInfo, GapSize);
	
	cerr << "\nFill gaps inside scaffold sequence done" << endl;



	//output singlet sequences not been incorpated into scafftig sequences
	string singlet_seq_file = Output_prefix + ".supertig_repeat.seq.fa";
	ofstream SingletFile ( singlet_seq_file.c_str() );
	if ( ! SingletFile )
	{   cerr << "fail to open file" << singlet_seq_file << endl;
	}
	string singlet_pos_file = Output_prefix + ".supertig_repeat.pos.tab";
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
	{	supertig_id += 2;
		SingletFile << ">spt_" + boost::lexical_cast<string>(supertig_id) + "   fragment_num:" + boost::lexical_cast<string>(1) + "   length:" + boost::lexical_cast<string>(LenSeqVec[i].len) + "   lenwogap:" + boost::lexical_cast<string>(LenSeqVec[i].len) + "   RepeatNode\n" << LenSeqVec[i].seq << "\n";
		SingletPosFile << ">spt_" + boost::lexical_cast<string>(supertig_id) + "\n" << LenSeqVec[i].pos;
	}


	cerr << "\nTotal supertig number:          " << total_supertig_num << endl;
	cerr << "Total supertig length[WithGap]: " << total_supertig_len << endl;
	cerr << "Total supertig length[NoGap]:   " << total_supertig_lenwogap << endl;

	cerr << "\nIncluded contig number: " << contig_included_num << "  "<< (float)contig_included_num / contig_total_num  << endl;
	cerr << "Included contig length: " << contig_included_len << "  "<< (float)contig_included_len / contig_total_len  << endl;
	cerr << "Excluded repeat contig number: " << contig_excluded_num << "  "<< (float)contig_excluded_num / contig_total_num  << endl;
	cerr << "Excluded repeat contig length: " << contig_excluded_len << "  "<< (float)contig_excluded_len / contig_total_len  << endl;

	
	cerr << "\nProgram finished !" << endl;
	time_end = clock();
	cerr << "Run time: " << double(time_end - time_start) / CLOCKS_PER_SEC << endl;
	
}



//fill the gaps by pacbio consensus sequence
void fill_gaps_inside_scaffold(vector <string> &contig_seqs, vector<string> &contig_ids, vector < vector<int> > &ScafInfo, map <string,string> &ReadsInfo, map <string, vector<vector<string> > > &MapInfo, map <string,vector<int> > &GapSize)
{	

	string scaff_seq_file = Output_prefix + ".supertig.seq.fa";
	string scaff_pos_file = Output_prefix + ".supertig.pos.tab";
	string scaff_gap_file = Output_prefix + ".supertig.gap.data";

	ofstream ScafPosFile ( scaff_pos_file.c_str() );
	if ( ! ScafPosFile )
	{   cerr << "fail to open file" << scaff_pos_file << endl;
	}

	ofstream ScafGapFile ( scaff_gap_file.c_str() );
	if ( ! ScafPosFile )
	{   cerr << "fail to open file" << scaff_gap_file << endl;
	}

	ofstream ScafSeqFile (scaff_seq_file.c_str());
	if (! ScafSeqFile)
	{   cerr << "fail to open file" << scaff_seq_file << endl;
	}
	

	int gap_id = 1;
	int gap_reads_id = 1;
	
	//ScafPosFile << "#supertig_id\tblock_id\tblock_start\tblock_end\tblock_size\tdirection\tgapsize_mode_freq\tgapsize_total_freq\tgapsize_variance\tgapseq_identity" << endl;
	
	vector<LenAndSeq> LenSeqVec;
	
	//close gaps inside scaffold
	for ( int i = 0; i < ScafInfo.size(); i ++)
	{	
		string output_scafpos_res;
		string output_scafseq_res;

		string scaff_seq; 
		int scaf_len = 0;
		int scaf_lenwogap = 0;
		int scaf_ctg_num = 0;
		int block_start = 0;
		int block_size = 0;

		vector<int> &combined_vec = ScafInfo[i];
		for( int j = 0; j < combined_vec.size(); j += 2)
		{	
			int ctg_id;
			string ctg_direction;
			string ctg_seq;
			if (combined_vec[j] % 2 == 1)
			{	ctg_id = combined_vec[j];
				ctg_direction = "F";
				ctg_seq = contig_seqs[ctg_id];
			}else
			{	ctg_id = combined_vec[j] - 1;
				ctg_direction = "R";
				ctg_seq = contig_seqs[ctg_id];
				rev_com_seq(ctg_seq);
			}
			
			scaf_ctg_num ++;

			if (j + 2 >= combined_vec.size()) // the last contig in scaffold
			{	scaff_seq += ctg_seq;
				block_start = scaf_len + 1;
				block_size = ctg_seq.size();
				scaf_len += block_size;
				output_scafpos_res +=  "\t" + contig_ids[ctg_id] +  "\t" + boost::lexical_cast<string>(block_start)  +  "\t" + boost::lexical_cast<string>(scaf_len)
 +  "\t" + boost::lexical_cast<string>(block_size)  +  "\t" + ctg_direction + "\t" + ctg_seq + "\n";

				contig_included_num ++;
				contig_included_len += block_size;
				break;
			}

			//fill one gap
			int ctg2_id;
			string ctg2_direction;
			if (combined_vec[j+2] % 2 == 1)
			{	ctg2_id = combined_vec[j+2];
				ctg2_direction = "F";
			}else
			{	ctg2_id = combined_vec[j+2] - 1;
				ctg2_direction = "R";
			}
			
			string ctg_id_str = contig_ids[ctg_id];
			string ctg2_id_str = contig_ids[ctg2_id];
			
			string gap_key =  (ctg_id_str < ctg2_id_str) ? (ctg_id_str + ctg2_id_str) : (ctg2_id_str + ctg_id_str) ;
			
			int mean_gap_size = GapSize[gap_key][0];
			int min_gap_size = GapSize[gap_key][1];
			int max_gap_size = GapSize[gap_key][2];
			int total_gapsize_freq = GapSize[gap_key][3];
			int gapsize_variance = GapSize[gap_key][4];

			if (mean_gap_size <= 0) 
			{	mean_gap_size = 1;
				cerr << "Error may happens: mean_gap_size <= 0" << endl;
			}
			if ( mean_gap_size > 0)
			{	//generate the gap sequences by calling the consensus from all the mapped reads with correct[mode] insert sizes
				
				vector<vector<string> > &AllReadsVec = MapInfo[gap_key];
				

				vector<LenAndSeq> AllGapSeqs;
				string ConsensusGapSeq;
				
				for ( int k = 0; k < AllReadsVec.size(); k ++)
				{	vector<string> &OneReadVec = AllReadsVec[k];
					string read_id = OneReadVec[0];
					int read_len = atoi(OneReadVec[1].c_str());
					int read_align1_end = atoi(OneReadVec[2].c_str());
					int read_align2_start = atoi(OneReadVec[3].c_str());
					string mapped_ctg_id = OneReadVec[4];
					string mapped_ctg_direction = OneReadVec[5];
					
					int gap_size = (read_align2_start > read_align1_end) ? read_align2_start - read_align1_end - 1 : 0;
					int middle_position = (read_align1_end + read_align2_start) / 2;
					int side_extend_len = 250;
					
					string &read_seq = ReadsInfo[read_id];
					string gap_seq = read_seq.substr(middle_position - side_extend_len - gap_size/2 , gap_size + side_extend_len * 2);

					if ( (mapped_ctg_id == ctg_id_str && mapped_ctg_direction != ctg_direction) || (mapped_ctg_id == ctg2_id_str && mapped_ctg_direction != ctg2_direction)  )
					{	rev_com_seq(gap_seq);
					}
					//cerr << gap_seq << endl;
					LenAndSeq lenandseq;
					lenandseq.len = gap_seq.size();
					lenandseq.seq = gap_seq;
					AllGapSeqs.push_back(lenandseq);
				}

				sort(AllGapSeqs.begin(), AllGapSeqs.end(), cmpSeqByLen);
				int median_idx = AllGapSeqs.size() / 2;
				string &median_seq = AllGapSeqs[median_idx].seq;
				
				//>utg4 length=186621 nodes=795
				string output_gapheader;
				string output_gapcontent;
				output_gapheader = ">gap" + boost::lexical_cast<string>( gap_id ) + " length=" + boost::lexical_cast<string>( median_seq.size() ) + " nodes=";
				output_gapcontent = "Y\tS" + boost::lexical_cast<string>( gap_reads_id ) + "\t+\t0\t" + boost::lexical_cast<string>( median_seq.size() ) + "\t" + boost::lexical_cast<string>( median_seq ) + "\n";
				gap_reads_id++;
				int nodes_num = 1;
				
				for (int idx = 0; idx < AllGapSeqs.size(); idx ++)
				{	if (idx != median_idx)
					{	string &this_seq = AllGapSeqs[idx].seq;
						if (this_seq.size() > median_seq.size() * 0.75 && this_seq.size() < median_seq.size() * 1.25)
						{	
							output_gapcontent += "N\tS"  + boost::lexical_cast<string>( gap_reads_id ) + "\t+\t0\t" +  boost::lexical_cast<string>( this_seq.size() ) + "\t" + boost::lexical_cast<string>( this_seq ) + "\n";

							gap_reads_id ++;
							nodes_num ++;

						}else
						{	cerr << "Altert message:  gap_id " << gap_id << "  " << median_seq.size() << "\t" << this_seq.size() << endl;
						}
					}
				}
				output_gapheader += boost::lexical_cast<string>(nodes_num) + "\n";
				ScafGapFile << output_gapheader << output_gapcontent;

				generate_Nstr(mean_gap_size , ConsensusGapSeq );

				scaff_seq += ctg_seq + ConsensusGapSeq;
				block_start = scaf_len + 1;
				block_size = ctg_seq.size();
				scaf_len += block_size;

				output_scafpos_res +=  "\t" + contig_ids[ctg_id] +  "\t" + boost::lexical_cast<string>(block_start)  +  "\t" + boost::lexical_cast<string>(scaf_len)
 +  "\t" + boost::lexical_cast<string>(block_size)  +  "\t" + ctg_direction  + "\t" + ctg_seq + "\n";

				contig_included_num ++;
				contig_included_len += block_size;

				block_start = scaf_len + 1;
				block_size = ConsensusGapSeq.size();
				scaf_len += block_size;

				

				output_scafpos_res +=  "\tgap" +  boost::lexical_cast<string>(gap_id) + "\t" + boost::lexical_cast<string>(block_start)  +  "\t" + boost::lexical_cast<string>(scaf_len) +  "\t" + boost::lexical_cast<string>(block_size)  +  "\tN\t" + boost::lexical_cast<string>(min_gap_size)  +  "\t" + boost::lexical_cast<string>(max_gap_size)  +  "\t" + boost::lexical_cast<string>(total_gapsize_freq) + "\t" + boost::lexical_cast<string>( gapsize_variance)  + "\n"; 
				
				gap_id ++;
				//cout << "####\t"<< ConsensusGapSeq << "\t" << ConsensusSupportRate << endl;

			}

		}
		scaf_lenwogap = scaf_len;
		
		output_scafseq_res =  "   fragment_num:" + boost::lexical_cast<string>(scaf_ctg_num) + "   length:" + boost::lexical_cast<string>(scaf_len) + "   lenwogap:" + boost::lexical_cast<string>(scaf_lenwogap) + "\n" + scaff_seq + "\n";
		
		//prepare data for sorting
		LenAndSeq this_node;
		this_node.len = scaf_len;
		this_node.seq = output_scafseq_res;
		this_node.pos = output_scafpos_res;
		LenSeqVec.push_back(this_node);

		total_supertig_num ++;
		total_supertig_len += scaf_len;
		total_supertig_lenwogap += scaf_lenwogap;

	}

	//sort the scafftig results by length
	sort(LenSeqVec.begin(), LenSeqVec.end(), cmpSeqByLen);

	for(int i = 0; i < LenSeqVec.size(); i ++)
	{	supertig_id += 2;
		ScafSeqFile << ">spt_" + boost::lexical_cast<string>(supertig_id) << LenSeqVec[i].seq;
		ScafPosFile << ">spt_" + boost::lexical_cast<string>(supertig_id) + "\n" << LenSeqVec[i].pos;
	}


	ScafPosFile.close();
	ScafGapFile.close();
	ScafSeqFile.close();


}

//decide the gap size information [mode_size, mode_freq, total_freq, average_variance] from all the mapped reads for each gap
void decide_gap_size(map <string, vector<vector<string> > > &MapInfo, map <string,vector<int> > &GapSize)
{
	map <string, vector<vector<string> > >::iterator  myiter;
	for (myiter = MapInfo.begin(); myiter != MapInfo.end(); myiter++)
	{	string gapkey = myiter->first;
		
		int mean_size = 0;
		int min_size = 1000000000;
		int max_size = -1000000000;
		int total_freq = 0;
		int average_variance = 0;

		vector<vector<string> > &gapvec = myiter->second;
		for (int i = 0; i < gapvec.size(); i++ )
		{	vector<string> &linevec = gapvec[i];
			int gap_size = atoi(linevec[3].c_str()) - atoi(linevec[2].c_str()) - 1; 
			total_freq ++;
			mean_size += gap_size;
			if (gap_size > max_size)
			{	max_size = gap_size;
			}
			if (gap_size < min_size)
			{	min_size = gap_size;
			}
		}
		mean_size /= total_freq;
		
		for (int i = 0; i < gapvec.size(); i++ )
		{	vector<string> &linevec = gapvec[i];
			int gap_size = atoi(linevec[3].c_str()) - atoi(linevec[2].c_str()) - 1; 
			average_variance += abs(mean_size - gap_size);
		}
		average_variance /= total_freq;

		vector<int> gapsize_vec;
		gapsize_vec.push_back(mean_size);
		gapsize_vec.push_back(min_size);
		gapsize_vec.push_back(max_size);
		gapsize_vec.push_back(total_freq);
		gapsize_vec.push_back(average_variance);
		
		GapSize[gapkey] = gapsize_vec;

	}
}
	


//load the read-end mapping file, store the mapped reads information for each gap between two contigs
//read_400_3/1    250     1       138     ctg_763 706     569     706     F       100%    read_400_3/1    250     109     250     ctg_563 7685    1       142     F       10
void load_map_twoctg_file(string &map_file, map <string, vector<vector<string> > >   &MapInfo)
{
	igzstream infile ( map_file.c_str() );
	if ( ! infile )
	{       cerr << "fail to open input file " << map_file << endl;
	}
	
	string line;
	while ( getline( infile, line, '\n' ) )
	{	
		if ( line[0] == '#' )
		{	continue;
		}

		vector<string> vec_line;
		split(line, vec_line, " \t\n");
		
		vector<string> line_info;
		line_info.push_back(vec_line[0]);    //read_id, 0
		line_info.push_back(vec_line[1]);    //read_len, 1
		line_info.push_back(vec_line[3]);    //read_align1_end, 2
		line_info.push_back(vec_line[12]);   //read_align2_start, 3
		line_info.push_back(vec_line[4]);    //contig1, 4
		line_info.push_back(vec_line[8]);    //direct1, 5

		
		string gap_key =  (vec_line[4] < vec_line[14]) ? (vec_line[4] + vec_line[14]) : (vec_line[14] + vec_line[4]) ;
		
		MapInfo[gap_key].push_back(line_info);

	}

}

//load the reads mapped to two different contigs, ready for gap filling
void load_reads_fa_file(string &mapped_twoctg_reads_file, map <string,string> &ReadsInfo)
{
	igzstream infile ( mapped_twoctg_reads_file.c_str() );
	if ( ! infile )
	{       cerr << "fail to open input file " << mapped_twoctg_reads_file << endl;
	}
	
	string line, read_id, read_seq;
	while ( getline( infile, line, '\n' ) )
	{	
		if ( line[0] == '>' )
		{	
			vector<string> vec_line;
			split(line, vec_line, "> \t\n");
			read_id = vec_line[0];

			getline( infile, read_seq, '\n' );
			ReadsInfo[read_id] = read_seq;
		}
	}

}


//get the scafftig related contig id and gap information
void read_out_scaffinfo( vector <string> &contig_seqs, CtgLink **ctgLink,  LinkStat *linkStat,  int contig_num, vector < vector<int> > &ScafInfo)
{	
		
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

		ScafInfo.push_back(combined_ids_and_gaps);

	}


}






