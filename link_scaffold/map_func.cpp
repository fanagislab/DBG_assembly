/*
This files contains functions for mapping pairs of reads onto the contig sequences.

(1) Basically, the alignment used seed-and-extension strategy, which is reallized by hash method:
   firstly build a hash, using contig kmers as keys, contig_id and contig_position as value.
   Then use a kmer as seed on read to locate the read, and then extend the alignment region and calculate identity.

(2) The alignment is not only used for scaffoling with pair-read relations, but also used for
   resolving tiny-repeat regions with linkage information inside a single read, which is equivalant to 
   the -R option in SOAPdenovo/contig.

Author: Fan Wei, fanweiagis@126.com
Date: 2015-11-30
*/

#include "map_func.h"

using namespace std;

int KmerSize = 31;                 //kmer size
double MinMapIdentity = 0.97;      //minimum mapping identity allowed to determine a true mapping
int SeedKmerNum = 5;               //number of kmers in the seed-alignment
int MinReadLen = 250;              //minimum read length, discard shorter reads
int MinCtgLen = 125;               //minimum contig length, ignore short contigs
int Input_file_format = 1;         //input file format: 1,fq; 2,fa
string Output_prefix = "./";   //output file prefix

uint64_t KmerHeadMaskVal = 0;          //used to get the kmer bit value on the forward strand from the previous one 
uint64_t KmerRCOrVal[4];               //used to get the kmer bit value on the reverse strand from the previous one 


//split a string into a vector
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


//read file_list into a vector
void reading_lib_file(string &lib_file, vector<string> &read_files)
{	
	ifstream infile ( lib_file.c_str() );
	if ( ! infile )
	{	cerr << "fail to open input file" << lib_file << endl;
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
		read_files.push_back(vec_line[0]);
	}
}


//read multiple-fasta format contig file
void read_contig_file(string &contig_seq_file, vector<string> &contig_ids, vector <string> &contig_seqs)
{
	ifstream infile ( contig_seq_file.c_str() );
	if ( ! infile )
	{       cerr << "fail to open input file " << contig_seq_file << endl;
	}

	string id, seq_str, line;
	
	while ( getline( infile, line, '\n' ) )
	{       
		if (line[0] == '>') 
		{	
			if(seq_str.size() > 0){
				
				contig_ids.push_back(id);
				contig_seqs.push_back(seq_str);
			}
			vector<string> vec_line;
			split(line, vec_line, "> \t");
			id = vec_line[0];
			seq_str.clear();
			
		}
		else if(line.size() > 0)
		{	seq_str += line;     
		}
	}

	if(seq_str.size() >= 0)
	{
		contig_ids.push_back(id);
		contig_seqs.push_back(seq_str);
	}

}

//chop contig sequences into kmers, and add the kmers into the hash
void chop_contig_to_kmerset(KmerSet *kset, vector <string> &contig_seqs)
{	
	KmerHeadMaskVal = pow_integer(2, KmerSize*2) - 1;

	KmerRCOrVal[3] = 0;
	KmerRCOrVal[1] = pow_integer(2,KmerSize*2-1);
	KmerRCOrVal[2] = pow_integer(2,KmerSize*2-1-1);
	KmerRCOrVal[0] = KmerRCOrVal[1] + KmerRCOrVal[2];
	


	int seq_num = contig_seqs.size();		
	for (int i=0; i<seq_num; i++)
	{	
		//divide sequnence into blocks seperated by Ns
		vector <string> contig_vec;
		vector <int> start_vec;
		scaffold_to_contig(contig_seqs[i], contig_vec, start_vec);
		
		for (int k=0; k<contig_vec.size(); k++)
		{
			string &contig_str = contig_vec[k];
			string kseq;
			uint64_t kbit = 0;
			uint64_t rc_kbit =  0;
			int base_bit = 4;

			for (int j=0; j<contig_str.size()-KmerSize+1; j++)
			{					
				if ( j == 0 ) 
				{	kseq = contig_str.substr(j, KmerSize);
					kbit = seq2bit(kseq);
					rc_kbit =  get_rev_com_kbit(kbit, KmerSize);
				}else
				{	base_bit = alphabet[contig_str[j+KmerSize-1]]; 
					kbit = ((kbit<<2)|base_bit) & KmerHeadMaskVal;
					rc_kbit = (rc_kbit>>2)|KmerRCOrVal[base_bit]; 
				}

				uint64_t kmer = 0;
				int direct = 0;
				if (kbit < rc_kbit)
				{	kmer = kbit;
					direct = 1;
				}else
				{	kmer = rc_kbit;
					direct = 0;
				}
				add_kmerset(kset, kmer, i, start_vec[k]+j, direct);

			}
		}
		
	}
}




//find the exact-match alignment of a specified region between read and contig, and determine the align direction
//search the seed alignment in a specified read region with given start and end positions
//use two kmers to locate the seed-alignment, use the first kmer to decide location and direction, use the second kmer to confirm
void get_align_seed(KmerSet *kset, string &read, int search_start, int search_end, int &contig_id_index, int &seed_contig_start, int &seed_contig_end, int &seed_read_start, int &seed_read_end, char &seed_direct)
{	
	KmerNode *array = kset->array;
	
	for (int i=search_start-1; i<=search_end-KmerSize-SeedKmerNum; i++)
	{	
		//get the left border kmer in a seed
		string kseq = read.substr(i,KmerSize);
		uint64_t kbit = seq2bit(kseq);
		uint64_t rc_kbit = get_rev_com_kbit(kbit, KmerSize);
		uint64_t kmer = 0;
		int direct = 0;
		if (kbit < rc_kbit)
		{	kmer = kbit;
			direct = 1;
		}else
		{	kmer = rc_kbit;
			direct = 0;
		}
		
		uint64_t idx = exist_kmerset(kset, kmer);
		if ( idx != kset->size && array[idx].freq == 1 )
		{	
			//get the right border kmer in a seed
			string kseq2 = read.substr(i+SeedKmerNum,KmerSize);
			uint64_t kbit2 = seq2bit(kseq2);
			uint64_t rc_kbit2 = get_rev_com_kbit(kbit2, KmerSize);
			uint64_t kmer2 = 0;
			uint64_t direct2 = 0;
			if (kbit2 < rc_kbit2)
			{	kmer2 = kbit2;
				direct2 = 1;
			}else
			{	kmer2 = rc_kbit2;
				direct2 = 0;
			}
			uint64_t idx2 = exist_kmerset(kset, kmer2);
			if ( idx2 != kset->size && array[idx2].freq == 1 && array[idx2].id == array[idx].id && abs(int(array[idx2].pos - array[idx].pos)) == SeedKmerNum )
			{	if (direct == array[idx].direct) //forward aligning
				{	seed_contig_start = array[idx].pos + 1;
					seed_contig_end = array[idx2].pos + KmerSize;
					seed_direct = 'F';
				}else  //reverse aligning
				{	seed_contig_start = array[idx2].pos + 1;
					seed_contig_end = array[idx].pos + KmerSize;
					seed_direct = 'R';
				}
				seed_read_start = i + 1;
				seed_read_end = i + SeedKmerNum + KmerSize;
				contig_id_index = array[idx].id;

				break;
			}
		}
	}
	
}


//based on the seed alignment, extend the alignment regions until the read ends or contig ends, and calculate the base identity value
void extend_align_region(string &read_original, string &contig, float &identity, int &align_contig_start, int &align_contig_end, int &align_read_start, int &align_read_end, char &direct)
{	
	int align_len = align_read_end - align_read_start + 1;
	int mis_match = 0;
	
	string read = read_original;

	if (direct == 'R')
	{	
		rev_com_seq(read);
		int rc_align_read_start = read.size() - align_read_start + 1;
		int rc_align_read_end = read.size() - align_read_end + 1;
		align_read_start = rc_align_read_end;
		align_read_end = rc_align_read_start;
	}
	

	//compare leftward
	while (align_read_start > 1)
	{	
		if (align_contig_start-1 < 1)
		{	break;
		}

		align_read_start --;
		align_contig_start --;
		
		align_len ++;
		if (read[align_read_start-1] != contig[align_contig_start-1])
		{	mis_match ++;
		}
	}

	//compare rightward
	while (align_read_end < read.size())
	{	
		if (align_contig_end-1 >= contig.size()-1)
		{	break;
		}

		align_read_end ++;
		align_contig_end ++;
		
		align_len ++;
		if (read[align_read_end-1] != contig[align_contig_end-1])
		{	mis_match ++;
		}
	}
	
	if (direct == 'R')
	{	
		int rc_align_read_start = read.size() - align_read_start + 1;
		int rc_align_read_end = read.size() - align_read_end + 1;
		align_read_start = rc_align_read_end;
		align_read_end = rc_align_read_start;
	}
	
	identity = 1.0 - (float)mis_match / align_len;
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


