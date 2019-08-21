/* This program is used to build contigs by kmer de bruijn graph method.


(1) All the processes are performed on the kmer de bruijn graphs, this is 
similar to APE; but different from velvet and grape, which firstly simplified
the kmer graph into edge graph, and then mainly processed on the edge graph. 

(2) One remained question in short-read assembly method, is how to use all the available
information. Previous program such as velvet and grape mainly focused on the 
topology structure, such as tip and bubble, but do not fully use the coverage depth
information. This program considers both graph topology and coverage depth (stored
from 0-255), which will benifit in distinguishing repeatitive and heterozygous sequences. 

(3) You are recommended to use the "-D" option, to remove the extrem low-frequency kmer-links
on the graph, which are mostly caused by sequencing errors. This action can also break down
many bubble structures into tips, which make it easy for error cleaning in later steps.

(4) By a combination of removing tips, bubbles, and low-coverage linear edges, most of the 
sequencing errors will be removed, and then the linear parts of the graph grows significantly. 
Note that all these graph-simplifying steps have coverage depth and length cutoff,
in order to avoid removing the genine sequences. 

(5) The bubble processing functions can be used to remove sequencing errors, they can also 
be used to merge heterozygous SNPs and Indels. The repeatitive and heterozgyous sequences will
both form bubble structures, however, repeatitive sequences will have higher coverage depth
than heterozygous sequences. 

(6) The final output is contig sequnces(length larger than a cutoff), which is defined as combination of linear kmer nodes
on the graph, i.e, continous unique Kmer regions in the genome. Contig will break down on 
the branching kmer nodes, which are formed by repeat structures. Those small contigs with
length smaller the the cutoff are output into another result file (*.small.fa).

(7) Along with the contig sequences, the coverage depth for each base is also output, the coverage
depth can be viewed as the qualily evaluation of a base, given that high coverage depth reflects
high accuracy. The leftward and rightward ending kmers were also output as the contig sequences,
all the contig squences that starting or ending on a branching kmer node, will share the ending
kmers. You can use this information to know which contigs are neighbors.

pass polyA and polyT kmers, i.e. polyA and polyT are not used in build DBG

Author: Fan Wei, fanweiagis@126.com
Version: 1.0
Date: 2015-11-9

*/


#include "seqKmer.h"
#include "kmerSet.h"
#include "DBGgraph.h"
#include "contig.h"


using namespace std;



void usage()
{
	
	cout << "Function introduction:\
\n\ndebrujin_contig constructs a kmer de bruijn graph storing in a kmer hash,\ 
Each kmer unit in the kmer hash occupies 146 bits [18.25 bytes] memory, the total\
memory is equal to hash_size * 18.25 bytes. The hash size can be set by -i parameter,\
if the kmer species number in the reads data exceeds hash_size * loading_factor,\
then the hash will expand its size automatically, doubling each time, the maximum\
expansion times can be set by -e parameter. This is helpful in control the final memory \
usage in order not to happen out of memory error. When it meets the maximum expansion times,\
the program will stop parsing any more reads, and the current kmer debruijn graph\
will be used for assembling. This means that some parts of reads data are not used.\
\n\nAll the processes are performed on the kmer de bruijn graphs. This is\
similar to APE; but different from Velvet and SOAPdenovo, which firstly simplified\
the kmer graph into edge graph, and then mainly processed on the edge graph.\
One remained question in short-read assembly method, is how to use all the available\
information. Previous program such as velvet and grape mainly focused on the\
topology structure, such as tip and bubble, but do not fully use the coverage depth\
information. This program considers both graph topology and coverage depth (stored\
from 0-255), which will be helpful in distinguishing repetitive and heterozygous sequences.\
\n\nYou are recommended to use the -d option, to remove the extrem low frequency kmers\
on the graph, which are mostly caused by sequencing errors. This action can also break down\
many bubble structures into tips, which make it easy for error cleaning in later steps.\
By a combination of removing tips, bubbles, and low-coverage linear edges, most of the\
sequencing errors will be removed, and then the linear parts of the graph grows significantly.\
Whether to remove tips, whether to remove bubbles, and whether to remove low coverage edges\
can be set by 3 paramters, along for the maximum length and coverage depth of these items.\
Note that all these graph-simplifying steps have coverage depth and length cutoff,\
in order to avoid removing the genine sequences.\
\n\nThe final output is contig sequnces, which is defined as combination of linear kmer nodes\
on the graph, i.e, continous unique Kmer regions in the genome.\
Along with the contig sequences, the coverage depth for each base is also output\
 The leftward and rightward ending kmers were also output as the contig sequences,\
all the contig squences that starting or ending on a branching kmer node, will share the ending\
kmers. Besides, the sequence and coverage depth information of tips, bubbles, and low coverage edges\
are also output to files. Almost all the information on the debruijn graph have been saved to files\n";

	
	cout << "\ndebruijn_contig   <reads_file.lib>\n" 
		 << "   \nFunction: This program is used to build contigs by kmer de bruijn graph method\n" << endl
		 << "   Verion: 1.0\n" << endl
	 	 << "   -k <int>   set kmer size, max 31, default=" << KmerSize << endl
		 << "   -r <int>   set maximum allowed read length, trimmed if longer, default=" << maxReadLen << endl
		 << "   -f <int>   set the input file format: 1: fq|gz(one-line), 2: fa|gz(one-line), default=" << Input_file_format << endl
		 << "   -o <str>   set the output file prefix, default = " << Output_prefix << endl
         << "   -t <int>   thread number to run in parallel, default=" << threadNum << endl
         << "   -i <float>  set initialization size (uint:G) of the kmer-hash, memory consumption ( * 16 G bytes ), default=" << initHashSize << endl
	     << "   -l <float>  set loading factor of the kmer hash, default=" << hashLoadFactor  << endl
	     << "   -e <int>  max doubling times of hash size allowed to enlarge memory consumption, default=" << maxDoubleHashTimes  << endl
	     << "   -b <int>  buffer size: number of reads loading into the buffer memory, default=" << BufferNum << endl

		 << "   -D <int>   delete kmer-links with frequency no larger than, default=" << KmerFreqCutoff << endl
		 << "   -T <int>   whether cut off tip-branch, 1:yes; 0:no; default=" << is_remove_tip << endl
		 << "   -I <int>   set the max allowed tip-branch length, default=" << Tip_len_cutoff << endl
		 << "   -P <float>  set the max allowed tip-branch depth, default=" << Tip_depth_cutoff << endl
		 << "   -W <int>   wheter cut off low-coverage branch between two branching nodes, 1:yes; 0:no; default=" << is_remove_lowedge << endl
		 << "   -C <int>    set the max allowed length for low-coverage branch, default=" << LowCovEdge_len_cutoff << endl
		 << "   -G <float>  set the max allowed depth for low-coverage branch, default=" << LowCovEdge_depth_cutoff << endl
		 << "   -B <int>   whether cut off the low-coverage branch for pairs of bubble branches, 1:yes; 0:no; default=" << is_remove_bubble << endl
		 << "   -U <int>   set the max allowed bubble-branch length, default=" << Bubble_len_cutoff << endl
		 << "   -L <float>   set the max allowed length difference rate between the two bubble-branchess, default=" << Bubble_len_diff_rate_cutoff << endl
		 << "   -E <float>  set the max allowed base difference rate between the two bubble-branches, default=" << Bubble_base_diff_rate_cutoff << endl
		 << "   -M <int>    set the minimum length for contig to output, default=" << Contig_len_cutoff << endl
		 << "   -h         get the help information" << endl
		 << endl
		 <<  "\nExample: \ndebruijn_contig  -k 31 -r 250  -t 10  -i 0.1  -M 125 -o Ecoli reads_files.lib   2> reads_files.debruijn_contig.log \n" << endl;
	exit(0);
}



void output_program_parameters()
{	
	cerr << "\nProgram parameters setting:" << endl;
	
	cerr << "   -k <int>   set kmer size, default=" << KmerSize << endl
		 << "   -r <int>   set maximum allowed read length, default=" << maxReadLen << endl
		 << "   -f <int>   set the input file format: 1: fq|gz(one-line), 2: fa|gz(one-line), default=" << Input_file_format << endl
		 << "   -o <str>   set the output prefix, default = " << Output_prefix  << endl
         << "   -t <int>   run the program in multiple thread mode, default=" << threadNum << endl
         << "   -i <float>  set initialization size (uint:G) of kmer-hash, memory consumption(* 16 G bytes), default=" << initHashSize <<"G" << endl
	     << "   -l <float>  set loading factor of the hash, default=" << hashLoadFactor  << endl
	     << "   -e <int>  max doubling times of hash size allowed to enlarge memory consumption, default=" << maxDoubleHashTimes  << endl
	     << "   -b <int>  buff size: number of reads loading into the buffer memory, default=" << BufferNum << endl

		 << "   -D <int>   delete kmer-links with frequency no larger than, default=" << KmerFreqCutoff << endl
		 << "   -T <int>   wheter cut off tip branches, 1:yes; 0:no; default=" << is_remove_tip << endl
		 << "   -I <int>   set the max allowed tip length, default=" << Tip_len_cutoff << endl
		 << "   -P <float>  set the max allowed tip depth, default=" << Tip_depth_cutoff << endl
		 << "   -W <int>   wheter cut off low coverage edges between two branching nodes, 1:yes; 0:no; default=" << is_remove_lowedge << endl
		 << "   -C <int>    set the max allowed length for low coverage edges, default=" << LowCovEdge_len_cutoff << endl
		 << "   -G <float>  set the max allowed depth for low coverage edges, default=" << LowCovEdge_depth_cutoff << endl
		 << "   -B <int>   wheter cut off bubble branches, 1:yes; 0:no; default=" << is_remove_bubble << endl
		 << "   -U <int>   set the max allowed bubble length, default=" << Bubble_len_cutoff << endl
		 << "   -L <float>   set the max allowed length difference rate between the two bubbles, default=" << Bubble_len_diff_rate_cutoff << endl
		 << "   -E <float>  set the max allowed base difference rate between the two bubbles, default=" << Bubble_base_diff_rate_cutoff << endl
		 << "   -M <int>    set the minimum length for contig to output, default=" << Contig_len_cutoff << endl
		 << endl;

}



int main(int argc, char *argv[])
{	
	//get options from command line
	int c;
	while((c=getopt(argc, argv, "k:r:f:o:t:i:l:e:b:D:T:I:P:W:C:G:B:U:L:E:M:h")) !=-1) {
		switch(c) {
			case 'k': KmerSize=atoi(optarg); break;
			case 'r': maxReadLen=atoi(optarg); break;
			case 'f': Input_file_format=atoi(optarg); break;
			case 'o': Output_prefix=optarg; break;
			case 't': threadNum=atoi(optarg); break;
			case 'i': initHashSize=atof(optarg); break;
			case 'l': hashLoadFactor=atof(optarg); break;
			case 'e': maxDoubleHashTimes=atoi(optarg); break;
			case 'b': BufferNum=atoi(optarg); break;

			case 'D': KmerFreqCutoff=atoi(optarg); break;
			case 'T': is_remove_tip=atoi(optarg); break;
			case 'I': Tip_len_cutoff=atoi(optarg); break;
			case 'P': Tip_depth_cutoff=atof(optarg); break;
			case 'W': is_remove_lowedge=atoi(optarg); break;
			case 'C': LowCovEdge_len_cutoff=atoi(optarg); break;
			case 'G': LowCovEdge_depth_cutoff=atof(optarg); break;
			case 'B': is_remove_bubble=atoi(optarg); break;
			case 'U': Bubble_len_cutoff=atoi(optarg); break;
			case 'L': Bubble_len_diff_rate_cutoff=atof(optarg); break;
			case 'E': Bubble_base_diff_rate_cutoff=atof(optarg); break;
			case 'M': Contig_len_cutoff=atof(optarg); break;
			case 'h': usage(); break;
			default: usage();
		}
	}
	
	if (argc < 2) usage();
	
	output_program_parameters();
	
	string reads_lib_file = argv[optind++];

	vector <string> reads_files;
	reading_file_list(reads_lib_file, reads_files);

	build_debruijn_graph(reads_files);
	cerr << "\nLoad reads, chop kmer, build kmer graph finished !" << endl;
	
	build_contig_sequence();
	cerr << "\nRemove tips, merge bubbles, output contig sequence finished !" << endl;

	cerr << "\nAssembly completely finished!" << endl;

}


