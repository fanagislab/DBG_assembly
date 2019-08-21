//This file include variables and routines to build the kmer graph which is stored in the KmerSet
//Author: Wei Fan, fanweiagis@126.com
//Date: 2015-11-4


#include "DBGgraph.h"


//parameters for building the kmer graph
int KmerSize = 31;                     //size of kmer
int maxReadLen = 250;                  //max read length allowed, accosiate with memory allocation, reads lengther longer than this value will be trimmed off
int KmerNumInRead = 0;                 //kmer number in each reads with max allowed read length, equal to: maxReadLen - KmerSize + 1;  
int Input_file_format = 1;             //input file format: fq, 1; fa,2;
string Output_prefix = "output";       //prefix of all the output files
int threadNum = 10;                     //thread number used for parallel computation
KmerSet *kset;                         //the major data structure that store the kmer deburijn graph
double initHashSize = 1.0;             //1 G array size in the hash, as the the unit size is 16 byte, so the hash takes 16G memory
uint64_t maxDoubleHashTimes = 10;      //2^10 times larger, 1T array size, 16T memory
uint64_t doubleHashTimes = 0;          //current doubling times of the Kmerset
float hashLoadFactor = 0.7;            //loading factor of the hash
int BufferNum = 10000;                 //number of reads stored each time in the memeory buffer 
string *RawReads;                      //buffer memory to store the reads sequence
uint64_t *StoreKmer;                   //buffer memory to store the chopped kmers
uint8_t *StoreLeftBase;                //buffer memory to store the left base of the kmer
uint8_t *StoreRightBase;               //buffer memory to store the right base of the kmer
uint8_t *Signal;                       //buffer memory to store the processing status: 0, empty; 1, filled; 2, end of file;
uint64_t Kmer_total_num = 0;           //the total number of processed kmers from input files
uint64_t Total_reads_num = 0;          //the total number of processed reads from input files
uint64_t KmerHeadMaskVal = 0;          //used to get the kmer bit value on the forward strand from the previous one 
uint64_t KmerRCOrVal[4];               //used to get the kmer bit value on the reverse strand from the previous one 
KmerNode *PolyA;                       //store the information for polyA and polyT kmers

clock_t time_start;
clock_t time_end;

//this is the thread routine to convert reads to kmers and get the leftbase and rightbase parallely
//each reads is assigned to a specific thread, i.e. different threads parse different reads
void *thread_parseBlock(void* threadId_p)
{
	int threadId = *((int*)threadId_p);

	for (int i = threadId; i < BufferNum; i += threadNum)
	{		
		while (1)
		{
			if (Signal[i] == 1)
			{  
				//get the kmers from a read sequence
				string &readi = RawReads[i];

				if(readi.size() < KmerSize) 
				{	break;  //when the read is too short, skip it
				}
				
				string kseq;
				uint64_t kbit = 0;
				uint64_t rc_kbit =  0;
				uint64_t kmer = 0;
				uint8_t left_base = 4;
				uint8_t right_base = 4;
				int base_bit = 4;
				uint64_t ary_pos = 0;
				int readlen = (readi.size() > maxReadLen) ?  maxReadLen : readi.size();
				for (int j=0; j < readlen - KmerSize + 1; j++)
				{
					if ( j == 0 ) 
					{	kseq = readi.substr(0,KmerSize);
						kbit = seq2bit(kseq);
						rc_kbit =  get_rev_com_kbit(kbit, KmerSize);
					}else
					{	base_bit = alphabet[readi[j+KmerSize-1]]; 
						kbit = ((kbit<<2)|base_bit) & KmerHeadMaskVal;
						rc_kbit = (rc_kbit>>2)|KmerRCOrVal[base_bit]; 
					}

					left_base = 4;
					right_base = 4;


					if ( kbit <= rc_kbit )
					{	kmer = kbit;
						if(j > 0)  { left_base = alphabet[readi[j-1]]; }
						if(j < readlen-KmerSize) { right_base = alphabet[readi[j+KmerSize]]; }

					}else if ( kbit > rc_kbit  )
					{	kmer = rc_kbit;
						if(j > 0)  { right_base =  3 - alphabet[readi[j-1]]; }
						if(j < readlen-KmerSize) { left_base = 3 - alphabet[readi[j+KmerSize]]; }
					}
					


					ary_pos = i * KmerNumInRead + j;
					StoreKmer[ary_pos] = kmer;
					StoreLeftBase[ary_pos] = left_base;
					StoreRightBase[ary_pos] = right_base;
				
				}

				//add the kmers into KmerFreq table, use the gcc built-in sync operations
				__sync_add_and_fetch(&Kmer_total_num, readi.size()-KmerSize+1);
				
				
				break;  //finished parsing kmers from one reads, and break loop to parse the next reads
		
			}else
			{
				if (Signal[i] == 0)
				{	usleep(1); // 休息百万分之一秒，不占计算资源
				}else  
				{	break;   //when signal[i] equals 2, reach the end of file
				}
			}
		}


	}		
			
	return NULL;
}



//this is the thread routine to convert reads to kmers parallely
//each specific kmer species is assigned to a specific thread, i.e. different threads parse different kmer species
void *thread_updatekmers(void* threadId_p)
{
	int threadId = *((int*)threadId_p);
	
	for (int i = 0; i < BufferNum; i ++)
	{
		
		if ( Signal[i] == 2 )
		{	break;
		}

		int readlen = (RawReads[i].size() > maxReadLen) ?  maxReadLen : RawReads[i].size();

		for(int j=0; j<readlen-KmerSize+1; j++)
		{	
			uint64_t ary_pos = i * KmerNumInRead + j;
			uint64_t kmer = StoreKmer[ary_pos];
			uint8_t left_base = StoreLeftBase[ary_pos];
			uint8_t right_base = StoreRightBase[ary_pos];
			
			
			//choose kmers belong to this thread
			if ( kmer % threadNum != threadId)
			{	continue;
			}

			//polyA and polyT are not added into kmerset at first, kmer[kbit==0] was processed by only the first child thread[threadId==0], so no need lock
			if (kmer == 0) 
			{
				if (left_base != 4 &&  get_next_kmer_depth(PolyA->l_link, left_base) < 255 )
				{	PolyA->l_link += BitAddVal[left_base];
				}
				
				if (right_base != 4 && get_next_kmer_depth(PolyA->r_link, right_base) < 255 )
				{	PolyA->r_link += BitAddVal[right_base];
				}

				continue;
			} 

					
			uint64_t hc = hash_code(kmer) % kset->size;

			//////////////////add a kmer-freq into the KmerSet///////////////////////
			do
			{   
				KmerNode *e = kset->array + hc;

				if( e->kmer == 0 && __sync_bool_compare_and_swap(&e->kmer, 0, kmer) )
				{ //put a new entity, no thread competition, so not need to lock
					
					e->l_link = (left_base != 4) ? BitAddVal[left_base] : 0;
					e->r_link = (right_base != 4) ? BitAddVal[right_base] : 0;
					
					__sync_or_and_fetch(kset->nul_flag+hc/8, BitOrVal[hc%8]); // need lock, set the nul_flag value to 1
					__sync_add_and_fetch(&kset->count ,1); // need lock, update kset->count
					break;
				}else
				{
					if ( e->kmer == kmer )
					{ // update an existing entity, no thread competition, so no need to lock
				
						if (left_base != 4 &&  get_next_kmer_depth(e->l_link, left_base) < 255 )
						{	e->l_link += BitAddVal[left_base];
						}
				
						if (right_base != 4 && get_next_kmer_depth(e->r_link, right_base) < 255 )
						{	e->r_link += BitAddVal[right_base];
						}
	
						break;
		
					}
				}
						
				__sync_add_and_fetch(&kset->count_conflict ,1);

				if(hc + 1 == kset->size) hc = 0;
				else hc = hc + 1;
			} while(1);
			//////////////////add a kmer-freq into the KmerSet///////////////////////

		}
	}


	return NULL;
}


//parse one reads file in parallel, controlling the children threads
void parse_one_reads_file(string &reads_file)
{
	string LineStr;
	igzstream currentFile; //Caution: gzstream.cpp has a bug, which can be overcomed by gcc/g++ optimation (-O1, -O2, -O3).
	currentFile.open(reads_file.c_str());
	

	uint64_t  ReadsNum = 0;

	while (1)
	{
		ReadsNum = 0;
		memset(Signal, 0, BufferNum);

		
		//create and run children threads to parse the reads that have loaded into the memory buffer
		pthread_t *pthread = new pthread_t[threadNum];
		int *pthreadId = new int[threadNum];
		for (int i=0; i<threadNum; i++)
		{
			pthreadId[i] = i;
			pthread_create((pthread+i), NULL, thread_parseBlock, (void*)(pthreadId+i));

		}
		cerr << "\n" << threadNum << " children threads created!" << endl;
		
		//load the reads data into memory buffer by the father thread
		if (Input_file_format == 1)
		{	//读取fq格式文件, support one-line fastq format
			while ( ReadsNum < BufferNum && getline( currentFile, LineStr, '\n' ) )
			{	
				if (LineStr[0] == '@') 
				{	
					getline( currentFile, RawReads[ReadsNum], '\n');


					getline( currentFile, LineStr, '\n');
					getline( currentFile, LineStr, '\n');
					Signal[ReadsNum] = 1;
					ReadsNum ++;
				}	
			}
		}
		else
		{	//读取fa格式文件, support one-line fasta format
			while ( ReadsNum < BufferNum && getline( currentFile, LineStr, '\n' )  )
			{	
				if (LineStr[0] == '>') 
				{	getline( currentFile, RawReads[ReadsNum], '\n');	


					Signal[ReadsNum] = 1;	
					ReadsNum ++;
				}	
			}
		}
		
		Total_reads_num += ReadsNum;

		cerr << "Load reads block " << Total_reads_num << endl;

		
		//judge the end of file, and make signal for the children threads
		if ( ReadsNum < BufferNum )
		{
			for (int i = ReadsNum; i < BufferNum; i ++) 
			{	
				Signal[i] = 2;
			}
			
			cerr << "this block has reach the end of file " << endl;
		}

		

		//当父线程读取数据结束后，等待全部子线程分析数据结束
		for (int i=0; i<threadNum; i++)
		{
			pthread_join(pthread[i], NULL);
		}
		
		delete [] pthread;
		delete [] pthreadId;
		
		cerr << "chop reads to kmers done" << endl;
		
		
		//add kmers to the kmerset by multiple threads
		/////////////////////////////////////////////////////////////////////////////////////////
		pthread = new pthread_t[threadNum];
		pthreadId = new int[threadNum];
		for (int i=0; i<threadNum; i++)
		{
			pthreadId[i] = i;
			pthread_create((pthread+i), NULL, thread_updatekmers, (void*)(pthreadId+i));

		}

		for (int i=0; i<threadNum; i++)
		{
			pthread_join(pthread[i], NULL);
		}

		delete [] pthread;
		delete [] pthreadId;

		cerr << "add kmers to hash done" << endl;

		/////////////////////////////////////////////////////////////////////////////////////////


		//when reach the end of file, and all children threads finished
		if ( ReadsNum < BufferNum )
		{	break;
		}
		
		//cerr << kset->count << "\t" << kset->max << "\t" << kset->size << "\t" << kset->count_conflict << "\t" << kset->load_factor << endl;
		//cerr << "Parsed reads in current block\n" << endl;

		//Judge whether need to enlarge the hash size, each enlarge double the original hash size
		if( kset->count > kset->max  ) 
		{
			if ( doubleHashTimes < maxDoubleHashTimes ) 
			{
				enlarge_kmerset_parallel(kset, 1, threadNum); //enlarge the memory space when the entity number exceeds the preset cutoff
				doubleHashTimes ++;
				cerr << "Enlarge hash array size to be: " << kset->size << endl;
				cerr << "The expanded memory used now:  "	<<  (double)kset->size/1000000000 * 16 << " G" << endl;
				
			}else
			{	
				cerr << "\nAlert message: Memory reach the maximum allowed, program have loaded " << Total_reads_num << " reads, the left others are ignored\n" << endl;
				break; 
			}
		}
			
	}


	currentFile.close();


}



//initialization, load each input files and construct the kmer de bruijn graph
void build_debruijn_graph( vector<string> &reads_files)
{	
	
	time_start = clock();

	

	KmerHeadMaskVal = pow_integer(2, KmerSize*2) - 1;

	KmerRCOrVal[3] = 0;
	KmerRCOrVal[1] = pow_integer(2,KmerSize*2-1);
	KmerRCOrVal[2] = pow_integer(2,KmerSize*2-1-1);;
	KmerRCOrVal[0] = KmerRCOrVal[1] + KmerRCOrVal[2];

	
	//cerr << initHashSize*1000000000 << "\t" << hashLoadFactor << "\t" << threadNum << endl;
	cerr << "Start to initialize the kmerset hash" << endl;
	kset = init_kmerset_parallel(initHashSize*1000000000, hashLoadFactor, threadNum); 
	cerr << "Hash initialization array size:  " << initHashSize <<  " G" << endl; 
	cerr << "The initialization memory used:  "	<<  initHashSize * 16 << " G" << endl;
	time_end = clock();
	cerr << "Finished! Run time: " << double(time_end - time_start) / CLOCKS_PER_SEC << endl;

	
	//assign buffer memory
	KmerNumInRead = maxReadLen - KmerSize + 1;
	

	RawReads = new string[BufferNum];  ////from here can get read length information
	Signal = new uint8_t[BufferNum];   ////from here can get whether this read reached the end of file
	StoreKmer = new uint64_t[BufferNum * KmerNumInRead];
	StoreLeftBase = new uint8_t[BufferNum * KmerNumInRead];
	StoreRightBase = new uint8_t[BufferNum * KmerNumInRead];
	
	
	PolyA = new KmerNode;
	PolyA->kmer = 0;
	PolyA->l_link = 0;
	PolyA->r_link = 0;

	//read data into memory, construct kmer de bruijn graph
	cerr << "\nparse input reads files: " << endl; 
	for (int i = 0; i < reads_files.size(); i++)
	{
		cerr << "\nStart to parse reads file: " << reads_files[i] << endl;
		parse_one_reads_file(reads_files[i]);
		cerr << "\nTotal number of reads loaded into memory: " << Total_reads_num  << endl;
		cerr << "Total number of kmers loaded into memory: " << Kmer_total_num << endl;
		time_end = clock();
		cerr << "Finished! Run time: " << double(time_end - time_start) / CLOCKS_PER_SEC << endl;
	}


	//add polyA and polyT [kmer: 0] to the kmerset
	add_node_to_kmerset(kset, PolyA);

	
	delete [] RawReads;
	delete [] Signal;
	delete [] StoreKmer;
	delete [] StoreLeftBase;
	delete [] StoreRightBase;
	delete PolyA;
	
	print_kmerset_parameter(kset);

}
