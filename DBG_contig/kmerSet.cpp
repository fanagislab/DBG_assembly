/******************************************************************************
We create a hash here, refer to APE/Grape's hash.

Originally, our purpose is to store the kmer (kmer) and frequency (value) data.
However, this is a generic hash, and can be used in any purpose. Take this matter 
as example, we define a kmer and its frequency as an entity, all the entities were
in fact stored in an array, the storing position for each entity is calculated
by the hash function (hash_code) taken the kmer as input. We also use two other
arrays to store the null (nul_flag) and delete (del_flag) status for each entity.

For each entity, the array suffix (or called index, idx) are same in the 3 arrays.
When an entity is existing in the hash, the nul_flag should be 1 and the del_flag 
should be 0, else either the nul_flag is 0 or the del_flag is 1, the entity was 
thought non-existed in the hash.

The arrays can't be full-filled by entities, often the occupying rate (load_factor)
should be low than 0.7 (0.5 best), when collision happens, which means two different
kmers have the same hash value (array suffix), we just find the next null position 
for the second entity. Together with 3 arrays, the load_factor, initial array size 
(number), entity size (bytes), and some other paraters, composed a hash data structure.

There are 3 basic functions to deal with hash: hash_code, hash_equal, and free_hash,
and several other special funcitons. The array size can be dynamically increased 
(x2,two times) during the run-time, however, it is better that we pre-give a proper 
initial array size, because dynamically increasing memory need to re-distribute all
the entities in the array, which will consum much computing time.

The code here are originally designed for calculating the kmer frequency on a genome.
To adapt the code to other applications, you should modify the KmerNode (struct) in the 
"hashSet.h", and the update_kmerset() function in "hashSet.cpp".

This is the multiple-thread version, the controlling parameters of hash:
(1) nul_flag[hc]     0 for empty;  1 for filled;
(2) array[hc].key    0 for empty;  1 for filled; [ignore polyA and polyT kmers at first]

Note:
Assign 0 values to all the array[].key and array[].val values, this is essential for
the multiple-threads function using atomic-processing. This is the major difference to 
the single-thread version.

Author: 
Fan Wei (fanweiagis@126.com), 
Li Zhenyu (lizhenyu@genomics.org.cn);  

Date: 2015-10-31

******************************************************************************/

#include "seqKmer.h"
#include "kmerSet.h"

//used to set the value of nul_flag and del_flag of KmerSet
uint8_t BitOrVal[8] = {128,64,32,16,8,4,2,1};   

//used to add frequency value to the l_links and r_links of KmerNode
uint32_t BitAddVal[4] = {0x1000000u, 0x10000u, 0x100u, 0x1u}; 




//free the memory of KmerSet data structure
void free_hash(KmerSet *set)
{ 
	free(set->array);
	free(set->nul_flag);
	free(set->del_flag);
	free(set); 
}


//judge whether a number is a prime
int is_prime(uint64_t num)
{
	uint64_t i, max;
	if(num < 4) return 1;
	if(num % 2 == 0) return 0;
	max = (uint64_t)sqrt((float)num);
	for(i=3;i<max;i+=2){ 
		if(num % i == 0) return 0; 
	}
	return 1;
}

//find the next prime number for a given number
//to decide the array size, used as the divisor to hash value, which is better to be a prime number
uint64_t find_next_prime(uint64_t num)
{
	if(num % 2 == 0) num ++;
	while(1){ 
		if(is_prime(num)) {
			return num; 
		}
		num += 2; 
	}
}

//Initialization of the KmerSet data structure, parallel version
KmerSet *init_kmerset_parallel(uint64_t init_size, float load_factor, int threadNum)
{
	KmerSet *set = new KmerSet;
	
	if(init_size < 3) init_size = 3;
	else init_size = find_next_prime(init_size);

	set->e_size = sizeof(KmerNode);
	set->size   = init_size;
	set->count  = 0;
	set->count_conflict = 0;

	if(load_factor <= 0) load_factor = 0.25f;
	else if(load_factor >= 1) load_factor = 0.75f;	
	
	set->load_factor = load_factor;
	set->max    = (uint64_t) (set->size * load_factor);
	set->iter_ptr    = 0;
	
	set->array = (KmerNode *)malloc(set->size * set->e_size);  //array uninitialized
	//multiple-threads 
	memset_parallel( set->array, 0, set->size * set->e_size, threadNum );


	set->nul_flag = (uint8_t *)malloc(set->size / 8 + 1);
	memset(set->nul_flag,0,set->size / 8 + 1);               //nul_flag initialized, each bit reprensents an entity status
	set->del_flag = (uint8_t *)malloc(set->size / 8 + 1);
	memset(set->del_flag,0,set->size / 8 + 1);               //del_flag initialized, each bit reprensents an entity status
	return set;
}


//enlarge the KmerSet memory space (*2 each time) when more than the allowed number of 
//entities were added into the hashset, this was automatically performed, parallel version
void enlarge_kmerset_parallel (KmerSet *set, uint64_t num,  int threadNum) 
{	
	uint64_t old_size = set->size;
	uint64_t new_size = set->size;
	do{ new_size = find_next_prime(new_size * 2); } while(new_size * set->load_factor < set->count + num);
	
	set->size = new_size;
	
	set->array = (KmerNode*) realloc(set->array, new_size*set->e_size); //realloc can be co-used with malloc/calloc, but can't with new(C++)
	//multiple-threads
	memset_parallel( set->array + old_size, 0, (new_size - old_size) * set->e_size, threadNum );


	set->max = (uint64_t)(new_size * set->load_factor);

	uint8_t *nul_flag, *del_flag;
	nul_flag = set->nul_flag;
	del_flag = set->del_flag;
	set->nul_flag = (uint8_t*) malloc(new_size/8 + 1);
	memset(set->nul_flag, 0, new_size/8+1);
	set->del_flag = (uint8_t*) malloc(new_size/8 + 1);
	memset(set->del_flag, 0, new_size/8+1);

	KmerNode *tmp1, *tmp2;    //tempary variables for exchange purpose
	tmp1 = (KmerNode*) malloc(set->e_size);
	tmp2 = (KmerNode*) malloc(set->e_size);

	for (uint64_t i=0; i<old_size; i++)
	{
		if (is_entity_null(nul_flag, i) || is_entity_delete(del_flag, i)) continue;
		memcpy(tmp1, set->array+i, set->e_size);
		memset(set->array+i,0,set->e_size); //assign 0 value
		set_entity_delete(del_flag, i);
		while(1)
		{
			uint64_t hc = hash_code(tmp1->kmer) % set->size;
			while (!is_entity_null(set->nul_flag, hc)) { hc = (hc + 1) % set->size; }
			set_entity_fill(set->nul_flag, hc);
			if ((hc < old_size) && (!is_entity_null(nul_flag, hc)) && (!is_entity_delete(del_flag, hc)) )
			{
				memcpy(tmp2, set->array+hc, set->e_size);
				memcpy(set->array+hc, tmp1, set->e_size);
				memcpy(tmp1, tmp2, set->e_size);
				set_entity_delete(del_flag, hc);
			}
			else
			{
				memcpy(set->array+hc, tmp1, set->e_size);
				break;
			}
		}
	}

	free(nul_flag);
	free(del_flag);
	free(tmp1);
	free(tmp2);
}



/*
//This is the single-thread routine for adding a kmerNode into the KmerSet
//There are two ways to add an entity to the hashset:
//(1) When the kmer is not existed, insert it in an null(nul_flag is 0) position.
//(2) When the kmer is already existed, just update the value of entity, but do not change the kmer  
int add_kmerset(KmerSet *set, uint64_t kmer, uint8_t left_base, uint8_t right_base)
{	
	//enlarge the memory space when the entity number exceeds the preset cutoff
	//After enlarge, all the delete tags were set to 0, the null tags were set to 0 or 1 as what they should be.
	if (set->count >= set->max)
	{	enlarge_kmerset(set, 1); 
	}
	
	uint64_t hc = hash_code(kmer) % set->size;
	do{
		if(is_entity_null(set->nul_flag, hc)){ //put a new entity
			
			//****************create a new node**************************
			KmerNode *e = new KmerNode;
			e->kmer = kmer;
			e->l_link = (left_base != 4) ? BitAddVal[left_base] : 0;
			e->r_link = (right_base != 4) ? BitAddVal[right_base] : 0;
			
			//****************create a new node**************************
			
			memcpy(set->array+hc,e,set->e_size);
			set_entity_fill(set->nul_flag, hc);
			set->count ++;
			delete e;
			return 1;
		
		}else{
			if (hash_equal(kmer,set->array+hc)){ // update an existing entity
				
				KmerNode *e = set->array+hc;
				
				if (left_base != 4 &&  get_next_kmer_depth(e->l_link, left_base) < 255 )
				{	e->l_link += BitAddVal[left_base];
				}
				
				if (right_base != 4 && get_next_kmer_depth(e->r_link, right_base) < 255)
				{	e->r_link += BitAddVal[right_base];
				}

				return 1;
			}
		}
		
		set->count_conflict ++;
		if(hc + 1 == set->size) hc = 0;
		else hc = hc + 1;
	} while(1);

	return 0;
}
*/


//copy one kmernode to the kmerset at an empty position
//currently, only used to add the polyA kmer at last step
int add_node_to_kmerset(KmerSet *set, KmerNode *e)
{	
	
	uint64_t hc = hash_code(e->kmer) % set->size;
	do{
		if(is_entity_null(set->nul_flag, hc))
		{ //put a new entity	
			memcpy(set->array+hc,e,set->e_size);
			set_entity_fill(set->nul_flag, hc);
			set->count ++;
			return 1;
		
		}
				
		set->count_conflict ++;
		if(hc + 1 == set->size) hc = 0;
		else hc = hc + 1;
	} while(1);

	return 0;
}



//return the array index, if the entity is existed (not null, and not delete)
//return the array_size, if the entity is not existed; we can't use 0 to mean un-exist, because 0 is an normal hc value.
//Usage: if( (uint64_t idx = exist_kmerset(set,entity)) != set->size ) { cout << idx << "\t" << set->arrary[idx].kmer << endl; }
uint64_t exist_kmerset(KmerSet *set, uint64_t kmer)
{
	uint64_t hc = hash_code(kmer) % set->size;
	do{
		if(is_entity_null(set->nul_flag, hc)){
			return set->size;
		}else { 
			if (hash_equal(kmer,set->array+hc)){
				
				if (!is_entity_delete(set->del_flag, hc))
				{	return  hc;	
				}else {
					return set->size;
				}
			}
		}

		if(hc + 1 == set->size) hc = 0;
		else hc = hc + 1;
	} while(1);

	return set->size;
}


//delete an existing entity from the KmerSet, by set the del_flag to 1
int delete_kmerset(KmerSet *set, uint64_t kmer)
{
	uint64_t idx = exist_kmerset(set, kmer);
	if (idx != set->size){
		set_entity_delete(set->del_flag, idx);
		set->count --;
		return 1;
	}
}

//print the entities, include: array id, hash kmer and value
void print_kmerset_entity(KmerSet *set)
{
	cout << "\narray_id\thash_kmer\thash_val\n";
	KmerNode *array = set->array;
	for (uint64_t i=0; i<set->size; i++)
	{	if (!is_entity_null(set->nul_flag, i) && !is_entity_delete(set->del_flag, i))
		{	cout << i << "\t" << array[i].kmer << "\n";
		}
	}
}

//print the hashset parameters
void print_kmerset_parameter(KmerSet *set)
{
	cerr << "\nKmerset hash parameters:" << endl;
	cerr << "element_size:\t" << set->e_size << endl;
	cerr << "array_size:\t" << set->size << endl;
	cerr << "load_factor:\t" << set->load_factor << endl;
	cerr << "max_cutoff:\t" << set->max << endl;
	cerr << "iter_ptr:\t" << set->iter_ptr << endl;
	cerr << "count:\t" << set->count << endl;
	cerr << "conflict:\t" << set->count_conflict << endl;
}

uint8_t get_next_kmer_depth(uint32_t link, uint8_t base)
{
	return (link >> ((3-base)*8)) & 0xFFu;
}

//this is the thread routine to do memset parallely
void *thread_memset(void* paras)
{

	THREAD threpara = *((THREAD*)paras);
	
	memset(threpara.pointer, threpara.value, threpara.memsize); 

}

//the parallel version of memset routine, using pthread functions inside
//void pointer is equivalent to char pointer, the unit is one byte, 以字节为计算单位，void指针就是按照字节计算的
void *memset_parallel( void *pointer, int value, uint64_t memsize, int threadNum )
{
	//Note the memsize uinte is a byte
	uint64_t block_size = memsize / threadNum;
	uint64_t tail_size = memsize % threadNum;

	threadNum += 1;  //needed threads to parse all the divided blocks
	

	pthread_t *pthread = new pthread_t[threadNum];
	for (int i=0; i<threadNum; i++)
	{
		THREAD *threpara = new THREAD;
		threpara->pointer = (char*)pointer + i * block_size;
		threpara->memsize = (i != threadNum - 1) ? block_size : tail_size;
		threpara->value = 0;

		pthread_create( (pthread+i), NULL, thread_memset, (void*)threpara );
	}
	
	//当父线程读取数据结束后，等待全部子线程分析数据结束
	for (int i=0; i<threadNum; i++)
	{
		pthread_join(pthread[i], NULL);
	}

	delete [] pthread;

}


