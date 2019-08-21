/******************************************************************************
We create a hash here, refer to APE/Grape's hash.

Originally, our purpose is to store the kmer (key) and frequency (value) data.
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
"hashSet.h", and the update_hashset() function in "hashSet.cpp".

Author: 
Fan Wei (fanw@genomics.org.cn), 
Li Zhenyu (lizhenyu@genomics.org.cn);  

Date: 2010-1-26
Modified for clear and speed, on 2010-6-9
******************************************************************************/


#ifndef __KMERSET_H_
#define __KMERSET_H_

#include <iostream>
#include <cmath>
#include <string.h>
#include <inttypes.h>
using namespace std;


extern uint8_t BitOrVal[8];

//define the entity data structure (elements of the hashset-array)
typedef struct {
        uint64_t kmer;
        uint64_t id:32;
		uint64_t pos:30;
		uint64_t freq:1;    //1 unique, 0 non-unique
		uint64_t direct:1;  //1 forward; 0 reverse;

} KmerNode;


//define the hashSet data structure(include data and control parameters)
typedef struct {
	uint32_t e_size;      //the byte size for each entity 
	uint64_t size;        //the space size of the array
	uint64_t count;       //the number of current stored entities
	uint64_t count_conflict; //the number of current stored entities
	uint64_t max;         //the max number of entities allowed to store
	float load_factor;    //the ratio (allowed entity number / total array size) 
	uint64_t iter_ptr;    //the pointer to the current processing entity 
	KmerNode *array;        //the arrary to store the body data
	uint8_t *nul_flag;    //the array to store null flags, one bit for an entity
	uint8_t *del_flag;    //the array to store delete flags, one bit for an entity
} KmerSet;


//Jenkins' hash function for 64-bit integers, used in APE, convert Kmer into 64-bit integers
inline uint64_t hash_code(uint64_t kmer)
{
        kmer += ~(kmer << 32);
        kmer ^= (kmer >> 22);
        kmer += ~(kmer << 13);
        kmer ^= (kmer >> 8);
        kmer += (kmer << 3);
        kmer ^= (kmer >> 15);
        kmer += ~(kmer << 27);
        kmer ^= (kmer >> 31);
        return kmer;
}


//The hash kmer comparison function to deal with conflicts, return 1 when equals
inline int hash_equal(uint64_t kmer, KmerNode *b)
{ 
	return (kmer == b->kmer); 
}

//free the memory of KmerSet data structure
void free_hash(KmerSet *set);

//judge whether a number is a prime
int is_prime(uint64_t num);

//find the next prime number for a given number
//to decide the array size, used as the divisor to hash value, which is better to be a prime number
uint64_t find_next_prime(uint64_t num);

//Initialization of the KmerSet data structure
KmerSet *init_kmerset(uint64_t init_size, float load_factor);

//judge whether there is an entity put on position of array
//one bit represents two status: 0 for null, and 1 for fill
inline int is_entity_null(uint8_t *nul_flag, uint64_t idx) 
{
	return 1 - ( ( nul_flag[idx/8] >> (7-idx%8) ) & 0x1u );
}

//set the nul_flag to 1, mean an entity is put on position of array
inline void set_entity_fill(uint8_t *nul_flag, uint64_t idx) 
{
	nul_flag[idx/8] |= BitOrVal[idx%8];
}

//judge whether an entity has been deleted, from the del_flag.
//one bit represents two status: 0 for no, and 1 for deleted
inline int is_entity_delete(uint8_t *del_flag, uint64_t idx) 
{
	return ( del_flag[idx/8] >> (7-idx%8) ) & 0x1u;
}

//set the del_flag to 1, mean an entity was deleted on position of array
inline void set_entity_delete(uint8_t *del_flag, uint64_t idx) 
{
	del_flag[idx/8] |= BitOrVal[idx%8];
}


//enlarge the KmerSet memory space (*2 each time) when more than the allowed number of 
//entities were added into the hashset, this was automatically performed.
void enlarge_kmerset (KmerSet *set, uint64_t num);

//There are two ways to add an entity to the hashset:
//(1) When the kmer is not existed, insert it in an null(nul_flag is 0) position.
//(2) When the kmer is already existed, just update the value of entity, but do not change the kmer  
int add_kmerset(KmerSet *set, uint64_t kmer, uint32_t seq_id, uint32_t seq_pos, int kmer_direct);

//return the array index, if the entity is existed (not null, and not delete)
//return the array_size, if the entity is not existed; we can't use 0 to mean un-exist, because 0 is an normal hc value.
//Usage: if( (uint64_t idx = exist_kmerset(set,entity)) != set->size ) { cout << idx << "\t" << set->arrary[idx].kmer << endl; }
uint64_t exist_kmerset(KmerSet *set, uint64_t kmer);

//delete an existing entity from the KmerSet, by set the del_flag to 1
int delete_kmerset(KmerSet *set, uint64_t kmer);

//print the entities, include: array id, hash kmer and value
void print_kmerset_entity(KmerSet *set);

//print the hashset parameters
void print_kmerset_parameter(KmerSet *set);



#endif


