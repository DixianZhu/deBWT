#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <time.h>
#include <inttypes.h>
#include <string.h>
#include <math.h>
#include <pthread.h>
#define MOD32 31
#define MOD64 63
#define BUFFERSIZE (1<<20)
#define BLACKLEN 16
#define BLACKCAPACITY ((uint64_t)1<<(BLACKLEN<<1))
#define REDLEN (KMER_LENGTH-BLACKLEN)
#define REDEXTRACT (((unsigned int)1<<(REDLEN<<1))-1)
#define KMER_LENGTH_PlusOne 32
#define KMER_LENGTH (KMER_LENGTH_PlusOne-1)
#define THREAD_NUM 32
#define PATH "/dupa-filer/bo/ManyBWT_Datasets/project1.24"
#define PATH_LEN 100
uint64_t convert(uint64_t mark);
int cmp(const void *a,const void *b);
uint64_t BinarySearch(uint64_t mk, uint64_t *target, int64_t up);//Search the # and $ point
uint64_t minusDimer(uint64_t cv,uint64_t index,int flag);//minus the index dimer 
char (*seeKMER(uint64_t *SA))[KMER_LENGTH+1];
void divideKmer(char (*pK)[KMER_LENGTH+1],uint64_t *SA);
int compare(const void *a,const void *b);
int compareI(uint64_t cva,uint64_t cvb,uint64_t cka,uint64_t ckb,uint64_t checka,uint64_t checkb);
char *getPath(char *name);
void *generateSpecialSA(void *arg);
#define ELIMINATE (((uint64_t)1<<((32-KMER_LENGTH)<<1))-1)
#define STKSIZ 40
uint64_t trans[256];
#define CUTOFF 8