//collect all kmers whose precursor or successor is # and store out into a file
//collect last k suffixes and other suffixes whose initial k-mers contain # 
#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <time.h>
#include <inttypes.h>
#include <string.h>
#include"kseq.h" 
#include"collect#$.h"
#define MOD32 31
#define MOD64 63
#define KMER_LENGTH 30
KSEQ_INIT(gzFile, gzread);
int collect(char *path)
{
   gzFile fp = gzopen(path,"r");
   kseq_t *seq=kseq_init(fp);
   uint64_t countRead=0;
   uint64_t i,count=0,j;
   while(kseq_read(seq)>=0)
   {
      printf("name: %s\n", seq->name.s);
      printf("counting the seq length:%ld\n",seq->seq.l);
      countRead++;
   }
   kseq_destroy(seq);
   gzclose(fp);
   fp = gzopen(path,"r");
   if(fp==NULL) {printf("can not open ref file\n");exit(0);}
   seq=kseq_init(fp);
   uint64_t *specialHead=(uint64_t*)calloc(countRead-1,sizeof(uint64_t));
   uint64_t *specialTail=(uint64_t*)calloc(countRead,sizeof(uint64_t));
   while(kseq_read(seq)>=0)
   {
   	printf("name: %s\n", seq->name.s);
   	printf("counting the seq length:%ld\n",seq->seq.l);
   	uint64_t temp=0;
   	if(count>0)
   	{
   		for(i=0,j=31;i<KMER_LENGTH;i++,j--)
    	{
    		//store the head characters of each read
    		temp=temp|trans[seq->seq.s[i]]<<(j<<1);
    	}
    	specialHead[count]=temp;
    	printf("%lu Head is:%#lx\n",count,specialHead[count]);
    	temp=0;
    }
    for(i=seq->seq.l-KMER_LENGTH,j=31;i<seq->seq.l;i++,j--)
    {
    	//store the tail characters of the first read
    	temp=temp|trans[seq->seq.s[i]]<<(j<<1);
    }
    specialTail[count]=temp;
    char *tail=&(seq->seq.s[seq->seq.l-KMER_LENGTH]);
    printf("%lu Tail is:%#lx\nthe string is %s\n",count,specialTail[count],tail);
    count++;
   }
   //kseq_destroy(seq); segment fault!
   gzclose(fp);
   
   return 1;
}