#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <time.h>
#include <inttypes.h>
#include <string.h>
#include "kseq.h"
#include "main.h"
#include "collect#$.h"
KSEQ_INIT(gzFile, gzread);
uint64_t trans[256];
int main()
{
	trans['A']=trans['a']=0;
    trans['C']=trans['c']=1;
    trans['N']=trans['G']=trans['g']=2;
    trans['T']=trans['t']=3;
    trans['#']=4;// patch for special characters
    trans['$']=5;
    /* Test the kmercounting is whether right or wrong
    FILE *fp = fopen("/media/dell/681CF0681CF03324/TempData/out","r");
    if(fp==NULL) {printf("can not open out!\n");exit(0);}
    uint64_t sum=0;
    char temp[32];
    while(!feof(fp))
    {
    	unsigned int count;
    	fscanf(fp,"%s%u\n",temp,&count);
    	sum+=count;
    	//printf("count = %u \n",count );
    }
    printf("The total counts are %lu\n",sum);
    */
    clock_t collect_t=clock();
    char *refPath=getPath("/bigGenome10FB.fa");
    if(collect(refPath)==1)
    {
    	printf("success test!\n");
    }
    else
    {
    	printf("no reads\n");
    }
    free(refPath);
    clock_t generateBlock_t=clock();
    double Collect=(double)(generateBlock_t-collect_t)/CLOCKS_PER_SEC;
    printf("collect#$ time is %lf\n",Collect);
    if(generateBlocks()==1)
    {
    	printf("success generate blocks!\n");
    }
    else
    {
    	printf("failed to generate blocks!\n");
    }
    clock_t generateSP_t=clock();
    time_t genSP_start=time(0);
    double GenerateBlock=(double)(generateSP_t-generateBlock_t)/CLOCKS_PER_SEC;
    printf("generate blocks time is %lf\n",GenerateBlock);
    if(generateSP()==1)
    {
        printf("success generate SP code\n");
    }
    else
    {
        printf("failed to generate SP code\n");
    }
    time_t genSP_end=time(0);
    long genSP=(genSP_end-genSP_start);
    clock_t sortBlue_t=clock();
    time_t sortBlue_start=time(0);
    double GenerateSP=(double)(sortBlue_t-generateSP_t)/CLOCKS_PER_SEC;
    printf("generate SP time is %lf, (%lu)\n",GenerateSP,genSP);
    if(sortBlue()==1)
    {
        printf("success sort blueTable\n");
    }
    else
    {
        printf("failed to sort blueTable\n");   
    }
    clock_t insertCase3_t=clock();
    time_t sortBlue_end=time(0);
    double SortBlue=(double)(insertCase3_t-sortBlue_t)/CLOCKS_PER_SEC;
    long soBlue=sortBlue_end-sortBlue_start;
    printf("Sort blue time is %lf, (%lu)\n",SortBlue,soBlue);
    if(insertCase3()==1)
    {
        printf("success merge bwt\n");
    }
    else
    {
        printf("failed to merge bwt\n");      
    }
    clock_t LFsearch_t=clock();
    double InsertCase3=(double)(LFsearch_t-insertCase3_t)/CLOCKS_PER_SEC;
    printf("InsertCase3 Time is %lf\n",InsertCase3);
    if(LFsearch()==1)
    {
        printf("success LFsearch\n");
    }
    else
    {
        printf("failed to LFsearch\n");
    }
    clock_t finish=clock();
    double LF=(double)(finish-insertCase3_t)/CLOCKS_PER_SEC;
    printf("LF search time is %lf\n",LF);
	return 0;
}
