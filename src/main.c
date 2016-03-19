#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <time.h>
#include <inttypes.h>
#include <string.h>
//#include "kseq.h"
#include "main.h"
#include "collect#$.h"
#define VERSION "1.0.1"
uint64_t trans[256];
int KMER_LENGTH_PlusOne=32;
int KMER_LENFGTH=31;
char *get_bin_dir(char *bin);
void usage(void);
int main(int argc, char *argv[])
{
    trans['A']=trans['a']=0;
    trans['C']=trans['c']=1;
    trans['G']=trans['g']=2;
    trans['T']=trans['t']=3;
    trans['#']=4;// patch for special characters
    trans['$']=5;
    ////////////////////////////GUI control///////////////////////////////
    if(argc>10||argc<6||(argc&1)==1) usage(),exit(1);
    char *source=argv[argc-1];
    char *obj=NULL;
    char *bin=get_bin_dir(argv[0]);
    char *jRoot=NULL;
    uint64_t THREAD_NUM=8;
    int i;
    for(i=1;i<argc-1;i+=2)
    {
        if(strcmp(argv[i],"-o")==0) obj=argv[i+1];
        if(strcmp(argv[i],"-t")==0) 
        {
            THREAD_NUM=atoi(argv[i+1]);
            if(THREAD_NUM==0) fprintf(stderr,"thread number must be a number!\n"), exit(1);
        }
        if(strcmp(argv[i],"-j")==0) jRoot=argv[i+1];
		if(strcmp(argv[i],"-k")==0)
		{
			KMER_LENGTH_PlusOne=atoi(argv[i+1]);
			KMER_LENGTH=KMER_LENGTH_PlusOne-1;
			if(KMER_LENGTH_PlusOne<12||KMER_LENGTH_PlusOne>32) 
				fprintf(stderr, "-k: k-mer length (from 12 to 32, default 32)\n" ),exit(1);
		}
    }
    if(jRoot==NULL)
    {
        fprintf(stderr, "Must specify your jellyfish root to use deBWT\n" );
        exit(1);
    }
    //////////////////////////if obj bwt cannot create, we tell it in advance///////
    FILE *fpbwt=fopen(obj,"wb");
    if(fpbwt==NULL) fprintf(stderr,"cannot create %s!\n",obj), exit(1);
    fclose(fpbwt);
    remove(obj);
    ////////////////////////////////////////////////////////////////////////////////
    char cmd[1024];
    fprintf(stderr,"---------------------------------------------------------------------------\n");
    fprintf(stderr,"run deBWT %s:\n",VERSION);
    fprintf(stderr,"sequence file: %s (shouldn't contain any uncertain char like 'N'...)\n", source);
    fprintf(stderr,"output bwt file: %s\n",obj);
    fprintf(stderr,"take %s as the bin to store some temporary files, which will be removed in the end\n",bin );
    fprintf(stderr,"use %lu-thread\n",THREAD_NUM);
	fprintf(stderr,"k-mer length: %d\n",KMER_LENGTH_PlusOne);
    fprintf(stderr,"try to run jellyfish in %s\n", jRoot);
    /////////////////////////////////////////jellyfish///////////////////////////////////////////////////
    sprintf(cmd, "bash %s/src/kmercounting.sh %s %lu %s %s %d", bin, source, THREAD_NUM, bin, jRoot, KMER_LENGTH_PlusOne);
    if (system(cmd) != 0 ) 
    { 
        fprintf(stderr, "\njellyfish exit abnormally.\n"); 
        exit(1);
    }
    fprintf(stderr,"done!\n");
    fprintf(stderr,"use deBWT construct BWT string ...\n");
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    void *arg_collect[3];
    arg_collect[0]=(void *)bin;
    arg_collect[1]=(void *)THREAD_NUM;
    arg_collect[2]=(void *)source;
    if(mySort((void**)arg_collect)==0) fprintf(stderr,"success sort kmers\n");
    else fprintf(stderr,"failed sort kmers, some error happened!\n"),exit(1);
    /////////////////////////////////////////////////////////////////////////////////////////////////////
    clock_t collect_t=clock();
    time_t mix_start=time(0);
    pthread_t myThread[2];
    
    int check=pthread_create( &myThread[0], NULL, collect, (void*)arg_collect);//need bin, THREAD_NUM and source
    if(check)
    {
        printf("fail to create thread for 'collect'\n");
        exit(EXIT_FAILURE);
    }
    check=pthread_create( &myThread[1], NULL, getKmer, (void*)bin);
    if(check)
    {
        printf("fail to create thread for 'getKmer'\n");
        exit(EXIT_FAILURE);
    }
    pthread_join( myThread[0], NULL);
    pthread_join( myThread[1], NULL);
    clock_t generateBlock_t=clock();
    time_t mix_end=time(0);
    long mix_dur=(mix_end-mix_start);
    double Collect=(double)(generateBlock_t-collect_t)/CLOCKS_PER_SEC;
    printf("collect#$ and getKmer time is %lf, (%lu)\n",Collect,mix_dur);
    if(generateBlocks(bin)==1)
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
    void **arg_generateSP=arg_collect;
    if(generateSP(arg_generateSP)==1)
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
    if(sortBlue(arg_generateSP)==1)
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
    if(insertCase3(obj, bin)==1)
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
    if(LFsearch(bin)==1)
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

void usage(void)
{
    fprintf(stderr,"usage:\n");
    fprintf(stderr,"deBWT [options] reference\n");
    fprintf(stderr,"Please make sure your sequence don't contain any uncertain characters like 'N'\n");
    fprintf(stderr, "options:\n" );
    fprintf(stderr, "-o: output bwt file(binary)\n");
    fprintf(stderr, "-t (optional): maximum thread number(default 8)\n" );
    fprintf(stderr, "-k (optional): k-mer length (from 12 to 32, default 32)\n" );
    fprintf(stderr, "-j: jellyfish directory\n" );
    fprintf(stderr, "reference: sequence in fasta or fastq format\n");
}

char *get_bin_dir(char *bin)
{
    char *end = strrchr(bin, '/');
    bin[end-bin] = '\0';
    return bin;
}
