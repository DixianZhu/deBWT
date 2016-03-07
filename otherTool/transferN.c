//transfer N into ACGT
//split the chromosome 1 into 10 segments
#include <stdio.h>
#include <zlib.h>
#include <time.h>
#include <inttypes.h>
#include "kseq.h" 
char randTable[15][4]={{},{},{},{},{'A','C','G','T'},{'A','C','G'},{'A','T','G'},{'T','C','G'},{'A','T','C'},{'A','T'},{'C','G'},{'T','G'},{'A','C'},
		{'C','T'},{'A','G'}};
unsigned int trans[256];
int randLen[15]={0,0,0,0,4,3,3,3,3,2,2,2,2,2,2};
KSEQ_INIT(gzFile, gzread);
void usage(void);
int main(int argc, char *argv[])
{
    trans['N']=trans['n']=4;
    trans['V']=trans['v']=5;
    trans['D']=trans['d']=6;
    trans['B']=trans['b']=7;
    trans['H']=trans['h']=8;
    trans['W']=trans['w']=9;
    trans['S']=trans['s']=10;
    trans['K']=trans['k']=11;
    trans['M']=trans['m']=12;
    trans['Y']=trans['y']=13;
    trans['R']=trans['r']=14;
    trans['A']=trans['a']=0;
    trans['C']=trans['c']=1;
    trans['G']=trans['g']=2;
    trans['T']=trans['t']=3;
	if(argc!=3) usage(),exit(1);
   uint64_t i,seg=0;
   char *source=NULL,*output=NULL;
   source=argv[1];
   output=argv[2]; 
   gzFile fp = gzopen(source,"r");
   //free(source);
   kseq_t *seq=kseq_init(fp);
   FILE *plf=fopen(output,"w");
   //free(output);
   char Buf_LF[71];
   srand(time(NULL));
   uint64_t minLen=30000000000;
   while(kseq_read(seq)>=0)
   {
      uint64_t count=0;
      unsigned int num;
      fprintf(plf,">%s\n", seq->name.s);
      for(i=0;seq->seq.s[i]!='\0';i++,seg++,count++)
      {
         if(count==70) 
         {
            Buf_LF[count]='\n';
            fputs(Buf_LF,plf);
            count=0;
         }
         if(trans[(int)seq->seq.s[i]]>3)
	 	 {
            seq->seq.s[i] = randTable[trans[(int)seq->seq.s[i]]][rand()%randLen[trans[(int)seq->seq.s[i]]]];
	 	 }  
         Buf_LF[count]=seq->seq.s[i];
      }
	  if(i<minLen) minLen=i;
      for(i=0;i<count;i++)
      {
        fprintf(plf, "%c",Buf_LF[i] );
      }
      fprintf(plf, "\n" );
      printf("seg=%lu\n",seg);
   }
   kseq_destroy(seq);
   gzclose(fp);
   fclose(plf);
   printf("minLen=%lu\n",minLen);
   return 0;
}

void usage(void)
{
  fprintf(stderr,"usage:\n");
  fprintf(stderr,"./transferN [sourcefile] [outputfile]\n");
  return;
}
