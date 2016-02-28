//transfer N into ACGT
//split the chromosome 1 into 10 segments
#include "collect#$.h"
#include "kseq.h" 
char randTable[15][4]={{},{},{},{},{'A','C','G','T'},{'A','C','G'},{'A','T','G'},{'T','C','G'},{'A','T','C'},{'A','T'},{'C','G'},{'T','G'},{'A','C'},
		{'C','T'},{'A','G'}};
int randLen[15]={0,0,0,0,4,3,3,3,3,2,2,2,2,2,2};
KSEQ_INIT(gzFile, gzread);
int main()
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
   char *bigGenome10Path=getPath("/bigGenome10.fa");
   //char *bigGenome10Path="/home/dxzhu/1000genome_individuals/bigGenome10.fa";
   gzFile fp = gzopen(bigGenome10Path,"r");
   free(bigGenome10Path);
   kseq_t *seq=kseq_init(fp);
   uint64_t i,seg=0;
   char *bigGenome10FBPath=getPath("/bigGenome10FB.fa");
   FILE *plf=fopen(bigGenome10FBPath,"w");
   free(bigGenome10FBPath);
   char Buf_LF[71];
   srand(time(NULL));
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
         if(trans[seq->seq.s[i]]>3)
	 {
            seq->seq.s[i] = randTable[trans[seq->seq.s[i]]][rand()%randLen[trans[seq->seq.s[i]]]];
	 }  
         Buf_LF[count]=seq->seq.s[i];
      }
      for(i=0;i<count;i++)
      {
        fprintf(plf, "%c",Buf_LF[i] );
      }
      fprintf(plf, "\n" );
      printf("i=%lu\n",i);
   }
   kseq_destroy(seq);
   gzclose(fp);
   fclose(plf);
   return 0;
}
char *getPath(char *name)
{
  char *path=(char *)calloc(PATH_LEN,sizeof(char));
  strcpy(path,PATH);
  strcat(path,name);
  return path;
}
