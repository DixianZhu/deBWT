//transfer N into ACGT
//split the chromosome 1 into 10 segments
#include "collect#$.h"
#include "kseq.h" 
KSEQ_INIT(gzFile, gzread);
int main()
{
   char *bigGenome10Path=getPath("/bigGenome10.fa");
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
         if(seq->seq.s[i]=='N'||seq->seq.s[i]=='n')
            seq->seq.s[i] = "ACGT"[rand()%4];
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
