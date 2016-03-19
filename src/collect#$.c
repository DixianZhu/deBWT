//collect all kmers whose precursor or successor is # and store out into a file    :head#   tail#
//generate the branch of the special kmer    : specialBranch
//store out reference     : reference
//store out specialKmer   :  specialKmer (may not be used)
//store out the special point SA : specialSA  
#include "collect#$.h"
#include "kseq.h" 
#include "INandOut.h"
#include "generateSP.h"
uint64_t *special,*reference;//special stand for the index of special character
uint64_t countRead=0;
uint64_t BWTLEN,ref_length=0,compress_length;
uint64_t specialBranchNum=0;
uint64_t ACGT[4]={0,0,0,0};
uint64_t *SA=NULL;
uint64_t multiFlag=1;
uint64_t branchNum=0;
uint64_t *specialHash=NULL;
uint64_t *invHash=NULL;
pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
//char *dbgBWT=NULL;
char cis[4]={'A','C','G','T'};
//char dbgStr1040[30];
#define MOD32 31
#define MOD64 63
KSEQ_INIT(gzFile, gzread);
void *collect(void *arg_collect)
{
   void **arg=(void **) arg_collect;
   char *path=(char *)arg[2];
   uint64_t THREAD_NUM=(uint64_t)arg[1];
   char *bin=(char *)arg[0];
   printf("path: %s\n",path );
   gzFile fp = gzopen(path,"r");
   kseq_t *seq=kseq_init(fp);
   uint64_t i=0,count; 
   while(kseq_read(seq)>=0)
   {
      //printf("name: %s (%lu)\n", seq->name.s, i++);
      //printf("counting the seq length:%ld\n",seq->seq.l);
      if(seq->seq.l<=32) 
      {
  	     fprintf(stderr,"Length <= 32!\n");
	       exit(1);
      }
      ref_length+=seq->seq.l;
      countRead++;
   }
   //printf("countRead=%lu\n",countRead);
   //exit(0);
   kseq_destroy(seq);
   gzclose(fp);
   fp = gzopen(path,"r");
   if(fp==NULL) {fprintf(stderr,"can not open ref file\n");exit(0);}
   seq=kseq_init(fp);
   ref_length+=countRead;//add # and $ to the ref
   BWTLEN=ref_length;//store the bwt length
   ref_length+=32;//add 32 'T' to the tail
   printf("BWTLEN=%lu\n",BWTLEN);      
   compress_length=ref_length>>5;
   if(ref_length&MOD32) compress_length++;
   reference=(uint64_t*)calloc(compress_length,sizeof(uint64_t));//success free
   special=(uint64_t*)calloc(countRead+1,sizeof(uint64_t));//success free
   uint64_t seg=0;
   count=0;
   while(kseq_read(seq)>=0)
   {
      //printf("name: %s\n", seq->name.s);
      for(i=0;seq->seq.s[i]!='\0';i++,seg++)
      {
	 /*
         if(seq->seq.s[i]!='A'&&seq->seq.s[i]!='C'&&seq->seq.s[i]!='G'&&seq->seq.s[i]!='T') 
 	 {
	   printf("%c!\n",seq->seq.s[i]);
	   exit(1);
	 }
	 */
         reference[seg>>5]=reference[seg>>5]|(trans[(int)(seq->seq.s[i])]<<((31-(seg&MOD32))<<1));
         ACGT[trans[(int)(seq->seq.s[i])]]++;
      }
      //printf("i=%lu\n",i);
      reference[seg>>5]=reference[seg>>5]|(trans['T']<<((31-(seg&MOD32))<<1));
      special[count]=seg;
      seg++;
      count++;
   }
   for(i=0;i<32;i++,seg++)
   {
    reference[seg>>5]=reference[seg>>5]|(trans['T']<<((31-(seg&MOD32))<<1));//add 32 'T' to the tail
   }
   printf("A:%lu, C:%lu, G:%lu, T:%lu\n",ACGT[0],ACGT[1],ACGT[2],ACGT[3]);
   for(i=1;i<4;i++)
   {
    ACGT[i]+=ACGT[i-1];
   }
   for(i=3;i>0;i--)
   {
    ACGT[i]=ACGT[i-1];
   }
   ACGT[0]=0;
   kseq_destroy(seq);
   gzclose(fp);
   char *refPath=getPath(bin,"/reference");
   FILE *pref=fopen(refPath,"wb");
   free(refPath);
   if(pref==NULL)
   {
    printf("fail to open reference!\n");
    exit(1);
   }
   fwrite(reference,sizeof(uint64_t),compress_length,pref);
   fclose(pref);
   char *specialSAPath=getPath(bin,"/specialSA");
   pref=fopen(specialSAPath,"wb");
   free(specialSAPath);
   fwrite(special,sizeof(uint64_t),countRead,pref);
   fclose(pref);
   SA=(uint64_t*)calloc(countRead*(KMER_LENGTH+1),sizeof(uint64_t));//success free
   /*
   ///////////////////////////collect special SA single thread//////////////////////////////////////////
   count=0;
   for(i=0;i<countRead;i++)//This part can be changed to multiThread, while it's fast enough
   {
    for(j=-KMER_LENGTH;j<=0;j++,count++)
    {
      SA[count]=special[i]+j;
    }
   }
   */
   ///////////////////////////collect special SA multi-thread///////////////////////////////////////////
   uint64_t boundLen=countRead/THREAD_NUM;
   pthread_t myThread[THREAD_NUM];
   for(i=0;i<THREAD_NUM;i++)
   {
    uint64_t *bound=(uint64_t *)calloc(2,sizeof(uint64_t));
    bound[0]=i*(boundLen);
    if(i<THREAD_NUM-1)
    {
      bound[1]=(i+1)*boundLen;
    }
    else
    {
      bound[1]=countRead;
    }
    int check=pthread_create( &myThread[i], NULL, generateSpecialSA, (void*) bound);
    if(check)
    {
        fprintf(stderr,"threadNum:%lu, Error - pthread_create() return code: %d\n",i,check);
        exit(EXIT_FAILURE);
    }
   }
   for(i=0;i<THREAD_NUM;i++)
   {
    pthread_join( myThread[i], NULL);
   }
   uint64_t temp_capa=countRead*(KMER_LENGTH+1);
   qsort(SA,temp_capa,sizeof(uint64_t),cmp);
   /*
   time_t hash_start=time(0);
   //////////////////////////////////////////A dbg module for LF-back search/////////////////////////////
   specialHash=(uint64_t *)calloc(countRead,sizeof(uint64_t));
   invHash=(uint64_t *)calloc(countRead,sizeof(uint64_t));
   uint64_t special_p=0;
   for(i=0;i<temp_capa;i++)
   {
      uint64_t j;
      for(j=0;j<countRead;j++)
      {
         if(special[j]==SA[i])
         {
            invHash[j]=special_p;
            specialHash[special_p++]=j;
         }
      }
   }
   time_t hash_end=time(0);
   printf("hash_dbg time cost: %lu\n",hash_end-hash_start);
   //////////////////////////////////////////////////////////////////////////////////////////////////////
   */
   /*print the kmer of SA*/
   divideKmer(seeKMER(SA,bin),SA,bin);
   /*
   ///////////////////////////////////////////RestoreRef/////////////////////////////////////////////////
   clock_t start=clock();
   char restoreBuf[71];
   int pointBuf=0;
   int64_t pointSpecial=countRead-2;
   char *restorePath=getPath("/restoreRef");
   FILE *fpRestore=fopen(restorePath,"w");
   free(restorePath);
   restoreBuf[70]='\0';
   int64_t restore_index;
   for(restore_index=BWTLEN-2;restore_index>=0;restore_index--)
   {
     if(pointSpecial>=0 && restore_index==special[pointSpecial])
     {
       restoreBuf[pointBuf++]='#';  
       pointSpecial--;
     }
     else
     {
       restoreBuf[pointBuf++]=getCharacter(restore_index);
     }
     if(pointBuf>=70)
     {
       fprintf(fpRestore, "%s\n",restoreBuf );
       pointBuf=0;
     } 
   }
   restoreBuf[pointBuf++]='$';
   for(i=0;i<pointBuf;i++)
   {
     fprintf(fpRestore, "%c", restoreBuf[i]);
   }
   fprintf(fpRestore, "\n" );
   fclose(fpRestore);
   clock_t finish=clock();
   double duration=(double)(finish-start)/CLOCKS_PER_SEC;
   printf("Restore_Ref time %lf\n",duration);
   //////////////////////////////////////////////////////////////////////////////////////////////////////
   */
   printf("branchNum=%lu\n",branchNum);
   free(reference);
   free(SA);
   free(special);
   return (void *)1;
}
void *generateSpecialSA(void *arg)
{
  uint64_t *bound=(uint64_t *)arg,i;
  int j;
  uint64_t start=bound[0]*(KMER_LENGTH+1);
  for(i=bound[0];i<bound[1];i++)
  {
   for(j=-KMER_LENGTH;j<=0;j++,start++)
   {
     SA[start]=special[i]+j;
   }
  }
  free(bound);
  return (void*)NULL;
}
uint64_t convert(uint64_t mark)//get the 32 mer of the ref
{
   uint64_t ref_mark=mark>>5;
   uint64_t move=(mark&MOD32)<<1;
   uint64_t res=reference[ref_mark]<<move;
   uint64_t temp=((uint64_t)0==move?(uint64_t)0:(reference[ref_mark+1]>>(64-move)));
   res=res|temp;
   return res;
}

int cmp(const void *a,const void *b)//the length of each read must longer than 32, or the special will have bug
{
   uint64_t mka=*(uint64_t*)a,mkb=*(uint64_t*)b;
   uint64_t checka=BinarySearch(mka,special,countRead-1),checkb=BinarySearch(mkb,special,countRead-1);
   uint64_t cka,ckb;
   uint64_t cva,cvb;
   while((cva=convert(mka))==(cvb=convert(mkb)))
   {
    cka=special[checka]-mka,ckb=special[checkb]-mkb;
    while((cka<=31)|(ckb<=31))
    {
      //printf("cka=%ld,ckb=%ld\n",cka,ckb );
      if(cka<ckb)//a string is bigger than b
      {
        return 1;
      }
      else if(cka>ckb)//a string is smaller than b
      {
        return -1;
      }
      else if(cka==ckb)
      {
        if(checka==countRead-1) return 1;//a string is bigger than b
        else if(checkb==countRead-1) return -1;//b string is bigger than a
      }
      if(cka<=31)
      {
          checka++;
          cka=special[checka]-mka;
      }
      if(ckb<=31)
      {
          checkb++;
          ckb=special[checkb]-mkb;
      }
    }
    mka+=32;
    mkb+=32;
   }
   cka=special[checka]-mka,ckb=special[checkb]-mkb;
   while(cka<=31||ckb<=31)
   {
     if(cka<=31)
     {
		checka++; //move the index of #
		if(checka==countRead) cvb=minusDimer(cvb,cka,1);
        else cvb=minusDimer(cvb,cka,0); 
        cka=special[checka]-mka;
     }          
     if(ckb<=31)
     {
		checkb++;
        if(checkb==countRead) cva=minusDimer(cva,ckb,1);
        else cva=minusDimer(cva,ckb,0);
        ckb=special[checkb]-mkb;
     }
   }
   return cva>cvb?1:-1;
}

uint64_t BinarySearch(uint64_t mk, uint64_t *target, int64_t up)//Search the # and $ point
{
  int64_t low=0,mid=0;
  while(low<=up)
  {
     mid=(low+up)>>1;
     if(mk<target[mid])   up=mid-1;
     else if(mk>target[mid])   low=mid+1;
     if(mk==target[mid])    return mid;
  }
  return low;
}

uint64_t minusDimer(uint64_t cv,uint64_t index,int flag)//minus the index dimer 
{
  unsigned int move=(31-index)<<1;
  uint64_t dimer=(cv>>move)&3;
  uint64_t res=0;
  if(dimer>=3)
  {
    if(flag==0)//# minus 1
    {
      dimer-=1;
    }
    else if(flag==1)//$ minus 2
    {
      dimer-=2;
    }
  }
  else return cv;
  res=cv&(~((uint64_t)3<<move));
  res=res|(dimer<<move);
  return res;
}

char **seeKMER(uint64_t *SA, char *bin)
{
  uint64_t len=countRead*(KMER_LENGTH+1);
  uint64_t i;
  int j;
  char *head$Path=getPath(bin,"/head$");
  FILE *fp$h=fopen(head$Path,"w");
  free(head$Path);
  uint64_t $h=convert(0);
  for(j=31;j>=31-KMER_LENGTH+1;j--)
  {
    unsigned int move=j<<1;
    uint64_t dimer=($h>>move)&3;
    fprintf(fp$h, "%c",cis[dimer]);
  }
  fprintf(fp$h, "\n" );
  fclose(fp$h);
  char (*res)[KMER_LENGTH+1];
  res=(char (*)[KMER_LENGTH+1])calloc(len*(KMER_LENGTH+1),sizeof(char));//success free
  char *bwt=(char *)calloc(countRead*KMER_LENGTH,sizeof(char));
  uint64_t *bwtSA=(uint64_t*)calloc(countRead*KMER_LENGTH,sizeof(uint64_t));//success free
  uint64_t bwtIndex=0;
  char *kmer=(char *)calloc(KMER_LENGTH+1,sizeof(char));//success free
  /* multiSeeKmer process
  void **para=NULL;
  pthread_t myThread[THREAD_NUM];
  for(i=0;i<THREAD_NUM;i++)
  {
    para=calloc(6,sizeof(void *));
    uint64_t *pi=calloc(1,sizeof(uint64_t));
    *pi=i;
    para[0]=(void *)pi;
    para[1]=(void *)len;
    para[2]=(void *)&bwtIndex;
    para[3]=(void *)res;
    para[4]=(void *)bwt;
    para[5]=(void *)bwtSA;
    printf("i(out)=%lu\n", i);
    int check=pthread_create( &myThread[i], NULL, multiSeeKmer, (void *)para);
    if(check)
    {
      fprintf(stderr,"threadNum:%lu, Error - pthread_create() return code: %d\n",i,check);
      exit(EXIT_FAILURE);
    }
  }
  for(i=0;i<THREAD_NUM;i++)
  {
    multiFlag=THREAD_NUM;
    pthread_join( myThread[i], NULL);
  }
  */
  /*
  uint64_t dbgSeq=0;
  dbgSeq=dbgSeq|(((uint64_t)1<<62)-1)-3;
  printf("target:\n");
  printf("A#");
  for(j=0;j<29;j++)
  {
    printf("%c",dbgStr1040[j]);
  }
  putchar('\n');
  */
  for(i=0;i<len;i++)
  {
    //printf("%-10lu",SA[i]);
    uint64_t temp=convert(SA[i]);
    uint64_t precursor=(convert(SA[i]-1)>>62)&3;
    uint64_t check=BinarySearch(SA[i],special,countRead-1);
    for(j=31;j>=31-KMER_LENGTH;j--)
    {
      int move=j<<1;
      uint64_t dimer=(temp>>move)&3;
      if(SA[i]+31-j==special[check]) 
      {
        if(check<countRead-1) kmer[31-j]='#';
        else kmer[31-j]='$';
      }
      else kmer[31-j]=cis[dimer];
      res[i][31-j]=kmer[31-j];
    }
    if(kmer[KMER_LENGTH]!='#'&&kmer[KMER_LENGTH]!='$') //if #$ in the end, we needn't them as patch
    {
      bwt[bwtIndex]=cis[precursor];
      unsigned int k,flag=0;
      if(kmer[0]=='#'||kmer[0]=='$') flag=1;
      //printf("flag=%d\n",flag);
      for(j=0,k=31;j<KMER_LENGTH;j++,k--)
      {
        unsigned int move=k<<1;
        if(flag)
        {
          bwtSA[bwtIndex]=bwtSA[bwtIndex]|((uint64_t)3<<move); //'T'
        }
        else 
        {
          bwtSA[bwtIndex]=bwtSA[bwtIndex]|(trans[(int)kmer[j]]<<move);
          if(kmer[j+1]=='#'||kmer[j+1]=='$') flag=1;
        }
      }
      /*
      if(bwtSA[bwtIndex]==dbgSeq)
      {
        printf("dbgSeq: %c|%s\n",bwt[bwtIndex],kmer);
      }
      */
      //decode(bwtSA[bwtIndex]);
      bwtIndex++;
    }
  }
  free(kmer);
  char *bwtPath=getPath(bin,"/specialBwt");
  FILE *fpbwt=fopen(bwtPath,"wb");
  free(bwtPath);
  fwrite(bwt,sizeof(char),countRead*KMER_LENGTH,fpbwt);
  fwrite(bwtSA,sizeof(uint64_t),countRead*KMER_LENGTH,fpbwt);
  free(bwtSA);
  free(bwt);
  fclose(fpbwt);
  return (char **)res;
}
void divideKmer(char **pk, uint64_t *SA, char *bin)
{
  char (*pK)[KMER_LENGTH+1]=(char (*)[KMER_LENGTH+1])pk;
  char *headSharpPath=getPath(bin,"/head#"); // change to binary file
  FILE *fph=fopen(headSharpPath,"wb");
  free(headSharpPath);
  char *tailSharpPath=getPath(bin,"/tail#"); // change to binary file
  FILE *fpt=fopen(tailSharpPath,"wb");
  free(tailSharpPath);
  uint64_t len=countRead*(KMER_LENGTH+1);
  uint64_t i,j;
  uint64_t *headBuffer=(uint64_t *)calloc(BUFFERSIZE,sizeof(uint64_t)),headI=0;
  uint64_t *tailBuffer=(uint64_t *)calloc(BUFFERSIZE,sizeof(uint64_t)),tailI=0;
  int k;
  unsigned int move;
  for(i=0;i<len;i++)
  {
    if(pK[i][0]=='#')
    {
      for(j=1,k=31;j<=KMER_LENGTH;j++,k--)
      {
        move=k<<1;
        headBuffer[headI]=headBuffer[headI]|(trans[(int)pK[i][j]]<<move);
        //fprintf(fph, "%c", pK[i][j]);
      }
      headI++;
      if(headI>=BUFFERSIZE)
      {
        fwrite(headBuffer,sizeof(uint64_t),BUFFERSIZE,fph);
        headI=0;
        free(headBuffer);
        headBuffer=(uint64_t *)calloc(BUFFERSIZE,sizeof(uint64_t));
      }
      //fprintf(fph, "\n");
    }
    else if(pK[i][KMER_LENGTH]=='#'||pK[i][KMER_LENGTH]=='$')
    {
      for(j=0,k=31;j<=KMER_LENGTH-1;j++,k--)
      {
        move=k<<1;
        tailBuffer[tailI]=tailBuffer[tailI]|(trans[(int)pK[i][j]]<<move);
        //fprintf(fpt, "%c", pK[i][j]);
      }
      tailI++;
      if(tailI>=BUFFERSIZE)
      {
        fwrite(tailBuffer,sizeof(uint64_t),BUFFERSIZE,fpt);
        tailI=0;
        free(tailBuffer);
        tailBuffer=(uint64_t *)calloc(BUFFERSIZE,sizeof(uint64_t));
      }
      //fprintf(fpt, "\n");
    }
  }
  if(headI>0)
  {
    fwrite(headBuffer,sizeof(uint64_t),headI,fph);
  }
  if(tailI>0)
  {
    fwrite(tailBuffer,sizeof(uint64_t),tailI,fpt);
  }
  fclose(fph);
  fclose(fpt);
  free(headBuffer);
  free(tailBuffer);
  char *specialBranchPath=getPath(bin,"/specialBranch");
  FILE *specialBranch=fopen(specialBranchPath,"wb");
  free(specialBranchPath);
  uint64_t low=0,up=0;
  int flag=0;
  uint64_t *branchBuffer=(uint64_t *)calloc(BUFFERSIZE,sizeof(uint64_t)),branchI=0;
  for(i=1;i<len;i++)//store the branch node mark
  {
    if(pK[i][KMER_LENGTH]!='#'&&pK[i][KMER_LENGTH]!='$')//because the $ kmer will not have muti-branch, 
    //it needn't any other judge
    {
      if(compare(&SA[low],&SA[i]))
      {
        up=i;
        if(pK[i][KMER_LENGTH]!=pK[i-1][KMER_LENGTH])
        {
          flag=1;
        }
      }
      else
      {
        if(up>low&&flag)
        {
          for(j=low;j<=up;j++)
          {
            //fprintf(specialBranch, "%lu\n", SA[j]); // change to binary file
            branchBuffer[branchI++]=SA[j];
            specialBranchNum++;
            if(branchI>=BUFFERSIZE)
            {
              fwrite(branchBuffer,sizeof(uint64_t),BUFFERSIZE,specialBranch);
              branchI=0;
            }
          }
          flag=0;
        }
        low=i,up=i;
      }
    }
    else
    {
      if(up>low&&flag) 
      {
         for(j=low;j<=up;j++)
         {
           //fprintf(specialBranch, "%lu\n", SA[j]); // change to binary file
           branchBuffer[branchI++]=SA[j];
           specialBranchNum++;
           if(branchI>=BUFFERSIZE)
           {
             fwrite(branchBuffer,sizeof(uint64_t),BUFFERSIZE,specialBranch);
	         branchNum+=BUFFERSIZE;
             branchI=0;
           }
         }
         flag=0;
      }
      low=i,up=i;
    }
  }
  if(branchI>0)
  {
    fwrite(branchBuffer,sizeof(uint64_t),branchI,specialBranch);
    branchNum+=branchI;
  }
  free(pK);
  fclose(specialBranch);
  free(branchBuffer);
}
int compare(const void *a,const void *b)
{
  uint64_t mka=*(uint64_t*)a,mkb=*(uint64_t*)b;
  uint64_t checka=BinarySearch(mka,special,countRead-1),checkb=BinarySearch(mkb,special,countRead-1);
  uint64_t cka,ckb;
  uint64_t cva,cvb;  
  cva=convert(mka);
  cvb=convert(mkb);
  cka=special[checka]-mka,ckb=special[checkb]-mkb;
  if(checka==countRead-1) checka=1;
  else checka=0;
  if(checkb==countRead-1) checkb=1;
  else checkb=0;
  return compareI(cva,cvb,cka,ckb,checka,checkb);
}
int compareI(uint64_t cva,uint64_t cvb,uint64_t cka,uint64_t ckb,uint64_t checka,uint64_t checkb)
{
  if(cka<=KMER_LENGTH-1)//check the # or $ point
  {
    uint64_t index=cka;
    cvb=minusDimer(cvb,index,checka);
  }
  if(ckb<=KMER_LENGTH-1)
  {
    uint64_t index=ckb;
    cva=minusDimer(cva,index,checkb);
  }
  uint64_t eliminate=~(ELIMINATE);
  cva=cva&eliminate;
  cvb=cvb&eliminate;
  return cva==cvb?1:0;
}

char *getPath(char *DIR, char *name)
{
  char *path=(char *)calloc(PATH_LEN,sizeof(char));
  strcpy(path,DIR);
  strcat(path,name);
  return path;
}
