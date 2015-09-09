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
pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
char *dbgBWT=NULL;
char cis[4]={'A','C','G','T'};
#define MOD32 31
#define MOD64 63
KSEQ_INIT(gzFile, gzread);
int collect(char *path)
{
   gzFile fp = gzopen(path,"r");
   kseq_t *seq=kseq_init(fp);
   uint64_t i,count;
   int j;
   while(kseq_read(seq)>=0)
   {
      printf("name: %s\n", seq->name.s);
      printf("counting the seq length:%ld\n",seq->seq.l);
      ref_length+=seq->seq.l;
      countRead++;
   }
   kseq_destroy(seq);
   gzclose(fp);
   fp = gzopen(path,"r");
   if(fp==NULL) {printf("can not open ref file\n");exit(0);}
   seq=kseq_init(fp);
   ref_length+=countRead;//add # and $ to the ref
   BWTLEN=ref_length;//store the bwt length
   ref_length+=32;//add 32 'T' to the tail
   printf("BWTLEN=%lu\n",BWTLEN);      
   compress_length=ref_length>>5;
   if(ref_length&MOD32) compress_length++;
   reference=(uint64_t*)calloc(compress_length,sizeof(uint64_t));//success free
   special=(uint64_t*)calloc(countRead,sizeof(uint64_t));//success free
   uint64_t seg=0;
   count=0;
   while(kseq_read(seq)>=0)
   {
      printf("name: %s\n", seq->name.s);
      for(i=0;seq->seq.s[i]!='\0';i++,seg++)
      {
         //if(seq->seq.s[i]=='n') printf("exist n!\n");
         reference[seg>>5]=reference[seg>>5]|(trans[seq->seq.s[i]]<<((31-(seg&MOD32))<<1));
         ACGT[trans[seq->seq.s[i]]]++;
      }
      printf("i=%lu\n",i);
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
   char *refPath=getPath("/reference");
   FILE *pref=fopen(refPath,"wb");
   free(refPath);
   if(pref==NULL)
   {
    printf("fail to open reference!\n");
    exit(1);
   }
   fwrite(reference,sizeof(uint64_t),compress_length,pref);
   fclose(pref);
   char *specialSAPath=getPath("/specialModule/specialSA");
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

   qsort(SA,countRead*(KMER_LENGTH+1),sizeof(uint64_t),cmp);
   //////////////////////////////////////////dbg/////////////////////////////////////////////////// 
   /*
   uint64_t *dbgSA=(uint64_t *)calloc(BWTLEN,sizeof(uint64_t));
   for(i=0;i<BWTLEN;i++) dbgSA[i]=i;
   qsort(dbgSA,BWTLEN,sizeof(uint64_t),cmp);
   dbgBWT=(char *)calloc(BWTLEN,sizeof(char));
   for(i=0;i<BWTLEN;i++)
   {
      if(dbgSA[i]>0)
      {
        dbgBWT[i]=getCharacter(dbgSA[i]-1);
        if(dbgSA[i]==special[0]) dbgBWT[i]='#';
      }
      else
      {
        dbgBWT[i]='$';
      }
      printf("%lu: %c\n",i,dbgBWT[i] );
   }
   free(dbgSA);
   free(dbgBWT);
   //////////////////////////////////////////////////////////////////////////////////////////////////////
   /*print the kmer of SA*/
   divideKmer(seeKMER(SA),SA);
   /*
   ///////////////////////////////////////////RestoreRef///////////////////////////////////////////////////////////
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
   free(reference);
   free(SA);
   free(special);
   return 1;
}
void *generateSpecialSA(void *arg)
{
  uint64_t *bound=(uint64_t *)arg,i;
  int j;
  uint64_t start=bound[0]*(KMER_LENGTH+1);
  for(i=bound[0];i<bound[1];i++)//This part can be changed to multiThread, while it's fast enough
  {
   for(j=-KMER_LENGTH;j<=0;j++,start++)
   {
     SA[start]=special[i]+j;
   }
  }
  free(bound);
}
uint64_t convert(uint64_t mark)//get the 32mer of the ref
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
   int64_t cka,ckb;
   uint64_t cva,cvb;
   while((cva=convert(mka))==(cvb=convert(mkb)))
   {
    if((cka=special[checka]-mka)<=31)//check the # or $ point
    {
      checka++; //move the index of #
    }
    else
    {
      cka=-1;
    }
    if((ckb=special[checkb]-mkb)<=31)
    {
      checkb++;
    }
    else
    {
      ckb=-1;
    }
    if(cka>=0||ckb>=0)
    {
      printf("cka=%ld,ckb=%ld\n",cka,ckb );
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
        if(checka==countRead) return 1;//a string is bigger than b
        else if(checkb==countRead) return -1;//b string is bigger than a
      }
    }
    mka+=32;
    mkb+=32;
   }
   cka=-1,ckb=-1;
   if((cka=special[checka]-mka)<=31)//check the # or $ point
   {
     checka++; //move the index of #
     uint64_t index=cka;
     if(checka==countRead) cvb=minusDimer(cvb,index,1);
     else cvb=minusDimer(cvb,index,0); 
   }
   if((ckb=special[checkb]-mkb)<=31)
   {
     checkb++;
     uint64_t index=ckb;
     if(checkb==countRead) cva=minusDimer(cva,index,1);
     else cva=minusDimer(cva,index,0);
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
  if(flag==0&&dimer>=1)//# minus 1
  {
    dimer-=1;
  }
  else if(flag==1&&dimer>=2)//$ minus 2
  {
    dimer-=2;
  }
  else return cv;
  res=cv&(~((uint64_t)3<<move));
  res=res|(dimer<<move);
  return res;
}

char (*seeKMER(uint64_t *SA))[KMER_LENGTH+1]
{
  uint64_t len=countRead*(KMER_LENGTH+1);
  uint64_t i;
  int j;
  char *head$Path=getPath("/specialModule/head$");
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
  
  for(i=0;i<len;i++)
  {
    char *kmer=(char *)calloc(KMER_LENGTH+1,sizeof(char));//success free
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
          bwtSA[bwtIndex]=bwtSA[bwtIndex]|(trans[kmer[j]]<<move);
          if(kmer[j+1]=='#'||kmer[j+1]=='$') flag=1;
        }
      }
      //decode(bwtSA[bwtIndex]);
      bwtIndex++;
    }
    free(kmer);
  }
  
  char *bwtPath=getPath("/specialModule/specialBwt");
  FILE *fpbwt=fopen(bwtPath,"wb");
  free(bwtPath);
  fwrite(bwt,sizeof(char),countRead*KMER_LENGTH,fpbwt);
  fwrite(bwtSA,sizeof(uint64_t),countRead*KMER_LENGTH,fpbwt);
  free(bwtSA);
  free(bwt);
  fclose(fpbwt);
  return res;
}
/*
void *multiSeeKmer(void *arg)//This method is wrong, for it cannot use multiThread, desolate it.
{
  void **para=(void **)arg;
  uint64_t *pi=(uint64_t *)para[0];
  uint64_t num=*pi;
  free(pi);
  uint64_t len=(uint64_t)para[1];
  uint64_t *bwtIndex=(uint64_t *)para[2];
  char (*res)[KMER_LENGTH+1];
  res=(char (*)[KMER_LENGTH+1])para[3];
  char *bwt=(char *)para[4];
  uint64_t *bwtSA=(uint64_t *)para[5];
  char kmer[KMER_LENGTH+1];
  uint64_t segLen=len/THREAD_NUM;
  uint64_t start=num*segLen,end;
  if(num<THREAD_NUM-1)
  {
    end=(num+1)*segLen;
  }
  else
  {
    end=len;
  }
  uint64_t i;
  int j,k;
  //printf("num=%lu, [%lu,%lu]\n",num,start,end );
  for(i=start;i<end;i++)
  {
    //printf("%-10lu",SA[i]);
    uint64_t temp=convert(SA[i]);
    uint64_t precursor=(convert(SA[i]-1)>>62)&3;
    uint64_t check=BinarySearch(SA[i],special,countRead-1);
    for(j=31,k=0;k<=KMER_LENGTH;j--,k++)
    {
      int move=j<<1;
      uint64_t dimer=(temp>>move)&3;
      if(SA[i]+k==special[check]) 
      {
        if(check<countRead-1) kmer[k]='#';
        else kmer[k]='$';
      }
      else kmer[k]=cis[dimer];
      res[i][k]=kmer[k];
    }
    //printf("i=%lu\n",i );  
    if(kmer[KMER_LENGTH]!='#'&&kmer[KMER_LENGTH]!='$') //check this 
    {
      pthread_mutex_lock( &mutex );
      bwt[*bwtIndex]=cis[precursor];
      unsigned int flag=0;
      if(kmer[0]=='#'||kmer[0]=='$') flag=1;
      //printf("flag=%d\n",flag);
      for(j=0,k=31;j<KMER_LENGTH;j++,k--)
      {
        unsigned int move=k<<1;
        if(flag)
        {
          bwtSA[*bwtIndex]=bwtSA[*bwtIndex]|((uint64_t)3<<move); //'T'
        }
        else 
        {
          bwtSA[*bwtIndex]=bwtSA[*bwtIndex]|(trans[kmer[j]]<<move);
          if(kmer[j+1]=='#'||kmer[j+1]=='$') flag=1;
        }
      }
      decode(bwtSA[*bwtIndex]);
      printf("bwtIndex=%lu\n",*bwtIndex );
      (*bwtIndex)++;
      pthread_mutex_unlock( &mutex );
    }
  }
  free(para);
}
*/
void divideKmer(char (*pK)[KMER_LENGTH+1],uint64_t *SA)//if reads shorter than 32, there will be some mistakes
{
  char *headSharpPath=getPath("/specialModule/head#"); // change to binary file
  FILE *fph=fopen(headSharpPath,"wb");
  free(headSharpPath);
  char *tailSharpPath=getPath("/specialModule/tail#"); // change to binary file
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
        headBuffer[headI]=headBuffer[headI]|(trans[pK[i][j]]<<move);
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
        tailBuffer[tailI]=tailBuffer[tailI]|(trans[pK[i][j]]<<move);
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
  char *specialBranchPath=getPath("/specialModule/specialBranch");
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

char *getPath(char *name)
{
  char *path=(char *)calloc(PATH_LEN,sizeof(char));
  strcpy(path,PATH);
  strcat(path,name);
  return path;
}