// we need to sort the blueTable with the help of SPcode in order to get case3's BWT
// blueTable, blueBound, spCode, spSpecial are available
#include "collect#$.h"
#include "sortBlue.h"
#include "INandOut.h"
uint64_t *blueBound=NULL;
uint64_t *spSpecialIndex=NULL;
const uint64_t MAXU64=-1;
pthread_rwlock_t rwlockPoint;
uint64_t bluePoint=0;
uint64_t tempOffset[THREAD_NUM];
/*
uint64_t tempOffsetLow[THREAD_NUM],tempOffsetHigh[THREAD_NUM];
uint64_t chooseBWT[THREAD_NUM];
*/
/*
uint64_t tempOffset,tempOffsetLow,tempOffsetHigh;
uint64_t chooseBWT;
*/
uint64_t blueSplit[THREAD_NUM+1];
int sortBlue(void)
{
    char *blueBoundPath=getPath("/blueBound");
    FILE *fpBlueBound=fopen(blueBoundPath,"rb");
    free(blueBoundPath);
    blueBound=(uint64_t *)calloc(blueBoundNum,sizeof(uint64_t)); 
    fread(blueBound,sizeof(uint64_t),blueBoundNum,fpBlueBound);
    fclose(fpBlueBound);
    printf("MAXU64=%lu\n",MAXU64 );
    /*
    char *spPath=getPath("/SPcode");
    FILE *fpSP=fopen(spPath,"rb");
    free(spPath);
    spCodeLen+=32;
    uint64_t spCodeSpace=spCodeLen>>5;
    uint64_t spLenMod=spCodeLen&MOD32;
    if(spLenMod) spCodeSpace++;
    //free(spCode);
    spCode=(uint64_t *)calloc(spCodeSpace,sizeof(uint64_t)); 
    fread(spCode,sizeof(uint64_t),spCodeSpace,fpSP);
    fclose(fpSP);
    */
    char *spSpecialIndexPath=getPath("/spSpecialIndex");
    FILE *fpspSpecial=fopen(spSpecialIndexPath,"rb");
    free(spSpecialIndexPath);
    spSpecialIndex=calloc(countRead,sizeof(uint64_t)); 
    fread(spSpecialIndex,sizeof(uint64_t),countRead,fpspSpecial); 
    fclose(fpspSpecial);
    uint64_t i;
    for(i=0;i<countRead;i++)
    {
        printf("spSpecialIndex[%lu]=%lu\n",i,spSpecialIndex[i] );
    }
    
    blueSplit[0]=0,blueSplit[THREAD_NUM]=blueBoundNum;
    pthread_t myThread[THREAD_NUM];
    /*
    clock_t temp_start=clock();
    for(i=1;i<THREAD_NUM;i++)
    {
        int *temp=(int *)calloc(1,sizeof(int));
        *temp=i;
        int check=pthread_create( &myThread[i], NULL, multiSplitBlue, (void*)temp);
        if(check)
        {
            fprintf(stderr,"threadNum:%lu, Error - pthread_create() return code: %d\n",i,check);
            exit(EXIT_FAILURE);
        }
    }
    for(i=1;i<THREAD_NUM;i++)
    {
        pthread_join( myThread[i], NULL);
    }
    clock_t temp_end=clock();
    double temp_time=(double)(temp_end-temp_start)/CLOCKS_PER_SEC;
    printf("BinarySearch time: %lf\n",temp_time );
    for(i=0;i<=THREAD_NUM;i++)
    {
        printf("split blue:%lu\n",blueSplit[i] );
    }
    */
    if(pthread_rwlock_init(&rwlockPoint, NULL))
    {
        printf("fail to create rwlock %lu\n",i);
        exit(1);
    }
    for(i=0;i<THREAD_NUM;i++)
    {
        int *temp=(int *)calloc(1,sizeof(int));
        *temp=i;
        int check=pthread_create( &myThread[i], NULL, multiQuickSort, (void*)temp);
        if(check)
        {
            fprintf(stderr,"threadNum:%lu, Error - pthread_create() return code: %d\n",i,check);
            exit(EXIT_FAILURE);
        }   
        //pthread_join( myThread[i], NULL);
    }
    
    for(i=0;i<THREAD_NUM;i++)
    {
        pthread_join( myThread[i], NULL);
    }
    
    /*
    for(i=0;i<blueBoundNum;i++)//sort each case3 in blueTable
    {
        uint64_t start,end=blueBound[i];
        if(i==0)
        {
            start=0;
        }
        else
        {
            start=blueBound[i-1]+1;
        }
        //printf("block: %lu  [%lu,%lu]\n", i,start,end);
        quickSort(blueTable,start,end,0);  
    }
    */
    free(blueBound);
    free(spCode);
    free(spSpecialIndex);
    /*
    for(i=0;i<blueCapacity;i++)
    {
        uint64_t dimer=blueTable[i]&15;
        if(dimer<4) dbgACGT[dimer]++;
        //printf("blue %c\n","ACGT#$"[dimer] );
    }
    */
    //printf("After sort blue, let's check ACGT:\n");
    //printf("A:%lu, C:%lu, G:%lu, T:%lu\n",dbgACGT[0],dbgACGT[1],dbgACGT[2],dbgACGT[3]); 
    return 1;
}

void swap(uint64_t *a,uint64_t *b)
{
    uint64_t swap=*(a);
    *(a)=*(b);
    *(b)=swap;
    return ;
}
uint64_t convertSP(uint64_t mark)//get the 32mer of the spCode
{
   uint64_t ref_mark=mark>>5;
   uint64_t move=(mark&MOD32)<<1;
   uint64_t res=spCode[ref_mark]<<move;
   uint64_t temp=((uint64_t)0==move?(uint64_t)0:(spCode[ref_mark+1]>>(64-move)));
   res=res|temp;
   return res;
}


void *multiSplitBlue(void *arg)
{
    int num=*(int *)arg;
    free((int *)arg);
    uint64_t segment=blueCapacity/THREAD_NUM;
    uint64_t start=segment*num;
    blueSplit[num]=BinarySearch(start,blueBound,blueBoundNum-1);
}

void *multiQuickSort(void *arg)
{
    clock_t temp_start=clock();
    int num=*(int *)arg;
    free((int *)arg);
    uint64_t i,control;
    /*
    printf("num=%d [%lu,%lu]\n",num,blueSplit[num],blueSplit[num+1] );
    uint64_t start_blue=(blueBoundNum/THREAD_NUM)*num,end_blue=(blueBoundNum/THREAD_NUM)*(num+1);
    if(num+1==THREAD_NUM) end_blue=blueBoundNum;
    */
    //for(i=start_blue;i<end_blue;i++)
    while(1)
    //for(i=num;i<blueBoundNum;i+=THREAD_NUM)//sort each case3 in blueTable
    {
        pthread_rwlock_wrlock(&rwlockPoint);
        i=bluePoint;
        bluePoint++;
        pthread_rwlock_unlock(&rwlockPoint);
        if(i>=blueBoundNum) break;
        uint64_t start,end=blueBound[i];
        if(i==0)
        {
            start=0;
        }
        else
        {
            start=blueBound[i-1]+1;
        }
        //printf("block: [%lu,%lu]\n",start,end );
        //quickSort(blueTable,start,end,0,num);  
        for(control=0;control<1;control++)
            myQsort(blueTable,start,end,num,cmpSP);
        //qsort(&blueTable[start],end-start+1,sizeof(uint64_t),cmpQsort);
    }
    clock_t temp_end=clock();
    double temp_time=(double)(temp_end-temp_start)/CLOCKS_PER_SEC;
    printf("%d thread time: %lf\n",num,temp_time);
}

int cmpSP(const void *a, const void *b, int num)
{
    uint64_t mka=*(uint64_t*)a, mkb=*(uint64_t*)b;
    mka=mka>>4,mkb=mkb>>4;
    uint64_t checka=BinarySearch(mka,spSpecialIndex,countRead-1);
    uint64_t checkb=BinarySearch(mkb,spSpecialIndex,countRead-1);
    int64_t cka,ckb;
    uint64_t cva,cvb;
    //record the comparing lenth in the 'tempOffset'
    while((cva=convertSP(mka))==(cvb=convertSP(mkb)))
    {
        cka=spSpecialIndex[checka]-mka,ckb=spSpecialIndex[checkb]-mkb;
        while(cka<=31||ckb<=31)
        {
            if(cka<ckb) return 1; //a string is bigger than b
            else if(cka>ckb) return -1; //b string is bigger than a
            else if(cka==ckb)
            {
                if(checka==countRead) return 1; //a string is bigger than b
                else if(checkb==countRead) return -1; //b string is bigger than a
            }
            if(cka<=31)
            {
                checka++;
                cka=spSpecialIndex[checka]-mka;
            }
            if(ckb<=31)
            {
                checkb++;
                ckb=spSpecialIndex[checkb]-mkb;
            }
        }
        mka+=32;
        mkb+=32;
    }
    cka=spSpecialIndex[checka]-mka,ckb=spSpecialIndex[checkb]-mkb;
    while(cka<=31||ckb<=31)
    {
        if(cka<=31)
        {
            checka++;
            if(checka==countRead)
            {
              cvb=minusDimer(cvb,cka,1);  
              break;
            } 
            else cvb=minusDimer(cvb,cka,0);
            cka=spSpecialIndex[checka]-mka;
        }
        if(ckb<=31)
        {
            checkb++;
            if(checkb==countRead) 
            {
              cva=minusDimer(cva,ckb,1);   
              break;
            }
            else cva=minusDimer(cva,ckb,0);
            ckb=spSpecialIndex[checkb]-mkb;
        }
    }
    //printf("mka=%lu  mkb=%lu\n",mka,mkb );
    //printf("cva=%lu  cvb=%lu   cva>cvb? %d\n",cva,cvb,cva>cvb );
    return cva>cvb;
}
/*
void quickSort(uint64_t *A, int64_t low, int64_t high, uint64_t offset, int num)
{
    uint64_t lostk[STKSIZ], histk[STKSIZ], offstk[STKSIZ];  
    int stkptr=0;                 
    while(1)
    {
        uint64_t countBWT[6];
        if(low<high)
        {
            int countZero=0;
            uint64_t i;
            for(i=0;i<=5;i++) countBWT[i]=0; //initialize countBWT
            for(i=low;i<=high;i++)
            {
                countBWT[A[i]&15]++;
            }
            uint64_t min=MAXU64;
            int minIndex=-1;
            for(i=0;i<=5;i++)
            {
                if(countBWT[i]==0) countZero++;
                else if(countBWT[i]<min)
                {
                    min=countBWT[i];
                    minIndex=i;
                }
            }
            if(countZero>=5) 
            {
                --stkptr;  
                if (stkptr >= 0) 
                {  
                  low  = lostk[stkptr];  
                  high = histk[stkptr];  
                  offset = offstk[stkptr];
                  continue;           // pop subarray from stack  
                }  
                else return;
            }
            chooseBWT[num]=minIndex;// need not to be stored in stack
            //end
            uint64_t offsetLow,offsetHigh;
            tempOffsetLow[num]=MAXU64,tempOffsetHigh[num]=MAXU64;
            int64_t p=partition(A,low,high,offset,num);

            //printf("low:%lu  p:%lu  high:%lu\n",low,p,high );
            offsetLow=tempOffsetLow[num]+offset;
            offsetHigh=tempOffsetHigh[num]+offset;
            //printf("offsetLow:%lu  offsetHigh:%lu\n",offsetLow,offsetHigh );
            if(p-low>=high-p)
            {
                if(low<p-1)
                {
                    lostk[stkptr] = low;  
                    histk[stkptr] = p-1;  
                    offstk[stkptr] = offsetLow;
                    ++stkptr;  
                }
                if(p+1<high)
                {
                    low=p+1;
                    offset=offsetHigh;
                    continue;
                }
            }
            else
            {
                if(p+1<high)
                {
                    lostk[stkptr] = p+1;  
                    histk[stkptr] = high;  
                    offstk[stkptr] = offsetHigh;
                    ++stkptr;
                }
                if(low<p-1)
                {
                    high=p-1;
                    offset=offsetLow;
                    continue;
                }  
            }
            --stkptr;  
            if (stkptr >= 0) 
            {  
              low  = lostk[stkptr];  
              high = histk[stkptr];  
              offset = offstk[stkptr];
              continue;           // pop subarray from stack  
            }  
            else  return; 
        }
        else return;
    }
    return ;
}
int64_t partition(uint64_t *A, uint64_t low, uint64_t high, uint64_t offset, int num)
{
    uint64_t pivot=choosePivot(A,low,high,num);
    uint64_t i;
    //if(pivot==((uint64_t)-1)) printf("error!\n");
    uint64_t pivotValue=A[pivot];
    //printf("choose character:%c\n","ACGT#$"[pivotValue&15] );
    swap(&A[pivot],&A[high]);
    uint64_t storeIndex=low;
    for(i=low;i<high;i++)
    {
        //if(cmpSP(&pivotValue,&A[i],offset))//if pivotValue > A[i]
        if(cmpSP(&pivotValue,&A[i],offset,num))
        {
            swap(&A[i],&A[storeIndex]);
            storeIndex++;
            //printf("storeIndex=%lu\n",storeIndex );
            if(tempOffset[num]<tempOffsetLow[num]) tempOffsetLow[num]=tempOffset[num];
        }
        else
        {
            if(tempOffset[num]<tempOffsetHigh[num]) tempOffsetHigh[num]=tempOffset[num];
        }
    }
    swap(&A[high],&A[storeIndex]);
    return storeIndex;
}
uint64_t choosePivot(uint64_t *A, uint64_t low, uint64_t high, int num)// can not choose the right one now!!
{
    uint64_t i;
    for(i=low;i<=high;i++)
    {
        if((A[i]&15)==chooseBWT[num])
        {
            return i;
        }
    }
    return -1;
}
*/
/*
int cmpSP(const void *a, const void *b, uint64_t offset)
{
    uint64_t mka=*(uint64_t*)a, mkb=*(uint64_t*)b;
    mka=mka>>4,mkb=mkb>>4;
    uint64_t offsetBonus=offset<<5;
    mka+=offsetBonus;
    mkb+=offsetBonus;
    uint64_t checka=BinarySearch(mka,spSpecialIndex,countRead-1);
    uint64_t checkb=BinarySearch(mkb,spSpecialIndex,countRead-1);
    int64_t cka,ckb;
    uint64_t cva,cvb;
    tempOffset=0;
    //record the comparing lenth in the 'tempOffset'
    while((cva=convertSP(mka))==(cvb=convertSP(mkb)))
    {
        tempOffset++; //offset is added for they are the same
        cka=spSpecialIndex[checka]-mka,ckb=spSpecialIndex[checkb]-mkb;
        while(cka<=31||ckb<=31)
        {
            if(cka<ckb) return 1; //a string is bigger than b
            else if(cka>ckb) return -1; //b string is bigger than a
            else if(cka==ckb)
            {
                if(checka==countRead) return 1; //a string is bigger than b
                else if(checkb==countRead) return -1; //b string is bigger than a
            }
            if(cka<=31)
            {
                checka++;
                cka=spSpecialIndex[checka]-mka;
            }
            if(ckb<=31)
            {
                checkb++;
                ckb=spSpecialIndex[checkb]-mkb;
            }
        }
        mka+=32;
        mkb+=32;
    }
    cka=spSpecialIndex[checka]-mka,ckb=spSpecialIndex[checkb]-mkb;
    while(cka<=31||ckb<=31)
    {
        if(cka<=31)
        {
            checka++;
            if(checka==countRead)
            {
              cvb=minusDimer(cvb,cka,1);  
              break;
            } 
            else cvb=minusDimer(cvb,cka,0);
            cka=spSpecialIndex[checka]-mka;
        }
        if(ckb<=31)
        {
            checkb++;
            if(checkb==countRead) 
            {
              cva=minusDimer(cva,ckb,1);   
              break;
            }
            else cva=minusDimer(cva,ckb,0);
            ckb=spSpecialIndex[checkb]-mkb;
        }
    }
    //printf("mka=%lu  mkb=%lu\n",mka,mkb );
    //printf("cva=%lu  cvb=%lu   cva>cvb? %d\n",cva,cvb,cva>cvb );
    return cva>cvb;
}

void quickSort(uint64_t *A, int64_t low, int64_t high, uint64_t offset)
{
    uint64_t countBWT[6];
    if(low<high)
    {
        int countZero=0;
        uint64_t i;
        for(i=0;i<=5;i++) countBWT[i]=0;
        for(i=low;i<=high;i++)
        {
            countBWT[A[i]&15]++;
        }
        uint64_t min=MAXU64;
        int minIndex=-1;
        for(i=0;i<=5;i++)
        {
            if(countBWT[i]==0) countZero++;
            else if(countBWT[i]<min)
            {
                min=countBWT[i];
                minIndex=i;
            }
        }
        if(countZero>=5) return ; // the BWT are all the same
        chooseBWT=minIndex;// need not to be stored in stack
        //end
        uint64_t offsetLow,offsetHigh;
        tempOffsetLow=MAXU64,tempOffsetHigh=MAXU64;
        int64_t p=partition(A,low,high,offset);
        //printf("low:%lu  p:%lu  high:%lu\n",low,p,high );
        offsetLow=tempOffsetLow+offset;
        offsetHigh=tempOffsetHigh+offset;
        //printf("offsetLow:%lu  offsetHigh:%lu\n",offsetLow,offsetHigh );
        quickSort(A,low,p-1,offsetLow);
        quickSort(A,p+1,high,offsetHigh);
    }
    return ;
}
int64_t partition(uint64_t *A, uint64_t low, uint64_t high, uint64_t offset)
{
    uint64_t pivot=choosePivot(A,low,high);
    //if(pivot==((uint64_t)-1)) printf("error!\n");
    uint64_t pivotValue=A[pivot];
    //printf("choose character:%c\n","ACGT#$"[pivotValue&15] );
    swap(&A[pivot],&A[high]);
    uint64_t storeIndex=low;
    uint64_t i;
    for(i=low;i<high;i++)
    {
        //if(cmpSP(&pivotValue,&A[i],offset))//if pivotValue > A[i]
        if(cmpSP(&pivotValue,&A[i],offset))
        {
            swap(&A[i],&A[storeIndex]);
            storeIndex++;
            //printf("storeIndex=%lu\n",storeIndex );
            if(tempOffset<tempOffsetLow) tempOffsetLow=tempOffset;
        }
        else
        {
            if(tempOffset<tempOffsetHigh) tempOffsetHigh=tempOffset;
        }
    }
    swap(&A[high],&A[storeIndex]);
    return storeIndex;
}
uint64_t choosePivot(uint64_t *A, uint64_t low, uint64_t high)// can not choose the right one now!!
{
    uint64_t i;
    for(i=low;i<=high;i++)
    {
        if((A[i]&15)==chooseBWT)
        {
            return i;
        }
    }
    return -1;
}
*/
int cmpQsort(const void *a, const void *b)
{
    uint64_t mka=*(uint64_t*)a, mkb=*(uint64_t*)b;
    mka=mka>>4,mkb=mkb>>4;
    uint64_t checka=BinarySearch(mka,spSpecialIndex,countRead-1);
    uint64_t checkb=BinarySearch(mkb,spSpecialIndex,countRead-1);
    int64_t cka,ckb;
    uint64_t cva,cvb;
    //record the comparing lenth in the 'tempOffset'
    while((cva=convertSP(mka))==(cvb=convertSP(mkb)))
    {
        cka=spSpecialIndex[checka]-mka,ckb=spSpecialIndex[checkb]-mkb;
        while(cka<=31||ckb<=31)
        {
            if(cka<ckb) return 1; //a string is bigger than b
            else if(cka>ckb) return -1; //b string is bigger than a
            else if(cka==ckb)
            {
                if(checka==countRead) return 1; //a string is bigger than b
                else if(checkb==countRead) return -1; //b string is bigger than a
            }
            if(cka<=31)
            {
                checka++;
                cka=spSpecialIndex[checka]-mka;
            }
            if(ckb<=31)
            {
                checkb++;
                ckb=spSpecialIndex[checkb]-mkb;
            }
        }
        mka+=32;
        mkb+=32;
    }
    cka=spSpecialIndex[checka]-mka,ckb=spSpecialIndex[checkb]-mkb;
    while(cka<=31||ckb<=31)
    {
        if(cka<=31)
        {
            checka++;
            if(checka==countRead)
            {
              cvb=minusDimer(cvb,cka,1);  
              break;
            } 
            else cvb=minusDimer(cvb,cka,0);
            cka=spSpecialIndex[checka]-mka;
        }
        if(ckb<=31)
        {
            checkb++;
            if(checkb==countRead) 
            {
              cva=minusDimer(cva,ckb,1);   
              break;
            }
            else cva=minusDimer(cva,ckb,0);
            ckb=spSpecialIndex[checkb]-mkb;
        }
    }
    //printf("mka=%lu  mkb=%lu\n",mka,mkb );
    //printf("cva=%lu  cvb=%lu   cva>cvb? %d\n",cva,cvb,cva>cvb );
    return cva>cvb;
}

void myQsort(uint64_t *A, int64_t lo, int64_t hi, int num,
    int comp(const void *, const void *, int num))
{
    //char *lo, *hi;
    //char *mid;
    //char *loguy, *higuy;
    size_t size;
    uint64_t lostk[STKSIZ], histk[STKSIZ];  
    int stkptr;                 
    stkptr = 0;                 //initialize stack  
    uint64_t countBWT[6],i;
    int chooseBWT;
    uint64_t pivot;
    while(lo<hi) //fake recurse (stack not empty)
    {
        //if(num==0) printf("num=%d,stkptr=%d,lo=%lu,hi=%lu\n",num,stkptr,lo,hi );
        int countZero=0;
        for(i=0;i<=5;i++) countBWT[i]=0; //initialize countBWT
        for(i=lo;i<=hi;i++)
        {
            countBWT[A[i]&15]++;
        }
        uint64_t min=MAXU64;
        int minIndex=-1;
        for(i=0;i<=5;i++)
        {
            if(countBWT[i]==0) countZero++;
            else if(countBWT[i]<min)
            {
                min=countBWT[i];
                minIndex=i;
            }
        }
        if(countZero>=5) 
        {
            --stkptr;  
            if (stkptr >= 0) 
            {  
              lo  = lostk[stkptr];  
              hi = histk[stkptr];  
              continue;           // pop subarray from stack  
            }  
            else return;
        }
        chooseBWT=minIndex;// need not to be stored in stack
        for(i=lo;i<=hi;i++)
        {
            if((A[i]&15)==chooseBWT)
            {
                pivot=i;
                break;
            }
        }
        uint64_t pivotValue=A[pivot];
        //printf("choose character:%c\n","ACGT#$"[pivotValue&15] );
        swap(&A[pivot],&A[hi]);
        int64_t storeIndex=lo;
        for(i=lo;i<hi;i++)
        {
            //if(cmpSP(&pivotValue,&A[i],offset))//if pivotValue > A[i]
            if(cmpSP(&pivotValue,&A[i],num))
            {
                swap(&A[i],&A[storeIndex]);
                storeIndex++;
            }
        }
        swap(&A[hi],&A[storeIndex]);
        if ( storeIndex - lo >= hi - storeIndex )//deal with the smaller block first
        {
            if (lo < storeIndex-1)
            {  
              lostk[stkptr] = lo;  
              histk[stkptr] = storeIndex-1;  
              ++stkptr;  
            }                           // save big recursion for later   
            if (storeIndex + 1 < hi) 
            {  
              lo = storeIndex + 1;  
              continue;           // do small recursion  
            }
        }  
        else
        {
            if (storeIndex + 1 < hi) 
            {  
              lostk[stkptr] = storeIndex + 1;  
              histk[stkptr] = hi;  
              ++stkptr;               // save big recursion for later   
            } 
            if (lo < storeIndex-1) 
            {  
              hi = storeIndex-1;  
              continue;           // do small recursion  
            } 
        }
        --stkptr;
        if (stkptr >= 0) 
        {
          lo = lostk[stkptr];  
          hi = histk[stkptr];  
          continue;           // pop subarray from stack  
        }
        else  return; 
    }
}