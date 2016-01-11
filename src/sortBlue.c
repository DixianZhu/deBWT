// we need to sort the blueTable with the help of SPcode in order to get case3's BWT
// blueTable, blueBound, spCode, spSpecial are available
#include "collect#$.h"
#include "sortBlue.h"
#include "INandOut.h"
uint64_t *blueBound=NULL;
uint64_t *spSpecialIndex=NULL;
pthread_rwlock_t rwlockPoint;
uint64_t bluePoint=0;
int sortBlue(void **arg)
{
    char *bin=(char *)arg[0];
    uint64_t THREAD_NUM=(uint64_t )arg[1];
    char *blueBoundPath=getPath(bin,"/blueBound");
    FILE *fpBlueBound=fopen(blueBoundPath,"rb");
    blueBound=(uint64_t *)calloc(blueBoundNum,sizeof(uint64_t)); 
    fread(blueBound,sizeof(uint64_t),blueBoundNum,fpBlueBound);
    fclose(fpBlueBound);
    remove(blueBoundPath);
    free(blueBoundPath);
    char *spSpecialIndexPath=getPath(bin,"/spSpecialIndex");
    FILE *fpspSpecial=fopen(spSpecialIndexPath,"rb");
    spSpecialIndex=calloc(countRead+1,sizeof(uint64_t)); 
    fread(spSpecialIndex,sizeof(uint64_t),countRead,fpspSpecial); 
    fclose(fpspSpecial);
    remove(spSpecialIndexPath);
    free(spSpecialIndexPath);
    uint64_t i;
    pthread_t myThread[THREAD_NUM];
    
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
    free(blueBound);
    free(spCode);
    free(spSpecialIndex);
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

void *multiQuickSort(void *arg)
{
    clock_t temp_start=clock();
    int num=*(int *)arg;
    free((int *)arg);
    uint64_t i;
    while(1)
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
        myQsort(blueTable,start,end,cmpSP);
        //qsort(&blueTable[start],end-start+1,sizeof(uint64_t),cmpSP);
    }
    clock_t temp_end=clock();
    double temp_time=(double)(temp_end-temp_start)/CLOCKS_PER_SEC;
    //printf("%d thread time: %lf\n",num,temp_time);
    return (void *)NULL;
}

int cmpSP(const void *a, const void *b)
{
    uint64_t mka=*(uint64_t*)a, mkb=*(uint64_t*)b;
    mka=mka>>4,mkb=mkb>>4;
    uint64_t checka=BinarySearch(mka,spSpecialIndex,countRead-1);
    uint64_t checkb=BinarySearch(mkb,spSpecialIndex,countRead-1);
    uint64_t cka,ckb;
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
                if(checka==countRead-1) return 1; //a string is bigger than b
                else if(checkb==countRead-1) return -1; //b string is bigger than a
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
              //break;
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
              //break;
            }
            else cva=minusDimer(cva,ckb,0);
            ckb=spSpecialIndex[checkb]-mkb;
        }
    }
    //printf("mka=%lu  mkb=%lu\n",mka,mkb );
    //printf("cva=%lu  cvb=%lu   cva>cvb? %d\n",cva,cvb,cva>cvb );
    return cva>cvb?1:-1;
}

void myQsort(uint64_t *A, int64_t lo, int64_t hi,
    int comp(const void *, const void *))
{
    //char *lo, *hi;
    //char *mid;
    //char *loguy, *higuy;
    size_t size;
    uint64_t lostk[STKSIZ], histk[STKSIZ];  
    int stkptr;                 
    stkptr = 0;                 //initialize stack  
    int64_t countBWT[6],i;
    int chooseBWT;
    uint64_t pivot=0;
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
            if(cmpSP(&pivotValue,&A[i])>0)
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
