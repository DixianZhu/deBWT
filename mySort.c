#include "collect#$.h"
#define BUCKET_LENGTH 12 //need less than KMER_LENGTH_PlusOne
#define BUCKET_CAPACITY (1<<(BUCKET_LENGTH<<1))
/*
void swap(char *a,char *b, size_t size);
int cmp(const void *a, const void *b);
static void shortsort (char *lo, char *hi, size_t width, int comp(const void *, const void *));
void quickSort(void *base, uint64_t num, size_t width, int comp(const void *, const void *));
*/
void *multiThreadSort(void *arg);
void decode(uint64_t obj);

void quickSort(uint64_t (*base)[2], uint64_t num, int comp(const void *, const void *));
void swap(const void *a, const void *b);
int cmp(const void *a, const void *b);

uint64_t *countKmer=NULL;
uint64_t (*hashKmer)[2]=NULL;
uint64_t trans[256];
uint64_t totalKmerNum;
uint64_t segCount[THREAD_NUM];
uint64_t multiFlag=1;
int main()
{
	trans['A']=trans['a']=0;
    trans['C']=trans['c']=1;
    trans['N']=trans['G']=trans['g']=2;
    trans['T']=trans['t']=3;
    trans['#']=4;// patch for special characters
    trans['$']=5;
    /*
    uint64_t testA=0,testB=100;
    printf("testA=%lu, testB=%lu\n",testA,testB );
    swap((char *)&testA,(char *)&testB,sizeof(uint64_t));
    printf("testA=%lu, testB=%lu\n",testA,testB );
    */
   	uint64_t i;
    uint64_t demo[5][2]={{3,3},{2,2},{1,1},{4,4},{5,5}};
    quickSort(demo,5,cmp);
    for(i=0;i<5;++i)
    {
    	printf("%lu,%lu\n",demo[i][0],demo[i][1] );
    }
    //////////////////////////////////////kmercounting///////////////////////////////////////////
	countKmer=(uint64_t *)calloc(BUCKET_CAPACITY,sizeof(uint64_t));
	char seqBuf[KMER_LENGTH_PlusOne];
	uint64_t bufferSize=BUFFERSIZE<<1;
	uint64_t occBuf;
	uint64_t *writeBuf=(uint64_t *)calloc(bufferSize,sizeof(uint64_t));
	uint64_t bufPoint=0;
	FILE *fpKmer=NULL, *fp1=NULL;
	char *outPath=getPath("/out");
	char *kmerPath=getPath("/kmerInfo");
	fp1=fopen(outPath,"r");
	fpKmer=fopen(kmerPath,"wb");
	free(outPath);
	uint64_t tempMove=(32-BUCKET_LENGTH)<<1;
	printf("start kmercounting\n");
	while(fscanf(fp1,"%s%lu",seqBuf,&occBuf)!=EOF)
	{
		int i,j;
		uint64_t seq=0,dimer;
		unsigned int move;
		for(i=0,j=31;i<BUCKET_LENGTH;i++,j--)
		{
			move=j<<1;
			dimer=trans[seqBuf[i]];
			seq=seq|(dimer<<move);
		}
		uint64_t tempSeq=seq>>(tempMove);
		countKmer[tempSeq]++;
		//Go on writing out the seq
		for(;i<KMER_LENGTH_PlusOne;i++,j--)
		{
			move=j<<1;
			dimer=trans[seqBuf[i]];
			seq=seq|(dimer<<move);
		}
		writeBuf[bufPoint++]=seq;
		writeBuf[bufPoint++]=occBuf;
		if(bufPoint>=bufferSize)
		{
			fwrite(writeBuf,sizeof(uint64_t),bufferSize,fpKmer);
			bufPoint=0;
		}
	}
	if(bufPoint)
	{
		fwrite(writeBuf,sizeof(uint64_t),bufPoint,fpKmer);
	}
	free(writeBuf);
	fclose(fpKmer);
	fclose(fp1);
	/////////////////////////////////////////kmer distributing///////////////////////////////////////////////
	//printf("countKmer[0]=%lu\n",countKmer[0]);
	for(i=1;i<BUCKET_CAPACITY;i++)
	{
		countKmer[i]=countKmer[i-1]+countKmer[i];//upper bound
		//printf("countKmer[%lu]=%lu\n",i,countKmer[i]);
	}
	totalKmerNum=countKmer[BUCKET_CAPACITY-1];
	segCount[0]=0;
	for(i=1;i<THREAD_NUM;i++)
	{
		uint64_t locateNum=i*(totalKmerNum/THREAD_NUM);
		uint64_t locate=BinarySearch(locateNum,countKmer,BUCKET_CAPACITY-1);
		segCount[i]=locate;
	}

	uint64_t *readBuf=(uint64_t *)calloc(bufferSize,sizeof(uint64_t));
	uint64_t readNum;
	hashKmer=(uint64_t (*)[2])calloc(totalKmerNum<<1,sizeof(uint64_t));
	fpKmer=fopen(kmerPath,"rb");
	printf("start kmer distributing\n");
	while((readNum=fread(readBuf,sizeof(uint64_t),bufferSize,fpKmer))>0)
	{
		//printf("readNum=%lu\n",readNum );
		for(i=0;i<readNum;i+=2)//distribute the kmer
		{
			uint64_t tempSeq=readBuf[i]>>(tempMove);
			uint64_t tempOcc=readBuf[i+1];
			//countKmer[tempSeq] stands for the position
			countKmer[tempSeq]--;//because it is upper bound(cannot reach)
			hashKmer[countKmer[tempSeq]][0]=readBuf[i];
			hashKmer[countKmer[tempSeq]][1]=readBuf[i+1];
		}
	}
	fclose(fpKmer);
	printf("start sorting\n");
	time_t start=time(0);
	clock_t mul_start=clock();
	//////////////////////////////////////multiThread//////////////////////////////////////////////////
	pthread_t myThread[THREAD_NUM];
	//pthread_attr_t attr[THREAD_NUM];
	for(i=0;i<THREAD_NUM;i++)
	{
		printf("i=%lu\n",i );
		//pthread_attr_init(&attr[i]);
		//pthread_attr_setscope(&attr[i], PTHREAD_SCOPE_SYSTEM);
		int check=pthread_create( &myThread[i], NULL, multiThreadSort, (void*) i);
		if(check)
    	{
        	fprintf(stderr,"threadNum:%lu, Error - pthread_create() return code: %d\n",i,check);
        	exit(EXIT_FAILURE);
    	}
    	//pthread_join( myThread[i], NULL);	
	}
	
	for(i=0;i<THREAD_NUM;i++)
	{
		multiFlag=THREAD_NUM;
		pthread_join( myThread[i], NULL);
	}
	
	time_t finish=time(0);
	clock_t mul_end=clock();
	long duration=(finish-start);
	double mul_dur=(double)(mul_end-mul_start)/CLOCKS_PER_SEC;
	printf("sort time is %lf (%ld)\n",mul_dur,duration );
	/////////////////////////////////////end///////////////////////////////////////////////////////////
	/*
	for (i = 0; i < totalKmerNum-1; ++i)
	{
		if(hashKmer[i][0]>hashKmer[i+1][0]) 
		{
			printf("hashKmer[%lu]=%lu, hashKmer[%lu]=%lu\n",i,hashKmer[i][0],i+1,hashKmer[i+1][0]);
		}
	}
	*/
	fpKmer=fopen(kmerPath,"wb");
	fwrite(hashKmer,sizeof(uint64_t)*2,totalKmerNum,fpKmer);
	fclose(fpKmer);
	free(kmerPath);
	free(readBuf);
	free(hashKmer);
	free(countKmer);
	return 0;
}

void *multiThreadSort(void *arg)
{
	clock_t start=clock();
	uint64_t num=(uint64_t )arg;
	uint64_t i,j;
	uint64_t low,up;
	low=segCount[num];
	if(num<THREAD_NUM-1)
	{
		up=segCount[num+1];
	}
	else//num==THREAD_NUM-1
	{
		up=BUCKET_CAPACITY-1;
	}
	printf("thread:%lu  [%lu,%lu]\n",num,low,up);
	for(i=low;i<up;i++)
	{
		quickSort(&hashKmer[countKmer[i]],countKmer[i+1]-countKmer[i],cmp);
		//quickSort((void *)&hashKmer[countKmer[i]],countKmer[i+1]-countKmer[i],2*sizeof(uint64_t),cmp);
		//qsort(hashKmer[countKmer[i]],(countKmer[i+1]-countKmer[i]),2*sizeof(uint64_t),cmp);
	}
	if(i==BUCKET_CAPACITY-1)
	{
		quickSort(&hashKmer[countKmer[i]],totalKmerNum-countKmer[i],cmp);
		//quickSort((void *)&hashKmer[countKmer[i]],totalKmerNum-countKmer[i],2*sizeof(uint64_t),cmp);
		//qsort(hashKmer[countKmer[i]],(totalKmerNum-countKmer[i]),2*sizeof(uint64_t),cmp);
	}
	
	clock_t finish=clock();
	double duration=(double)(finish-start)/CLOCKS_PER_SEC;
	duration=duration/multiFlag;
	printf("thread:%lu  multiFlag=%lu cost:%lf\n",num,multiFlag,duration );
	return ;
}
///////////////////////////////////not used yet////////////////////////////////////////////////

void quickSort(uint64_t (*base)[2], uint64_t num, int comp(const void *, const void *))
{
	uint64_t lo=0, hi=num-1, mid;
	int64_t loguy, higuy;
	uint64_t lostk[STKSIZ], histk[STKSIZ];
	int stkptr=0;
	size_t size;
	if(num<2) return;
	while(1)
	{
		size=hi-lo+1;
		mid=lo+(size>>1);
		if(comp((void *)base[lo],(void *)base[mid]) > 0)
		{
			swap((void *)base[lo],(void *)base[mid]);
		}
		if(comp((void *)base[lo],(void *)base[hi]) > 0)
		{
			swap((void *)base[lo],(void *)base[hi]);
		}
		if(comp((void *)base[mid],(void *)base[hi]) > 0)
		{
			swap((void *)base[mid],(void *)base[hi]);	
		}
		loguy=lo;
		higuy=hi;
		while(1)
		{
			do
			{
				loguy++;
			}while(loguy<=hi && comp((void *)base[loguy],(void *)base[mid])<=0);
			do
			{
				higuy--;
			}while(higuy>mid && comp((void *)base[higuy],(void *)base[mid])>0);
			if(higuy<loguy) 
			{
				break;
			}
			swap((void *)base[loguy],(void *)base[higuy]);
			if(mid==higuy) mid=loguy;
		}
		//swap((void *)base[mid],(void *)base[higuy]);
		while(higuy>lo&&comp((void *)base[higuy],(void *)base[mid])==0)
		{
			higuy--;
		}
		if ( higuy - lo >= hi - loguy )//deal with the smaller block first
        {
        	if (lo < higuy)
        	{
              lostk[stkptr] = lo;  
              histk[stkptr] = higuy;  
              ++stkptr;  
            }                           // save big recursion for later   
            if (loguy < hi) 
            {  
              lo = loguy;  
              continue;           // do small recursion  
            }
        }  
        else
        {
        	if (loguy < hi) 
        	{  
              lostk[stkptr] = loguy;  
              histk[stkptr] = hi;  
              ++stkptr;               // save big recursion for later   
        	} 
        	if (lo < higuy) 
           	{  
              hi = higuy;  
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
void swap(const void *a, const void *b)
{
	uint64_t temp;
	uint64_t (*va)[2]=((uint64_t (*)[2])a),(*vb)[2]=((uint64_t (*)[2])b);
	temp=va[0][0];
	va[0][0]=vb[0][0];
	vb[0][0]=temp;
	temp=va[0][1];
	va[0][1]=vb[0][1];
	vb[0][1]=temp;
}

/*
static void shortsort (char *lo, char *hi, size_t width, int comp(const void *, const void *))  
{  
    char *p, *max;  
    while (hi > lo) 
    { 
      max = lo;  
      for (p = lo+width; p <= hi; p += width) 
      {  
        if (comp(p, max) > 0) 
        {  
            max = p;  
        }  
      }  
      swap(max, hi, width);  
      hi -= width;  
    }  
}
void quickSort(void *base, uint64_t num, size_t width, int comp(const void *, const void *))
{
	char *lo, *hi;
	char *mid;
	char *loguy, *higuy;
	size_t size;
	char *lostk[STKSIZ], *histk[STKSIZ];  
    int stkptr;                 
    if (num < 2 || width == 0)  
      return; 
    stkptr = 0;                 //initialize stack  
    lo = base;  
    hi = lo+ width * (num-1);        // initialize limits   
    while(1) //fake recurse (stack not empty)
    {
    	size=(hi-lo)/width+1; 
    	if(size<=CUTOFF)
    	{
    		shortsort(lo, hi, width, comp);
    	}
    	else
    	{
    		mid = lo + (size>>1) * width;      // find middle element 
    		if (comp(lo, mid) > 0) 
    		{  
            	swap(lo, mid, width);  
	        }  
	        if (comp(lo, hi) > 0) 
	        {  
            	swap(lo, hi, width);  
	        }  
	        if (comp(mid, hi) > 0) 
	        {  
            	swap(mid, hi, width);  
	        }
	        loguy = lo; // traveling pointers for partition step  循环中的游动指针
       		higuy = hi; // traveling pointers for partition step  循环中的游动指针 
       		while(1)
       		{
       			do 
                {
                  loguy += width;  
	            } while (loguy <= hi && comp(loguy, mid) <= 0);     
	            do  
	            {
                  higuy -= width;  
                } while (higuy > mid && comp(higuy, mid) > 0);                 
                if (higuy < loguy)  break;
                swap(loguy, higuy, width);
                if (mid == higuy)  mid = loguy;
       		}
       		//swap(mid,higuy,width);
        	while (higuy > lo && comp(higuy, mid) == 0)
        	{
				higuy -= width;
        	}
        	if ( higuy - lo >= hi - loguy )//deal with the smaller block first
        	{
        		if (lo < higuy)
        		{  
        	      lostk[stkptr] = lo;  
        	      histk[stkptr] = higuy;  
        	      ++stkptr;  
        	    }                           // save big recursion for later   
        	    if (loguy < hi) 
        	    {  
        	      lo = loguy;  
        	      continue;           // do small recursion  
        	    }
        	}  
        	else
        	{
        		if (loguy < hi) 
        		{  
        	      lostk[stkptr] = loguy;  
        	      histk[stkptr] = hi;  
        	      ++stkptr;               // save big recursion for later   
        		} 
        		if (lo < higuy) 
        	   	{  
        	      hi = higuy;  
        	      continue;           // do small recursion  
        	   	} 
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
void swap(char *a,char *b, size_t size)
{
    char temp;
    if(a!=b)
    {
    	while(size--)
    	{
    		temp=*a;
    		*a++=*b;
    		*b++=temp;
    	}
    }
    return ;
}
*/

int cmp(const void *a, const void *b)
{
	uint64_t va=((uint64_t (*)[2])a)[0][0],vb=((uint64_t (*)[2])b)[0][0];
	//printf("va=%lu; vb=%lu\n",va,vb );
	if(va<vb) return -1;
	else if(va==vb) return 0;
	else return 1;
}

char *getPath(char *name)
{
  char *path=(char *)calloc(PATH_LEN,sizeof(char));
  strcpy(path,PATH);
  strcat(path,name);
  return path;
}
void decode(uint64_t obj)
{
	int i,j;
	uint64_t eliminate=3;
	for(i=0,j=31;i<KMER_LENGTH;i++,j--)
	{
		unsigned int move=j<<1;
		printf("%c","ACGT"[(obj>>move)&eliminate]);
	}
	putchar('\n');
}
uint64_t BinarySearch(uint64_t mk, uint64_t *target, int64_t up)//return the upper bound
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


