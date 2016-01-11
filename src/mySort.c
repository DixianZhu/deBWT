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
//void decode(uint64_t obj);

void quickSort(uint64_t (*base)[2], uint64_t num, int comp(const void *, const void *));
void swapTwo(const void *a, const void *b);
int cmpKmer(const void *a, const void *b);
void *multiDistri(void *arg);
pthread_rwlock_t *rwlockBkt;
uint64_t *countKmer=NULL;
uint64_t (*hashKmer)[2]=NULL;
uint64_t trans[256];
uint64_t totalKmerNum;
uint64_t *segCount=NULL;
uint64_t *readBuf=NULL;
uint64_t readNum;
uint64_t tempMove=(32-BUCKET_LENGTH)<<1;
int mySort(void **arg)
{
	char *bin=(char *)arg[0];
	uint64_t THREAD_NUM=(uint64_t )arg[1];
	segCount=(uint64_t *)calloc(THREAD_NUM,sizeof(uint64_t));
	trans['A']=trans['a']=0;
    trans['C']=trans['c']=1;
    trans['N']=trans['G']=trans['g']=2;
    trans['T']=trans['t']=3;
    trans['#']=4;// patch for special characters
    trans['$']=5;
    uint64_t i;
    //////////////////////////////////////kmercounting///////////////////////////////////////////
	countKmer=(uint64_t *)calloc(BUCKET_CAPACITY,sizeof(uint64_t));
	char *seqBuf=(char *)calloc(KMER_LENGTH_PlusOne,sizeof(char));
	uint64_t bufferSize=BUFFERSIZE<<1;
	uint64_t occBuf;
	uint64_t *writeBuf=(uint64_t *)calloc(bufferSize,sizeof(uint64_t));
	uint64_t bufPoint=0;
	FILE *fpKmer=NULL, *fp1=NULL;
	char *outPath=getPath(bin,"/out");
	char *kmerPath=getPath(bin,"/kmerInfo");
	fp1=fopen(outPath,"r");
    if(fp1==NULL) printf("failed open out!\n"), exit(1);
	fpKmer=fopen(kmerPath,"wb");
    if(fpKmer==NULL) printf("failed open kmerInfo!\n"), exit(1);
	free(outPath);
	time_t countTime_start=time(0);
	printf("start kmercounting\n");
	while(fscanf(fp1,"%s%lu",seqBuf,&occBuf)!=EOF)
	{
                //if(dbgIO%1000==0) printf("dbgIO=%lu\n",dbgIO);
                //if(dbgIO>=3431775000) printf("%s: %lu\n",seqBuf,occBuf);
		int i,j;
		uint64_t seq=0,dimer;
		unsigned int move;
		for(i=0,j=31;i<BUCKET_LENGTH;i++,j--)
		{
			move=j<<1;
			dimer=trans[(int)seqBuf[i]];
			seq=seq|(dimer<<move);
		}
		uint64_t tempSeq=seq>>(tempMove);
      	countKmer[tempSeq]++;
		//Go on writing out the seq
		for(;i<KMER_LENGTH_PlusOne;i++,j--)
		{
			move=j<<1;
			dimer=trans[(int)seqBuf[i]];
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
    free(seqBuf);
	free(writeBuf);
    fclose(fpKmer);
	fclose(fp1);
	time_t countTime_end=time(0);
	fprintf(stderr, "read kmercounting results time (txt transfer): %ld\n",countTime_end-countTime_start );
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
	clock_t dis_start=clock();
	readBuf=(uint64_t *)calloc(bufferSize,sizeof(uint64_t));
	hashKmer=(uint64_t (*)[2])calloc(totalKmerNum<<1,sizeof(uint64_t));
	fpKmer=fopen(kmerPath,"rb");
	rwlockBkt=(pthread_rwlock_t*)calloc(BUCKET_CAPACITY,sizeof(pthread_rwlock_t));
	for(i=0;i<BUCKET_CAPACITY;i++)
	{
		if(pthread_rwlock_init(&rwlockBkt[i], NULL))
		{
			printf("fail to create rwlock %lu\n",i);
			exit(1);
		}	
	}
	pthread_t myThread[THREAD_NUM];
	printf("start kmer distributing\n");
	time_t multiDistri_start=time(0);
	while((readNum=fread(readBuf,sizeof(uint64_t),bufferSize,fpKmer))>0)
	{           
		for(i=0;i<THREAD_NUM;i++)
		{
			void **tt=(void **)calloc(2,sizeof(void *));
			tt[0]=(void *)i;
			tt[1]=(void *)THREAD_NUM;
			int check=pthread_create( &myThread[i], NULL, multiDistri, (void*)tt);
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
                
	}
	fclose(fpKmer);
	clock_t dis_end=clock();
	double dis_duration=(double)(dis_end-dis_start)/CLOCKS_PER_SEC;
	time_t multiDistri_end=time(0);
	long multiDistri_dur=(multiDistri_end-multiDistri_start);
	printf("distributing time = %lf (%ld)\n",dis_duration,multiDistri_dur);
	free(rwlockBkt);
	time_t start=time(0);
	printf("start sorting\n");
	clock_t mul_start=clock();
	//////////////////////////////////////multiThread//////////////////////////////////////////////////
	//pthread_attr_t attr[THREAD_NUM];
	for(i=0;i<THREAD_NUM;i++)
	{
		void **tt=(void **)calloc(2,sizeof(void *));
		tt[0]=(void *)i;
		tt[1]=(void *)THREAD_NUM;
		int check=pthread_create( &myThread[i], NULL, multiThreadSort, (void*) tt);
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
	void **argT=(void **)arg;
	uint64_t num=(uint64_t )argT[0];
	uint64_t THREAD_NUM=(uint64_t)argT[1];
	free(argT);
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
	//printf("thread:%lu  [%lu,%lu]\n",num,low,up);
	for(i=low;i<up;i++)
	{
		//quickSort(&hashKmer[countKmer[i]],countKmer[i+1]-countKmer[i],cmp);
		//quickSort((void *)&hashKmer[countKmer[i]],countKmer[i+1]-countKmer[i],2*sizeof(uint64_t),cmp);
		qsort(hashKmer[countKmer[i]],(countKmer[i+1]-countKmer[i]),2*sizeof(uint64_t),cmpKmer);
	}
	if(i==BUCKET_CAPACITY-1)
	{
		//quickSort(&hashKmer[countKmer[i]],totalKmerNum-countKmer[i],cmp);
		//quickSort((void *)&hashKmer[countKmer[i]],totalKmerNum-countKmer[i],2*sizeof(uint64_t),cmp);
		qsort(hashKmer[countKmer[i]],(totalKmerNum-countKmer[i]),2*sizeof(uint64_t),cmpKmer);
	}
	
	clock_t finish=clock();
	double duration=(double)(finish-start)/CLOCKS_PER_SEC;
	return (void *)NULL;
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
			swapTwo((void *)base[lo],(void *)base[mid]);
		}
		if(comp((void *)base[lo],(void *)base[hi]) > 0)
		{
			swapTwo((void *)base[lo],(void *)base[hi]);
		}
		if(comp((void *)base[mid],(void *)base[hi]) > 0)
		{
			swapTwo((void *)base[mid],(void *)base[hi]);	
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
			swapTwo((void *)base[loguy],(void *)base[higuy]);
			if(mid==higuy) mid=loguy;
		}
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
void swapTwo(const void *a, const void *b)
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

int cmpKmer(const void *a, const void *b)
{
	uint64_t va=((uint64_t (*)[2])a)[0][0],vb=((uint64_t (*)[2])b)[0][0];
	//printf("va=%lu; vb=%lu\n",va,vb );
	if(va<vb) return -1;
	else if(va==vb) return 0;
	else return 1;
}
/*
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
*/
void *multiDistri(void *arg)
{
	void **argT=(void **)arg;
	uint64_t num=(uint64_t )argT[0];
	uint64_t THREAD_NUM=(uint64_t)argT[1];
	free(argT);
	uint64_t segment=(readNum>>1)/THREAD_NUM;
	uint64_t start=(num*segment)<<1,i;
	uint64_t end;
	if(num<THREAD_NUM-1)
	{
		end=((num+1)*segment)<<1;
	}
	else//num==THREAD_NUM-1
	{
		end=readNum;
	}
        //printf("thread:%d  (%lu, %lu)\n",num,start,end);
	for(i=start;i<end;i+=2)
	{
		uint64_t tempSeq=readBuf[i]>>(tempMove);
		uint64_t tempOcc=readBuf[i+1];
		//countKmer[tempSeq] stands for the position
		pthread_rwlock_wrlock(&rwlockBkt[tempSeq]);
		countKmer[tempSeq]--;//because it is upper bound(cannot reach)
		hashKmer[countKmer[tempSeq]][0]=readBuf[i];
		hashKmer[countKmer[tempSeq]][1]=readBuf[i+1];
		pthread_rwlock_unlock(&rwlockBkt[tempSeq]);
	}
    return (void *)NULL;
}
