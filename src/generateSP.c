#include "collect#$.h"
#include "generateSP.h"
#include "INandOut.h"
uint64_t *redSeq=NULL;
uint64_t *redPoint=NULL;
uint64_t *blueTable=NULL;//point to the SP code's index
uint64_t spCodeLen=0;
uint64_t *specialBranch=NULL;
uint64_t *specialSA=NULL;
uint64_t dbgQue=0;
uint64_t *splitIndex=NULL;
uint64_t *spSplit,*spSplitCapacity;
uint64_t *spCode=NULL;
pthread_rwlock_t *rwlockRed;
uint64_t *fusionNode;
uint64_t *fusionMod;
uint64_t **tempSP;
//uint64_t spIndex=0;
int generateSP(void **arg)
{
	char *bin=(char *)arg[0];
	uint64_t THREAD_NUM=(uint64_t )arg[1];
	spSplit=(uint64_t *)calloc(THREAD_NUM,sizeof(uint64_t));
	spSplitCapacity=(uint64_t *)calloc(THREAD_NUM,sizeof(uint64_t));
	fusionNode=(uint64_t *)calloc(THREAD_NUM,sizeof(uint64_t));
	fusionMod=(uint64_t *)calloc(THREAD_NUM,sizeof(uint64_t));
	tempSP=(uint64_t **)calloc(THREAD_NUM,sizeof(uint64_t*));
	printf("the case3num is %lu\n", case3num);
	printf("the blueBoundNum is %lu\n", blueBoundNum);
	printf("the redCapacity is %lu\n", redCapacity);
	printf("the blueCapacity is %lu\n", blueCapacity);
	redSeq=(uint64_t *)calloc(redCapacity,sizeof(uint64_t));
	redPoint=(uint64_t *)calloc(redCapacity,sizeof(uint64_t));
	blueTable=(uint64_t *)calloc(blueCapacity,sizeof(uint64_t));
	char *redSeqPath=getPath(bin,"/redSeq");
	FILE *fpredSeq=fopen(redSeqPath,"rb");
	char *redPointPath=getPath(bin,"/redPoint");
	FILE *fpredPoint=fopen(redPointPath,"rb");
	if(fpredSeq==NULL||fpredPoint==NULL)
	{
		printf("fail to open the red file!\n");
		exit(1);
	}
	fread(redSeq,sizeof(uint64_t),redCapacity,fpredSeq);
	fread(redPoint,sizeof(uint64_t),redCapacity,fpredPoint);
	fclose(fpredSeq);
	fclose(fpredPoint);
	remove(redPointPath);
	remove(redSeqPath);
	free(redPointPath);
	free(redSeqPath);
	uint64_t i;
	for(i=1;i<BLACKCAPACITY;i++)
	{
		if(blackTable[i]==0)
		{
			blackTable[i]=blackTable[i-1];
		}
	}	
	//////////////////////////////go through reference////////////////////////////////////
	char *refPath=getPath(bin,"/reference");
	FILE *fpRef=fopen(refPath,"rb");
	if(fpRef==NULL)
	{
		printf("fail to open the ref file\n");
		exit(1);	
	}
	reference=(uint64_t *)calloc(compress_length,sizeof(uint64_t)); 
	fread(reference,sizeof(uint64_t),compress_length,fpRef);
	fclose(fpRef);
	remove(refPath);
	free(refPath);
	////////////////////////////////////////////////////////////////////////////////////////////////
	char *specialBranchPath=getPath(bin,"/specialBranch");
	FILE *fpSB=fopen(specialBranchPath,"rb");
	specialBranch=(uint64_t *)calloc(specialBranchNum,sizeof(uint64_t));
	if(fpSB==NULL)
	{
		printf("fail to open the specialBranch file\n");
		exit(1);	
	}
	/* replaced by binary file
	for(i=0;i<specialBranchNum;i++)
	{
		fscanf(fpSB,"%lu",&specialBranch[i]);
	}
	*/
	fread(specialBranch,sizeof(uint64_t),specialBranchNum,fpSB);
	fclose(fpSB);
	remove(specialBranchPath);
	free(specialBranchPath);
	qsort(specialBranch,specialBranchNum,sizeof(uint64_t),ascend);
	specialSA=(uint64_t *)calloc(countRead,sizeof(uint64_t));
	char *specialSAPath=getPath(bin,"/specialSA");
	FILE *fpSpecialSA=fopen(specialSAPath,"rb");
	if(fpSpecialSA==NULL)
	{
		printf("fail to open the specialBranch file\n");
		exit(1);	
	}
	fread(specialSA,sizeof(uint64_t),countRead,fpSpecialSA);
	fclose(fpSpecialSA);
	remove(specialSAPath);
	free(specialSAPath);
	///////////////////////////////////split the ref first//////////////////////////////////////////
	splitIndex=(uint64_t *)calloc(THREAD_NUM+1,sizeof(uint64_t));
	splitIndex[0]=0;
	splitIndex[THREAD_NUM]=BWTLEN;
	pthread_t myThread[THREAD_NUM];
	for(i=1;i<THREAD_NUM;i++)
	{
		void **tt=(void **)calloc(2,sizeof(void *));
		tt[0]=(void *)i;
		tt[1]=(void *)THREAD_NUM;
		int check=pthread_create( &myThread[i-1], NULL, multiGenerateSplit, (void*)tt);
		if(check)
    	{
        	fprintf(stderr,"threadNum:%lu, Error - pthread_create() return code: %d\n",i,check);
        	exit(EXIT_FAILURE);
    	}
	}
	for(i=1;i<THREAD_NUM;i++)
	{
		pthread_join( myThread[i-1], NULL);
	}
	//////////////////////////////////multi-thread generate spcode////////////////////////////////////
	rwlockRed=(pthread_rwlock_t*)calloc(redCapacity,sizeof(pthread_rwlock_t));
	for(i=0;i<redCapacity;i++)
	{
		if(pthread_rwlock_init(&rwlockRed[i], NULL))
		{
			printf("fail to create rwlock %lu\n",i);
			exit(1);
		}	
	}
	for(i=0;i<THREAD_NUM;i++)
	{
		void **tt=(void **)calloc(2,sizeof(void *));
		tt[0]=(void *)i;
		tt[1]=(void *)bin;
		int check=pthread_create( &myThread[i], NULL, multiGenerateSP, (void*)tt);
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
	spSplitCapacity[0]=spSplit[0];
	for(i=1;i<THREAD_NUM;i++)
	{
		spSplitCapacity[i]=spSplit[i];
		spSplit[i]=spSplit[i]+spSplit[i-1];
	}
	free(rwlockRed);
	free(reference);
	free(specialSA);
	free(specialBranch);
	free(redSeq);
	free(redPoint);
	free(blackTable);
	/////////////////////////////////add seg to the blue table and special sp index////////////////////
	for(i=1;i<THREAD_NUM;i++)
	{
		void **tt=(void **)calloc(2,sizeof(void *));
		tt[0]=(void *)i;
		tt[1]=(void *)bin;
		int check=pthread_create( &myThread[i], NULL, multiAddSeg, (void*)tt);
		if(check)
    	{
        	fprintf(stderr,"threadNum:%lu, Error - pthread_create() return code: %d\n",i,check);
        	exit(EXIT_FAILURE);
    	}	
    	//pthread_join( myThread[i], NULL);
	}
	for(i=1;i<THREAD_NUM;i++)
	{
		pthread_join( myThread[i], NULL);
	}
	char *spSpecialIndexPath=getPath(bin,"/spSpecialIndex");
	FILE *fpspSpecial=fopen(spSpecialIndexPath,"wb");
	free(spSpecialIndexPath);
	spSpecialIndexPath=getPath(bin,"/spSpecialIndex0");
	FILE *fpspSpecial0=fopen(spSpecialIndexPath,"rb");
	uint64_t *readBuf=(uint64_t *)calloc(BUFFERSIZE,sizeof(uint64_t)),bufReadNum;
	while((bufReadNum=fread(readBuf,sizeof(uint64_t),BUFFERSIZE,fpspSpecial0))>0)
	{
		fwrite(readBuf,sizeof(uint64_t),bufReadNum,fpspSpecial);
	}
	remove(spSpecialIndexPath);
	free(spSpecialIndexPath);
	fclose(fpspSpecial0);
	for(i=1;i<THREAD_NUM;i++)
	{
		char cNum[4];
		sprintf(cNum,"%lu",i);
		char spSpecialName[30]="/spSpecialIndexAdded";
		strcat(spSpecialName,cNum);
		spSpecialIndexPath=getPath(bin,spSpecialName);
		fpspSpecial0=fopen(spSpecialIndexPath,"rb");
		while((bufReadNum=fread(readBuf,sizeof(uint64_t),BUFFERSIZE,fpspSpecial0))>0)
		{
			fwrite(readBuf,sizeof(uint64_t),bufReadNum,fpspSpecial);
		}
		remove(spSpecialIndexPath);
		free(spSpecialIndexPath);
		fclose(fpspSpecial0);
	}
	fclose(fpspSpecial);
	free(readBuf);
       	//////////////////////////connect the sp codes///////////////////////////////////////////////
	spCodeLen=spSplit[THREAD_NUM-1];
	spCodeLen+=32;
    uint64_t spCodeSpace=spCodeLen>>5;
    uint64_t spLenMod=spCodeLen&MOD32;
    if(spLenMod) spCodeSpace++;
    spCode=(uint64_t *)calloc(spCodeSpace,sizeof(uint64_t));
	spSpecialIndexPath=getPath(bin,"/SPcode0");
	fpspSpecial=fopen(spSpecialIndexPath,"rb");//just used as temp paramater.
	uint64_t capacity0=spSplitCapacity[0];
	uint64_t space0=capacity0>>5;
	if((capacity0&MOD32)!=0) space0++;
	bufReadNum=fread(spCode,sizeof(uint64_t),space0,fpspSpecial);
	remove(spSpecialIndexPath);
	free(spSpecialIndexPath);
	if(bufReadNum!=space0) fprintf(stderr,"alert! bufReadNum=%lu, space0=%lu\n",bufReadNum,space0 ),exit(1);
	fclose(fpspSpecial);
	fusionNode[0]=spCode[space0-1];
	fusionMod[0]=capacity0&MOD32;
	for(i=1;i<THREAD_NUM;i++)
	{
		void **tt=(void **)calloc(2,sizeof(void *));
		tt[0]=(void *)i;
		tt[1]=(void *)bin;
		int check=pthread_create( &myThread[i], NULL, multiCatSP, (void*)tt);
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
	for(i=1;i<THREAD_NUM;i++)
	{
		void **tt=(void **)calloc(2,sizeof(void *));
		tt[0]=(void *)i;
		tt[1]=(void *)bin;
		int check=pthread_create( &myThread[i], NULL, multiConnect, (void*)tt);
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
	uint64_t space=spSplit[THREAD_NUM-1]>>5;
	spCode[space]=fusionNode[THREAD_NUM-1];	
	free(spSplit);
	free(spSplitCapacity);
	free(fusionNode);
	free(fusionMod);
	return 1;
}
void *multiGenerateSplit(void *arg)
{
	void **argT=(void **)arg;
	uint64_t num=(uint64_t )argT[0];
	uint64_t THREAD_NUM=(uint64_t)argT[1];
	free(argT);
	uint64_t segment=BWTLEN/THREAD_NUM;
	uint64_t start=num*segment,i;
	uint64_t specialSApoint,branchPoint;
	specialSApoint=BinarySearch(start,specialSA,countRead-1);
	branchPoint=BinarySearch(start,specialBranch,specialBranchNum-1);
	for(i=start;i<BWTLEN;i++)
	{
		unsigned int multioutFlag=0;
		//check whether is multiout or not
		uint64_t distance=specialSA[specialSApoint]-i;
		if(distance>=KMER_LENGTH)// use the main module, for kmer info was stored in it
		{
			uint64_t checkSeq=convert(i)>>((32-KMER_LENGTH)<<1);
			unsigned int formerSeq=checkSeq>>(REDLEN<<1);
			uint64_t latterSeq=(checkSeq&REDEXTRACT);
			uint64_t redLen;
			//operation to check multiout, use hash
			if(formerSeq>0)
			{
				redLen=blackTable[formerSeq]-blackTable[formerSeq-1];
			}
			else
			{
				redLen=blackTable[0];
			}
			if(redLen>0)
			{
				//binary search
				uint64_t startPoint=blackTable[formerSeq]-redLen;
				uint64_t *searchRedSeq=&redSeq[startPoint];
				uint64_t candidate=BinarySearch_red(latterSeq,searchRedSeq,redLen-1);
				if(candidate<redLen)
				{
					if((searchRedSeq[candidate]>>2)==latterSeq)
					{
						// Then we can judge whether it is multiout
						unsigned int tempFlag=searchRedSeq[candidate]&3;
						multioutFlag=tempFlag&1;
					}
				}
			}
		}
		else// use the special module
		{
			if(distance==0)
			{
				specialSApoint++;
			}
			//operation to check multiout, use specialBranch
			if(branchPoint<specialBranchNum&&specialBranch[branchPoint]==i)
			{
				multioutFlag=1;
				branchPoint++;
			}
		}
		if(multioutFlag)
		{
			splitIndex[num]=i+1;
			break;
		}
	}
	return (void *)NULL;
}
void *multiCatSP(void *arg)
{
	void **argT=(void **)arg;
	uint64_t num=(uint64_t )(argT[0]);
	char *bin=(char *)(argT[1]);
	free(argT);
	char cNum[4];
	sprintf(cNum,"%lu",num);	
	uint64_t startPoint=spSplit[num-1];
	uint64_t space=spSplitCapacity[num]>>5;
	space++;
	tempSP[num]=(uint64_t *)calloc(space,sizeof(uint64_t));
	char spName[12]="/SPcode";
	strcat(spName,cNum);
	char *spPath=getPath(bin,spName);
	FILE *fpSP=fopen(spPath,"rb");
	free(spPath);
	fread(tempSP[num],sizeof(uint64_t),space,fpSP);
	fusionMod[num]=spSplit[num]&MOD32;
	fclose(fpSP);
	uint64_t tempMod=spSplitCapacity[num]&MOD32;
	uint64_t tempNode=0;
	if(fusionMod[num]>tempMod)//we move last 2 one to the left and add last one.
	{                         
		unsigned int pledge=fusionMod[num]-tempMod;
		tempNode=tempSP[num][space-2]<<((32-pledge)<<1);
		tempNode=tempNode|(tempSP[num][space-1]>>(pledge<<1));
	}
	else//we move last one to left
	{
		unsigned int pledge=tempMod-fusionMod[num];//Decide whether there are something left
		tempNode=tempSP[num][space-1]<<(pledge<<1);
	}
	fusionNode[num]=tempNode;
	return (void *)NULL;
}

void *multiConnect(void *arg)
{
	void **argT=(void **)arg;
	uint64_t num=(uint64_t )(argT[0]);
	char *bin=(char *)(argT[1]);
	free(argT);
	uint64_t start=spSplit[num-1]>>5, end=spSplit[num]>>5;
	uint64_t i,pointSP=0,recurTemp=fusionNode[num-1];//use the node we prepared before
	unsigned int move=fusionMod[num-1]<<1;
	unsigned int recurMove=(32-fusionMod[num-1])<<1;
	if(move!=0)
	{
		for(i=start;i<end;i++,pointSP++)
		{	
			uint64_t temp=tempSP[num][pointSP];
			uint64_t latter=temp>>move;
			spCode[i]=recurTemp|latter;
			recurTemp=temp<<recurMove;
		}	
	}
	else
	{
		for(i=start;i<end;i++,pointSP++)
		{	
			spCode[i]=tempSP[num][pointSP];
		}
	}
	char spName[12]="/SPcode";
	char cNum[4];
	sprintf(cNum,"%lu",num);
	strcat(spName,cNum);
	char *spPath=getPath(bin,spName);
	remove(spPath);
	free(spPath);
	free(tempSP[num]);
	return (void *)NULL;
}

void *multiAddSeg(void *arg)
{
	void **argT=(void **)arg;
	uint64_t num=(uint64_t )(argT[0]);
	char *bin=(char *)(argT[1]);
	free(argT);
	char cNum[4];
	sprintf(cNum,"%lu",num);
	uint64_t bufferSize=BUFFERSIZE;
	char spSpecialNameAdded[30]="/spSpecialIndexAdded";
	strcat(spSpecialNameAdded,cNum);
	char *spSpecialIndexPath=getPath(bin,spSpecialNameAdded);
	FILE *fpspSpecialAdded=fopen(spSpecialIndexPath,"wb");
	free(spSpecialIndexPath);
	char spSpecialName[30]="/spSpecialIndex";
	strcat(spSpecialName,cNum);
	spSpecialIndexPath=getPath(bin,spSpecialName);
	FILE *fpspSpecial=fopen(spSpecialIndexPath,"rb");
	uint64_t *readBuf=(uint64_t *)calloc(bufferSize,sizeof(uint64_t));
	uint64_t bufReadNum,i;
	while((bufReadNum=fread(readBuf,sizeof(uint64_t),bufferSize,fpspSpecial))>0)
	{
		for(i=0;i<bufReadNum;i++)
		{
			readBuf[i]+=spSplit[num-1];
		}
		fwrite(readBuf,sizeof(uint64_t),bufReadNum,fpspSpecialAdded);
	}
	fclose(fpspSpecial);
	fclose(fpspSpecialAdded);
	remove(spSpecialIndexPath);
	free(spSpecialIndexPath);
	char addBlueSegName[30]="/addBlueSeg";//the indexes has been stored out
	strcat(addBlueSegName,cNum);
	char *addBlueSegPath=getPath(bin,addBlueSegName);
	FILE *fpAddBlueSeg=fopen(addBlueSegPath,"rb");
	while((bufReadNum=fread(readBuf,sizeof(uint64_t),bufferSize,fpAddBlueSeg))>0)
	{
		for(i=0;i<bufReadNum;i++)
		{
			uint64_t temp=blueTable[readBuf[i]];
			uint64_t tempSeg=temp>>4;
			tempSeg+=spSplit[num-1];
			uint64_t tempEliminate=15;
			blueTable[readBuf[i]]=(tempSeg<<4)|(temp&tempEliminate);
		}
	}
	fclose(fpAddBlueSeg);
	free(readBuf);
	remove(addBlueSegPath);
	free(addBlueSegPath);
	return (void *)NULL;
}

void *multiGenerateSP(void *arg)
{
	void **argT=(void **)arg;
	uint64_t num=(uint64_t )(argT[0]);
	char *bin=(char *)(argT[1]);
	free(argT);
	uint64_t bufferSize=BUFFERSIZE;
	uint64_t queueSize=bufferSize;
	uint64_t start=splitIndex[num],end=splitIndex[num+1];
	uint64_t i;
	uint64_t spIndex=0;
	char cNum[4];
	sprintf(cNum,"%lu",num);
	char spName[12]="/SPcode";
	strcat(spName,cNum);
	char *spPath=getPath(bin,spName);
	FILE *fpSP=fopen(spPath,"wb");
	free(spPath);
	char spSpecialName[30]="/spSpecialIndex";
	strcat(spSpecialName,cNum);
	char *spSpecialIndexPath=getPath(bin,spSpecialName);
	FILE *fpspSpecial=fopen(spSpecialIndexPath,"wb");
	free(spSpecialIndexPath);
	char addBlueSegName[30]="/addBlueSeg";//the indexes has been stored out
	strcat(addBlueSegName,cNum);
	char *addBlueSegPath=getPath(bin,addBlueSegName);
	FILE *fpAddBlueSeg=fopen(addBlueSegPath,"wb");
	uint64_t *addBlueSegBuf=NULL,addBlueSegIndex=0; //used to store seg bonus
	if(num!=0)
	{
		addBlueSegBuf=(uint64_t *)calloc(BUFFERSIZE,sizeof(uint64_t));
	}
	uint64_t specialSApoint,branchPoint,spSpecialPoint=0;
	uint64_t spBufPoint=0;
	uint64_t qup=0;
	specialSApoint=BinarySearch(start,specialSA,countRead-1);
	branchPoint=BinarySearch(start,specialBranch,specialBranchNum-1);
	//printf("specialSApoint=%lu,branchPoint=%lu\n",specialSApoint,branchPoint);
	uint64_t *myQueue=(uint64_t *)calloc(bufferSize,sizeof(uint64_t));//stand for the red index
	if(NULL==myQueue)
	{
		printf("fail to alloc myQueue\n");
		exit(1);
	}
	char *bwtQueue=(char *)calloc(bufferSize,sizeof(char));//store in the blue table
	if(NULL==bwtQueue)
	{
		printf("fail to alloc bwtQueue\n");
		exit(1);
	}
	uint64_t *spCodeBuf=(uint64_t *)calloc(bufferSize,sizeof(uint64_t));
	if(NULL==spCodeBuf)
	{
		printf("fail to alloc spCodeBuf\n");
		exit(1);
	}
	uint64_t *spSpecialbuf=(uint64_t *)calloc(bufferSize,sizeof(uint64_t));
	if(NULL==spSpecialbuf)
	{
		printf("fail to alloc spSpecialbuf\n");
		exit(1);
	}
	//printf("start:%lu, end:%lu\n",start,end );
	for(i=start;i<end;i++)
	{
		unsigned int multioutFlag=0;
		unsigned int multiinFlag=0;
		//check whether is multiout or not
		uint64_t distance=specialSA[specialSApoint]-i;
		if(distance>=KMER_LENGTH)// use the main module, for kmer info was stored in it
		{
			uint64_t checkSeq=convert(i)>>((32-KMER_LENGTH)<<1);
			unsigned int formerSeq=checkSeq>>(REDLEN<<1);
			uint64_t latterSeq=(checkSeq&REDEXTRACT);
			uint64_t redLen;
			//operation to check multiout, use hash
			if(formerSeq>0)
			{
				redLen=blackTable[formerSeq]-blackTable[formerSeq-1];
			}
			else
			{
				redLen=blackTable[0];
			}
			if(redLen>0)
			{
				//binary search
				uint64_t startPoint=blackTable[formerSeq]-redLen;
				uint64_t *searchRedSeq=&redSeq[startPoint];
				uint64_t candidate=BinarySearch_red(latterSeq,searchRedSeq,redLen-1);
				if(candidate<redLen)
				{
					if((searchRedSeq[candidate]>>2)==latterSeq)
					{
						// Then we can judge whether it is multiout
						unsigned int tempFlag=searchRedSeq[candidate]&3;
						multioutFlag=tempFlag&1;
						multiinFlag=tempFlag>>1;
						if(multiinFlag)
						{
							//put the index of red table into the queue
							//put the bwt of red table into bwt queue
							if(qup>=queueSize)
							{
								queueSize=(queueSize<<1);
								printf("realloc queue...\n");
								myQueue=realloc(myQueue,queueSize*sizeof(uint64_t));
								//printf("success realloc queue\n");
								bwtQueue=realloc(bwtQueue,queueSize*sizeof(char));
								//printf("success realloc bwtQueue\n");
							}
							myQueue[qup]=startPoint+candidate;
							//record the previous character
							if(specialSApoint>0)
							{
								if(specialSA[specialSApoint-1]+1==i)
								{
									bwtQueue[qup]='#';//stand for #
								}
								else
								{
									bwtQueue[qup]=getCharacter(i-1);//get the precursor of i
								}
							}
							else
							{
								if(0==i)
								{
									bwtQueue[qup]='$';//stand for $
								}
								else
								{
									bwtQueue[qup]=getCharacter(i-1);//get the precursor of i
								}
							}
							qup++;
						}
					}
				}
			}
		}
		else// use the special module
		{
			if(distance==0)
			{
				specialSApoint++;
			}
			//operation to check multiout, use specialBranch
			if(branchPoint<specialBranchNum&&specialBranch[branchPoint]==i)
			{
				multioutFlag=1;
				branchPoint++;
			}
		}
		////////////////////////////end check multiout/////////////////////////////////
		if(multioutFlag)
		{
			//get the SP code's characters
			char tempSP;
			if(distance==KMER_LENGTH)
			{
				//sp code only use 2 bits.
				tempSP='T';
				spSpecialbuf[spSpecialPoint++]=spIndex;
				//printf("# in spCode:%lu\n",spIndex );
				if(spSpecialPoint>=bufferSize)
				{
					//fwrite the spSpecial into the file
					fwrite(spSpecialbuf,sizeof(uint64_t),bufferSize,fpspSpecial);
					spSpecialPoint=0;
				}
			}
			else
			{
				tempSP=getCharacter(i+KMER_LENGTH);
			}
			//printf("%lu: SP:%c\n",spBufPoint,tempSP );
			uint64_t lowBound=spBufPoint>>5;
			uint64_t lowMod=spBufPoint&MOD32;
			unsigned int move=(31-lowMod)<<1;
			spCodeBuf[lowBound]=spCodeBuf[lowBound]|(trans[(int)tempSP]<<move);//store the sp code
			spBufPoint++;
			if((spBufPoint>>5)==bufferSize)
			{
				//fwrite the sp code into the file
				fwrite(spCodeBuf,sizeof(uint64_t),bufferSize,fpSP);
				free(spCodeBuf);
				spCodeBuf=(uint64_t *)calloc(bufferSize,sizeof(uint64_t));
				spBufPoint=0;
			}
			//pop the queue and generate the SP code. Fill the blue table.
			while(qup)
			{
				qup--;
				//pop
				uint64_t tempBlueChar=trans[(int)bwtQueue[qup]];//precursor character store in low position
				uint64_t redIndex=myQueue[qup];
                                uint64_t blueTemp=tempBlueChar|(spIndex<<4);
				pthread_rwlock_wrlock(&rwlockRed[redIndex]);
				blueTable[redPoint[redIndex]]=blueTemp;
				if(num!=0) addBlueSegBuf[addBlueSegIndex++]=redPoint[redIndex];
				redPoint[redIndex]--;
				pthread_rwlock_unlock(&rwlockRed[redIndex]);
				if(addBlueSegIndex==bufferSize)
				{
					fwrite(addBlueSegBuf,sizeof(uint64_t),addBlueSegIndex,fpAddBlueSeg);
					addBlueSegIndex=0;
				}

			}
			spIndex++;
		}
	}
	if(addBlueSegIndex>0)
	{
		fwrite(addBlueSegBuf,sizeof(uint64_t),addBlueSegIndex,fpAddBlueSeg);
	}
	if(spSpecialPoint>0)
	{
		fwrite(spSpecialbuf,sizeof(uint64_t),spSpecialPoint,fpspSpecial);
	}
	if(spBufPoint>0)
	{
		uint64_t lowBound=spBufPoint>>5;
		uint64_t lowMod=spBufPoint&MOD32;
		if(lowMod) lowBound++;
		fwrite(spCodeBuf,sizeof(uint64_t),lowBound,fpSP);
	}
	spSplit[num]=spIndex;
	if(spSpecialbuf!=NULL)	free(spSpecialbuf);
	if(spCodeBuf!=NULL)	free(spCodeBuf);
	if(bwtQueue!=NULL)	free(bwtQueue);
	if(myQueue!=NULL)	free(myQueue);
	if(num==0) remove(addBlueSegPath);
	free(addBlueSegPath);
	free(addBlueSegBuf);
	fclose(fpSP);
	fclose(fpspSpecial);
	fclose(fpAddBlueSeg);
	return (void *)NULL;
}

int ascend(const void *a, const void *b)
{
	uint64_t mka=*(uint64_t *)a, mkb=*(uint64_t *)b;
	return mkb<mka;
}
char getCharacter(uint64_t index)
{
	uint64_t ref_mark=index>>5;
	uint64_t move=(31-(index&MOD32))<<1;
	uint64_t res=(reference[ref_mark]>>move)&3;
	return "ACGT"[res];
}
uint64_t BinarySearch_red(uint64_t mk, uint64_t *target, int64_t up)
{
  int64_t low=0,mid=0;
  while(low<=up)
  {
     mid=(low+up)>>1;
     uint64_t temp=target[mid]>>2;
     if(mk<temp)   up=mid-1;
     else if(mk>temp)   low=mid+1;
     if(mk==temp)    return mid;
  }
  return low;
}
