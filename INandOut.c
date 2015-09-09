#include "collect#$.h"
#include "INandOut.h"
uint64_t case3num=0,blueBoundNum=0;
uint64_t redCapacity=0,mergeCapacity=0,blueCapacity=0;
uint64_t *blackTable=NULL;
uint64_t tempBWTLen=0;
//uint64_t dbgACGT[4]={0,0,0,0};
//uint64_t dbgACGT3[4]={0,0,0,0};
int generateBlocks(void)
{
	char *outPath=getPath("/out");
	char *tailSharpPath=getPath("/specialModule/tail#");
	if(multiOut(outPath,
		tailSharpPath)==1)
	{
		printf("success generate multiOut!\n");
	}
	else
	{
		printf("Failed to generate multiOut!\n");
	}
	free(outPath);
	free(tailSharpPath);
	return 1;
}
int multiOut(char* path1, char* path2)
{
	FILE *fp1=fopen(path1,"r");//The main module
	FILE *fp2=fopen(path2,"rb");//The anciliary module whose last character is # or &
	uint64_t i;
	uint64_t bufferSize=BUFFERSIZE<<1;
	uint64_t *tailSharp=(uint64_t*)calloc(countRead+1,sizeof(uint64_t));//extra 0 as sentry
	//printf("Decode tailSharp:\n");
	/* replaced by binary file
	for(i=0;i<countRead;i++)
	{
		int j,k;
		char temp[KMER_LENGTH];
		fscanf(fp2,"%s",temp);
		uint64_t tempU=0;
		for(j=0,k=31;j<KMER_LENGTH;j++,k--)
		{
			unsigned int move=k<<1;
			tempU=tempU|trans[temp[j]]<<move;
		}
		tailSharp[i]=tempU;//|AGCTAGC...|#
		decode(tailSharp[i]);
	}
	*/
	fread(tailSharp,sizeof(uint64_t),countRead,fp2);
	fclose(fp2);
	uint64_t *tailBuffer=(uint64_t*)calloc(bufferSize,sizeof(uint64_t));//create the buffer for tail 
	uint64_t *headBuffer=(uint64_t*)calloc(bufferSize,sizeof(uint64_t));//create the buffer for head 
	uint64_t *readBuf=(uint64_t*)calloc(bufferSize,sizeof(uint64_t));
	unsigned int headFlag=0;//distinguish the AGCT
	char *multiInAPath=getPath("/multiInA");
	FILE *fph=fopen(multiInAPath,"wb");
	char *multiOutPath=getPath("/multiOut");
	FILE *fpt=fopen(multiOutPath,"wb");
	FILE *fphSharp,*fphDollar;
	char *kmerInfoPath=getPath("/kmerInfo");
	FILE *fpKmerInfo=fopen(kmerInfoPath,"rb");
	free(kmerInfoPath);
	uint64_t bufReadNum;
	////////////////////////////////prepare for merge//////////////////////////////////////
	while((bufReadNum=fread(readBuf,sizeof(uint64_t),bufferSize,fpKmerInfo))>0)
	{
		/* check the order
		for(i=1;i<(bufReadNum>>1);i++)
		{
			uint64_t latter=i<<1,former=(i-1)<<1;
			if(readBuf[latter]<readBuf[former])
			{
				printf("sort error! index=%lu\n",latter);
				decode(readBuf[former]);
				decode(readBuf[latter]);
				exit(1);
			}
		}
		*/
		bufferIO(bufReadNum,&headFlag,readBuf,tailBuffer,headBuffer,fph,fpt);
	}
	free(readBuf);
	fclose(fpKmerInfo);
	printf("T prefix is completed!\n");
	free(tailBuffer);
	free(headBuffer);
	fclose(fpt);
	fclose(fph);
	fpt=fopen(multiOutPath,"rb");
	if(getKmer(fpt)==1) printf("Success get kmer\n");
	fclose(fpt);
	fpt=fopen(multiOutPath,"rb");
	free(multiOutPath);
	char *headSharpPath=getPath("/specialModule/head#");
	char *head$Path=getPath("/specialModule/head$");
	fphSharp=fopen(headSharpPath,"rb");
	fphDollar=fopen(head$Path,"r");
	free(headSharpPath);
	free(head$Path);
	uint64_t *headSharp=(uint64_t*)calloc(countRead+1,sizeof(uint64_t));
	/* replaced by binary file
	for(i=0;i<countRead-1;i++)
	{
		int j,k;
		char temp[KMER_LENGTH];
		fscanf(fphSharp,"%s",temp);
		uint64_t tempU=0;
		for(j=0,k=31;j<KMER_LENGTH;j++,k--)
		{
			unsigned int move=k<<1;
			tempU=tempU|trans[temp[j]]<<move;
		}
		headSharp[i]=tempU;//#|AGCTAGC...|
	}
	*/ 
	fread(headSharp,sizeof(uint64_t),countRead-1,fphSharp);
	uint64_t mk;
	for(i=countRead-1;i<countRead;i++)
	{
		int j,k;
		char temp[KMER_LENGTH];
		fscanf(fphDollar,"%s",temp);
		uint64_t tempU=0;
		for(j=0,k=31;j<KMER_LENGTH;j++,k--)
		{
			unsigned int move=k<<1;
			tempU=tempU|trans[temp[j]]<<move;
		}
		headSharp[i]=tempU;//$|AGCTAGC...|	
		mk=BinarySearch(headSharp[countRead-1],headSharp,countRead-2);
		printf("mk=%lu\n",mk );
		if(insert(headSharp,mk,headSharp[countRead-1],countRead-1)) printf("successfully insert\n");
	}
	fclose(fphSharp);
	fclose(fphDollar);
	printf("Decode headSharp:\n");
	for(i=0;i<countRead;i++)
	{
		decode(headSharp[i]);
	}
	printf("end\n");
	/////////////////////////////////////////do merge/////////////////////////////////////////////
	uint64_t mkt=0,mktSharp=0,mkhSharp=0,
		mkhA=0,mkhC=0,mkhG=0,mkhT=0;
	uint64_t *buft=(uint64_t *)calloc(bufferSize,sizeof(uint64_t));
	uint64_t *bufhA=(uint64_t *)calloc(bufferSize,sizeof(uint64_t));
	uint64_t *bufhC=(uint64_t *)calloc(bufferSize,sizeof(uint64_t));
	uint64_t *bufhG=(uint64_t *)calloc(bufferSize,sizeof(uint64_t));
	uint64_t *bufhT=(uint64_t *)calloc(bufferSize,sizeof(uint64_t));
	uint64_t eliminate=~(ELIMINATE);
	uint64_t extract=ELIMINATE;
	FILE *fphA,*fphC,*fphG,*fphT;
	char *multiInCPath=getPath("/multiInC");
	char *multiInGPath=getPath("/multiInG");
	char *multiInTPath=getPath("/multiInT");
	fphA=fopen(multiInAPath,"rb");
	fphC=fopen(multiInCPath,"rb");
	fphG=fopen(multiInGPath,"rb");
	fphT=fopen(multiInTPath,"rb");
	free(multiInAPath);
	free(multiInCPath);
	free(multiInGPath);
	free(multiInTPath);
	fread(buft,sizeof(uint64_t),bufferSize,fpt);
	fread(bufhA,sizeof(uint64_t),bufferSize,fphA);
	fread(bufhC,sizeof(uint64_t),bufferSize,fphC);
	fread(bufhG,sizeof(uint64_t),bufferSize,fphG);
	fread(bufhT,sizeof(uint64_t),bufferSize,fphT);
	uint64_t storage=0,bwtIndex=0;          //storage stand for bound level 
	uint64_t *case3buf=(uint64_t *)calloc(bufferSize,sizeof(uint64_t)),case3count=0;
	uint64_t *bwtBuf=(uint64_t *)calloc(bufferSize,sizeof(uint64_t));
	uint64_t *blueBuf=(uint64_t *)calloc(bufferSize,sizeof(uint64_t));
	char *case3boundPath=getPath("/case3bound");
	char *case2bwtPath=getPath("/case2bwt");
	FILE *fpcase3bound=fopen(case3boundPath,"wb");
	FILE *fpbwt=fopen(case2bwtPath,"wb");
	free(case3boundPath);
	free(case2bwtPath);
	///////////////////////////////read the special module/////////////////////////////////////
	char *specialBwtPath=getPath("/specialModule/specialBwt");
	FILE *fpSpecialBwt=fopen(specialBwtPath,"rb");
	free(specialBwtPath);
	char *specialBwt=(char *)calloc(countRead*KMER_LENGTH,sizeof(char));//insert into bwt
  	uint64_t *specialBwtSA=(uint64_t*)calloc(countRead*KMER_LENGTH,sizeof(uint64_t));//insert into bwt
	fread(specialBwt,sizeof(char),countRead*KMER_LENGTH,fpSpecialBwt);
	fread(specialBwtSA,sizeof(uint64_t),countRead*KMER_LENGTH,fpSpecialBwt);
	/*
	for(i=0;i<countRead*KMER_LENGTH;i++)
	{
		printf("%c   ", specialBwt[i]);
		//decode(specialBwtSA[i]);
	}
	*/
	fclose(fpSpecialBwt);
	char *getKmerPath=getPath("/getKmer");
	FILE *fpkmer=fopen(getKmerPath,"rb");	
	free(getKmerPath);
	mergeKmer(fpkmer,deleteSame(tailSharp,countRead),deleteSame(specialBwtSA,countRead*KMER_LENGTH));
	fclose(fpkmer);
	char *kmerMergePath=getPath("/kmerMerge");
	fpkmer=fopen(kmerMergePath,"rb");
	free(kmerMergePath);
	//////////////////////////////////////////////////////////////////////////////////////////
	char *blueBoundPath=getPath("/blueBound");
	FILE *fpBlueBound=fopen(blueBoundPath,"wb");
	free(blueBoundPath);
	uint64_t specialIndex=0;
	uint64_t blueIndex=0;
	uint64_t lowBound,lowMod;
	uint64_t blueBound=0;
	uint64_t *redPointBuf=(uint64_t *)calloc(bufferSize,sizeof(uint64_t));
	uint64_t redIndex=0;
	unsigned int *redSeqBuf=(unsigned int *)calloc(bufferSize,sizeof(unsigned int));
	char *redPointPath=getPath("/redPoint");
	FILE *fpRedPoint=fopen(redPointPath,"wb");
	free(redPointPath);
	char *redSeqPath=getPath("/redSeq");
	FILE *fpRedSeq=fopen(redSeqPath,"wb");
	free(redSeqPath);
	uint64_t kmerCapacity=0;
	uint64_t *kmerBuffer=(uint64_t *)calloc(bufferSize,sizeof(uint64_t));
	uint64_t blackPoint=1;
	blackTable=(uint64_t *)calloc(BLACKCAPACITY,sizeof(uint64_t));
	uint64_t dbgSP=0;
	while((kmerCapacity=fread(kmerBuffer,sizeof(uint64_t),bufferSize,fpkmer))>0)
	{
		for(i=0;i<kmerCapacity;i++)
		{
			uint64_t transI=kmerBuffer[i],temp=0;
			//decode(transI);
			//if(i+1==kmerCapacity) decode(transI);
			int multioutFlag=0,countOut=-1,multiinFlag=0,countIn=-1;
			uint64_t bwtSingle;
			if(tailSharp[mktSharp]==transI)//That kmer is multiOut
			{								
				multioutFlag=1;
				temp++;
				mktSharp++;
				while(tailSharp[mktSharp]==transI) mktSharp++,temp++;
			}
			while((buft[mkt]&eliminate)==transI)
			{
				countOut++;
				temp+=buft[mkt+1];
				mkt+=2;
				if(mkt>=bufferSize)
				{
					fread(buft,sizeof(uint64_t),bufferSize,fpt);
					mkt=0;
				}
			}
			if(countOut>0)//That kmer is multiOut
			{
				multioutFlag=1;
			}
			if(headSharp[mkhSharp]==transI)//That kmer is multiIn 
			{
				multiinFlag=1;
				mkhSharp++;
				while(headSharp[mkhSharp]==transI) 
				{
					mkhSharp++;
				}
				//multiIn needn't to take care for store the #and$ by now.
			}
			if(bufhA[mkhA]==transI)
			{
				countIn++;
				mkhA+=2;
				bwtSingle=0;
				if(mkhA>=bufferSize)
				{
					fread(bufhA,sizeof(uint64_t),bufferSize,fphA);
					mkhA=0;
				}
			}
			if(bufhC[mkhC]==transI)
			{
				countIn++;
				mkhC+=2;
				bwtSingle=1;
				if(mkhC>=bufferSize)
				{
					fread(bufhC,sizeof(uint64_t),bufferSize,fphC);
					mkhC=0;
				}
			}
			if(bufhG[mkhG]==transI)
			{
				countIn++;
				mkhG+=2;
				bwtSingle=2;
				if(mkhG>=bufferSize)
				{
					fread(bufhG,sizeof(uint64_t),bufferSize,fphG);
					mkhG=0;
				}
			}
			if(bufhT[mkhT]==transI)
			{
				countIn++;
				mkhT+=2;
				bwtSingle=3;
				if(mkhT>=bufferSize)
				{
					fread(bufhT,sizeof(uint64_t),bufferSize,fphT);
					mkhT=0;
				}
			}
			if(countIn>0)//That kmer is multiIn
			{
				multiinFlag=1;
			}
			uint64_t case3low=storage;
			storage+=temp;    //temp stand for the each kmer's frequences
			uint64_t case3up=storage-1;
			if(multioutFlag==1) dbgSP+=temp;
			if(multiinFlag==1)//record the case3 region should be inserted latter
			{
				case3buf[case3count]=case3low;
				case3buf[case3count+1]=case3up;
				case3count+=2;
				if(case3count>=bufferSize)
				{
				 	fwrite(case3buf,sizeof(uint64_t),bufferSize,fpcase3bound);
					case3count=0;
					case3num+=bufferSize;
				}
				blueBound+=temp;
				blueBuf[blueIndex++]=blueBound-1;
				if(blueIndex>=bufferSize)
				{
				 	fwrite(blueBuf,sizeof(uint64_t),bufferSize,fpBlueBound);//store the up bound
					blueIndex=0;
					blueBoundNum+=bufferSize;
		 		} 
			}	
			else if(temp)//the kmer is single in, then write the bwt out(case2),which use bwtSingle 
			{
		  		uint64_t j,upIndex=bwtIndex+temp; 
				uint64_t upBound=upIndex>>5;
				while(upBound>=bufferSize)//stimulate the renewing buffer
				{
                    uint64_t tempSeg=bufferSize<<5;
                    for(j=bwtIndex;j<tempSeg;j++)
                    {
                        uint64_t block=j>>5;
                        unsigned int move=(31-(j&MOD32))<<1;
                        bwtBuf[block]=bwtBuf[block]|(bwtSingle<<move);
                    }
                    tempBWTLen+=tempSeg;
                    fwrite(bwtBuf,sizeof(uint64_t),bufferSize,fpbwt);
                    free(bwtBuf);
                    bwtBuf=(uint64_t *)calloc(bufferSize,sizeof(uint64_t));
                    upIndex=upIndex-tempSeg;
                    upBound=upIndex>>5;
                    bwtIndex=0;
                    /*
		 			lowBound=bwtIndex>>5;
					lowMod=bwtIndex&MOD32;
					uint64_t tempStore=0;
					if(lowMod)
					{
						tempStore=bwtBuf[lowBound+1];
					}
					fwrite(bwtBuf,sizeof(uint64_t),lowBound,fpbwt);
                    */
                    /*
                    uint64_t ACGTi=0,ACGTSeg=lowBound<<5,ACGTI,ACGTmove,ACGTdimer;
                    for(ACGTi=0;ACGTi<ACGTSeg;ACGTi++)
                    {
                        ACGTI=ACGTi>>5;
                        ACGTmove=(31-(ACGTi&MOD32))<<1;
                        ACGTdimer=(bwtBuf[ACGTI]>>ACGTmove)&3;
                        dbgACGT[ACGTdimer]++;
                    }
                    */
                    /*
					tempBWTLen+=(bwtIndex-lowMod);
					free(bwtBuf);
					bwtBuf=(uint64_t *)calloc(bufferSize,sizeof(uint64_t));
					bwtIndex=lowMod;
					bwtBuf[0]=tempStore;
		 			upIndex=bwtIndex+temp; 
                    */
		 		}
				for(j=bwtIndex;j<upIndex;j++)//There exists a bug of buf writing
				{
					uint64_t block=j>>5; //2^5=32 each uint64_t for AGCT
					unsigned int move=(31-(j&MOD32))<<1;//62,60,...,2,0.
					bwtBuf[block]=bwtBuf[block]|(bwtSingle<<move);
                    //dbgACGT3[bwtSingle]++; 
				}
				bwtIndex=upIndex;
			}
			if(multioutFlag==1||multiinFlag==1)
			{
				//record the bound of the blue table in red table
				unsigned int tempFlag=0;
				tempFlag=(0|multioutFlag)|(multiinFlag<<1);
				uint64_t tempTransI=transI>>((32-KMER_LENGTH)<<1);
				unsigned int tempSeq=(tempTransI)&REDEXTRACT;
				redSeqBuf[redIndex]=(tempSeq<<2)|tempFlag;
				redPointBuf[redIndex]=blueBound-1; //This is upper bound
				redIndex++;
				if(redIndex>=bufferSize)
				{
					fwrite(redSeqBuf,sizeof(unsigned int),bufferSize,fpRedSeq);
					fwrite(redPointBuf,sizeof(uint64_t),bufferSize,fpRedPoint);
					redIndex=0;
					redCapacity+=bufferSize;
				}
				//record the black point 
				unsigned int formerSeq=tempTransI>>(REDLEN<<1);
				blackTable[formerSeq]=blackPoint;	//This is upper bound
				blackPoint++;
			}
			//Do not forget the special module's BWT (AGCT)|(......#......)
			while(specialBwtSA[specialIndex]==transI)//forget to sort specialBwtSA again
			{
				lowBound=bwtIndex>>5;
				//printf("%lu, block:%lu : %c   \n",storage,lowBound,specialBwt[specialIndex]);
				//decode(specialBwtSA[specialIndex]);
				if(lowBound>=bufferSize)
				{
					fwrite(bwtBuf,sizeof(uint64_t),bufferSize,fpbwt);
                    /*
                    uint64_t ACGTi,ACGTSeg=bwtIndex,ACGTI,ACGTmove,ACGTdimer;
                    for(ACGTi=0;ACGTi<ACGTSeg;ACGTi++)
                    {
                        ACGTI=ACGTi>>5 ; 
                        ACGTmove=(31-(ACGTi&MOD32))<<1;
                        ACGTdimer=(bwtBuf[ACGTI]>>ACGTmove)&3;
                        dbgACGT[ACGTdimer]++;
                    }
                    */
					free(bwtBuf);
					bwtBuf=(uint64_t *)calloc(bufferSize,sizeof(uint64_t));
					tempBWTLen+=bwtIndex;
					bwtIndex=0;
					lowBound=0;
				}   
				lowMod=bwtIndex&MOD32;
				uint64_t move=(31-lowMod)<<1;
				bwtBuf[lowBound]=bwtBuf[lowBound]|(trans[specialBwt[specialIndex]]<<move);
                //dbgACGT3[trans[specialBwt[specialIndex]]]++;
				//printf("specialIndex=%lu\n",specialIndex );
				specialIndex++;
				storage++;
				bwtIndex++;
		 	}
		} 
	}
	printf("check storage:%lu\n",storage );
	if(case3count>0)
	{
		fwrite(case3buf,sizeof(uint64_t),case3count,fpcase3bound);
		case3num+=case3count;
	}
	if(bwtIndex>0) 
	{
		tempBWTLen+=bwtIndex;
		lowBound=bwtIndex>>5;
		lowMod=bwtIndex&MOD32;
		if(lowMod) lowBound++;
        /*
        uint64_t ACGTi,ACGTI,ACGTmove,ACGTSeg,ACGTdimer;
        for(ACGTi=0;ACGTi<bwtIndex;ACGTi++)
        {
            ACGTI=ACGTi>>5;
            ACGTmove=(31-(ACGTi&MOD32))<<1;
            ACGTdimer=(bwtBuf[ACGTI]>>ACGTmove)&3;
            dbgACGT[ACGTdimer]++;
        }
        */
		fwrite(bwtBuf,sizeof(uint64_t),lowBound,fpbwt);// store the single in bwt and special's bwt
		for(i=bwtIndex-20;i<bwtIndex;i++)
		{
			uint64_t I=i>>5;
			uint64_t modI=i&MOD32;
			uint64_t moveI=(31-modI)<<1;
			uint64_t dimer=(bwtBuf[I]>>moveI)&3;
			printf("block:%lu :%c\n",I,"ACGT"[dimer] );
		}
	}
	if(blueIndex>0)
	{
		fwrite(blueBuf,sizeof(uint64_t),blueIndex,fpBlueBound);
		blueBoundNum+=blueIndex;
	}
	if(redIndex>0)
	{
		fwrite(redSeqBuf,sizeof(unsigned int),redIndex,fpRedSeq);
		fwrite(redPointBuf,sizeof(uint64_t),redIndex,fpRedPoint);
		redCapacity+=redIndex;
	}
	blueCapacity=blueBound;
	printf("dbgSP=%lu\n",dbgSP );
	printf("blackPoint=%lu\n",blackPoint );
	printf("bwtlen=%lu\n",tempBWTLen+blueCapacity );
	free(case3buf);
	free(bwtBuf);
	free(blueBuf);
	free(specialBwtSA);
	free(specialBwt);
	free(headSharp);
	free(tailSharp);
	free(buft);
	free(bufhA);
	free(bufhC);
	free(bufhG);
	free(bufhT);
	free(redPointBuf);
	free(redSeqBuf);
	fclose(fpRedSeq);
	fclose(fpRedPoint);
	fclose(fpbwt);
	fclose(fpcase3bound);
	fclose(fpBlueBound);
	fclose(fpt);
	fclose(fphA);
	fclose(fphG);
	fclose(fphC);
	fclose(fphT);
	return 1;
}
void bufferIO (uint64_t bufferSize, unsigned int *headFlag, uint64_t *kmerInfo, uint64_t *tailBuffer,
	uint64_t *headBuffer, FILE *fph, FILE *fpt)
{
	uint64_t j;
	uint64_t blockMark=0;//stand for the rest block size
	for(j=0;j<bufferSize;j+=2)
	{
		uint64_t temp=0;
		temp=kmerInfo[j];
		tailBuffer[j]=temp;
		tailBuffer[j+1]=kmerInfo[j+1];
		if((temp>>62)==*headFlag)//take care for low letter
		{
			headBuffer[j]=temp<<2;
			headBuffer[j+1]=kmerInfo[j+1];
		}
		else//write the head buffer out 
		{
			fwrite(headBuffer+blockMark,sizeof(uint64_t),(j-blockMark),fph);
			blockMark=j;
			if(*headFlag==0)
			{
				fclose(fph);
				char *multiInCPath=getPath("/multiInC");
				fph=fopen(multiInCPath,"wb");
				free(multiInCPath);
				printf("A prefix is completed! trans is %lu\n",(temp>>62));
			}
			else if(*headFlag==1)
			{
				fclose(fph);
				char *multiInGPath=getPath("/multiInG");
				fph=fopen(multiInGPath,"wb");
				free(multiInGPath);
				printf("C prefix is completed! trans is %lu\n",(temp>>62));
			}
			else if(*headFlag==2)
			{
				fclose(fph);
				char *multiInTPath=getPath("/multiInT");
				fph=fopen(multiInTPath,"wb");
				free(multiInTPath);
				printf("G prefix is completed! trans is %lu\n",(temp>>62));
			}
			(*headFlag)++;
			if((temp>>62)==*headFlag)
			{
				headBuffer[j]=temp<<2;
				headBuffer[j+1]=kmerInfo[j+1];
			}
			else
			{
				j-=2;
			}
		}
	}
	fwrite(tailBuffer,sizeof(uint64_t),bufferSize,fpt);
	fwrite(headBuffer+blockMark,sizeof(uint64_t),(bufferSize-blockMark),fph);
}
int insert (uint64_t *target, uint64_t mk, uint64_t value, uint64_t len)
{
	if(len>mk)
	{
		uint64_t i;
		for(i=len;i>mk;i--)
		{
			target[i]=target[i-1];
		}
		target[mk]=value;
		return 1;
	}
	else
	{
		printf("mk out of boundary\n");
		return -1;
	}
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
uint64_t *deleteSame(uint64_t *origin,uint64_t len)
{
	uint64_t bufferSize=BUFFERSIZE<<1;
	uint64_t *kmerBuffer=(uint64_t *)calloc(bufferSize,sizeof(uint64_t));
	char *deleteSamePath=getPath("/deleteSame");
	FILE *fpkmer=fopen(deleteSamePath,"wb");
	uint64_t i,tempKmer,kmer,kbIndex=0,capacity=0;
	tempKmer=origin[0];
	kmerBuffer[kbIndex++]=tempKmer;
	for(i=1;i<len;i++)
	{
		if(origin[i]>tempKmer) 
		{
			tempKmer=origin[i];
			kmerBuffer[kbIndex++]=tempKmer;
			if(kbIndex>=bufferSize)
			{
				fwrite(kmerBuffer,sizeof(uint64_t),bufferSize,fpkmer);
				kbIndex=0;
				capacity+=bufferSize;
			}
		}
	}
	if(kbIndex>0)
	{
		fwrite(kmerBuffer,sizeof(uint64_t),kbIndex,fpkmer);
		capacity+=kbIndex;
	}
	free(kmerBuffer);
	fclose(fpkmer);
	uint64_t *res=(uint64_t *)calloc(capacity+1,sizeof(uint64_t));
	fpkmer=fopen(deleteSamePath,"rb");
	free(deleteSamePath);
	res[0]=capacity;
	fread(&res[1],sizeof(uint64_t),capacity,fpkmer);
	fclose(fpkmer);
	return res;
}
int getKmer(FILE *fpt)
{
	uint64_t bufferSize=BUFFERSIZE<<1;
	uint64_t *tailBuffer=(uint64_t *)calloc(bufferSize,sizeof(uint64_t));
	uint64_t *kmerBuffer=(uint64_t *)calloc(bufferSize,sizeof(uint64_t));
	char *getKmerPath=getPath("/getKmer");
	FILE *fpkmer=fopen(getKmerPath,"wb");
	free(getKmerPath);
	uint64_t capacity,kbIndex=0,tempKmer,countRound=0;
	uint64_t i;
	uint64_t eliminate=~(ELIMINATE);
	while((capacity=fread(tailBuffer,sizeof(uint64_t),bufferSize,fpt))>0)
	{
		if(countRound==0)
		{
			tempKmer=tailBuffer[0]&eliminate;
			kmerBuffer[kbIndex++]=tempKmer;
		}
		for(i=0;i<capacity;i+=2)
		{
			uint64_t kmer=tailBuffer[i]&eliminate;
			if(kmer>tempKmer)
			{
				tempKmer=kmer;
				kmerBuffer[kbIndex++]=tempKmer;		
				if(kbIndex>=bufferSize)
				{
					fwrite(kmerBuffer,sizeof(uint64_t),bufferSize,fpkmer);
					kbIndex=0;
				}
			}
		}
		countRound++;
	}
	if(kbIndex>0)
	{
		fwrite(kmerBuffer,sizeof(uint64_t),kbIndex,fpkmer);
	}
	fclose(fpkmer);
	free(tailBuffer);
	free(kmerBuffer);
	return 1;
}
int mergeKmer(FILE *fpt,uint64_t *tailSharp,uint64_t *specialBwtSA)
{
	uint64_t bufferSize=BUFFERSIZE<<1;
	uint64_t *tailBuffer=(uint64_t *)calloc(bufferSize,sizeof(uint64_t));
	uint64_t *kmerBuffer=(uint64_t *)calloc(bufferSize,sizeof(uint64_t));
	char *kmerMergePath=getPath("/kmerMerge");
	FILE *fpkmer=fopen(kmerMergePath,"wb");
	free(kmerMergePath);
	uint64_t kbIndex=0;
	uint64_t point[3]={0,0,0};
	uint64_t capacityTS=tailSharp[0],capacitySB=specialBwtSA[0];
	uint64_t CP[3];
	CP[1]=capacityTS,CP[2]=capacitySB;
	uint64_t *seq[3];
	seq[1]=&tailSharp[1],seq[2]=&specialBwtSA[1];
	unsigned int i;
	uint64_t smallest,markSmallest;
	while((CP[0]=fread(tailBuffer,sizeof(uint64_t),bufferSize,fpt))>0)
	{
		seq[0]=tailBuffer;
		while(point[0]<CP[0])
		{
			markSmallest=0;
			smallest=seq[0][point[0]];
			for(i=1;i<=2;i++)
			{
				if(point[i]<CP[i])
				{
					if(seq[i][point[i]]<smallest)
					{
						markSmallest=i;
						smallest=seq[i][point[i]];
					}
					else if(seq[i][point[i]]==smallest)
					{
						point[i]++;
					}
				}
			}
			kmerBuffer[kbIndex++]=smallest;
			point[markSmallest]++;
			if(kbIndex>=bufferSize)
			{
				fwrite(kmerBuffer,sizeof(uint64_t),bufferSize,fpkmer);
				kbIndex=0;
				mergeCapacity+=bufferSize;
				printf("mergeCapacity=%lu\n", mergeCapacity);
			}
		}
		point[0]=0;
	}
	uint64_t times=0;
	for(i=1;i<=2;i++) times+=CP[i]-point[i];
	uint64_t MAX64=((((uint64_t)1<<63)-1)<<1)|1;
	while(times--)
	{
		markSmallest=0;
		smallest=MAX64;
		for(i=1;i<=2;i++)
		{
			if(point[i]<CP[i])
			{
				if(seq[i][point[i]]<smallest)
				{
					markSmallest=i;
					smallest=seq[i][point[i]];
				}
				else if(seq[i][point[i]]==smallest)
				{
					point[i]++;
				}
			}
		}
		kmerBuffer[kbIndex++]=smallest;
		point[markSmallest]++;
		if(kbIndex>=bufferSize)
		{
			fwrite(kmerBuffer,sizeof(uint64_t),bufferSize,fpkmer);
			kbIndex=0;
			mergeCapacity+=bufferSize;
			printf("mergeCapacity=%lu\n", mergeCapacity);
		}
	}
	if(kbIndex>0) fwrite(kmerBuffer,sizeof(uint64_t),kbIndex,fpkmer);
	mergeCapacity+=kbIndex;
	printf("mergeCapacity=%lu\n", mergeCapacity);
	free(kmerBuffer);
	free(tailBuffer);
	free(tailSharp);
	free(specialBwtSA);
	fclose(fpkmer);
	return 1;
}
