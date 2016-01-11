#include "collect#$.h"
#include "INandOut.h"
uint64_t case3num=0,blueBoundNum=0;
uint64_t redCapacity=0,mergeCapacity=0,blueCapacity=0;
uint64_t *blackTable=NULL;
uint64_t tempBWTLen=0;
char *specialBwt=NULL;
uint64_t *specialBwtSA=NULL;
//uint64_t check_single[4]={0,0,0,0};
//uint64_t dbgACGT[4]={0,0,0,0};
//uint64_t dbgACGT3[4]={0,0,0,0};
////////////////////////////////////////////////////////////////////////
int generateBlocks(char *bin)
{
	uint64_t i;
	///////////////////////////////////////////get kmer////////////////////////////////////////////
	/*
	kmerInfoPath=getPath("/kmerInfo");
	FILE *fpt=fopen(kmerInfoPath,"rb");
	if(getKmer(fpt)==1) 
	    printf("Success get kmer\n");
	fclose(fpt);
	*/
	///////////////////////////////read the |...#$...|'s SA and bwt/////////////////////////////////////
	char *specialBwtPath=getPath(bin,"/specialBwt");
	FILE *fpSpecialBwt=fopen(specialBwtPath,"rb");
	specialBwt=(char *)calloc(countRead*KMER_LENGTH,sizeof(char));//insert into bwt
  	specialBwtSA=(uint64_t*)calloc(countRead*KMER_LENGTH,sizeof(uint64_t));//insert into bwt
	if(fread(specialBwt,sizeof(char),countRead*KMER_LENGTH,fpSpecialBwt)<countRead*KMER_LENGTH) 
		fprintf(stderr, "failed to read specialBwt!\n" ),exit(1);
	if(fread(specialBwtSA,sizeof(uint64_t),countRead*KMER_LENGTH,fpSpecialBwt)<countRead*KMER_LENGTH)
		fprintf(stderr, "failed to read specialBwtSA!\n" ),exit(1);
	fclose(fpSpecialBwt);
	remove(specialBwtPath);
	free(specialBwtPath);
	///////////////////////////////////////////deal with head#$////////////////////////////////////
	char *headSharpPath=getPath(bin,"/head#");
	char *head$Path=getPath(bin,"/head$");
	FILE *fphSharp,*fphDollar;
	fphSharp=fopen(headSharpPath,"rb");
	fphDollar=fopen(head$Path,"r");
	uint64_t *headSharp=(uint64_t*)calloc(countRead+1,sizeof(uint64_t));
	if(fread(headSharp,sizeof(uint64_t),countRead-1,fphSharp)<countRead-1)
		fprintf(stderr, "failed to read headSharp!\n" ),exit(1);
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
			tempU=tempU|trans[(int)temp[j]]<<move;
		}
		headSharp[i]=tempU;//$|AGCTAGC...|	
		mk=BinarySearch(headSharp[countRead-1],headSharp,countRead-2);
		//printf("mk=%lu\n",mk );
		if(insert(headSharp,mk,headSharp[countRead-1],countRead-1)) ;
		else exit(1);
	}
	fclose(fphSharp);
	fclose(fphDollar);
	remove(headSharpPath);
	remove(head$Path);
	free(headSharpPath);
	free(head$Path);
	/////////////////////////////////////////get tailSharp////////////////////////////////////////////
	char *tailSharpPath=getPath(bin,"/tail#");
	FILE *fptSharp=fopen(tailSharpPath,"rb");//The anciliary module whose last character is # or &
	uint64_t *tailSharp=(uint64_t*)calloc(countRead+1,sizeof(uint64_t));//extra 0 as sentry
	if(fread(tailSharp,sizeof(uint64_t),countRead,fptSharp)<countRead)
		fprintf(stderr, "failed to read tailSharp\n" ),exit(1);
	fclose(fptSharp);
	remove(tailSharpPath);
	free(tailSharpPath);
	////////////////////////////////////merge kmer line///////////////////////////////////////////////////
	char *getKmerPath=getPath(bin,"/getKmer");
	FILE *fpkmer=fopen(getKmerPath,"rb");	
	mergeKmer(fpkmer, deleteSame(tailSharp,countRead,bin), deleteSame(specialBwtSA,countRead*KMER_LENGTH,bin), 
		headSharp, tailSharp, bin);
	fclose(fpkmer);
	remove(getKmerPath);
	free(getKmerPath);
	/////////////////////////////////////////start generate tags/////////////////////////////////////////////////
	//exit(1);
	return 1;
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
		fprintf(stderr,"insert out of boundary!\n");
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
uint64_t *deleteSame(uint64_t *origin,uint64_t len,char *bin)
{
	uint64_t bufferSize=BUFFERSIZE<<1;
	uint64_t *kmerBuffer=(uint64_t *)calloc(bufferSize,sizeof(uint64_t));
	char *deleteSamePath=getPath(bin,"/deleteSame");
	FILE *fpDel=fopen(deleteSamePath,"wb");
	uint64_t i,tempKmer,kbIndex=0,capacity=0;
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
				fwrite(kmerBuffer,sizeof(uint64_t),bufferSize,fpDel);
				kbIndex=0;
				capacity+=bufferSize;
			}
		}
	}
	if(kbIndex>0)
	{
		fwrite(kmerBuffer,sizeof(uint64_t),kbIndex,fpDel);
		capacity+=kbIndex;
	}
	free(kmerBuffer);
	fclose(fpDel);
	uint64_t *res=(uint64_t *)calloc(capacity+1,sizeof(uint64_t));
	fpDel=fopen(deleteSamePath,"rb");
	res[0]=capacity;
	fread(&res[1],sizeof(uint64_t),capacity,fpDel);
	fclose(fpDel);
	remove(deleteSamePath);
	free(deleteSamePath);
	return res;
}
int mergeKmer(FILE *fpkmer,uint64_t *tailSharpKmer, uint64_t *specialKmer, 
	uint64_t *headSharp, uint64_t *tailSharp, char *bin)
{
	uint64_t bufferSize=BUFFERSIZE<<1;
	uint64_t *tailBuffer=(uint64_t *)calloc(bufferSize,sizeof(uint64_t));
	uint64_t *kmerBuffer=(uint64_t *)calloc(bufferSize,sizeof(uint64_t));
	uint64_t kbIndex=0;
	uint64_t point[3]={0,0,0};
	uint64_t capacityTS=tailSharpKmer[0],capacitySB=specialKmer[0];
	uint64_t CP[3];
	CP[1]=capacityTS,CP[2]=capacitySB;
	uint64_t *seq[3];
	seq[1]=&tailSharpKmer[1],seq[2]=&specialKmer[1];
	unsigned int i;
	uint64_t smallest,markSmallest;
	/////////////////////////////////////////buffers&files&pointers for get tags///////////////////////////
	uint64_t mkt=0,mktSharp=0,mkhSharp=0,mkhA=0,mkhC=0,mkhG=0,mkhT=0;
	uint64_t *buft=(uint64_t *)calloc(bufferSize,sizeof(uint64_t));
	uint64_t *bufhA=(uint64_t *)calloc(bufferSize,sizeof(uint64_t));
	uint64_t *bufhC=(uint64_t *)calloc(bufferSize,sizeof(uint64_t));
	uint64_t *bufhG=(uint64_t *)calloc(bufferSize,sizeof(uint64_t));
	uint64_t *bufhT=(uint64_t *)calloc(bufferSize,sizeof(uint64_t));
	uint64_t eliminate=~(ELIMINATE);
	FILE *fphA,*fphC,*fphG,*fphT;
	char *multiInAPath=getPath(bin,"/multiInA");
	char *multiInCPath=getPath(bin,"/multiInC");
	char *multiInGPath=getPath(bin,"/multiInG");
	char *multiInTPath=getPath(bin,"/multiInT");
	fphA=fopen(multiInAPath,"rb");
	fphC=fopen(multiInCPath,"rb");
	fphG=fopen(multiInGPath,"rb");
	fphT=fopen(multiInTPath,"rb");
	char *kmerInfoPath=getPath(bin,"/kmerInfo");
	FILE *fpTail=fopen(kmerInfoPath,"rb");
	fread(buft,sizeof(uint64_t),bufferSize,fpTail);
	fread(bufhA,sizeof(uint64_t),bufferSize,fphA);
	fread(bufhC,sizeof(uint64_t),bufferSize,fphC);
	fread(bufhG,sizeof(uint64_t),bufferSize,fphG);
	fread(bufhT,sizeof(uint64_t),bufferSize,fphT);
	uint64_t storage=0,bwtIndex=0;          //storage stand for bound level 
	uint64_t *case3buf=(uint64_t *)calloc(bufferSize,sizeof(uint64_t)),case3count=0;
	uint64_t *bwtBuf=(uint64_t *)calloc(bufferSize,sizeof(uint64_t));
	uint64_t *blueBuf=(uint64_t *)calloc(bufferSize,sizeof(uint64_t));
	char *case3boundPath=getPath(bin,"/case3bound");
	char *case2bwtPath=getPath(bin,"/case2bwt");
	FILE *fpcase3bound=fopen(case3boundPath,"wb");
	FILE *fpbwt=fopen(case2bwtPath,"wb");
	free(case3boundPath);
	free(case2bwtPath);
	char *blueBoundPath=getPath(bin,"/blueBound");
	FILE *fpBlueBound=fopen(blueBoundPath,"wb");
	free(blueBoundPath);
	uint64_t specialIndex=0;
	uint64_t blueIndex=0;
	uint64_t lowBound,lowMod;
	uint64_t blueBound=0;
	uint64_t *redPointBuf=(uint64_t *)calloc(bufferSize,sizeof(uint64_t));
	uint64_t redIndex=0;
	uint64_t *redSeqBuf=(uint64_t *)calloc(bufferSize,sizeof(uint64_t));
	char *redPointPath=getPath(bin,"/redPoint");
	FILE *fpRedPoint=fopen(redPointPath,"wb");
	free(redPointPath);
	char *redSeqPath=getPath(bin,"/redSeq");
	FILE *fpRedSeq=fopen(redSeqPath,"wb");
	free(redSeqPath);
	uint64_t blackPoint=1;
	blackTable=(uint64_t *)calloc(BLACKCAPACITY,sizeof(uint64_t));
	while((CP[0]=fread(tailBuffer,sizeof(uint64_t),bufferSize,fpkmer))>0)
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
				/////////////////////////////////get tags/////////////////////////////////////
				for(i=0;i<kbIndex;i++)
				{
					uint64_t transI=kmerBuffer[i],temp=0;
					//decode(transI);
					//if(i+1==kmerCapacity) decode(transI);
					int multioutFlag=0,countOut=-1,multiinFlag=0,countIn=-1;
					uint64_t bwtSingle=0;
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
							fread(buft,sizeof(uint64_t),bufferSize,fpTail);
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
						bwtSingle=0;
						//check_single[bwtSingle]+=bufhA[mkhA+1];
						mkhA+=2;
						if(mkhA>=bufferSize)
						{
							fread(bufhA,sizeof(uint64_t),bufferSize,fphA);
							mkhA=0;
						}
					}
					if(bufhC[mkhC]==transI)
					{
						countIn++;
						bwtSingle=1;
						//check_single[bwtSingle]+=bufhC[mkhC+1];
						mkhC+=2;
						if(mkhC>=bufferSize)
						{
							fread(bufhC,sizeof(uint64_t),bufferSize,fphC);
							mkhC=0;
						}
					}
					if(bufhG[mkhG]==transI)
					{
						countIn++;
						bwtSingle=2;
						//check_single[bwtSingle]+=bufhG[mkhG+1];
						mkhG+=2;
						if(mkhG>=bufferSize)
						{
							fread(bufhG,sizeof(uint64_t),bufferSize,fphG);
							mkhG=0;
						}
					}
					if(bufhT[mkhT]==transI)
					{
						countIn++;
						bwtSingle=3;
						//check_single[bwtSingle]+=bufhT[mkhT+1];
						mkhT+=2;
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
				 		}
						for(j=bwtIndex;j<upIndex;j++)//There exists a bug of buf writing
						{
							uint64_t block=j>>5; //2^5=32 each uint64_t for AGCT
							unsigned int move=(31-(j&MOD32))<<1;//62,60,...,2,0.
							bwtBuf[block]=bwtBuf[block]|(bwtSingle<<move); 
						}
						bwtIndex=upIndex;
					}
					if(multioutFlag==1||multiinFlag==1)
					{
						//record the bound of the blue table in red table
						unsigned int tempFlag=0;
						tempFlag=(0|multioutFlag)|(multiinFlag<<1);
						uint64_t tempTransI=transI>>((32-KMER_LENGTH)<<1);
						uint64_t tempSeq=(tempTransI)&REDEXTRACT;
						redSeqBuf[redIndex]=(tempSeq<<2)|tempFlag;
						redPointBuf[redIndex]=blueBound-1; //This is upper bound
						redIndex++;
						if(redIndex>=bufferSize)
						{
							fwrite(redSeqBuf,sizeof(uint64_t),bufferSize,fpRedSeq);
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
							free(bwtBuf);
							bwtBuf=(uint64_t *)calloc(bufferSize,sizeof(uint64_t));
							tempBWTLen+=bwtIndex;
							bwtIndex=0;
							lowBound=0;
						}   
						lowMod=bwtIndex&MOD32;
						uint64_t move=(31-lowMod)<<1;
						bwtBuf[lowBound]=bwtBuf[lowBound]|(trans[(int)specialBwt[specialIndex]]<<move);
						specialIndex++;
						storage++;
						bwtIndex++;
					}
				}
				//////////////////////////////////////////////////////////////////////////////////
				kbIndex=0;
				mergeCapacity+=bufferSize;
				//printf("mergeCapacity=%lu\n", mergeCapacity);
			}
		}
		point[0]=0;
	}
	uint64_t times=0;
	for(i=1;i<=2;i++) times+=(CP[i]-point[i]);	
	while(times--)
	{
		markSmallest=0;
		smallest=MAXU64;
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
			/////////////////////////////////get tags/////////////////////////////////////
			for(i=0;i<kbIndex;i++)
			{
				uint64_t transI=kmerBuffer[i],temp=0;
				//decode(transI);
				//if(i+1==kmerCapacity) decode(transI);
				int multioutFlag=0,countOut=-1,multiinFlag=0,countIn=-1;
				uint64_t bwtSingle=0;
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
						fread(buft,sizeof(uint64_t),bufferSize,fpTail);
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
					bwtSingle=0;
					//check_single[bwtSingle]+=bufhA[mkhA+1];
					mkhA+=2;
					if(mkhA>=bufferSize)
					{
						fread(bufhA,sizeof(uint64_t),bufferSize,fphA);
						mkhA=0;
					}
				}
				if(bufhC[mkhC]==transI)
				{
					countIn++;
					bwtSingle=1;
					//check_single[bwtSingle]+=bufhC[mkhC+1];
					mkhC+=2;
					if(mkhC>=bufferSize)
					{
						fread(bufhC,sizeof(uint64_t),bufferSize,fphC);
						mkhC=0;
					}
				}
				if(bufhG[mkhG]==transI)
				{
					countIn++;
					bwtSingle=2;
					//check_single[bwtSingle]+=bufhG[mkhG+1];
					mkhG+=2;
					if(mkhG>=bufferSize)
					{
						fread(bufhG,sizeof(uint64_t),bufferSize,fphG);
						mkhG=0;
					}
				}
				if(bufhT[mkhT]==transI)
				{
					countIn++;
					bwtSingle=3;
					//check_single[bwtSingle]+=bufhT[mkhT+1];
					mkhT+=2;
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
			 		}
					for(j=bwtIndex;j<upIndex;j++)//There exists a bug of buf writing
					{
						uint64_t block=j>>5; //2^5=32 each uint64_t for AGCT
						unsigned int move=(31-(j&MOD32))<<1;//62,60,...,2,0.
						bwtBuf[block]=bwtBuf[block]|(bwtSingle<<move); 
					}
					bwtIndex=upIndex;
				}
				if(multioutFlag==1||multiinFlag==1)
				{
					//record the bound of the blue table in red table
					unsigned int tempFlag=0;
					tempFlag=(0|multioutFlag)|(multiinFlag<<1);
					uint64_t tempTransI=transI>>((32-KMER_LENGTH)<<1);
					uint64_t tempSeq=(tempTransI)&REDEXTRACT;
					redSeqBuf[redIndex]=(tempSeq<<2)|tempFlag;
					redPointBuf[redIndex]=blueBound-1; //This is upper bound
					redIndex++;
					if(redIndex>=bufferSize)
					{
						fwrite(redSeqBuf,sizeof(uint64_t),bufferSize,fpRedSeq);
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
    					free(bwtBuf);
						bwtBuf=(uint64_t *)calloc(bufferSize,sizeof(uint64_t));
						tempBWTLen+=bwtIndex;
						bwtIndex=0;
						lowBound=0;
					}   
					lowMod=bwtIndex&MOD32;
					uint64_t move=(31-lowMod)<<1;
					bwtBuf[lowBound]=bwtBuf[lowBound]|(trans[(int)specialBwt[specialIndex]]<<move);
					specialIndex++;
					storage++;
					bwtIndex++;
				}
			}
			/////////////////////////////////////////////////////////////////////////////////
			kbIndex=0;
			mergeCapacity+=bufferSize;
			//printf("mergeCapacity=%lu\n", mergeCapacity);
		}
	}
	if(kbIndex>0) 
	{
		/////////////////////////////////get tags/////////////////////////////////////
		for(i=0;i<kbIndex;i++)
		{
			uint64_t transI=kmerBuffer[i],temp=0;
			//decode(transI);
			//if(i+1==kmerCapacity) decode(transI);
			int multioutFlag=0,countOut=-1,multiinFlag=0,countIn=-1;
			uint64_t bwtSingle=0;
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
					fread(buft,sizeof(uint64_t),bufferSize,fpTail);
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
				bwtSingle=0;
				//check_single[bwtSingle]+=bufhA[mkhA+1];
				mkhA+=2;
				if(mkhA>=bufferSize)
				{
					fread(bufhA,sizeof(uint64_t),bufferSize,fphA);
					mkhA=0;
				}
			}
			if(bufhC[mkhC]==transI)
			{
				countIn++;
				bwtSingle=1;
				//check_single[bwtSingle]+=bufhC[mkhC+1];
				mkhC+=2;
				if(mkhC>=bufferSize)
				{
					fread(bufhC,sizeof(uint64_t),bufferSize,fphC);
					mkhC=0;
				}
			}
			if(bufhG[mkhG]==transI)
			{
				countIn++;
				bwtSingle=2;
				//check_single[bwtSingle]+=bufhG[mkhG+1];
				mkhG+=2;
				if(mkhG>=bufferSize)
				{
					fread(bufhG,sizeof(uint64_t),bufferSize,fphG);
					mkhG=0;
				}
			}
			if(bufhT[mkhT]==transI)
			{
				countIn++;
				bwtSingle=3;
				//check_single[bwtSingle]+=bufhT[mkhT+1];
				mkhT+=2;
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
		 		}
				for(j=bwtIndex;j<upIndex;j++)//There exists a bug of buf writing
				{
					uint64_t block=j>>5; //2^5=32 each uint64_t for AGCT
					unsigned int move=(31-(j&MOD32))<<1;//62,60,...,2,0.
					bwtBuf[block]=bwtBuf[block]|(bwtSingle<<move); 
				}
				bwtIndex=upIndex;
			}
			if(multioutFlag==1||multiinFlag==1)
			{
				//record the bound of the blue table in red table
				unsigned int tempFlag=0;
				tempFlag=(0|multioutFlag)|(multiinFlag<<1);
				uint64_t tempTransI=transI>>((32-KMER_LENGTH)<<1);
				uint64_t tempSeq=(tempTransI)&REDEXTRACT;
				redSeqBuf[redIndex]=(tempSeq<<2)|tempFlag;
				redPointBuf[redIndex]=blueBound-1; //This is upper bound
				redIndex++;
				if(redIndex>=bufferSize)
				{
					fwrite(redSeqBuf,sizeof(uint64_t),bufferSize,fpRedSeq);
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
    				free(bwtBuf);
					bwtBuf=(uint64_t *)calloc(bufferSize,sizeof(uint64_t));
					tempBWTLen+=bwtIndex;
					bwtIndex=0;
					lowBound=0;
				}   
				lowMod=bwtIndex&MOD32;
				uint64_t move=(31-lowMod)<<1;
				bwtBuf[lowBound]=bwtBuf[lowBound]|(trans[(int)specialBwt[specialIndex]]<<move);
				specialIndex++;
				storage++;
				bwtIndex++;
			}
		}
		////////////////////////////////////////////////////////////////////
	}
	mergeCapacity+=kbIndex;
	//printf("mergeCapacity=%lu\n", mergeCapacity);
	free(kmerBuffer);
	free(tailBuffer);
	free(tailSharpKmer);
	free(specialKmer);
	printf("end merge\n");
	///////////////////////////////////////add tags free module/////////////////////////////////////////
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
	}
	if(blueIndex>0)
	{
		fwrite(blueBuf,sizeof(uint64_t),blueIndex,fpBlueBound);
		blueBoundNum+=blueIndex;
	}
	if(redIndex>0)
	{
		fwrite(redSeqBuf,sizeof(uint64_t),redIndex,fpRedSeq);
		fwrite(redPointBuf,sizeof(uint64_t),redIndex,fpRedPoint);
		redCapacity+=redIndex;
	}
	blueCapacity=blueBound;
	printf("blackPoint=%lu\n",blackPoint );
	printf("bwtlen=%lu\n",tempBWTLen+blueCapacity );
    printf("mkhSharp=%lu\n",mkhSharp);
	free(case3buf);
	free(bwtBuf);
	free(blueBuf);
	free(specialBwtSA);
	free(specialBwt);
	free(headSharp);
	free(tailSharp);
	remove(multiInAPath);
	remove(multiInCPath);
	remove(multiInGPath);
	remove(multiInTPath);
	free(multiInAPath);
	free(multiInCPath);
	free(multiInGPath);
	free(multiInTPath);
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
	fclose(fpTail);
	fclose(fphA);
	fclose(fphG);
	fclose(fphC);
	fclose(fphT);
	remove(kmerInfoPath);
	free(kmerInfoPath);
	return 1;
}
