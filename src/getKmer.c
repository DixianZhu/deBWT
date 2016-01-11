#include "collect#$.h"
void bufferIO (uint64_t bufferSize, unsigned int *headFlag, uint64_t *kmerInfo, 
	uint64_t *headBuffer, FILE *fph, char *bin);
///////////////////////////////////////////////////////////////////////
uint64_t countRound=0; //combine getKmer&getIn
uint64_t *kmerBuffer=NULL;
uint64_t kmerBufferSize=BUFFERSIZE<<1;
uint64_t kbIndex=0,tempKmer;
FILE *fpkmer=NULL;
////////////////////////////////////////////////////////////////////////

void *getKmer(void *str)
{
	char *bin=(char *)str;
	uint64_t bufferSize=BUFFERSIZE<<1;
	char *getKmerPath=getPath(bin, "/getKmer");
	fpkmer=fopen(getKmerPath,"wb");
	if(fpkmer==NULL) printf("failed to open write file\n"),exit(1);
	kmerBuffer=(uint64_t *)calloc(kmerBufferSize,sizeof(uint64_t));
    /////////////////////////////////////////generate In ACGT////////////////////////////////////////////
	uint64_t *headBuffer=(uint64_t*)calloc(bufferSize,sizeof(uint64_t));//create the buffer for head ACGT
	uint64_t *readBuf=(uint64_t*)calloc(bufferSize,sizeof(uint64_t));
	unsigned int headFlag=0;//distinguish the AGCT
	char *multiInAPath=getPath(bin, "/multiInA");
	FILE *fph=fopen(multiInAPath,"wb");
	free(multiInAPath);
	char *kmerInfoPath=getPath(bin, "/kmerInfo");
	FILE *fpKmerInfo=fopen(kmerInfoPath,"rb");
	if(fpKmerInfo==NULL) printf("failed to open kmerInfo!\n"),exit(1);
	uint64_t bufReadNum;
	while((bufReadNum=fread(readBuf,sizeof(uint64_t),bufferSize,fpKmerInfo))>0)
	{
		bufferIO(bufReadNum,&headFlag,readBuf,headBuffer,fph,bin);
	}
	free(readBuf);
	fclose(fpKmerInfo);
	free(kmerInfoPath);
	//printf("T prefix is completed!\n");
	free(headBuffer);	
	fclose(fph);
	if(kbIndex>0)
	{
		fwrite(kmerBuffer,sizeof(uint64_t),kbIndex,fpkmer);
	}
	fclose(fpkmer);
	free(getKmerPath);
	free(kmerBuffer);
	return (void *)1;
}

void bufferIO (uint64_t bufferSize, unsigned int *headFlag, uint64_t *kmerInfo, 
	uint64_t *headBuffer, FILE *fph, char *bin)
{
	uint64_t j;
	uint64_t blockMark=0;//stand for the rest block size
	uint64_t eliminate=~(ELIMINATE);
	if(countRound==0)
	{
		tempKmer=kmerInfo[0]&eliminate;
		kmerBuffer[kbIndex++]=tempKmer;
	}
	for(j=0;j<bufferSize;j+=2)
	{
		uint64_t temp=0;
		temp=kmerInfo[j];
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
				char *multiInCPath=getPath(bin, "/multiInC");
				fph=fopen(multiInCPath,"wb");
				free(multiInCPath);
				//printf("A prefix is completed! trans is %lu\n",(temp>>62));
			}
			else if(*headFlag==1)
			{
				fclose(fph);
				char *multiInGPath=getPath(bin, "/multiInG");
				fph=fopen(multiInGPath,"wb");
				free(multiInGPath);
				//printf("C prefix is completed! trans is %lu\n",(temp>>62));
			}
			else if(*headFlag==2)
			{
				fclose(fph);
				char *multiInTPath=getPath(bin, "/multiInT");
				fph=fopen(multiInTPath,"wb");
				free(multiInTPath);
				//printf("G prefix is completed! trans is %lu\n",(temp>>62));
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
		uint64_t kmer=kmerInfo[j]&eliminate;
		if(kmer>tempKmer)
		{
			tempKmer=kmer;
			kmerBuffer[kbIndex++]=tempKmer;		
			if(kbIndex>=kmerBufferSize)
			{
				fwrite(kmerBuffer,sizeof(uint64_t),kmerBufferSize,fpkmer);
				kbIndex=0;
			}
		}
	}
	countRound++;
	//fwrite(tailBuffer,sizeof(uint64_t),bufferSize,fpt);
	fwrite(headBuffer+blockMark,sizeof(uint64_t),(bufferSize-blockMark),fph);
	return;
}
