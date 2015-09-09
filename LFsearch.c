#include "collect#$.h"
#include "sortBlue.h"
#include "INandOut.h"
#include "insertCase3.h"
#include "LFsearch.h"
#define FASTALEN 70
uint64_t occLen;
uint64_t specialSeg;
char *writeBuf=NULL;
uint64_t writePoint;
FILE *fpRestore=NULL;
uint64_t countBack=0;
int LFsearch(void)
{
	occLen=(bwtLen>>5)+1;
	specialSeg=bwtLen-countRead;
	writeBuf=(char*)calloc(FASTALEN+1,sizeof(char));
	writeBuf[FASTALEN]='\0';
	writePoint=0;
	char *restorePath=getPath("/restore");
	fpRestore=fopen(restorePath,"w");
	free(restorePath);
	uint64_t i;
	for(i=0;i<countRead-1;i++) 
	{
		printf("%lu#: %lu\n",i,specialPos[i] );
	}
	for(i=0;i<4;i++) printf("%lu\n",ACGT[i] );
    printf("countRead=%lu\n",countRead);
	printf("specialSeg: %lu\n",specialSeg );
	printf("$: %lu\n",dollarPos);
	i=bwtLen-1;
	uint64_t I=i>>5;
	uint64_t modI=i&MOD32;
	uint64_t moveI=(31-modI)<<1;
	uint64_t dimer=(bwt[I]>>moveI)&3;
	printf("bwt[%lu]=%c\n",bwtLen-1,"ACGT"[dimer] );
	backSearch(bwtLen-1);
	fclose(fpRestore);
	free(bwt);
	free(specialPos);
	free(occ);
	free(writeBuf);
	return 1;
}
void backSearch(uint64_t num) // stack is too deep! change the method!
{
	while(1&&countBack<31000000000)
	{
		//printf("%lu %lu\n",countBack,writePoint );
		uint64_t dimer;
		uint64_t indexI,modI,moveI;
		indexI=num>>5;
		modI=num&MOD32;
		moveI=(31-modI)<<1;
		dimer=(bwt[indexI]>>moveI)&3;
		uint64_t check;
		char base;
		base="ACGT"[dimer];
		//printf("arrive at %lu\n",num );
		if(dimer==3)//check whether is #$ or not
		{
			if(num==dollarPos)//$
			{
				base='$';
				writeBuf[writePoint++]=base;
				printf("get the $\n");
				break;
			}
			else//#
			{
				check=BinarySearch(num,specialPos,countRead-2);
				if(check<countRead-1&&specialPos[check]==num)
				{
					base='#';
					writeBuf[writePoint++]=base;
					printf("%lu:get an #\n",countBack );
					if(writePoint>=FASTALEN)
					{
						fprintf(fpRestore, "%s\n",writeBuf);
						//putchar('\n');
						writePoint=0;
					}
					uint64_t nextNum=specialSeg+check;
					//printf("nextNum(#): %lu\n",nextNum );
					num=nextNum;
					countBack++;
					continue;
				}
			}
		}
		//printf("%c", base);
		writeBuf[writePoint++]=base;
		if(writePoint>=FASTALEN)
		{
			fprintf(fpRestore, "%s\n",writeBuf);
			//putchar('\n');
			writePoint=0;
		}
		uint64_t nextNum=findSeg(dimer,num);
		//printf("nextNum: %lu\n",nextNum );
		num=nextNum;
		countBack++;
	}
	if(writePoint>0)
	{
		int temp;
		for(temp=0;temp<writePoint;temp++)
		{
			fprintf(fpRestore, "%c", writeBuf[temp]);
		}
		fprintf(fpRestore, "\n" );
	}
	return;
}
uint64_t findSeg(uint64_t type, uint64_t num)//use occ of type and ACGT seg to find next
{
	uint64_t indexI,modI,seg;
	indexI=num>>5;
	modI=num&MOD32;
	uint64_t indexJ,modJ,moveJ,dimer,tempNum;
	uint64_t pSharp;
	tempNum=(indexI<<5);
	if(type==3)
	{
		pSharp=BinarySearch(tempNum,specialPos,countRead-2);
	}
	/*
	if(modI<=15||indexI==(occLen-1))
	{
		seg=occ[indexI][type];
		for(indexJ=indexI,modJ=0;modJ<=modI;modJ++)
		{
			moveJ=(31-modJ)<<1;
			dimer=(bwt[indexJ]>>moveJ)&3;
			if(dimer==type)
			{
				seg++;
			}
		}
	}
	else
	{
		seg=occ[indexI+1][type];
		for(indexJ=indexI+1,modJ=31;modJ>modI;modJ--)
		{
			moveJ=(31-modJ)<<1;
			dimer=(bwt[indexJ]>>moveJ)&3;
			if(dimer==type)
			{
				seg--;
			}
		}
	}
	*/
	seg=occ[indexI][type];
    if(countBack==0) printf("seg=%lu\n",seg);
	for(indexJ=indexI,modJ=0;modJ<=modI;modJ++,tempNum++)
	{
		moveJ=(31-modJ)<<1;
		dimer=(bwt[indexJ]>>moveJ)&3;
		if(type==3)
		{
			if(pSharp<countRead-1&&tempNum==specialPos[pSharp])
			{
				pSharp++;
				continue;	
			}
			else if(tempNum==dollarPos)
			{
				continue;
			}
		}
		if(dimer==type)
		{
			seg++;
            if(countBack==0) printf("seg=%lu\n",seg);
		}
    }
    if(countBack==0) printf("Total G num is %lu\n",seg);
	seg--;//if seg is the first one, then it should be 0+ACGT
	seg+=ACGT[type];
	return seg;
}
