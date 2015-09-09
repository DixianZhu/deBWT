#include "collect#$.h"
#include "sortBlue.h"
#include "INandOut.h"
#include "insertCase3.h"
uint64_t *case3bound=NULL;//freed
uint64_t *case2bwt=NULL;//freed
uint64_t *bwt=NULL;
uint64_t bwtLen;
uint64_t *specialPos=NULL;
uint64_t dollarPos;
uint64_t (*occ)[4]=NULL;
uint64_t dbgACGT2[4]={0,0,0,0};
int insertCase3(void)
{
	//read bound and case2bwt into the RAM
	char *case3boundPath=getPath("/case3bound");
	FILE *fpCase3Bound=fopen(case3boundPath,"rb");
	free(case3boundPath);
	case3bound=(uint64_t *)calloc(case3num,sizeof(uint64_t));
	fread(case3bound,sizeof(uint64_t),case3num,fpCase3Bound);
	fclose(fpCase3Bound);
	char *case2bwtPath=getPath("/case2bwt");
	FILE *fpbwt=fopen(case2bwtPath,"rb");
	free(case2bwtPath);
	uint64_t case2len=tempBWTLen;
	uint64_t case2space=case2len>>5; 
	if(case2len&MOD32) case2space++;
	case2bwt=(uint64_t *)calloc(case2space,sizeof(uint64_t));
	fread(case2bwt,sizeof(uint64_t),case2space,fpbwt);
	fclose(fpbwt);
	//end
	//merge the bwt file
	bwtLen=case2len+blueCapacity;
	printf("bwtLen is %lu\n",bwtLen );
	uint64_t bwtSpace=bwtLen>>5;
	if(bwtLen&MOD32) bwtSpace++;
	bwt=(uint64_t *)calloc(bwtSpace,sizeof(uint64_t));
	uint64_t i,segBlue=0,j,segCase2=0;
	printf("let's check case2bwt:\n");
	for(i=tempBWTLen-11;i<tempBWTLen;i++)
	{
		uint64_t I=i>>5;
		uint64_t modI=i&MOD32;
		uint64_t moveI=(31-modI)<<1;
		uint64_t dimer=(case2bwt[I]>>moveI)&3;
		printf("%c\n","ACGT"[dimer] );
	}
	uint64_t modJ,indexJ,moveJ,dimer;//target bwt string
	specialPos=(uint64_t *)calloc(countRead-1,sizeof(uint64_t));
	uint64_t specialPoint=0;
	for(i=0;i<case3num;i+=2) //fill it block by block
	{
		//printf("No. %lu [%lu,%lu]\n",i>>1,case3bound[i],case3bound[i+1]);
		if(i==0)
		{
			if(0<case3bound[i])
			{
				//{0,case3bound[i]-1} inserted the case2
				//printf("segCase2=%lu\n", segCase2);
				segCase2=insertCase2(0,case3bound[i]-1,segCase2);
				//printf("segCase2=%lu\n", segCase2);
			}
		}
		else
		{
			//{case3bound[i-1]+1, case3bound[i]-1} inserted the case2
			segCase2=insertCase2(case3bound[i-1]+1,case3bound[i]-1,segCase2);
			//printf("segCase2=%lu\n", segCase2);
		}
		for(j=case3bound[i];j<=case3bound[i+1];j++,segBlue++)
		{
			// copy the bwt of blue Table using segBlue
			indexJ=j>>5;
			modJ=j&MOD32;
			moveJ=(31-modJ)<<1;
			dimer=blueTable[segBlue]&15;
            //if(dimer<4) dbgACGT2[dimer]++;
			// remember # and $
			if(dimer==4)//record all the bwt position of #
			{
				specialPos[specialPoint++]=j;
				printf("#:%lu\n",j);
				dimer=3;
			}
			if(dimer==5)//record the bwt position of $
			{
				dollarPos=j;
				printf("$:%lu\n",j);
				dimer=3;
			}
			bwt[indexJ]=bwt[indexJ]|(dimer<<moveJ);
		}
		//printf("bwtIndex=%lu\n" ,j );
	}
	if(case3bound[case3num-1]<bwtLen-1)
	{
		//{case3bound[case3num-1]+1,bwtLen-1} inserted the case2
		segCase2=insertCase2(case3bound[case3num-1]+1,bwtLen-1,segCase2);
	}
    /*
    printf("After sort blue, let's check ACGT:(case2bwt verify)\n");
    printf("A:%lu, C:%lu, G:%lu, T:%lu\n",dbgACGT[0],dbgACGT[1],dbgACGT[2],dbgACGT[3]);
    printf("dbgACGT2(check whether insert error):\n");
    printf("A:%lu, C:%lu, G:%lu, T:%lu\n",dbgACGT2[0],dbgACGT2[1],dbgACGT2[2],dbgACGT2[3]);
    */
	printf("segCase2=%lu\n",segCase2 );
	printf("segBlue=%lu\n",segBlue );
    printf("blueCapacity=%lu\n",blueCapacity);
	free(case3bound);
	free(case2bwt);
	free(blueTable);
	// get the occ
	uint64_t occLen=(bwtLen>>5)+1;
	occ=(uint64_t (*)[4])calloc(4*occLen,sizeof(uint64_t));
	uint64_t tempOcc[4]={0,0,0,0};
	uint64_t modI,indexI,moveI;//target bwt string
	uint64_t check=0;
	printf("bwt:\n");
	for(i=0;i<bwtLen;i++)
	{
		indexI=i>>5;
		modI=i&MOD32;
		moveI=(31-modI)<<1;	
		dimer=(bwt[indexI]>>moveI)&3;
		if(dimer==3)
		{
			//check=BinarySearch(i,specialPos,countRead-2);
			if(check<countRead-1&&specialPos[check]==i)
			{
			 	if(modI==31)
				{	
					//printf("smapling %lu\n", indexI);
					int j;
					for(j=0;j<4;j++) occ[indexI+1][j]=tempOcc[j];// 0 ... 30 31(sample) 0 ... 30 31(sample)
					//printf("I:%lu, A:%lu, C:%lu, G:%lu, T:%lu\n",
					//	(indexI<<5)+31,tempOcc[0],tempOcc[1],tempOcc[2],tempOcc[3] );
				 	// >>5 map the sample point ahead of this point.
				}
                check++;
				continue;// it is a #
			}
			else if(i==dollarPos)
			{
				if(modI==31)
				{	
					//printf("smapling %lu\n", indexI);
					int j;
					for(j=0;j<4;j++) occ[indexI+1][j]=tempOcc[j];// 0 ... 30 31(sample) 0 ... 30 31(sample)
					//printf("I:%lu, A:%lu, C:%lu, G:%lu, T:%lu\n",
					//	(indexI<<5)+31,tempOcc[0],tempOcc[1],tempOcc[2],tempOcc[3] );
					// >>5 map the sample point ahead of this point.
				}
				continue;// it is a $
			} 
		}
		tempOcc[dimer]++;
		if(modI==31)
		{
			//printf("smapling %lu\n", indexI);
			int j;
			for(j=0;j<4;j++) occ[indexI+1][j]=tempOcc[j];// 0 ... 30 31(sample) 0 ... 30 31(sample)
			//printf("I:%lu, A:%lu, C:%lu, G:%lu, T:%lu\n",
			//	(indexI<<5)+31,tempOcc[0],tempOcc[1],tempOcc[2],tempOcc[3] );
			// >>5 map the sample point ahead of this point.
		}
	}
    printf("I:%lu, A:%lu, C:%lu, G:%lu, T:%lu\n",
      (indexI<<5)+31,tempOcc[0],tempOcc[1],tempOcc[2],tempOcc[3] );
	// write out the bwt information
	
	for(i=(indexI<<5)+32;i<bwtLen;i++)
	{
		indexI=i>>5;
    	modI=i&MOD32;
		moveI=(31-modI)<<1;
		dimer=(bwt[indexI]>>moveI)&3;
	 	printf("%lu: %c \n",i,"ACGT"[dimer]);
	}
	char *bwtPath=getPath("/bwt");
	fpbwt=fopen(bwtPath,"wb");
	free(bwtPath);
	fwrite(bwt,sizeof(uint64_t),bwtSpace,fpbwt);
	fclose(fpbwt);
	return 1;
}
uint64_t insertCase2(uint64_t low, uint64_t up, uint64_t segCase2)
{
	uint64_t i,j;
	uint64_t modI,indexI,moveI;//target bwt string
	uint64_t modJ,indexJ,moveJ;//source case2bwt string
	for(i=low,j=segCase2;i<=up;i++,j++)
	{
		indexI =i>>5;
		modI=i&MOD32;
		moveI=(31-modI)<<1;
		indexJ=j>>5;
		modJ=j&MOD32;
		moveJ=(31-modJ)<<1;
		uint64_t dimer=(case2bwt[indexJ]>>moveJ)&3;
        //dbgACGT2[dimer]++;
		bwt[indexI]=bwt[indexI]|(dimer<<moveI);
        //dbgACGT[dimer]++;
	}
	//printf("bwtIndex=%lu\n",i );
	return j;  
}
