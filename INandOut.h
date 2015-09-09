int multiOut(char* path1, char* path2);
void bufferIO (uint64_t bufferSize, unsigned int *headFlag, uint64_t  *kmerInfo, uint64_t *tailBuffer,
	uint64_t *headBuffer, FILE *fph, FILE *fpt);
uint64_t BinarySearch(uint64_t mk, uint64_t *target, int64_t up);//Search the # and $ point
void decode(uint64_t obj);
int getKmer(FILE *fpt);
int mergeKmer(FILE *fpkmer,uint64_t *tailSharp,uint64_t *specialBwtSA);
uint64_t *deleteSame(uint64_t *origin,uint64_t len);
int insert (uint64_t *target, uint64_t mk, uint64_t value, uint64_t len);
uint64_t countRead;
