uint64_t BinarySearch(uint64_t mk, uint64_t *target, int64_t up);//Search the # and $ point
void decode(uint64_t obj);
int mergeKmer(FILE *fpkmer,uint64_t *tailSharpKmer, uint64_t *specialKmer, 
	uint64_t *headSharp, uint64_t *tailSharp, char *bin);
uint64_t *deleteSame(uint64_t *origin,uint64_t len,char *bin);
int insert (uint64_t *target, uint64_t mk, uint64_t value, uint64_t len);
uint64_t countRead;
