uint64_t *bwt;
uint64_t bwtLen;
uint64_t *specialPos;
uint64_t dollarPos;
uint64_t (*occ)[4];
uint64_t ACGT[4];
uint64_t findSeg(uint64_t type, uint64_t num);//use occ of type and ACGT seg to find next
void backSearch(uint64_t num);
uint64_t *specialHash;
uint64_t *invHash;
//void reportError(int index);
