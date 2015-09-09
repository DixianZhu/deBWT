uint64_t *blueTable;//point to the SP code's index, can be used instantly this time.
uint64_t spCodeLen;
uint64_t blueBoundNum;
uint64_t blueCapacity;
uint64_t *spCode;
/* This is only for single thread
int cmpSP(const void *a, const void *b, uint64_t offset);
void quickSort(uint64_t *A, int64_t low, int64_t high, uint64_t offset);
int64_t partition(uint64_t *A, uint64_t low, uint64_t high, uint64_t offset);
uint64_t choosePivot(uint64_t *A, uint64_t low, uint64_t high);
*/
void swap(uint64_t *a,uint64_t *b);
uint64_t convertSP(uint64_t mark);//get the 32mer of the spCode
void *multiSplitBlue(void *arg);
void *multiQuickSort(void *arg);

/*
void quickSort(uint64_t *A, int64_t low, int64_t high, uint64_t offset, int num);
int64_t partition(uint64_t *A, uint64_t low, uint64_t high, uint64_t offset, int num);
uint64_t choosePivot(uint64_t *A, uint64_t low, uint64_t high, int num);
*/

int cmpQsort(const void *a, const void *b);
void myQsort(uint64_t *A, int64_t lo, int64_t hi, int num,
    int comp(const void *, const void *, int num));
int cmpSP(const void *a, const void *b, int num);