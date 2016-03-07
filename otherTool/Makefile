CC = gcc
CFLAGS = -g -Wall -O2 -Wno-unused-variable -Wno-unused-result -Wno-unused-function
LIB = -lz 

BIN_DIR = .
SRC_DIR = .

SOURCE = $(wildcard ${SRC_DIR}/*.c)
OBJS = $(SOURCE:.c=.o)

transferN : $(OBJS)
	$(CC) $(CFLAGS) -o transferN $(OBJS) $(LIB)

clean :
	rm transferN $(OBJS)

