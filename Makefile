CC = gcc
CFLAGS = -g -Wall -O2 -Wno-unused-variable -Wno-unused-result -Wno-unused-function
LIB = -lz -lpthread

BIN_DIR = .
SRC_DIR = ./src

SOURCE = $(wildcard ${SRC_DIR}/*.c) 

OBJS = $(SOURCE:.c=.o)

deBWT : $(OBJS) 
	$(CC) $(CFLAGS) -o deBWT $(OBJS) $(LIB)

clean : 
		rm deBWT $(OBJECTS) 