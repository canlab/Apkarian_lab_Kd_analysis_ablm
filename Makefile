# For Linux or any machines with gcc compiler
CC = gcc
CFLAGS = -g -ansi -Wall -pedantic

# For SunOS
#CFLAGS = -Aa

all: ablm

clean:
	/bin/rm *.o ablm

OBJ = allocate_cao.o

ablm: ablm.o $(OBJ) 
	$(CC) $(CFLAGS) -o ablm ablm.o $(OBJ) -lm

