#!/usr/bin/make -f

SHELL=/bin/sh
CC=g++
C=gcc

CFLAGS=-O -Wno-unused-result


# If OpenMP is NOT installed, then uncomment the following line
LDFLAGS= 

# If OpenMP is installed, then uncomment the following line
# LDFLAGS=-fopenmp

EXEC = isosel fastdist

# Set the install directory

DESTDIR=/usr/local/bin/

all: $(EXEC)

isosel: main.cpp fonction.cpp fonction.h memory.h memory.cpp
	$(CC) $(CFLAGS) main.cpp fonction.cpp memory.cpp -o ../bin/isosel $(LDFLAGS)

fastdist: fastdist.cpp fastdist.h
	$(CC) $(CFLAGS) -I eigen fastdist.cpp -o ../bin/fastdist $(LDFLAGS)

# bioNJ: BioNJ.c
# 	$(C) $(CFLAGS) BioNJ.c -o ../bin/bioNJ

install:
	cp -f ../bin/isosel $(DESTDIR) 2> /dev/null
	cp -f ../bin/clustalo $(DESTDIR) 2> /dev/null
	cp -f ../bin/fastdist $(DESTDIR) 2> /dev/null
	cp -f ../bin/mafft $(DESTDIR) 2> /dev/null
	cp -f ../bin/muscle $(DESTDIR) 2> /dev/null
	cp -f ../bin/run_isosel $(DESTDIR) 2> /dev/null
	cp -f ../bin/seaview $(DESTDIR) 2> /dev/null
	cp -R ../bin/mafftdir $(DESTDIR) 2> /dev/null
clean:
	rm -f *.o

