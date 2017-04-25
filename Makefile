CC=gcc
CFLAGS=-g -Wall -O2 -pthread

HTSDIR=htslib
INCLUDES=-I$(HTSDIR)
LFLAGS=-L$(HTSDIR)
LIBS=$(HTSDIR)/libhts.a
#LDLIBS=-lm -lz -lpthread # older version of htslib
LDLIBS=-lm -lz -llzma -lbz2 -lpthread

MAIN = eagle
AUX = util.o vector.o

all: UTIL HTSLIB READCLASSIFY
	$(CC) $(CFLAGS) $(LFLAGS) $(INCLUDES) $(MAIN).c -o $(MAIN) $(AUX) $(LIBS) $(LDLIBS)

HTSLIB:
	$(MAKE) -C $(HTSDIR)/

UTIL:
	$(CC) $(CFLAGS) -c util.c vector.c $(LDLIBS)

READCLASSIFY:
	$(CC) $(CFLAGS) $(LFLAGS) $(INCLUDES) readclassify.c -o readclassify $(AUX) $(LIBS) $(LDLIBS)

clean:
	rm -f eagle readclassify *.o

# DO NOT DELETE THIS LINE -- make depend needs it
