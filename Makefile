CC = gcc
CFLAGS = -g -Wall -O2 -pthread

HTSDIR = htslib-1.3.2
INCLUDES = -I$(HTSDIR)
LFLAGS = -L$(HTSDIR)
LIBS = $(HTSDIR)/libhts.a -lm -lz -lpthread

MAIN = eagle
AUX = util.o vector.o

all: UTIL HTSLIB READCLASSIFY
	$(CC) $(CFLAGS) $(LFLAGS) $(INCLUDES) $(MAIN).c -o $(MAIN) $(LIBS) $(AUX)

HTSLIB:
	$(MAKE) -C $(HTSDIR)/

UTIL:
	$(CC) $(CFLAGS) $(LDLAGS) -c util.c vector.c

READCLASSIFY:
	$(CC) $(CFLAGS) $(LFLAGS) $(INCLUDES) readclassify.c -o readclassify $(LIBS) $(AUX)

clean:
	rm -f eagle readclassify *.o

# DO NOT DELETE THIS LINE -- make depend needs it
