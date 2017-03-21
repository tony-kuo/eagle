CC = gcc
CFLAGS = -g -Wall -O2 -pthread

HTSDIR = htslib-1.3.2
INCLUDES = -I$(HTSDIR)
LFLAGS = -L$(HTSDIR)
LIBS = $(HTSDIR)/libhts.a -lm -lz -lpthread

MAIN = eagle

all: HTSLIB
	$(CC) $(CFLAGS) $(LFLAGS) $(INCLUDES) $(MAIN).c -o $(MAIN) $(LIBS)
	$(CC) $(CFLAGS) $(LFLAGS) $(INCLUDES) readclassify.c -o readclassify $(LIBS)

HTSLIB:
	$(MAKE) -C $(HTSDIR)/

clean:
	rm -f eagle

# DO NOT DELETE THIS LINE -- make depend needs it
