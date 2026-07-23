prefix ?= /usr/local/bin #This can be changed
CC ?= gcc
LIBS ?=  # e.g., -L$PREFIX/lib, or where ever htslib is
LIBBIGWIG ?=
CFLAGS ?= -Wall -g -O3 -pthread

.PHONY: all clean install

.SUFFIXES:.c .o

all: MethylSEA

OBJS = common.o bed.o svg.o overlaps.o extract.o MBias.o mergeContext.o perRead.o
VERSION = 0.6.1

version.h:
	echo '#define VERSION "$(VERSION)"' > $@

.c.o:
	$(CC) -c $(CFLAGS) $(LIBS) -IlibBigWig $< -o $@

libMethylSEA.a: version.h $(OBJS)
	-@rm -f $@
	$(AR) -rcs $@ $(OBJS)

lib: libMethylSEA.a

MethylSEA: libMethylSEA.a version.h $(OBJS)
	$(CC) $(CFLAGS) $(LIBS) -o MethylSEA $(OBJS) main.c libMethylSEA.a $(LIBBIGWIG) -lm -lz -lpthread -lhts -lcurl

test: MethylSEA
	python tests/test.py

clean:
	rm -f *.o MethylSEA libMethylSEA.a

install: MethylSEA
	install MethylSEA $(prefix)
