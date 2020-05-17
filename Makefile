CC=			gcc
#CC=			clang --analyze
CFLAGS=		-g -Wall -Wno-unused-function -O2
WRAP_MALLOC=-DUSE_MALLOC_WRAPPERS
DFLAGS=		-DHAVE_PTHREAD $(WRAP_MALLOC)
AOBJS=		bwamez.o
DEPOBJ=		bwa/bwashm.o bwa/kopen.o
PROG=		zs-mem
INCLUDES=	-Ibwa/
LIBS=		-lm -lz -lpthread

ifeq ($(shell uname -s),Linux)
	LIBS += -lrt
endif

.SUFFIXES:.c .o

.c.o:
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

deps:
		make -C bwa

all:$(PROG)

$(PROG):$(AOBJS) main.o
		$(CC) $(CFLAGS) $(DFLAGS) $(AOBJS) $(DEPOBJ) main.o -o $@ -Lbwa/ -lbwa $(LIBS)

clean-zips:
		rm -f *.o *.out $(PROG)

clean:
		make -C bwa clean
		rm -f *.o *.out $(PROG)

depend:
	( LC_ALL=C ; export LC_ALL; makedepend -Y -- $(CFLAGS) $(DFLAGS) -- *.c )
