#################################################################
# WGA-Power Makfile
# Modified from RSEARCH's
# CVS $Id: Makefile.in,v 1.7 2003/01/10 21:33:31 rjklein Exp $
#################################################################
SHELL  = /bin/sh

# On most Unices, you can build the package without modifying anything
#   by just typing "./configure; make".
#
# You may want to modify the following make variables:
#   BINDIR  - where the executables will be installed by a 'make install'
#   MANDIR  - where the man pages will be installed by a 'make install'
#   CC      - which compiler to use
#   CFLAGS  - compiler flags to use

# where you want things installed
# Sort of uses GNU coding standards. ${prefix} might be /usr/local.
# ${exec_prefix} gives you some flexibility for installing architecture
# dependent files (e.g. the programs): an example ${exec_prefix} might be
# /nfs/share/irix64/
#
prefix      = /usr/local
exec_prefix = ${prefix}
BINDIR      = ${exec_prefix}/bin
MANDIR      = ${prefix}/man
DATADIR     = ${prefix}/share

## your compiler and compiler flags
CC      = gcc
#CC = mpicc
#CFLAGS = -O3 -DUSE_MPI 
#CFLAGS = -O3 -Wall 
CFLAGS = -g -Wall

## other defined flags. 
#   contains stuff that autoconf 
#  decides on.  contains stuff that we added to
#  the configure script tests.  -lm -lccmalloc -ldl contains system
#  libraries that the configure script decides we need.
#
MDEFS = 
LIBS = -lm -lRmath -lgsl

# Where my libraries/includes (distribured with program) are
MYLIBS   = -lsquid
MYLIBDIR = -L${prefix}/lib -L/sw/lib
MYINCDIR = -I${prefix}/include -I/sw/include


PROGS = wga-power power-from-lambda-alpha simple-power compute-params

OBJS  = calcpower.o findmarkers.o powerfuncs.o mpifuncs.o
HDRS  = 

.c.o: 
	$(CC) $(CFLAGS) $(MDEFS) $(MYINCDIR) -c $<

#################################################################
## Targets defining how to make RSEARCH executables.
##
all: 	
	make progs

progs:	$(PROGS)

$(PROGS): %: %.o $(OBJS) 
	$(CC) $(CFLAGS) $(MDEFS) $(MYLIBDIR) -o $@ $@.o $(OBJS) $(MYLIBS) $(LIBS)

clean:
	-rm -f *.o *~ Makefile.bak core TAGS gmon.out
