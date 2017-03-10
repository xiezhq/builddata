# SGI C/C++ compiler
#CC = cc

# SGI C++ compiler
#CC = CC

# Intel C/C++ compiler
#CC = icc

# Gnu C/C++ compiler
CC = gcc

# debug flag for SGI r12k cpu
#CFLAGS = -g -fullwarn -n32 -mips4 -r12000 #-apo

# normal optimization for SGI r12k multiprocessor machine
#CFLAGS = -O2 -n32 -mips4 -r12000 #-apo

# debug flag for intel icc compiler
CFLAGS = -g -O2 #-Wall

# normal optimization for intel icc compiler
#CFLAGS = -O2 -parallel -Wall

# debug flag for gcc
#CFLAGS = -g -O2 -Wall


MAKEFILE = makefile

PROGS = main interf asa2pep asa2prot searchpdb

OBJECTS4MAIN	= main.o pdb.o memory.o
OBJECTS4INTERF	= interf.o pdb.o asa2pdb.o open_file.o cnctarea2.o check.o sortag.o calzscore.o memory.o property2res.o
OBJECTS4ASA2PEP	= asa2pep.o pdb.o cnctarea2.o check.o sortag.o asa2pdb.o open_file.o memory.o
OBJECTS4ASA2PROT = asa2prot.o pdb.o cnctarea2.o check.o sortag.o asa2pdb.o open_file.o memory.o

OBJECTS4SEARCHPDB = searchpdb.o open_file.o pdb.o memory.o


all: $(PROGS)

# main
main: $(OBJECTS4MAIN)
	$(CC) $(CFLAGS) -o $@ $(OBJECTS4MAIN)

main.o: pdb.h memory.h $(MAKEFILE)

# interf
interf: $(OBJECTS4INTERF)
	$(CC) $(CFLAGS) -o $@ $(OBJECTS4INTERF) -lm

interf.o: pdb.h memory.h max.h interf.h asa2pdb.h property2res.h $(MAKEFILE)

# asa2pep
asa2pep: $(OBJECTS4ASA2PEP)
	$(CC) $(CFLAGS) -o $@ $(OBJECTS4ASA2PEP) -lm

asa2pep.o: pdb.h asa2pdb.h max.h memory.h $(MAKEFILE)

# asa2prot
asa2prot: $(OBJECTS4ASA2PROT)
	$(CC) $(CFLAGS) -o $@ $(OBJECTS4ASA2PROT) -lm

asa2prot.o: pdb.h asa2pdb.h max.h memory.h $(MAKEFILE)

# searchpdb
searchpdb: $(OBJECTS4SEARCHPDB)
	$(CC) $(CFLAGS) -o $@ $(OBJECTS4SEARCHPDB)

searchpdb.o: searchpdb.h open_file.h pdb.h $(MAKEFILE)


# shared files
pdb.o: pdb.h memory.h max.h property4res.h $(MAKEFILE)
asa2pdb.o: pdb.h cnctarea2.h asa2pdb.h $(MAKEFILE)
cnctarea2.o: $(MAKEFILE)
calzscore.o: memory.h $(MAKEFILE)
property.o: property4res.h $(MAKEFILE)

.c.o:
	$(CC) $(CFLAGS) -c $<


#remove the unnecessary binary files
.PHONY: clean
clean:
	rm -f $(OBJECTS4MAIN) $(OBJECTS4INTERF) $(OBJECTS4ASA2PEP) $(OBJECTS4ASA2PROT) $(OBJECTS4SEARCHPDB) $(PROGS) core a.out
