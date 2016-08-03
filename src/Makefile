# Makefile created by Anna Lobley (a.lobley@cs.ucl.ac.uk) Jan 2009
# for compiling pGenTHREADER and mGenTHREADER

CC=gcc
MAKE=make
CFLAGS= -O3
LIBS= -lm

all:	svm_prob pseudo_bas chkparse mtx2con

svm_prob:
	$(CC) $(CFLAGS) -o svm_prob svm_prob.c -lm
	$(CC) $(CFLAGS) -o svm_prob_dom svm_prob_dom.c -lm

pseudo_bas:
	$(CC) $(CFLAGS) -o pseudo_bas pseudo_bas.c potential.c -lm
	$(CC) $(CFLAGS) -o pseudo_bas_dom pseudo_bas_dom.c potential.c -lm

chkparse:
	$(CC) $(CFLAGS) chkparse.c $(LIBS) -o chkparse

mtx2con:
	$(CC) $(CFLAGS) mtx2con.c $(LIBS) -o mtx2con

clean:
	-rm -f ../bin/svm_prob ../bin/pseudo_bas ../bin/pseudo_bas_dom ../bin/svm_prob_dom ../bin/chkparse ../bin/mtx2con
	-rm -f svm_prob pseudo_bas pseudo_bas_dom svm_prob_dom chkparse mtx2con

install:
	-cp mtx2con chkparse svm_prob pseudo_bas pseudo_bas_dom svm_prob_dom ../bin/
