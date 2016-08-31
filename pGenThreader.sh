#!/bin/bash

#------------------------------------------------------------------------------
# pGenTHREADER
#
# simple shell script for running pGenTHREADER code
#
# Jan 2009 alobley@cs.ucl.ac.uk
#------------------------------------------------------------------------------

#JOB name
JOB=$1
#input file
FSA=$2
#path to psiblast
PSIB=/Users/dbuchan/Code/blast-2.2.26/bin
#path to database
DB=/Users/dbuchan/Code/uniref_test_db/uniref_test.fasta
#path to PSIPRED data files
PDATA=/scratch0/NOT_BACKED_UP/dbuchan/psipred/data/
#path to pGenTHREADER data files
DATA=./data/
#path to PGenThreader binary directory
PGT=./bin
#path to PSIPRED binary directory
PSIP=/Users/dbuchan/Code/psipred/bin
#path to fold library
TDB=/Users/dbuchan/Code/pGenTHREADER/
#TDB=/scratch0/NOT_BACKED_UP/dbuchan/tdb


#make a masked copy of input file
#$PSIP/pfilt -f $FSA > $JOB.fsa
cp %FSA $JOB.fsa

export TDB_DIR=$TDB
export THREAD_DIR=./data

echo started `date` $HOST > $JOB.pgt.log

#Run 3 iterations of PSI-BLAST
$PSIB/blastpgp -a 4 -F T -t 1 -j 3 -v 5000 -b 0 -h 0.001 -i $JOB.fsa -d $DB -C $JOB.chk -F T > /dev/null


 echo "Finished 3 iterations PSI-BLAST"

 echo $JOB.fsa > $JOB.sn
 echo $JOB.chk > $JOB.pn

#make a profile matrix
 $PSIB/makemat -P $JOB

#rename checkpoint file
  mv $JOB.chk $JOB.iter3.chk
  mv $JOB.mtx $JOB.iter3.mtx

# Run PSI-PRED for PGT
$PSIP/psipred $JOB.iter3.mtx $PDATA/weights.dat $PDATA/weights.dat2 $PDATA/weights.dat3 > $JOB.pgen.ss
#$PSIP/psipred $JOB.iter3.mtx $PDATA/weights.dat $PDATA/weights.dat2 $PDATA/weights.dat3 $PDATA/weights.dat4 > $JOB.pgen.ss

$PSIP/psipass2 $PDATA/weights_p2.dat 1 1.0 1.0 $JOB.pgen.ss2 $JOB.pgen.ss > $JOB.horiz


if [ ! -s "$JOB.pgen.ss2" ]
then
   echo "PSIPRED failed ... exiting early"
   echo "PSIPRED failed ... exiting early" >> $JOB.pgt.log
   exit;
fi

echo "Finished PSI-PRED"

#Run further 3 iterations of psi-blast

$PSIB/blastpgp -a 4 -F T -t 1 -i $JOB.fsa -R $JOB.iter3.chk -d $DB -j 3 -v 5000 -b 0 -h 0.001 -C $JOB.chk > /dev/null

echo "Finished 6 iterations of PSI BLAST"

#make PSSM from 6th iteration

 echo $JOB.fsa > $JOB.sn
 echo $JOB.chk > $JOB.pn

#make a profile matrix
 $PSIB/makemat -P $JOB

#rename checkpoint file
  mv $JOB.chk $JOB.iter6.chk
  mv $JOB.mtx $JOB.iter6.mtx

# Run PGenThreader process

$PGT/pseudo_bas -c11.0 -C20 -h0.2 -F$JOB.pgen.ss2 $JOB.iter6.mtx $JOB.pgen.pseudo $DATA/psichain.lst

if [ ! -s "$JOB.pgen.pseudo" ]
then
    echo "pseudo_bas failed"
    echo "pseudo_bas failed" >> $JOB.pgt.log
    exit;
fi

$PGT/svm_prob $JOB.pgen.pseudo | sort -k 2,2rn -k 6,6rn -k 5,5g  > $JOB.pgen.presults

if [ ! -s "$JOB.pgen.presults" ]
then
    echo "sortprob failed"
    echo "sortprob failed" >> $JOB.pgt.log
    exit;
fi

$PGT/pseudo_bas -S -p -c11.0 -C20 -h0.2 -F$JOB.pgen.ss2 $JOB.iter6.mtx $JOB.pgen.pseudo $JOB.pgen.presults > $JOB.pgen.align

echo Finished `date` >> $JOB.pgt.log
