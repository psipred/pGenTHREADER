#!/bin/bash

#------------------------------------------------------------------------------
# pDomTHREADER
#
# simple shell script for running pDomTHREADER code
#
# Jan 2009 alobley@cs.ucl.ac.uk
#------------------------------------------------------------------------------



#JOB name 
JOB=$1
#input file
FSA=$2
#path to psiblast
PSIB=/usr/local/bin/
#name of database
DB=uniref90
#prefix to PDomThreader directory
HOME=/cs/research/bioinf/software0/green/software/pDomThreader
#path to PSIPRED data files
PDATA=/usr/local/lib/psipred/
#path to pDomTHREADER data files
DATA=./data
#path to pDomThreader binary directory
PDT=./bin
#path to PSIPRED binary directory
PSIP=/usr/local/bin
#path to fold library
TDB=./tdb


#make a masked copy of input file
$PSIP/pfilt -f $FSA > $JOB.fsa

export TDB_DIR=$TDB
export THREAD_DIR=./data

echo started `date` $HOST > $JOB.pdt.log

#Run 3 iterations of PSI-BLAST
 $PSIB/blastpgp -a 4 -F T -i $JOB.fsa -d $DB -j 3 -v 5000 -b 0 -h 0.001 -C $JOB.chk -F T > /dev/null


 echo "Finished 3 iterations PSI-BLAST"

 echo $JOB.fsa > $JOB.sn
 echo $JOB.chk > $JOB.pn

#make a profile matrix
 $PSIB/makemat -P $JOB

#rename checkpoint file
  mv $JOB.chk $JOB.iter3.chk
  mv $JOB.mtx $JOB.iter3.mtx

# Run PSI-PRED for PDT

$PSIP/psipred $JOB.iter3.mtx $PDATA/weights.dat $PDATA/weights.dat2 $PDATA/weights.dat3 > $JOB.pdom.ss

$PSIP/psipass2 $PDATA/weights_p2.dat 1 1.0 1.0 $JOB.pdom.ss2 $JOB.pdom.ss > $JOB.horiz


if [ ! -s "$JOB.pdom.ss2" ] 
then
   echo "PSIPRED failed ... exiting early"
   echo "PSIPRED failed ... exiting early" >> $JOB.pdt.log
   exit;
fi

echo "Finished PSI-PRED"

# Run PDomThreader

$PDT/pseudo_bas_dom -c11.0 -C20 -h0.2 -F$JOB.pdom.ss2 $JOB.iter3.mtx $JOB.pdom.pseudo $DATA/cath.lst > $JOB.pdom.output

if [ ! -s "$JOB.pdom.pseudo" ]
then
    echo "pseudobas failed" 
    echo "pseudobas failed" >> $JOB.pdt.log
    exit;
fi

$PDT/svm_prob_dom $JOB.pdom.pseudo | sort -k 2,2rn -k 6,6rn -k 5,5g  > $JOB.pdom.presults

if [ ! -s "$JOB.pdom.presults" ]
then
    echo "sortprob failed or no hits returned"
    echo "sortprob failed" >> $JOB.pdt.log
    exit;
fi

$PDT/pseudo_bas_dom -S -p -c11.0 -C20 -h0.2 -F$JOB.pdom.ss2 $JOB.iter3.mtx $JOB.pdom.pseudo $JOB.pdom.presults > $JOB.pdom.align

echo Finished `date` >> $JOB.pdt.log
