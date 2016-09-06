#!/bin/bash
#------------------------------------------------------------------------------
# pGenTHREADER
#
# new shell script for running pGenTHREADER, pDOmTHREADER code with BLAST+
# and the new MSA options. Now takes command line options!
#
################################################################################
#
# Usage:
#
# Genthreader usage:
# genthreader.sh -i fasta_file -j job_title
#
# DomTHREADER usage
# genthreader.sh -i fasta_file -j job_title -d
#
################################################################################
#
# Changes:
#
# Jan 2009	alobley@cs.ucl.ac.uk
#			pGenTHREADER and pDomTHREADER scripts written
# Feb 2011	dbuchan@cs.ucl.ac.uk
#			Add BLAST+ support and multiple alignment and contact alignments
# Sep 2011	dbuchan@cs.ucl.ac.uk
#			United pGenTHREADER.sh and pDomTHEADER.sh to a single GenThreader.sh
# Jan 2012	dbuchan@cs.ucl.ac.uk
#			Added ligand binding prediction step
#-------------------------------------------------------------------------------

#JOB name
JOB=''
#input file
FSA=''
#TYpe of job to run
TYPE=0
#msa flag
MSA=0
#msa has all close relatives
MSA_COMPLETE=0
#msa has all possible relatives
RELATIVES_COMPLETE=0
#Whether to build C-alpha structural models
BUILD_MODELS=0
#flag to toggle building multiple alignments, requires MAFFT install
BUILD_MULTI=0
#flat to toggle building of contact consensus mapping
BUILD_CONTACTS=0
#Name of the Catalytic site atlas file
CSA=CSA_2_2_12.dat
#path to psiblast/ncbi tools
#PSIB=/usr/local/bin/
PSIB=/scratch0/NOT_BACKED_UP/dbuchan/Applications/ncbi-blast-2.2.31+/bin/
#path to database
DB=/scratch0/NOT_BACKED_UP/dbuchan/uniref/uniref90.fasta
#path to PSIPRED data files
PDATA=/usr/local/lib/psipred/
#PDATA=/scratch0/NOT_BACKED_UP/dbuchan/psipred/data/
#path to pGenTHREADER data files
DATA=./data/
#path to PGenThreader binary directory
PGT=./bin
#path to PSIPRED binary directory
PSIP=/cs/research/bioinf/home1/green/dbuchan/Code/psipred40/bin/
#path to fold library
TDB=/scratch0/NOT_BACKED_UP/dbuchan/GenTHREADER/corrected_tdb/
#location of experimental scripts
BLAST=./BLAST+
#Name of your PSICHAIN file
PSICHAIN=psichain.lst
#Name of your dom chain file
DOMCHAIN=cath_dom.lst
#Toggle whether or not to run the ligand binding prediction at the end
LIGAND_PRED=0
#Directory that contains the svm-predict executable from libsvm
LIBSVM=/scratch0/NOT_BACKED_UP/dbuchan/libsvm-3.11/
ERROR=0;

while getopts i:j:b:u:p:P:t:C:L:ldcmrRhM o
do case "$o" in
	i) FSA="$OPTARG";;
	j) JOB="$OPTARG";;
	c) BUILD_CONTACTS=1;;
	C) CSA="$OPTARG";;
	d) TYPE=1;;
	m) MSA=1;;
	r) MSA_COMPLETE=1;;
	R) RELATIVES_COMPLETE=1;;
	b) PSIB="$OPTARG";;
	u) DB="$OPTARG";;
	p) PDATA="$OPTARG";;
	P) PSIP="$OPTARG";;
	t) TDB="$OPTARG";;
	M) BUILD_MULTI=1;;
	l) LIGAND_PRED=1;;
	L) LIBSVM="$OPTARG";;
	s) BUILD_MODELS=1;;
	h) echo >&2 "Usage: $0 -i file -j 'job_name' [-d] [-m] [-M] [-r] [-R] [-c] [-l] [-s] [-L /libsvm/] [-C CSA_file] [-b /blast_bin_path] [-u /uniref90_path] [-p /psipred_data_path] [-P /psipred_bin_path] [-t /tdb_path]"
		echo "  -i Input fasta file"
		echo "  -j Job name"
		echo "  -d toggle between genthreader and domthreader code (ensure you set -t correctly)"
		echo "  -m flag for input is multiple structural alignment"
		echo "  -r flag for : MSA contains all known close relatives (surpress running initial psi-blast)"
		echo "  -R flag for : MSA additionally contains all known distant relatives (surpress running second psi-blast during pGenTHREADER jobs)"
		echo "  -b location of BLAST+ bin directory"
		echo "  -u location of uniref90 database (prepared with createblastdb)"
		echo "  -p location of PSIPRED data/ directory"
		echo "  -P location of PSIPRED bin/ directory"
		echo "  -t location of tdb files (psichain.lst and cath_dom.lst should be in the genthreader data/ directory)"
		echo "  -M build multiple alignment of CERT, GOOD and MED hits, requires MAFFT and that the mafft_spliced_alignment.pl script is correctly configured "
		echo "  -c build map of residue contacts; requires internet connection. Will not work with -d option"
		echo "  -C specify the name of a different CSA zip file availalble from http://www.ebi.ac.uk/thornton-srv/databases/CSA/archive/"
		echo "  -l flag for running ligand binding prediction"
		echo "  -L specify directory location of libsvm's svm-predict"
		echo "  -s control is structural models are produced"
		echo ""
		echo "pGenTHREADER Example: genthreader.sh -i test.fa -j test"
		echo "pDomTHREADER Example: genthreader.sh -i test.fa -j test -d"
		echo ""
		exit 1;;
	[?]) ERROR=1;;
	esac
done

if [ $ERROR == 1 ]
then
	echo >&2 "Usage: $0 -i file -j 'job_name' [-d] [-m] [-r] [-R] [-c] [-l] [-C CSA_file] [-L /libsvm/] [-b /blast_bin_path] [-u /uniref90_path] [-p /psipred_data_path] [-P /psipred_bin_path] [-t /tdb_path]"
	exit 1
fi

#override some settings if some clashing settings are set
if [ $TYPE == 1 ]
then
	echo "pDomTHREADER job running so -c set to 0"
	BUILD_CONTACTS=0
fi

if [ $MSA == 0 ]
then
	MSA_COMPLETE=0
	RELATIVES_COMPLETE=0
fi

if [ $LIGAND_PRED ]
then
	if [ $TYPE != 1 ]
	then
	echo "Ligand binding prediction selected so -c is set to 1"
	BUILD_CONTACTS=1
	fi
fi
#make a masked copy of input file
$PSIP/pfilt -f $FSA > $JOB.fsa

export TDB_DIR=$TDB
export THREAD_DIR=./data

echo started `date` $HOST > $JOB.pgt.log

if [ $MSA == 1 ]
then
 #Input file is an MSA
	echo "Running initial MSA psiblast"
 	if [ $MSA_COMPLETE == 1 ]
	then
		cp $FSA ./tmp_db.fa
		echo "$PSIB/makeblastdb -in $FSA -out tmp_db"
		$PSIB/makeblastdb -in $FSA -out tmp_db
		echo "$PSIB/psiblast -db ./tmp_db -in_msa $FSA -inclusion_ethresh 0.001 -out_pssm $JOB.chk -num_iterations 2 -num_alignments 0 >& $JOB.blast"
		$PSIB/psiblast -db ./tmp_db -in_msa $FSA -inclusion_ethresh 0.001 -out_pssm $JOB.chk -num_iterations 2 -num_alignments 0 >& $JOB.blast
		rm tmp_db*
		echo "chk file prepared"

	else
		echo "$PSIB/psiblast -db $DB -in_msa $FSA -inclusion_ethresh 0.001 -out_pssm $JOB.chk -num_iterations 3 -num_alignments 0 >& $JOB.blast"
		$PSIB/psiblast -db $DB -in_msa $FSA -inclusion_ethresh 0.001 -out_pssm $JOB.chk -num_iterations 3 -num_alignments 0 >& $JOB.blast
		echo "Finished 3 iterations PSI-BLAST with MSA"
	fi
else
 #Input file is a single fasta sequence echo "$PSIB/psiblast -db $DB -query $FSA -inclusion_ethresh 0.001 -out_pssm $JOB.chk -in_pssm $JOB.iter3.chk -num_iterations 3 -num_alignments 0 >& $JOB.blast"
 echo "Running initial fasta psiblast"
 echo "$PSIB/psiblast -db $DB -query $JOB.fsa -inclusion_ethresh 0.001 -out_pssm $JOB.chk -num_iterations 3 -num_alignments 0 >& $JOB.blast"
 $PSIB/psiblast -db $DB -query $JOB.fsa -inclusion_ethresh 0.001 -out_pssm $JOB.chk -num_iterations 3 -num_alignments 0 >& $JOB.blast
 echo "Finished 3 iterations PSI-BLAST"
fi

#make a profile matrix
$PGT/chkparse $JOB.chk > $JOB.mtx
#exit
#rename checkpoint file
mv $JOB.chk $JOB.iter3.chk
mv $JOB.mtx $JOB.iter3.mtx

echo "Running PSIPRED"

$PSIP/psipred $JOB.iter3.mtx $PDATA/weights.dat $PDATA/weights.dat2 $PDATA/weights.dat3 > $JOB.pgen.ss
$PSIP/psipass2 $PDATA/weights_p2.dat 1 1.0 1.0 $JOB.pgen.ss2 $JOB.pgen.ss > $JOB.horiz

if [ ! -s "$JOB.pgen.ss2" ]
then
   echo "PSIPRED failed ... exiting early"
   echo "PSIPRED failed ... exiting early" >> $JOB.pgt.log
   exit;
fi

echo "Finished PSIPRED"


if [ $TYPE == 1 ]
then
	#We're going to run domthreader
	echo "Matching domTHREADER tdb files"
	$PGT/pseudo_bas_dom -c11.0 -C20 -h0.2 -F$JOB.pgen.ss2 $JOB.iter3.mtx $JOB.pdom.pseudo $DATA/$DOMCHAIN > $JOB.pdom.output

	if [ ! -s "$JOB.pdom.pseudo" ]
	then
   		echo "pseudobas failed"
    	echo "pseudobas failed" >> $JOB.pdt.log
    	exit;
	fi
echo "svm_prob_dom"
	$PGT/svm_prob_dom $JOB.pdom.pseudo | sort -k 2,2rn -k 6,6rn -k 5,5g  > $JOB.pdom.presults

	if [ ! -s "$JOB.pdom.presults" ]
	then
    	echo "sortprob failed or no hits returned"
    	echo "sortprob failed" >> $JOB.pdt.log
    	exit;
	fi
echo "pseudo bas dom 2"
	#No multiple alignments yet, watch this space.
	$PGT/pseudo_bas_dom -S -p -c11.0 -C20 -h0.2 -F$JOB.pgen.ss2 $JOB.iter3.mtx $JOB.pdom.pseudo $JOB.pdom.presults > $JOB.pdom.align

else
	#We're going to run genthreader
	if [ $RELATIVES_COMPLETE == 0 ]
	then
		#$PSIB/blastpgp -a 4 -F T -t 1 -i $JOB.fsa -R $JOB.iter3.chk -d $DB -j 3 -v 5000 -b 0 -h 0.001 -C $JOB.chk > /dev/null
		echo "$PSIB/psiblast -db $DB -inclusion_ethresh 0.001 -out_pssm $JOB.chk -in_pssm $JOB.iter3.chk -num_iterations 3 -num_alignments 0 >& $JOB.blast"
		$PSIB/psiblast -db $DB -inclusion_ethresh 0.001 -out_pssm $JOB.chk -in_pssm $JOB.iter3.chk -num_iterations 3 -num_alignments 0 >& $JOB.blast

		echo "Finished 6 iterations of PSI BLAST"

		#make PSSM from 6th iteration
		$PGT/chkparse $JOB.chk > $JOB.mtx
		# $PSIB/makemat -P $JOB

		#rename checkpoint file
 		mv $JOB.chk $JOB.iter6.chk
  		mv $JOB.mtx $JOB.iter6.mtx
	fi

	if [ $RELATIVES_COMPLETE == 1 ]
	then
		#rename checkpoint file
		cp $JOB.iter3.chk $JOB.iter6.chk
  		cp $JOB.iter3.mtx $JOB.iter6.mtx
	fi

	# Run PGenThreader process
	echo "Matching GenTHREADER tdb files"
	$PGT/pseudo_bas  -c11.0 -C20 -h0.2 -F$JOB.pgen.ss2 $JOB.iter6.mtx $JOB.pgen.pseudo $DATA/$PSICHAIN

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

	#Ok we still don't have anything to do here about returning alignments. Answers on a postcard to the usual address
	if [ $BUILD_MODELS == 1]
  then
	  $PGT/pseudo_bas -m$JOB -S -p -c11.0 -C20 -h0.2 -F$JOB.pgen.ss2 $JOB.iter6.mtx $JOB.pgen.pseudo $JOB.pgen.presults > $JOB.pgen.align
  else
		$PGT/pseudo_bas -S -p -c11.0 -C20 -h0.2 -F$JOB.pgen.ss2 $JOB.iter6.mtx $JOB.pgen.pseudo $JOB.pgen.presults > $JOB.pgen.align
	fi
fi

#Here's a whole load of ancilliary bonus processes that you can calculate
#"for free" if you want
if [ $BUILD_MULTI == 1 ]
then
	echo "$BLAST/build_multi_alignment.pl $JOB $FSA"
	$BLAST/build_multi_alignment.pl $JOB $FSA
fi

if [ $BUILD_CONTACTS == 1 ]
then
	echo "$PGT/GenAlignmentHandler/main.rb ./$JOB.pgen ./$CSA"
	$PGT/GenAlignmentHandler/main.rb ./$JOB.pgen ./$CSA
fi

if [ $LIGAND_PRED == 1 ]
then
	echo "$PGT/predict_contacts.pl $LIBSVM $DATA $JOB > output"
	$PGT/predict_contacts.pl $LIBSVM $DATA $JOB > output
	echo "mv output $JOB.pgen.contactcons"
	mv output $JOB.pgen.contactcons
	#echo "$PGT/process_contacts.pl"
fi
echo Finished `date` >> $JOB.pdt.log
