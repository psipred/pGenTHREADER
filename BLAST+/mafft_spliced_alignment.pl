#!/usr/bin/perl -w
#
# This script takes an alignment a coordinate series and cuts out the initial
# and trailing segments and then MAFFT aligns the target sequence to the msa sub
# segment. Finally it rebuilds these into a finished alignment.
#
# User needs to set the $MAFFT_DIR to the correct location to MAFFT
#
# useage: ./mafft_spliced_alignment.pl <ALIGNMENT FASTA> <NEW SEQUENCE> <START> <STOP>
#
# ./mafft_spliced_alignment.pl multiple_alignment.fa new_sequnce.fa 103 156
#
################################################################################
#
# Author: Daniel Buchan
# Date: 19th of May 2011
# Version: 1.1
# This code is distrubted under GPL v3.0
# http://www.gnu.org/licenses/gpl.html
#
################################################################################


use strict;
use English;
use Data::Dumper;
use FileHandle;

#my $MAFFT_DIR = "/webdata/binaries/current/mafft/bin";
my $MAFFT_DIR = "/scratch0/NOT_BACKED_UP/dbuchan/mafft-6.849-without-extensions/bin";
#chdir "/webdata/tmp/NewPredServer/";

#Setting these as global variables!
my $alignment = $ARGV[0];
my $seq = $ARGV[1];
my $start= $ARGV[2];
my $stop = $ARGV[3];
my $residue_count = 0;

if(!defined $ARGV[0] || !defined $ARGV[1] || !defined $ARGV[2] || !defined $ARGV[2] )
{
	print_useage();
}

#Open the alignment and divide it into three files
#Strictly you could just hold the sequences in memory for the whole process
#but if you've got a v.v.large alignment doing it on the disk is unlikely
#to run into memory problems
divide_alignment();

#run mafft on our alignment_segment.fa
run_mafft();

#Now we open our 3 segment files, concatentate the sequences and pad the
#newly aligned segment.
recombine_segments();

#here we clear up the temporary files.
clean();

sub print_useage
{
	print "\nuseage:\nmafft_spliced_alignment.pl <ALIGNMENT FASTA> <NEW SEQUENCE> <START COORD> <STOP COORD>\n";
	print "\nexample:\n./mafft_spliced_alignment.pl multiple_alignment.fa new_sequnce.fa 103 156\n";
	exit;
}

sub clean
{
	`rm leading_segment.fa`;
	`rm new_insert.fa`;
	`rm trailing_segment.fa`;
	`rm alignment_segment.fa`;
	#`rm realigned.fasta`;
}

sub recombine_segments
{
	tidy_insert();
	my $fhOut = new FileHandle("realigned.fasta","w");
	my $fhLeading_Segment = new FileHandle("leading_segment.fa","r");
	my $fhAlignment_Segment = new FileHandle("new_insert.fa","r");
	my $fhTrailing_Segment = new FileHandle("trailing_segment.fa","r");
	
	my $seqs = 0;
	my $full_length = 0;
	my $max_length = 0;
	
	while(my $line = $fhAlignment_Segment->getline)
	{
		my $lead_line = $fhLeading_Segment->getline;
		my $trail_line = $fhTrailing_Segment->getline;
		chomp $lead_line;
		if($line =~ /^>/)
		{
			print $fhOut $line;
		}
		else
		{
			chomp $line;
			print $fhOut $lead_line.$line.$trail_line;
		}
	}
	
	$fhLeading_Segment->close;
	$fhAlignment_Segment->close;
	$fhTrailing_Segment->close;
}

sub tidy_insert
{
	my $fhInsert = new FileHandle("new_insert.fa","r");
	my $fhOut = new FileHandle("tmp_insert.fa","w");
	my $seq_count=0;
	while(my $line = $fhInsert->getline)
	{
		if($line =~ /^>/)
		{
			if($seq_count != 0)
		    {
				print $fhOut "\n";
			}
			print $fhOut $line;
			$seq_count++;
		}
		else
		{
			chomp $line;
			print $fhOut $line;
		}
	}
	
	`mv tmp_insert.fa new_insert.fa`;
	$fhInsert->close;
	$fhOut->close;
}

sub divide_alignment
{
	my $fhIn = new FileHandle($alignment,"r");
	
	my $fhLeading_Segment = new FileHandle("leading_segment.fa","w");
	my $fhAlignment_Segment = new FileHandle("alignment_segment.fa","w");
	my $fhTrailing_Segment = new FileHandle("trailing_segment.fa","w");
	
	my $alignment_count=0;
	my $leading_length = 0;
	my $trailing_length = 0;
	while(my $line = $fhIn->getline)
	{
		#print $line;
		if($line =~ /^>/)
		{
			$leading_length = 0;
			$trailing_length = 0;
		   if($alignment_count != 0)
		   {
			   print $fhLeading_Segment "\n";
			   print $fhAlignment_Segment "\n";
			   print $fhTrailing_Segment "\n";
		   }
		   $alignment_count++;
		   
		   print $fhLeading_Segment $line;
		   print $fhAlignment_Segment $line;
		   print $fhTrailing_Segment $line;
		   $residue_count = 0;
		   
		}
		else
		{
		   chomp $line;
		   
		   my $aResidues = [];
		   @$aResidues = split //, $line;
		   foreach my $residue (@$aResidues)
		   {
			   $residue_count++;
			   if($residue_count < $start)
			   {
				   print $fhLeading_Segment $residue;
				   $leading_length++;
			   }
			   elsif($residue_count > $stop)
			   {
				   print $fhTrailing_Segment $residue;
				   $trailing_length++;
			   }
			   else
			   {
				   print $fhAlignment_Segment $residue;
			   }
		   }
		}
	}
	
	print $fhLeading_Segment "\n>\n";
	for(my $i = 0; $i<$leading_length; $i++)
	{
		print $fhLeading_Segment "-";
	}
	print $fhLeading_Segment "\n";
	
	print $fhTrailing_Segment "\n>\n";
	for(my $i = 0; $i<$trailing_length; $i++)
	{
		print $fhTrailing_Segment "-";
	}
	print $fhTrailing_Segment "\n";
	
	
	$fhIn->close;
	$fhLeading_Segment->close;
	$fhAlignment_Segment->close;
	$fhTrailing_Segment->close;
}

sub run_mafft
{
	my $cmd = $MAFFT_DIR."/mafft --quiet --add $seq alignment_segment.fa > new_insert.fa";
	print $cmd."\n";
	`$cmd`;
}
