#!/usr/bin/perl -w

use strict;
use FileHandle;
use DirHandle;
use English;
use Data::Dumper;

my $blastDir = "/home/dbuchan/genome3d/blast_out/";
my $fastaDir = "/home/dbuchan/genome3d/fasta/";
my $domDir = "/home/dbuchan/genome3d/domthreader_out/";

my $drBlast = new DirHandle($blastDir);
my $fhBlastAlignOut = new FileHandle("blast_aligns.txt","w");
my $hBlastData ={};
my $length = 0;
my $ID = '';
my $hPDomData = {};
my $fhSSF = new FileHandle("genome3d.ssf","w");
my $fhAligns = new FileHandle("genome3d_aligns.txt","w");
while(my $file = $drBlast->read)
{
	if($file =~ /^\./){next;}
	#if($file !~ /^85\./){next;}
	my $pgen_file = $file;
	$pgen_file =~ s/bls/pdom.presults/;
	my $align_file = $pgen_file;
	$align_file =~ s/presults/align/;
	
	if(-e $blastDir.$file && -e $domDir.$pgen_file && -e $domDir.$align_file)
	{
	
	
	print $fhBlastAlignOut $file."\n";
	$length = read_fasta($file);
	$hBlastData = read_blast_data($file);
	print $fhBlastAlignOut "-------\n";
	$hPDomData = read_pdom_data($file);
	#remove bad overlaps
	remove_low_overlaps();
	print_ssf();
	print_alignments();
	}
}

sub print_alignments
{
	foreach my $id (keys %$hBlastData)
	{
		print $fhAligns $ID." ".$id."\n";
		print $fhAligns $hBlastData->{$id}{ALIGNMENT_HEADER}."\n";
		print $fhAligns $hBlastData->{$id}{ALIGNMENT}."\n";
	}
	foreach my $id (keys %$hPDomData)
	{
		print $fhAligns $ID." ".$id."\n";
		print $fhAligns $hPDomData->{$id}{ALIGNMENT_HEADER}."\n";
		print $fhAligns $hPDomData->{$id}{ALIGNMENT}."\n";
	
	}
}

sub print_ssf
{
	foreach my $id (keys %$hBlastData)
	{
		print $fhSSF $ID." ".$id." 0 0 0 0 0 0 0 0 ".$hBlastData->{$id}{EVAL}." 0.00 0.00 1 ".$hBlastData->{$id}{START}.":".$hBlastData->{$id}{STOP}."\n";
	}
	foreach my $id (keys %$hPDomData)
	{
		print $fhSSF $ID." ".$id." 0 0 0 0 0 0 0 0 ".$hPDomData->{$id}{PVAL}." 0.00 0.00 1 ".$hPDomData->{$id}{START}.":".$hPDomData->{$id}{STOP}."\n";
	}
}

sub remove_low_overlaps
{
	foreach my $id (keys %$hBlastData)
	{
		my $align_length = $hBlastData->{$id}{STOP}-$hBlastData->{$id}{START};
		my $ratio = $align_length/$length;
		
		if($ratio < 0.4)
		{
			delete $hBlastData->{$id};
		}
	}
	
	foreach my $id (keys %$hPDomData)
	{
		my $align_length =$hPDomData->{$id}{STOP}-$hPDomData->{$id}{START};
		my $ratio = $align_length/$length;
		
		if($ratio < 0.4)
		{
			delete $hPDomData->{$id};
		}
	}
}



$fhBlastAlignOut->close;

sub read_fasta
{
	my ($ffile) = @ARG;
	$ffile =~ s/bls/fasta/;
	my $fhIn = new FileHandle($fastaDir.$ffile,"r");
	my $seq = '';
	while(my $line = $fhIn->getline)
	{
		if($line =~ /^>(.+)/){$ID = $1;next;}
		chomp $line;
		$seq.=$line;
		
	}
	my $length = length $seq;
	
	return($length);
}

sub read_blast_data
{
	my ($file) = @ARG;
	my $fhIn = new FileHandle($blastDir.$file,"r");
	
	my $passed_count = 0;
	my $found_align = 0;
	my $hData = {};
	my $current_id = '';
	while(my $line = $fhIn->getline)
	{
		chomp $line;
		if($line =~ /^\s\s(\Sdb\|.{7}).+?\s+(\d+|\d+\.\d+)\s+(.+)/)
		{
			my $id = $1;
			my $score = $2;
			my $eval = $3;
			#print $id." ".$score." ".$eval."\n";
			if($eval == 0.00 || $eval < 0.00005)
			{
				$passed_count++;
				$current_id = "CATHDOM|".$id.":".$passed_count;
				$hData->{$current_id}{SCORE} = $score;
				$hData->{$current_id}{EVAL} = $eval;
				$hData->{$current_id}{START} = 0;
				$hData->{$current_id}{STOP} = 0;
				$hData->{$current_id}{ALIGNMENT} = '';
				$hData->{$current_id}{DOMAINID} = $id;
			}
		}
		if($line =~ /^>\s(.{11})\s/)
		{
			$current_id = "CATHDOM|".$1.":".($found_align+1);
			if($found_align < $passed_count)
			{
				$hData->{$current_id}{ALIGNMENT_HEADER} = $line;
			}
			$found_align++;
						
		}
		if($line =~ /^(Query|Sbjct)\s/ && $found_align <= $passed_count)
		{
			$hData->{$current_id}{ALIGNMENT} .= $line."\n";
			if($line =~ /^Query\s+(\d+)\s.+\s(\d+)/)
			{
				my $tmp_start = $1;
				my $tmp_stop = $2;
				if($hData->{$current_id}{START} == 0)
				{
					$hData->{$current_id}{START} =$tmp_start;
				}
				if($tmp_stop > $hData->{$current_id}{STOP})
				{
					$hData->{$current_id}{STOP} = $tmp_stop;
				}
				
			}
		}
		
	}
	#print Dumper $hData;
	return($hData);
}

sub read_pdom_data
{
	my ($pfile) = @ARG;
	
	$pfile =~ s/bls/pdom.presults/;
	my $fhIn = new FileHandle($domDir.$pfile,"r");
	my $hit_count = 0;
	my $hData = {};
	while(my $line = $fhIn->getline)
	{
		if($line =~ /^CERT|^HIGH/)
		{
			$hit_count++;
			chomp $line;
			my $aEntries = [];
			@$aEntries = split /\s+/, $line;
			my $current_id = "PDOM|".@$aEntries[11].":".$hit_count;
			$hData->{$current_id}{SCORE} = @$aEntries[1];
			$hData->{$current_id}{PVAL} = @$aEntries[2];
			$hData->{$current_id}{START} = @$aEntries[9];
			$hData->{$current_id}{STOP} = @$aEntries[10];
			$hData->{$current_id}{DOMAINID} = @$aEntries[11];
			
		}
	}
	
	$pfile =~ s/presults/align/;
	$fhIn = new FileHandle($domDir.$pfile,"r");
	my $align_count = 0;
	my $current_id = '';
	my $id = '';
	while(my $line = $fhIn->getline)
	{
		if($line =~ />>>\sAlignment\swith\s(.+):/)
		{
			$align_count++;
			$id = $1;
			$current_id = "PDOM|".$id.":".$align_count;
			if(exists $hData->{$current_id})
			{
				$hData->{$current_id}{ALIGNMENT_HEADER} = ">".$current_id;
				$hData->{$current_id}{ALIGNMENT} = '';
			}
		}
		if(exists $hData->{$current_id})
		{
			if($line =~ /^Query/)
			{
				$hData->{$current_id}{ALIGNMENT}.=$line;
			}
			if($line =~ /^$id/)
			{
				$hData->{$current_id}{ALIGNMENT}.=$line;
			}
		}
	}
	return($hData);
}
