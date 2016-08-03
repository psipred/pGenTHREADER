#!/usr/bin/perl -w
#
# This script takes the output of genthreader and GenAlignment handler
# and outputs the input files for libsvm for contact prediction.
#
######################################################################
#
# Version 1.0
# D.Buchan
#
# This software is distributed under the GPLv3
# http://www.gnu.org/licenses/gpl.html
#
######################################################################

use strict;
use FileHandle;
use DirHandle;
use Data::Dumper;
use English;

my $LIBSVM = $ARGV[0];
my $DATA = $ARGV[1];
my $FILE = $ARGV[2];

if(-e "./".$FILE.".iter6.mtx")
{
	`./bin/mtx2con $FILE.iter6.mtx > $FILE.cons`;
}

if(-e "./".$FILE.".pgen.contactcons")
{
	my($aResidues,$aConsensus,$hContact_con) = get_data("./", $FILE);
	if($aResidues =~ /NO/ || $aConsensus =~ /NO/){print STDERR "No residue or contact data\n";exit;}
	my $hPred_residues = print_data($aResidues,$aConsensus,$hContact_con,$FILE.".libsvm");
	run_libsvm($LIBSVM, $DATA, $FILE.".libsvm", $FILE);
	my $hPredictions = read_libsvm($FILE.".libsvmoutput");
	output_predictions($hPred_residues,$hPredictions, "./", $FILE);
	#print Dumper $hPredictions;
	#print Dumper $aPred_residues;
}

sub output_predictions
{
	my($hPred_residues, $hPredictions, $dir, $file) = @ARG;

	my $fhIn = new FileHandle($dir.$file.".pgen.contactcons","r");
	my $seq_length = 0;
	while(my $line = $fhIn->getline)
	{
		if($line !~ /^TYPE/)
		{
			$seq_length = 0;
			my $aEntries = [];
			@$aEntries = split /\s+/, $line; 
			my $aRes = [];
			@$aRes = split //, @$aEntries[4];
			foreach my $res (@$aRes)
			{
				#print $res;
				$seq_length++;
			}
		}
		print $line;
	}
	print "---------\n";
	#9 tabs
	print "CONTACTS\t\t\t\t\t\t\t";
	my $pred_count = 0;
	for(my $i = 0; $i < $seq_length; $i++)
	{
		if(exists $hPred_residues->{$i})
		{
			print $hPredictions->{$pred_count}{INITIAL};
			$pred_count++;
		}
		else
		{
			print "N";
		}
	}
	print "\n";
	print "PROB\t\t\t\t\t\t\t\t";
	$pred_count = 0;
	for(my $i = 0; $i < $seq_length; $i++)
	{
		if(exists $hPred_residues->{$i})
		{
			if($hPredictions->{$pred_count}{PROB} =~ /0\.(\d)/)
			{
				print $1;
			}
			else
			{
				print STDERR "CATASTROPHIC ERROR: Probability of 1\n";
				exit;
			}
			$pred_count++;	
		}
		else
		{
			print "N";
		}
	}
	print "\n";
	print "CONTACTS\t\t\t\t\t\t\t\t";
	$pred_count = 0;
	for(my $i = 0; $i < $seq_length; $i++)
	{
		if(exists $hPred_residues->{$i})
		{
			print $hPredictions->{$pred_count}{PRED};
			$pred_count++;
		}
		else
		{
			print "N";
		}
	}
	print "\n";
	
}

sub read_libsvm
{
	my($file) = @ARG;

	my $fhIn = new FileHandle($file,"r");
	my $count = 0;
	my $hPredictions = {};
	while(my $line = $fhIn->getline)
	{
		if($line =~ /^-1\s.+\s(.+)/)
		{
			$hPredictions->{$count}{INITIAL} = 0;
			$hPredictions->{$count}{PRED} = 0;
			$hPredictions->{$count}{PROB} = $1;
			$count++;
		}
		if($line =~ /^1\s(.+)\s.+/)
		{
			$hPredictions->{$count}{PROB} = $1;
			$hPredictions->{$count}{INITIAL} = 1;
			if($1 >= 0.6)
			{
				$hPredictions->{$count}{PRED} = 1;	
			}
			else
			{
				$hPredictions->{$count}{PRED} = 0;
			}
			
			$count++;
		}
		
	}
	return($hPredictions);
}

sub run_libsvm
{
	my($libsvm, $data, $file, $name) = @ARG;
	my $model = $data."/probability_model";
	print "$libsvm/svm-predict -b 1 ./$file $model ./$name.libsvmoutput\n";
	`$libsvm/svm-predict -b 1 ./$file $model ./$name.libsvmoutput`;
}

#go through data outputting and output
sub print_data
{
	my($aResidues,$aConsensus,$hContact_con,$file) = @ARG;
	my $fhOut = new FileHandle($file,"w");
		
	my $seq_length = scalar @$aResidues;
	my $hPredictions = {};
	#print $seq_length."\n";
	for(my $i = 0; $i < $seq_length; $i++)
	{
		if($i < 2 || $i > ($seq_length-3)){next;}#skip the ones we can't calculate for due to window size
		my $line  = '0 ';
		
		my $consensus_target = @$aConsensus[$i];
		my $average_cons_score = (@$aConsensus[($i-2)]+@$aConsensus[($i-1)]+@$aConsensus[$i]+@$aConsensus[($i+1)]+@$aConsensus[($i+2)])/5;
		my $best_cons_score = 0;
		for(my $it = ($i-2); $it <= ($i+2); $it++)
		{
			if(@$aConsensus[$it] > $best_cons_score)
			{
				$best_cons_score = @$aConsensus[$it];
			}
		}
		
		my $ContactCon = $hContact_con->{$i+1};
		my $Ave_ContactCon = ($hContact_con->{$i-1}+$hContact_con->{$i}+$hContact_con->{$i+1}+$hContact_con->{$i+2}+$hContact_con->{$i+3})/5;
		my $Best_ContactCon = 0.00000000000;
		my $a = $hContact_con->{$i-1};
		my $b = $hContact_con->{$i};
		my $c = $hContact_con->{$i+1};
		my $d = $hContact_con->{$i+2};
		my $e = $hContact_con->{$i+3};
		#so, this ought to be a for loop but I can't for the life of me work out why it wasn't working, see commented out section below
		if($a > $Best_ContactCon){$Best_ContactCon = $a;}
		if($b > $Best_ContactCon){$Best_ContactCon = $b;}
		if($c > $Best_ContactCon){$Best_ContactCon = $c;}
		if($d > $Best_ContactCon){$Best_ContactCon = $d;}
		if($e > $Best_ContactCon){$Best_ContactCon = $e;}
		
		#for(my $it = ($i-1); $it <= ($i+3); $it++)
		#{
		#	my $test_value_a = sprintf("%.4f",$hContact_con->{$i+1});
		#	my $test_value_b = sprintf("%.4f",$Best_ContactCon);
		#	
		#	#if($hContact_con->{$i+1} > $Best_ContactCon)
		#	if($test_value_a > $test_value_b)
		#	{
		#		$Best_ContactCon = $hContact_con->{$i+1};
		#	}
		#	print $file." ".$hContact_con->{$i+1}." ".$test_value_a." ".$test_value_b." ".$Best_ContactCon."\n";
		#}
		if($Best_ContactCon > 0.0 && $ContactCon != 0)
		{
			$hPredictions->{$i} = 1;
			print $fhOut $line."1:".$consensus_target." 2:".$average_cons_score." 3:".$best_cons_score." 4:".$ContactCon." 5:".$Ave_ContactCon." 6:".$Best_ContactCon."\n";
		}
	}
	$fhOut->close;
	return($hPredictions);
}

#Here we read in the bit list of true contacts (doesn't include the ignored ones)
sub get_true_contacts
{
	my($file) = @ARG;
	
	my $fhIn = new FileHandle($file,"r");
	my $hData = {};
	while(my $line = $fhIn->getline)
	{
		chomp $line;
		my $aEntries = [];
		@$aEntries = split /\s+/, $line;
		my $pdb = (uc @$aEntries[0]).":".@$aEntries[1];
		@$aEntries[4] =~ s/[A-Z]//;
		$hData->{$pdb}{@$aEntries[4]} = 1;
	}
	return($hData);
}


#we read in the data that we're going to process. The Fasta file of the query sequence, the consensus scores and the contact_consensus ratio
sub get_data
{
	my($dir, $file_name) = @ARG;
	
	my $aResidues = read_fsa($dir."/".$file_name.".fsa");
	if($aResidues =~ /^NO/){return("NO")}
	my $aConsensus = read_cons($dir."/".$file_name.".cons");
	if($aConsensus =~ /^NO/){return("NO")}	
	my $hContact_con = read_contactcons($dir."/".$file_name.".pgen.contactcons");
	if($hContact_con =~ /^NO/){return("NO")}	
	return($aResidues,$aConsensus,$hContact_con)
}

#produce and array of residues for the query sequence
sub read_fsa
{
	my($file) = @ARG;
	my $fhIn;
	if(-e $file)
	{
		$fhIn = new FileHandle($file,"r");
	}
	else
	{
		return("NO");
	}
	my $seq = '';
	while(my $line = $fhIn->getline)
	{	
		chomp $line;
		if($line =~ /^>/){next;}
		$seq.=$line;
	}
	
	my $aResidues = [];
	@$aResidues = split //, $seq;
	return($aResidues)
}

#read the consensus string
sub read_cons
{
	my($file) = @ARG;
	
	my $fhIn = new FileHandle($file,"r");
	my $seq = '';
	my $line = $fhIn->getline;
	$line = $fhIn->getline;
	chomp $line;
	
	my $aConsensus = [];
	@$aConsensus = split //, $line;
	return($aConsensus);
}

#read contact consensus alignments, skipping the invalid ligands
sub read_contactcons
{
	my($file) = @ARG;
	my $hInvalid = {"SO4" => 1, 
	"GOL" => 1, 
	"CL" => 1, 
	"PO4" => 1, 
	"NAG" => 1, 
	"MG" => 1, 
	"ACT" => 1, 
	"EDO" => 1, 
	"PLP" => 1, 
	"HEM" => 1, 
	"MPD" => 1, 
	"CA" => 1, 
	"ACY" => 1, 
	"BME" => 1, 
	"MES" => 1, 
	"DMS" => 1, 
	"HG" => 1, 
	"CD" => 1, 
	"TRS" => 1, 
	"NAG NAG" => 1, 
	"FMT" => 1, 
	"NAG MAN" => 1, 
	"PG4" => 1, };
	if(! -e $file)
	{
		return("NO");
	}
	
	my $fhIn = new FileHandle($file,"r");
	
	my $aResults = [];
	my $line_count = 0;
	my $seq_length = 0;
	while(my $line = $fhIn->getline)
	{
		chomp $line;
		#if($line =~ /^TYPE/){next;}
		my $aEntries = [];
		@$aEntries = split /\s+/, $line;
		if($line =~ /^TYPE/)
		{
			$seq_length = length @$aEntries[4];
			next;
		}
		if(length @$aEntries[4] != $seq_length){next;}
		if(exists $hInvalid->{@$aEntries[1]}){next;}
		if(@$aEntries[3] =~ /CERT|HIGH/)
		{
			$line_count++;
			push @$aResults,  @$aEntries[4];
		}
	}
	if($line_count == 0)
	{
		return("NO");
	}
	
	my $hResults = {};
	$line_count = 0;
	#print Dumper @$aResults;
	foreach my $result (@$aResults)
	{
		$line_count++;
		my $aResidues = [];
		@$aResidues = split //, $result;
		my $res_count = 0;
		foreach my $res (@$aResidues)
		{
			$res_count++;
			$hResults->{$res_count} += $res;
		}
				
	}
	foreach my $res_index (keys %$hResults)
	{
		$hResults->{$res_index} = $hResults->{$res_index}/$line_count;
	}
	return($hResults);
}
