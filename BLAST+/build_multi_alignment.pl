#!/usr/bin/perl -w

#
# This script parses the alignments for the segments of the pdb files to
# realign into a multiple alignment and the coordinates in the template 
# sequence.
#

use strict;
use English;
use Data::Dumper;
use FileHandle;

#chdir "/webdata/tmp/NewPredServer/";

my $name = $ARGV[0];
my $alignment = $ARGV[1];

#`cp $alignment tmp_input_align.fa`;
#`echo ">Query" | cat > tmp.fa`;
#`cat tmp.fa $alignment > tmp_input_align.fa`;

#exit;
my $cert_count = 0;
my $high_count = 0;
my $med_count = 0;
my $hResults = {};
get_list();
my $med_master = $med_count;
my $high_master = $high_count;
my $cert_master = $cert_count;
my $hAlignments={};
get_data();
#print Dumper $hAlignments;
process_alignments();
clean();

sub clean
{
	#`rm tmp.fa`;
	if(-e "tmp_seq.fa")
	{
		`rm tmp_seq.fa`;
	}
	if(-e "tmp_input_align.fa")
	{
		`rm tmp_input_align.fa`;
	}
	if(-e "realigned.fasta")
	{
		`rm realigned.fasta`;
	}
}

sub get_list
{
	my $fhInput;
	
	if(-e $name.".pgen.presults")
	{
		$fhInput = new FileHandle($name.".pgen.presults","r");
	}
	if(-e $name.".pdom.presults")
	{
		$fhInput = new FileHandle($name.".pdom.presults","r");
	}
	if(-e $name.".presults")
	{#print "hello\n";
		$fhInput = new FileHandle($name.".presults","r");
	}
	
	#print $name.".presults\n";
	
	while(my $line = $fhInput->getline)
	{
		if($line =~ /^(CERT|HIGH|MEDIUM).+\s+(.+)$/)
		{
			my $type = $1;
			my $pdb = $2;
			if($type =~ /CERT/)
			{
				$cert_count++;
			}
			if($type =~ /HIGH/)
			{
				$high_count++;
			}
			if($type =~ /MEDIUM/)
			{
				$med_count++;
			}
			
			if($hResults->{$pdb})
			{
				$hResults->{$pdb}++;
			}
			else
			{
				$hResults->{$pdb}=1;
			}
		}
	}
}

sub process_alignments
{
	my $count = 0;
	my $cert_high = $cert_master+$high_master;
	my $cert_high_med = $cert_master+$high_master+$med_master;
	foreach my $alignment_count (sort {$a <=> $b} keys %$hAlignments)
	{
		#print $alignment_count."\n";
		if($alignment_count == 1)
		{
			output_first_alignment("tmp_input_align.fa");
			`cp tmp_input_align.fa realigned.fasta`;
			$count++;
		}
		else
		{
			my $pdb_start = 0;
			my $pdb_stop = 0;
			my $query_start = 0;
			my $query_stop = 0;
			my $alignment_start = 0;
			my $alignment_stop = 0;
			my $pdb_subseq = [];
			($pdb_start, $pdb_stop, $query_start, $query_stop, $pdb_subseq) = find_coords($hAlignments->{$alignment_count}{HIT_SEQ},$hAlignments->{$alignment_count}{QUERY_SEQ});
			#exit;
			#print $pdb_start."  ".$pdb_stop." : ".$query_start." ".$query_stop."\n";
			($alignment_start,$alignment_stop) = map_coords($query_start,$query_stop);
			#print $pdb_start."  ".$pdb_stop." : ".$alignment_start." ".$alignment_stop."\n";
			print_sub_seq($alignment_count."_".$hAlignments->{$alignment_count}{PDB},$pdb_subseq);
			#Ok here we do some dialign specific things
			print "mafft_spliced_alignment.pl tmp_input_align.fa tmp_seq.fa $alignment_start $alignment_stop\n";
			`./BLAST+/mafft_spliced_alignment.pl tmp_input_align.fa tmp_seq.fa $alignment_start $alignment_stop`;
			
			#`/webdata/binaries/current/psipred/bin/mafft_spliced_alignment.pl tmp_input_align.fa tmp_seq.fa $start $stop`;
			#exit;
			#exit;
			`cp realigned.fasta tmp_input_align.fa`;
			
			$count++;
		}
		
		if($count == $cert_master)
		{
			`cp realigned.fasta $name.pgen.cert.multialign`;
		}
		if($count == $cert_high && $high_master != 0)
		{
			`cp realigned.fasta $name.pgen.certhigh.multialign`;
		}
		if($count == $cert_high_med && $high_master != 0 && $med_master != 0)
		{
			`cp realigned.fasta $name.pgen.all.multialign`;
		}
	}
	
	
}

sub get_data
{
	my $fhInput;
	if(-e $name.".pgen.align")
	{
		$fhInput = new FileHandle($name.".pgen.align","r");
	}
	if(-e $name.".pdom.align")
	{
		$fhInput = new FileHandle($name.".pdom.align","r");
	}
	if(-e $name.".align")
	{
		$fhInput = new FileHandle($name.".align","r");
	}
	my $alignment_count = 0;
	my $alignment_total = $cert_master+$high_master+$med_master;
	my $pdb_id="";
	while(my $line = $fhInput->getline)
	{
		if($line =~ />>>\sAlignment\swith\s(.+):/)
		{
			$pdb_id = $1;
			$alignment_count++;
			if($alignment_count > $alignment_total)
			{
				last;
			}
		}
		
		if($line =~ /^$pdb_id\s+(.+)/)
		{
			my $aTmp = [];
			@$aTmp = split //, $1;
			foreach my $residue (@$aTmp)
			{
			   push @{$hAlignments->{$alignment_count}{HIT_SEQ}}, $residue;
			}
			$hAlignments->{$alignment_count}{PDB}=$pdb_id;
		}
		if($line =~ /^Query\s+(.+)/)
		{
			my $aTmp = [];
			@$aTmp = split //, $1;
			foreach my $residue (@$aTmp)
			{
			   push @{$hAlignments->{$alignment_count}{QUERY_SEQ}}, $residue;
			}
			}
		
		
	}
	
	return($hAlignments);
}

sub output_first_alignment
{
	my($file) = @ARG;
	
	my $fhOut = new FileHandle($file,"w");
	
	my $start = 0;
	my $stop = 0;
	my $count = 0;
	foreach my $residue (@{$hAlignments->{1}{QUERY_SEQ}})
	{
		$count++;
		if($residue !~ /-/)
		{
			if($start == 0)
			{
				$start = $count;
			}
			$stop = $count;
		}
	}
	
	print $fhOut ">Query\n";
	my $seq_count = 0;
	foreach my $residue (@{$hAlignments->{1}{QUERY_SEQ}})
	{
		$seq_count++;
		if($seq_count < $start){next;}
		if($seq_count > $stop){next;}
		print $fhOut $residue;
	}
	print $fhOut "\n";
	
	print $fhOut ">".$hAlignments->{1}{PDB}."\n";
	 $seq_count = 0;
	foreach my $residue (@{$hAlignments->{1}{HIT_SEQ}})
	{	
		$seq_count++;
		if($seq_count < $start){next;}
		if($seq_count > $stop){next;}
		print $fhOut $residue;
	}
	print $fhOut "\n";
	$fhOut->close;
}

sub write_anc
{
	my($start, $stop,$sub_seq_length) = @ARG;
	
	my $fhInput = new FileHandle("dialign_input.fasta","r");
	my $seq_count =0;
	while(my $line = $fhInput->getline)
	{
		if($line =~ /^>/){$seq_count++;}
	}
	
	my $fhOut = new FileHandle("dialign_input.anc","w");
	print $fhOut "1 ".$seq_count." ".$start." 1 ".$sub_seq_length." 10.00000\n";
	$fhOut->close;
}
sub print_sub_seq
{
	my($pdb_id,$pdb_subseq) = @ARG;
	my $fhOut = new FileHandle("tmp_seq.fa","w");
	print $fhOut ">".$pdb_id."\n";
	foreach my $residue (@$pdb_subseq)
	{
		print $fhOut $residue;
	}
	print $fhOut "\n";
	$fhOut->close;
}
sub find_coords
{
	my($pdb_seq, $query_seq) = @ARG;
	
	my $alignment_count = 0;
	my $pdb_count = 0;
	my $pdb_start = 0;
	my $pdb_stop = 0;
	
	my $query_start = 0;
	my $query_stop = 0;
	my $pdb_subseq = [];
	
	my $length = @$pdb_seq;
	#print Dumper $query_seq;
	foreach my $residue (@$pdb_seq)
	{
		if($residue !~ /-/ && @$query_seq[$alignment_count] !~ /-/)
		{
			push @$pdb_subseq,  $residue;
			
			if($pdb_start == 0)
			{
				$pdb_start=$alignment_count;
			}
			$pdb_stop= $alignment_count;
			
			#push @$pdb_subseq, @$pdb_seq[$alignment_count];
		}
		$alignment_count++;
	}
	
	$alignment_count = 0;
	foreach my $residue (@$query_seq)
	{
		if($residue !~ /-/)
		{
			if($query_start == 0)
			{
				$query_start=$alignment_count;
			}
			$query_stop=$alignment_count;
			
			#push @$pdb_subseq, @$pdb_seq[$alignment_count];
			$alignment_count++;
		}
	}
	
	return(($pdb_start+1), $pdb_stop, $query_start, ($query_stop+1), $pdb_subseq);
}

sub map_coords
{
	my ($start,$stop)  = @ARG;
	
	my $fhInput = new FileHandle("tmp_input_align.fa","r");
	my $seq_count = 0;
	my $seq = [];
	while(my $line = $fhInput->getline)
	{
		if($line =~ /^>/)
		{
			$seq_count++;
		}
		if($seq_count == 1)
		{
			if($line !~ /^>/)
			{
				chomp $line;
				my $aTmp = [];
				@$aTmp = split //, $line;
				foreach my $residue (@$aTmp)
				{
					#print $residue."\n";
					push @$seq, $residue;
				}
			}
		}
	}
	#print Dumper $seq;
	my $align_count = 0;
	my $residue_count = 0;
	my $start_changed = 0;
	my $stop_changed = 0;
	foreach my $residue (@$seq)
	{
		$align_count++;
		if($residue !~ /-/)
		{
			$residue_count++;
			if($residue_count == $start && $start_changed==0)
			{
				$start = $align_count;
				$start_changed=1;
			}
			if($residue_count == $stop  && $stop_changed==0)
			{
				$stop = $align_count;
				$stop_changed=1;
			}
		}
	}
	return($start, $stop);
}
