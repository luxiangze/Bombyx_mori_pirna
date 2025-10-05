#!/usr/bin/perl
use Getopt::Long;
$|=1;

$script=$0;
$script=~s/^.+[\\\/]//;
$info_text=
"
============================== NGS TOOLBOX ==============================
Program: collapse
Version: 2.1
LAST MODIFIED: 11. March 2015

Usage: perl $script -i input -o output

Input files can be in FASTA or FASTQ format. Avoid linebreaks within
sequences! Default output is 'collapse.out'. Sequence headers in the
output file will refer to sequence read counts, e.g.:
>18
TTAGCCGATACGTCGCATAGCTCGACTG

When using FASTQ input files the qualities will be set to 'B' indicating
unknown quality, e.g.:
@18
TTAGCCGATACGTCGCATAGCTCGACTG
+
BBBBBBBBBBBBBBBBBBBBBBBBBBBB

-h or -help will print this information.

Contact:
David Rosenkranz
Institute of Anthropology
Johannes Gutenberg University Mainz, Germany
email: rosenkranz\@uni-mainz.de

";

$output_file="collapse.out";
GetOptions
	(
	"help"=>\$print_info,
	"h"=>\$print_info,
	"i=s"=>\$input_file,
	"o=s"=>\$output_file,
	);

# print info
if($print_info)
	{
	print$info_text;
	exit;
	}

# check input format (FASTA/FASTQ)
print"\nChecking input file format...";
open(IN,$input_file)||die print"\nCould not open input file $input_file.\n$!\n\n";
$i=0;
$fasta_heads=0;
$fastq_heads=0;
while(<IN>)
	{
	$i++;
	last if($i==1000);
	if($_=~/^>/)
		{
		$fasta_heads++;
		}
	elsif($_=~/^@/||$_=~/^\+/)
		{
		$fastq_heads++;
		}
	}
close IN;
if($fasta_heads>$fastq_heads)
	{
	$format='FASTA';
	}
elsif($fastq_heads>$fasta_heads)
	{
	$format='FASTQ';
	}
else
	{
	die print"\nUnknown file format of input file $input_file.\n\n";
	}
print" done. -> $format";

print"
============================== NGS TOOLBOX ==============================
Program: collapse
Version: 2.1
last modified: 11. March 2015
=========================================================================
";

# process input file
print"\nProcess sequences in $input_file...";
open(IN,$input_file)||die print"\nCould not open input file $input_file.\n$!\n\n";
open(OUT,">$output_file")||die print"\nCould not create output file $output_file.\n$!\n\n";

$i=0;
$reads=0;
%seqs=();
@seq_data=();
while(<IN>)
	{
	$_=~s/\s*$//;
	push(@seq_data,$_);
	if($format eq'FASTQ')
		{
		$i+=0.25;
		}
	elsif($format eq'FASTA')
		{
		$i+=0.5;
		}
	if($i==1)
		{
		$i=0;
		$reads++;
		$seqs{$seq_data[1]}++;
		@seq_data=();
		}
	}
close IN;
$non_id=keys%seqs;
if($format eq'FASTQ')
	{
	foreach$seq(keys%seqs)
		{
		print OUT"\@$seqs{$seq}\n$seq\n+\n";
		$seq=~s/./B/g;
		print OUT"$seq\n";
		}
	}
elsif($format eq'FASTA')
	{
	foreach$seq(keys%seqs)
		{
		print OUT"\>$seqs{$seq}\n$seq\n";
		}
	}
close OUT;

print"\nFound $reads reads representing $non_id non-identical sequences.\n\n";
exit;