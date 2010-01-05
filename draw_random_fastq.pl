#!/usr/bin/perl

# PROGRAM  : draw_random_reads_better.pl
# AUTHOR   : Dave
# Last Revision Date: July 22, 2009
# Last Revised By: Dave

# Pull random fastq reads
# Selection done WITHOUT replacement
# Usage: perl ~draw_random_reads.pl [filename] (numreads) (total lines in file)

# NOTE: If you do not enter the total lines in the file,
#       it will count them for you
#

use strict;
use FileHandle;
use POSIX qw(ceil floor);

my $filename = $ARGV[0];
my $totalreads = $ARGV[1];
my $totallines = $ARGV[2];

chomp $filename;
my $seq;

# Automatic line counting, unless the user entered a value

unless ($totallines > 0) {
	$totallines = 0;
	open (SEQFILE, $filename);	
	while ($seq = <SEQFILE>) {
		$totallines += 1;
	}
	close SEQFILE;
}

#print "There are $totallines total lines\n"; 

# Pick random line numers (divisible by 4);

my $cnt = 1;
my $randline;

my %randhash;

while ($cnt <= $totalreads) {

	$randline = int(rand($totallines - 4));
	
	# Make sure it is divisible by 4
	if ($randline % 4 == 0) {
	
		# Check to see it hasn't already been chosen
		
		unless (exists($randhash{$randline})) {
		
			$randhash{$randline} = $randline;
			# Also add the next three lines for 
			# read sequence, quality string header, quality string
			# This can be changed for fasta sequence
			#$randhash{$randline + 1} = $randline;
			#$randhash{$randline + 2} = $randline;
			#$randhash{$randline + 3} = $randline;
			$cnt += 1;
	
		}
	}
}



# Now print the reads you chose

my $linenum = 0;

open (SEQFILE, $filename);

while ($seq = <SEQFILE>) {
	
	if (exists($randhash{$linenum}) | exists($randhash{$linenum - 1}) | exists($randhash{$linenum - 2}) | exists($randhash{$linenum - 3} )) {
		print $seq;
	} 

	$linenum += 1;

}
close SEQFILE;






=pod 


open (SEQFILE, $filename);

my $seq;
my %seqhash;
my $txid = 0;
my $tx;
my $previd;
my $header;
my $prevheader = "";
my $refcount = 0;
my @keys;

=cut