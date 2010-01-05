#!/usr/bin/perl

# PROGRAM  : simrnaseq_generate_reads.pl
# AUTHOR   : Dave Sturgill

use strict;
use FileHandle;
use Getopt::Long;
use POSIX qw(ceil floor);


my $cov;
my $bp;
my $ends;
my $modeldir;
my $outdir;
my $type;
my $fast;
my $total = 0;
my $verbose;
my $varcov;

my $seqfilename;



# Specify cmdline options and process command line...
my $options_okay = GetOptions (
    # Application-specific options...
    'c=f'     => \$cov,     # --c option expects a floating point num
    'b=i'     => \$bp,     # --b option expects an integer
    'ends=i'    => \$ends,    # --ends option expects an integer
    'd=s'		=> \$modeldir, # --models directory is a string
    'o=s'		=> \$outdir, # --models directory is a string
    't=s'		=> \$type, # --tx, intron, intergenic
    'fast=s'		=> \$fast, # --fast, (fasta or fastq), is a string (a or q)
    'varcov=s'		=> \$varcov, # --verbose, is a string (T or F)
    'verbose=s'		=> \$verbose, # --verbose, is a string (T or F)
    'i=s'     => \$seqfilename     # --in option expects a string for sequence file name


);


chomp $seqfilename;
mkdir $outdir;

usage () if (! $options_okay);
	
sub usage
{
	print "Unknown option: @_\n" if ( @_ );
	print "usage: program [-c=coverage] [-b=base pairs ][-e 1=single ends 2=paired ends] [-d=directory where models are] [-o=output dir] [-t=tx,intron,intergenic] [-fast={a for fasta, q for fastq}] <filename>\n";
	exit;
}

unless (($ends < 3) && ($ends > 0)) {
	print "ends must be 1 or 2\n";
	exit;
}

unless ($cov > 0) {die "Don't forget to enter a valid coverage you want\n"};
unless ($bp > 0) {die "Don't forget to enter a valid read length\n"};

my $rtxfilename;

if ($type eq "intron") {
	#$seqfilename =  "annotation/dmel-all-intron-r5.21.fasta";
	#$rtxfilename = "$modeldir"."/rand_transcripts.bed"; ##!
	print "Intron reads selected\n";
	print "NOTE! You must have a file 'rand_transcripts.bed' in your directory\n";
	print "this file is from your random transcript selection, and this script\n";
	print "will attempt to only generated reads for introns that were in your\n";
	print "random transcript selection\n";
	sleep 2;
	$rtxfilename = "rand_transcripts.bed"; ##!
	print "seqfilename is $seqfilename\n";
}
if ($type eq "intergenic") {
	#$seqfilename =  "annotation/dmel-all-intergenic-r5.21.fasta";
	print "seqfilename is $seqfilename\n";
}

unless (($fast eq "a") || ($fast eq "q")) {
	print "'fast' must be 'a' or 'q'\n";
	exit;
}

unless (($verbose eq "T") || ($verbose eq "F")) {
	print "verbose must be 'T' or 'F'\n";
	exit;
}

unless (($varcov eq "T") || ($varcov eq "F")) {
	print "varcov must be 'T' or 'F'\n";
	exit;
}

if ($varcov eq "T") { 
	$varcov = "var";
} else {
	$varcov = "";
}

my $originalcov = $cov;

my $SIMREAD;
my $SIMREAD1;
my $SIMREAD2;
my $ANSWERS = new FileHandle ">"."$outdir"."/answers_tab_$varcov$cov"."x.txt" or die "can't open answers_tab.txt"; #OUTPUT: table of transcripts
$ANSWERS->printf("FB_txid\tFB_id\tChr\tTx_length\tNumreads\tRPK\n");  # Print header for my "answers" table

if ($ends == 1) {
	if ($fast eq "q") { $SIMREAD = new FileHandle ">"."$outdir"."/sim_"."$ends"."_"."$varcov$cov"."x.fastq" or die "can't open read output file"; } #OUTPUT: sim. reads
	if ($fast eq "a") { $SIMREAD = new FileHandle ">"."$outdir"."/sim_"."$ends"."_"."$varcov$cov"."x.fa" or die "can't open read output file"; } #OUTPUT: sim. reads
}
if ($ends == 2) {
	if ($fast eq "q") { $SIMREAD1 = new FileHandle ">"."$outdir"."/sim_"."$ends"."_"."$varcov$cov"."x_1.fastq" or die "can't open read output file"; } #OUTPUT: sim. reads
	if ($fast eq "q") { $SIMREAD2 = new FileHandle ">"."$outdir"."/sim_"."$ends"."_"."$varcov$cov"."x_2.fastq" or die "can't open read output file"; } #OUTPUT: sim. reads
	if ($fast eq "a") { $SIMREAD1 = new FileHandle ">"."$outdir"."/sim_"."$ends"."_"."$varcov$cov"."x_1.fa" or die "can't open read output file"; } #OUTPUT: sim. reads
	if ($fast eq "a") { $SIMREAD2 = new FileHandle ">"."$outdir"."/sim_"."$ends"."_"."$varcov$cov"."x_2.fa" or die "can't open read output file"; } #OUTPUT: sim. reads

}



my $revcom;
my $randrev;
my $read;

my @f;

################################################
################################################
#
# 				Subroutines
#
################################################
################################################

my @mmprobs;
my %mmbase; 
my $j;
my %qualpwm;
my %indist_weights;


#~~~~~~~~~~~~~~~~~~~~~~~
# pretty print
#~~~~~~~~~~~~~~~~~~~~~~~


sub print_sequence {
	my($sequence, $length) = @_;
	for (my $pos=0 ; $pos < length($sequence) ; $pos += $length) {
		#$RANDCDS->printf(substr($sequence, $pos, $length));
		#$RANDCDS->printf("\n");
	}
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~
# For converting mate pair distance weights
# Adapted from perl cookbook
# http://docstore.mik.ua/orelly/perl/cookbook/ch02_11.htm
#~~~~~~~~~~~~~~~~~~~~~~~~~~~

# weight_to_dist: takes a hash mapping key to weight and returns
# a hash mapping key to probability
sub weight_to_dist {
    my %weights = @_;
    my %dist    = ();
    my $total   = 0;
    my ($key, $weight);
    local $_;

    foreach (values %weights) {
        $total += $_;
    }

    while ( ($key, $weight) = each %weights ) {
        $dist{$key} = $weight/$total;
    }

    return %dist;
}

# weighted_rand: takes a hash mapping key to probability, and
# returns the corresponding element
sub weighted_rand {
    my %dist = @_;
    my ($key, $weight);

    while (1) {                     # to avoid floating point inaccuracies
        my $rand = rand;
        while ( ($key, $weight) = each %dist ) {
            return $key if ($rand -= $weight) < 0;
        }
    }
}

# This is to generate random coverage, with gaussian distribution
sub gaussian_rand {
    my ($u1, $u2);  # uniformly distributed random numbers
    my $w;          # variance, then a weight
    my ($g1, $g2);  # gaussian-distributed numbers

    do {
        $u1 = 2 * rand() - 1;
        $u2 = 2 * rand() - 1;
        $w = $u1*$u1 + $u2*$u2;
    } while ( $w >= 1 );

    $w = sqrt( (-2 * log($w))  / $w );
    $g2 = $u1 * $w;
    $g1 = $u2 * $w;
    # return both if wanted, else just one
    return wantarray ? ($g1, $g2) : $g1;
}


#~~~~~~~~~~~~~~~~~~~~~~~
# round
#~~~~~~~~~~~~~~~~~~~~~~~

 sub round
 {
     my($number) = shift;
     return int($number + .5 * ($number <=> 0));
 }
 

#~~~~~~~~~~~~~~~~~~~~~~~
# mismatch generator
#~~~~~~~~~~~~~~~~~~~~~~~


my %testmmhash; # For printing later

sub mmgen {

	my ($sequence) = @_;
	my $position;
	my $randtemp;
	my $mmatch;
	my $newbase;
	my $oldbase;
	my $randnum;
	my $checkstring = "";
	my @nucs = ('A','C','G','T');
	my $checkstring;
	## Iterate over each position in the sequence
	## Use probabilities to make a mismatch
	for ($position = 0; $position < length($sequence); ++$position) {
		$randtemp = rand();
		$mmatch = 0;
		if ($randtemp < (@mmprobs[$position]/100)) {
			$mmatch += 1;
		}
		$oldbase = substr($sequence,$position, 1);
		if ($mmatch > 0) {  ## If a mismatch is called, what to change base to
		
			## $mmbase{$j}{$k}
			## this object is a hash of hash. j is ref, k is read
			
			#For refbase = A
			#T: 0.186246603700349 N: 0.00187605123560616 C: 0.411502134816923 G: 0.400375210247121 
			
			# Convert to a hash of just probs for old base
			my %temphash;
			
			if ($oldbase ne "N") {
			for my $k ( keys %{ $mmbase{$oldbase} } ) {
				$temphash{$k} = $mmbase{$oldbase}{$k};
			}
			
			$newbase = weighted_rand(%temphash);
			} else {
			$newbase = $oldbase;
			}
			
			$checkstring = $checkstring . $newbase;
			
			$testmmhash{"$oldbase".">"."$newbase"} += 1;
		} else {
			$checkstring = $checkstring . $oldbase;
		
		}
		
		
		}
	return($checkstring);
}


#~~~~~~~~~~~~~~~~~~~~~~~
# Quality score generator
#~~~~~~~~~~~~~~~~~~~~~~~

sub qualgen {
	my ($mer) = @_;
	my $position;
	my $maxprob;
	my $total;
	my $randqual;
	my $randtemp;
		
	my $qualsum = 0;
	my $cumulativeprob = 0;
	my $qualfound;
	
	my $qualstring = "";
	my $mostlikely_bypos = "B"; # This is to plug in most common quality score if an error occurs
								# A warning will print if it is used.
	my $maxprob = 0;
	
	my %tempqualhash;
	
	for ($position = 0; $position < $mer; ++$position) {		
		$maxprob = 0;
		foreach $j (keys %qualpwm) {	
			$tempqualhash{$j} = $qualpwm{$j}[$position];
			if ($qualpwm{$j}[$position] > $maxprob) {
				$maxprob = $qualpwm{$j}[$position];
				$mostlikely_bypos = $j;
				
			}
		}	
		$randqual = weighted_rand(%tempqualhash);
		if (($randqual =~ /\s+/) || (length($randqual) < 1)) {		
			print "\n*Warning - had to substitute $mostlikely_bypos for $randqual\n";
			$randqual = $mostlikely_bypos;		
		}		
		$qualstring = $qualstring.$randqual;	
	}	
	return($qualstring);		
}



##################################################
##################################################
#~~~~~~~~~~~~~~~~~~~~~~~
# read generator
#~~~~~~~~~~~~~~~~~~~~~~~
##################################################
##################################################

sub generate_reads {
	my($sequence, $id, $mer, $numreads, $fast) = @_;
		if ($verbose eq "T") {print "d\n";}	
	my $count2;
	my $position;
	my $randtemp;
	my $randqual;
	my $randstart;
	my $qualstring;
	my $seqlen = length($sequence);
	#print "*c\n";
	my $readid;	
	
	$sequence =~ tr/acgtn/ACGTN/;
	
	for ($count2=1; $count2<$numreads; $count2++) {
		if ($verbose eq "T") {print "e\n";}	
		if ($fast eq "q") {	$readid = "@"."$id"."_"."$count2" };
		if ($fast eq "a") {	$readid = ">"."$id"."_"."$count2" };
		#print "$readid\n";
		$SIMREAD->printf($readid);
		$SIMREAD->printf("\n");
		
		$randstart = int(rand($seqlen-$mer));
		$read = substr($sequence, $randstart, $mer);
		if ($verbose eq "T") {print "f\n";}	
		if (length($read) != $mer) {
			print "Read length is wrong!\n";
			print $read, "\n";
			sleep 2;
		}
		
		#---------------------------------
		# Randomly get forward or reverse
		$randrev = rand();

		if ($randrev > 0.5) { 

			$revcom = reverse $read;
			$revcom =~ tr/ACGTacgtnN/TGCATGCANN/;
			$read = $revcom;
		}

		if ($verbose eq "T") {print "g\n";}	
		$SIMREAD->printf(mmgen($read));
		$SIMREAD->printf("\n");
		
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# Now make a quality string
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		if ($fast eq "q") {	
			$readid = "+"."$id"."_"."$count2";
			$SIMREAD->printf($readid);
			$SIMREAD->printf("\n");	
			$qualstring = qualgen($mer);
			if (length($qualstring) != $mer) {
				print "Qualstring length is wrong!\n";
				print $qualstring, "\n";
				sleep 2;
			}
			print { $SIMREAD } $qualstring;
			$SIMREAD->printf("\n");
		}
	}
}

# Sample Bowtie map file

#HWUSI-EASXXX:4:1:0:1368#0/1 Ê Ê + Ê Ê Ê 2R Ê Ê Ê16938356 NAGGGCTCCTCGCCGTACATGTCGGCCAATTTGCGGAAGCGCGGTCCAAAGTTGGATCGGTAGTCGAAGTTGAGAT Ê 0 Ê Ê Ê 0:G>N,56:C>T,57:A>C 
#HWUSI-EASXXX:4:1:0:1368#0/2 Ê Ê - Ê Ê Ê 2R Ê Ê Ê16938454 ÊCGGACGAGAGAGGCTGCCATCGGAGTTGCCGTCACCTTCGTACGCGTAATGCCGCACATCGTCCACGGTTGTGGC Ê Ê0 Ê Ê Ê 69:G>A,70:T>G,74:A>G
sub generate_pe_reads {
	my($sequence, $id, $mer, $numreads, $fast) = @_;
	my $count2;
	my $position;
	my $randtemp;
	my $randqual;
	my $qualstring;
	my $innerd;
	my $randstart;
	my $randstrand;
	my $seqlen = length($sequence);
	my $readid;	
	
	$sequence =~ tr/acgtn/ACGTN/;	
	
	for ($count2=1; $count2<$numreads; $count2++) {
		
		$innerd = weighted_rand(%indist_weights);
		while ($seqlen < (($mer * 2) + $innerd)) { 
		$innerd = weighted_rand(%indist_weights);
		}
		if ($fast eq "q") {	$readid = "@"."$id"."_"."$count2"."#0/1" };
		if ($fast eq "a") {	$readid = ">"."$id"."_"."$count2"."#0/1" };
		$SIMREAD1->printf($readid);
		$SIMREAD1->printf("\n");
		
		$randstart = int(rand($seqlen - $mer - $mer - $innerd));
		$read = substr($sequence,$randstart, $mer);
		if (length($read) != $mer) {
			print "Read length is wrong!\n";
			print $read, "\n";
			sleep 2;
		}
		
		
		#---------------------------------
		$randstrand = "+";
		$SIMREAD1->printf(mmgen($read));
		$SIMREAD1->printf("\n");
		
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# Now make a quality string
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		if ($fast eq "q") {
			$readid = "+"."$id"."_"."$count2"."#0/1";
			$SIMREAD1->printf($readid);
			$SIMREAD1->printf("\n");
	
			$qualstring = qualgen($mer);
			
			if (length($qualstring) != $mer) {
				print "Qualstring length is wrong!\n";
				print $qualstring, "\n";
				sleep 2;
			}
			print { $SIMREAD1 } $qualstring;
			$SIMREAD1->printf("\n");
		
		}
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# Now print the 2nd read
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		if ($fast eq "q") {	$readid = "@"."$id"."_"."$count2"."#0/2" };
		if ($fast eq "a") {	$readid = ">"."$id"."_"."$count2"."#0/2" };
		$SIMREAD2->printf($readid);
		$SIMREAD2->printf("\n");
		
		$read = substr($sequence,$randstart + $mer + $innerd, $mer);
		# Reverse the read 
		$revcom = reverse $read;
		$revcom =~ tr/ACGTacgtnN/TGCATGCANN/;
		$read = $revcom;
		
		$SIMREAD2->printf(mmgen($read));
		$SIMREAD2->printf("\n");

		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# Now make a quality string
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		if ($fast eq "q") {
			$readid = "+"."$id"."_"."$count2"."#0/2";
			$SIMREAD2->printf($readid);
			$SIMREAD2->printf("\n");
	
			$qualstring = qualgen($mer);
			
			if (length($qualstring) != $mer) {
				print "Qualstring length is wrong!\n";
				print $qualstring, "\n";
				sleep 2;
			}
			print { $SIMREAD2 } $qualstring;
			$SIMREAD2->printf("\n");

		}
		
	}
}


################################################
################################################
#
# 			Read in Models
#		(Qualities, mismatches, inner d)
#
################################################
################################################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read in quality PWM
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

my $pwmfile = "$modeldir"."/qualitiespwm.txt";
my $pwmline;

my $qualkey;

my $i;
my $j;
my $k;



open (PWMFILE, $pwmfile);

while ($pwmline = <PWMFILE>) {
	if ($verbose eq "T") {print "Reading in quality pwm\n";}
	@f = split /\t+/,$pwmline;
	$qualkey = shift @f;
	$i = 0;
	for (@f) {
	$qualpwm{$qualkey}[$i] = $_;
	$i += 1;
	
	}


}

close (PWMFILE);


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read in mismatch prob PWM
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


my $mmfile = "$modeldir"."/mmprobs.txt";
my $mmline;


open (MMFILE, $mmfile);

while ($mmline = <MMFILE>) {
	if ($verbose eq "T") {print "Reading in mismatch pwm\n";}
	@f = split /\t+/,$mmline;
	@mmprobs[$f[0]] = $f[1];
	
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read in mismatch type PWM
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

my %mmtypes;

my $mmfile = "$modeldir"."/mmtypes.txt";
my $mmline;


open (MMFILE, $mmfile);

while ($mmline = <MMFILE>) {
	if ($verbose eq "T") {print "Reading in mismatch type probs\n";}
	@f = split /\s+/,$mmline;
	$mmtypes{$f[0]} = $f[1];
	
}

sleep 1;
my %mmbasetotals;
my $refbase;
my $readbase;

foreach $j (keys %mmtypes) {
	if ($verbose eq "T") {print "Converting mmtype probs\n";}	
	#print "$j $mmtypes{$j}\n";
	$refbase = substr($j,0,1);
	$readbase = substr($j,2,1);
	$mmbase{$refbase}{$readbase} = $mmtypes{$j};
	$mmbasetotals{$refbase} += $mmtypes{$j};
	
	
}

### Convert to probabilities for each base
# Code is 'base in ref' > 'base in read'
# So the base in ref = my read generated from transcript
# 'base in read' will be what it changes to after mismatch generation

my $mmtypecount = 0;

for $j (keys %mmbase) {
    for $k ( keys %{ $mmbase{$j} } ) {
         $mmbase{$j}{$k} = $mmbase{$j}{$k}/$mmbasetotals{$j};
    }
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read in mate inner distance distribution
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if ($ends == 2) {

my $indistfile = "$modeldir"."/mateinnerdists.txt";
my $indistline;

open (INDIST, $indistfile);

while ($indistline = <INDIST>) {

	@f = split /\t+/,$indistline;
	$indist_weights{$f[0]} = $f[1];
	
}
}

close (INDIST);

my $weightsum = 0;

## Convert to prob. distribution
if ($ends == 2) {
	%indist_weights = weight_to_dist(%indist_weights);
}

################################################
################################################
#
# 		Read in Transcipt Files
#	(and call the read generator subs)
#
################################################
################################################


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# If you are doing introns, make a hash of FBtr ids, so that only introns
# from the selected transcripts are used.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

my %rtxs;
my $rtxline;
my $fbtr;

if ($type eq "intron") {
	print "Trying to load rtxfile $rtxfilename\n";
	open (RTXFILE, $rtxfilename);	
	while ($rtxline = <RTXFILE>) {
		@f = split /\t+/,$rtxline;		
		if ($f[3] =~ /FBtr[0-9]+/) {			
			$fbtr = $&;
			$rtxs{$fbtr} = 1;
			print "adding $fbtr to rtxfile hash\n";
		} else {	
			print "Can't find FBtr number $f[3]\n $rtxline\n";
			exit;	
		}		
	}
}


close (RTXFILE);

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load in sequence file
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

print "Starting to load the sequence file $seqfilename \n";
open (SEQFILE, $seqfilename);

my $seq;
my %seqhash;
my $txid = 0;
my $tx;
my $previd;
my $header;
my $prevheader = "";
my $fbtr;
while ($seq = <SEQFILE>) {
	if ($seq =~ /^>/) {
		$header = $seq;
		if ($txid > 0) {
			if ($prevheader =~ /FBtr[0-9]+/) {$fbtr = $&;}  ## This is only used if doing intron file
			
			if (($type ne "intron") || (($type eq "intron") && (exists($rtxs{$fbtr})))) {
			print "adding transcript # $txid to hash\n";
			$seqhash{$txid}[0] = $previd;
			$seqhash{$txid}[1] = $tx;
			$seqhash{$txid}[2] = $prevheader;
			$tx = "";
			}

		}
	$txid += 1;
	$previd = $seq;
	$prevheader = $header;
	
	} else {
	
	$seq =~ s/\n//;
	$tx = $tx.$seq;
		
	}
	
}

close (SEQFILE);


if ($prevheader =~ /FBtr[0-9]+/) {$fbtr = $&;}  ## This is only used if doing intron file

if (($type ne "intron") || (($type eq "intron") && (exists($rtxs{$fbtr})))) {
	$seqhash{$txid}[0] = $previd;
	$seqhash{$txid}[1] = $tx;
	$seqhash{$txid}[2] = $prevheader;
	print "added transcript # $txid to hash\n";
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Iterate over sequences and generate reads
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


my $k;
my $txstring;
my $hstring;

my $range;
my $minimum;
my $coverage;

my $numreads;

my $FB;
my $fbtr;
my $chr;
my $covtimeslength;
my $txlength;
my $RPK;

my $mean;
my $sdev;
my $newcov;

my $innerd = 0;

my $lengthcutoff;

if ($ends == 1) {
	$lengthcutoff = 75;
	
} else {
	$lengthcutoff = 350;
}

#Coverage is the average number of reads representing a given nucleotide in the reconstructed sequence. It can be calculated from the length of the original genome (G), the number of reads(N), and the average read length(L) as NL / G. For example, a hypothetical genome with 2,000 base pairs reconstructed from 8 reads with an average length of 500 nucleotides will have 2x redundancy.

# To calculate reads from the other data points:
# N = (Cov * G) / L
# coverage * G = NL
# cov = NL/G
my $kcount = 0;

foreach $k (sort keys %seqhash) {
	if ($verbose eq "T") {print "a\n";}	
	$kcount += 1;
	$txstring = $seqhash{$k}[1];
	$hstring = $seqhash{$k}[0];
	 if ($hstring =~ /(FBtr[0-9]+).+loc\=([a-zA-Z0-9\_]+)\:([a-zA-Z]*)\(*([0-9\.\,]+).+(FBgn[0-9]{7})/) {

		 $fbtr = $1;	
		 $chr = $2;
		 $FB = $5;  # This is the parent gene

	} else {
	
		$hstring =~ s/\>//;
		$hstring =~ s/\n//;
		$fbtr = $hstring;
	
	}
	
	unless ($type eq "tx") {
			if ($hstring =~ /(FBgn[0-9]{7})/) {
	 
			 $fbtr = $type."_".$kcount;
			
			 
			 $FB = $&;  # This is the parent gene

			}

	
	}
	
	
	if ($varcov eq "var") {
	
		$mean = $originalcov;
		$sdev = 1;
		$newcov = gaussian_rand() * $sdev + $mean;
		printf("Random coverage is %.2f\n", $newcov);
		$cov = $newcov;
	
	}
	
	
	$txlength = length($txstring);
	$covtimeslength = $cov * $txlength;
	if ($ends == 1) {
	$numreads = round($covtimeslength / $bp);
	} else {
	$numreads = round($covtimeslength / ($bp * 2));	
	}
	$RPK = $numreads / $txlength * 1000;
	$RPK = sprintf("%.4f", $RPK);
	#print "If it was 36bp, it would be ",int(($cov * length($txstring)) / 36), " reads\n";

		if ($verbose eq "T") {print "b\n";}	

	
	if ($txlength > $lengthcutoff) {
		$ANSWERS->printf("$fbtr\t$FB\t$chr\t$txlength\t$numreads\t$RPK\n");  # Print my "answers" table
		if ($verbose eq "T") {print "c\n";}	
	
			if ($ends == 1) {
			#print "*b\t$k\n";
			if ($verbose eq "T") {print "Starting 'generate reads' sub $fbtr $numreads $txstring\n";}	
			generate_reads($txstring, $fbtr, $bp, $numreads, $fast);
			} else {
			#print "*b\t$k\n";
			generate_pe_reads($txstring, $fbtr, $bp, $numreads, $fast);
			}
			
			#generate_qualities($txstring, substr($hstring,0,length($hstring) - 1), $bp, $numreads);
	} else {
		$ANSWERS->printf("$fbtr\t$FB\t$chr\t$txlength\tskipped\ttoo short\n");  # Print my "answers" table
	}
	
	
}		

print "\n";

#== Below prints out the models
=pod
foreach $j (keys %mmtypes) {
	$refbase = substr($j,0,1);
	$readbase = substr($j,2,1);
	$mmbase{$refbase}{$readbase} = $mmtypes{$j};
	$mmbasetotals{$refbase} += $mmtypes{$j};	
}

for my $k ( keys (%testmmhash) ) {

	print "$k:",  $testmmhash{$k}, "\n";
}
=cut
