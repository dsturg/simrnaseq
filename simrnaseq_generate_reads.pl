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

my $sdev;

my %wighash;

my $mmtf;
my $qualspeed;
my @randquals = ();

my @coordstring;
my $optid;

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
    'sdev=f'		=> \$sdev, # st. dev of converage, floating
    'verbose=s'		=> \$verbose, # --verbose, is a string (T or F)
    'mmtf=s'		=> \$mmtf, # --mmtf, is a string (T or F)
    'qualspeed=s'   => \$qualspeed, # --qualspeed is a string (T or F)
    'optid=s'   => \$optid, # --id is an optional identifier for the run (string)
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

unless (($mmtf eq "T") || ($mmtf eq "F")) {
	print "mmtf must be 'T' or 'F'\n";
	exit;
}


unless (($qualspeed eq "T") || ($qualspeed eq "F")) {
	print "qualspeed must be 'T' or 'F'\n";
	exit;
}



my $covstring;

if ($varcov eq "T") { 
	$varcov = "var";
	if ($sdev eq "") { 
		print "No sdev specified, using 1 as default\n";
		sleep 1;
		$sdev = 1;
	}
	$covstring = "$type"."_"."$varcov$cov"."x_sd$sdev";

} else {
	$varcov = "";
	$covstring = "$type"."_"."$varcov$cov"."x_fixed";
}



my $originalcov = $cov;
my %junchash;

my $e;
if ($ends == 1) {$e = "se"}
if ($ends == 2) {$e = "pe"}


my $SIMREAD;
my $SIMREAD1;
my $SIMREAD2;
my $ANSWERS = new FileHandle ">"."$outdir"."/answers_tab$optid"."_"."$e$covstring".".txt" or die "can't open answers_tab.txt"; #OUTPUT: table of transcripts
$ANSWERS->printf("Feature_id\tFB_id\tChr\tTx_length\tNumreads\tRPK\tcov\n");  # Print header for my "answers" table


my $WIG = new FileHandle ">"."$outdir"."/basecoverage$optid"."_$e$covstring".".wig" or die "can't open wigfile"; #OUTPUT: wigle file

my $READBED = new FileHandle ">"."$outdir"."/readbed$optid"."_$e$covstring".".bed" or die "can't open readbed"; #OUTPUT: bed file

print { $READBED } "track name=simulatedReads description='Simulated Reads'\n";

my $JUNCTAB;

if ($type eq "tx") {
$JUNCTAB = new FileHandle ">"."$outdir"."/junctab$optid"."_$covstring".".txt" or die "can't open junction table"; #OUTPUT: junction table
}



if ($ends == 1) {
	if ($fast eq "q") { $SIMREAD = new FileHandle ">"."$outdir"."/sim$optid"."_"."$e"."_"."$covstring".".fastq" or die "can't open read output file"; } #OUTPUT: sim. reads
	if ($fast eq "a") { $SIMREAD = new FileHandle ">"."$outdir"."/sim$optid"."_"."$e"."_"."$covstring".".fa" or die "can't open read output file"; } #OUTPUT: sim. reads
}
if ($ends == 2) {
	if ($fast eq "q") { $SIMREAD1 = new FileHandle ">"."$outdir"."/sim$optid"."_"."$e"."_"."$covstring"."_1.fastq" or die "can't open read output file"; } #OUTPUT: sim. reads
	if ($fast eq "q") { $SIMREAD2 = new FileHandle ">"."$outdir"."/sim$optid"."_"."$e"."_"."$covstring"."_2.fastq" or die "can't open read output file"; } #OUTPUT: sim. reads
	if ($fast eq "a") { $SIMREAD1 = new FileHandle ">"."$outdir"."/sim$optid"."_"."$e"."_"."$covstring"."_1.fa" or die "can't open read output file"; } #OUTPUT: sim. reads
	if ($fast eq "a") { $SIMREAD2 = new FileHandle ">"."$outdir"."/sim$optid"."_"."$e"."_"."$covstring"."_2.fa" or die "can't open read output file"; } #OUTPUT: sim. reads

}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Read in GFF file (This is for getting the strand of single-exon genes)

my $gfffile = "/users/labuser/data/annotation/dmel-all-no-analysis-r5.17.gff";
my %gffhash;
my $gffline;
my @g = ();

open (GFFFILE, $gfffile);

while ($gffline = <GFFFILE>) {

	@g = split /\t+/,$gffline;
	#if (@g[2] eq "mRNA") {
	if (@g[2] eq "gene") {
		if (@g[8] =~ /ID\=([a-zA-Z0-9]+)/) {
			$gffhash{$1}[0] = @g[0];  #chr
			$gffhash{$1}[1] = @g[3];  #start
			$gffhash{$1}[2] = @g[4];  #stop
			$gffhash{$1}[3] = @g[6];  #strand
			
			#print "$1\t$gffhash{$1}[0]\t$gffhash{$1}[1]\t$gffhash{$1}[2]\t$gffhash{$1}[3]\n";
		}
	}
	


}


close (GFFFILE);
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



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
# histogram
#~~~~~~~~~~~~~~~~~~~~~~~

 sub histogram
 {
   my ($bin_width, @list) = @_;
 
   # This calculates the frequencies for all available bins in the data set
   my %histogram;
   $histogram{ceil(($_ + 1) / $bin_width) -1}++ for @list;
 
   my $max;
   my $min;
 
   # Calculate min and max
   while ( my ($key, $value) = each(%histogram) )
   {
     $max = $key if !defined($min) || $key > $max;
     $min = $key if !defined($min) || $key < $min;
   }
 
 
   for (my $i = $min; $i <= $max; $i++)
   {
     my $bin       = sprintf("% 10d", ($i) * $bin_width);
     my $frequency = $histogram{$i} || 0;
 
     $frequency = "#" x $frequency;
 
     print $bin." ".$frequency."\n";
   }
 
   print "===============================\n\n";
   print "    Width: ".$bin_width."\n";
   print "    Range: ".$min."-".$max."\n\n";
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
# explode coordinate string
#~~~~~~~~~~~~~~~~~~~~~~~

 sub coordexplode
 {
     my($join) = @_;
     my @coordstring;
     my $num1;
     
     my $i;
     my @matches;
     
     #$string = "12 34 56 78 90 98 76 54 32 10";
	(@matches) = ($join =~ /[0-9]+\.\.[0-9]+/g);
     
     $i = 0;
    foreach (@matches) {
 	  #print "* $_\n";
 	  if ($_ =~ /([0-9]+)\.\.([0-9]+)/) {
 	      	push(@coordstring, $1);
 	      	
     		for ($i = $2 - $1; $i >= 1; $i--) {
     			push(@coordstring, $coordstring[$#coordstring] + 1);
     			
     	}

 	  }
	}

     
     
     #if ($join =~ /([0-9]+)\.\.([0-9]+)/) {
     
     #	$i = 0;
     #	@coordstring[$i] = $1;
     #	while ($i < ($2 - $1)) {
     #		$i += 1;
     #		@coordstring[$i] = @coordstring[$i - 1] + 1
     #	}
     
    # }
     
     return @coordstring;
     
 }


#(4999006..5000837,5000900..5001458)

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
	my $i;
	## Iterate over each position in the sequence
	## Use probabilities to make a mismatch
	if ($verbose eq "T") {print "..1 *\n";}
	for ($position = 0; $position < length($sequence); ++$position) {
		$randtemp = rand();
		$mmatch = 0;
		if ($randtemp < (@mmprobs[$position]/100)) {
			$mmatch += 1;
		}
		$oldbase = substr($sequence,$position, 1);
		#print "..2 *\n";
		if ($mmatch > 0) {  ## If a mismatch is called, what to change base to
			if ($verbose eq "T") {print "..3 *\n";}
			## $mmbase{$j}{$k}
			## this object is a hash of hash. j is ref, k is read
			
			#For refbase = A
			#T: 0.186246603700349 N: 0.00187605123560616 C: 0.411502134816923 G: 0.400375210247121 
			
			# Convert to a hash of just probs for old base
			my %temphash;
			
			if ($oldbase =~ /[ACGT]/) {
				for my $k ( keys %{ $mmbase{$oldbase} } ) {
				if ($verbose eq "T") {print "..4 *\n";}
					$temphash{$k} = $mmbase{$oldbase}{$k}; #make a new hash for mms for this base
					
				}
				if ($verbose eq "T") {print "..5 *\n";}
				if ($verbose eq "T") {print "position is $position oldbase is $oldbase\n";}
				if ($verbose eq "T") {print "sequence is  $sequence\n";}
				$i = 0;
				for my $k (keys %temphash) {
					#print "$k\t$temphash{$k}\n";
					$i += $temphash{$k};
				}
				if ($verbose eq "T") {print "total is $i\n";}
				$newbase = weighted_rand(%temphash);
				if ($verbose eq "T") {print "..6 *\n";}
			} else {
				$newbase = "N";
				print "An $oldbase was mmgen'ed to N\n";
			}
			
			$checkstring = $checkstring . $newbase;
			if ($verbose eq "T") {print "..7 *\n";}
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

#track name=pairedReads description="Clone Paired Reads" useScore=1
#chr22 1000 5000 cloneA 960 + 1000 5000 0 2 567,488, 0,3512
#chr22 2000 6000 cloneB 900 - 2000 6000 0 2 433,399, 0,3601

sub generate_reads {
	my($sequence, $id, $mer, $numreads, $fast, $chr, $strand,$direction,$mmtf,$qualspeed) = @_;
		if ($verbose eq "T") {print "d\n";}
		#print " coordstring 0 is $coordstring[0]\n"; 
	my $count2;
	my $position;
	my $randtemp;
	my $randqual;
	my $randstart;
	my $qualstring;
	my $seqlen = length($sequence);
	my $readid;	
	my $count;
	my $genomecoord;
	
	my @junctions;
	my $numjuncs;
	
	my @readcoords = ();
	my $seqstrand;
	
	my $blockcoordstring;
	my $blocksizestring;
	my $tempstart;
	my $tempend;
	
	my $tempjunc;
	
	my $origstrand = $strand;
	
	my $j;
	my $prev;
	my $junckey;
	my $offset;
	
	$sequence =~ tr/acgtn/ACGTN/;
	
	for ($count2=1; $count2<$numreads; $count2++) { # iterate based on number of reads you want
		
		$strand = $origstrand;
		
		$blockcoordstring = "";
		$blocksizestring = "";
				
		$randstart = int(rand($seqlen-$mer)); # Note that zeros are possible, zero-based start
		
		
		
		if ($verbose eq "T") {print "$seqlen\t$mer\t$randstart\n";}

		@readcoords = (); # This will store each coordinate for this read
						  # It will be used to create a base-level wiggle track (+1 based)
						  # Note these coords are already in +1 format
				
		for ($count=0; $count<$mer; $count++) { # Get genomic coordinates for this read
			$genomecoord = @coordstring[$randstart + $count];
			# Add coordinates to hash
			$wighash{$chr}[$genomecoord] += 1; # These coords are independent of strand, so OK to instantiate this hash now			
			push(@readcoords, $genomecoord); #readcoords are still in 1-based system
			
		}
		
		my $gstart = @readcoords[0];
		my $rangestring = "$gstart..";
		@junctions = (); # This will store all the junctions spanned by this read
		foreach (@readcoords) { # Generate a rangestring, find junctions
			if ($_ > ($gstart + 1)) {
				push(@junctions,"$gstart".".."."$_");
				$rangestring = $rangestring."$gstart".",";
				$gstart = $_;
				$rangestring = $rangestring.$gstart."..";
				
			}
		$gstart = $_;
		}
		$rangestring = $rangestring."$readcoords[$#readcoords]";
				
		#---------------------------------
		# Randomly get forward or reverse
		# Figure out the strand
		
		$randrev = rand();

		if ($randrev > 0.5) { $strand = $strand * -1 } # Farther down, the read will be revcom'ed
		
		if ($strand < 0) {
			$rangestring = "_$chr"."\:minus(".$rangestring.")";
			$seqstrand = "-";
		} elsif ($strand > 0) {
			$rangestring = "_$chr"."\:plus(".$rangestring.")";	
			$seqstrand = "+";
		} else {
		print "I'm confused\n";
		sleep 2;
		}
		
		
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		## Print full readid to filehandle

		
		if ($fast eq "q") {	$readid = "@"."$id"."_"."$count2"; };
		if ($fast eq "a") {	$readid = ">"."$id"."_"."$count2"; };
		$readid = $readid.$rangestring;
		#print "readid here is $readid\n";
		#sleep 1;
		$SIMREAD->printf($readid); 
		$SIMREAD->printf("\n");
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		
		# Get read sequence
		
		$read = substr($sequence, $randstart, $mer);
		if ($verbose eq "T") {print "f\n";}	
		if (length($read) != $mer) {
			print "Read length is wrong!\n";
			print $read, "\n";
			sleep 2;
		}
		
		#---------------------------------
		# Reverse the read, based on random selector above

		if ($randrev > 0.5) { 
			$revcom = reverse $read;
			$revcom =~ tr/ACGTacgtnN/TGCATGCANN/;
			$read = $revcom;
		}

		#!  Printing to BED in zero-based system
		$readid =~ s/\>//;
		$readid =~ s/\@//;
				
		if (@junctions) { # Junctions array is NOT empty

			print { $READBED } "chr$chr\t",$readcoords[0] - 1,"\t$readcoords[$#readcoords]\t$readid\t"."1"."\t$seqstrand\t",$readcoords[0] - 1,"\t$readcoords[$#readcoords]\t255,0,0\t";
			$numjuncs = @junctions;
			print { $READBED } $numjuncs + 1,"\t";
			$tempstart = $readcoords[0];
			my $sizetemp;
			my $coordtemp;
			$blockcoordstring = "0";
			$blocksizestring = "";
			#print "\n** $tempstart to $readcoords[$#readcoords]\t$blockcoordstring\t$blocksizestring **\n";
			foreach (@junctions) {
				#$j += 1;
				$tempjunc = $_;
				$junckey = "chr$chr"."_".$_."_$direction";
				if ($tempjunc =~ /([0-9]+)\.\.([0-9]+)/) {$gstart = $1};
	
				$offset = $gstart - @readcoords[0] + 1;

				if (exists($junchash{$junckey})) {
					$junchash{$junckey}[0] += 1;
					$prev = pop @{$junchash{$junckey}};
					#print "prev for $junckey is $prev\n";
					
					push(@{$junchash{$junckey}}, "$prev\,".$offset);

				} else {
					$junchash{$junckey}[0] += 1;
					#print "added ", $offset, " to hash with key $junckey\n";
					push(@{$junchash{$junckey}},$offset);
					#sleep 1;
				}
				
				if ($verbose eq "T") {print $_."\n";}
					if ($_ = /([0-9]+)\.\.([0-9]+)/) {
						#print "1 is $1 and 2 is $2\n"; 
						$sizetemp = $1 - $tempstart + 1;
						$blocksizestring = $blocksizestring.$sizetemp."\,";
						$coordtemp = $2 - $readcoords[0];
						$blockcoordstring = $blockcoordstring."\,".$coordtemp;
						if ($verbose eq "T") {print "\n$tempstart\t$blockcoordstring\t$blocksizestring";}
						$tempstart = $2;
					}
				
			}
			$sizetemp = $readcoords[$#readcoords] - $tempstart + 1;
			$blocksizestring = $blocksizestring.$sizetemp;
			#$blockcoordstring = $blockcoordstring.",".$readcoords[$#readcoords] - $tempend;
			if ($verbose eq "T") {print "\n Final: $blockcoordstring\t$blocksizestring";}
			
			print { $READBED } "$blocksizestring\t$blockcoordstring\n";

		} else { # If there are no junctions, the line is simple
			print { $READBED } "chr$chr\t",$readcoords[0] - 1,"\t$readcoords[$#readcoords]\t$readid\t"."1"."\t$seqstrand\t",$readcoords[0] - 1,"\t$readcoords[$#readcoords]\t255,0,0\t1\t75\t0\n";



		}

		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# Mismatch generator 
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#print "readid *here* is $readid\n";
		#sleep 1;

		if ($verbose eq "T") {print "\nstart $randstart, strand $strand, $readid g\n$read\n";}	
		
		if ($mmtf eq "T") {
		$SIMREAD->printf(mmgen($read));
		} elsif ($mmtf eq "F") {
		$SIMREAD->printf($read);
		} else {
		print "I dont get it \n";
		sleep 3;
		}
		$SIMREAD->printf("\n");
		if ($verbose eq "T") {print "h\n";}	
		
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# Now make a quality string
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		if ($fast eq "q") {	
			#$readid = "+"."$id"."_"."$count2";
			$readid = "+".$readid;
			$SIMREAD->printf($readid);
			$SIMREAD->printf("\n");			
			if ($qualspeed eq "T") {			
				$qualstring = @randquals[int(rand(999))];
				#print "qualstring is $qualstring\n";
				print { $SIMREAD } $qualstring;
				$SIMREAD->printf("\n");		
			} else {
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
	

}




sub generate_pe_reads {
	my($sequence, $id, $mer, $numreads, $fast, $chr, $strand,$direction,$mmtf,$qualspeed) = @_;
		if ($verbose eq "T") {print "d\n";}
		#print " coordstring 0 is $coordstring[0]\n"; 
	my $count2;
	my $position;
	my $randtemp;
	my $randqual;
	my $randstart;
	my $qualstring;
	my $seqlen = length($sequence);
	my $readid;	
	my $count;
	my $genomecoord;
	
	my @junctions;
	my $numjuncs;
	
	my @readcoords = ();
	my $seqstrand;
	
	my $blockcoordstring;
	my $blocksizestring;
	my $tempstart;
	my $tempend;
	
	my $tempjunc;
	
	my $origstrand = $strand;
	
	my $j;
	my $prev;
	my $junckey;
	my $offset;
	
	my $innerd;
	
	$sequence =~ tr/acgtn/ACGTN/;
	
	for ($count2=1; $count2<$numreads; $count2++) { # iterate based on number of reads you want
		$innerd = weighted_rand(%indist_weights);
		while ($seqlen < (($mer * 2) + $innerd)) { 
		$innerd = weighted_rand(%indist_weights);
		}
		
		$strand = $origstrand;
		
		$blockcoordstring = "";
		$blocksizestring = "";
				
		#$randstart = int(rand($seqlen-$mer)); # Note that zeros are possible, zero-based start
		$randstart = int(rand($seqlen - $mer - $mer - $innerd));

		
		
		if ($verbose eq "T") {print "$seqlen\t$mer\t$randstart\n";}

		@readcoords = (); # This will store each coordinate for this read
						  # It will be used to create a base-level wiggle track (+1 based)
						  # Note these coords are already in +1 format
				
		for ($count=0; $count<$mer; $count++) { # Get genomic coordinates for this read
			$genomecoord = @coordstring[$randstart + $count];
			# Add coordinates to hash
			$wighash{$chr}[$genomecoord] += 1; # These coords are independent of strand, so OK to instantiate this hash now			
			push(@readcoords, $genomecoord); #readcoords are still in 1-based system
			
		}
		
		my $gstart = @readcoords[0];
		my $rangestring = "$gstart..";
		@junctions = (); # This will store all the junctions spanned by this read
		foreach (@readcoords) { # Generate a rangestring, find junctions
			if ($_ > ($gstart + 1)) {
				push(@junctions,"$gstart".".."."$_");
				$rangestring = $rangestring."$gstart".",";
				$gstart = $_;
				$rangestring = $rangestring.$gstart."..";
				
			}
		$gstart = $_;
		}
		$rangestring = $rangestring."$readcoords[$#readcoords]";
				
		#---------------------------------
		# Randomly get forward or reverse
		# Figure out the strand
		
		if ($strand < 0) {
			$rangestring = "_$chr"."\:minus(".$rangestring.")";
			$seqstrand = "-";
		} elsif ($strand > 0) {
			$rangestring = "_$chr"."\:plus(".$rangestring.")";	
			$seqstrand = "+";
		} else {
		print "I'm confused 1\n";
		sleep 2;
		}
		
		
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		## Print full readid to filehandle

		
		if ($fast eq "q") {	$readid = "@"."$id"."_"."$count2"; };
		if ($fast eq "a") {	$readid = ">"."$id"."_"."$count2"; };
		$readid = $readid.$rangestring."#0/1";
		#print "readid here is $readid\n";
		#sleep 1;
		$SIMREAD1->printf($readid); 
		$SIMREAD1->printf("\n");
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		
		# Get read sequence
		
		$read = substr($sequence, $randstart, $mer);
		if ($verbose eq "T") {print "f\n";}	
		if (length($read) != $mer) {
			print "Read length is wrong!\n";
			print $read, "\n";
			sleep 2;
		}
		
		#---------------------------------

		#!  Printing to BED in zero-based system
		$readid =~ s/\>//;
		$readid =~ s/\@//;
				
		if (@junctions) { # Junctions array is NOT empty

			print { $READBED } "chr$chr\t",$readcoords[0] - 1,"\t$readcoords[$#readcoords]\t$readid\t"."1"."\t$seqstrand\t",$readcoords[0] - 1,"\t$readcoords[$#readcoords]\t255,0,0\t";
			$numjuncs = @junctions;
			print { $READBED } $numjuncs + 1,"\t";
			$tempstart = $readcoords[0];
			my $sizetemp;
			my $coordtemp;
			$blockcoordstring = "0";
			$blocksizestring = "";
			#print "\n** $tempstart to $readcoords[$#readcoords]\t$blockcoordstring\t$blocksizestring **\n";
			foreach (@junctions) {
				#$j += 1;
				$tempjunc = $_;
				$junckey = "chr$chr"."_".$_."_$direction";
				if ($tempjunc =~ /([0-9]+)\.\.([0-9]+)/) {$gstart = $1};
	
				$offset = $gstart - @readcoords[0] + 1;

				if (exists($junchash{$junckey})) {
					$junchash{$junckey}[0] += 1;
					$prev = pop @{$junchash{$junckey}};
					#print "prev for $junckey is $prev\n";
					
					push(@{$junchash{$junckey}}, "$prev\,".$offset);

				} else {
					$junchash{$junckey}[0] += 1;
					#print "added ", $offset, " to hash with key $junckey\n";
					push(@{$junchash{$junckey}},$offset);
					#sleep 1;
				}
				
				if ($verbose eq "T") {print $_."\n";}
					if ($_ = /([0-9]+)\.\.([0-9]+)/) {
						#print "1 is $1 and 2 is $2\n"; 
						$sizetemp = $1 - $tempstart + 1;
						$blocksizestring = $blocksizestring.$sizetemp."\,";
						$coordtemp = $2 - $readcoords[0];
						$blockcoordstring = $blockcoordstring."\,".$coordtemp;
						if ($verbose eq "T") {print "\n$tempstart\t$blockcoordstring\t$blocksizestring";}
						$tempstart = $2;
					}
				
			}
			$sizetemp = $readcoords[$#readcoords] - $tempstart + 1;
			$blocksizestring = $blocksizestring.$sizetemp;
			#$blockcoordstring = $blockcoordstring.",".$readcoords[$#readcoords] - $tempend;
			if ($verbose eq "T") {print "\n Final: $blockcoordstring\t$blocksizestring";}
			
			print { $READBED } "$blocksizestring\t$blockcoordstring\n";

		} else { # If there are no junctions, the line is simple
			print { $READBED } "chr$chr\t",$readcoords[0] - 1,"\t$readcoords[$#readcoords]\t$readid\t"."1"."\t$seqstrand\t",$readcoords[0] - 1,"\t$readcoords[$#readcoords]\t255,0,0\t1\t75\t0\n";



		}

		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# Mismatch generator 
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#print "readid *here* is $readid\n";
		#sleep 1;

		if ($verbose eq "T") {print "\nstart $randstart, strand $strand, $readid g\n$read\n";}	
		
		if ($mmtf eq "T") {
		$SIMREAD1->printf(mmgen($read));
		} elsif ($mmtf eq "F") {
		$SIMREAD1->printf($read);
		} else {
		print "I dont get it 1 \n";
		sleep 3;
		}
		$SIMREAD1->printf("\n");
		if ($verbose eq "T") {print "h\n";}	
		
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# Now make a quality string
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		if ($fast eq "q") {	
			#$readid = "+"."$id"."_"."$count2";
			$readid = "+".$readid;
			$SIMREAD1->printf($readid);
			$SIMREAD1->printf("\n");			
			if ($qualspeed eq "T") {			
				$qualstring = @randquals[int(rand(999))];
				#print "qualstring is $qualstring\n";
				print { $SIMREAD1 } $qualstring;
				$SIMREAD1->printf("\n");		
			} else {
				$qualstring = qualgen($mer);
				if (length($qualstring) != $mer) {
					print "Qualstring length is wrong!\n";
					print $qualstring, "\n";
					sleep 2;
				}
				print { $SIMREAD1 } $qualstring;
				$SIMREAD1->printf("\n");
			}
			
		
		}
######################################################################		
######################################################################	
######################################################################		
######################################################################	
		#-----------------------------------
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# Now print the 2nd read
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
######################################################################		
######################################################################	
######################################################################		
######################################################################	
		#print "** I am doing the 2nd read\n";
		$strand = $strand * -1;
		
		$blockcoordstring = "";
		$blocksizestring = "";
		
		@readcoords = ();
		$randstart = $randstart + $mer + $innerd;

		for ($count=0; $count<$mer; $count++) { # Get genomic coordinates for this read
			$genomecoord = @coordstring[$randstart + $count];
			# Add coordinates to hash
			$wighash{$chr}[$genomecoord] += 1; # These coords are independent of strand, so OK to instantiate this hash now			
			push(@readcoords, $genomecoord); #readcoords are still in 1-based system
			
		}
		
		my $gstart = @readcoords[0];
		my $rangestring = "$gstart..";
		@junctions = (); # This will store all the junctions spanned by this read
		foreach (@readcoords) { # Generate a rangestring, find junctions
			if ($_ > ($gstart + 1)) {
				push(@junctions,"$gstart".".."."$_");
				$rangestring = $rangestring."$gstart".",";
				$gstart = $_;
				$rangestring = $rangestring.$gstart."..";
				
			}
		$gstart = $_;
		}
		$rangestring = $rangestring."$readcoords[$#readcoords]";
				
		#---------------------------------
		# Figure out the strand
		
		if ($strand < 0) {
			$rangestring = "_$chr"."\:minus(".$rangestring.")";
			$seqstrand = "-";
		} elsif ($strand > 0) {
			$rangestring = "_$chr"."\:plus(".$rangestring.")";	
			$seqstrand = "+";
		} else {
		print "I'm confused strand is $strand 2\n";
		sleep 2;
		}

		
		if ($fast eq "q") {	$readid = "@"."$id"."_"."$count2"; };
		if ($fast eq "a") {	$readid = ">"."$id"."_"."$count2"; };
		$readid = $readid.$rangestring."#0/2";
		#print "readid here is $readid\n";
		#sleep 1;
		$SIMREAD2->printf($readid); 
		$SIMREAD2->printf("\n");
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#!  Printing to BED in zero-based system
		$readid =~ s/\>//;
		$readid =~ s/\@//;
				
		if (@junctions) { # Junctions array is NOT empty

			print { $READBED } "chr$chr\t",$readcoords[0] - 1,"\t$readcoords[$#readcoords]\t$readid\t"."1"."\t$seqstrand\t",$readcoords[0] - 1,"\t$readcoords[$#readcoords]\t255,0,0\t";
			$numjuncs = @junctions;
			print { $READBED } $numjuncs + 1,"\t";
			$tempstart = $readcoords[0];
			my $sizetemp;
			my $coordtemp;
			$blockcoordstring = "0";
			$blocksizestring = "";
			#print "\n** $tempstart to $readcoords[$#readcoords]\t$blockcoordstring\t$blocksizestring **\n";
			foreach (@junctions) {
				#$j += 1;
				$tempjunc = $_;
				$junckey = "chr$chr"."_".$_."_$direction";
				if ($tempjunc =~ /([0-9]+)\.\.([0-9]+)/) {$gstart = $1};
	
				$offset = $gstart - @readcoords[0] + 1;

				if (exists($junchash{$junckey})) {
					$junchash{$junckey}[0] += 1;
					$prev = pop @{$junchash{$junckey}};
					#print "prev for $junckey is $prev\n";
					
					push(@{$junchash{$junckey}}, "$prev\,".$offset);

				} else {
					$junchash{$junckey}[0] += 1;
					#print "added ", $offset, " to hash with key $junckey\n";
					push(@{$junchash{$junckey}},$offset);
					#sleep 1;
				}
				
				if ($verbose eq "T") {print $_."\n";}
					if ($_ = /([0-9]+)\.\.([0-9]+)/) {
						#print "1 is $1 and 2 is $2\n"; 
						$sizetemp = $1 - $tempstart + 1;
						$blocksizestring = $blocksizestring.$sizetemp."\,";
						$coordtemp = $2 - $readcoords[0];
						$blockcoordstring = $blockcoordstring."\,".$coordtemp;
						if ($verbose eq "T") {print "\n$tempstart\t$blockcoordstring\t$blocksizestring";}
						$tempstart = $2;
					}
				
			}
			$sizetemp = $readcoords[$#readcoords] - $tempstart + 1;
			$blocksizestring = $blocksizestring.$sizetemp;
			#$blockcoordstring = $blockcoordstring.",".$readcoords[$#readcoords] - $tempend;
			if ($verbose eq "T") {print "\n Final: $blockcoordstring\t$blocksizestring";}
			
			print { $READBED } "$blocksizestring\t$blockcoordstring\n";

		} else { # If there are no junctions, the line is simple
			print { $READBED } "chr$chr\t",$readcoords[0] - 1,"\t$readcoords[$#readcoords]\t$readid\t"."1"."\t$seqstrand\t",$readcoords[0] - 1,"\t$readcoords[$#readcoords]\t255,0,0\t1\t75\t0\n";



		}

	#-------------------------		
		# Get read sequence
		$read = substr($sequence,$randstart, $mer);
		# Reverse the read 
		$revcom = reverse $read;
		$revcom =~ tr/ACGTacgtnN/TGCATGCANN/;
		$read = $revcom;
		if (length($read) != $mer) {
			print "Read length is wrong!\n";
			print $read, "\n";
			sleep 2;
		}
		
		#---------------------------------

				
		if ($mmtf eq "T") {
		$SIMREAD2->printf(mmgen($read));
		} elsif ($mmtf eq "F") {
		$SIMREAD2->printf($read);
		} else {
		print "I dont get it mmtf is $mmtf 2\n";
		sleep 3;
		}
		$SIMREAD2->printf("\n");

		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# Now make a quality string
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		# Now make a quality string
		#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		if ($fast eq "q") {	
			#$readid = "+"."$id"."_"."$count2";
			$readid = "+".$readid;
			$SIMREAD2->printf($readid);
			$SIMREAD2->printf("\n");			
			if ($qualspeed eq "T") {			
				$qualstring = @randquals[int(rand(999))];
				#print "qualstring is $qualstring\n";
				print { $SIMREAD2 } $qualstring;
				$SIMREAD2->printf("\n");		
			} else {
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
	

}


# Sample Bowtie map file

#HWUSI-EASXXX:4:1:0:1368#0/1 Ê Ê + Ê Ê Ê 2R Ê Ê Ê16938356 NAGGGCTCCTCGCCGTACATGTCGGCCAATTTGCGGAAGCGCGGTCCAAAGTTGGATCGGTAGTCGAAGTTGAGAT Ê 0 Ê Ê Ê 0:G>N,56:C>T,57:A>C 
#HWUSI-EASXXX:4:1:0:1368#0/2 Ê Ê - Ê Ê Ê 2R Ê Ê Ê16938454 ÊCGGACGAGAGAGGCTGCCATCGGAGTTGCCGTCACCTTCGTACGCGTAATGCCGCACATCGTCCACGGTTGTGGC Ê Ê0 Ê Ê Ê 69:G>A,70:T>G,74:A>G
sub generate_pe_reads_old {
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
			#print "adding $fbtr to rtxfile hash\n";
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
			#print "putting $txid in seqhash, $fbtr was found $previd\n";
			#print "length of seq is ". length($tx);
			#sleep 1;
			$seqhash{$txid}[0] = $previd;
			$seqhash{$txid}[1] = $tx;
			$seqhash{$txid}[2] = $prevheader;
			$tx = "";
			
			} elsif (($type eq "intron") && (! exists($rtxs{$fbtr}))) {
				$tx = ""  # This is an important change.. hope it is OK
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
			print "putting $txid in seqhash\n";
			print length($tx);
			#sleep 1;

	$seqhash{$txid}[0] = $previd;
	$seqhash{$txid}[1] = $tx;
	$seqhash{$txid}[2] = $prevheader;
	#print "added transcript # $txid to hash\n";
}

print "Loaded $txid transcripts to hash\n";
print "Working on generating reads\n";




#xxxxxxxxxxxxxxxxxxxxxxxxxx
# Quality score generation speedup
# Create 1000 qualscores by model, then randomly pick from them.

if (($qualspeed eq "T") & ($fast eq "q")) {
	
	my %tempqualhash;
	my $position;
	my $maxprob;
	my $j;
	my $qualstring;

	
	my $i;
	my $randqual;
	#for ($count = 10; $count >= 1; $count--)
	for ($i = 0; $i < 1000; ++$i) {

		$qualstring = "";
		for ($position = 0; $position < $bp; ++$position) {		

			$maxprob = 0;
			foreach $j (keys %qualpwm) {	

				$tempqualhash{$j} = $qualpwm{$j}[$position];
				if ($qualpwm{$j}[$position] > $maxprob) {$maxprob = $qualpwm{$j}[$position];}
			}	
			$randqual = weighted_rand(%tempqualhash);
			$qualstring = $qualstring.$randqual;	
		}	
			
			
		push(@randquals,$qualstring);	
		#print length($qualstring), "\t$i\n";
	}
	
	
	if ($verbose eq "T") {print "There were ". @randquals . "random qualities generated\n"};
#sleep 3;
if ($verbose eq "T") {
	foreach (@randquals) {print "$_ \n";}
}#sleep 3;
} 


#xxxxxxxxxxxxxxxxxxxxxxxxxxx


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

my $newcov;

my $innerd = 0;

my $lengthcutoff;

my $mincovtx = 0;


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

my @randcovs;
my $headerstrand;
my $headercoords;
my $strand = 0;

my $transcriptid;

foreach $k (sort keys %seqhash) {
	if ($verbose eq "T") {print "a\n";}	
	$kcount += 1;
	$txstring = $seqhash{$k}[1];
	$hstring = $seqhash{$k}[0];
	
	#"loc=2R:complement(4999006..5000837,5000900..5001458);"
	#print $hstring;
	 if ($hstring =~ /(FBtr[0-9]+).+loc\=([a-zA-Z0-9\_]+)\:([a-zA-Z]*)\(*([0-9\.\,]+).+(FBgn[0-9]{7})/) {

		 $fbtr = $1;
		 $transcriptid = $1;
		 $chr = $2;
		 $headerstrand = $3;
		 $headercoords = $4;
		 $FB = $5;  # This is the parent gene

	} elsif ($hstring =~ /\>([A-Za-z0-9\_\:]+).+loc\=([a-zA-Z0-9\_]+)\:([a-zA-Z]*)\(*([0-9\.\,]+).+(FBgn[0-9]{7}).+(FBtr[0-9]{7})/) { 
		#print "* hstring is $hstring\n";	
		 $fbtr = $1;	
		 $chr = $2;
		 $headerstrand = $3;
		 $headercoords = $4;
		 $FB = $5;  # This is the parent gene
		 $transcriptid = $6;
		# print "headerstrand is $headerstrand\n";
		#print "headerccords is $headercoords\n";
		#print "fb is $FB\n";
		#print "chr is $chr\n";
		#print "fbtrr is $fbtr\n";
	} elsif ($hstring =~ /\>([A-Za-z0-9\_\:]+).+loc\=([a-zA-Z0-9\_]+)\:([a-zA-Z]*)\(*([0-9\.\,]+)/) { 
		#print "* hstring is $hstring\n";	
		 $fbtr = $1;	
		 $chr = $2;
		 $headerstrand = $3;
		 $headercoords = $4;
		 if ($fbtr =~ /(FBgn[0-9]{7})/) {$FB = $&};  # This is the parent gene
		 #print "headerstrand is $headerstrand\n";
		#print "headerccords is $headercoords\n";
		#print "fb is $FB\n";
		#print "chr is $chr\n";
		#print "fbtrr is $fbtr\n";
		#sleep $1
	} else {
		$hstring =~ s/\>//;
		$hstring =~ s/\n//;
		$fbtr = $hstring;	
	
	}
	
	#print "$fbtr\n";
	#print "$chr\n";
	#print "$headerstrand\n";
	#print "$headercoords\n";
	#print $headercoords, "\n";
	@coordstring = coordexplode($headercoords);
	#print "headercoords was $headercoords\n";
	#print "before generating reads, coord 0 is $coordstring[0]\n";
	#sleep 2;
	#foreach (@coordstring) {
 	#  print "$_\n";
	#}
	#sleep 2;
	my $gffstrand;
	my $direction;
	
	if ($type eq "intergenic") {$headerstrand = "join"}
	
	if ($headerstrand eq "join") {
		$strand = 1;
		$direction = "f";
	} elsif ($headerstrand eq "complement") {
		$strand = -1;
		$direction = "r";
	} else {
		#$gffstrand = $gffhash{$transcriptid}[3];
		$gffstrand = $gffhash{$FB}[3];
		#print "Couldn't determine strand from header line for $fbtr, got $gffstrand from hash\n";	
		#print $hstring, "\n";
		if ($gffstrand eq "+") {
			$strand = 1;
			$direction = "f";
			} elsif ($gffstrand eq "-") {
			$strand = -1;
			$direction = "r";
			} else { 
			print "Couldn't determine strand for $fbtr, so picked +\n";
			print "transcript is was $transcriptid\n";
			if (exists($gffhash{$transcriptid})) {print "hash entry exists $gffhash{$transcriptid}[0]\t $gffhash{$transcriptid}[3]\n";}
			$strand = 1;
			$direction = "f";

			}
	}
	
	
	
	unless ($type eq "tx") {
			#if ($hstring =~ /(FBgn[0-9]{7})/) {
	 
			 #$fbtr = $type."_".$kcount;
			
			 
			 #$FB = $&;  # This is the parent gene

			#}

	
	}
	
	
	if ($varcov eq "var") {
	
		$mean = $originalcov;
		#$sdev = 1;
		$newcov = (gaussian_rand() * $sdev) + $mean;
		# Sanity check, if rand coverage is negative or too low
		if ($newcov < 0.5) {
			$newcov = 0.5;
			$mincovtx += 1;
		}
		#printf("Random coverage is %.2f\n", $newcov);
		push (@randcovs, $newcov);
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


	# Convert to revcom if on minus strand.  This
	# is important to getting the coordinates right
	my $revcom;
		if ($strand < 0) { 
			$revcom = reverse $txstring;
			$revcom =~ tr/ACGTacgtnN/TGCATGCANN/;
			$txstring = $revcom;
			$strand = 1;
		}


	
	if ($txlength > $lengthcutoff) {
		$ANSWERS->printf("$fbtr\t$FB\t$chr\t$txlength\t$numreads\t$RPK\t$cov\n");  # Print my "answers" table
		if ($verbose eq "T") {print "c\n";}	
	
			if ($ends == 1) {
			#print "*b\t$k\n";
			if ($verbose eq "T") {print "Starting 'generate reads' sub $fbtr $numreads strand: $strand $direction $txstring $chr\n";}	
			generate_reads($txstring, $fbtr, $bp, $numreads, $fast, $chr, $strand, $direction, $mmtf, $qualspeed);
			} else {
			#print "*b\t$k\n";
			generate_pe_reads($txstring, $fbtr, $bp, $numreads, $fast, $chr, $strand, $direction, $mmtf, $qualspeed);

			#generate_pe_reads($txstring, $fbtr, $bp, $numreads, $fast);
			}
			
			#generate_qualities($txstring, substr($hstring,0,length($hstring) - 1), $bp, $numreads);
	} else {
		$ANSWERS->printf("$fbtr\t$FB\t$chr\t$txlength\tskipped\ttoo short\n");  # Print my "answers" table
	}
	
	
}		

print "\n";

if ($varcov eq "var") {
print "Here is a histogram of random coverages (rounded)\n";

histogram(1,@randcovs);

}

if ($mincovtx > 0) {
print "\n---------------------\n";
print "Warning:  For $mincovtx transcripts, the random coverage was less than 0.5, and was set to equal 0.5\n";
}


#== Below prints out the models
=pod
foreach $j (keys %mmtypes) {
	$refbase = substr($j,0,1);
	$readbase = substr($j,2,1);
	$mmbase{$refbase}{$readbase} = $mmtypes{$j};
	$mmbasetotals{$refbase} += $mmtypes{$j};	
}
=cut


=pod
for $k ( keys %wighash ) {
    print "k is $k\n";
    sleep 2;
    print "variableStep chrom=chr$k\n";
    for $i ( 0 .. $#{ $wighash{$k} } ) {
        
        if ($wighash{$k}[$i] < 1) {$wighash{$k}[$i] = 0}
        print "$i\t$wighash{$k}[$i]\n";
    }
    #print "\n";
}
=cut


#for $family ( keys %HoA ) {
#    print "$family: ";
#    for $i ( 0 .. $#{ $HoA{$family} } ) {
#        print " $i = $HoA{$family}[$i]";
#    }
#    print "\n";
#}



my $k;
my $wigcoord;

# Print Wig from bedhash.  Note that WIG files are +1 based, BED files are zero-based.
print "Generating wiggle file\n";
print { $WIG } "track type=wiggle_0\n";
#print "track type=bedGraph name='simulated read coverage'\n";
for $k ( keys %wighash ) {
    print { $WIG }  "variableStep chrom=chr$k\n";
    for $i ( 0 .. $#{ $wighash{$k} } ) {  #Bedhash is a hash of arrays, with chr as key
        
        $wigcoord = $wighash{$k}[$i] + 1;
        if ($wighash{$k}[$i] >= 1) {
        
        print { $WIG }  "$i\t$wighash{$k}[$i]\n";
		}

    }
    #print "\n";
}

# Here is junction hash
=pod
for $k ( keys %junchash ) {
   print "$k\t$junchash{$k}\n";
}


for $k ( keys %junchash ) {
    print "$k: ";
    for $i ( 0 .. $#{ $junchash{$i} } ) {
        print " $i = $junchash{$i}[$i]";
    }
    print "\n";
}



foreach my $k (keys %junchash) {
    print "The members of $k are\n";
    foreach (@{$junchash{$k}}) {
        print "\t$_\n";
    }
}

=cut
if ($type eq "tx") {
	print { $JUNCTAB } "junction(chr_coord_direct)\tdensity\toffsets\n";
	foreach my $k (keys %junchash) {
		print { $JUNCTAB } "$k\t";
		foreach (@{$junchash{$k}}) {
			print { $JUNCTAB } "$_\t";
		}
		print { $JUNCTAB } "\n";
	}
}


# Below is a BED-like format

=pod
print "track type=bedGraph name='simulated read coverage'\n";
for $k ( keys %wighash ) {
#    print "k is $k\n";
#    sleep 2;
#    print "variableStep chrom=chr$k\n";
    for $i ( 0 .. $#{ $wighash{$k} } ) {
        if ($wighash{$k}[$i] >= 1) {
        
        #if ($wighash{$k}[$i] < 1) {$wighash{$k}[$i] = 0}
        #print "$i\t$wighash{$k}[$i]\n";
        print "chr$k\t$i\t", $i + 1, "\t$wighash{$k}[$i]\n";
		}

    }
    #print "\n";
}

=cut