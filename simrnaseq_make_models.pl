#!/usr/bin/perl
#
# simrnaseq_make_models.pl

use strict;
use FileHandle;
use Data::Dumper;
use Getopt::Long;


my $FILENAME;
my $ends;

my @f;
my $line;

my $bp;
my $modeldir;

# Specify cmdline options and process command line...
my $options_okay = GetOptions (
    # Application-specific options...
    'i=s'     => \$FILENAME,     # --in option expects a string
    'ends=i'    => \$ends,    # --ends option expects an integer
    'length=i' => \$bp,     # --length option expects an integer
    'd=s'		=> \$modeldir # --models directory is a string

);


usage () if (! $options_okay);
	
sub usage
{
print "Unknown option: @_\n" if ( @_ );
print "usage: program [-i INPUT MAP FILE] [-e 1=single ends 2=paired ends] [-length read length] [-d model directory]\n";
exit;
}

print "my mapfile is $FILENAME\n";
print "my number of ends is $ends\n";

unless (($ends < 3) && ($ends > 0)) {
print "ends must be 1 or 2\n";
exit;
}

# Make directory to put models in

mkdir($modeldir, 0777) || print $!, ", data in your output directory will be overwritten\n";



my $MM_pos_PROBS;
my $MM_pos_counts;
my $MM_type_PROBS;
my $QUAL_PWM;
my $QUAL_counts;

my $LOGFILE = new FileHandle ">$modeldir/logfile.txt" or die "can't open $modeldir/logfile.txt"; 

my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
print { $LOGFILE } "Model building started on:\n";
$LOGFILE -> printf("%4d-%02d-%02d %02d:%02d:%02d\n",$year+1900,$mon+1,$mday,$hour,$min,$sec);
print { $LOGFILE } "Options for this run:\n";
print { $LOGFILE } "filename = $FILENAME\n";
print { $LOGFILE } "ends = $ends\n";
print { $LOGFILE } "length = $bp\n\n";

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# For single end reads
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if ($ends == 1) {
$MM_pos_PROBS = new FileHandle ">$modeldir/mmprobs.txt" or die "can't open $modeldir/mmprobs.txt"; #OUTPUT: mm probabilities
$MM_pos_counts = new FileHandle ">$modeldir/mmcounts.txt" or die "can't open $modeldir/mmprobs.txt"; #OUTPUT: mm probabilities
$MM_type_PROBS = new FileHandle ">$modeldir/mmtypes.txt" or die "can't open $modeldir/mmtypes.txt"; #OUTPUT: mm type frequencies
$QUAL_PWM = new FileHandle ">$modeldir/qualitiespwm.txt" or die "can't open $modeldir/qualitiespwm.txt"; #OUTPUT: PWM of qual.scores
$QUAL_counts = new FileHandle ">$modeldir/qualitiescounts.txt" or die "can't open $modeldir/qualitiespwm.txt"; #OUTPUT: PWM of qual.scores
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# For paired end reads
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
my $MM_pos_PROBS1;
my $MM_pos_counts1 ;
my $MM_type_PROBS1;
my $QUAL_PWM1;
my $QUAL_counts1; 
my $MM_pos_PROBS2;
my $MM_pos_counts2;
my $MM_type_PROBS2;
my $QUAL_PWM2 ;
my $QUAL_counts2 ;
my $INNERD ;

if ($ends == 2) {
$MM_pos_PROBS1 = new FileHandle ">$modeldir/mmprobs_1.txt" or die "can't open mmprobs.txt"; #OUTPUT: mm probabilities
$MM_pos_counts1 = new FileHandle ">$modeldir/mmcounts_1.txt" or die "can't open mmprobs.txt"; #OUTPUT: mm probabilities
$MM_type_PROBS1 = new FileHandle ">$modeldir/mmtypes_1.txt" or die "can't open mmtypes.txt"; #OUTPUT: mm type frequencies
$QUAL_PWM1 = new FileHandle ">$modeldir/qualitiespwm_1.txt" or die "can't open qualitiespwm.txt"; #OUTPUT: PWM of qual.scores
$QUAL_counts1 = new FileHandle ">$modeldir/qualitiescounts_1.txt" or die "can't open qualitiespwm.txt"; #OUTPUT: PWM of qual.scores

$MM_pos_PROBS2 = new FileHandle ">$modeldir/mmprobs_2.txt" or die "can't open $modeldir/mmprobs.txt"; #OUTPUT: mm probabilities
$MM_pos_counts2 = new FileHandle ">$modeldir/mmcounts_2.txt" or die "can't open $modeldir/mmprobs.txt"; #OUTPUT: mm probabilities
$MM_type_PROBS2 = new FileHandle ">$modeldir/mmtypes_2.txt" or die "can't open $modeldir/mmtypes.txt"; #OUTPUT: mm type frequencies
$QUAL_PWM2 = new FileHandle ">$modeldir/qualitiespwm_2.txt" or die "can't open $modeldir/qualitiespwm.txt"; #OUTPUT: PWM of qual.scores
$QUAL_counts2 = new FileHandle ">$modeldir/qualitiescounts_2.txt" or die "can't open $modeldir/qualitiespwm.txt"; #OUTPUT: PWM of qual.scores

$INNERD = new FileHandle ">$modeldir/mateinnerdists.txt" or die "can't open output file"; #OUTPUT: PWM of qual.scores

# Also output combined models (both ends put in same model)
$MM_pos_PROBS = new FileHandle ">$modeldir/mmprobs.txt" or die "can't open $modeldir/mmprobs.txt"; #OUTPUT: mm probabilities
$MM_pos_counts = new FileHandle ">$modeldir/mmcounts.txt" or die "can't open $modeldir/mmprobs.txt"; #OUTPUT: mm probabilities
$MM_type_PROBS = new FileHandle ">$modeldir/mmtypes.txt" or die "can't open $modeldir/mmtypes.txt"; #OUTPUT: mm type frequencies
$QUAL_PWM = new FileHandle ">$modeldir/qualitiespwm.txt" or die "can't open $modeldir/qualitiespwm.txt"; #OUTPUT: PWM of qual.scores
$QUAL_counts = new FileHandle ">$modeldir/qualitiescounts.txt" or die "can't open $modeldir/qualitiespwm.txt"; #OUTPUT: PWM of qual.scores



}

################################
################################
# Load in map file
################################
################################

my $read_name;
my $bowtie_read_count=0;
my $start;
my $end;
my $strand;
my $readlength;

my %mappings;
my $mapseq;

my $totalreads = 0;

my %uniqueseq;
   
my @mmprobs;
my @mmprobs1; ## For paired ends;
my @mmprobs2; ## For paired ends;

my @mms;
my $mm;
my $mmstring;

my $feature;

my $mmpos;

my %mmtype;
my %mmtype1;  ## For paired ends;
my %mmtype2;  ## For paired ends;

my %pwm;
my %pwm1;  ## For paired ends;
my %pwm2;  ## For paired ends;

my $i = 0;
while ($i < $bp) {  # If read is longer than the length entered, don't add all positions to model
	@mmprobs[$i] = 0;
	$i += 1;
}

my $qualstring;
my $sitepos;
my $mate;
my $read_name_save;
my @mate1starts;
my @mate2starts;

my $first = 0;

open(BOWTIE_OUT, "$FILENAME") or die("can't open $FILENAME"); 

while ($line = <BOWTIE_OUT>){       # read in each line of map file
	#print $line;
	if ($first < 1) {
	print { $LOGFILE } "Here is the first line of the map file\n";
	print { $LOGFILE } "$line\n";
	}
	$first += 1;
	@f = split /\s+/,$line;
	$read_name = $f[0];
	$read_name_save = $read_name;
	$mate = chop $read_name;
	$strand = $f[1];
	$feature = $f[2];
	$start = $f[3];
	$mapseq = $f[4];
	$readlength = length($f[4]);
	$end = $start+$readlength;
	$qualstring = $f[5];

	$mmstring = $f[7];
	@mms = split(',', $mmstring);

	# Truncate the read and qual strings to equal user-input length
	$qualstring = substr($qualstring,0,$bp);
	$mapseq = substr($mapseq,0,$bp);
	
	if ($strand eq "-") {
		$qualstring = reverse $qualstring;
	} elsif ($strand ne "+") {
		print "Can't determine strand\n";
		print "strand is $strand\n";
		exit;
	}


      
	foreach $mm (@mms) {   # Instantiate mismatch array
	   if ($mm =~ /([0-9]+)\:([A-Z]+)\>([A-Z]+)/) {
		   $mmpos = $1;
		   unless ($mmpos > ($bp - 1)) {
		   		$mmtype{"$2>$3"} += 1;
		   		@mmprobs[$mmpos] += 1;
		   }
	   }
	}
   
   
   if ($ends == 2) {  ## If paired ends, also do a mm tab by mate
		if ($mate eq 1) {
			push (@mate1starts, $start);
			foreach $mm (@mms) {   # Instantiate mismatch array
		   		if ($mm =~ /([0-9]+)\:([A-Z]+)\>([A-Z]+)/) {
			   		$mmpos = $1;
			   		unless ($mmpos > ($bp - 1)) {
						$mmtype1{"$2>$3"} += 1;
						@mmprobs1[$mmpos] += 1;
			  		}
		  		 }
	   		}	   
	   		while ($qualstring =~ /\S/gi) {
				$sitepos = pos($qualstring) - 1;
				++$pwm1
				{$&}[$sitepos];		
			}
		} elsif ($mate eq 2) {
			push (@mate2starts, $start);
	   		foreach $mm (@mms) {   # Instantiate mismatch array
		   		if ($mm =~ /([0-9]+)\:([A-Z]+)\>([A-Z]+)/) {
			  		$mmpos = $1;
			  		unless ($mmpos > ($bp - 1)) {
						$mmtype2{"$2>$3"} += 1;
						@mmprobs2[$mmpos] += 1;
			   		}
		   		}
	   		}
	   		while ($qualstring =~ /\S/gi) {
				$sitepos = pos($qualstring) - 1;
				++$pwm2{$&}[$sitepos];		
			}	   
		} else {
			print "Could not identify which mate this is\n";
			exit;
		} 
   }
   
  $totalreads += 1;

	while ($qualstring =~ /\S/gi) {
		$sitepos = pos($qualstring) - 1;
		++$pwm{$&}[$sitepos];		
	}
}

close BOWTIE_OUT;  

if ($readlength ne $bp) {
	print "NOTE!  Read length $readlength is different than input length $bp bp\n";
	print "Truncating the reads to $bp\n";
}


my $k;
my $j;
my $qualprob;

print "There are $totalreads reads mapping\n";
print { $LOGFILE } "$totalreads mapped reads were used in this model\n";

my $prob;
my $mmprob;
my $i = 0;


my $temp = @mmprobs;
print "mmprobs is $temp\n";

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# print the probability models
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# In this version, print the combined model if ends is 1 OR 
#if ($ends == 1) {

foreach $mmprob(@mmprobs) {
	$prob =  $mmprob/$totalreads*100;
	$MM_pos_PROBS->printf("$i\t$prob\n");
	$MM_pos_counts->printf("$i\t$mmprob\n");
	$i += 1;
}


foreach $j (keys %mmtype) {
	$MM_type_PROBS->printf("$j: $mmtype{$j}\n");
}
	
$i = 0;

foreach $j (keys %pwm) {
	$QUAL_PWM->printf("$j\t");
	$QUAL_counts->printf("$j\t");	
	while ($i < $bp) { 
		if ($pwm{$j}[$i] >= 1) {
			$qualprob = $pwm{$j}[$i] / $totalreads;
			$qualprob = sprintf("%1.5f",$qualprob);
			$QUAL_PWM->printf("$qualprob\t");
			$QUAL_counts->printf("$pwm{$j}[$i]\t");	
		} else {
			$QUAL_PWM->printf("0\t");
			$QUAL_counts->printf("0\t");	
		}
		$i += 1;		
	}
	$QUAL_PWM->printf("\n");
	$QUAL_counts->printf("\n");
	$i = 0;
}

#}

if ($ends == 2) {

	$i = 0;
	#######
	foreach $mmprob(@mmprobs1) {
		$prob =  $mmprob/$totalreads/2*100;
		$MM_pos_PROBS1->printf("$i\t$prob\n");
		$MM_pos_counts1->printf("$i\t$mmprob\n");
		$i += 1;
	}
	
	
	foreach $j (keys %mmtype1) {
		$MM_type_PROBS1->printf("$j: $mmtype1{$j}\n");
	}
	
	$i = 0;
	
	foreach $j (keys %pwm1) {
	$QUAL_PWM1->printf("$j\t");
	$QUAL_counts1->printf("$j\t");
	
		while ($i < $bp) {
			if ($pwm1{$j}[$i] >= 1) {
				$qualprob = $pwm1{$j}[$i] / ($totalreads / 2);
				$qualprob = sprintf("%1.5f",$qualprob);
				$QUAL_PWM1->printf("$qualprob\t");
				$QUAL_counts1->printf("$pwm1{$j}[$i]\t");
			} else {
				$QUAL_PWM1->printf("0\t");
				$QUAL_counts1->printf("0\t");
			}
			$i += 1;		
		}
		$QUAL_PWM1->printf("\n");	
		$QUAL_counts1->printf("\n");		
		$i = 0;
	}
	########
	$i = 0;
	foreach $mmprob(@mmprobs2) {
		$prob =  $mmprob/$totalreads/2*100;
		$MM_pos_PROBS2->printf("$i\t$prob\n");
		$MM_pos_counts2->printf("$i\t$mmprob\n");
		$i += 1;
	}
	
	
	foreach $j (keys %mmtype2) {
		$MM_type_PROBS2->printf("$j: $mmtype2{$j}\n");
	}
	
	$i = 0;
	
	foreach $j (keys %pwm2) {
	$QUAL_PWM2->printf("$j\t");
	$QUAL_counts2->printf("$j\t");
		while ($i < $bp) {
			if ($pwm2{$j}[$i] >= 1) {
				$qualprob = $pwm2{$j}[$i] / ($totalreads / 2);
				$qualprob = sprintf("%1.5f",$qualprob);
				$QUAL_PWM2->printf("$qualprob\t");
				$QUAL_counts2->printf("$pwm2{$j}[$i]\t");	
			} else {
				$QUAL_PWM2->printf("0\t");
				$QUAL_counts2->printf("0\t");
			}
			$i += 1;		
		}
		$QUAL_PWM2->printf("\n");
		$QUAL_counts2->printf("\n");
		$i = 0;
	}
	########
	
	# Now get inner distances:
	
	my %innerdistances;
	my $innerd;
	
	my $arraysize = scalar @mate1starts;
	#print "There are $arraysize starts\n";
	
	my $pos;
	
	for ($pos=0 ; $pos < $arraysize ; $pos += 1) {
	$innerd = abs(@mate1starts[$pos] - @mate2starts[$pos]) - $readlength;
	#print "abs of @mate1starts[$pos] - @mate2starts[$pos] - $readlength is $innerd\n";
	$innerdistances{$innerd} += 1
	}


	my $cnt = 0;
	
	foreach $k (sort keys(%innerdistances)) {

	#print "$k\t$innerdistances{$k}\n";
	$INNERD->printf("$k\t$innerdistances{$k}\n");
	$cnt += $innerdistances{$k};

	}

	#print "There are $cnt total innerdistances\n";
	
}


($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
print { $LOGFILE } "Model building finished at:\n";
$LOGFILE -> printf("%4d-%02d-%02d %02d:%02d:%02d\n",$year+1900,$mon+1,$mday,$hour,$min,$sec);





