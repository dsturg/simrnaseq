#!/usr/bin/perl

# PROGRAM  : simrnaseq_create_random_seqlibraries.pl
# AUTHOR   : Dave
# Last Revision Date: Dec 3, 2009
# Last Revised By: Dave

# Usage: perl simrnaseq_create_random_seqlibraries.pl <ref. transcript file> <number> <options> 
# Example: perl simrnaseq_create_random_seqlibraries.pl ~/data/indexes/dmel-transcript-r5.20.fa 1000 2


use strict;
use FileHandle;

my $filename = $ARGV[0];
chomp $filename;

my $size = $ARGV[1];

my $option = 0;
$option = $ARGV[2];

die "Don't forget to enter how many transcripts you want\n" unless ($size > 0) ;

unless ($option > 0) {
print "Choose 0 pick any random transcript\n";
print "Choose 1 to reject mitochondrial genome transcripts\n";
print "Choose 2 to reject mitochondrial and chromosome U genes\n";
print "Choose 3 to reject mitochondrial and chromosome U genes AND to limit to one isoform per gene\n";
print "NOTE: random selection is always performed WITHOUT replacement\n";
die "** Don't forget to enter your options\n" ;

}

unless (($option == 0) | ($option == 1) | ($option == 2) | ($option == 3)) {
print "Choose 0 pick any random transcript\n";
print "Choose 1 to reject mitochondrial genome transcripts\n";
print "Choose 2 to reject mitochondrial and chromosome U genes\n";
print "Choose 3 to reject mitochondrial and chromosome U genes AND to limit to one isoform per gene\n";

die "** Invalid option specified\n" ;

}


my $RANDCDS = new FileHandle ">rand_transcripts.fa" or die "can't open rand_transcripts.fa"; #OUTPUT: random CDS
my $BED = new FileHandle ">rand_transcripts.bed" or die "can't open randcds.bed"; #OUTPUT: random CDS
my $JUNCBED = new FileHandle ">rand_transcripts.juncs.bed" or die "can't open randtranscripts.juncs.bed"; #OUTPUT: random transcripts

sub print_sequence {
	my($sequence, $length) = @_;
	for (my $pos=0 ; $pos < length($sequence) ; $pos += $length) {
		$RANDCDS->printf(substr($sequence, $pos, $length));
		$RANDCDS->printf("\n");
	}
}

sub median { 
@_ == 1 or die ('Sub usage: $median = median(\@array);'); 
my ($array_ref) = @_; 
my $count = scalar @$array_ref; 
# Sort a COPY of the array, leaving the original untouched 
my @array = sort { $a <=> $b } @$array_ref; 
if ($count % 2) { 
return $array[int($count/2)]; 
} else { 
return ($array[$count/2] + $array[$count/2 - 1]) / 2; 
} 
} 

my $median;
my $k;


###############################
# The part below is for reading
# in a fasta file and instantiating
# a hash
###############################


open (SEQFILE, $filename);

my $seq;
my %seqhash;
my $txid = 0;
my $tx;
my $previd;
my $header;
my $prevheader = "";


while ($seq = <SEQFILE>) {   ## First load all transcripts into a hash
	if ($seq =~ /^>/) {
		$header = $seq;
		if ($txid > 0) {
			$seqhash{$txid}[0] = $previd;
			$seqhash{$txid}[1] = $tx;
			$seqhash{$txid}[2] = $prevheader;
			$tx = "";

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

$seqhash{$txid}[0] = $previd;
$seqhash{$txid}[1] = $tx;
$seqhash{$txid}[2] = $prevheader;
print "\n ..putting $txid transcripts in hash\n";
###############################
###############################
###############################


##########################
##########################

#~~~~~~~~~~~~~~~~~~~~~~~~~
# For parsing header line:
#~~~~~~~~~~~~~~~~~~~~~~~~~

my %FBseen;
my $FB;
my $coords;
my @coordstring;
my $chr;
my $strand;

#~~~~~~~~~~~~~~~~~~~~~~~~~

print "\n";
print "NOTE:  This program expects a fasta file (from flybase) with headers like this:\n";
print ">FBtr0071764 type=mRNA; loc=2R:join(18024938..18025756,18039159..18039200,18050410..18051199,18052282..18052494,18056749..18058222,18058283..18059490,18059587..18059757,18059821..18059938,18060002..18060346); ID=FBtr0071764; name=a-RB; dbxref=FlyBase_Annotation_IDs:CG6741-RB,FlyBase:FBtr0071764,REFSEQ:NM_079902; score=7; score_text=Moderately Supported; MD5=8ef39c97ee0a0abc31900fefe3fbce8f; length=5180; parent=FBgn0000008; release=r5.20; species=Dmel; \n";
print "\n";


# Sample header line
#>abo-PA type=CDS; loc=2L:complement(10973536..10974210,10974268..10975212); name=abo-RA; dbxref=FlyBase:FBpp0079757,FlyBase_Annotation_IDs:CG6093-PA,REFSEQ:NP_524784,GB_protein:AAF53040; MD5=5d20f439f7478b3302f76a317d7fe9c9; length=1620; parent=FBgn0000018,FBtr0080168; release=r5.17; species=Dmel; 
# Sample header line for transcripts file																																																			  FBgn0000008
#>FBtr0071764 type=mRNA; loc=2R:join(18024938..18025756,18039159..18039200,18050410..18051199,18052282..18052494,18056749..18058222,18058283..18059490,18059587..18059757,18059821..18059938,18060002..18060346); ID=FBtr0071764; name=a-RB; dbxref=FlyBase_Annotation_IDs:CG6741-RB,FlyBase:FBtr0071764,REFSEQ:NM_079902; score=7; score_text=Moderately Supported; MD5=8ef39c97ee0a0abc31900fefe3fbce8f; length=5180; parent=FBgn0000008; release=r5.18; species=Dmel; 

# Sample for a single-exon gene:
#>FBtr0080316 type=mRNA; loc=2L:complement(11988735..11989861); ID=FBtr0080316; name=zuc-RA; dbxref=FlyBase:FBtr0080316,FlyBase_Annotation_IDs:CG12314-RA,REFSEQ:NM_135686; score=11; score_text=Strongly Supported; MD5=0e270a2540dd9faa1dbb971b

#~~~~~~~~~~~~~~~~~~~~~~~~~
# For choosing random transcripts
#~~~~~~~~~~~~~~~~~~~~~~~~~

my $randnum;
my $count;

my $keptcount = 0;
my $rejectedcount = 0;
my $totalcount = 0;

my $exoncount;
my @exoncount;

my $reject;
#~~~~~~~~~~~~~~~~~~~~~~~~~

#### For the JUNCBED file:
my $leftedge;
my $rightedge;
my $intronstart;
my $intronend;
my $junccount = 0;
my $juncid;

my $singlecheck;

my $fbtr;

my %seentx;

my $mitocount = 0;
my $Ucount = 0;

my $isoformrej = 0;

for ($count=1; $count < $size+1; $count++) {
	
	$totalcount += 1;
	$randnum = int(rand($txid));
	 
	 # This is where I ensure it is random WITHOUT replacement	 
	while (exists($seentx{$randnum})) {$randnum = int(rand($txid))};
	 	 
	$seentx{$randnum} = 1;
	 	 
	# Parsing the header line
	$header = $seqhash{$randnum}[2];
	$reject = 0;
		 
	if ($header =~ /(FBtr[0-9]+).+loc\=([a-zA-Z0-9\_]+)\:([a-zA-Z]*)\(*([0-9\.\,]+).+(FBgn[0-9]{7})/) {
	 
		 $fbtr = $1;			
		 $chr = $2;
		 
		 if ($chr =~ /ito/) {
			 $mitocount += 1;
			if ($option > 0) {$reject += 1};
		 }
		 
		 if ($chr =~ /U/) {
			$Ucount += 1;
			if ($option > 1) {$reject += 1};
		 }
			 
		
		$strand = $3;
		$coords = $4;
		$FB = $5;  # This is the parent gene
		
		if ($option == 3) {
			if (exists ($FBseen{$FB})) {
				$reject += 1;
				$isoformrej += 1;			
			}
		}
		
		$coords =~ s/\(//;
		$coords =~ s/\)//;
		
		if ($strand eq "complement") {$strand = "-"};
		if ($strand eq "join") {$strand = "+"};
		if ($strand eq "") {$strand = "+"};	

		@coordstring = ();
		if ($coords =~ /\,{1,}?/) {
			@coordstring = split(",",$coords);
		} else {
			@coordstring[0] = $coords;
			$singlecheck += 1;
		}

		$exoncount = 0;
			
		if ($reject < 1) {
			foreach my $coordstring (@coordstring) {  # Print junction coordinates
				$BED->printf("$chr\t");
				if ($exoncount > 0) {
					$junccount += 1;
					$juncid = "JUNC".sprintf("%05d",$junccount)."_".$fbtr."_";
					$JUNCBED->printf("$chr\t");
					$intronstart = $rightedge + 1;
					$JUNCBED->printf("$intronstart\t")
				}
				if ($coordstring =~ /([0-9]+)\.\.([0-9]+)/) {					
					$BED->printf($1 - 1);
					$BED->printf("\t$2\t");
					$leftedge = $1;
					$rightedge = $2 + 1;
					
					if ($exoncount > 0) {
						$intronend = $leftedge - 1;
						$JUNCBED->printf("$intronend\t");
						$JUNCBED->printf("$juncid$FB\t");
						$JUNCBED->printf("1\t");
						$JUNCBED->printf("$strand\n");							
					}					
				}
				$exoncount += 1;
				$BED->printf("$fbtr"."_"."$FB"."_"."$exoncount"."\t");
				$BED->printf("1\t");
				$BED->printf("$strand\n");
				$FBseen{$FB} = $FB;	
			}			
			$keptcount += 1;
			push(@exoncount, $exoncount);
			$RANDCDS->printf($seqhash{$randnum}[0]);
			print_sequence($seqhash{$randnum}[1], 50);
		
		} else {
			$rejectedcount += 1;
			$count = $count - 1;
		}
	
	 } else {	 
		# Here, the transcript is rejected if the header is malformed
		$rejectedcount += 1;
		$count = $count - 1; 
	 }
  
 
}

print "$keptcount random transcripts were selected\n";
print "While doing the random selection, $rejectedcount transcripts were rejected based on your options:\n";
print "\t$mitocount mapping to mitochondrial genome\n";
print "\t$Ucount mapping to 'chromosome u'\n";
print "\t$isoformrej isoforms from genes already selected\n\n";
print "--------------------------------------------------\n";
print "Characteristics of your random set of transcripts:\n";
print "--------------------------------------------------\n";

$median = median(\@exoncount);
print "Median exon count is $median\n";

my $singles;
my $min;
my $max;
my $i;
$max = $exoncount[0];
$min = $exoncount[0];

foreach $i (@exoncount[1..$#exoncount])
{
	if ($i > $max)
	{
		$max = $i;
	}
	elsif ($i < $min)
	{
		$min = $i;
	}
	
	if ($i == 1) {
	$singles += 1;
	}
	
}
print "The maximum # of exons is " . $max . "\n";
print "The minimum # of exons is " . $min . "\n";
print "There are $singles single exon transcripts\n";

