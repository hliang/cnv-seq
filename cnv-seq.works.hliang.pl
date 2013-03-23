#!/usr/bin/perl
use strict;
use Getopt::Long;
use List::Util qw(min max sum);

my $usage = <<'USAGE';

################ CNV-seq ################

	usage: cnv.seq.pl [options]

		--test = test.hits.file
		--ref = ref.hits.file

		--log2-threshold = number
		  (default=0.6)
		--p-value = number
		  (default=0.001)

		--bigger-window = number
		  (default=1.5)
		--window-size
		  (default: determined by log2 and p-value;
		   if set, then log2 and p-value are not used)

		--genome = (human, chicken, chrom1, autosome, sex, chromX, chromY)
		  (higher priority than --genome-size)
		--genome-size = number
		  (no use when --genome is avaliable)

		--annotate
		--no-annotate
		  (default: do annotation, ie --annotate)
		--minimum-windows-required = number 
		  (default=4; only for annotation of CNV)
		--global-normalization
		  (if used, normalization on whole genome,
		  instead of oneach chromosome)

		--Rexe = path to your R program
		  (default: R)

		--help
		--debug


#########################################

USAGE

my ($test_hit, $ref_hit, $genome_size);
my $log2 = 0.6;
my $pvalue = 0.001;
my $genome;
my $min_window=4;
my $chromosomal_normalization = 'TRUE';
my $annotation='TRUE';
my $window_size;
my $bigger = 1.5;
my $debug = 0;
my $Rexe = 'R';

my $result = GetOptions(
	"test=s" => \$test_hit,
	"ref=s" => \$ref_hit,
	"log2-threshold=f" => \$log2,
	"p-value=f" => \$pvalue,
	"bigger-window=f" => \$bigger,
	'window-size=i' => \$window_size,
	"genome=s" => \$genome,
	"genome-size=o" => \$genome_size,
	"global-normalization" => sub{$chromosomal_normalization='FALSE'},
	"no-annotate" => sub{$annotation='FALSE'},
	"annotate" => sub{$annotation='TRUE'},
	"Rexe=s" => \$Rexe,
	"debug" => sub{$debug='TRUE'},
	"minimum-windows-required=i" => \$min_window,
	"help|?" => sub{print $usage; exit}
);

die $usage unless($test_hit && $ref_hit);
open(TEST,$test_hit) or die "can not find $test_hit file\n";
open(REF, $ref_hit) or die "can not find $ref_hit file\n";
die "must set --genome or --genome-size or --window-size\n" unless($genome_size||$genome||$window_size);

my $out;
my $temp = $test_hit;
$temp =~ s/.+\///;
$out = $temp;
$temp = $ref_hit;
$temp =~ s/.+\///;
$out .= "-vs-$temp";
my ($total_test, $total_ref);
unless($window_size)
{
	$genome_size = 3253037807 if $genome=~m/human/i;
	$genome_size = 1050947331 if $genome=~m/chicken/i;
	$genome_size = 247249719  if $genome=~m/chrom1|chr1|chromosome1/i;
	$genome_size = 2867732772 if $genome=~m/auto/i;
	$genome_size = 212686708  if $genome=~m/sex|xy/i;
	$genome_size = 154913754  if $genome=~m/chrX|chromX|chromosomeX/i;
	$genome_size = 57772954   if $genome=~m/chrY|chromY|chromosomeY/i;
	print "genome size used for calculation is $genome_size\n";
	die "$usage\n ERROR: check genome size pls\n" unless $genome_size; 

	$total_test++ while(<TEST>);
	close TEST;
	print "$test_hit: $total_test reads\n";

	$total_ref++ while(<REF>);
	close REF;
	print "$ref_hit: $total_ref reads\n";

	# there was an typo in the CNV-seq paper:
	# Equation 6 should be qnorm(1 - 0.5*p), instead of qnorm(0.5*(1-p))
	my $bt = `echo 'options(digits=16); qnorm(1-0.5*$pvalue)' | $Rexe --vanilla --slave`;
	die "\n Error:\ncan not find program $Rexe" unless $bt;
	my $st = `echo 'options(digits=16); qnorm(0.5*$pvalue)' | $Rexe --vanilla --slave`;
	die "\n Error:\ncan not find program $Rexe" unless $st;
	$bt = $1 if $bt =~ m/^\[1\]\s+([\d.e\+\-]+)/m;
	$st = $1 if $st =~ m/^\[1\]\s+([\d.e\+\-]+)/m;
	$log2 = abs($log2);
	my $brp = 2**$log2;
	my $srp = 1/(2**$log2);
	my $bw = ($total_test * $brp**2 + $total_ref) * $genome_size * $bt**2 / ((1-$brp)**2 * $total_test * $total_ref);
	my $sw = ($total_test * $srp**2 + $total_ref) * $genome_size * $st**2 / ((1-$srp)**2 * $total_test * $total_ref);
	$window_size = max($bw, $sw);
	print "The minimum window size for detecting log2>= $log2 should be $bw\n";
	print "The minimum window size for detecting log2<=-$log2 should be $sw\n";
	#print "window size to use is $window_size x $bigger = ";
	print "window size (even number) to use is $window_size x $bigger ~= ";
	$window_size *= $bigger;
	$window_size = sprintf("%.0f", $window_size);
	#change window size to an even number. will help bypass a bug in this script.
	$window_size++ if $window_size%2;
	print "$window_size\n";
	$out .= ".log2-$log2.pvalue-$pvalue";
}else
{
	#change window size to an even number. will help bypass a bug in this script.
	$window_size++ if $window_size%2;
	$out .= ".window-$window_size";
}
my $step = sprintf("%.0f",$window_size/2);
#print "window size to be used: $window_size\n";
print "window size (even number) to be used: $window_size\n";

my $cnvout = $out;
$cnvout .= ".minw-$min_window" if $annotation eq 'TRUE';

print "start counting test hits ... \n" if $debug;
my (%count, %test_total, %ref_total);
my $pat = qr/^\s*(\S+)\s+([\d.eE\+-]+)\s*$/;
my %chrom;
my $n=0;
open(TEST,$test_hit) or die;
open(REF, $ref_hit) or die;
while(<TEST>)
{
	if(m/$pat/)
	{
		$n++;
		if($debug)
		{
			print "$n of $total_test ...\n" if $n%100000==0;
		}
		$chrom{$1}=1;
		my $window_pos = int($2/$step);
		$count{$1}->[0][$window_pos]++;
		$count{$1}->[0][$window_pos-1]++ if $window_pos>0;
		print "\thit at window $window_pos ... \n" if $debug && $n<10;
		$count{$1}->[1][$window_pos]+=0;
		$count{$1}->[1][$window_pos-1]+=0 if $window_pos>0;
	}
}
print "read $n test reads, out of $total_test lines\n";
print "start counting ref hits ... \n" if $debug;
$n=0;
while(<REF>)
{
	if(m/$pat/)
	{
		$n++;
		if($debug)
		{
			print "$n of $total_ref ...\n" if $n%100000==0;
		}
		$chrom{$1}=1;
		my $window_pos = int($2/$step);
		$count{$1}->[1][$window_pos]++;
		$count{$1}->[1][$window_pos-1]++ if $window_pos>0;
		print "\thit at window $window_pos ... \n" if $debug && $n<10;
		$count{$1}->[0][$window_pos]+=0;
		$count{$1}->[0][$window_pos-1]+=0 if $window_pos>0;
	}
}
print "read $n ref reads, out of $total_ref lines\n";

if($debug)
{
	print "done counting:\n";
	use Data::Dumper;
	my $chrom = Dumper(%chrom);
	my $count = Dumper(%count);
	print '%chrom: ', substr($chrom, 0, 100), " ...\n";
	print '%count: ', substr($count, 0, 100), " ...\n";
}

open(OUT, ">$cnvout.count") or die;
print "write read-counts into file: $out.count\n";
print OUT "chromosome\tstart\tend\ttest\tref\n";
#foreach my $chr(sort{$a<=>$b} keys %chrom)
foreach my $chr(sort{$a cmp $b} keys %chrom)
{
	my @test = @{$count{$chr}->[0]};
	my @ref = @{$count{$chr}->[1]};
	for my $id (0..min($#test,$#ref))
	{
		my $start = $step*$id+1;
		my $end = $start+$window_size-1;
		my $test_raw = $test[$id]+0;
		my $ref_raw = $ref[$id]+0;
		my ($test, $ref)=(0,0);
		print OUT "$chr\t$start\t$end\t$test_raw\t$ref_raw\n";
	}
}
close OUT;

print "$Rexe package cnv output: $cnvout.cnv\n";
system(qq`echo 'library(cnv);\ncnv<-cnv.cal("$cnvout.count", log2=$log2, min=$min_window, chromosomal=$chromosomal_normalization, annotate=$annotation);\nwrite.table(cnv,"$cnvout.cnv",row.names=FALSE,sep="\\t");\n' | $Rexe --vanilla --slave`);
