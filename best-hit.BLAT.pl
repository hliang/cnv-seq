#!/usr/bin/perl
use strict;

die "usage: $0 blat.psl [blat.2.psl, ...]\n" unless @ARGV;

my %blat;
while(<>)
{
	if(m/^(\d+)\t(\d+\t){7}[+-]\t(\S+)\t(\d+\t){3}(\S+)\t\d+\t(\d+)\t(\d+)/)
	{
		my ($match, $other, $query, $other2, $chrom, $chromstart, $chromend) = ($1,$2,$3,$4,$5,$6,$7);
		if($match > ${$blat{$query}}[0])
		{
			$blat{$query} = [$match, $chrom, $chromstart, $chromend];
		}
	}
}

foreach(keys %blat)
{
	my @temp = @{$blat{$_}};
	my $mid = int(($temp[2]+$temp[3])/2);
	print "$temp[1]\t$mid\n";
}
