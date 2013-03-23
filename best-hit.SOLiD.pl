#!/usr/bin/perl
use strict;
use List::Util qw(min shuffle);

die "usage: $0 SOLiD.csfasta [more SOLiD match files]\n" unless @ARGV;

while(<>)
{
	next unless m/^>/;
	my %map;
	push @{$map{$3}}, "$1\t$2\n" while(m/\b(\d+)_-?(\d+)\.(\d+)/g);
	my @map = shuffle(@{$map{min(keys(%map))}});
	print shift @map;
}
