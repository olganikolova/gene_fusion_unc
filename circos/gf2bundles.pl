#!/usr/local/bin/perl

use strict;
use Data::Dumper;
use POSIX;

my $fin = shift();
my $fout = $fin;
$fout =~ s/\.txt/\.circos\.txt/;

open INF, "<$fin" or die "Cannot open input file $fin...\n";
open OUTF, ">>$fout" or die "Cannot open output file $fout...\n";

my @coord = ();
my $count = 0;
my %db = ();
# threshold for frequency across samples
my $thresh = 0.05;
my $nsamples = <INF>;
chomp($nsamples);
$nsamples =~ s/\s+//g;

while (my $line=<INF>) {
  unless( $line =~ /donorChromosome/ ){ # skip the header
	chomp;
	$line =~ s/chr/hs/g;
	
	if(	exists($db{$line}) ) {
		$db{$line}++;
	}
	else {
		$db{$line} = 1;
	}
	}
}

foreach( keys (%db) ) {
	my $line = $_;
	my $freq = $db{$line};
	my $thick = floor(100*($freq/$nsamples));
	my $z = 0;
	
# 	if($db{$line} > 100 ){
# 		$th = 100;
# 		$z = 12
# 	}
# 	elsif($db{$line} > 50 ){
# 		$th = 50;
# 		$z = 6;
# 	}
# 	elsif($db{$line} > 25 ) {
# 		$th = 25;
# 		$z = 3;
# 	}
# 	else {
# 		$th = 2;
# 	}
	
	@coord = split(/\s/, $line);
	if( $coord[0] =~ /hs[\d+|X|Y]/ && 
		$coord[3] =~ /hs[\d+|X|Y]/ &&
		$db{$line}/$nsamples >= $thresh) {
        #if( $coord[0] =~ /hs\d+/ && $coord[3] =~ /hs\d+/) {
		print OUTF join(" ", "link$count", @coord[0..2], "stroke_thickness=$thick,z=$z")."\n";
		print OUTF join(" ", "link$count", @coord[3..5], "stroke_thickness=$thick,z=$z")."\n";
		$count++;
	}
}
close (INF); 
close (OUTF); 
