#!/usr/local/bin/perl

use strict;
use Data::Dumper;
use POSIX;

my $fin = shift();
my $fout = $fin;
$fout =~ s/\.output/\.hist/;

open INF, "<$fin" or die "Cannot open input file $fin...\n";
open OUTF, ">>$fout" or die "Cannot open output file $fout...\n";

my $count = 0;
my %db = ();

while (my $line=<INF>) {
  if( $line !~ /"donorGeneID/ && $line !~ /NA/ && $line !~ ""){
	chomp($line);
	$line =~ s/"//g;
	#print Dumper(\$line);
	my @parts = split(/\t/, $line);
	my @ids1 = split(/;/, $parts[0]);
	my @genes1 = split(/,/, $ids1[1]);
	
	my @ids2 = split(/;/, $parts[1]);
	my @genes2 = split(/,/, $ids2[1]);
	
	my $fus = $genes1[0]."--".$genes2[0];
	chomp($fus);
	
	if(	exists($db{$fus}) ) {
		$db{$fus}++;
	}
	else {
		$db{$fus} = 1;
	}
	}
	$count++;
}

#print Dumper(\%db);

foreach( keys (%db) ) {
	my $freq = $db{$_}/$count;
		print OUTF join("\t", "$_", $freq)."\n";
}
close (INF); 
close (OUTF); 
