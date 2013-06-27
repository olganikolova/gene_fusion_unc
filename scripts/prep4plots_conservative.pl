#!/usr/bin/perl

use strict;
#use warnings;
use Data::Dumper;
use POSIX;
use List::Util qw[min max];

my $fin = shift();
my $base = $fin;
$base =~ s/\.tmp//;	
my $scatGeneFout = $base.".scat.gene";		# for scatterplot per gene file out
my $scatFusFout = $base.".scat.fus";		# for scatterplot per fusion file out
my $circFout = $base.".circ";				# for circos map file out

open INF, "<$fin" or die "Cannot open input file $fin...\n";
open OUTSCATGENE, ">>$scatGeneFout" or die "Cannot open 
					output file $scatGeneFout...\n";
open OUTSCATFUS, ">>$scatFusFout" or die "Cannot open 
					output file $scatFusFout...\n";
open OUTCIRC, ">>$circFout" or die "Cannot open output file $circFout...\n";

# print headers
print OUTSCATGENE "GeneID\tFrequency\n";
print OUTSCATFUS "FusionID\tFrequency\n";
# no header for circos maps

my %patientDB = ();

my %Gene = ();
my %GenePatientIDs = ();
my %Fusion = ();
my %FusionPatientIDs = ();
my %Fusion2Circos = ();

# <SAGE_sampleID> <SAGE_donorGeneID> <SAGE_acceptorGeneID> <SAGE_donorGeneLength> <SAGE_acceptorGeneLength> <SAGE_donorChromosome> <SAGE_acceptorChromosome> <donorEnd> <acceptorEnd> <strand>
# Gene ID format: <E1,E2,..;H1,H2..>

while (my $line=<INF>) {
	chomp($line);
	
	# skip header
	if( $line =~ /^SAGE_sampleID/ && $line !~ /NA/ && $line !~ ""){
		next;
	}
	
	$line =~ s/"//g;
	my ($tcgaid, $idsD, $idsA, $lengthsD, $lengthsA, $chrD, $chrA, $posD, $posA, $strand) = split(/\t/, $line);
	
	my ($ensemblidsD, $hugoidsD) = split(/;/, $idsD);
	my @hugoidsDarr = split(/,/, $hugoidsD);
	my ($ensemblidsA, $hugoidsA) = split(/;/, $idsA);
	my @hugoidsAarr = split(/,/, $hugoidsA);
	
	for(my $i=0; $i < max(length(scalar(@hugoidsDarr)), length(scalar(@hugoidsAarr))); $i++) {

		# here we only analyze fusions where 
		# both D and A map to different genes
		my $dg = defined($hugoidsDarr[$i]) ? $hugoidsDarr[$i] : "";
		my $ag = defined($hugoidsAarr[$i]) ? $hugoidsAarr[$i] : "";
		#my $ag = $hugoidsAarr[$i];
		
		# skip anything that doesn't map to a 
		# tuple of mapped distinctive genes
		if( $dg !~ /[a-zA-Z0-9]+/ or 
			$ag !~ /[a-zA-Z0-9]+/ or
			$dg eq $ag ) {
			next;
		} 
		
		##########################################	
		# count patients
		##########################################	
		$patientDB{$tcgaid} = 1 unless(exists($patientDB{$tcgaid})); 

		##########################################	
		# count genes
		##########################################
		my @genes = ($dg, $ag);
		
		foreach my $g ( @genes ) {	
			if(!exists($Gene{$g})) {
				$GenePatientIDs{$g} = "$tcgaid;"; 
				$Gene{$g} = 1;
			}
			elsif (  $GenePatientIDs{$g} !~ /^$tcgaid;/ &&
					 $GenePatientIDs{$g} !~ /;$tcgaid;/ ){
				$GenePatientIDs{$g} .= "$tcgaid;";
				$Gene{$g}++;	
			}# end count genes
		}
		
		##########################################	
		# count fusions (direct output for circos maps)
		##########################################	
		my $fus = "$dg---$ag";
			
		if(!exists($Fusion{$fus})) {
			$FusionPatientIDs{$fus} = "$tcgaid;"; 
			$Fusion{$fus} = 1;
			
			# we only add this string once for each fusion
			# frequency parameter and z are added later
			$Fusion2Circos{$fus} = "$chrD $posD $posD;$chrA $posA $posA";
		} #end_if_fusion does not exist
		elsif (  $FusionPatientIDs{$fus} !~ /^$tcgaid;/ &&
				 $FusionPatientIDs{$fus} !~ /;$tcgaid;/ ){
			$FusionPatientIDs{$fus} .= "$tcgaid;";
			$Fusion{$fus}++;	

			#print "Fusion $fus exists for patient $tcgaid..\n";
		}# end_if_exists_fusion
	} #_for_i_hugo_ids		
} # while each line in file

# foreach my $gene (keys %Gene){
#  	print "$gene\t".$Gene{$gene}."\n";
# }

# write to files
my $tpat = scalar keys %patientDB;
print '"'.$tpat.'",';
print "\n";

# scatterplot per fusion
# circos map
my $cnt = 0;

# order hashes by value in desc order
my @sortedF = sort { $Fusion{$b} <=> $Fusion{$a} } keys %Fusion; 
my @sortedG = sort { $Gene{$b} <=> $Gene{$a} } keys %Gene; 


#foreach my $fus (keys %Fusion){
foreach my $fus (@sortedF) {
	my $freqF = $Fusion{$fus}/$tpat; # should be equal to the sum of the above
	# write scatter plot per fusion
	print OUTSCATFUS "$fus\t$freqF\n";	
	
	# write circos map input
	my $str = $Fusion2Circos{$fus}; #."stroke_thickness=$freqF,z=$freqF"; 
	my($line1, $line2) = split(/;/, $str);
	
	my $thick = 100*$freqF;
	
	print OUTCIRC "link$cnt $line1 stroke_thickness=$thick,z=$thick\n";
	print OUTCIRC "link$cnt $line2 stroke_thickness=$thick,z=$thick\n";
	$cnt++; 
}

#foreach my $gene (keys %Gene){
foreach my $gene (@sortedG){
	my $freq = $Gene{$gene}/$tpat;
	if($freq > 1){die"$gene\n";}
 	print OUTSCATGENE "$gene\t$freq\n";
}


close (INF); 
close (OUTSCATGENE);
close (OUTSCATFUS);
close (OUTCIRC);
