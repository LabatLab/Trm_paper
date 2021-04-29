#!usr/bin/perl

use strict;
use warnings;

opendir(DIR, ".") or die "Cannot open directory.";

my $output_file = './results_concat.tsv';
open(OUTPUT, '>', $output_file) or die $!;

my @docs = grep(/\.filtered\.condensed\.ranked\.tsv$/, readdir(DIR));

if (scalar(@docs) == 0) {
	die "No expected result files found.";
}

# Define variables storing intermediate results
my %HGVSc_summary;
my %gene_summary;
my %neoantigen_sample_count_summary;
my %shared_private_summary;
my %samples_per_patient;

##########################################

# Go through each sample result file
foreach my $file (@docs) {
    open (RES, $file) or die "Could not open result file: $file\n";
    
    my $sample = $file =~ s/\.filtered\.condensed\.ranked\.tsv//r;
    my $patient = $sample =~ s/(.*)_.*/$1/r;
    $samples_per_patient{$patient}++;
    
    # Make sure a zero found of neoantigens is initialized for each sample
	if (!defined($neoantigen_sample_count_summary{$sample})) {
		$neoantigen_sample_count_summary{$sample} = 0;
	}
   
   # Go through each line in the result file
    while (<RES>){
    	# Ignore the header
    	if ($_ =~ /^Gene/) {
    		next;
    	}
    	
    	my @columns = split(/\t/, $_);
    	
    	# Output to the concatenated results file with the sample name in the first column
        print OUTPUT $sample."\t$_";
        
        # Save the HGVSc summary results (samples for each unique neoantigen)
        $HGVSc_summary{$columns[3]}{"counts"}++;
        $HGVSc_summary{$columns[3]}{"samples"} .= $sample."\t";
        
        # Save the gene summary results (number of times a gene contains neoantigens)
        $gene_summary{$columns[0]}{"counts"}++;
        $gene_summary{$columns[0]}{"samples"} .= $sample."\t";
        
        # Save a per-patient summary of how many times each unique neoantigen was seen
        $shared_private_summary{$patient}{$columns[3]}++;
        
        # Save the counts of neoantigens per sample
       	$neoantigen_sample_count_summary{$sample}++;
    }
}

close(OUTPUT);

##########################################

$output_file = './results_samples_per_unique_neoantigen.tsv';
open(OUTPUT, '>', $output_file) or die $!;

print OUTPUT "HGVSc\tCount\tSamples\n";

# Output the HGVSc summary results
foreach my $key (keys %HGVSc_summary) {
	print OUTPUT $key."\t".$HGVSc_summary{$key}{"counts"}."\t".$HGVSc_summary{$key}{"samples"}."\n";
}

close OUTPUT;

##########################################

$output_file = './results_samples_per_unique_neoantigen.tsv';
open(OUTPUT, '>', $output_file) or die $!;

print OUTPUT "HGVSc\tCount\tSamples\n";

# Output the HGVSc summary results
foreach my $key (keys %HGVSc_summary) {
	print OUTPUT $key."\t".$HGVSc_summary{$key}{"counts"}."\t".$HGVSc_summary{$key}{"samples"}."\n";
}

close OUTPUT;

##########################################

$output_file = './results_shared_private_neoantigen_summary.tsv';
open(OUTPUT, '>', $output_file) or die $!;

print OUTPUT "Patient\tPrivate\tShared\tTrunk\n";

foreach my $patient (sort keys %samples_per_patient) {
	my $neoantigens_private_count = 0;
	my $neoantigens_shared_count = 0;
	my $neoantigens_trunk_count = 0;
	
	if (defined($shared_private_summary{$patient})) {
		foreach my $HGVSc (keys %{$shared_private_summary{$patient}}) {
			if ($shared_private_summary{$patient}{$HGVSc} == 1) {
				$neoantigens_private_count++;
			} elsif ($shared_private_summary{$patient}{$HGVSc} == $samples_per_patient{$patient}) {
				$neoantigens_trunk_count++;
			} else {
				$neoantigens_shared_count++;
			}
		}
	}
	
	print OUTPUT $patient."\t".$neoantigens_private_count."\t".$neoantigens_shared_count."\t".$neoantigens_trunk_count."\n";
}

close OUTPUT;

##########################################

$output_file = './results_unique_neoantigen_count_per_sample_combination.tsv';
open(OUTPUT, '>', $output_file) or die $!;

print OUTPUT "Samples\tCount\n";

my %unique_neoantigens_per_sample_combination;

# Pull out unique sample combinations for each unique neoantigen
foreach my $key (keys %HGVSc_summary) {
	$unique_neoantigens_per_sample_combination{$HGVSc_summary{$key}{"samples"}}++;
}

# Output the HGVSc summary results
foreach my $key (keys %unique_neoantigens_per_sample_combination) {
	my $samples = $key =~ s/\t/;/rg;
	chop($samples);
	print OUTPUT $samples."\t".$unique_neoantigens_per_sample_combination{$key}."\n";
}

close OUTPUT;

##########################################

$output_file = './results_neoantigens_per_gene.tsv';
open(OUTPUT, '>', $output_file) or die $!;

print OUTPUT "Gene\tCount\tSamples\n";

# Output the gene summary results
foreach my $key (sort keys %gene_summary) {
	print OUTPUT $key."\t".$gene_summary{$key}{"counts"}."\t".$gene_summary{$key}{"samples"}."\n";
}

close OUTPUT;

##########################################

$output_file = './results_neoantigens_per_sample.tsv';
open(OUTPUT, '>', $output_file) or die $!;

print OUTPUT "Sample\tCount\n";

# Output the HGVSc summary results
foreach my $key (sort keys %neoantigen_sample_count_summary) {
	print OUTPUT $key."\t".$neoantigen_sample_count_summary{$key}."\n";
}

close OUTPUT;

##########################################

exit;