#!/usr/bin/perl
use warnings;
use strict;
use JSON;
use Getopt::Long;

my $inDir = '';

GetOptions(
   'dir=s' => \$inDir
);

my $count = 0;
select(STDOUT);
$|++;

opendir DATA, "$inDir" or die "Unable to open directory: $!";
my @data = readdir DATA;
chomp @data;
closedir DATA;

open my $file, ">vcfTracks.json" or die "Could not create vcfTracks.json file: $!";
print $file "[";

my %check;

foreach my $t (@data) {  
    if ($t =~ /\.vcf.gz$/) {
        $t =~ s/\.chr.+?\./\.{refSeq}\./;
        my $chrom = $1;
        my $jsonOut = {
            storeClass => "JBrowse/Store/SeqFeature/VCFTabix",
            urlTemplate => "http://54.235.177.158/1000genomes/phase1/analysis_results/integrated_call_sets/$t",
            tbiUrlTemplate => "http://54.235.177.158/1000genomes/phase1/analysis_results/integrated_call_sets/$t.tbi",
            label => $t,
            type => "JBrowse/View/Track/HTMLVariants",
            key => $t,
            metadata => { source => "1000 genomes variants" }
        };
        print $file encode_json($jsonOut),",\n" unless $check{$t};
        $check{$t} = 1;
    }    
    print "\r",$count++, " genomes processed";    
}

print $file "]";

close $file;
