#!/usr/bin/perl;
use warnings;
use strict;
use JSON;

my $count = 0;
select(STDOUT);
$|++;

opendir DATA, "/jbrowse/www/jbrowse/genomes/phase1/analysis_results/integrated_call_sets" or die "Unable to open directory: $!";
my @data = readdir DATA;
chomp @data;
closedir DATA;

open my $file, ">vcfTracks.json" or die "Could not create vcfTracks.json file: $!";
print $file "[";

foreach my $t (@data) {  
    if ($t =~ /\.vcf.gz$/) {
        my $jsonOut = {
            storeClass => "JBrowse/Store/SeqFeature/VCFTabix",
            urlTemplate => "http://54.235.177.158/1000genomes/phase1/analysis_results/integrated_call_sets/$t",
            label => $t,
            type => "JBrowse/View/Track/HTMLFeatures",
            key => $t,
            style => {
                className => "features",
                arrowheadClass => "arrowhead",
                labelScale => 100
            }
        };
        print $file encode_json($jsonOut),",\n";
    }    
    print "\r",$count++, " genomes processed";    
}

print $file "]";

close $file;