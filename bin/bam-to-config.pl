#!/usr/bin/perl;
use warnings;
use strict;

my $count = 0;
select(STDOUT);
$|++;

opendir DATA, "/jbrowse/www/jbrowse/1000genomes/data" or die "Unable to open directory: $!";
my @data = readdir DATA;
chomp @data;
closedir DATA;

open my $file, ">bamTracks.json" or die "Could not create bamTracks.json file: $!";
print $file "[";

foreach my $t (@data) {
    if (-d "/jbrowse/www/jbrowse/1000genomes/data/$t/alignment") {
        opendir HG, "/jbrowse/www/jbrowse/1000genomes/data/$t/alignment";
        my @chrom = readdir HG;
        closedir HG;
        foreach (@chrom) {
            if ($_ =~ /\.bam$/) {
                print $file "{\n\t\"storeClass\" : \"/jbrowse/www/jbrowse/1000genomes/data/$t/$_\",";
                print $file "\n\t\"urlTemplate\" : \"http://s3.amazonaws.com/1000genomes/data/$t/alignment/$_\",";
                print $file "\n\t\"label\" : \"$_\",";
                print $file "\n\t\"type\" : \"JBrowse/View/Track/Alignments2\",";
                print $file "\n\t\"key\" : \"$_\",";
                print $file "\n\t\"style\" : {\n\t\t\"className\" : \"alignment\",\n\t\t\"arrowheadClass\" \"arrowhead\",\n\t\t\"labelScale\" : 100\n\t}\n},";
            }
        }
        print "\r",$count++, " genomes processed";
    }
}

print $file "]";

close $file;
