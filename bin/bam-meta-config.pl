#!/usr/bin/perl;
use warnings;
use strict;
use Getopt::Long;
use JSON;

my %population;
my %genome_centre;
my %individuals;
my %extras;
my $inDir = "";
my $metaDir = "1000genomes_metadata1.csv";

GetOptions(
    'dir=s' => \$inDir,
    'meta=s' => \$metaDir
);

$inDir .= "/" unless $inDir =~ /\/$/;

open META, "$metaDir" or die "Unable to open metadata file: $!";
my @meta = <META>;

#create a hash containing population acronyms and their names and another hash containing the genome centres and their acronyms
my @general = splice @meta, 0, 21;
foreach (@general) {
    $_ =~ s/\r//g;
    my @info = split /;/, $_;
    my $pop_acr = $info[1];
    my $pop = $info[2];
    $population{$pop_acr} = $pop;
    my $centre_acr = $info[3];
    my $centre = $info[4];
    $genome_centre{$centre_acr} = $centre;
    my $web_label = $info[6];
    my $link = $info[8];
    $extras{$web_label} = $link;
}

my @samples = splice @meta, 7;
foreach (@samples) {
    my @indiv = split /;/, $_;
    my $popu = $indiv[0];
    my $popul = $popu . "-" . $population{$popu};
    my $acnum = $indiv[1];
    my $coriell_ID = $indiv[2];
    my $family = $indiv[3];
    my $gender = $indiv[4];
    my $relation = $indiv[5];
    my $type = $indiv[6];
    my $genome_cen = $indiv[7];
    #my $genome_centr = $genome_cen . "-" . $genome_centre{$genome_cen};
    my $exome_cen = $indiv[8];
    #my $exome_centr = $exome_cen . "-" . $genome_centre{$exome_cen};
    my $pilot1 = $indiv[9];
    my $pilot2 = $indiv[10];
    my $pilot3 = $indiv[11];
    my $omni_gen = $indiv[12];
    my $axiom_gen = $indiv[13];
    my $aligned = $indiv[14];
    chomp(my $exome_targets_covered = $indiv[15]);
    #create a hash of arrays... key = coriell id and value contains an array with all the associated data
    $individuals{$coriell_ID} = [$popul, $acnum, $coriell_ID, $family, $gender, $relation, $type, $genome_cen, $exome_cen, $pilot1, $pilot2, $pilot3, $omni_gen, $axiom_gen, $aligned, $exome_targets_covered];
}

my $count = 0;
select(STDOUT);
$|++;

opendir DATA, "$inDir" or die "Unable to open directory: $!";
my @data = readdir DATA;
chomp @data;
closedir DATA;

open my $file, ">bamTracks.json" or die "Could not create bamTracks.json file: $!";
print $file "[";

foreach my $t (@data) {
    if (-d "$inDir$t/alignment") {
        opendir HG, "$inDir$t/alignment";
        my @chrom = readdir HG;
        closedir HG;
        foreach (@chrom) {
            if ($_ =~ /\.bam$/) {
                my %style = ( className => "alignment",
                             arrowheadClass => "arrowhead",
                             labelScale => 100);
                my %meta = ( population => $individuals{$t}->[0],
                            sampleAcNum => $individuals{$t}->[1],
                            coriellID => $individuals{$t}->[2],
                            family => $individuals{$t}->[3],
                            gender => $individuals{$t}->[4],
                            relationship => $individuals{$t}->[5],
                            type => $individuals{$t}->[6],
                            genomeCentre => $individuals{$t}->[7],
                            exomeCentre => $individuals{$t}->[8],
                            pilot1 => $individuals{$t}->[9],
                            pilot2 => $individuals{$t}->[10],
                            pilot3 => $individuals{$t}->[11],
                            omniGenome => $individuals{$t}->[12],
                            axiomGenome => $individuals{$t}->[13],
                            aligned => $individuals{$t}->[14],
                            exomeTargetsCovered => $individuals{$t}->[15]
                            );
                my $json = { urlTemplate => "http://54.235.177.158//1000genomes/data/$t/alignment/$_",
                            label => $_,
                            type => "Alignments2",
                            key => $_,
                            style => \%style,
                            metadata => \%meta
                            };
                my $jsonText = to_json($json);
                print $file "$jsonText,\n";
            }
        }
        print "\r",$count++, " genomes processed";
    }
}

print $file "]";

close $file;