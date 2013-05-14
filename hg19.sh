#!/bin/bash

echo "This script will create the following directory structure for hg19 data\n"
echo "goldenPath/hg19/chromosomes and goldenPath/hg19/database\n"
echo "hg19 database files will now be downloaded, this may take a while\n\n"

mkdir -p goldenPath/hg19/chromosomes
mkdir -p goldenPath/hg19/database

echo "Downloading reference sequences...\n"
rsync -ahvzP rsync://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/ goldenPath/hg19/chromosomes/
echo "Reference sequence download finished\n\n"

echo "Downloading database files\n"
rsync -ahvzP --exclude="chain*" --exclude="net*" rsync://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/ goldenPath/hg19/database/
echo "Database download finished\n\n"

echo "Adding reference sequence track\n\n"
for x in goldenPath/hg19/chromosomes/*.fa
do
    bin/prepare-refseqs.pl --fasta $x
done

tracks=( "wgRna:sno and miRNA"
        "gad:Genetic Association Studies of Complex Diseases and Disorders(GAD)"
        "cosmic:Catalogue of Somatic Mutations(COSMIC)"
        "gwasCatalog:NHGRI Genome-Wide Association Studies(GWAS)"
        "pseudoYale60:Yale Pseudogenes based on Ensembl 60"
        "wgEncodeGencodeCompV14:Gene models from GENCODE project"
        "lincRNAsTranscripts:lincRNA and TUCP transcripts"
        "lincRNAsCTAdipose:lincRNAs from adipose"
        "lincRNAsCTAdrenal:lincRNAs from adrenal"
        "lincRNAsCTBrain:lincRNAs from brain"
        "lincRNAsCTBreast:lincRNAs from breast"
        "lincRNAsCTColon:lincRNAs from colon"
        "lincRNAsCTHeart:lincRNAs from heart"
        "lincRNAsCTKidney:lincRNAs from kidney"
        "lincRNAsCTLiver:lincRNAs from liver"
        "lincRNAsCTLung:lincRNAs from lung"
        "lincRNAsCTLymphNode:lincRNAs from lymph node"
        "lincRNAsCTOvary:lincRNAs from ovary"
        "lincRNAsCTProstate:lincRNAs from prostate"
        "lincRNAsCTSkeletalMuscle:lincRNAs from skeletal muscle"
        "lincRNAsCTTestes:lincRNAs from testes"
        "lincRNAsCTThyroid:lincRNAs from thyroid"
        "lincRNAsCTWhiteBloodCell:lincRNAs from white blood cell" )

echo "Adding default tracks\n\n"        
for i in ${tracks[@]}
do
    in=(${i//:/ })
    bin/ucsc-to-json.pl --in goldenPath/hg19/database --track "${in[0]}" --key "${in[1]}" --clientConfig '{ "metadata" : { "source" : "UCSC genome browser tracks" }}'
done

echo "Adding UCSC Gene and RGD QTL tracks\n"
bin/ucscsql-to-gff3.pl --primaryTable knownGene --secondaryTable kgXref --link link kgID --primaryName geneSymbol --out .
bin/flatfile-to-json.pl --gff knownGene.gff3 --key 'UCSC Gene'

bin/ucscsql-to-gff3.pl --primaryTable rgdQtl --secondaryTable rgdQtlLink --link name name --out .
bin/flatfile-to-json.pl --gff rgdQtl.gff3 --key 'Human RGD QTL'

if [ -e /*bin/nginx || -e /usr/*bin/nginx ]
then
    echo "Please setup nginx as a reverse proxy for 1000genomes at http://s3.amazonaws.com/1000genomes\n"
    echo "Then run bam-meta-config.pl and vcf-config.pl to generate config files\n\n"
else
    echo "Please install nginx first then\n"
    echo "Setup nginx as a reverse proxy for 1000genomes at http://s3.amazonaws.com/1000genomes\n"
    echo "Then run bam-meta-config.pl and vcf-config.pl to generate config files\n\n"
fi

echo "Enabling faceted track selection\n"
bin/add-json '{ "trackSelector" : { "type" : "faceted" }, "dataset_id" : "hg19" }' data/trackList.json