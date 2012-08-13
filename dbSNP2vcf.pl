#!/usr/bin/perl -w
# A script to convert a dbSNP database into a tabix-ed bgzip-ed vcf file.
# Author: Peter Hickey (hickey@wehi.edu.au)
# Date: 02/11/2011

##########################################################################
# WARNING: This script has only been tested on dbSNP128 for the mm9 genome
##########################################################################

# VCF file format is #CHROM POS ID REF ALT QUAL FILTER INFO (however the ID, QUAL, and INFO fields can be ".").
# A header should be included in the vcf file, see example bleow (also see GATK human vcf for more examples of information that could be included)
##fileformat=VCFv4.0
##dbSNP_BUILD_ID=128 
##fileDate=20110721
##phasing=unknown
##reference=mm9
##source=dbSNP(http://hgdownload.cse.ucsc.edu/goldenPath/mm9/database/snp128.txt.gz)
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO

# Note: This script uses the refNCBI (dbSNP) reference genomic allele (see http://www.ncrna.org/glocal/cgi-bin/hgTables?db=mm9&hgta_group=varRep&hgta_track=snp128&hgta_table=snp128&hgta_doSchema=describe+table+schema)

# It is probably worthwhile running vcf-validator -d from the vcftools package to check for duplicate entries and other anomolies

use strict;
use warnings;
use Getopt::Long;

# Global variables
my $input;
my $output;
my $errout;
my $line_cnt;
# UPDATE $header as appropriate
my $header="##fileformat=VCFv4.0\n##dbSNP_BUILD_ID=128\n##fileDate=20110721\n##phasing=unknown\n##reference=mm9\n##source=dbSNP(http://hgdownload.cse.ucsc.edu/goldenPath/mm9/database/snp128.txt.gz)\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
my @values;
my $CHROM;
my $POS;
my $ID;
my $REF;
my $ALT;
my $QUAL;
my $FILTER;
my $INFO;
my $molType;
my @alleles;
my $isSNPval;

# Get input, output and error file names
GetOptions("i=s"=>\$input,"o=s"=>\$output,"e=s"=>\$errout);

if(!defined($input)||!defined($output)){
  print "Takes in a dbSNP file and converts it to a tabix-ed, bgzip-ed vcf file. Only SNPs can be converted, i.e. ignores indels and complex variants. \n\n";
  print "Usage:dbSNP2vcf.pl -i <dbSNP.txt> -o <output.vcf> -e <failed2covert.txt> \n\n";
  print "Peter Hickey (hickey\@wehi.edu.au), Bioinformatics Division, The Walter and Eliza Hall Institute of Medical Research. Last updated Nov 2011.\n\n";
  print "\t<dbSNP.txt> - unzipped dnSNP file e.g. http://hgdownload.cse.ucsc.edu/goldenPath/mm9/database/snp128.txt.gz \n\n";
  print "\t<output.vcf> - the output file of successfully converted variants\n\n";
  print "\t<failed2convert.txt> - a file with a list of variants that can't be converted to vcf format by this script\n\n";
  exit(1);
}

# Open file handles
open(IN1,"<$input") || die print "Can't open file $input\n";
print "Opening $input ...\n";
open(OUT1, ">$output") || die print "Can't open file $output\n";
print "Opening $output ...\n";
open(ERR1, ">$errout") || die print "Can't open file $errout\n";
print "Opening $errout ...\n";

print "Processing...\n";

$line_cnt=0;
print OUT1 "$header";

while (<IN1>) {
  chomp;
  $line_cnt++;
  if($line_cnt%1000000==0){
    print "Processed $line_cnt lines\n";
  }
  # Split each line by whitespace
  @values = split(/\s+/,$_);
  $CHROM=$values[1];
  $POS=$values[2]+1; # Genomic coordinates in dbSNP are 0-based, therefore we use start position + 1 to comply with vcf specification.
  $ID=$values[4];
  $REF=$values[7]; # refNCBI (Reference allele from dbSNP).
  @alleles = split(/\//,$values[9]); # Array of possible alleles
  @alleles = grep { $_ ne $REF } @alleles; # Remove the reference allele from the set of possible alleles
  $ALT=join(",", @alleles); # Join the alternate alleles by commas
  $QUAL="."; # Set to unknown
  $FILTER="PASS"; # As per Eric Banks' advice (http://getsatisfaction.com/gsa/topics/how_to_make_dnsnp_vcf_file_for_mouse_and_or_organisms)
  $INFO="."; # Set to unknown
  $molType=$values[10];
  $isSNPval=&isSNP(); # Check whether variant is a SNP or a more complex variant - this script only deals with SNPs
  if(($molType eq "genomic") & (length($REF)==1) & ($REF ne "-") & ($ALT !~ m/-/) & ($isSNPval)){ # Only using variants that are genomic, not "-" (deletions) and with 1bp reference allele (i.e. no long INDELS).
    print OUT1 "$CHROM\t$POS\t$ID\t$REF\t$ALT\t$QUAL\t$FILTER\t$INFO\n";
  } else{
    print ERR1 "@values\n";
  }
}

# bgzip and tabix the vcf file
# Comment out the following two lines to skip the bgzip and tabix steps
print "bgzip-ing and tabix-ing the vcf file\n";
system("bgzip $output; tabix -p vcf $output.gz");
  
# Close files
close(IN1);
close(OUT1);
close(ERR1);
print "File converted!\n";

# isSNV reports 1 if *all* alternate alleles are 1bp in length, 0 otherwise.
sub isSNP
  {
    my $alt_length;
    # print "@alleles\n";
    # foreach(@alleles){
    #   $alt_length=length($_);
    #   print "$alt_length\n";
    # }
      foreach(@alleles){
        $alt_length=length($_);
        if($alt_length>1){
        return(0);
        }
     }
    return(1);
  }
