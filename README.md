NOTE: There are now better tools available to do this conversion process, such as `vcfutils ucscsnp2vcf` in the SAMtools 
package and `VariantsToVCF` function in the Genome Analysis Toolkit (GATK), and I recommend these tools over dbSNP2vcf.pl.

dbSNP2vcf.pl is a (simplistic) Perl script to convert a dbSNP flat database into a tabix-ed bgzip-ed vcf file.
This is useful for generating VCF files for non-human organisms to use in conjunction with, for example, GATK (http://www.broadinstitute.org/gatk/).
