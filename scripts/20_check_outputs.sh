####Check Stats file from ipyrad output files####
##Remove bad samples if there are any, after going through the stats file##
##NB; I removed sample ***** from all outputs except .nex .seq.hdf5 .snps.hdf5 .str##

##Missingness (optional)
vcftools --vcf ipyrad_AMBYSE.vcf --missing-indv
##This generates an .imiss file for each individual sample

##Check Snp count (optional)
vcftools --vcf ipyrad_AMBYSE.vcf --freq --out vcf_qc
