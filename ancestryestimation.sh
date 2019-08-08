
#Prune study data
plink --vcf ../data/HumanCoreExome_PhLiPS_ipsc_impute2_reheader.vcf.gz --recode --out ipsc
plink --bfile ipsc --exclude range high-LD-regions.txt --indep-pairwise 50 5 0.2 --out ipsc 
plink --bfile ipsc --extract ipsc.prune.in --make-bed --out ipsc.pruned


#Get reference data from 1000 genomes
pgen=https://www.dropbox.com/s/afvvf1e15gqzsqo/all_phase3.pgen.zst?dl=1
pvar=https://www.dropbox.com/s/0nz9ey756xfocjm/all_phase3.pvar.zst?dl=1
sample=https://www.dropbox.com/s/yozrzsdrwqej63q/phase3_corrected.psam?dl=1

wget $pgen
mv 'all_phase3.pgen.zst?dl=1' all_phase3.pgen.zst
plink2 --zst-decompress all_phase3.pgen.zst > all_phase3.pgen

wget $pvar
mv 'all_phase3.pvar.zst?dl=1' all_phase3.pvar.zst
plink2 --zst-decompress all_phase3.pvar.zst > all_phase3.pvar

wget $sample
mv 'phase3_corrected.psam\?dl\=1' all_phase3.psam

#Convert pfiles to vcf to filter out multiallelic variants
plink2 --pfile all_phase3 --recode vcf
bgzip -c plink2.vcf > all_phase3.vcf.gz
tabix -p vcf all_phase3.vcf.gz

#Filter out multiallelic variants
bcftools view --max-alleles 2 --exclude-types indels all_phase3.vcf.gz

#Convert vcf back to plink files
plink --vcf all_phase3.vcf.gz --recode --out 1000genomes

#Reannotate  plink files
plink --bfile ipsc --update-map /project/labprojects/projects/Morrisey/dave/Teachey.RNA-seq/annovar/humandb/snplist4.txt --update-name --make-bed --out ipsc2
plink --bfile ipsc2 --exclude range high-LD-regions.txt --indep-pairwise 50 5 0.2 --out ipsc_new
plink --bfile ipsc2 --extract ipsc_new.prune.in --make-bed --out ipsc_new.pruned

#Filter reference data for the same SNP set as in study
#plink --bfile /project/labprojects/Genotype/1kg_plink/1kg_phase1_all --extract ipsc.prune.in --make-bed --out 1000genomes_phase1.pruned
plink --bfile reference/1000genomes_p3 --extract ipsc_new.prune.in --make-bed --out 1000genomes_new.pruned

#Check and correct mismatch between reference and pruned dataset
awk 'BEGIN {OFS="\t"} FNR==NR {a[$2]=$1; next} \
($2 in a && a[$2] != $1) {print a[$2],$2}' \
ipsc_new.pruned.bim 1000genomes_new.pruned.bim | \
sed -n '/^[XY]/!p' > 1000genomes.toUpdateChr

plink --bfile 1000genomes_new.pruned \
--update-chr 1000genomes.toUpdateChr 1 2 \
--make-bed \
--out 1000genomes_new.updateChr

name = 'ipsc_new'
refname = '1000genomes_new'
#check for position mismatch
awk 'BEGIN {OFS="\t"} FNR==NR {a[$2]=$4; next} \
($2 in a && a[$2] != $4) {print a[$2],$2}' \
ipsc_new.pruned.bim 1000genomes_new.pruned.bim > \
1000genomes_new.toUpdatePos

#check for allele flips
awk 'BEGIN {OFS="\t"} FNR==NR {a[$1$2$4]=$5$6; next} \
($1$2$4 in a && a[$1$2$4] != $5$6 && a[$1$2$4] != $6$5) {print $2}' \
ipsc_new.pruned.bim 1000genomes_new.pruned.bim > \
1000genomes_new.toFlip

#correct for position mismatch and allele flip
plink --bfile 1000genomes_new.updateChr \
--update-map 1000genomes_new.toUpdatePos 1 2 \
--flip 1000genomes_new.toFlip \
--make-bed \
--out 1000genomes_new.flipped

#Remove any mismatches
awk 'BEGIN {OFS="\t"} FNR==NR {a[$1$2$4]=$5$6; next} \
($1$2$4 in a && a[$1$2$4] != $5$6 && a[$1$2$4] != $6$5) {print $2}' \
ipsc_new.pruned.bim 1000genomes_new.flipped.bim > \
1000genomes_new.mismatch

plink --bfile 1000genomes_new.flipped \
--exclude 1000genomes_new.mismatch \
--make-bed \
--out 1000genomes_new.clean


#MErge
plink --bfile ipsc_new.pruned \
--bmerge 1000genomes_new.clean.bed 1000genomes_new.clean.bim \
1000genomes_new.clean.fam \
--make-bed \
--out ipsc_new.merge.1000genomes_new

#Run PCA
plink --bfile ipsc_new.merge.1000genomes_new \
--pca \
--out ipsc_new.1000genomes_new
