# Code for The Maintenance of deleterious variation in wild Chinese rhesus macaques
### Camille Steux and Zachary A. Szpiech


### Filtering data
***
```bash
# Removing CR1 and CR2 individuals
bcftools view -s ^CR1,CR2 81wildChineseRhesus2.vcf.gz -Oz -o 79wildChineseRhesus.vcf.gz
# Filtering multi-allelic snps and indels, and keeping only passed sites
bcftools view -m2 -M2 -v snps -i 'FILTER="PASS"' 79wildChineseRhesus.vcf.gz -Oz -o 79wildChineseRhesus_passed.vcf.gz
# Setting GT to missing with GQ<20 or lowDP flag
bcftools filter 79wildChineseRhesus_passed.vcf.gz -i 'FMT/GQ>19 & FMT/FT!="lowDP"' -S . -Oz -o 79wildChineseRhesus_passed_missing.vcf.gz
# Filtering missing GT
bcftools view -i 'F_MISSING<0.2' 79wildChineseRhesus_passed_missing.vcf.gz -Oz -o 79wildChineseRhesus_passed_missing_filt.vcf.gz
bcftools stats 79wildChineseRhesus_passed_missing_filt.vcf.gz > 79wildChineseRhesus_passed_missing_filt.STAT.txt
```

### Lifting over
***
```bash
liftOver ~/Liftover/liftover-3/Variants_for_liftover_5pop.txt ~/Liftover/liftover-3/rheMac8ToRheMac10.over.chain.gz ~/Liftover/liftover-3/lifted_variants_5pop.txt ~/Liftover/liftover-3/unlifted_variants_5pop.txt
```

### Annotating variants with SIFT 4G
***
Building SIFT database (chromosome per chromosome to save time). The config file is available in the repository (config_SIFTdb_chr_i.txt).
```bash
perl make-SIFT-db-all.pl -config ~/SIFT_database/Mmul10_SIFT_chromo/chr${i}/config_SIFTdb_chr${i}.txt
```

Running SIFT4G annotator using the previously created database.
```bash
java -jar ~/Code/SIFT4G_Annotator.jar -c -t -i ~/Liftover/liftover-3/vcf/5pop_lifted_sorted.vcf -d ~/SIFT_database/Mmul10_SIFT_db/ -r ~/Variants_annotation/Mmul10_3/
```

Parsing the output.
```bash
# Keeping only the annotated variants in the vcf
bcftools view -i 'INFO/SIFTINFO ~ "|"' 5pop_lifted_sorted_SIFTpredictions.vcf -Ov -o 5pop_SIFTpredictions_only.vcf
# Removing the annotations with low confidence
bcftools view -e 'INFO/SIFTINFO ~ "WARNING"' 5pop_SIFTpredictions_only.vcf -Ov -o 5pop_SIFTpredictions_only_filt.vcf
# Extracting the annotations
bcftools query -f '%CHROM %POS %INFO/SIFTINFO \n' 5pop_SIFTpredictions_only_filt.vcf > 5pop_SIFTannotations_extract.txt
# Sorting the annotated variants in different files
awk -v OFS='\t' '$3 ~ "STOP" || $3 ~ "SYNONYMOUS" {print}' 5pop_SIFTannotations_extract.txt > 5pop_SIFTannotations_extract_cod.txt
awk -v OFS='\t' '{printf $1"\t"$2"\t"; split ($3,trans,","); for (t=1;t<=split($3,trans,",");t++) {split (trans[t],a,"|"); printf a[6]","a[13]"|"}; print ""}' 5pop_SIFTannotations_extract_cod.txt | sed "s/|$//g" > 5pop_SIFTannotations_extract_cod_red.txt
sed "s/STOP-GAIN,NA/STOP-GAIN,STOP-GAIN/g" 5pop_SIFTannotations_extract_cod_red.txt | sed "s/STOP-LOSS,NA/STOP-LOSS,STOP-LOSS/g" > 5pop_SIFTannotations_extract_cod_red_stop.txt
awk '$3 !~ "SYNONYMOUS,NA"' 5pop_SIFTannotations_extract_cod_red_stop.txt  > 5pop_SIFTannotations_extract_cod_red_stop_filt.txt
awk '{if ($3 ~ "STOP-GAIN") $4="STOP-GAIN,STOP-GAIN"; else if ($3 ~ "STOP-LOSS") $4="STOP-LOSS,STOP-LOSS"; else if ($3 ~ "NONSYNONYMOUS,DELETERIOUS") $4="NONSYNONYMOUS,DELETERIOUS"; else if ($3 ~ "SYNONYMOUS,DELETERIOUS") $4="SYNONYMOUS,DELETERIOUS"; else if ($3 ~ "NONSYNONYMOUS,TOLERATED") $4="NONSYNONYMOUS,TOLERATED"; else if ($3 ~ "SYNONYMOUS,TOLERATED") $4="SYNONYMOUS,TOLERATED"; print}' 5pop_SIFTannotations_extract_cod_red_stop_filt.txt > 5pop_SIFTannotations_extract_cod_red_stop_filt_uniq.txt
awk '{split ($4,a,","); print $1,$2,a[1],a[2]}' 5pop_SIFTannotations_extract_cod_red_stop_filt_uniq.txt > 5pop_SIFTannotations_uniq.txt
grep STOP 5pop_SIFTannotations_uniq.txt > 5pop_SIFTannotations_uniq_STOP.txt
grep SYNONYMOUS 5pop_SIFTannotations_uniq.txt > 5pop_SIFTannotations_uniq_SYN.txt
```


### Calling ROH with garlic
***
Creating the required input files to run garlic.
```bash
# Creating .tgls files containing genotype quality information for each site (with personal python script)
./create_tgls_from_vcf.py ~/Liftover/liftover-3/vcf/5pop_lifted_sorted.vcf.gz > ~/ROH_calling/roh-calling-3/tped_tfam_tgls/5pop_lifted_sorted.tgls
# Transforming vcf files into tped and tfam files to run garlic
plink --vcf ~/Liftover/liftover-3/vcf/5pop_lifted_sorted.vcf.gz --recode transpose --double-id --chr-set 20 --out ~/ROH_calling/roh-calling-3/tped_tfam_tgls/5pop_lifted_sorted
# Adding the population name as the first column of the tfam file
awk -v OFS="\ " '{print "pan",$2,$3,$4,$5,$6}' 5pop_lifted_sorted.tfam > 5pop_lifted_sorted_mod.tfam
```

Running garlic.
```bash
garlic --auto-winsize --auto-overlap-frac --winsize 100 --tped ~/ROH_calling/roh-calling-3/tped_tfam_tgls/5pop_lifted_sorted.tped --tfam ~/ROH_calling/roh-calling-3/tped_tfam_tgls/5pop_lifted_sorted_mod.tfam --tgls ~/ROH_calling/roh-calling-3/tped_tfam_tgls/5pop_lifted_sorted.tgls --gl-type GQ --resample 40 --centromere ~/ROH_calling/roh-calling-3/centromeres.txt --out ~/ROH_calling/roh-calling-3/roh-calls/5pop
```

## Computing the distribution of variants in ROH
***

Code for deleterious variants (identical for tolerated and LOF variants).

```bash
# Creating a file with the positions of the deleterious variants.
awk -v OFS='\t' '$4=="DELETERIOUS" {print $1, $2}' ../5pop_SIFTannotations_uniq_SYN.txt > deleterious_positions_5pop.txt
# Extracting the genotypes from the vcf
bcftools query -f '%CHROM %POS [ %GT]\n' -H -R deleterious_positions_5pop.txt ../5pop_lifted_sorted_SIFTpredictions.vcf.gz > deleterious_genotypes_5pop.txt
awk '{if ($1 ~ "#") {printf "CHROM POS"; for (t=4;t<=NF;t++) {split ($t,a,/[\] :]/); printf " "a[2]}; print ""} else {print $0}}' deleterious_genotypes_5pop.txt | sed 's/\ /\t/g' | sed 's/\t\t/\t/g' > deleterious_genotypes_5pop2.txt
# Extracting the genotypes for each individual seperately
for (( i=3; i<=81; i++ )); do ind=$(grep "CHROM" deleterious_genotypes_5pop2.txt | awk -v i=$i '{print $i}'); awk -v i=$i '{print $1,$2-1,$2,$i}' deleterious_genotypes_5pop2.txt > ../../../bedtools/bedtools-3/genotypes/deleterious_genotype_5pop_${ind}.bed; done
# Extracting only alternate homozygotes
for i in *; do arrname=(${i//./ }); name=${arrname[0]}; grep "1/1" $i > ${name}_homoalt.bed ; done

# Intersecting the deleterious homozygotes and the ROH with bedtools to get the deleterious sites falling within ROH
for i in *.allroh.bed ; do arrname=(${i//./ }); ind=${arrname[1]}; bedtools intersect -a 5pop.${ind}.allroh.bed -b ../genotypes/deleterious_genotype_5pop_${ind}_homoalt.bed > ../output_homozygotes/5pop_deleterious_intersect.${ind}.allroh.txt; done
# Counting the number of deleterious homozygotes falling inside and outside ROH
echo -e "ind\tnA\tnB\tnC\tnallroh\tntotal" > ../summary/deleterious_distribution_5pop.txt
for i in 5pop_deleterious_intersect*.txt; do arrname=(${i//./ }); ind=${arrname[1]}; echo -e "${ind}\t `grep 'A' $i | wc -l` \t `grep 'B' $i | wc -l` \t `grep 'C' $i | wc -l` \t `cat $i | wc -l` \t `cat ../genotypes/deleterious_genotype_5pop_${ind}_homoalt.bed | wc -l`"  >> ../summary/deleterious_distribution_5pop.txt; done
```

### Running shapeit
***
```bash
# For chromosome i:
shapeit2 --input-vcf ~/shapeit/vcf/5pop_stripped.chr$i.vcf --force --thread 8 --rho 0.00017043 -O ~/shapeit/output/phased.chr$i
# Converting the output into a vcf
shapeit2 -convert --input-haps ~/shapeit/output/phased.chr$i --output-vcf ~/shapeit/output/phased.converted.chr$i.vcf
```

### Running selscan
***
```bash
# For chromosome i:
for pop in brevicaudus lasiotis littoralis mulatta tcheliensis; do selscan --nsl --vcf ~/shapeit/output/phased.converted.${pop}.chr$i.vcf --threads 8 --out ~/selscan/output2/${pop}.chr$i; done
# Normalizing
for pop in brevicaudus lasiotis littoralis mulatta tcheliensis; do norm --nsl --files ${pop}.chr*.nsl.out --bp-win; done
```

### Constructing the SFS
***
```bash
## Intersecting deleterious variants with regions under selection with Bedtools
bedtools intersect -a  ~/bedtools/bedtools-selscan/input/regions_under_selection.bed -b ~/bedtools/bedtools-selscan/input/deleterious_positions_5pop.bed > ~/bedtools/bedtools-selscan/output/ss_delpos.bed
# Extracting deleterious variants from vcf
for pop in brevicaudus lasiotis littoralis mulatta tcheliensis; do bcftools view -R ~/Variants_annotation/Mmul10_3/variants_annot_summary/deleterious_positions_5pop.txt ~/shapeit/output/phased.converted.${pop}.vcf.gz -Oz -o ~/sfs/vcf/phased.converted.${pop}.del.vcf.gz ; done
for pop in brevicaudus lasiotis littoralis mulatta tcheliensis; do bcftools index ~/sfs/vcf/phased.converted.${pop}.del.vcf.gz; done
# Extracting deleterious variants in selected regions
for pop in brevicaudus lasiotis littoralis mulatta tcheliensis; do bcftools view -R ~/selscan/output2/${pop}.allchromo.extr.nsl.bed ~/sfs/vcf/phased.converted.${pop}.del.vcf.gz -Ov -o ~/sfs/vcf/phased.converted.${pop}.del.sel.vcf ; done
# Extracting deleterious variants not in selected regions
for pop in brevicaudus lasiotis littoralis mulatta tcheliensis; do bcftools view -R ~/selscan/output2/${pop}.allchromo.notextr.nsl.bed ~/sfs/vcf/phased.converted.${pop}.del.vcf.gz -Ov -o ~/sfs/vcf/phased.converted.${pop}.del.notsel.vcf ; done
```

Computing and resampling the SFS using python scripts (available in the repository, sfs2.py and subsampling_sfs2.py).
```bash
# Computing the SFS
for pop in brevicaudus lasiotis littoralis mulatta tcheliensis; do  python3 sfs2.py vcf/phased.converted.${pop}.del.sel.vcf ${pop}.sel; done
for pop in brevicaudus lasiotis littoralis mulatta tcheliensis; do  python3 sfs2.py vcf/phased.converted.${pop}.del.notsel.vcf ${pop}.notsel; done
# Resampling the sfs
for pop in lasiotis littoralis mulatta; do  python3 subsampling_sfs2.py output/${pop}.sel.sfs.txt  10 ${pop}.sel.sfs; done
for pop in lasiotis littoralis mulatta; do  python3 subsampling_sfs2.py output/${pop}.notsel.sfs.txt  10 ${pop}.notsel.sfs; done
```
