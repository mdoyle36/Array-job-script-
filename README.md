# Script for array job

This page contains the script which: 
- Extracted the vcf file for each participant 
- Subset vcf file for recessive EDS genes 
- Ran VEP 111 on the subset vcf files 
- Ran filter_vep on these annotated files 
-filtered for heterozygous variants that passed quality control 

```
#!/bin/bash                                                                                                                                                 
#BSUB -q short                                                                                                                                             
#BSUB -P re_gecip_paediatrics                                                                                                                               
#BSUB -cwd "/re_gecip/paediatrics/mdoyle1/110/array_job_6"
#BSUB -o logs/array_job_6_V6.stdout                                                                                                        
#BSUB -e logs/array_job_6_V6.stderr                                                                                                        
#BSUB -n 1                                                                                                                                                  
#BSUB -R "rusage[mem=2000]"                                                                                                                                
#BSUB -M 2000
#BSUB -J "myArray[1-503]%50"

mkdir -p ./logs

#The file EDS_case_not_solved_GRCh38.txt contains a list of file paths for each participant, which lead to a directory where all genomic files are stored. The next steps allow the extraction of each paricipant's vcf file

cp /re_gecip/paediatrics/mdoyle1/110/array_job_3/EDS_case_not_solved_GRCh38.txt . 
paths=$(cut -f 2 EDS_case_not_solved_GRCh38.txt)
ids=$(cut -f 1 EDS_case_not_solved_GRCh38.txt)

patharr=($paths)
idarr=($ids)

id=${idarr[$LSB_JOBINDEX]}

path=${patharr[$LSB_JOBINDEX]}
path2=$(echo "$path" | tr -d '\r')


vcffile=$(find "$path2" -type f -name "*.vcf.gz" | grep -v -e "SV" -e "genome.vcf" -e "Genotyping" -e "ExpansionHunter")

# The next step was to intersect the files for the coordinates of recessive EDS genes, contained in the recessive_EDS_genes.bed file

module load bedtools/2.31.0 

bedtools intersect -wa -a "$vcffile" -b recessive_EDS_genes.bed > "$id".vcf

# Singularity was then loaded, to allow VEP and filter_vep to run
                                                                                                                                               
module purge
module load singularity/3.8.3

source /re_gecip/paediatrics/mdoyle1/110/array_job_6/vep_111_but_109.conf

export PERL5LIB=$PERL5LIB:/plugins/loftee-GRCh38:/plugins/loftee-GRCh37:/plugins

input="$id".vcf
output="$(basename "${input}" .vcf)"
echo "${output}"

# VEP 111 was then run on the participants' vcf files                                                                                                                                                                   
singularity run --bind "${MOUNT_WD}","${MOUNT_GENOMES}","${MOUNT_GEL_DATA_RESOURCES}","${MOUNT_PUBLIC_DATA_RESOURCES}","${MOUNT_SCRATCH}" "${IMG}" vep \
    --cache \
    --species homo_sapiens \
    --format vcf \
    --dir_plugins /plugins/ \
    --offline --vcf \
    --cache \
    --force_overwrite \
    --assembly GRCh38 \
    --cache_version 111 \
    --dir_cache "${CACHE}" \
    --fasta "${REFFASTA}" \
    --input_file "${input}" \
    --warning_file "${output}"_errors.vcf \
    --output_file "${output}"_annotated.vcf \
    --plugin dbNSFP,"${DB_NSFP_DB}","${DB_NSFP_REPLACEMENT_LOGIC}",ALL \
    --plugin LoF,loftee_path:/plugins/loftee-GRCh38,human_ancestor_fa:"${LOFTEE38HA}",gerp_bigwig:"${LOFTEE38GERP}",conservation_file:"${LOFTEE38SQL}" \
    --plugin SpliceAI,snv="${SPLICEAIRAW38}",indel="${SPLICEAIINDEL38}" \
    --plugin CADD,"${CADD16}" \
    --plugin REVEL,"/public_data_resources/vep_resources/REVEL/revel_v1.3_GRCh38.tsv.gz" \
    --plugin mutfunc,db="${MUTFUNC_DB}"


# The command filter_vep was then used to filter for high or moderate impact terms, low allele frequencies, and high scores for either CADD Phred or REVEL                                                                                                                                                
singularity exec --bind "${MOUNT_WD}","${MOUNT_GENOMES}","${MOUNT_GEL_DATA_RESOURCES}","${MOUNT_PUBLIC_DATA_RESOURCES}","${MOUNT_SCRATCH}" "${IMG}" filter_vep -i "${output}"_annotated.vcf --filter "(CSQT matches transcript_ablation or CSQT matches splice_acceptor_variant or CSQT matches splice_donor_variant or CSQT matches stop_gained or CSQT matches frameshift_variant or CSQT matches stop_lost or CSQT matches start_lost or CSQT matches transcript_amplification or CSQT matches feature_elongation or CSQT matches feature_truncation or CSQT matches inframe_insertion or CSQT matches inframe_deletion or CSQT matches missense_variant or CSQT matches protein_altering_variant) and (AF < 0.01 or not AF)" --filter "(CADD_PHRED > 30 or REVEL > 0.80)" > "$id"_filtered_variants.vcf

# Variants were then filtered so that only those that were heterozygous and PASS the filter were retained

module load vcftools/0.1.16
module load bcftools/1.16
module load tabix/1.18

for file in "$id"_filtered_variants.vcf; do

    # Extract the filename without the extension                                                                                                                                  
    filename=$(basename "$file" .vcf)

    # Filtering for variants that pass filter and are heterozygous                                                                                                                
    vcftools --vcf "$file" --out "$filename" --recode --recode-INFO-all --keep-filtered PASS
    bgzip "$filename".recode.vcf
    tabix -p vcf "$filename".recode.vcf.gz
    bcftools view -g het "$filename".recode.vcf.gz -o "$filename"_pass_filter_heterozygous.vcf
    rm -f *recode.vcf.gz *.recode.vcf.gz.tbi
    
done

# Remove empty files

rm -f "$id"_annotated.vcf -f "$id"_annotated.vcf_summary.html -f "$id"_errors.vcf -f "$id".vcf

if ! grep -q -v "^#" "$id"_filtered_variants.vcf; then
            # If the file does not contain any variants, remove it                                                                                                                
            rm -f "$id"_filtered_variants.vcf
        fi
done

if ! grep -q -v "^#" "$filename"_pass_filter_heterozygous.vcf; then
            # If the file does not contain any variants, remove it                                                                                                                                                                                                                                                                                               
            rm -f "$filename"_pass_filter_heterozygous.vcf
        fi
done                                                                                                                                                                      
        
