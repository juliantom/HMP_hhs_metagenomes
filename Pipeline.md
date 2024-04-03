### Obtaining the necessary files
1. Access the (HMP website)[https://portal.hmpdacc.org/] and download the **manifest** and **metadata** files. <br>
Criteria (all sites):
    * Click on *Data*
    * On *Samples* TAB
        * Select *HHS* (Healthy Human Study) from *Studies* category.
    * On *Files* TAB
        * Select *FASTQ* (fastQ format - raw sequencing) from *Format* category.
        * Select *wgs_raw_seq_set* (metagenomes) from *Type* category.
    * Click *Add all files to the Cart*
    * Go to Cart (top right corner)
    * Download *File Manifest* and *Sample Metadata* locally.
- Note: Human sequencing reads (data) have bee removed from the metagenomic samples prior to upload.
2. Prepare working folder
```bash
# Create a working directory (mine looks as)
my_date=$(date +'%Y_%m_%d')
mkdir HMP_HHS_${my_date} && cd HMP_HHS_${my_date}
```
3. Edit manifest and metadata file
```bash
# Set variables
my_date=2023_07_07 # date - should have the date the metadata and manifest files were downloaded
# Folders
dir_original=01-Orig_manifest_metadata
dir_new_files=02-Edit_manifest_metadata
# Files
original_manifest=hmp_manifest_35005f35fb.tsv # ajust as needed
original_metadata=hmp_manifest_metadata_a0061aaf77.tsv # ajust as needed
new_manifest="manifest-$my_date"
new_metadata="metadata-$my_date"

mkdir $dir_new_files

# Copy files from local to server 
scp ~/Downloads/hmp_manifest_35005f35fb.tsv server_name@ip:/path/to/HMP_HHS_${my_date}/$dir_new_files
scp ~/Downloads/hmp_manifest_metadata_a0061aaf77.tsv server_name@ip:/path/to/HMP_HHS_${my_date}/$dir_new_files

# Copy files and edit name
cp $dir_original/$original_manifest $dir_new_files/$new_manifest.tsv
cp $dir_original/$original_metadata $dir_new_files/$new_metadata.tsv

# Remove 454 sequencing entires and duplicated entries
### Manifest
grep -v '_454' $dir_new_files/$new_manifest.tsv | awk 'BEGIN{FS=OFS="\t"}NR==1{print $0}NR>1{print $0 | "sort | uniq "}' | sed -e 's/ /_/g' > $dir_new_files/$new_manifest-sort_unique.tsv
### Metadata
awk 'BEGIN{FS=OFS="\t"}NR==1{print $0}NR>1{print $0 | "sort | uniq "}' $dir_new_files/$new_metadata.tsv | sed -e 's/ /_/g' > $dir_new_files/$new_metadata-sort_unique.tsv
```
4. Download metagenomes (subset of 50 metagenomes with most data - biggest file size - for external naris and feces)
```bash
# Copy file
cp $dir_new_files/$new_manifest-sort_unique.tsv 98-DATA/02-dl_metag-temp_1.txt
# Generate ID for each sample
cat 98-DATA/02-dl_metag-temp_1.txt | awk 'BEGIN{FS=OFS="\t"}NR>1{print $4}' | awk -F',' '{print $1}' | sed -e 's/\.tar\.bz2//' |awk -F'/' 'NR==1{print "Internal_ID"}{print $12"_"$11"\.tar\.bz2"}' | sed -e 's/\.bz2//' > 98-DATA/02-dl_metag-temp_2.txt
# Combine original manifest with new Internal IDs
paste 98-DATA/02-dl_metag-temp_1.txt 98-DATA/02-dl_metag-temp_2.txt > 98-DATA/04-dl_metag.txt
# Generate a two-sample test file
head -n3 98-DATA/04-dl_metag.txt > 98-DATA/03-dl_metag-test.txt
# Run test
./99-SCRIPTS/01-download_hmp-v2.py -i 98-DATA/03-dl_metag-test.txt -d 03-Download_test/ -c -r 2 --report-file 98-DATA/03-dl_metag-test-report.txt
```
4. Check download and retry download if INCOMPLETE or MISSING <br>
Total samples: 3183 <br>
454 seq samples (removed): -15 <br>
Total Illumina samples: 31 <br>
Samples not found in https: 4<br>
Total samples: 1364<br>
Missing samples: SRS013216_v2,SRS016990_v2,SRS022093_v2,SRS022093_v2 <br>
```bash
nohup ./99-SCRIPTS/01-download_hmp-v2.py -i 98-DATA/04-dl_metag.txt -d 04-Download_mg/ -c -r 10 --report-file 98-DATA/04-dl_metag-report.txt &
mv nohup.out nohup.out.v1.txt
################# Check for errors
# Remove all lines before SRS013216 (already downloaded)
cat 98-DATA/04-dl_metag.txt | awk 'BEGIN{FS=OFS="\t"}NR==1{print $0}NR>1{print $0 | "sed '1,/SRS013216/d' "}' > 98-DATA/04-dl_metag-v2.txt
nohup ./99-SCRIPTS/01-download_hmp-v2.py -i 98-DATA/04-dl_metag-v2.txt -d 04-Download_mg/ -c -r 10 --report-file 98-DATA/04-dl_metag-report-v2.txt &
#################
################# Check for errors
# Remove all lines before SRS016990 (already downloaded)
mv nohup.out nohup.out.v2.txt
cat 98-DATA/04-dl_metag.txt | awk 'BEGIN{FS=OFS="\t"}NR==1{print $0}NR>1{print $0 | "sed '1,/SRS016990/d' "}' > 98-DATA/04-dl_metag-v3.txt
nohup ./99-SCRIPTS/01-download_hmp-v2.py -i 98-DATA/04-dl_metag-v3.txt -d 04-Download_mg/ -c -r 10 --report-file 98-DATA/04-dl_metag-report-v3.txt &
#################

################# Check for errors
# Remove all lines before SRS022093 (already downloaded)
mv nohup.out nohup.out.v3.txt
cat 98-DATA/04-dl_metag.txt | awk 'BEGIN{FS=OFS="\t"}NR==1{print $0}NR>1{print $0 | "sed '1,/SRS022093/d' "}' > 98-DATA/04-dl_metag-v4.txt
nohup ./99-SCRIPTS/01-download_hmp-v2.py -i 98-DATA/04-dl_metag-v4.txt -d 04-Download_mg/ -c -r 10 --report-file 98-DATA/04-dl_metag-report-v4.txt &
#################

################# Check for errors
# Remove all lines before SRS078182 (already downloaded)
mv nohup.out nohup.out.v4.txt
cat 98-DATA/04-dl_metag.txt | awk 'BEGIN{FS=OFS="\t"}NR==1{print $0}NR>1{print $0 | "sed '1,/SRS078182/d' "}' > 98-DATA/04-dl_metag-v5.txt
nohup ./99-SCRIPTS/01-download_hmp-v2.py -i 98-DATA/04-dl_metag-v5.txt -d 04-Download_mg/ -c -r 10 --report-file 98-DATA/04-dl_metag-report-v5.txt &
mv nohup.out nohup.out.v5.txt
#################
```
5. Rename files. Substitute HMP version with visit number.
```bash
# Create a list with names of files excluding extensions and path
ls 04-Download_mg/ | sed -e 's/\.tar\.bz2//' > 98-DATA/05-metagenomeID-hmp_to_visit-1.txt

# Create new IDs for all samples (visit number, HMP version, and Biosample ID)
# Visit_number (v1,v2,v3)
# HMP version (hmp1, hmp2) - Resampling?
# Biosample ID (SRS...)
for sample_name in `cat 98-DATA/05-metagenomeID-hmp_to_visit-1.txt `
do
my_sample_id="$sample_name"
sample_search=$( echo "$my_sample_id" | awk -F'_' '{print $2"/"$1}' )
original_sampl_id=$( grep "$sample_search" 02-Edit_manifest_metadata/manifest-2023_07_07.tsv | awk -F'\t' '{print $NF}' )
visit_number=$( grep "$original_sampl_id" 02-Edit_manifest_metadata/metadata-2023_07_07-sort_unique.tsv | awk -F'\t' '{print $5}' )
visit_id=$( echo "v$visit_number" )
biosample_id=$( echo "$my_sample_id" | awk -F'_' '{print $1}' )
hmp_version=$( echo "$my_sample_id" | awk -F'_' '{print "hmp"$2}' | sed -e 's/v//' )
echo -e "${my_sample_id}\t${visit_id}_${hmp_version}_${biosample_id}" >> 98-DATA/05-metagenomeID-hmp_to_visit-2.txt
done

```
6. Uncompress files
```bash
# Create directory for uncompressed files
mkdir 05-Uncompressed_metagenomes

# Create directory for renamed files
mkdir 06-metagenomes_raw

# Extract files (TEST)
./99-SCRIPTS/02-extract_hmp_data.py -i 98-DATA/05-metagenomeID-hmp_to_visit-2-test.txt --in-dir 04-TEST/ --out-dir 05-Uncompressed_metagenomes/ --report report.tsv
# Extract files
nohup ./99-SCRIPTS/02-extract_hmp_data.py -i 98-DATA/05-metagenomeID-hmp_to_visit-2.txt --in-dir 04-Download_mg/ --out-dir 05-Uncompressed_metagenomes/ --report 98-DATA/05-extraction_report.tsv &
mv nohup.out nohup.out.v6.txt
```

6. Make copy of files and add isolation site
```bash
mkdir 06-metagenomes_raw/
# Store path to FASTQ files for each sample (only R1 & R2)
cat 98-DATA/05-extraction_report.tsv | awk 'NR>1{print $2"|"$3}' | rev | sed -e 's/,/|/' | rev | awk -F'|' -v OFS="\t" '{print $1,$NF}' | sed -e 's/\.1\.fastq//' | sed -e 's/\.2\.fastq//' | sed -e 's/\.singleton\.fastq//' | awk 'BEGIN{FS=OFS="\t"}{print $1,$2}' | sed -e 's/\.\///g' > 98-DATA/06-path_to_raw_fastq.txt

# Rapid test of the follwing script
#head -n3 98-DATA/06-path_to_raw_fastq.txt > 98-DATA/06-path_to_raw_fastq-test.txt

# Move FASTQ files to new folder (this will exclude singletons)
while IFS= read -r line
do
metagenome_id=$( echo "$line" | awk '{print $1}' )
path_original_fastq=$( echo "$line" | awk '{print $2}' )
mv 05-Uncompressed_metagenomes/$metagenome_id/$path_original_fastq.1.fastq 06-metagenomes_raw/${metagenome_id}_R1.fastq
mv 05-Uncompressed_metagenomes/$metagenome_id/$path_original_fastq.2.fastq 06-metagenomes_raw/${metagenome_id}_R2.fastq
done < 98-DATA/06-path_to_raw_fastq.txt

# Verify only SINGLETONS are remaining in the extracted files
# Number of row should be the same
find 05-Uncompressed_metagenomes -type f -name "*.fastq" > 98-DATA/06-find_singletonts.txt
wc -l 98-DATA/06-find_singletonts.txt
grep 'singleton' 98-DATA/06-find_singletonts.txt | wc -l 
# Delet all subfolders and singleton files
rm -r 05-Uncompressed_metagenomes/*

```

7. Get (fetch) subsite from NCBI Biosample
```bash
# Obtain isolation site for each sample
# 1. Create a list with metagenome IDs
cat 98-DATA/06-path_to_raw_fastq.txt |  awk '{print $1}' > 98-DATA/07-list_metagenomeID-v1.txt
# 2. Find isolation site from BioSample repository using NCBI E_utilities (E-ENTREZ)
for sample_id in `cat 98-DATA/07-list_metagenomeID-v1.txt`
do
metagenome=$( echo "$sample_id" | awk -F'_' '{print $NF}' )
hmp_code=$( echo "$sample_id" | awk 'BEGIN{FS=OFS="_"}{print $1,$2}' )
isolation_source=$( esearch -db biosample -query "$metagenome" | efetch -db biosample | grep "isolation source" | sed 's/.*\=//' | sed -e 's/"//g' | sed -e 's/ /_/' | grep 'G_DNA' | sed -s 's/G_DNA_//' )
# If variable is empty try second option
if [ -z $isolation_source ]
then
isolation_source=$( esearch -db biosample -query "$metagenome" | efetch -db biosample | grep 'Identifiers' | tr ' ' '\n' | grep 'G_DNA' | sed -e 's/;//' | sed -e 's/.*G_DNA_//' )
fi
# Define an associative array to store word abbreviations
declare -A word_abbreviations=(
    ["Anterior_nares"]="S_AN"
    ["Attached/Keratinized_gingiva"]="O_KG"
    ["Buccal_mucosa"]="O_BM"
    ["Hard_palate"]="O_HP"
    ["L_Retroauricular_crease"]="S_RCL"
    ["Mid_vagina"]="V_MV"
    ["Palatine_Tonsils"]="O_PT"
    ["Posterior_fornix"]="V_PF"
    ["R_Antecubital_fossa"]="S_AFR"
    ["R_Retroauricular_crease"]="S_RCR"
    ["Saliva"]="O_SV"
    ["Stool"]="G_ST"
    ["Subgingival_plaque"]="O_SUBP"
    ["Supragingival_plaque"]="O_SUPP"
    ["Throat"]="O_TH"
    ["Tongue_dorsum"]="O_TD"
    ["Vaginal_introitus"]="V_VI"
)
# Check if empty
if [ -z $isolation_source ]
then
    isolation_source=EMPTY
    site_abbreviation=EMPTY
else 
    # Find abbreviation
    if [ -n "${word_abbreviations[$isolation_source]}" ]
    then
        site_abbreviation="${word_abbreviations[$isolation_source]}"
    else
        site_abbreviation=NOT_FOUND
    fi
fi
# Print 
echo -e "$sample_id\t$metagenome\t$isolation_source\t$site_abbreviation\t${site_abbreviation}_${metagenome}_${hmp_code}" >> 98-DATA/07-list_metagenomeID-v2-site.txt
done

# Make copy of original result (only FIRST round)
mv 98-DATA/07-list_metagenomeID-v1.txt 98-DATA/07-list_metagenomeID-v0.txt

# Check for empty fields/errors during http fetching
cat 98-DATA/07-list_metagenomeID-v2-site.txt | grep "EMPTY" | cut -f1 > 98-DATA/07-list_metagenomeID-v1.txt
cat 98-DATA/07-list_metagenomeID-v2-site.txt | grep -v 'EMPTY' > 98-DATA/07-list_metagenomeID-v2-site-temp.txt
mv 98-DATA/07-list_metagenomeID-v2-site-temp.txt 98-DATA/07-list_metagenomeID-v2-site.txt

# Redo fetching until all samples have a sampling site

# Rename file to full list of samples
mv 98-DATA/07-list_metagenomeID-v0.txt 98-DATA/07-list_metagenomeID-v1.txt
```
8. Prepare files for QC (illumina-utils + minoche det al. recommendations)
```bash
# Create directory for renaming metagenome
mkdir 07-metagenomes_raw_site

 98-DATA/07-list_metagenomeID-v2-site.txt

# Rename file adding the human subsite
while IFS= read -r line
do
old_metagenome_id=$( echo $line | awk '{print $1}' )
new_metagenome_id=$( echo $line | awk '{print $5}' )
path_old_metagenomes=06-metagenomes_raw/$old_metagenome_id
path_new_metagenomes=07-metagenomes_raw_site/$new_metagenome_id
mv ${path_old_metagenomes}_R1.fastq ${path_new_metagenomes}_R1.fastq
mv ${path_old_metagenomes}_R2.fastq ${path_new_metagenomes}_R2.fastq
done < 98-DATA/07-list_metagenomeID-v2-site.txt
```

9. Perform QC (Minoche et al. implementation with illumina-utils)
```bash
mkdir 08-qc_metagenomes
# activate anvio (docker)
docker exec -it anvio7_1_main bash
cd /data/JULIAN/DATABASES/HMP_HHS_2023_07_07
# make list of samples
cat 98-DATA/07-list_metagenomeID-v2-site.txt | cut -f5 > metagenomes_list-v1-3164.txt
# Create config ini file for each metagenome
for sample in ` cat metagenomes_list-v1-3164.txt `
do
raw_dir="$PWD/07-metagenomes_raw_site"
qc_dir="$PWD/08-qc_metagenomes"
cp 98-DATA/08-ini-template.txt $qc_dir/$sample.ini
sed -i "s/METAGENOME/$sample/" $qc_dir/$sample.ini
sed -i "s|RAW_DIR|$raw_dir|" $qc_dir/$sample.ini
sed -i "s|QC_DIR|$qc_dir|" $qc_dir/$sample.ini
done

# Run QC on metagenomes

./99-SCRIPTS/s-qc_compress_metagenomes.sh
```

10. Rename files and create quality matrix
```bash
# Rename INI files
./99-SCRIPTS/s-rename_files.sh 08-qc_metagenomes/ 08-qc_metagenomes/ ini
# Rename RAW metagenomes
./99-SCRIPTS/s-rename_files.sh 07-metagenomes_raw_site/ 07-metagenomes_raw_site/ fastq
# Rename QC metagenomes
./99-SCRIPTS/s-rename_files.sh 08-qc_metagenomes/ 09-qc_metagenomes_renamed/ gz
# Rename STATS files
./99-SCRIPTS/s-rename_files.sh 08-qc_metagenomes/ 09-qc_metagenomes_renamed/ txt
# Build QC matrix
# Create an updated list of metagenomes IDs
ls 09-qc_metagenomes_renamed/*.gz | sed -e 's/09-qc_metagenomes_renamed\///' | sed -e 's/-QUALITY.*//'| sort -u > metagenomes_list-v2-3164.txt

# header
echo -e "Internal_id\tBioSample\tSite\tSubsite\tSubject_gender\tTotal_reads_raw\tTotal_reads_QC\tFraction_pairs_pass_QC\tMin_1M_paired_end_read_compliance\tMin_1K_paired_end_read_compliance" > hmp_qc_info-m_3164-231016.tsv
# Add qc info
for metagenome in `cat metagenomes_list-v2-3164.txt`
do
internal_id=$metagenome
biosample_id=$( echo "$internal_id" | awk -F'_' '{print $3}' )
site_id=$( echo "$internal_id" | awk -F'_' '{print $1}' )
subregion_id=$( echo "$internal_id" | awk -F'_' '{print $2}' )
hmp_samp_id=$( grep "$biosample_id" 02-Edit_manifest_metadata/manifest-2023_07_07-sort_unique.tsv | awk -F '\t' '{print $5}' | sort | uniq )
subject_gender=$( grep "$hmp_samp_id" 02-Edit_manifest_metadata/metadata-2023_07_07-sort_unique.tsv | awk -F '\t' '{print $6}' )
raw_reads_total=$( cat 09-qc_metagenomes_renamed/${internal_id}-STATS.txt | grep 'number of pairs analyzed' | awk -F':' '{print $2}' | sed -e 's/^ //' | awk '{print $1}' )
qc_reads_total=$( cat 09-qc_metagenomes_renamed/${internal_id}-STATS.txt | grep 'total pairs passed' | awk -F':' '{print $2}' | sed -e 's/^ //' | awk '{print $1}' )
fraction_pairs_QC=$(awk -v num="$qc_reads_total" -v den="$raw_reads_total" 'BEGIN {if (den == 0 || den == "") print "ERROR"; else printf "%.1f\n", (num / den) * 100}')
pass_min_total_reads_1M=$( echo "$qc_reads_total" | awk '{if ($1>=1000000)print "PASS"; else print "FAIL"}' )
pass_min_total_reads_1k=$( echo "$qc_reads_total" | awk '{if ($1>=100000)print "PASS"; else print "FAIL"}' )
echo -e "$internal_id\t$biosample_id\t$site_id\t$subregion_id\t$subject_gender\t$raw_reads_total\t$qc_reads_total\t$fraction_pairs_QC\t$pass_min_total_reads_1M\t$pass_min_total_reads_1k" >> hmp_qc_info-m_3164-231016.tsv
done
```
10. Extract oral metagenomes with at least 1M paired-end reads (0.5M R1 + 0.5M R2)
```bash
# List of HMP QC oral metagenomes (1M paired-end reads)
# A 1M read sample should have enough reads (10,000 of 100 bases long) to cover 1X at least 50% of a genome (2Mbp) that is 1% abundant.
cat hmp_qc_info-m_3164-231016.tsv | grep "^O_" | grep 'hmp2'| awk '{if($9=="PASS") print $1 | "sort -u"}' > metagenome-hmp_qc-oral_1218-list-231016.txt

head -n1 hmp_qc_info-m_3164-231016.tsv > metagenome-hmp_qc-oral_1218-metadata-231016.txt
cat hmp_qc_info-m_3164-231016.tsv | grep "^O_" | grep 'hmp2'| awk '{if($9=="PASS") print $0 | "sort -nk1,1"}' >> metagenome-hmp_qc-oral_1218-metadata-231016.txt
```

```bash
path_to_perio_files="/data/JULIAN/CALIFF_PERIO/califf_perio_metagenomes_QC_PERIO"
ls $path_to_perio_files | sed -e 's/_PERIO.*//' | sort | uniq > metagenomes-perio_calif-list.txt

for perio_sample in `cat metagenomes-perio_calif-list.txt`
do
cp ${path_to_perio_files}/${perio_sample}_PERIO_R1.fastq ${path_to_perio_files}/O_PERIO_${perio_sample}_R1.fastq
cp ${path_to_perio_files}/${perio_sample}_PERIO_R2.fastq ${path_to_perio_files}/O_PERIO_${perio_sample}_R2.fastq
done

for perio_sample in `cat metagenomes-perio_calif-list.txt`
do
rm ${path_to_perio_files}/${perio_sample}_PERIO_R1.fastq
rm ${path_to_perio_files}/${perio_sample}_PERIO_R2.fastq
done

```
Clearn
```bash
# remove unnecessary folders
rm -r 03-Download_test/ 04-Download_mg/ 04-TEST/ 05-Uncompressed_metagenomes/ 06-metagenomes_raw/

```
```bash
for perio_samp in `cat /data/JULIAN/CALIFF_PERIO/samples_id-QC-CALIFF.txt`
do
int_id="O_PERIO_${perio_samp}"
bio_samp="${perio_samp}"
site_samp="O"
subsite_samp="PERIO"
sub_gend="NA"
reads_raw=$( cat /data/JULIAN/CALIFF_PERIO/califf_perio_metagenomes_QC/${perio_samp}-STATS.txt | grep "number of pairs analyzed" | awk -F":" '{print $2}' | sed -e 's/ //g' )
reads_qc=$( cat /data/JULIAN/CALIFF_PERIO/califf_perio_metagenomes_QC/${perio_samp}-STATS.txt | grep "number of pairs analyzed" | awk -F":" '{print $2}' | sed -e 's/ //g' | sed -e 's/(.*//' )
frac_reads_qc=$(awk -v reads_raw="$reads_raw" -v reads_qc="$reads_qc" 'BEGIN { result = (reads_qc / reads_raw) * 100; printf "%.1f\n", result }')
test_1M=$( awk -v value="$reads_qc" 'BEGIN { if (value >= 1000000) print "PASS"; else print "FAIL" }' )
test_1K=$( awk -v value="$reads_qc" 'BEGIN { if (value >= 100000) print "PASS"; else print "FAIL" }' )
echo -e "$int_id\t$bio_samp\t$site_samp\t$subsite_samp\t$sub_gend\t$reads_raw\t$reads_qc\t$frac_reads_qc\t$test_1M\t$test_1K" >> metagenome-hmp_qc-oral_1218-metadata-231016-copy.txt
done
```