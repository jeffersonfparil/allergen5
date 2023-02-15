# allergen5
Trying to detect signatures of selection in allergen 5 of Lolium perenne and Lolium rigidum

## Load conda environment
```shell
# conda install -y wget r-base
# conda install -y -c conda-forge julia
# conda install -y -c bioconda mafft

```

## Download the genome annotations and coding DNA sequences of Lolium perenne and Lolium rigidum
```shell
echo 'Lolium_rigidum.gff,https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/022/539/505/GCF_022539505.1_APGP_CSIRO_Lrig_0.1/GCF_022539505.1_APGP_CSIRO_Lrig_0.1_genomic.gff.gz
Lolium_perenne.gff,https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/019/359/855/GCF_019359855.1_MPB_Lper_Kyuss_1697/GCF_019359855.1_MPB_Lper_Kyuss_1697_genomic.gff.gz
Lolium_rigidum.cds,https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/022/539/505/GCF_022539505.1_APGP_CSIRO_Lrig_0.1/GCF_022539505.1_APGP_CSIRO_Lrig_0.1_cds_from_genomic.fna.gz
Lolium_perenne.cds,https://ftp.ncbi.nlm.nih.gov/genomes/refseq/plant/Lolium_perenne/latest_assembly_versions/GCF_019359855.1_MPB_Lper_Kyuss_1697/GCF_019359855.1_MPB_Lper_Kyuss_1697_cds_from_genomic.fna.gz
' > urls.txt

for line in $(cat urls.txt)
do
    spe=$(echo ${line} | cut -d, -f1)
    url=$(echo ${line} | cut -d, -f2)
    wget ${url} -O - | gunzip -c - > ${spe}
done

```

## Extract allergen 5a sequences and merge them into a single fasta file
```shell
grep -i "allergen" *.gff | \
    grep "mRNA" | \
    grep -v "exon" > \
    allergens_mRNA.gff

grep "5a" allergens_mRNA.gff > \
    allergens_mRNA_5a.gff

touch IDs.txt
for i in $(seq 1 $(cat allergens_mRNA_5a.gff | wc -l))
do
    line=$(head -n${i} allergens_mRNA_5a.gff | tail -n1)
    ID=$(echo "$line" | cut -f9 | sed 's/;/\n/g' | grep "gene=LOC" | sed 's/gene=//g')
    species=$(echo ${line} | cut -d':' -f1 | sed 's/.gff//g')
    cds=${species}.cds
    julia src/extract_sequence_using_name_query.jl \
        ${cds} \
        ${ID} \
        ${species}-${ID}.fasta \
        "" \
        true
    echo "$species,$ID" >> IDs.txt
done
### Merge
cat *.fasta > allergens_5a.cds
### Clean-up
rm allergens_mRNA*.gff
rm *.fasta
```

## Align the orthologs and paralogs of allergen 5a
```shell
mafft allergens_5a.cds > allergens_5a.aln
```

## Divide the alignements into pairs
```shell
aln="allergens_5a.aln"
gene="LOC124688403"
julia src/extract_sequence_using_name_query.jl \
        ${aln} \
        ${gene} \
        focal.tmp \
        Lolium_rigidum-${gene} \
        false
for line in $(grep -v ${gene} IDs.txt)
do
    s=$(echo "$line" | cut -d',' -f1)
    g=$(echo "$line" | cut -d',' -f2)
    julia src/extract_sequence_using_name_query.jl \
            ${aln} \
            ${g} \
            other.tmp \
            ${s}-${g} \
            false
    cat focal.tmp other.tmp > ${gene}-${g}.pw.aln
    rm other.tmp
done
rm *.tmp
```

## Divide the pairs into windows
```shell
for f in $(ls *.pw.aln)
do
    # f=$(ls *.pw.aln | head -n1)
    window_size=15
    slide_size=15
    julia src/split_alignment_pairs.jl \
        ${f} \
        ${window_size} \
        ${slide_size} \
        ${f}.split
done

```

## Estimate non-synonymous/synonymous mutaton rates
```shell
# cd src/
# git clone https://github.com/kullrich/kakscalculator2.git
# export PREFIX="."
# cd kakscalculator2;make clean;make
for f in $(ls *.pw.aln.split)
do
    # f=$(ls *.pw.aln.split | head -n1)
    src/KaKs_Calculator \
        -m MS \
        -i ${f} \
        -o ${f}.kaks
    Rscript src/plot_KaKs_across_windows.R \
        ${f}.kaks \
        0.001 \
        allergen-5a
done
```