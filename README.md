# Refinement EVEs annotation in reference genome

Pipeline developped by Elisa Pischedda, while in the Bonizzoni Lab at the University of Pavia (Italy).

## Purpose

The pipeline used in Whitfield et al. 2017, allow to find EVEs within reference genomes. However, output files of the pipeline include results that need to be manually curated. We decide to extend the pipeline of Whitefield with a further reverse blastx step against the database of NR and/or RefSeq and to sostitute the manual curation with an automatic method, avoiding untraceable errors. Optionally the user can manually check the result and modify the scripts to obtaine the desired result.

### Reverse blastx
The input of the pipeline is based on a blastx obtained with one of the following commands. It is necessary to set the output file as descibed in the format parameter:

```sh
blastx -query tophits.fasta \
-db NR/RefSeq_protein \
-evalue 1e-06 \
-outfmt '6 qseqid qstart qend salltitles saccver evalue qframe pident qcovs sstart send slen staxid' \
-out TopHits.blastx
```
```sh
diamond blastx -d nr_diamond.dmnd \
-e 1e-06 \
-f 6 qseqid qstart qend salltitles evalue qframe pident qcovhsp sstart send slen staxids \
-q TopHits.fasta -o TopHits.blastx &
```

- - - -

## Requirements
The pipeline uses the following programs.

* PYTHON 2.7 https://www.python.org/download/releases/
* PYTHON 3 https://www.python.org/download/releases/
* Virus-Host Classifier https://github.com/Kzra/VHost-Classifier
  It requires ETE3 toolkit http://etetoolkit.org/download/
  Note: the first time that it is used, it download the NCBI taxonomy database and put it in /home/username/.etetoolkit , thus the user have to install it in own account
        to update the database please do the following command:

```sh
python3 /ngs-data/Pipelines/Refinement_EVE_annotation_080620/Update_taxon_localdb.py

```

* check the latest version of Taxonkit https://bioinf.shenwei.me/taxonkit/
* update the NCBI taxonomt database for taxonkit with the following command and then uncompress it

```sh
cd /home/username/.taxonkit

wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz

```

* download the file virushostdb.tsv from virus host database in VHost-Classifier-master directory

```sh
cd /ngs-data/Pipelines/Refinement_EVE_annotation_080620/VHost-Classifier-master

wget ftp://ftp.genome.jp/pub/db/virushostdb/virushostdb.tsv

```

## Command line

```sh
bash /ngs-data/Pipelines/Refinement_EVE_annotation_080620/Refinement_EVE_Annotation.sh \
-pipeline_directory /ngs-data/Pipelines/Refinement_EVE_annotation_080620 \
-tool diamond \
-file_blastx /pathTo/TopHits.blastx \
-file_bed_tophit /pathTo/TopHits.bed \
-output_directory /pathTo/Output_directory \
-taxonkit_exe /ngs-data/share_tools/tools/Taxonkit_0.6/taxonkit
```

### Paramters

| parameter name                          | description                                                                             | default |
|-----------------------------------------|-----------------------------------------------------------------------------------------|---------|
| pipeline_directory                      | absolute path of the pipeline directory                                                 |         |
| tool                                    | reverse blastx tool, options: blastx, diamond.                                          | blastx  |
| file_blastx                             | absolute path of the output file of the reverse blast                                   |         |
| file_bed_tophit                         | absolute path of the bed file of the tophits from the pipeline of Whitfield             |         |
| output_directory                        | absolute path of the output directory                                                   |         |
| taxonkit_exe                            | bsolute path of blastn executable                                                       |         |

<br>

## References
[Whitfield et al., 2017]('https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5698160' "Whitfield et al., 2017")

[Kitson et al., 2019]('https://academic.oup.com/bioinformatics/article-abstract/35/19/3867/5368011?redirectedFrom=fulltext' "Kitson et al., 2019")

<br>
```sh
if re.search(r'.*AGAP.*-PA.*', k_hit_db) is None and \
    re.search(r'.*AAEL.*-P.*', k_hit_db) is None and \
    re.search(r'.*ncharacterized.*', k_hit_db) is None and \
    re.search(r'.*PREDICTED.*', k_hit_db) is None and \
    re.search(r'.*hypothetical.*', k_hit_db) is None and \
    re.search(r'.*CLUMA_CG.*', k_hit_db) is None and \
    re.search(r'.*baculoviral.*', k_hit_db) is None and \
    re.search(r'.*LOW QUALITY PROTEIN:.*', k_hit_db) is None and \
    re.search(r'.*unnamed protein product.*', k_hit_db) is None:
    selected_nonviral_evalues.append(v_hit_db.evalue)
    selected_nonviral_elem.append(k_hit_db)
```
<br>
<br>
<br>
