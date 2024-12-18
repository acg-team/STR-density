# STR-density
Scripts and background to determine short tandem repeat (STR) densities in genes and genomes

## Background 
### The structure of a (eukaryotic) gene
Genes contain the information necessary for living cells to create proteins.
Read [this wikiversity article](https://en.wikiversity.org/wiki/WikiJournal_of_Medicine/Eukaryotic_and_prokaryotic_gene_structure) to find out more about how genes are made up of exons, introns, and more.

![](images/gene_structure.png)

There are around 20'000 protein coding genes in the human genome.
What may be surprising is that there are many more proteins than there are genes, with estimates of the number of human proteins ranging between 80'000 and 400'000.
This means that a single gene has to be able to produce several different protein products.
This is possible through a process known as alternative splicing.
The different protein products generated from the same gene are known as isoforms.
[Wikipedia](https://en.wikipedia.org/wiki/Alternative_splicing) has a good overview of alternative splicing.

### Short tandem repeats
Short tandem repeats (STRs) are consecutive repetitions of 1-6 basepair (bp) motifs. 
They are estimated to make up around 3% of the human genome. Below is an example STR locus:  

![](images/str_example.png)  

STRs are a rich source of genomic variation, with some loci having mutation rates up to a 10000 times higher than point mutations. 
STR mutations are typically the result of DNA polymerase slippage during replication, where strand misalignment after polymerase detachment results in the insertion or deletion of one or more repeat units at an STR locus:  

![](images/str_slippage_example.png)

For more background on STRs, the way they mutate, and why they are relevant, you can have a look at our recent review: [Mutation and selection processes regulating short tandem repeats give rise to genetic and phenotypic diversity across species](https://onlinelibrary.wiley.com/doi/full/10.1111/jeb.14106) (there's also a pdf in the [literature](literature) folder).

## Project goal
We are interested in determining whether the evolutionary age of genes (where they first emerged in the tree of life) has a relationship to the number of STR loci. 
Specifically, the goal of this project is to write software that calculates the number and density of STRs in genes. 
This would provide clues to understanding the mechanisms behind how genes emerge.
Ideally, the software should provide specific statistics for intronic, exonic, and/or coding regions.

In order to write software to do this, it is important to understand how genetic data is represented computationally, and how we can work with such data. 
The most important formats to be aware of initially are [FASTA](https://www.ncbi.nlm.nih.gov/genbank/fastaformat/), [BED](https://samtools.github.io/hts-specs/BEDv1.pdf), and [GFF/GTF](https://www.gencodegenes.org/pages/data_format.html). 
For each of these, different libraries exist that can be used to parse them.

## Useful tools and libraries
### General
* [GitHub](https://github.com/): 
Source control is critical in any non-trivial project.
* [miniconda](https://docs.anaconda.com/miniconda/): 
Super useful for managing and separating work environments. 
Is used mainly for setting up Python environments, but can also be used for R.
* [bedtools](https://bedtools.readthedocs.io/en/latest/index.html): 
"...a swiss-army knife of tools for a wide-range of genomics analysis tasks.". 
Useful command line tools for working with bed files, fasta files, etc.
### Python
* [Biopython](https://biopython.org/): 
Library for biological computation in Python. 
The [SeqIO](https://biopython.org/wiki/SeqIO) module can be useful to read and write sequences (fasta, fastq, etc.), although I tend to prefer pysam these days (see below).
* [pysam](https://niyunyun-pysam-fork.readthedocs.io/en/latest/api.html):
Python wrapper for the original [htslib](https://github.com/samtools/htslib) C library. 
Focused on reading/writing alignment files, but also supports reading and writing of fasta, vcf, and more.
* [gtftools](http://www.genemine.org/gtftools.php):
Software that can parse GTF files and provides functionality useful for this project.
Importantly, gtftools can merge the exonic regions for different isoforms of a gene.
This will be useful when determining the overlap between STRs, genes, and their exons.
There is a paper describing gtftools in the [literature](literature) folder.

## 2024-10-23: taskset 1
Goal: use some test files located (in [data](./data/)) to get familiar with the problem we are trying to address.
A good starting point will be to determine the number of STR positions across all STRs specified in [this BED file](./data/test.bed), overlap the gene region of the SAMD11 gene.
The coordinates of the SAMD11 gene are specified in [this GTF file](./data/test.gtf).
According to my tests, there are 68 nucleotides in the entire SAMD11 gene sequence that consist of STRs specified in the `test.bed` file.
If you manage to get to the same value, you can try the next challenge, which will be to determine the amount of overlap between STRs in `test.bed` and the exonic regions of the SAMD11 gene.
According to my tests, the answer here should be 24 nucleotides.
I will leave it up to you to figure out how you want to tackle this, but it may be useful to merge all overlapping exons from the different SAMD11 transcripts using.
This could be done using [gtftools](http://www.genemine.org/gtftools.php), for example, although you will need to figure out how to install it.

The most important part: don't get too frustrated or discouraged if something doesn't work.
It is expected that it will not work perfectly the first time.
This will be an iterative process where we figure out the best way to do this.
Don't hestitate to reach out if you have any questions or are stuck.

## 2024-11-08: taskset 2
**Goal 1:** create a standalone python script that calculates STR densities in genes and has a command-line interface (CLI).
To implement the CLI, you can use python's [argparse](https://docs.python.org/3/library/argparse.html) library.
There are many tutorials online that show how to use it.
The CLI should take at least two arguments: one that specifies the filepath to a GTF file with gene annotations, and one that specifies the filepath to a BED file that contains STR regions.
The script should then calculate the STR density of each gene and return this to the user.

**Goal 2:** create a folder in this GitHub repository called `scripts` that we can use to store our program in. 
The script that you create in `Goal 1` should be stored here, for example.
To learn more about Git, take a look at [this page](https://github.com/git-guides).
For now, you should propably only use the [Add](https://github.com/git-guides/git-add), [Commit](https://github.com/git-guides/git-commit), and [Push](https://github.com/git-guides/git-push) commands.
You can use these commands either from the command line or through VSCode (again, there are many tutorials online for this).

## 2024-11-13: taskset 3
Good job on implementing the script!
I made a few small changes to the `calculate_density()` function so it runs a bit faster.
I also propose some changes to the format in which the results are returned by the script:

1)  Currently, the output is just written to stdout using `print()`, please also add some functionality to write the results to an output file specified via a `--output` commmand line argument.
2)  Please return the output in CSV format so it will be easier to work with in the future. 
    You can make a pandas DataFrame that contains the results and then use `pd.DataFrame.to_csv()` to write it to a file.
3)  Currently, the script only returns the `gene_name`. 
    Please also return the `gene_id` as well.

Once you've implemented these changes, you can run the script on the full dataset.
You can find the BED file we'll use [here](https://github.com/acg-team/ConSTRain/tree/main/resources/h_sapiens) (the file you need is called 'GRCh38_repeats.bed.gz').
You can download it and then uncompress it at the command line using `gunzip GRCh38_repeats.bed.gz`.
You can download a complete GTF file of the human genome directly from [GENCODE](https://www.gencodegenes.org/human/#) (download the 'Basic gene annotation' GTF file for the 'CHR' regions).
Get the files, run the script, and store the output.
Then we can get started on the analysis!

## 2024-11-19: taskset 4
By now you should have the STR density of all human genes in the GENCODE GTF file.
This is the point where we switch from pure coding to investigating biological questions.
The main motivation for this project was to determine if there is a link between the evolutionary age of a gene and it's STR content.
To find the evolutionary ages of the genes in our dataset, you can use [GenOrigen](http://chenzxlab.hzau.edu.cn/GenOrigin/#!/).
Specifically, [this CSV file](http://chenzxlab.hzau.edu.cn/static/GenOrigin/download/Age_Inference/Homo_sapiens.csv) contains the gene ages for human protein coding genes.
You can map the gene ages contained in this CSV to the genes you already analysed by matching on the `ensembl_gene_id`.

*Hint: you might have to remove the ENSEMBL version numbers from the gene ids in your dataset first. For example, ENSG00000290825.1 -> ENSG00000290825*