# ZS-MEM: Zip-Seeding implemented in BWA-MEM
Zip-seeding is a compressive seeding strategy based on reordering-based compression.
ZS-MEM is the incorporation of zip-seeding and BWA-MEM. It takes the reordered reads(FASTQ/FASTA) as input,
and outputs mapping alignments in SAM format.

## Download and install

	git clone https://github.com/icode-xiaohu/ZS-MEM.git
	cd ZS-MEM; make


## Reordering-based compression

In order to conveniently test multiple compressors, ZS-MEM takes 
the decompressed file as input. The state-of-the-art compressors: 
SPRING (https://github.com/shubhamchandak94/Spring), minicom (https://github.com/yuansliu/minicom) and PgRC (https://github.com/kowallus/PgRC) 
prove that reordering-based compression can  facilitate the read mapping. Users can choose any one of them to evaluate the ZS-MEM. 
We provide the example of SPRING because it's a full FASTQ compressor.


**Reordering**

If you download the *master* branch of SPRING, you have to compress and 
decompress to abtain the reordered FASTQ. For example:

	./spring -c -i data.fastq -r -o data.comp
	./spring -d -i data.comp -o spring.fastq
	
Or you can choose the another branch - *reorder-only* (recommended)
which only reorders reads without compressing and decompressing.  

	./spring-reorder -i data.fastq -o spring.fastq

## ZS-MEM mapping
Running ZS-MEM to map reads. Note that the reads file has been reordered by compressors.

	./bwa/bwa index ref.fasta
	./zs-mem ref.fasta spring.fastq > zs_mem.sam
	
Running BWA-MEM for comparison. We add minor modifications to ./bwa/bwa without affecting 
the mapping performance, so that it can output the mapping time profiling. 
Then the speedup of zip-seeding can be evaluated.

	./bwa/bwa mem ref.fasta spring.fatsq > mem.sam
	
There are many existing software to evaluate the mapping consistency of ZS-MEM and BWA-MEM, 
such as wgsim_eval.pl (https://github.com/lh3/wgsim).

## Citations
* Chandak, S., *et al*. (2019) SPRING: a next-generation compressor for FASTQ data, *Bioinformatics*, **35**, 
 2674-2676. [PMID:[30535063][1]].
 
* Liu, Y., *et al*. (2019) Index suffix–prefix overlaps by (w, k)-minimizer to generate long contigs for reads compression, *Bioinformatics*, **35**, 
 2066-2074. [PMID:[30407482][2]].
 
* Kowalski, T.M. and Grabowski, S. (2020) PgRC: pseudogenome-based read compressor, *Bioinformatics*, **36**, 
 2082-2089. [PMID:[31893286][3]].
 
* Li H. (2013) Aligning sequence reads, clone sequences and assembly contigs
 with BWA-MEM. [arXiv:1303.3997v2][4] [q-bio.GN].


[1]: https://www.ncbi.nlm.nih.gov/pubmed/30535063
[2]: https://www.ncbi.nlm.nih.gov/pubmed/30407482
[3]: https://www.ncbi.nlm.nih.gov/pubmed/31893286
[4]: https://arxiv.org/abs/1303.3997
