Prerequisites:
1.FLASH-1.2.6
2.bowtie2-2.1.0
3.PHRAP
Install the three tools and add them to the environment variable of Linux system.
Download dataset from NCBI BioProject PRJNA341910 (https://www.ncbi.nlm.nih.gov/sra?linkname=bioproject_sra_all&from_uid=341910) to this folder.

Pipeline to extract barcode pairs from the cloning vector pool of pBACode-1
#Merge mates of mated reads using FLASH
perl PEflash101.pl configure_pool.txt prepool
#Extract raw paired barcodes
perl barcodeRMvectorSite.pl configure_pool.txt
#Clean sequencing error and hybrid barcode pairs
perl barcodePrepool.pl configure_pool.txt
#Bootstrapping
perl bootstrap_1.pl configure_pool.txt
##Pipeline input: Yeast_BAC_barcodepair_bootstrap_only_<1/2>.fastq, the prefix of input file should be added to the file named configure_pool.txt on the line behind "prepool_reads". Other parameters are already set.
##Pipeline output: barcodePrepoolYeast_BAC_barcodepair_bootstrap_only.txt, format: <left barcode>:<rightbarcode>[TAB]<read count>

Pipeline to extract barcode pairs from BAC library using pBACode-1
#Merge mates of mated reads using FLASH
perl PEflash101.pl configure_pool.txt postpool
#Extract raw paired barcodes
perl barcodeRMvectorSitePost.pl configure_pool.txt
#Clean sequencing error and hybrid barcode pairs
perl barcodePostpool.pl configure_pool.txt
##Pipeline input: Yeast_BAC_barcodepair_<1/2>.fastq, the prefix of input file should be added to the file named configure_pool.txt on the line behind "postpool_reads". Other parameters are already set.
##Pipeline output: barcodePostpoolYeast_BAC_barcodepair.txt, format: <left barcode>:<rightbarcode>[TAB]<read count>

Pipeline to process yeast BAC-PE data
#Merge mates of mated reads using FLASH
perl flashBarcode.pl configure_BAC.txt left
perl flashBarcode.pl configure_BAC.txt right
perl flashBarcode2.pl configure_BAC.txt left
perl flashBarcode2.pl configure_BAC.txt right
#Extract raw barcodes linking end sequences
perl rmVectorBarcodeMergedL.pl configure_BAC.txt
perl rmVectorBarcodeMerged.pl configure_BAC.txt
perl rmVectorBarcode1mateL.pl configure_BAC.txt
perl rmVectorBarcode1mate.pl configure_BAC.txt
#Local assembly
perl groupAssemble.pl configure_BAC.txt left
perl groupAssemble.pl configure_BAC.txt right
#Bowtie index buiding and alignment of end sequences
perl revcomGenome.pl configure_BAC.txt
perl groupLocation.pl configure_BAC.txt left
perl groupLocation.pl configure_BAC.txt right
#Connect end sequences with pooled barcodes
perl poolizeBarcode.pl configure_BAC.txt left
perl poolizeBarcode.pl configure_BAC.txt right
#Remove grouped end sequences that aligned multiply and barcode used more than once
perl groupBarcode.pl configure_BAC.txt left
perl groupBarcode.pl configure_BAC.txt right
#Mark each barcode on genome
perl barcodeLoci.pl configure_BAC.txt left
perl barcodeLoci.pl configure_BAC.txt right
#Pair BAC-PEs
perl S288CchrLen.pl configure_BAC.txt
perl barcodepairDistance.pl configure_BAC.txt
##Pipeline input: Yeast_BAC_PE_<L/R>_<1/2>.fastq, the prefix of input file should be added to the file named configure_BAC.txt on the line behind "left_reads" and "right_reads". Other parameters are already set.
##Pipeline output: 1. barcodepairDistanceYeast_BAC_PE_LYeast_BAC_PE_RbarcodePostpoolYeast_BAC_barcodepair.txt, format: <left barcode>[space]<right barcode>[TAB]<chromosome>[TAB]<starting coordinate>[TAB]<ending coordinate>[TAB]<estimated BAC size>; 2. groupAssemblyYeast_BAC_PE_<L/R>CorrectedUniqueInpool.fa, local assembly result in fasta format, barcode sequence are used as assembly sequence id.

Pipeline to locate BAC by barcode in 3D
#To extract barcodes
perl system.pl filename.txt configure_pooling.txt
#To locate BAC by barcode
perl 3Dcross.pl <threshold, 0.01 or 0.05 recommended>
##Pipeline input: all fastq file names should be add to the file named filename.txt. Other parameters are already set.
##Pipeline output: Yeast_BAC_library_3Dcross_cutoff_<threshold>.out, format: <Clone location>[TAB]<left barcode>
