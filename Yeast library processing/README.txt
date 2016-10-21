Prerequisites:
1.FLASH-1.2.6
2.bowtie2-2.1.0
3.PHRAP
4.SOAPdenovo2
5.SSPACE
Install the five tools and add them to the environment variable of Linux system.
Download dataset from NCBI BioProject PRJNA341910 (https://www.ncbi.nlm.nih.gov/sra?linkname=bioproject_sra_all&from_uid=341910) to this folder.

Pipeline to extract barcode pairs from the cloning vector pool of pBACode-1
##Pipeline input: Yeast_BAC_barcodepair_bootstrap_only_<1/2>.fastq: raw read files.
                  configure_pool.txt: configeration file. the prefix of raw read files should be added to the line "prepool_reads". Other parameters have already been set.
##Pipeline output: barcodePrepoolYeast_BAC_barcodepair_bootstrap_only.txt, format: <left barcode>:<rightbarcode>[TAB]<read count>
#Merge mates of mated reads using FLASH
perl PEflash101.pl configure_pool.txt prepool
#Extract raw paired barcodes
perl barcodeRMvectorSite.pl configure_pool.txt
#Clean sequencing error and hybrid barcode pairs
perl barcodePrepool.pl configure_pool.txt
#Bootstrapping
perl bootstrap_1.pl configure_pool.txt

Pipeline to extract barcode pairs from BAC library using pBACode-1
##Pipeline input: Yeast_BAC_barcodepair_<1/2>.fastq: raw read files.
                  configure_pool.txt: configeration file. the prefix of raw read files should be added to the line "postpool_reads". Other parameters have already been set.
##Pipeline output: barcodePostpoolYeast_BAC_barcodepair.txt, format: <left barcode>:<rightbarcode>[TAB]<read count>
#Merge mates of mated reads using FLASH
perl PEflash101.pl configure_pool.txt postpool
#Extract raw paired barcodes
perl barcodeRMvectorSitePost.pl configure_pool.txt
#Clean sequencing error and hybrid barcode pairs
perl barcodePostpool.pl configure_pool.txt

SOAPdenovo2 assembly
#Configuration files for SOAPdenovo2 assembly: soapdenovo_config_file
SOAPdenovo-127mer all -s soapdenovo_config_file -o yeast -K 81 -R -p 8 -L 500

Pipeline to process yeast BAC-PE data
##Pipeline input: Yeast_BAC_PE_<L/R>_<1/2>.fastq: raw read files. 
                  barcodePostpoolYeast_BAC_barcodepair.txt： barcode pair file
                  configure_BAC.txt: configeration file. The prefix of raw read files should be added to the line "left_reads" and "right_reads". The prefix of barcode pair file should be added to the line "pool_file". The prefix of reference sequence file should be added to the line "genome". Other parameters have already been set.
##Pipeline output: barcodepairDistanceYeast_BAC_PE_LYeast_BAC_PE_RbarcodePostpoolYeast_BAC_barcodepair.txt, format: <left barcode>[space]<right barcode>[TAB]<chromosome>[TAB]<starting coordinate>[TAB]<ending coordinate>[TAB]<estimated BAC size>; 
                   groupAssemblyYeast_BAC_PE_<L/R>CorrectedUniqueInpool.fa, local assembly result in fasta format, barcode sequences are used as assembly sequence IDs.
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

SSPACE assembly
SSPACE_Standard_v3.0.pl -l library.txt -s contig_81.scafSeq.fa -T 8 -z 500 -b sspace_out

Pipeline to locate BAC by barcode in 3D
##Pipeline input:  filename.txt. All 36 pooled BAC fastq file names should be add to the file. Other parameters are already set.
##Pipeline output: Yeast_BAC_library_3Dcross_cutoff_<threshold>.out, format: <Clone location>[TAB]<left barcode>
#To extract barcodes
perl system.pl filename.txt configure_pooling.txt
#To locate BAC by barcode
perl 3Dcross.pl <threshold, 0.01 or 0.05 recommended>

Configeration file explanation
1. configure_pool.txt
prepool/postpool_reads: raw read file prefix
prepool/postpool_readlength: raw read readlength
prepool/postpool_fragmentlength: expected readlength after merging
barcode_length_standard: expected barcode length
barcode_length_min: minimal barcode length
barcode_length_max：maximal barcode length
SQ_threshold: minimal sequencing quality requirement of barcode sequence
prepool/postpool_left_upstream: vector sequence upstream of left barcode
prepool/postpool_left_downstream: vector sequence downstream of left barcode
prepool/postpool_right_upstream: vector sequence upstream of right barcode
prepool/postpool_right_downstream: vector sequence downstream of right barcode
site_prepool: endonuclease cut site used for inverse PCR
prepool/postpool_left_anchor: expected vector sequence length on the 5' end of merged reads
prepool/postpool_right_anchor: expected vector sequence length on the 3' end of merged reads
2. configure_BAC.txt
genome: reference sequence prefix
left/right_reads: raw read file prefix
pool_file: prefix of the file that lists barcode pair relation
