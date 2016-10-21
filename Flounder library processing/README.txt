Prerequisites:
1.FLASH-1.2.6
2.bowtie2-2.1.0
3.SEED
4.ALLPATHS-LG
5.SSPACE
Install the five tools and add them to the environment variable of Linux system.
Download dataset from NCBI BioProject (https://www.ncbi.nlm.nih.gov/sra?linkname=bioproject_sra_all&from_uid=344006) to this folder.

Pipeline to extract barcode pairs from the cloning vector pool of pBACode-2
##Pipeline input: flounder_BAC_barcodepair_bootstrap_only_<1/2>.fastq: raw read files.
                  configure_pool.txt: configeration file. the prefix of raw read files should be added to the line "prepool_reads". Other parameters have already been set.
##Pipeline output: barcodePrepoolflounder_BAC_barcodepair_bootstrap_only.txt, format: <left barcode>:<rightbarcode>[TAB]<read count>
#Merge mates of mated reads using FLASH
perl PEflash101.pl configure_pool.txt prepool
#Extract raw paired barcodes
perl barcodeRMvectorSite.pl configure_pool.txt
#Clean sequencing error and hybrid barcode pairs
perl barcodePrepool.pl configure_pool.txt

Pipeline to extract barcode pairs from BAC library using pBACode-2
##Pipeline input: flounder_BAC_barcodepair_part<1/2/3>_<1/2>.fastq: raw read files.
                  configure_pool.txt: configeration file. the prefix of raw read files should be added to the line "postpool_reads". Other parameters have already been set.
##Pipeline output: barcodePostpoolflounder_BAC_barcodepair_part<1/2/3>.txt, format: <left barcode>:<rightbarcode>[TAB]<read count> 
#Merge mates of mated reads using FLASH
perl PEflash101.pl configure_pool.txt postpool
#Extract raw paired barcodes
perl barcodeRMvectorSitePost.pl configure_pool.txt
#Clean sequencing error and hybrid barcode pairs
perl barcodePostpool.pl configure_pool.txt

ALLPATHS-LG assembly
#Configuration files for ALLPATHS-LG assembly: input_groups.csv and input_libraries.csv
ulimit -s 100000
mkdir -p flounder.genome/data
PrepareAllPathsInputs.pl DATA_DIR=$PWD/flounder.genome/data PLOIDY=1 IN_GROUPS_CSV=input_groups.csv IN_LIBS_CSV=input_libraries.csv OVERWRITE=True | tee prepare.out
ulimit -s 100000
RunAllPathsLG PRE=$PWD REFERENCE_NAME=flounder.genome DATA_SUBDIR=data RUN=run SUBDIR=test OVERWRITE=True MAXPAR=8 | tee -a assemble.out

Pipeline to process flounder BAC-PE data
##Pipeline input: flounder_BAC_PE_part<1/2/3>_<L/R>_<1/2>.fastq: raw read files.
                  barcodePostpoolflounder_BAC_barcodepair_part<1/2/3>.txt: barcode pair file
                  configure_BAC.txt: configeration file. The prefix of raw read files should be added to the line "left_reads" and "right_reads". The prefix of barcode pair file should be added to the line "pool_file". Other parameters have already been set.
##Pipeline output: barcodepairDistanceflounder_BAC_PE_part1_Lflounder_BAC_PE_part1_RbarcodePostpoolflounder_BAC_barcodepair_part1.txt, format: <left barcode>[space]<right barcode>[TAB]<chromosome>[TAB]<starting coordinate>[TAB]<ending coordinate>[TAB]<estimated BAC size>; 
                  barcodepairreadflounder_BAC_PE_part<1/2/3>_Lflounder_BAC_PE_part<1/2/3>_R.tab.txt, tab delimited file for SSPACE input, format: <left end scaffold ID>[TAB]<left end starting coordinate>[TAB]<left end ending coordinate>[TAB]<right end scaffold ID>[TAB]<right end starting coordinate>[TAB]<right end ending coordinate>
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
#Flip pair1 and pair2 reads, should run only once in the pipeline
perl rmVectorBarcodeRevcom.pl configure_BAC.txt
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
#Remove multiply-used barcodes
#perl barcodeRMNonunique.pl configure_BAC.txt left
#perl barcodeRMNonunique.pl configure_BAC.txt right
#perl barcodeRMNonuniquePair.pl configure_BAC.txt
#Pair BAC-PEs
perl flounderchrLen.pl configure_BAC.txt
perl barcodepairDistance.pl configure_BAC.txt
#To generate input for SSPACE assembler
perl barcoderead.pl configure_BAC.txt left
perl barcoderead.pl configure_BAC.txt right
perl barcodepairread_coordinate.pl configure_BAC.txt

Pipeline to locate BAC by barcode in 5D
##Pipeline input: filename.txt. All 36 pooled BAC fastq file names should be add to the file. Other parameters are already set.
##Pipeline output: Flatfish_BAC_library_3Dcross_cutoff_<threshold>.out, format: <Clone location>[TAB]<left barcode>
#To extract barcodes
perl system.pl filename.txt configure_pooling.txt
#To locate BAC by barcode
perl 5Dcross.pl <threshold, 0.01 or 0.05 recommended>

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
