#Tang\_stats version 20200615

### Software for calculating iES and Rsb statistics. 
#### Following Tang, Thornton & Stoneking, PloS Biology 2007.


##### Usage:
	Tang_stats [genotype filename (one chrom)] [number of SNPs] [number of indiv] [threshold value (eg=0.1)] [seed (eg=123456)] [number pops] [size pop1] [size pop2] ... [size popN] [name_pop1] ... [name_popN]

##### Output file is automatically generated with the input filename plus '_Results_Tang.txt'

Note: The analysis uses unphased data. That is, the results obtained using the unphosed data are based on the observed homozygosity per individual. The results can be very different that using haplotypes (phased data), which uses the expected homozygosity assuming panmixia (you can uses the R library "rehh" to perform analysis with phased data).

#### To compile, enter into the source folder and type:
	gcc ./source/*.c -lm -o ./bin/Tang_stats -O3 -Wall

#### Examples:

	../bin/Tang_stats ./Window_1.txt 10000 61 0.1 123456 6 3 8 18 33 45 60 BRG BRM CHCU CHCA CHFR TXL
	Rscript --vanilla ../bin/run_plots_Rsb.R ./Window_1.txt_Results_Tang.txt 20 6 BRG BRM CHCU CHCA CHFR TXL
	
	../bin/Tang_stats ./Window_100_rows.txt 100 61 0.1 123456 6 5 10 15 12 15 4 BRG BRM CHCU CHCA CHFR TXL
	Rscript --vanilla ../bin/run_plots_Rsb.R ./Window_100_rows.txt_Results_Tang.txt 20 6 BRG BRM CHCU CHCA CHFR TXL
