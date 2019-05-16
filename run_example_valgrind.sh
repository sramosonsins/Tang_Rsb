#run example:
set -x

cd ./source
gcc *.c -lm -o ../bin/Tang_stats -O3 -Wall -g

cd ../example
valgrind --leak-check=full --track-origins=yes -v  ../bin/Tang_stats ./Window_1.txt 10000 61 0.1 123456 6 3 8 18 33 45 60 BRG BRM CHCU CHCA CHFR TXL 2>./Window_1.txt_valgrind.txt
#Rscript --vanilla ../bin/run_plots_Rsb.R ./Window_1.txt_Results_Tang.txt 20 6 BRG BRM CHCU CHCA CHFR TXL

#valgrind --leak-check=full --track-origins=yes -v  ../bin/Tang_stats ./Window_100_rows.txt 100 61 0.1 123456 6 3 8 18 33 45 60 BRG BRM CHCU CHCA CHFR TXL 2>./Window_100_rows.txt_valgrind.txt
#Rscript --vanilla ../bin/run_plots_Rsb.R ./Window_100_rows.txt_Results_Tang.txt 20 6 BRG BRM CHCU CHCA CHFR TXL

cd ..
