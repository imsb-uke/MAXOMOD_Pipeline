FOLDER=datasets/consortium/00_download/
xargs -a $FOLDER/SRR_Acc_List.txt -I {} -n 1 -P 8 fastq-dump --split-3 --disable-multithreading --skip-technical -O $FOLDER/fastq $FOLDER/prefetch/{}
