FOLDER=datasets/consortium/00_download/
xargs -a $FOLDER/SRR_Acc_List.txt -n 1 -P 8 prefetch -O $FOLDER/prefetch {}
exit 0
