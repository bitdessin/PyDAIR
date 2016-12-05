```
sra=("SRR017328" "SRR017329" "SRR017330" "SRR017331" "SRR017332" "SRR017333" "SRR017334" \
     "SRR017335" "SRR017336" "SRR017337" "SRR017338" "SRR017339" "SRR017340" "SRR017341")

for sid in ${sra[@]}
do
    prefetch ${sid}
    fastq-dump ${sid} -O ./data/
done
```
