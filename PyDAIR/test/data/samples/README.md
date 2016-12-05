# Sample data

The scripts is used to create sample PYDAIR format files.

```
for (( i = 1; i < 4; ++i ))
do
pydair parse -q sample.${i}.fa \
             -v ../db/v.fa  -d ../db/d.fa -j ../db/j.fa \
             --v-blastdb ../db/v --d-blastdb ../db/d --j-blastdb ../db/j \
             --v-match-score 3 --v-mismatch-score -3         \
             --v-gap-open-penalty 6 --v-gap-extend-penalty 6 \
             --v-wordsize 21 --v-evalue-cutoff 1e-100         \
             --d-match-score 1 --d-mismatch-score -1         \
             --d-gap-open-penalty 0 --d-gap-extend-penalty 2 \
             --d-wordsize 4 --d-evalue-cutoff 1e-2           \
             --j-match-score 3 --j-mismatch-score -3         \
             --j-gap-open-penalty 6 --j-gap-extend-penalty 6 \
             --j-wordsize 7 --j-evalue-cutoff 1e-10           \
             -o sample.${i}
done
rm *.tsv
rm *.pydair.simple
rm *.unaligned.fa
rm *.vj.pydair

mv sample.1.vdj.pydair sample.1.pydair
mv sample.2.vdj.pydair sample.2.pydair
mv sample.3.vdj.pydair sample.3.pydair
```




