# Non-B DNA annotation and enrichment
code written by Kaivan Kamali, Linnéa Smeds and Edmundo Torres-Gonzalez

## Non-B DNA annotation
Non-B DNA annotation for A-phased repeats, short tandem repeats, direct repeats (slipped DNA), mirror repeats (triplex DNA), inverted repeats (cruciform DNA), and Z DNA was done with gfa ([https://github.com/abcsFrederick/non-B_gfa](https://github.com/abcsFrederick/non-B_gfa)); a detailed description of the file structure and how the code was run can be found in: [https://github.com/kxk302/non-B_gfa/blob/master/gfa_README.txt](https://github.com/kxk302/non-B_gfa/blob/master/gfa_README.txt).
G-quadruplexes (G4s) where annotated using Quadron. A dockerized version is found in [https://github.com/kxk302/Quadron_Docker](https://github.com/kxk302/Quadron_Docker).  

Quadron output format was converted to bed format using:
 ```
for i in {1..23} "X" "Y"
do
  awk -v chr=$i -v OFS="\t" '(/^DATA:/ && $5!="NA"){s=$2-1; e=s+$4; print "chr"chr,s,e,$3}' output/Quadron/chr${i}_out.txt |sort -k2,2n  >output/chr${i}_GQ.bed
done
 ```

## Identification of new sequence and non-B DNA enrichment
The code for this analyses can be found in the shell script [new_sequence_enrichment_commands.sh](https://github.com/makovalab-psu/T2T_primate_autosomes/blob/main/nonB_DNA/new_sequence_enrichment_commands.sh)), that requires that lastz and bedtools are installed.

#### Assembly versions used as 'old' (previous T2T) assemblies:
| Species | pre-T2T Version |
| -------- | ------- |
| Bonobo | panPan3.fa |
| Chimpanzee | panTro6.fa |
| Gorilla | gorGor6.fa |
| Sumatran orangutan | ponAbe3.fa |

*Note that there were no previous assemblies for Bornean orangutan nor siamang*

#### T2T assembly versions:
| Species | T2T Version |
| -------- | ------- |
| Bonobo | mPanPan1.pri.cur.20231122.fasta |
| Chimpanzee | mPanTro3.pri.cur.20231122.fasta |
| Gorilla | mGorGor1.pri.cur.20231122.fasta |
| Sumatran orangutan | mPonAbe1.pri.cur.20231205.fasta |
