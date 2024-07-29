################################################################################
# ALIGNMENT BETWEEN OLD AND NEW ASSEMBLIES
# Assuming a subdirectory for each species, with the old and new chromosomes
# named old.chr*.fa and t2t.chr*.fa. Note that the chromosomes are numbered
# according to the NEW assemblies, so for example gorilla/old.chr2.fa contains
# chr3 from gorGor6.


for sp in "chimp" "bonobo" "gorilla" "sorang"
do
  for i in {1..23}
  do
    lastz $sp/t2t.chr${i}.fa $sp/old.chr${i}.fa \
    --format=general:name1,size1,zstart1,end1,name2,strand2,zstart2+,end2+,nmatch,nmismatch,cgap,score,id%,blastid%,con% \
    --progress --ambiguous=iupac \
    > $sp/lastz.chr${i}.tab
  done
done

################################################################################
# MAKE BED FILES WITH ALIGNED SEQUENCE AND UNALIGNED SEQUENCE


# Make genome files for T2T assemblies
cut -f1,2 mPanTro3.pri.cur.20231122.fasta.fai |grep -v "random" >chimp/T2Tv2.genome
cut -f1,2 mPanPan1.pri.cur.20231122.fasta.fai |grep -v "random" >bonobo/T2Tv2.genome
cut -f1,2 mGorGor1.pri.cur.20231122.fasta.fai |grep -v "random" >gorilla/T2Tv2.genome
cut -f1,2 mPonAbe1.pri.cur.20231205.fasta.fai |grep -v "random" >sorang/T2Tv2.genome

# Extract aligned sequence (assume bed format), sort, merge and take complement
for sp in "chimp" "bonobo" "gorilla" "sorang"
do
  for i in {1..23}
  do
    cat $sp/lastz.chr${i}.*tab | cut -f 1,3,4 |sort -k2,2n -k3,3n |mergeBed -i - >$sp/chr${i}.aligned.bed
    grep "^chr${i}_" $sp/T2Tv2.genome | complementBed -i $sp/chr${i}.aligned.bed -g - |grep "^chr${i}_" >$sp/chr${i}.unaligned.bed
  done
done

# To compare lengths, make genome files of old sequences
cat ponAbe3.fa.fai |grep -v "random" | grep -v "Un" |sed 's/chr/hsa/' |sed 's/A/a/' |sed 's/B/b/' |cut -f1,2 >sorang/old.genome
cat panTro6.fa.fai |grep -v "random" | grep -v "Un" |sed 's/chr/hsa/' |sed 's/A/a/' |sed 's/B/b/' |cut -f1,2 >chimp/old.genome
cat gorGor6.fa.fai  |grep -v "random" | grep -v "Un" |sed 's/chr/hsa/' |sed 's/A/a/' |sed 's/B/b/' |cut -f1,2 >gorilla/old.genome
cat panPan3.fa.fai   |grep -v "random" | grep -v "Un" |sed 's/chr/hsa/' |sed 's/A/a/' |sed 's/B/b/' |cut -f1,2 >bonobo/old.genome

# Make sum of each unaligned segments, and compare to the size difference of
# new and old
for sp in "chimp" "bonobo" "gorilla" "sorang"
do
  for i in {1..23}
  do
  cat $sp/chr${i}.unaligned.bed |awk '{sum+=$3-$2; c=$1}END{print c, sum}'
  done |cut -f3 -d"_" |sort -k1,1 >$sp/sum.unaliged.txt
 cut -f3 -d"_" $sp/T2Tv2.genome | sort | join -1 1 -2 1 - <(sort $sp/old.genome) |awk '{d=$2-$3; print $0,d}' |sort -k1,1 | join -1 1 -2 1 -  <(sort -k1,1 $sp/sum.unaliged.txt) >$sp/comb.stats.txt
done


################################################################################
# GET FOLD DIFFERENCE FROM NON-B DNA DENSITY IN "NEW" VS "OLD" SEQUENCE

# Make merged files for the aligned (old) sequence
for sp in "chimp" "bonobo" "gorilla" "sorang"
do
  # Merged bedfile
  rm -f $sp/merged_aliged.bed
  for i in {1..23}
  do
    cat $sp/chr${i}.aligned.bed >>$sp/merged_aliged.bed
  done
done

# Intersect with non-B DNA annotation and calculate fold difference 
cat species_list.txt | while read -r sp latin;
do
  echo "============== $sp =============="
  for nb in "all" "GQ" "MR" "STR" "APR" "IR" "DR" "Z"
  do
    new=`cat non-B-DNA-Annotations/output/$latin/chr*_${nb}_merged.bed |intersectBed -a <(awk -v OFS="\t" '{split($1,s,"_"); print s[1],$2,$3}' $sp/merged_unaliged.bed) -b - -wao |cut -f1,2,3,7| awk -v OFS="\t" '{if(NR==0){chr=$1; s=$2; e=$3; sum=$4}else{if($1==chr && $2==s){sum+=$4}else{print chr,s,e,sum; chr=$1; s=$2; e=$3; sum=$4}}}END{print chr,s,e,sum}' | sed '/^\s*$/d' |awk -v nb=$nb '{sum_l+=$3-$2; sum_nb+=$4}END{d=sum_nb/sum_l; print nb,sum_l,sum_nb,d}'`
    old=`cat non-B-DNA-Annotations/output/$latin/chr*_${nb}_merged.bed |intersectBed -a <(awk -v OFS="\t" '{split($1,s,"_"); print s[1],$2,$3}' $sp/merged_aliged.bed) -b - -wao |cut -f1,2,3,7| awk -v OFS="\t" '{if(NR==0){chr=$1; s=$2; e=$3; sum=$4}else{if($1==chr && $2==s){sum+=$4}else{print chr,s,e,sum; chr=$1; s=$2; e=$3; sum=$4}}}END{print chr,s,e,sum}' | sed '/^\s*$/d' |awk '{sum_l+=$3-$2; sum_nb+=$4}END{d=sum_nb/sum_l; print sum_l,sum_nb,d}'`
    echo $new" "$old
 done
done
