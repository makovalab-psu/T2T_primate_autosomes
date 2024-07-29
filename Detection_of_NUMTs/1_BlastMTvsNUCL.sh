##########################################################################
#   Identify NUMTs: Pairwise alignment of MT genome to Nuclear genome   #
##########################################################################

# Mitochondrial assemblies (REF_MT).
GORILLA=mGorGor1.MT.cur.20231122.fasta
BONOBO=mPanPan1.MT.cur.20231122.fasta
CHIMP=mPanTro3.MT.cur.20231122.fasta
SORANG=mPonAbe1.MT.cur.20231122.fasta
BORANG=mPonPyg2.MT.cur.20231122.fasta
SIAMANG=mSymSyn1.MT.cur.20231122.fasta
HUMAN=NC_012920.1.fa

# Nuclear T2T assemblies (REF_NUCL).
GORILLA=mGorGor1.pri.cur.20231122.fasta
BONOBO=mPanPan1.pri.cur.20231122.fasta
CHIMP=mPanTro3.pri.cur.20231122.fasta
SORANG=mPonAbe1.pri.cur.20231122.fasta
BORANG=mPonPyg2.pri.cur.20231122.fasta
SIAMANG=mSymSyn1.pri.cur.20231122.fasta
HUMAN=GCF_009914755.1_T2T-CHM13v2.0_genomic.fa

# Nuclear non-T2T assemblies (REF_NUCL).
GORILLA=mGorGor1.MT.cur.20231122.fasta
GORILLA_Y=GCA_015021865.1_gorGor.msY.makovalab.ver3_genomic.fna
BONOBO=panPan3.fa
BONOBO_Y=GCA_015021855.1_panPan.msY.makovalab.ver1_genomic.fna
CHIMP=panTro6.fa
SORANG=ponAbe3.fa
SORANG_Y=GCA_015021835.1_ponAbe.msY.makovalab.ver3_genomic.fna
HUMAN=GCF_000001405.40_GRCh38.p14_genomic.fna



# Following Tao et al. Genes 2023 to identify NUMTs in Human T2T.
blastn -evalue 0.0001 -gapopen 5 -gapextend 2 -penalty -3 -reward 2 -task blastn \
       -query $REF_MT -subject $REF_NUCL -outfmt '7' > MTvsNUCL.${ASSEMBLY}.tab

# Repeat the alignment with a relinearized MT genome (to account for breakpoint in a linear assembly of a circular genome).
blastn -evalue 0.0001 -gapopen 5 -gapextend 2 -penalty -3 -reward 2 -task blastn \
       -query $RELIN_MT -subject $REF_NUCL -outfmt '7' > relinMTvsNUCL.${ASSEMBLY}.tab
