##########################################################################
#   Identify NUMTs: Pairwise alignment of MT genome to Nuclear genome   #
##########################################################################

# Default settings.
#blastn -num_alignments 1000 -query $REF_MT -subject $REF_NUCL -outfmt '7' > blast.MTvsNUCL.${ASSEMBLY}.tab
#blastn -num_alignments 1000 -query $RELIN_MT -subject $REF_NUCL -outfmt '7' > blast.relinMTvsNUCL.${ASSEMBLY}.tab

# Following Tao et al. Genes 2023 to identify NUMTs in Human T2T.
blastn -evalue 0.0001 -gapopen 5 -gapextend 2 -penalty -3 -reward 2 -task blastn \
       -query $REF_MT -subject $REF_NUCL -outfmt '7' > MTvsNUCL.${ASSEMBLY}.tab

# Repeat the alignment with a relinearized MT genome (to account for breakpoint in assembly of a circular genome).
blastn -evalue 0.0001 -gapopen 5 -gapextend 2 -penalty -3 -reward 2 -task blastn \
       -query $RELIN_MT -subject $REF_NUCL -outfmt '7' > relinMTvsNUCL.${ASSEMBLY}.tab
