# %% [markdown]
# # NUMTs: Analyze the BLAST alignments of MT and NUCL genomes for Primates

# %%
import pandas as pd
import re
import sys
import subprocess
import matplotlib.pyplot as plt
import seaborn as sns

# %%
#!ls ../../results/blast/*/blast*

# %%
#!tail -n2 ../../results/blast/*/blast*

# %% [markdown]
# ## Utils

# %%
def change_species_names(df):
    '''
    Shortens orangutan common names. 
    '''
    dict_names = {
        'Gorilla':'Gorilla', 'Bonobo':'Bonobo', 'Chimpanzee':'Chimpanzee', 'Human':'Human', 'Siamang':'Siamang',
        'Sumatran orangutan':'S. orangutan', 'Bornean orangutan':'B. orangutan'
    }
    df['Species'] = [ dict_names[row] for row in df['Species'] ]
    return df

def hsa_to_species_chr(hsa,genus_or_species):
    '''
    Takes the HSA and returns the species-specific chromosome.
    Uses Saswat's map for human analogs in primates (hsa).
    '''
    # Import Saswat's conversion table.
    df = pd.read_table("../../data/align.hsa.map.txt").fillna(0)
    df = df.rename(columns={'spe':'Species'})
    df.index = df['Species']
    # Process input species.
    if genus_or_species not in df['Species']:
        genus_or_species = {
            'Human':'hs1', 'Gorilla':'gor', 'Bonobo':'pan', 'Chimpanzee':'pan', 'Siamang':'sia', 
            'Sumatran orangutan':'pon', 'S. orangutan':'pon', 'Bornean orangutan':'pon', 'B. orangutan':'pon'
        }[genus_or_species]
    # Process input HSA.
    hsa = str(hsa)
    if hsa in ['hsaM','hsaUn_NW','hsaUn_NT']: #These rows do not have an HSA.
        return hsa.replace('hsa','chr')
    if hsa[:3] == 'chr':
        chr = hsa
        return chr
    # Special case for uppercase letters.
    if hsa[-2:] in ['2A','2B']:
        hsa = hsa.replace('2A','2a').replace('2B','2a')
    # Remove 'hsa' prefix.
    hsa = hsa.split('hsa')[-1]
    # Special case for humans.
    if genus_or_species == 'hs1':
        hsa = hsa.split('chr')[-1]
    # Special cases for gorilla hsa5 and hsa17.
    if hsa != 'x':
        hsa = ''.join(hsa.split('x'))
    if genus_or_species == 'gor':
        if hsa == '5':
            hsa = '517'
        elif hsa == '17':
            hsa = '175'
        else:
            pass
    # Find the corresponding CHR.
    chr = df[df['Species']==genus_or_species][hsa].iloc[0]
    if chr == 0:
        sys.exit(f'### This HSA <{hsa}> does not exist for the genus/species <{genus_or_species}> ###')
    chr = str(chr).split('.')[0]
    return f'chr{chr}'

def species_chr_to_hsa(chr,genus_or_species):
    '''
    Takes the HSA and returns the species-specific chromosome.
    Uses Saswat's map for human analogs in primates (hsa).
    '''
    # Import Saswat's conversion table.
    df = pd.read_table("../../data/align.hsa.map.txt").fillna(0)
    df = df.rename(columns={'spe':'Species'})
    df.index = df['Species']
    if chr[-2:] == '2A':
        sys.exit( f'### Expected a species-specific chromosome number, this is an HSA <{chr},{genus_or_species}> ###')
    # Process input species.
    if genus_or_species not in df['Species']:
        genus_or_species = {
            'Human':'hs1', 'Gorilla':'gor', 'Bonobo':'pan', 'Chimpanzee':'pan', 'Siamang':'sia', 
            'Sumatran orangutan':'pon', 'S. orangutan':'pon', 'Bornean orangutan':'pon', 'B. orangutan':'pon'
        }[genus_or_species]
    # Process input CHR.
    try:
        chr = int(chr)
    except:
        chr = chr.split('chr')[-1]
    try:
        chr = int(chr)
    except:
        pass
    # Find the corresponding HSA.
    df = df[df['Species']==genus_or_species]
    if df.columns[df.isin([chr]).any()].shape[0] > 0:
        hsa = df.columns[df.isin([chr]).any()][0]
    else:
        sys.exit(f'### This chr <{chr}> does not exist for the genus/species <{genus_or_species}> ###')
    if hsa == '517':
        hsa = '5x17'
    elif hsa == '175':
        hsa = '17x5'
    return f'hsa{hsa}'

def nc_to_human_chr(df_in):
    '''
    Change NC IDs to human-readable human chromosomes.
    '''
    df = df_in.copy()
    # Copy Subject to CHR column.
    df.loc[df['Species']=='Human','CHR'] = df['Subject']
    # Rename 'Unknown' chromosomes.
    df.loc[df['Subject'].str.contains('NT'),'CHR'] = 'chrUn_NT'
    df.loc[df['Subject'].str.contains('NW'),'CHR'] = 'chrUn_NW'
    # Change NC IDs to chromosome numbers.
    df['HSA'] = df['CHR'].map(get_dict_human_chr()).fillna(df['HSA'])
    df.loc[df['Species']=='Human','CHR'] = df['HSA']
    return df

def get_CHR_from_Subject(df):
    '''
    The subject column can be split into Chromosome, haplotype, and HSA for T2T primate samples.
    Other samples use different methods.
    '''
    # This properly annotates CHR and HSA for only non-human T2T assemblies.
    df['CHR'] = df['Subject'].str.split('_').str[0]
    df['HSA'] = df['Subject'].str.split('_').str[2].fillna(df['CHR'])
    # Annotate chromosomes of Human assemblies.
    df = nc_to_human_chr(df)
    # Annotate non-human pre-T2T assemblies. 
    ### These used "hsa" numbering as the standard numbering, 
    ### instead of species-specific numbering by chromosome size.
    assemblies_nonHuman_preT2T = filter( lambda s: not (s in ['CHM13','GRCh38'] or s[0]=="m" or s[-2:]=='_y'), 
                                        get_dict_species().keys()) #list of (non-human, pre-T2T) assemblies.
    for ASSEMBLY in list(assemblies_nonHuman_preT2T):
        # Subject_ID had the HSA chromosome number.
        df.loc[df['Assembly']==ASSEMBLY,'HSA'] = df['Subject'].str.split('_').str[0].str.replace('chr','hsa')
        # Highlight unfinished hsa numbers in CHR column.
        df.loc[df['Assembly']==ASSEMBLY,'CHR'] = df['CHR'].str.replace('chr','hsa')
    assemblies_chrY = filter( lambda s: s[-2:]=='_y', 
                                        get_dict_species().keys()) #list of Y chr assemblies.
    # Annotates Y chr.
    df.loc[df['Assembly'].str.contains("_Y"),'HSA'] = "hsaY"
    df.loc[df['Assembly'].str.contains("_Y"),'CHR'] = "chrY"
    # Annotate with species-specific chromosome.
    df["CHR"] = df.apply(lambda x: hsa_to_species_chr(x["CHR"], x["Species"]) , axis=1)
    return df

def exclude_unknown_chr(df_in):
    '''
    Excludes secondary contigs from non-T2T genomes.
    These are called 'chrUn' and '_random' chromosomes, or have 'NT' and 'NW' chromosome IDs.
    '''
    df = df_in.copy()
    # Exclude 'chrUn' and '_random' chromosomes from non-T2T genomes.
    df = df[~df['Subject'].str.contains("chrUn")]
    df = df[~df['Subject'].str.contains("_random")]
    # Exclude other Unknown and incomplete contigs from GRCh38 chromosomes.
    df = df[~df['Subject'].str.contains("NT")]
    df = df[~df['Subject'].str.contains("NW")]
    # Exclude chrM from BLAST results (mostly non-T2T genomes).
    df = df[(~df['Subject'].str.contains("chrM"))]
    df = df[(~df['Subject'].str.contains("NC_012920.1"))]
    return df

def exclude_sex_chr(df_in):
    '''
    Excludes X and Y chromosomes from all genomes.
    '''
    df = df_in.copy()
    # Exclude 'chrX' and 'chrY' chromosomes from all genomes.
    if 'Subject' in df.columns:
        df = df[~df['Subject'].str.contains("chrX")]
        df = df[~df['Subject'].str.contains("chrY")]
    df = df[~df['CHR'].str.contains("chrX")]
    df = df[~df['CHR'].str.contains("chrY")]
    return df

def get_category_assembly(ass):
    '''
    Takes assembly name as input and returns whether it is a T2T assambly.
    '''
    if ass[0] == 'm':
        categ = 'T2T'
    elif ass == 'CHM13':
        categ = 'T2T'
    elif ass[-2:] == '_Y':
        categ = 'Y-chr'
    else:
        categ = 'non-t2t'
    return categ

def get_dict_species():
    '''
    Links assembly names to scientific and common species names.
    '''
    dict_species = {
                'mPanPan1':['Pan paniscus','Bonobo'],
                'panPan3':['Pan paniscus','Bonobo'],
                'mPanTro3':['Pan troglodytes','Chimpanzee'],
                'panTro6':['Pan troglodytes','Chimpanzee'],
                'CHM13':['Homo sapiens','Human'],
                'GRCh38':['Homo sapiens','Human'],
                'mGorGor1':['Gorilla gorilla','Gorilla'],
                'gorGor6':['Gorilla gorilla','Gorilla'],
                'mPonAbe1':['Pongo abelii','Sumatran orangutan'],
                'mPonAbe1.hap1':['Pongo abelii','Sumatran orangutan'],
                'mPonAbe1.hap2':['Pongo abelii','Sumatran orangutan'],
                'ponAbe3':['Pongo abelii','Sumatran orangutan'],
                'mPonPyg2':['Pongo pygmaeus','Bornean orangutan'],
                'mPonPyg2.hap1':['Pongo pygmaeus','Bornean orangutan'],
                'mPonPyg2.hap2':['Pongo pygmaeus','Bornean orangutan'],
                'mSymSyn1':['Symphalangus syndactylus (gibbon)','Siamang'],
                'gorilla_Y':['Gorilla gorilla','Gorilla'],
                'bonobo_Y':['Pan paniscus','Bonobo'],
                'sorang_Y':['Pongo abelii','Sumatran orangutan'],
            }
    return dict_species

def get_dict_human_chr():
    '''
    Changes human chromosome names to common usage.
    Affects human CHM13 and GRCh38 assemblies.
    '''
    dict_chr = {
        'NC_060925.1':'chr1','NC_060926.1':'chr2','NC_060927.1':'chr3',
        'NC_060928.1':'chr4','NC_060929.1':'chr5','NC_060930.1':'chr6',
        'NC_060931.1':'chr7','NC_060932.1':'chr8','NC_060933.1':'chr9',
        'NC_060934.1':'chr10','NC_060935.1':'chr11','NC_060936.1':'chr12',
        'NC_060937.1':'chr13','NC_060938.1':'chr14','NC_060939.1':'chr15',
        'NC_060940.1':'chr16','NC_060941.1':'chr17','NC_060942.1':'chr18',
        'NC_060943.1':'chr19','NC_060944.1':'chr20','NC_060945.1':'chr21',
        'NC_060946.1':'chr22','NC_060947.1':'chrX','NC_060948.1':'chrY',
        'NC_000001.11':'chr1','NC_000002.12':'chr2','NC_000003.12':'chr3',
        'NC_000004.12':'chr4','NC_000005.10':'chr5','NC_000006.12':'chr6',
        'NC_000007.14':'chr7','NC_000008.11':'chr8','NC_000009.12':'chr9',
        'NC_000010.11':'chr10','NC_000011.10':'chr11','NC_000012.12':'chr12',
        'NC_000013.11':'chr13','NC_000014.9':'chr14','NC_000015.10':'chr15',
        'NC_000016.10':'chr16','NC_000017.11':'chr17','NC_000018.10':'chr18',
        'NC_000019.10':'chr19','NC_000020.11':'chr20','NC_000021.9':'chr21',
        'NC_000022.11':'chr22','NC_000023.11':'chrX','NC_000024.10':'chrY',
        'NC_012920.1':'chrM', 'chrUn_NT':'chrUn_NT', 'chrUn_NW':'chrUn_NW'
    }
    return dict_chr

def get_key(my_dict,val_in):
    '''
    Returns dictionary key when given the value.
    '''
    for key, value in my_dict.items():
        if val_in == value:
            return key
        else:
            sys.exit(f"### The key for value <{val_in}> doesn't exist ###")

#%tb
#hsa_to_species_chr('hsaX','gor')
#species_chr_to_hsa('chr4','gor')

# %% [markdown]
# ## Import BLAST results for original reference

# %%
def get_blast_results():
    '''
    Imports BLASTn results as pandas dataframes.
    '''
    dict_species = get_dict_species()
    
    cols = [ 
        'Query', 'Subject', 'Perc_ident', 'Align_length', 'Mismatches', 'Gap_opens', 'Qu_START', 'Qu_END', 'Su_START', 'Su_END', 'E-value', 'Bit-score',
        'Species', 'Assembly', 'Scientific_name'
    ]
    
    list_dfs = []
    for ASSEMBLY in dict_species.keys(): #iterate over all species.
        #print(ASSEMBLY)
        try:
            #df_iter = pd.read_table( f'../../results/blast/{ASSEMBLY}/blast.MTvsNUCL.{ASSEMBLY}.tab', comment='#', header=None, names = cols )
            #df_iter = pd.read_table( f'../../results/blast_new/blast.MTvsNUCL.{ASSEMBLY}.tab', comment='#', header=None, names = cols )
            df_iter = pd.read_table( f'../../results/blast_new_01/blast01.MTvsNUCL.{ASSEMBLY}.tab', comment='#', header=None, names = cols )
        except:
            sys.exit(f"### Did not find a file for assembly {ASSEMBLY} ###")
        '''
        # Fix human chromosome names.
        if ASSEMBLY in ['CHM13','GRCh38']:
            df_iter.loc[df_iter['Subject'].str.contains('NT'),'Subject'] = 'chrUn_NT'
            df_iter.loc[df_iter['Subject'].str.contains('NW'),'Subject'] = 'chrUn_NW'
            df_iter['Subject'] = [ get_dict_human_chr(ASSEMBLY)[NC_ID] for NC_ID in df_iter['Subject'] ]
        '''
        # Annotate species names (common, assembly, scientific).
        common_name = dict_species[ASSEMBLY][1]
        df_iter['Species'] = common_name
        df_iter['Assembly'] = ASSEMBLY
        df_iter['Scientific_name'] = dict_species[ASSEMBLY][0]
        
        #print( common_name, df_iter.shape[0] )

        # Add to list of dataframes.
        list_dfs.append(df_iter)

    # Concatenate all dataframes.
    df = pd.concat(list_dfs).reset_index(drop=True)
    
    # Exclude chrM from BLAST results (mostly non-T2T genomes).
    df = df[(~df['Subject'].str.contains("chrM"))]
    df = df[(~df['Subject'].str.contains("NC_012920.1"))]
    
    # Get the Subject chromosome info.
    #df['CHR'] = df['Subject'].str.split('_').str[0]#.str.split('r').str[-1]
    df = get_CHR_from_Subject(df)
    
    # Final merged dataset size.
    #print( "\n Total: ", df.shape )

    # Filter by E-value (≤0.005).
    df = df[ df['E-value'] <= 0.005 ]
    #print( "\n Total (E-value ≤0.005): ", df.shape )

    # Column data types.
    df['Su_START'] = df['Su_START'].astype('int')
    df['Su_END'] = df['Su_END'].astype('int')
    #print('\n',df.dtypes)

    return df

#%tb
#get_blast_results()

# %% [markdown]
# ## Import BLAST results for shifted breakpoint

# %%
def get_shifted_blast_results():
    '''
    Imports BLASTn results for the alternate (shifted breakpoint) mtDNA
    '''
    dict_species = get_dict_species()
    
    cols = [ 
        'Query', 'Subject', 'Perc_ident', 'Align_length', 'Mismatches', 'Gap_opens', 'Qu_START', 'Qu_END', 'Su_START', 'Su_END', 'E-value', 'Bit-score',
        'Species', 'Assembly', 'Scientific_name'
    ]
    
    list_dfs = []
    for ASSEMBLY in dict_species.keys():
        try:
            #df_iter = pd.read_table( '../../results/blast/'+ASSEMBLY+'/blast.relinMTvsNUCL.'+ASSEMBLY+'.tab', comment='#', header=None, names = cols )
            df_iter = pd.read_table( f'../../results/blast_new_01/blast01.relinMTvsNUCL.{ASSEMBLY}.tab', comment='#', header=None, names = cols )
        except:
            sys.exit(f"### Did not find a file for assembly {ASSEMBLY} ###")

        # Annotate species names (common, assembly,scientific).
        common_name = dict_species[ASSEMBLY][1]
        df_iter['Species'] = common_name
        df_iter['Assembly'] = ASSEMBLY
        df_iter['Scientific_name'] = dict_species[ASSEMBLY][0]

        # Update shifted breakpoint coordinates to match original reference coords.
        df_iter = fix_shifted_coords(common_name,df_iter).fillna('not_annot')

        #print( common_name, df_iter.shape )

        # Add to list of dataframes.
        list_dfs.append(df_iter)

    # Concatenate all dataframes.
    df = pd.concat(list_dfs).reset_index(drop=True)

    # Exclude chrM from BLAST results (mostly non-T2T genomes).
    df = df[(~df['Subject'].str.contains("chrM"))]
    df = df[(~df['Subject'].str.contains("NC_012920.1"))]
    
    # Get the Subject chromosome and HSA info.
    df = get_CHR_from_Subject(df)

    # Final merged dataset size.
    #print( "\n Total: ", df.shape )

    # Filter by E-value (≤0.005).
    df = df[ df['E-value'] <= 0.005 ]
    #print( "\n Total (E-value ≤0.005): ", df.shape )

    # Column data types.
    df['Su_START'] = df['Su_START'].astype('int')
    df['Su_END'] = df['Su_END'].astype('int')
    #print('\n',df.dtypes)

    # Exclude chrM from non-T2T genomes.
    df = df[~df['Subject'].str.contains("chrM")]
    
    return df


## Realigns the relinearized positions to match the original reference.
def fix_shifted_coords(species,df_relin):
    # Select genome size based on species.
    dict_genome_sizes = {
        'Gorilla':16407, 'Bonobo':16569, 'Chimpanzee':16560, 
        'Sumatran orangutan':16496, 'Bornean orangutan':16461, 'Siamang':16515, 'Human':16564
    }
    genome_size = dict_genome_sizes[species]
    # Fixing the position to match the original reference.
    df_relin['Qu_START'] += (8400 - genome_size - 1)
    df_relin.loc[df_relin['Qu_START'] <= 0, 'Qu_START'] += genome_size
    df_relin['Qu_END'] += (8400 - genome_size - 1)
    df_relin.loc[df_relin['Qu_END'] <= 0, 'Qu_END'] += genome_size
    ## Sort by new position.
    ##df_relin = df_relin.sort_values(['Subject','Qu_START'])
    return df_relin


#get_shifted_blast_results()#.head()

# %% [markdown]
# # Join alignment results for BLAST for original reference and the one with a shifted breakpoint

# %%
def join_blast_results(write=True):
    '''
    Joins BLASTn results from the two mtDNA alignments (excluding duplicates).
    '''
    df1 = get_blast_results()
    df2 = get_shifted_blast_results()
    # Drop repeated events.
    df_merged = pd.concat([df1,df2]).drop_duplicates(subset=['Assembly','Subject','Qu_START','Qu_END','Su_START','Su_END']).reset_index(drop=True)
    # Exclude chrM from BLAST results (mostly non-T2T genomes).
    df_merged = df_merged[(~df_merged['Subject'].str.contains("chrM"))]
    df_merged = df_merged[(~df_merged['Subject'].str.contains("NC_012920.1"))]
    # New columns to clarify relevant (nuclear coordinates).
    df_merged['START'] = df_merged['Su_START']
    df_merged['END'] = df_merged['Su_END']
    if write == True:
            df_merged.to_csv( f"../../results/blast_new_01/blast01.MTvsNUCL.correctedByRelin.ALL.tab", sep='\t', header=None, index=None )
    return df_merged

# Blast results corrected with shifted breakpoint.
#join_blast_results().head()

# %% [markdown]
# ## Export coordinates as BED files.

# %%
def export_BED_numt_coords(df_in=False, write=True):
    '''
    Export a BED file to use bedtools intersect.
    '''
    if isinstance(df_in, pd.DataFrame):
        df_blast = df_in.copy()
    else:
        df_blast = join_blast_results()
        
    list_assemblies = get_dict_species().keys()
    for ASSEMBLY in list_assemblies:
        df = df_blast[df_blast['Assembly']== ASSEMBLY]
        #df = df[['CHR','Su_START','Su_END']]
        df = df[['Subject','Su_START','Su_END','Perc_ident']]
        if df.shape[0] == 0: #blast for gorilla_Y did not return results.
            if write == True:
                subprocess.run( f"echo '' > ../../results/blast_new_01/blast.rawCoords.{ASSEMBLY}.bed")
                subprocess.run( f"echo '' > ../../results/blast_new_01/merged.clustered.{ASSEMBLY}.bed")
                continue
            continue
        # Determine stand direction.
        df.loc[ df['Su_END'] - df['Su_START'] > 0 , 'Su_Strand' ] = '+'
        df.loc[ df['Su_END'] - df['Su_START'] < 0 , 'Su_Strand' ] = '-'
        # Reorganize positions based on strand direction.
        df.loc[ df['Su_Strand'] == '+' , 'START' ] = df['Su_START']
        df.loc[ df['Su_Strand'] == '+' , 'END' ] = df['Su_END']
        df.loc[ df['Su_Strand'] == '-' , 'START' ] = df['Su_END']
        df.loc[ df['Su_Strand'] == '-' , 'END' ] = df['Su_START']
        # Change data types.
        df['START'] = df['START'].astype('int')
        df['END'] = df['END'].astype('int')
        # Add info to 'Name' column (useful for IGV tracks).
        lengths = abs(df['END'] - df['START']).astype('str')
        perc_ident = round(df['Perc_ident'],0).astype(str).str.split('.').str[0]
        df['Annotation'] = "NUMT_"+ lengths +"bp_"+ perc_ident +"%"
        # Empty columns.
        df['Empty'] = '.'
        # Presort your data by chromosome and then by start position.
        df = df[['Subject','START','END','Annotation','Empty','Su_Strand']].sort_values(['Subject','START','END'])
        #### Strand (+/-) must be 6th row in BED file ####
        # Export as BED file.
        if write == True:
            #df_out.to_csv( f"../../results/blast/{ASSEMBLY}/blast.rawCoords.{ASSEMBLY}.bed", sep='\t', header=None, index=None )
            df.to_csv( f"../../results/blast_new_01/blast.rawCoords.{ASSEMBLY}.bed", sep='\t', header=None, index=None )

    return df.tail()
    

#export_BED_numt_coords()

# %% [markdown]
# # Resolve overlapping coordinates.

# %% [markdown]
# ## Run bash script to merge/cluster NUMTs by overlap and proximity using bedtools.

# %%
'''
#!/bin/bash
set -ue

for ASSEMBLY in mGorGor1 mPanPan1 mPanTro3 mPonAbe1 mPonPyg2 mSymSyn1 gorGor6 panTro6 ponAbe3 panPan3 CHM13 GRCh38 gorilla_Y bonobo_Y sorang_Y
do

# Skips assembly if it did not have Blast coordinates.
if [[ $ASSEMBLY == 'gorilla_Y' ]]
then
 echo "Skipping: ${ASSEMBLY}"
 cat rawCoords.${ASSEMBLY}.bed > mergedCoords.${ASSEMBLY}.bed
 continue
fi

# Merge overlapping NUMT coordinates (regardless of strand info) and within 1000 bp flanks.
bedtools merge -d 1000 -c 4,5,6 -o distinct -i rawCoords.${ASSEMBLY}.bed \
 | awk -v OFS='\t' '{ TOTALBP = $3 - $2 } $4 ~ /,/ { $4 = "NUMT_" TOTALBP "bp_merged" }  1' \
 > mergedCoords.${ASSEMBLY}.bed

done

echo "DONE"
'''

# %%
#subprocess.run("bash mergeCoords.sh blast_new", shell=True, check=True)

# %% [markdown]
# ## Import raw coordinates (derived from Blast results).

# %%
def get_raw_coords():
    dict_species = get_dict_species()

    cols = ['Subject','START','END','Annotation','Empty','Su_Strand','Species']
    
    df = pd.DataFrame(columns=cols)
    
    for ASSEMBLY in dict_species.keys():
        #df_iter = pd.read_table( "../../results/blast/"+ASSEMBLY+"/blast.rawCoords."+ASSEMBLY+".bed", header=None, names = cols )
        df_iter = pd.read_table( f"../../results/blast_new_01/blast.rawCoords.{ASSEMBLY}.bed", header=None, names = cols )
        
        # Annotate species names (common, assembly,scientific).
        common_name = dict_species[ASSEMBLY][1]
        df_iter['Species'] = common_name
        df_iter['Assembly'] = ASSEMBLY

        print( common_name, df_iter.shape )
        
        # Merge results for this iteration.
        df = pd.concat( [df_iter, df] ).reset_index(drop=True)
    
    # Final merged dataset size.
    print( "\n Total: ", df.shape )

    return df


#get_raw_coords()#.tail()


# %% [markdown]
# ## Import merged coordinates (accounts for overlaps rergardless of strand direction).
# 
# Used `bedtools merge` on `numt.rawCoords.${ASSEMBLY}.bed` for each species. This joined overlapping coordinates but did not use strand information.

# %%
def get_merged_coords():
    '''
    Import BED file after merging overlapping coordinates.
    All coordinates should be on the same strand (+).
    '''
    dict_species = get_dict_species()

    cols = ['Subject','START','END','Annotation','Empty','Su_Strand']
    
    df = pd.DataFrame(columns=cols)
    
    for ASSEMBLY in dict_species.keys():
        try:
            #df_iter = pd.read_table( "../../results/blast/"+ASSEMBLY+"/noStrand.blast.mergedCoords."+ASSEMBLY+".bed", header=None, names = cols )
            #df_iter = pd.read_table( f"../../results/blast_new/noStrand.blast.mergedCoords.{ASSEMBLY}.bed", header=None, names = cols )
            df_iter = pd.read_table( f"../../results/blast_new_01/merged.clustered.{ASSEMBLY}.bed", header=None, names = cols )
        except:
            sys.exit(f"### Did not find a file for assembly {ASSEMBLY} ###")

        # Annotate species names (common, assembly,scientific).
        common_name = dict_species[ASSEMBLY][1]
        df_iter['Species'] = common_name
        df_iter['Assembly'] = ASSEMBLY

        #print( common_name, df_iter.shape )

        # Merge results for this iteration.
        df = pd.concat( [df_iter, df] ).reset_index(drop=True)
    
    # Final merged dataset size.
    #print( "\n Total: ", df.shape )

    return df


#get_merged_coords().head()


# %% [markdown]
# # What is the Percentage of the genome derived from NUMTs?

# %% [markdown]
# ## Import chromosome sizes for each species.

# %%
def get_chrom_sizes():
    '''
    Chromosome or contig sizes.
    '''
    dict_species = get_dict_species()

    cols = ['Subject','Chrom_size']
    
    df = pd.DataFrame(columns=cols)
    
    for ASSEMBLY in dict_species.keys():
        df_iter = pd.read_table( f"../../results/chromSizes/chromSizes.{ASSEMBLY}.tab", header=None, names = cols )
        #if ASSEMBLY == 'CHM13':
            #df_iter = pd.read_table( "../../results/blast/"+ASSEMBLY+"/chromSizes."+ASSEMBLY+".tab", header=None, names = cols )
        #elif ASSEMBLY != 'CHM13':
            #df_iter = pd.read_table( "../../results/blast/"+ASSEMBLY+"/chromSizes."+ASSEMBLY[:-1]+".tab", header=None, names = cols )
        
        # Annotate species names (common, assembly,scientific).
        common_name = dict_species[ASSEMBLY][1]
        df_iter['Species'] = common_name
        df_iter['Assembly'] = ASSEMBLY

        #print( common_name, df_iter.shape )

        # Merge results for this iteration.
        df = pd.concat( [df_iter, df] ).reset_index(drop=True)
    
    # Final merged dataset size.
    #print( "\n Total: ", df.shape )

    df = df.reset_index(drop=True)
    # Add CHR column.
    df = get_CHR_from_Subject(df)

    # Sum Y contigs.
    df = df.drop('Subject',axis=1)
    df = df.groupby(['Species','Assembly','CHR','HSA'])['Chrom_size'].sum().reset_index()

    return df


# Sum chromosome sizes.
def calc_genome_size(keep=False):
    '''
    Genome size for each genome assembly. Sum of all contigs.
    '''
    df = get_chrom_sizes()

    if keep == 'autosome':
        df = exclude_sex_chr(df)
    elif keep == 'x':
        df = df[df['CHR']=='chrX']
    elif keep == 'y':
        df = df[df['CHR']=='chrY']
        #df = df.value_counts()
        #df = df.drop('Subject', axis=1)
        #return df
    else:
        pass
    
    list_assemblies = df['Assembly'].drop_duplicates().tolist()
    genome_sizes = [ df[df['Assembly']==assembly]['Chrom_size'].sum() for assembly in list_assemblies ]

    # Turn two lists into a dictionary.
    tuples = [(key, value)
          for i, (key, value) in enumerate(zip(list_assemblies, genome_sizes))]

    dict_genome_sizes = dict(tuples)

    return dict_genome_sizes


#get_chrom_sizes()#.head()


# %% [markdown]
# ## Abundance, Total bp and Percentage of Genome covered by NUMTs

# %%
def compute_total_numt_bp( df_in=False, excl_second_contigs=True, write=True, exclude_sex=True):
    if isinstance(df_in, pd.DataFrame):
        df = df_in.copy()
    else:
        df = get_merged_coords()
    if 'CHR' not in df:
        df = get_CHR_from_Subject(df)
    df = df[['Species','Assembly','Subject','CHR','START','END']].drop_duplicates()
    # Removes secondary contigs from non-T2T assemblies.
    if excl_second_contigs == True:
        df = exclude_unknown_chr(df)
    # Removes Sex chromosomes.
    if exclude_sex == True:
        df = exclude_sex_chr(df)
    # Categorize assemblies.
    df['Version'] = 'non-T2T'
    df.loc[df['Assembly'].str[0]=='m','Version'] = 'T2T'
    df.loc[df['Assembly']=='CHM13','Version'] = 'T2T'
    df = df.sort_values(by=['Species','Version'])
    # Length of each NUMT.
    df['Length'] = abs(df['END'] - df['START'] + 1)
    # Abundance of NUMTs by species.
    df_abundance = df.value_counts(['Species','Version','Assembly']).reset_index().rename(columns={0:'Abundance','count':'Abundance'})
    # Sum of NUMT lengths.
    df_totalbp = df.groupby(['Species','Version','Assembly'])['Length'].sum().reset_index().rename(columns={'Length':'Total_bp'})
    df = pd.merge( df_abundance, df_totalbp ) #merge results.
    # Numts as percentage of genome.
    dict_genome_sizes = calc_genome_size(keep='autosomes') #genome sizes per species assembly.
    df['Genome_size'] = [ dict_genome_sizes[assembly] for assembly in df['Assembly'].drop_duplicates().tolist() ]
    df = df.assign(Perc_genome=100*df['Total_bp']/df['Genome_size']) #percentage of genome covered by NUMTs.
    # Export as table.
    if write == True:
        df.to_csv("tables/abundance_Numts_autosomes.table", index=None, sep='\t')
    return df


#compute_total_numt_bp()

# %% [markdown]
# ### Abundance by Chromosome

# %%
def compute_total_numt_bp_by_chr(df_in=False, excl_second_contigs=True, write=True):
    if isinstance(df_in, pd.DataFrame):
        df = df_in.copy()
    else:
        df = get_merged_coords()
    if 'CHR' not in df:
        df = get_CHR_from_Subject(df)
    df = df[['Species','Assembly','Subject','CHR','HSA','START','END']].drop_duplicates()
    # Exclude Y-chr.
    df = df[df['CHR']!='chrY']
    # Removes secondary contigs from non-T2T assemblies.
    if excl_second_contigs == True:
        df = exclude_unknown_chr(df)
    # Categorize assemblies.
    df['Version'] = 'non-T2T'
    df.loc[df['Assembly'].str[0]=='m','Version'] = 'T2T'
    df.loc[df['Assembly']=='CHM13','Version'] = 'T2T'
    # Length of each NUMT.
    df['Length'] = abs(df['END'] - df['START'] + 1)
    # Abundance of NUMTs by species.
    df_abundance = df.value_counts(['Species','Version','Assembly','CHR','HSA']).reset_index().rename(columns={0:'Abundance','count':'Abundance'})
    # Sum of NUMT lengths.
    df_totalbp = df.drop(['Subject'], axis=1).groupby(['Species','Version','Assembly','CHR','HSA'])['Length'].sum().reset_index().rename(columns={'Length':'Total_bp'})
    df_out = pd.merge( df_abundance, df_totalbp )
    # Numts as percentage of genome.
    df_chrom_sizes = get_chrom_sizes() #chrom sizes per species assembly.
    df_out = pd.merge( df_out, df_chrom_sizes )
    df_out = df_out.assign(Perc_chrom=100*df_out['Total_bp']/df_out['Chrom_size']) #percentage of genome covered by NUMTs.
    df = df_out.sort_values(['CHR','Species'])
    # Filter for X chr.
    df_x = df[df['CHR']=='chrX'].sort_values(by=['Species','Version'])
    df_x =  df_x.drop(['HSA'], axis=1)
    # Export as table.
    if write == True:
        df.to_csv("tables/abundance_Numts_byChromosome_sansY.table", index=None, sep='\t')
        df_x.to_csv("tables/abundance_Numts_Xchr.table", index=None, sep='\t')
    return df_x
    

#compute_total_numt_bp_by_chr()

# %% [markdown]
# ### Abundance in Y chr

# %%
def compute_total_numt_bp_chrY(df_in=False, excl_second_contigs=True, write=True):
    if isinstance(df_in, pd.DataFrame):
        df = df_in.copy()
    else:
        df = get_merged_coords()
    if 'CHR' not in df:
        df = get_CHR_from_Subject(df)
    df = df[['Species','Assembly','Subject','CHR','HSA','START','END']].drop_duplicates()
    # Keep just Y-chr.
    df = df[df['CHR']=='chrY']
    # Removes secondary contigs from non-T2T assemblies.
    if excl_second_contigs == True:
        df = exclude_unknown_chr(df)
    # Categorize assemblies.
    df['Version'] = 'non-T2T'
    df.loc[df['Assembly'].str[0]=='m','Version'] = 'T2T'
    df.loc[df['Assembly']=='CHM13','Version'] = 'T2T'
    # Length of each NUMT.
    df['Length'] = abs(df['END'] - df['START'] + 1)
    df = df[df['CHR']=='chrY'].drop('Subject', axis=1)
    # Abundance of NUMTs by species.
    df_abundance = df.value_counts(['Species','Version','Assembly','CHR']).reset_index().rename(columns={0:'Abundance','count':'Abundance'})
    # Collapse each non-T2T Y-chr assembly into one row.
    #df = df.drop_duplicates().groupby(['Species','CHR','Version','Assembly']).sum().reset_index()
    #df.loc[(df['Species']=='Bonobo')&(df['Version']=='non-T2T'), 'Assembly'] = 'bonobo_Y' #GCA_015021855.1
    #df.loc[(df['Species']=='Sumatran orangutan')&(df['Version']=='non-T2T'), 'Assembly'] = 'sorang_Y' #GCA_015021835.1
    #df.loc[(df['Species']=='Gorilla')&(df['Version']=='non-T2T'), 'Assembly'] = 'gorilla_Y' #GCA_015021865.1
    # Sum of NUMT lengths.
    df_totalbp = df.groupby(['Species','Version','Assembly','CHR'])['Length'].sum().reset_index().rename(columns={'Length':'Total_bp'})
    df = pd.merge( df_abundance, df_totalbp )
    # Numts as percentage of genome.
    df_chrom_sizes = get_chrom_sizes() #chrom sizes per species assembly.
    df = pd.merge( df, df_chrom_sizes )
    #df = df.drop('HSA',axis=1)
    df = df.assign(Perc_chrom=100*df['Total_bp']/df['Chrom_size']) #percentage of genome covered by NUMTs.
    df = df.sort_values(['CHR','Species'])
    if write == True:
        df.to_csv("tables/abundance_Numts_Ychr.table", index=None, sep='\t')
    return df

#compute_total_numt_bp_chrY()


# %% [markdown]
# ### Abundance in Autosomes + XY

# %%
def join_total_numt_bp_allchr(write=True):
    '''
    Joins the NUMT results for Autosomes + X, and Y.
    ''' 
    # Includes Autosomes and X.
    df_aut = pd.read_table( "tables/abundance_Numts_byChromosome_sansY.table" )
    # Includes Y.
    df_y = pd.read_table( "tables/abundance_Numts_Ychr.table" )
    df = pd.concat( [df_aut,df_y] )
    if write == True:
        df.to_csv("tables/abundance_Numts_byChromosome_all.table", index=None, sep='\t')
    return df

#join_total_numt_bp_allchr()

# %% [markdown]
# # How big are the NUMTs we identify? Should we set a minimum fragment size?

# %%
def sort_chrom(df_in):
    '''
    Simple sorting of chromosomes by (human) numerical order, plus sex chromosomes.
    '''
    df = df_in.copy()
    df['short_Subject'] = df['CHR'].str.split("r").str[1]
    df['short_Subject'] = df['short_Subject'].str.split("_").str[0]
    # Index.
    dict_chr = {
        '1':1, '2':2, '2A':2.3, '2B':2.7, '3':3, '4':4, '5':5, '6':6, '7':7, '8':8, '9':9, '10':10, '11':11, '12':12, '13':13, '14':14, '15':15, '16':16, '17':17, '18':18, '19':19, '20':20, '21':21, '22':22, '23':23, '24':24, 'X':30, 'Y':31
        }
    df['index_Subject'] = df['short_Subject'].map(dict_chr).fillna('MISSING')
    return df.sort_values('index_Subject')


def sort_hsa(df_in):
    '''
    Sort HSA column in plotting order.
    '''
    df = df_in.copy()
    # Annotate with CHR and HSA.
    df = get_CHR_from_Subject(df)
    #Exclude unknown chromosomes from non-T2T assemblies.
    df = exclude_unknown_chr(df)
    # Convert human NC chr IDs to hsa numbers.
    df = nc_to_human_chr(df)
    # Index for HSA.
    dict_hsa = {
        'hsa1':1, 'hsa2':2, 'hsa3':3, 'hsa4':4, 'hsa5':5, 'hsa6':6, 'hsa7':7, 'hsa8':8, 'hsa9':9, 'hsa10':10, 'hsa11':11, 'hsa12':12, 'hsa13':13, 'hsa14':14, 'hsa15':15, 'hsa16':16, 'hsa17':17, 'hsa18':18, 'hsa19':19, 'hsa20':20, 'hsa21':21, 'hsa22':22, 'hsa23':23, 'hsa24':24, 'hsaX':30, 'hsaY':31, 'hsa17x5':4, 'hsa5x17':19, 'hsa2a':2, 'hsa2b':2.5, 'hsa2A':2, 'hsa2B':2.5,
        'chr1':1, 'chr2':2, 'chr2A':2.3, 'chr2B':2.7, 'chr3':3, 'chr4':4, 'chr5':5, 'chr6':6, 'chr7':7, 'chr8':8, 'chr9':9, 'chr10':10, 'chr11':11, 'chr12':12, 'chr13':13, 'chr14':14, 'chr15':15, 'chr16':16, 'chr17':17, 'chr18':18, 'chr19':19, 'chr20':20, 'chr21':21, 'chr22':22, 'chr23':23, 'chr24':24, 'chrX':30, 'chrY':31
        }
    df['index_HSA'] = df['HSA'].map(dict_hsa).fillna('MISSING')
    if df[df['index_HSA']=='MISSING'].shape[0] > 0:
        #return df[df['index_HSA']=='MISSING']
        sys.exit(f'### Missing values in HSA column ###')
    return df.sort_values('index_HSA')


def get_species_order(species):
    '''
    Sort species in plotting order.
    '''
    order_species = {
            'Human':0, 'Chimpanzee':1, 'Bonobo':2, 'S. orangutan':3, 'B. orangutan':4, 'Gorilla':5, 'Siamang':6
        }
    return order_species[species]


def get_assembly_order(assembly):
    '''
    Sort assembly in plotting order.
    '''
    order_assemblies = {
            'CHM13':0, 'GRCh38':1, 'mPanTro3':2, 'panTro6':3, 'mPanPan1':4, 'panPan3':5, 'mPonAbe1':6, 'mPonAbe1.hap1':7, 'mPonAbe1.hap2':8, 'ponAbe3':9, 'mPonPyg2':10, 'mPonPyg2.hap1':11, 'mPonPyg2.hap2':12, 'mSymSyn1':13, 'mGorGor1':14, 'gorGor6':15, 'gorilla_Y':16, 'bonobo_Y':17, 'sorang_Y':18
        }
    return order_assemblies[assembly]

def get_assembly_order_t2t(assembly):
    '''
    Sort assembly in plotting order.
    '''
    order_assemblies = {
            'mPanPan1':0, 'mPanTro3':1, 'CHM13':2, 'mGorGor1':3, 'mPonAbe1':4, 'mPonPyg2':5, 'mSymSyn1':6
        }
    return order_assemblies[assembly]


#nc_to_human_chr(get_merged_coords_strand())

#sort_hsa(get_merged_coords())

#get_species_order('Human')
#get_assembly_order('CHM13')



# %% [markdown]
# # Define main function

# %%
def main( skip_tables=False ):
    # Import BLAST results.
    df_tab = join_blast_results()
    # Export BLAST results as BED files.
    export_BED_numt_coords(df_in=df_tab)
    # Merge overlapping reads.
    subprocess.run("bash mergeCoords.sh blast_new", shell=True, check=True)
    # Import the BED-formatted LASTZ results.
    df_bed = get_merged_coords()
    
    # Tables of summary info.
    compute_total_numt_bp(df_in=df_bed)
    print('table fini')
    compute_total_numt_bp_by_chr(df_in=df_bed)
    print('table byCHR fini')
    compute_total_numt_bp_chrY(df_in=df_bed)
    print('table chrY fini')
    join_total_numt_bp_allchr()
    print('table joinAll fini')

if __name__ == "__main__":
    main()

# %%



