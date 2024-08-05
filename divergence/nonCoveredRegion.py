import pandas as pd
import pybedtools
import sys

'''
This file takes in a part of maf file (lines marked with an 's')
with a single species with it's single chromsome. So for example
the file should look like:

s sorang.chr1 10 5 + 203644362 ATCCT
s sorang.chr1 20 7 - 203644362 GATCCTT

'''

mafDat = []
with open(sys.argv[1],"r") as data:
    for line in data:
        dat = line.rstrip().split()
        mafDat.append(dat[2:-2])

chrSpecies = dat[1].split(".")[1]
chrLen = int(dat[-2])

mafDat = pd.DataFrame(mafDat, columns=['MAF_START','LENGTH','STRAND'])
mafDat[["MAF_START", "LENGTH"]] = mafDat[["MAF_START", "LENGTH"]].apply(pd.to_numeric)
mafDat["CHR"] = chrSpecies

def calStart(maf_start, chrLen, strand):
    if strand == "-":
        start = chrLen - maf_start - 1
    elif strand == "+":
        start = maf_start
    return start

mafDat["START"] = mafDat.apply(lambda row: calStart(row["MAF_START"], chrLen, row["STRAND"]), axis=1)
mafDat["END"] = mafDat["START"] + mafDat["LENGTH"]
mafDat.drop(["MAF_START", "LENGTH","STRAND"], inplace=True, axis=1)

mafDat = mafDat.to_csv(sep='\t', header=False, index=False)

mafDatBed = pybedtools.BedTool(mafDat, from_string=True)
mafDatBed = mafDatBed.sort().merge()

mafDat = mafDatBed.to_dataframe()
mafDat["length"] = mafDat["end"] - mafDat["start"]

covered = mafDat["length"].sum()
uncovered = chrLen - covered

print(f"{dat[1]}\t{chrLen}\t{covered}\t{covered/chrLen*100:.2f}%\t{uncovered}\t{uncovered/chrLen*100:.2f}%")
