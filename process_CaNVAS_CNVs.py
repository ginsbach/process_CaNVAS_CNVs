#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import getopt
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'
from datetime import datetime
import argparse

def createOutputs(outDF, outputfile):
	print("Generating Outputs")
	for cnv in outDF['CNVID'].unique():
		outpath="output/"+cnv+"_"+outputfile+".txt"
		rel = outDF[outDF.CNVID ==cnv]
		rel[["CNVID","SubjectID","meanLRRIn","meanLRROutleft","meanLRROutright","medianLRRIn","medianLRROutleft","medianLRROutright","BAFHomozygote","BAFTrisomic","BAFDisomic","SNPsIn", "SNPsOutleft","SNPsOutright"]].to_csv(outpath, sep="\t", index=False)
		outPlotName="output/plots/"+cnv+"_"+outputfile+".png"
		x = rel["meanLRRIn"]
		y1 = (rel["meanLRROutleft"] + rel["meanLRROutright"])/2
		y2 = rel["BAFHomozygote"]
		y3 = rel["BAFTrisomic"]
		y4 = rel["BAFDisomic"]
		fig, (ax1, ax2, ax3, ax4) = plt.subplots(nrows=4, sharex=True)
		ax1.scatter(x, y1)
		ax1.grid()
		ax1.set_title('Mean LRR Out vs Mean LRR In',fontsize=8)
		ax2.scatter(x, y2)
		ax2.grid()
		ax2.set_title('BAF Homozygote vs Mean LRR In',fontsize=8)
		ax3.scatter(x, y3)
		ax3.grid()
		ax3.set_title('BAF Trisomic vs Mean LRR In',fontsize=8)
		ax4.scatter(x, y4)
		ax4.grid()
		ax4.set_title('BAF Disomic vs Mean LRR In',fontsize=8)
		fig.suptitle(cnv, fontsize=10)
		plt.tight_layout()
		fig.set_size_inches(7, 10)
		fig.savefig(outPlotName, dpi=200)
		plt.close()

def renameColumns(columnHeader):
    if "Log R Ratio" in columnHeader:
        return  "logRratio"
    elif "B Allele Freq" in columnHeader:
        return "BAF"
    elif "GType" in columnHeader:
        return "Genotype"
    else:
        return columnHeader

def readInputs(inputfile, inputpath):
	cnvInputFile=pd.read_csv(inputfile, delim_whitespace=True)
	cnvInputFile['CNVID']=cnvInputFile['Chr'].map(str)+"_"+cnvInputFile['CNVStart'].map(str)+"_"+cnvInputFile['CNVEnd'].map(str)
	outDF= pd.DataFrame()
	print("Processing CNVs")
	with os.scandir(inputpath) as entries:
	   for entry in entries:
	    	print(entry)
	    	genotypeInput=pd.read_csv(entry, sep='\t', low_memory=False)
	    	genotypeInput.rename(mapper = renameColumns, axis = 'columns', inplace = True)
	    	genotypeInput[["Chr","Position","BAF","logRratio"]] = genotypeInput[["Chr","Position","BAF","logRratio"]].apply(pd.to_numeric, errors='coerce')
	    	bothInput=pd.merge(genotypeInput,cnvInputFile,on='Chr')
	    	bothInput= bothInput.loc[bothInput['Genotype'] != "NC"]
	    	inCNV=bothInput[bothInput["Position"].between(bothInput["CNVStart"], bothInput["CNVEnd"])]
	    	leftCNV=bothInput[bothInput["Position"].between(bothInput["LeftStart"], bothInput["LeftEnd"])]
	    	rightCNV=bothInput[bothInput["Position"].between(bothInput["RightStart"], bothInput["RightEnd"])]
	    	SNPin_size=inCNV.groupby(['CNVID']).size()
	    	SNPleft_size=leftCNV.groupby(['CNVID']).size()
	    	SNPright_size=rightCNV.groupby(['CNVID']).size()
	    	inCNV_sum=inCNV.groupby(['CNVID'])['logRratio'].median()
	    	leftCNV_sum=leftCNV.groupby(['CNVID'])['logRratio'].median()
	    	leftCNV_sum.rename("medLRROutleft", inplace=True)
	    	rightCNV_sum=rightCNV.groupby(['CNVID'])['logRratio'].median()
	    	rightCNV_sum.rename("medLRROutright", inplace=True)
	    	inCNV_mean=inCNV.groupby(['CNVID'])['logRratio'].mean()
	    	inCNV_mean.rename("meanLRRIn", inplace=True)
	    	leftCNV_mean=leftCNV.groupby(['CNVID'])['logRratio'].mean()
	    	leftCNV_mean.rename("meanLRROutleft", inplace=True)
	    	rightCNV_mean=rightCNV.groupby(['CNVID'])['logRratio'].mean()
	    	rightCNV_mean.rename("meanLRROutright", inplace=True)
	    	inCNV['BAF_bins']=pd.cut(inCNV.BAF,[0,0.005,0.29,0.37,0.47,0.53,0.63,0.71,0.995,1], include_lowest=True, ordered=False, labels=pd.Categorical(['Homozygote',"None","Trisomic","None","Disomic","None","Trisomic","None","Homozygote"]))
	    	BAFbin_size=inCNV.groupby(['CNVID',"BAF_bins"], as_index=False).count()
	    	relBAFbin_size=BAFbin_size[['CNVID', 'BAF_bins',"Name"]]
	    	relBAFbin_size2=relBAFbin_size.pivot(index='CNVID', columns='BAF_bins', values='Name')
	    	aggregate_df=pd.concat([SNPin_size, SNPleft_size, SNPright_size, inCNV_sum, leftCNV_sum, rightCNV_sum, inCNV_mean, leftCNV_mean, rightCNV_mean,relBAFbin_size2], axis=1)
	    	aggregate_df["SubjectID"]=entry
	    	outDF=outDF.append(aggregate_df)
	   outDF.index.name = 'CNVID'
	   outDF.reset_index(inplace=True)
	   outDF.rename(columns={outDF.columns[1]: "SNPsIn",outDF.columns[2]: "SNPsOutleft",outDF.columns[3]: "SNPsOutright",outDF.columns[4]:"medianLRRIn", outDF.columns[5]:"medianLRROutleft", outDF.columns[6]:"medianLRROutright",outDF.columns[7]:"meanLRRIn", outDF.columns[8]:"meanLRROutleft", outDF.columns[9]:"meanLRROutright", outDF.columns[10]: "BAFDisomic",outDF.columns[11]: "BAFHomozygote",outDF.columns[12]: "None",outDF.columns[13]: "BAFTrisomic",outDF.columns[14]: "SubjectID"}, inplace=True)
	print("Finished CNV Processing")
	return outDF

def main(argv):
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', type=str)
	parser.add_argument('-p', type=str)
	parser.add_argument('-o', type=str)
	args = parser.parse_args()
	inputfile = args.i
	inputpath = args.p
	outputfile = args.o
	outDF=[]
	print('Input file is "'+inputfile+'"')
	print('Input path is "'+inputpath+'"')
	print('Output file is "'+outputfile+'"')
	if not os.path.exists('output/'):
		os.makedirs('output/')
	if not os.path.exists('output/plots/'):
		os.makedirs('output/plots/')
	outDF= readInputs(inputfile, inputpath)
	createOutputs(outDF, outputfile)



if __name__ == "__main__":
   main(sys.argv[1:])