import os,sys
from shutil import copyfile
import shutil
import uuid
import xlsxwriter
import pandas as pd
import math

sys.path.append("/home/nj8612/bin/XlsxWriter-1.1.2")

inputDir=sys.argv[1]
meta=sys.argv[2]
cmdfile=inputDir+"/Vcmd.txt"
outputExcel = inputDir+"/Results_"+meta+".xlsx"
writer = pd.ExcelWriter(outputExcel, engine='xlsxwriter')
workbook  = writer.book
cell_format = workbook.add_format()
cell_format.set_font_name('Courier New')
cell_format.set_font_size(10)
inputFa = inputDir+"/input.fa"
outputFa = inputDir+"/output.fa"

excelRow={'1':'A','2':'B','3':'C','4':'D','5':'E','6':'F','7':'G','8':'H','9':'I','10':'J','11':'K'}

def isNaN(num):
	return num != num
	
def addCellToExcel(cellnum,data,allcellnum):
	if isNaN(data):
		worksheet.write_blank(cellnum, None,cell_format)
		worksheet1.write_blank(allcellnum, None,cell_format)
	else:
		if type(data) == int or type(data) == float:
			worksheet.write_number(cellnum, float(data),cell_format)
			worksheet1.write_number(allcellnum, float(data),cell_format)
		else:
			worksheet.write(cellnum, data,cell_format)
			worksheet1.write(allcellnum, data,cell_format)

def colorHeader(ws,df):
	header_format = workbook.add_format({
	'bold': True,
	'fg_color': '#ffcccc',
	'border': 1})
	for col_num, value in enumerate(df.columns.values):
		ws.write(0, col_num, value, header_format)
		column_len = df[value].astype(str).str.len().max()
		column_len = max(column_len, len(value)) + 3
		ws.set_column(col_num, col_num, column_len)			

def colorMutation(seq,allpos,num,allrows):
	red = workbook.add_format({'color': 'red'})
	red.set_font_name('Courier New')
	red.set_font_size(10)
	
	seqlen = len(seq)
	start = allpos[0]
	part1 = seq[0:start]
	
	
	string_parts = []
	string_parts.append(cell_format)
	string_parts.append(part1)

	for p in range(0,len(allpos)):
		if p > 0:
			part = seq[allpos[p-1]+1:allpos[p]]
			if part:
				string_parts.append(cell_format)
				string_parts.append(part)

		part = seq[allpos[p]:allpos[p]+1]
		string_parts.append(red)
		string_parts.append(part)
		
	partn = seq[allpos[-1]+1:]
	string_parts.append(cell_format)
	string_parts.append(partn)
	
	cellnum = "E"+str(num)
	allcellnum="E"+str(allrows)
	worksheet.write_rich_string(cellnum, *string_parts)
	worksheet1.write_rich_string(allcellnum, *string_parts)

def colorTarget(seq,pos,length,num,colr,allrows):
	green = workbook.add_format({'color': colr, 'bold' : True })
	green.set_font_name('Courier New')
	green.set_font_size(10)
	
	seqlen = len(seq)
	part1 = seq[0:pos]
	part2 = seq[pos:pos+length]
	part3 = seq[pos+length:]

	string_parts = [ cell_format, part1, green, part2, cell_format, part3]
	cellnum = "E"+str(num)
	allcellnum="E"+str(allrows)
	worksheet.write_rich_string(cellnum, *string_parts)	
	worksheet1.write_rich_string(allcellnum, *string_parts)

def addColumnNames(df,allRows):
	lnum=1
	cols=list(df.columns.values)
	for i in range(0,11):
		data=cols[i]
		key = i+1
		cellnum = excelRow[str(key)]+str(lnum)
		addCellToExcel(cellnum,data,cellnum)


settingsFile = inputDir + "/Settings.txt"
df1 = pd.read_csv(settingsFile, sep="\t")
df1.to_excel(writer,sheet_name='0_CRIA-Settings',index=False)
worksheet = writer.sheets['0_CRIA-Settings']
colorHeader(worksheet,df1)

summaryFile = inputDir + "/AllStats.txt"
df1 = pd.read_csv(summaryFile, sep="\t")
df1 = df1.sort_values(by=['Biotracker', 'Experiment'])
df1.to_excel(writer,sheet_name='1_Exp-Summary',index=False)
worksheet = writer.sheets['1_Exp-Summary']
colorHeader(worksheet,df1)

statsFile = inputDir + "/RunStats.txt"
df2 = pd.read_csv(statsFile, sep="\t")
df2.sort_values(by=['SampleID', 'SampleName'])
df2.to_excel(writer,sheet_name='2_CRIA-Stats',index=False)
worksheet = writer.sheets['2_CRIA-Stats']
colorHeader(worksheet,df2)

worksheet1 = workbook.add_worksheet("3_All-Mutations")

totRows = 0
fileNum = 0
num = 0
with open(cmdfile,"r") as f:
	for line in f:
		fileNum = fileNum + 1
		arr = line.rstrip().split(' ')
		inputFile = arr[2].split('/')[-1]
		inputName = inputFile.split('_')[0]
		inputExcelName = arr[2].split('/')[-1].replace(".csv",".xlsx")
		inputExcel = inputDir + "/Filtered." + inputExcelName 
		edit = arr[3]
		ref=arr[4]
		target=arr[5]
		try:
			worksheet = workbook.add_worksheet(inputName)
		except:
			num = num + 1
			name = inputName+"-"+str(num)
			worksheet = workbook.add_worksheet(name)
		cell_format = workbook.add_format()
		cell_format.set_font_name('Courier New')
		cell_format.set_font_size(10)
		tar_loc=ref.index(target)
		tar_len=len(target)
		limit = int(tar_loc) + int(tar_len) + 10
		xlsxFile = pd.ExcelFile(inputExcel)
		lnum=1
		for sheet in xlsxFile.sheet_names:
			df = xlsxFile.parse(sheet)
			df = df.drop(['TargetSite','AlleleSequence'],axis=1)
			rows = len(df)
			cols = len(df.columns)
			for r in range(0,rows):
				allRows = totRows + r 
				if lnum == 1:
					addColumnNames(df,allRows+2)
				lnum=lnum+1
				for i in range(0,11):
					try:
						data=df.iloc[r].iloc[i]
					except:
						data=""
					key = i+1
					allCellNum = excelRow[str(key)]+str(allRows+2)
					cellnum = excelRow[str(key)]+str(lnum)
					addCellToExcel(cellnum,data,allCellNum)
				if df.iloc[r].iloc[3] == "Wild Type":
					ref=df.iloc[r].iloc[4]
					if target in ref:
						tar_loc=ref.index(target)
						tar_len=len(target)
						colorTarget(ref,tar_loc,tar_len,lnum,'green',allRows+2)
				elif (df.iloc[r].iloc[9] == "yes") and edit in df.iloc[r].iloc[4]:
					pos=df.iloc[r].iloc[4].index(edit)
					colorTarget(df.iloc[r].iloc[4],pos,len(edit),lnum,'blue',allRows+2)
				elif (df.iloc[r].iloc[9] == "yes"):
					seq = df.iloc[r].iloc[4]
					oname = str(lnum)+".fa"
					fasta= open(inputFa,"w+")
					fasta.write(">Ref\n%s\n>AlleleSeq\n%s\n" % (ref,seq))
					fasta.close()
					muscmd = "muscle -in " + inputFa + " -out " + outputFa + " -quiet -diags"
					os.system(muscmd)
					tname=str(lnum)+".fa"
					seq=""
					refseq=""
					with open(outputFa,"r") as fi:
						for faline in fi:
							faline=faline.rstrip()
							if faline.startswith(">Ref"):
								seq=""
							elif faline.startswith(">AlleleSeq"):
								refseq=seq
								seq=""
							else:
								seq = seq+faline
					fi.close()
					
					ref_list = list(refseq)
					allele_list = list(seq)
					ind = -1
					allpos = []
					
					for r,a in zip(ref_list, allele_list):
						ind = ind + 1
						
						if r != a and ind < limit:
							allpos.append(ind)
					
					colorMutation(seq,allpos,lnum,allRows+2)
		totRows = totRows + rows
		colorHeader(worksheet,df)
		colorHeader(worksheet1,df)

f.close()
workbook.worksheets_objs.sort(key=lambda x: x.name)
workbook.close()
os.remove(inputFa)
os.remove(outputFa)
