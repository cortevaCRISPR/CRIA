import os,sys
from shutil import copyfile
import uuid

sys.path.append("/home/nj8612/bin/XlsxWriter-1.1.2")
import xlsxwriter

inputFile = sys.argv[1]
edit = sys.argv[2]
ref=sys.argv[3]
targetS=sys.argv[4]
includeSNP=sys.argv[5]

tar_loc=ref.index(targetS)
tar_len=len(targetS)
limit = int(tar_loc) + int(tar_len) + 10

basename=os.path.basename(inputFile)
pathname=os.path.dirname(inputFile)
outname = basename.replace(".csv",".xlsx")
outputExcel = pathname+"/Filtered."+outname

workbook = xlsxwriter.Workbook(outputExcel, {'strings_to_numbers': True})
worksheet = workbook.add_worksheet()
cell_format = workbook.add_format()
cell_format.set_font_name('Courier New')
cell_format.set_font_size(10)

tmpname=str(uuid.uuid4())+".csv"
tmp=open(tmpname,"w+")
samplename = outname.split("_")[0]

def isNaN(num):
    return num != num
	
def addCellToExcel(cellnum,data):
	if isNaN(data):
		worksheet.write_blank(cellnum, None,cell_format)
	else:
		if type(data) == int or type(data) == float:
			worksheet.write_number(cellnum, float(data),cell_format)
		else:
			worksheet.write(cellnum, data,cell_format)

def addColumnNames(df):
	lnum=1
	cols=list(df.columns.values)
	for i in range(0,13):
		data=cols[i]
		key = i+1
		cellnum = excelRow[str(key)]+str(lnum)
		addCellToExcel(cellnum,data,i)
		
validLines=0
newtot=0
id=""
name=""
nhej=0
hdr=0
totMutations=0	

def colorMutation(seq,allpos,num):
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
	
	cellnum = "G"+str(num)
	worksheet.write_rich_string(cellnum, *string_parts)


with open(inputFile,"r") as f:
	for line in f:
		arr = line.rstrip().split(",")	
		if arr[0]:
			id=arr[0]
			#samplename=arr[0]
			name=arr[1]
		if len(arr) >= 12 and (arr[5] == "AlleleChange" or arr[5] == "Wild Type"):
			tmp.write(line)
			if arr[5] != "AlleleChange":
				validLines=validLines+1
				newtot=newtot+int(arr[7])
			id = 0
			#samplename=arr[0]
		elif len(arr) >= 12 and (arr[5] and arr[11] == "yes" and includeSNP == "0"):	
			if ";" in arr[5]:
				dat = arr[5].split(";")
				for i in range(0,len(dat)-1):
					if "SNP" in dat[i] or "outside" in dat[i] or "to" in dat[i]:
						a = 0
					else:
						tmp.write(line)
						validLines=validLines+1
						newtot=newtot+int(arr[7])
						a = 0
						break
			else:
				if "SNP" in arr[5] or "outside" in arr[5] or "to" in arr[5]:
					a = 0
				else:
					tmp.write(line)
					validLines=validLines+1
					newtot=newtot+int(arr[7])
					a = 0
		elif len(arr) >= 12 and (arr[5] and arr[11] == "yes" and includeSNP == "1"):
			tmp.write(line)
			validLines=validLines+1
			newtot=newtot+int(arr[7])
			a = 0
				
tmp.close()

excelRow={'1':'A','2':'B','3':'C','4':'D','5':'E','6':'F','7':'G','8':'H','9':'I','10':'J','11':'K','12':'L','13':'M'}

lnum=0
if validLines == 0:
	with open(tmpname,"r") as f:
		for line in f:
			lnum=lnum+1
			arr = line.rstrip().split(",")
			arr[8] = "100%"
			for i in range(0,12):
				data=arr[i]
				key = i+1
				cellnum = excelRow[str(key)]+str(lnum) 
				addCellToExcel(cellnum,data)
	print("%s\t%s\t%s\t%s\t%s" %(samplename,name,str(nhej),str(hdr),str(totMutations)))	
	os.remove(tmpname)
	workbook.close()
	sys.exit()


totMutations=0	
nhej=0
hdr=0
lnum=0
target=""
tar_loc=""
with open(tmpname,"r") as f:
	for line in f:
		lnum=lnum+1
		arr = line.rstrip().split(",")
		arr[8] = arr[8][:-1]
		if lnum == 2:
			target = arr[3]
		for i in range(0,13):
			try:
				data=arr[i]
			except:
				data=""
			if i == 8 and lnum != 1:
				data = round((float(arr[7])/float(newtot))*100,4)
				arr[8] = data
			elif i == 9 and lnum != 1:
				data = newtot
			key = i+1
			cellnum = excelRow[str(key)]+str(lnum) 
			addCellToExcel(cellnum,data)
		if arr[5] == "Wild Type":
			ref=arr[6]
			if target in ref:
				tar_loc=ref.index(target)
				tar_len=len(target)
			else:
				target = arr[3]
				if target in ref:
					tar_loc=ref.index(target)
					tar_len=len(target)
		elif (arr[11] == "yes") and edit in arr[6]:
			pos=arr[6].index(edit)
			hdr=hdr+float(arr[8])
			totMutations=totMutations+int(arr[7])
		elif (arr[11] == "yes"):
			seq = arr[6]
			oname = str(lnum)+".fa"
			fasta= open("input.fa","w+")
			fasta.write(">Ref\n%s\n>AlleleSeq\n%s\n" % (ref,seq))
			fasta.close()
			os.system("muscle -in input.fa -out output.fa -quiet -diags")
			tname=str(lnum)+".fa"
			seq=""
			refseq=""
			totMutations=totMutations+int(arr[7])
			nhej=nhej+float(arr[8])
			with open("output.fa","r") as fi:
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
			
			
			colorMutation(seq,allpos,lnum)

				
lnum=lnum+1
cellnum = "G" +str(lnum) 
addCellToExcel(cellnum,"Total Mutations")
cellnum = "H" +str(lnum) 
addCellToExcel(cellnum,totMutations)

lnum=lnum+1
cellnum = "G" +str(lnum) 
addCellToExcel(cellnum,"NHEJ Frequencies")
cellnum = "H" +str(lnum) 
data=str(nhej)+"%"
addCellToExcel(cellnum,data)

lnum=lnum+1
cellnum = "G" +str(lnum) 
addCellToExcel(cellnum,"HDR Frequencies")
cellnum = "H" +str(lnum) 
data=str(hdr)+"%"
addCellToExcel(cellnum,data)

print("%s\t%s\t%s\t%s\t%s" %(samplename,name,str(nhej),str(hdr),str(totMutations)))

lnum=0
with open(tmpname,"r") as f:
	for line in f:
		lnum=lnum+1
		arr = line.rstrip().split(",")
		if lnum > 2:
			val = round((float(arr[7])/float(totMutations))*100,4)
			cellnum = "M" +str(lnum) 
			addCellToExcel(cellnum,val)
		
f.close()
os.remove(tmpname)
workbook.close()
