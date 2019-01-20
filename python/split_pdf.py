import argparse
import os
from PyPDF2 import PdfFileReader, PdfFileWriter

#Parse arguments

parser = argparse.ArgumentParser(description='Split PDF with many figures into individual pdf')
parser.add_argument('-i', '--input',
                        help='Name of PDF file to split',
                        required='True')
parser.add_argument('-n', '--nMain',
                        help='Number of main figures in the pdf file',
                        required='True',
                        type=int)
parser.add_argument('-s', '--suffix',
                        help='Suffix of output file name',
                        required='True')
args = parser.parse_args()

#Read PDF

inputpdf = PdfFileReader(open(args.input, "rb"))

if not os.path.exists("Output"):
    os.makedirs("Output")

def getOneFigure(inputpdf,nFigures):
	for i in range(inputpdf.getNumPages()):
		output = PdfFileWriter()
		output.addPage(inputpdf.getPage(i))
		if(i<nFigures):
			with open("Output/" + args.suffix + "%s.pdf" % str(i+1), "wb") as outputStream:
				output.write(outputStream)
		else:
			with open("Output/" + args.suffix + "S%s.pdf" % str(i-(nFigures-1)), "wb") as outputStream:
				output.write(outputStream)
			
getOneFigure(inputpdf,args.nMain)