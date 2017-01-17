#!/usr/bin/python3.4

import cgitb; cgitb.enable()
import os, string, sys, cgi
import binascii
#import binascii, mkdir
sys.path.insert(0, '/nfshome/agd2q/local/lib/python3.4/site-packages/biopython-1.65-py3.4-linux-x86_64.egg/')
sys.path.insert(0, '/nfshome/hcarroll/public_html/apps/clustalOmega/bin/')
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import ClustalOmegaCommandline
import itertools
import colorsys
import argparse #import argaparse library for addArgv, --afa
sys.stderr = sys.stdout

UTR_EXON = -2

form = cgi.FieldStorage()
path = ""

#creating an object to hold all the information for each gene
class Gene(object):
    def __init__(self, name):
        #name of gene
        self.name = name

        #array that will hold exon lengths 
        self.ExonCount = []

        # length of AA (unaligned)
        self.AA_len = 0

        # string to store the translated (unaligned) gene
        self.AA = ''

        #string to store the CDS region and lengths
        self.CDS = ''
        self.CDS_len = []

        # string to store the translated (aligned) gene
        self.Aligned_str = ''

        #string to store the frame string
        self.Frame_str = ''

        #array to keep track of exon frames
        self.Frames = []

        #line color for this gene
        self.color = ''

dnaToProtD = {
    "TTT":'F', "TTC":'F', "TTA":'L', "TTG":'L',
    "CTT":'L', "CTC":'L', "CTA":'L', "CTG":'L',
    "ATT":'I', "ATC":'I', "ATA":'I', "ATG":'M',
    "GTT":'V', "GTC":'V', "GTA":'V', "GTG":'V',

    "TCT":'S', "TCC":'S', "TCA":'S', "TCG":'S',
    "CCT":'P', "CCC":'P', "CCA":'P', "CCG":'P',
    "ACT":'T', "ACC":'T', "ACA":'T', "ACG":'T',
    "GCT":'A', "GCC":'A', "GCA":'A', "GCG":'A',

    "TAT":'Y', "TAC":'Y', "TAA":'', "TAG":'',
    "CAT":'H', "CAC":'H', "CAA":'Q', "CAG":'Q',
    "AAT":'N', "AAC":'N', "AAA":'K', "AAG":'K',
    "GAT":'D', "GAC":'D', "GAA":'E', "GAG":'E',

    "TGT":'C', "TGC":'C', "TGA":'', "TGG":'W',
    "CGT":'R', "CGC":'R', "CGA":'R', "CGG":'R',
    "AGT":'S', "AGC":'S', "AGA":'R', "AGG":'R',
    "GGT":'G', "GGC":'G', "GGA":'G', "GGG":'G'
}

def dna_to_prot(strand):
    aminos = [ dnaToProtD[strand[i:i+3] ] for i in range(0, len(strand), 3) ]
    protein = "".join(aminos)
    return protein

# make a new directory with a random name (in the files directory)
def mkDir():
    if not os.path.exists("files"):
        os.makedirs("files")
        #os.chmod(path, 0o755)
        os.chmod("files", 0o777) # DEBUGGING
        
    newDir = binascii.hexlify(os.urandom(16)).decode()
    global path
    path = "files/" + newDir
    if not os.path.exists(path):
        os.makedirs(path)
        #os.chmod(path, 0o755)
        os.chmod(path, 0o777) # DEBUGGING
	
# Argument afa used for checking if there is a multiple sequence alignment supplied from the command line
def dataProcess(afa): 
    # get the files in the directory
    listing = os.listdir(path)

    # create a dictionary to hold the Gene objects
    gene_dic = {}

    # loop through every file in the directory
    for infile in listing:

        # open the file
        with open(os.path.join(path, infile), 'r') as genes:
            # initalize the counter to zero for each gene/file
            CDSExon = ''
            CDSExonCount = 0


            # get the name of the gene (from the filename) and create an empty object
            name = os.path.splitext(infile)[0]
            gene = Gene(name)


            #go through the file get the gene name and then
            #get the length of all exons and put into array
            for seq_record in SeqIO.parse(genes, "fasta"):
                startFrame = ''
                stopFrame = ''

                #add the exon length into ExonCount array
                gene.ExonCount.append(len(seq_record))

                #convert Bio.SeqRecord to fasta format to get the sequence
                #in string format to check letter casing
                exon = seq_record.format("fasta")
                #this removes the first line of the fasta file and all newlines
                exon = exon.split('\n', 1)[-1].replace('\n', '')

                exonLength = len(exon)
		# remove the UTR (lowercase letters)
                CDS = exon.strip(string.ascii_lowercase)
                CDSExon += CDS
                CDSlength = len(CDS)
                gene.CDS_len.append(CDSlength)
                if(CDSExonCount == 0):
                    startFrame = '^'
                else:
                    startFrame = gene.Frames[CDSExonCount-1][2]  # get the ORF from the previous exon
		    # adjust the coding length to account for the reading frame from the previous exon
                    if(startFrame == '<'):
                        CDSlength -= 2
                    elif(startFrame == '>'):
                        CDSlength -= 1
			
		# determine the reading frame
                if(CDSlength % 3 == 0):
                    stopFrame = '^'
                elif(CDSlength % 3 == 1):
                    stopFrame = '<'
                else:
                    stopFrame = '>'

                gene.Frames.append([CDSExonCount, startFrame, stopFrame, CDSlength, exonLength])
                CDSExonCount += 1

            # trim off 1-2 nucleotides if the coding region is not a multiple of 3
            if len(CDSExon) % 3 != 0:
                excessNucleotides = len(CDSExon) % 3
                CDSExon = CDSExon[:-excessNucleotides]
                
            gene.CDS = CDSExon
            gene.AA = dna_to_prot(CDSExon)
            gene.AA_len = len(gene.AA)

            #create a directory to store AA fasta files
            aa_path = path + '/AA/'
            if not os.path.exists(aa_path):
                os.makedirs(aa_path)
                os.chmod(aa_path, 0o777) # DEBUGGING

	    # FUTURE: open the file once, append to it in this loop, then close it once
            #create one giant fasta file with all AA sequences 
            #to pass to clustal omega
            with open(path + "/ALL.txt", "a") as ALL_AA:
                ALL_AA.write('>' + name + '\n')
                ALL_AA.write(gene.AA + '\n')
                ALL_AA.close()

            # #create AA fasta files
            # with open(os.path.join(path + '/AA/' + name + '_AA.txt'), 'w') as AA_file:
            #     AA_file.write('>' + name + '\n')
            #     AA_file.write(gene.AA + '\n')
            #     AA_file.close()

	    # FUTURE:
            #NEED TO ASK ABOUT THIS...IS IT 100 AA OR 100 NUCLEOTIDES
            #Check to make sure all AA sequences are at least 100 aa


        gene_dic[name] = gene

    # call clustal omega
    in_file = path + "/ALL.txt"

    # if there is --afa on the command line, use the MSA from user
    # else use the aligned from the path in file folder
    out_file = ""
    if afa:
        out_file = "./" + afa
    else:
        out_file= path + "/aligned.fasta"
        clustalo = ClustalOmegaCommandline(infile=in_file, outfile=out_file, auto=True, cmd="/nfshome/hcarroll/public_html/apps/clustalOmega/bin/clustalo")
        clustalo()

    #
    # read in the aligned sequences and store them in their respective objects
    #
    with open(out_file, 'r') as clustal:
        for line in clustal:
            line = line.rstrip()
            if(line[0] is '>'): # and count is 0):
                name = line[1:] # remove the > character at the beginning of the line
            else:
                gene_dic[name].Aligned_str += line

    if afa:
        # verify that the user-supplied alignment file matches exactly with the translated versions of the user-supplied input files
        # for each gene
        #     compare AA with Aligned_str (except for "-"s) (by making a temporary copy of Aligned_str that has all "-"s replaced with "")
        for name in gene_dic:
            gene = gene_dic[name]
            if gene.AA != gene.Aligned_str.replace("-",""):
                print( "Translated version of input files for", name, "does not match with the supplied alignment:")
                print( "Translated:", gene.AA)
                print( "Alignment: ", gene.Aligned_str.replace("-",""), "(any gaps were removed for this error message)")
                sys.exit(1)    

    for name, gene in gene_dic.items():
        count = 0  # alignment index
        exon = 0   # index of the current        exon for this gene
        start = False  # first time?
        
        if exon >= len( gene.CDS_len):
            # no valid exons, go to the next gene
            continue
            
        frame_len = int(gene.CDS_len[exon]/3)

	# process all of the 5' exons that are entirely UTR (i.e., have at most 2 nucleotides)
        while frame_len == 0:  # only UTR
            gene.Frames[exon].append(UTR_EXON)
            exon += 1
            if exon < len( gene.CDS_len):
                frame_len = int(gene.CDS_len[exon]/3)
            else:
                break

	# record the alignment index for
        # 1) the first AA and
        # 2) the last AA for each exon    
        for char in gene.Aligned_str:
            if char is '-':
                pass
            else:
                if start is False:
                    start = True
                    # FUTURE: can this be pulled out of the for loop (i.e., can the frame_len == 0 for the first AA)?
                    gene.Frames[exon].append( count )  # record the index of the alignment
                elif frame_len is 0:  # done processing this exon
                    exon += 1
                    gene.Frames[exon].append( count )  # record the index of the alignment
                    frame_len = int(gene.CDS_len[exon]/3)
                    while frame_len == 0: # only UTR 
                        gene.Frames[exon].append(UTR_EXON)
                        exon += 1
                        frame_len = int(gene.CDS_len[exon]/3)
                else:
                    frame_len -= 1
            count += 1
            
	# process all of the 3' exons that are entirely UTR
        while exon < len(gene.CDS_len)-1:
            exon += 1
            gene.Frames[exon].append(UTR_EXON)

    # Make Final 2D array
    CombinedLists = []  # 1st index: genes; 2nd index items are: name, CDSExonCount, startFrame, stopFrame, CDSlength, exonLength, then AA alignment indicies
    geneIndex = 0
    # each gene is a row in CombinedLists
    for name, gene in gene_dic.items():
        CombinedLists.append([name])
        for col in gene.Frames:
            CombinedLists[geneIndex].append(col)
        geneIndex += 1

    preFinal = []
    count = 0
    for frameInfoList in CombinedLists:
        name1 = frameInfoList[0]
        preFinal.append([name1])
        for col in frameInfoList[1:]:
            rowCount = count + 1
            rowOn = 0
            exonNum = col[0]  # coding exon
            startFrame = col[1]
            endFrame = col[2]
            # = col[3]
            exonSize = col[4] # coding length
            alignmentIndex_firstExon = col[5]
            
            if alignmentIndex_firstExon is UTR_EXON:  # if the first exon is only a UTR
                preFinal[count].append([alignmentIndex_firstExon, exonSize])
            elif count == 0:
                preFinal[count].append([alignmentIndex_firstExon, exonSize])
                for Nrow in CombinedLists[rowCount:]:
                    name2 = Nrow[0]
                    for Ncol in Nrow[1:]:
                        if startFrame == Ncol[1] and endFrame == Ncol[2]:  # start and end frames match
                            if (abs(exonSize - Ncol[4]) % 3) is 0: # lengths are the same ORF
                                if abs(alignmentIndex_firstExon - Ncol[5]) <= 10:
                                    # each of first alignment positions for each gene are within 10 of each other
                                    preFinal[count][exonNum+1].append([name2, Ncol[5], Ncol[4], Ncol[0]])
            else:
                found = False
                for row1 in preFinal:
                    if found is False:
                        lastName = row1[0]
                        colOn = 0
                        for col1 in row1[1:]:
                            if len(col1) > 2:
                                count3 = 0
                                while count3 < len(col1)-2 and found is False:
                                    if name1 == col1[count3+2][0] and exonNum == col1[count3+2][3]:
                                        found = True
                                        preFinal[count].append([-1, [rowOn, colOn, col1[0]]])
                                    count3 += 1
                            colOn += 1
                    rowOn += 1

                if found is False:
                    preFinal[count].append([alignmentIndex_firstExon, exonSize])
                    for Nrow in CombinedLists[rowCount:]:
                        name2 = Nrow[0]
                        for Ncol in Nrow[1:]:
                            if startFrame == Ncol[1] and endFrame == Ncol[2]:
                                if (abs( exonSize - Ncol[4]) % 3) is 0:
                                    if abs(alignmentIndex_firstExon - Ncol[5]) <= 10:
                                        preFinal[count][exonNum+1].append([name2, Ncol[5], Ncol[4], Ncol[0]])
            rowCount += 1
        count +=1

    # print("preFinal")
    # for row in preFinal:
    #     print(len(row))
    #     print(row)
    # print()

    ###DEBUGGING SECTION TO SEE THE 2D ARRAY OF EXON LENGTHS
    i = 1
    x = 100
    y = 20
    temp = '<svg width=300% height=300%>\n'

    #find out how many nonCDS exons there to find how far to shift SVG
    maxStart = 0
    for row in preFinal:
        startNum = 0
        for col in row[1:]:
            if col[0] == UTR_EXON:
                startNum += 1
        if startNum > maxStart:
            maxStart = startNum
    # print('max start is', maxStart)


    #create random colors for each gene using colorsys library
    numOfGenes = len(preFinal)
    HSV_tuples = [(x*1.0/numOfGenes, 1, 1) for x in range(numOfGenes)]
    RGB_tuples = map(lambda x: colorsys.hsv_to_rgb(*x), HSV_tuples)
    RGB = []
    for rgb in RGB_tuples:
        temp = list(rgb)
        temp = [int(color*255) for color in temp]
        temp = ','.join(str(x) for x in temp)
        RGB.append('rgb('+temp+')')

    with open("createSVGtemp.html", "w") as printSVG:
        printSVG.write("<html>\n")
        printSVG.write("<body>\n")
        printSVG.write('<svg width=\"1000%\" height=\"300%\" style=\"overflow-x: auto; overflow-y: auto;\">\n')

        colCount = 0
        xStart = 100
        yStart = 130
        xText = 0
        yText = 0
        xPrev = 0
        yLine = 0
        xLine = 0
        xLineStart = 0
        yLineStart = 0
        xLineStop = 0
        yLineStop = 0
        rowCount = 0
        prevExonCount = 0
        prevNumOfExons = 0
        exonCount = 0
        firstX = 0
        for row in preFinal:
            y += yStart
    #        if prevNumOfExons > 2:
    #            y += ((prevNumOfExons-1)*3)
            x = xStart
            # print('xStart is ', xStart)
            name = row[0]
            colCount = 0
            color = RGB[rowCount]
            printSVG.write('\n<!--'+name+'-->\n')
            # print(name)
            prevNumOfExons = -1
            prevExonCount = 0
            prevDy = -1
            for col in row[1:]:
                count = 0
                if colCount is 0:
                    if col[0] >= -1:
                        x = xStart + (maxStart * 75)
                        firstX = x + 50
                        # print('firstx is', firstX)
                    else:
                        firstX = xStart + (maxStart * 75) + 50
                elif colCount > 0:
                    if col[0] > -1:
                        x = 3*col[0]
                        # print('col[0] is', col[0], x, firstX)
                        if col[0] < firstX:
                            x = firstX + 60
                        if x < prevX+50:
                            x = prevX + 100
                    elif  col[0] is -1:
                        x = 3*col[1][2]
                        # print('col[1][2]', col[1][2], x, firstX)
                        if col[1][2] < firstX:
                            x = firstX + 60
                        if x < prevX+50:
                            x = prevX + 100
                    elif col[0] is UTR_EXON:
                        x = prevX + 60
                    printSVG.write('<text x=\"'+str(25)+'\" y=\"'+str((y+(y+50))/2)+'\" style=\"stroke:'+color+'\">\"'+name+'\"</text>')
     #               printSVG.write('<line x1=\"'+str((prevX+50))+'\" y1=\"'+str((yLineStart))+'\"x2=\"'+str((x))+'\" y2=\"'+str((yLineStop))+'\" style=\"stroke:'+color+'; stroke-width:2\"/>\n')
                numOfExons = len(col)-1
                dy = 50 + (50 * 0.25 * (numOfExons-1))
                if col[0] == UTR_EXON:
                    dasharray = 'stroke-dasharray: 10 5;'
                else:
                    dasharray = ''
                if col[0] != -1:
                    printSVG.write('<rect x=\"'+str(x)+'\" y=\"'+str(y)+'\" rx=\"10\" ry=\"10\" width=\"50\" height=\"'+str(dy)+'\" style=\"fill:white;stroke:rgb(0,0,0);stroke-width:2;'+dasharray+'\" />\n')
                    while count < numOfExons:
                        xText = x + 10
                        if numOfExons == 1:
                            yText = (y + (y+dy))/2
                        elif count is 0:
                            yText = y + dy/numOfExons
                        if count == 0:
                            printSVG.write('<text x =\"'+str(xText)+'\" y=\"'+str(yText)+'\">'+str(col[1])+'</text>\n')
                        else:
                            yText += 15
                            printSVG.write('<text x =\"'+str(xText)+'\" y=\"'+str(yText)+'\">'+str(col[count+1][2])+'</text>\n')
                        count += 1
                if col is not row[1]:
                    exonCount = 0
                    yLine = y
                    prev = ''
                    if col[0] is -1:
                        yLine = 20 + (yStart * (col[1][0] + 1))
                        prev = preFinal[col[1][0]][col[1][1] + 1]
                        numOfExons = len(prev) -1
                        dy =50 + (50 * 0.25 * (numOfExons-1))


                    if numOfExons > 0:
                        for i in range(len(prev)-2):
                            if prev[i+2][0] is name:
                                exonCount = i + 1
                    if colCount is 1:
                        prevExonCount = 0
                    yLineStart = prevY + ((prevDy/(prevNumOfExons+1)) * (prevExonCount+1))
                    yLineStop = yLine + ((dy/(numOfExons+1))) * (exonCount+1)
                    printSVG.write('<line x1=\"'+str((prevX+50))+'\" y1=\"'+str((yLineStart))+'\"x2=\"'+str((x))+'\" y2=\"'+str((yLineStop))+'\" style=\"stroke:'+color+'; stroke-width:2\"/>\n')

                colCount +=1
                if col[0] is not -1:
                    prevNumOfExons = numOfExons
                    prevDy = dy
                    prevY = y
                    prevX = x
                    prevExonCount = exonCount
                elif col[0] is -1:
                    prevNumOfExons = numOfExons
                    prevDy = dy
                    prevY = yLine
                    prevX = x
                    prevExonCount = exonCount


            rowCount += 1

        printSVG.write("\n</body>\n")
        printSVG.write("</html>")


def main():
    # os.system("chmod -R 777 files")
    # os.system("echo 'DEBUGGING:' > /tmp/hcarroll.tmp; chmod 777 /tmp/hcarroll.tmp")

    
    message = ""

    mkDir()
    #declaration for argparse
    useArgv = False
    parser = argparse.ArgumentParser()
    parser.add_argument("-a", "--afa", metavar="", help="provide aligned file for cluster mega")
    parser.add_argument("files",nargs="*", help="fasta file(s)")
    args = parser.parse_args()

    #
    # get input file(s) from either the CGI form or from the command-line
    #    
    if 'GATEWAY_INTERFACE' in os.environ:  # Called from a CGI form
        fileItems = form['filename[]']
        for fileItem in fileItems:
            if fileItem.file:
                fn = os.path.basename(fileItem.filename.replace("\\", "/"))  # change Windows filenames
                try:
		    # copy the file
                    open(path + '/' + fn, 'wb').write(fileItem.file.read())
                except:
		    # display error message
                    print("""Content-Type: text/html\n\n
                          <html>
                          <body>
                              <p>%s</p>
                          </body>
                          </html>
                          """ %(path+'/'+fn))
                    sys.exit()
    else:
        # get file(s) from the command-line
        useArgv = True
        fileItems = args.files          #if using the argument line
        for fileItem in fileItems:
            op = open(fileItem,'r').read()
            fn = os.path.basename(fileItem)
            try:
                # Need to change 'wb' to 'w' due to this error message
                # 'str' does not support the buffer interface
                open(path + '/' + fn, 'w').write(op)
            except:
		# display error message
                print("""Content-Type: text/html\n\n
                      <html>
                      <body>
                         <p>%s</p>
                      </body>
                      </html>
                      """ %(path+'/'+fn))
                sys.exit()
                     
    afa = args.afa      #get the aligned provide, or set it as False if user dont provide the aligned
    dataProcess(afa)

    redirectStr="""Content-Type: text/html\n\n
<html>
    <head>
        <meta http-equiv="refresh" content="0; url=createSVGtemp.html" />
    </head>
    <body>
        <p>Redirecting to <a href="createSVGtemp.html">createSVGtemp.html</a></p>
    </body>
</html>
"""
    if(not useArgv):
        os.system("echo 'DEBUGGING: redirectStr: " + redirectStr + "' > /tmp/hcarroll.tmp; chmod 777 /tmp/hcarroll.tmp")
        print(redirectStr)
    else:
        print("Successfully completed")

if __name__ == "__main__":
    main() 
    sys.exit(0)
