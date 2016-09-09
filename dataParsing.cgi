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

sys.stderr = sys.stdout

form = cgi.FieldStorage()
path = ""

#creating an object to hold all the information for each gene
class Gene(object):
    def __init__(self, name):
        #name of gene
        self.name = name

        #array that will hold exon lengths 
        self.ExonCount = []

        #length of AA
        self.AA_len = 0

        #string to store the translated gene
        self.AA = ''

        #string to store the CDS region and lengths
        self.CDS = ''
        self.CDS_len = []

        #string to store the aligned gene
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


#THIS SHOULD BE PASSED IN OR MAYBE JUST THE USER FILE 
#WILL BE PASSED IN

#def ParseData(path):
#path = 'files/'
#path = 'files/test/'
#path = 'files/smallTest/'
#path = 'files/testAll/'

def mkDir():
    newDir = binascii.hexlify(os.urandom(16)).decode()
    global path
    path = "files/" + newDir
    if not os.path.exists(path):
        os.makedirs(path)
        #os.chmod(path, 0o755)
        os.chmod(path, 0o777) # DEBUGGING

def dataProcess():
    #get the files in the directory
    listing = os.listdir(path)

    #initialize the 2d array as a 1d array for exon lengths 
    #geneExonCount = []

    #string to store the translated gene
    #AA = ""

    #create an empty object
    #gene = Gene()

    #create a dictionary to hold the Gene objects
    gene_dic = {}

    #counter for 2D array
    #row = 0

    #loop through every file in the directory
    for infile in listing:

        #open the file
        with open(os.path.join(path, infile), 'r') as genes:
            #initalize the counter to zero for each gene/file
            count = 0
            CDSExon = ''
            CDSExonCount = 0


            #get the name of the gene and create an empty object
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
                CDS = exon.strip(string.ascii_lowercase)
                CDSExon += CDS
                CDSlength = len(CDS)
                gene.CDS_len.append(CDSlength)
                if(CDSExonCount == 0):
                    startFrame = '^'
                    if(CDSlength % 3 == 0):
                        stopFrame = '^'
                    elif(CDSlength % 3 == 1):
                        stopFrame = '<'
                    else:
                        stopFrame = '>'
                else:
                    startFrame = gene.Frames[CDSExonCount-1][2]
                    if(startFrame == '<'):
                        CDSlength -= 2
                    elif(startFrame == '>'):
                        CDSlength -= 1
                    if(CDSlength % 3 == 0):
                        stopFrame = '^'
                    elif(CDSlength % 3 == 1):
                        stopFrame = '<'
                    else:
                        stopFrame = '>'

                gene.Frames.append([count, startFrame, stopFrame, CDSlength, exonLength])
                CDSExonCount += 1
                count += 1

            gene.CDS = CDSExon
            gene.AA = dna_to_prot(CDSExon)
            gene.AA_len = len(gene.AA)

            #create a directory to store AA fasta files
            aa_path = path + '/AA/'
            if not os.path.exists(aa_path):
                os.makedirs(aa_path)
                os.chmod(aa_path, 0o777) # DEBUGGING

            #create one giant fasta file with all AA sequences 
            #to pass to clustal omega
            with open(path + "/ALL.txt", "a") as ALL_AA:
                ALL_AA.write('>' + name + '\n')
                ALL_AA.write(gene.AA + '\n')
                ALL_AA.close()

            #create AA fasta files
            with open(os.path.join(path + '/AA/' + name + '_AA.txt'), 'w') as AA_file:
                AA_file.write('>' + name + '\n')
                AA_file.write(gene.AA + '\n')
                AA_file.close()

            #NEED TO ASK ABOUT THIS...IS IT 100 AA OR 100 NUCLEOTIDES
            #Check to make sure all AA sequences are at least 100 aa


        gene_dic[name] = gene

    #call clustal omega
    in_file = path + "/ALL.txt"
    out_file= path + "/aligned.fasta"
    clustalo = ClustalOmegaCommandline(infile=in_file, outfile=out_file, auto=True, cmd="/nfshome/hcarroll/public_html/apps/clustalOmega/bin/clustalo")
    clustalo()


    ##remove this line and uncomment above increment after done with debugging
    #        row += 1

    clustalDic = {}
    tempDic = gene_dic
    count = 0
    frame_count = {}
    with open(out_file, 'r') as clustal:
            for line in clustal:
                line = line.rstrip()
                if(line[0] is '>'): # and count is 0):
                    name = line.lstrip('>')
                else:
                    gene_dic[name].Aligned_str += line

    for name, gene in gene_dic.items():
        frame_seq = ''
        count = 0
        CDS_length = 0
        exon = 0
        start = False
        start_pos = 0
        frame_len = int(gene.CDS_len[exon]/3)
        while frame_len == 0:
            gene.Frames[exon].append(-2)
            exon += 1
            frame_len = int(gene.CDS_len[exon]/3)
        for char in gene.Aligned_str:
            if char is '-':
                frame_seq += '.'
            else:
                if start is False:
                    frame_seq += gene.Frames[exon][1]
                    start = True
                    start_pos = count
                    gene.Frames[exon].append(start_pos)
                elif frame_len is 0:
                    exon += 1
                    frame_seq += gene.Frames[exon][1]
                    start_pos = count
                    gene.Frames[exon].append(start_pos)
                    frame_len = int(gene.CDS_len[exon]/3)
                    while frame_len == 0:
                        gene.Frames[exon].append(-2)
                        exon += 1
                        frame_len = int(gene.CDS_len[exon]/3)
                else:
                    frame_seq += '.'
                    frame_len -= 1
            count += 1
            CDS_length += 1
        while exon < len(gene.CDS_len)-1:
            exon += 1
            gene.Frames[exon].append(-2)

    #Make Final 2D array
    CombinedLists = []
    finalList = []
    rowCount = 0
    for name, gene in gene_dic.items():
        CombinedLists.append([name])
        for col in gene.Frames:
            CombinedLists[rowCount].append(col)
        rowCount += 1

    preFinal = []
    temp = []
    rowCount = 0
    count = 1
    count2 = 0
    #for row in testList:
    for row in CombinedLists:
        name1 = row[0]
        preFinal.append([name1])
        for col in row[1:]:
            rowCount = count
            rowOn = 0
            exonNum = col[0]
            exonSize = col[4]
            MSALoc = col[5]
            if col[5] is -2:
                preFinal[count2].append([col[5], col[4]])
            elif count is 1:
                preFinal[count2].append([col[5], col[4]])
                for Nrow in CombinedLists[rowCount:]:
                    name2 = Nrow[0]
                    for Ncol in Nrow[1:]:
                        if col[1] == Ncol[1] and col[2] == Ncol[2]:
                            if (abs(col[4] - Ncol[4]) % 3) is 0:
                                if abs(col[5] - Ncol[5]) <= 10:
                                    preFinal[count2][col[0]+1].append([name2, Ncol[5], Ncol[4], Ncol[0]])
            else:
                found = False
                temp = preFinal
                for row1 in temp:
                    if found is False:
                        lastName = row1[0]
                        colOn = 0
                        for col1 in row1[1:]:
                            if len(col1) > 2:
                                count3 = 0
                                while count3 < len(col1)-2 and found is False:
                                    if name1 == col1[count3+2][0] and exonNum == col1[count3+2][3]:
                                        found = True
                                        preFinal[count2].append([-1, [rowOn, colOn, col1[0]]])
                                    count3 += 1
                            colOn += 1
                    rowOn += 1

                if found is False:
                    preFinal[count2].append([col[5], exonSize])
                    for Nrow in CombinedLists[rowCount:]:
                        name2 = Nrow[0]
                        for Ncol in Nrow[1:]:
                            if col[1] == Ncol[1] and col[2] == Ncol[2]:
                                if (abs(col[4] - Ncol[4]) % 3) is 0:
                                    if abs(col[5] - Ncol[5]) <= 10:
                                        preFinal[count2][col[0]+1].append([name2, Ncol[5], Ncol[4], Ncol[0]])
            rowCount += 1
        count +=1
        count2 +=1

    print("preFinal")
    for row in preFinal:
        print(len(row))
        print(row)
    print()

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
            if col[0] == -2:
                startNum += 1
        if startNum > maxStart:
            maxStart = startNum
    print('max start is', maxStart)


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
            print('xStart is ', xStart)
            name = row[0]
            colCount = 0
            color = RGB[rowCount]
            printSVG.write('\n<!--'+name+'-->\n')
            print(name)
            prevNumOfExons = -1
            prevExonCount = 0
            prevDy = -1
            for col in row[1:]:
                count = 0
                if colCount is 0:
                    if col[0] >= -1:
                        x = xStart + (maxStart * 75)
                        firstX = x + 50
                        print('firstx is', firstX)
                    else:
                        firstX = xStart + (maxStart * 75) + 50
                elif colCount > 0:
                    if col[0] > -1:
                        x = 3*col[0]
                        print('col[0] is', col[0], x, firstX)
                        if col[0] < firstX:
                            x = firstX + 60
                        if x < prevX+50:
                            x = prevX + 100
                    elif  col[0] is -1:
                        x = 3*col[1][2]
                        print('col[1][2]', col[1][2], x, firstX)
                        if col[1][2] < firstX:
                            x = firstX + 60
                        if x < prevX+50:
                            x = prevX + 100
                    elif col[0] is -2:
                        x = prevX + 60
                    printSVG.write('<text x=\"'+str(25)+'\" y=\"'+str((y+(y+50))/2)+'\" style=\"stroke:'+color+'\">\"'+name+'\"</text>')
     #               printSVG.write('<line x1=\"'+str((prevX+50))+'\" y1=\"'+str((yLineStart))+'\"x2=\"'+str((x))+'\" y2=\"'+str((yLineStop))+'\" style=\"stroke:'+color+'; stroke-width:2\"/>\n')
                numOfExons = len(col)-1
                dy = 50 + (50 * 0.25 * (numOfExons-1))
                if col[0] == -2:
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
    #os.system("chmod -R 777 files/*")
    os.system("echo 'DEBUGGING:' > /tmp/hcarroll.tmp; chmod 777 /tmp/hcarroll.tmp")
    fileItems = form['filename[]']
    
    message = ""

    mkDir()

    for fileItem in fileItems:
        if fileItem.file:
            fn = os.path.basename(fileItem.filename.replace("\\", "/"))
            try:
                open(path + '/' + fn, 'wb').write(fileItem.file.read())
            except:
                print("""Content-Type: text/html\n\n
                        <html>
                        <body>
                            <p>%s</p>
                        </body>
                        </html>
                        """ %(path+'/'+fn))
                exit()
#            message = message + 'The file "' + fn + '"was uploaded successfully with the path' + path + '\n'
#        else:
#            message = 'No file was uploaded'

    dataProcess()

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
    os.system("echo 'DEBUGGING: redirectStr: " + redirectStr + "' > /tmp/hcarroll.tmp; chmod 777 /tmp/hcarroll.tmp")
    print(redirectStr)

if __name__ == "__main__":
    main() 
    exit(0)
