#!/usr/bin/python3.4
import json


global overArchingArray
overArchingArray = {}
exons = []
genes = []
overArchingArray['Exons'] = exons
overArchingArray['GeneFamilies'] = genes



def addExon(x, y, size=1, scale=1, xscale=1):
    exon = {
    'xpos' : x,
    'ypos' : y,
    'size' : size,
    'scale' : scale,
    'xscale' : xscale
    }
    exons.insert(-1, exon)
def addGeneFamily(listOfExons, color, geneFamilyName):
    geneFamily = {
    'exonsInFamily' : listOfExons,
    'color' : color
    }
    genes.insert(-1, geneFamily)

addExon(1, 1)
addExon(1, 4)
addExon(1, 5)
addExon(1, 6)
addExon(2, 1)
addExon(2, 2)
addExon(2, 3)
addExon(2, 4)
addExon(2, 5)
addExon(2, 6)
addExon(2, 7)
addExon(4, 2, 1.3)
addExon(4, 4, 1.7)
addExon(4, 5.5, 1.9, xscale=3)
addExon(6, 2, 1.9)
addExon(6, 4, 2)
addExon(8, 4, 3)
addExon(10, 4, 3)
addExon(12, 4, 3)
addExon(14, 4, 3)
addExon(17, 3, 2)
addExon(18, 5.5, 2.5, 2)
addExon(19, 3, 2)
addExon(20, 7, 2)
addExon(22, 4, 3)
addExon(25, 4, 3)

addGeneFamily([0, 4, 11, 14, 16, 17, 18, 19, 20, 22, 24, 25], 'blue', 1)
addGeneFamily([5, 11, 15, 16, 17, 18, 19, 20, 22, 24, 25], 'red', 2)
#addGeneFamily([6, 12, 15, 16, 17, 18, 19, 20, 22, 24, 25], 'orange', 3)
#addGeneFamily([1, 7, 12, 15, 16, 17, 18, 19, 21, 24, 25], 'lightblue', 4)
#addGeneFamily([2, 8, 12, 15, 16, 17, 18, 19, 21, 23, 24, 25], 'green', 5)
#addGeneFamily([3, 9, 12, 15, 16, 17, 18, 19, 21, 23, 24, 25], 'brown', 6)
#addGeneFamily([10, 13, 16], 'black', 7)

with open('exonGeneFamilyArray.json', 'w') as test:
    json.dump(overArchingArray,test)
