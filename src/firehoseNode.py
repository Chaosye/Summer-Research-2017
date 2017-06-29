## Firehose Experimentation

## 6/15/17
## Uploads the Wnt pathway to Graphspace with information from fbget

#Citations:

#Notes:
#DB = Debugging

#####################################
## Nick Egan
## Advisor: Anna Ritz
#####################################

# Plan of Attack!  (\\ means complete)
#Write a program to:
##\\(1) get ONE TYPE of data from ONE TYPE of cancer from the API
## mRNAseq data from COAD
##(2) get the data for every NODE (scaled/normalized to [0,1] (p-value?))
### Normaliz mRNAseq data by their z-scores
##(3) "map" the data to the NetPath nodes (NP nodes are in the UniProtID)
### Use the normalized data to produce a color gradient from white to red of their significance, as well as possibly line thickness and shape size.
##(4) annotate graph with this info.
### Possibly for MAF (Mutation Annotation Format) data,
##(5) post to GS

##Start all functions with:
#Input:
#Processes:
#Output:
#Initialize variables:
#Variables initialized.

#Initializing packages...
import sys
from NetPathProcessor import fileReader, NetPath_pathwayGenerator, NetPath_nodeGenerator, NetPath_GraphSpace
#Pulling in packages
sys.path.insert(0, '../packages')
import firebrowse
from firebrowse import fbget

import graphspace_python
#Initialized packages.

def main():
    #print(fbget.mrnaseq(gene="wnt3", cohort="coad")) #DB
    ##Formatted as...
    ##tcga_participant_barcode, gene, expression_log2, z-score, cohort, sample_type, protocol, geneID
    ###expression_log2: For that person, the expression of the gene in log2 ratio of (cancer:healthy) in the form of 2**.  e.g. if a sample has a score of 3, that means it's 2**3 (8 times) what their healthy tissue is expressing.

    #Initialize variables
    pathways = [] #The pathways that graphs will be made of.
    nodedata = {} #Holds unprocessed and processed dat for nodes of the pathway.
    #Variables initialized

    '''
    ADD IN FUNCTION THAT PULLS FROM FBGET DATA ALL PATHWAYS ASSOCIATED WITH A CANCER, AND THEN CREATES AN ARRAY OF STRINGS FROM THEM
    Temporarily just do an individual pathway, such as Wnt.
    '''
    namedict = UniProtID_Common_Dictionary()
    pathways = ["Wnt"]

    for paths in range(len(pathways)):
        fileText = fileReader((pathways[paths] + "-nodes.txt"), True) #Reads in the -node.txt NetPath file for a particular pathway in pathways.
        fileNodes = NetPath_nodeGenerator(fileText) #Processes fileText into a format of [Uniprot ID, node type, Common name]
        fileText = fileReader((pathways[paths] + "-edges.txt"), True)
        filePathways = NetPath_pathwayGenerator(fileText)

        #print(fileNodes)
        nodedata.update(fbget_PerPathwayInfo(fileNodes)) #Produces the dictionary holding unprocessed and processed data for each node.

        NetPath_GraphSpace(filePathways, fileNodes, "(Node Version) "+ (pathways[paths] + "Pathway (NetPath with fbget data)"), pathways[paths] + " pathway, uploaded to GraphSpace and generated using NetPath data.", "Common", nodedata, namedict, "Node")
    #print(nodedata["BTRCnormalizedavg"]) #DB

    ##Other options for fbget analysis below
    #fbget.maf(gene="wnt3", cohort="coad") #Collapsible nodes? popup panel that links to other graphs.  Note that exporting graphs of laid out nodes saves location of nodes, and so you can copy graph layouts between each other.
    #fbget.smg(gene="wnt3", cohort="coad")
    #fbget.cn_genes_all()
    #GraphSpace Uploading

    return

def fbget_PerPathwayInfo(pathwayNodes):
    #Input: Nodes from NetPath for a particular pathway
    #Process: Produces fbget data for every single node in the pathway, puts the data into a dictionary for that node that is formatted as "[Node][output]" in a string.
    #Output: An array>array containing the information in the text files separated by column and line.

    #Initialize variables
    pathwaydict = {} #Holds all the fbget information per node in a dictionary.  Processed information is accessible by "[Node][output]" - which is the node, and the type of information in one string.
    #Variables initialized.

    for node in range(len(pathwayNodes)):
        print("The node currently being analyzed is: "+ pathwayNodes[node][2])
        fbget_Output = fbget_Reader(fbget.mrnaseq(gene=pathwayNodes[node][2], cohort="coad"), {3})
        #Nodes in a pathway, and then eventually nodes for an entire type of cancer on the pathways it has been seen to affect.
        normalizedvalue, avg , normalizedavg = fbget_Normalizing(fbget_Output)
        pathwaydict[pathwayNodes[node][2]] = fbget_Output
        pathwaydict[(pathwayNodes[node][2]) + "normalizedvalue"] = normalizedvalue
        pathwaydict[(pathwayNodes[node][2]) + "avg"] = avg
        pathwaydict[(pathwayNodes[node][2]) + "normalizedavg"] = normalizedavg
        pathwaydict[(pathwayNodes[node][2]) + "huevalue"] = rgb_to_hex(255, int(255-(255*normalizedavg)), int(255-(255*normalizedavg)))
        #e.g. To access WNT7's normalizedavg, use the value "pathwaydict["WNT7normalizedavg"]""

    #print(pathwaydict["WNT7Anormalizedavg"])  #DB

    return pathwaydict


def fbget_Reader(fbget_Output, chosencolumns):
    #Input: A unicode text input from fbget, the columns to be read(set), and if
    #Process: Processes the lines and columns into an array in array format.  Converts from unicode to string.  Isolates the desired columns.
    #Output: An array>array containing the information in the text files separated by column and line.

    #Initialize variables:
    lineList = [] # The array to contain the file information
    columns = [] # Columns of a row after being selected from lineList
    rows = [] # Rows of the above columns in an array of arrays
    #Variables initialized.

    lineString = fbget_Output.strip() ## remove extra whitespace

    lineList = lineString.split('\n') ## Split the string by line

    for pattern in range(len(lineList)):
        lineList[pattern] = lineList[pattern].split('\t')  ## Split the list by col

    lineList.pop(0) #Removes header

    #print(lineList) #DB
    #print(lineList[0]) #DB

    #The loop below isolates the desired columns, and turns the data in the arrays from unicode into strings.
    for row in range(len(lineList)):
        for col in range(7):
            if col in chosencolumns: # If in the columns containing the info we want...
                columns.append(lineList[row][col].encode('ascii')) # Converts from unicode to an ascii string, and then appends onto the overall array
        rows.append(columns)
        columns = []

    #print("The rows produced are:", rows) #DB
    #print(type(rows[0][0])) #DB for unicode encoding to ascii
    return rows

def fbget_Normalizing(fbget_Output):
    #Input: An array, containing the info of a text file separated by col internally and then line from fbget file.
    #Process: Turns strings into floats.
    #Ouput:

    #Initialize variables:
    valuesToBeNormalized = [] #Single array of values to be normalized in a range
    normalizedValues = [] #Array of values after being normalized in a range
    smallestValue = 0.0
    biggestValue = 0.0
    #Variables initialized.

    #The loop below pulls the values out of a row & col format and puts them into a single layer array.
    for row in range(len(fbget_Output)):
        for col in range(len(fbget_Output[0])):
            fbget_Output[row][col] = float(fbget_Output[row][col])
            valuesToBeNormalized.append(fbget_Output[row][col])
    #Done!
    #print(valuesToBeNormalized) #DB

    #Overall normalizing of values below from [0,1]
    smallestValue = min(valuesToBeNormalized)
    biggestValue = max(valuesToBeNormalized)
    difference = biggestValue - smallestValue
    avg = sum(valuesToBeNormalized)/len(valuesToBeNormalized)
    print("big = ", biggestValue, " small = ", smallestValue, " diff = ", difference, "avg = ", avg) #DB
    #Function where (max-min)/(max-min) = 1, (min-min) = 0, and so (x-min)/(max-min)is less than 1.  Therefore normalizing to a range of [0,1]
    for x in range(len(valuesToBeNormalized)):
        normalizedValues.append((valuesToBeNormalized[x]-smallestValue)/(biggestValue - smallestValue))
    normalizedavg = sum(normalizedValues)/len(normalizedValues)
    #Done!
    #print(normalizedValues) #DB

    return normalizedValues, avg, normalizedavg

def rgb_to_hex(red, green, blue):
    """Return color as #rrggbb for the given color values."""
    return '#%02x%02x%02x' % (red, green, blue)

def UniProtID_Common_Dictionary():
    #Input: Reads from "human-gene-map.txt"
    #Process: Processes the text file, pulling the first and last columns, and making entries of each other.
    #Ouput: Dictionary with keys of both UniProtID and Common names for genes, with entries of each other.

    namedict = {}

    namedatabase = fileReader("human-gene-map.txt", True)
    for names in range(len(namedatabase)):
        namedict[namedatabase[names][0]] = namedatabase[names][len(namedatabase[names])-1]
        namedict[namedatabase[names][len(namedatabase[names])-1]] = namedatabase[names][0]
    print(namedict)
    return namedict

if __name__ == '__main__':
    main()
