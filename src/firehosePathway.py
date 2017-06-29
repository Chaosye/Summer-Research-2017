## Firehose Experimentation

## 6/15/17
## Uploads the Wnt pathway to Graphspace with information from fbget

#Citations:
# hex_to_rgb and rgb_to_hex functions are from https://stackoverflow.com/questions/214359/converting-hex-color-to-rgb-and-vice-versa

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
import fnmatch
from utils import rgb_to_hex, UniProtID_Common_Dictionary
from NetPathProcessor import fileReader, NetPath_pathwayGenerator, NetPath_nodeGenerator, NetPath_GraphSpace
#Pulling in packages from the packages directory.
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
        print("filePathways: ", filePathways) #DB
        #print(fileNodes)
        nodedata.update(fbget_PerPathwayInfo(fileNodes)) #Produces the dictionary holding unprocessed and processed data for each node.
        nodedata.update(lineWeighter(filePathways, nodedata, "normalizedavg", namedict))

        NetPath_GraphSpace(filePathways, fileNodes, ("(Full Pathway Version) " + pathways[paths] + "Pathway (NetPath with fbget data)"), "The "+ pathways[paths] + " pathway, uploaded to GraphSpace and generated using NetPath data.", "Common", nodedata, namedict, "FullPathway")
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
    valuesToBeNormalized = []
    #Variables initialized.

    for node in range(len(pathwayNodes)):

        print("The node currently being read from fbget is: "+ pathwayNodes[node][2])
        fbget_Output = fbget_Reader(fbget.mrnaseq(gene=pathwayNodes[node][2], cohort="coad"), {3})
        #Nodes in a pathway, and then eventually nodes for an entire type of cancer on the pathways it has been seen to affect.
        value = fbget_PreNormalizing(fbget_Output)
        pathwaydict[pathwayNodes[node][2]] = fbget_Output
        pathwaydict[(pathwayNodes[node][2]) + "unchangedvalue"] = value

    smallBigValues = fbget_PathwayNormalizing(pathwaydict)
    for node in range(len(pathwayNodes)):
        print("The node currently being normalized is: "+ pathwayNodes[node][2])
        pathwaydict.update(fbget_NodeNormalizing(pathwayNodes[node][2], pathwaydict, smallBigValues))
    pathwaydict.update(hueNormalizer(pathwaydict))
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

def fbget_PreNormalizing(fbget_Output):
    #Input: An array, containing the info of a text file separated by col internally and then line from fbget file.
    #Process: Turns strings into floats.
    #Ouput:

    #Initialize variables:
    valuesToBeNormalized = [] #Single array of values to be normalized in a range
    normalizedValues = [] #Array of values after being normalized in a range
    #Variables initialized.

    #The loop below pulls the values out of a row & col format and puts them into a single layer array.
    for row in range(len(fbget_Output)):
        for col in range(len(fbget_Output[0])):
            fbget_Output[row][col] = float(fbget_Output[row][col])
            valuesToBeNormalized.append(fbget_Output[row][col])
            #print("PreNormalize's valuesToBeNormalized are: ", valuesToBeNormalized)
    #print(valuesToBeNormalized) #DB
    return valuesToBeNormalized

def fbget_PathwayNormalizing(pathwaydict):
    #Input: Every value frome every node as part of a dictionary
    #Process: Takes every single value from each node in an entire pathway by reading in the pathway's nodes from pathwaydict.keys() and finds the min/max of them so that each value can be scaled from the same [0,1] scale once returning to the node.
    #Ouput: Minimum value, maximum value.

    #Initialize variables:
    valuesToBeNormalized = [] #Single array of values to be normalized in a range for all nodes
    smallestValue = 0.0
    biggestValue = 0.0
    #Variables initialized.

    l = pathwaydict.keys()
    keysToBeNormalized = fnmatch.filter(l, "*unchangedvalue")
    print("The keys being looked at are: ", keysToBeNormalized)
    for keys in range(len(keysToBeNormalized)):
        for items in range(len(keysToBeNormalized[0])):
            valuesToBeNormalized.append(pathwaydict[keysToBeNormalized[keys]][items])

    #Overall normalizing of values below from [0,1]
    smallestValue = min(valuesToBeNormalized)
    biggestValue = max(valuesToBeNormalized)
    #print("big = ", biggestValue, " small = ", smallestValue, " diff = ", difference, "avg = ", avg) #DB

    #Function where (max-min)/(max-min) = 1, (min-min) = 0, and so (x-min)/(max-min)is less than 1.  Therefore normalizing to a range of [0,1]
    print("An example of a value to be normalized is: ", valuesToBeNormalized[1])
    print("Smallest value is: ", smallestValue)
    print("Biggest value is: ", biggestValue)

    smallBigValues = [smallestValue, biggestValue]

    print("Done!")

    return smallBigValues

def fbget_NodeNormalizing(nodeName, pathwaydict, smallBigValues):
    #Input:
    #Process:
    #Ouput:

    #Initialize variables:
    valuesToBeNormalized = [] #Single array of values to be normalized in a range for this specific node
    smallestValue = smallBigValues[0] #Smallest value of all nodes in the pathway
    biggestValue = smallBigValues[1] #Biggest value of all nodes in the pathway
    pathwaydict[nodeName + "normalizedvalue"] = [] #Contains the node's normalized values, crunched down to [0, 1]
    pathwaydict[nodeName + "normalizedavg"] = 0 #Creates an average number for values of the node, after being normalized across the full pathway's values.
    #e.g. To access WNT7's normalizedavg, use the value "pathwaydict["WNT7normalizedavg"]""

    #Variables initialized.

    valuesToBeNormalized = (pathwaydict[nodeName + "unchangedvalue"])
    #print("The valuesToBeNormalized are ", valuesToBeNormalized)
    for entry in range(len(valuesToBeNormalized)):
        #Takes for each node, scales down all values to fit within the pathway's min/max values using this function:
        # y= (x - min)/(max-min)
        # since (max-min)/(max-min) = 1, turning the first value into a smaller number pushes the result into a range of [0,1]
        pathwaydict[nodeName + "normalizedvalue"].append((valuesToBeNormalized[entry]-smallestValue)/(biggestValue - smallestValue))
    print("The normalized values are ", pathwaydict[nodeName + "normalizedvalue"])
    pathwaydict[nodeName + "normalizedavg"] = sum(pathwaydict[nodeName + "normalizedvalue"])/len(pathwaydict[nodeName + "normalizedvalue"])
    print("The normalized average is ", pathwaydict[nodeName + "normalizedavg"])

    return pathwaydict

def hueNormalizer(pathwaydict):
    #Input:
    #Processes:
    #Output:
    #Initialize variables:
    #Variables initialized.

    huesToBeNormalized = []
    l = pathwaydict.keys()
    keysToBeNormalized = fnmatch.filter(l, "*normalizedavg")
    print("The keys being looked at are: ", keysToBeNormalized)
    for keys in range(len(keysToBeNormalized)):
        huesToBeNormalized.append(pathwaydict[keysToBeNormalized[keys]])
    print("The huesToBeNormalized are: ", huesToBeNormalized)

    largestavg = max(huesToBeNormalized)
    smallestavg = min(huesToBeNormalized)

    print("The largest hue avg is: ", largestavg, " and the smallest is: ", smallestavg)

    for entry in range(len(huesToBeNormalized)):
        print("The node being looked at is: ", keysToBeNormalized[entry].replace("normalizedavg", "normalizedhueavg"))

        print("The huesToBeNormalizedp[entry] is: ", huesToBeNormalized[entry])

        pathwaydict[keysToBeNormalized[entry].replace("normalizedavg", "normalizedhueavg")] = (huesToBeNormalized[entry]-smallestavg)/(largestavg - smallestavg)

        print("The hue pre-processing values are: ", pathwaydict[keysToBeNormalized[entry].replace("normalizedavg","normalizedhueavg")], "and ints they are: ", int(255-(255*pathwaydict[keysToBeNormalized[entry].replace("normalizedavg","normalizedhueavg")])))

        pathwaydict[keysToBeNormalized[entry].replace("normalizedavg", "huevalue")] = rgb_to_hex(255, int(255-(255*pathwaydict[keysToBeNormalized[entry].replace("normalizedavg","normalizedhueavg")])), int(255-(255*pathwaydict[keysToBeNormalized[entry].replace("normalizedavg","normalizedhueavg")])))

    return pathwaydict

def lineWeighter(filePathways, pathwaydict, measuredVariable, namedict):
    #Input: The pathway dictionary, with all associated data.
    #Processes: Compares the values between a node and every other node, and assigns a value from [0, 100] based on their similarity.  0. is not the same at all, while 100 is is identical.  This information is passed on to a dictionary that is consulted by the GraphSpace portion of the program that uses this info to change line weight.
    #Output: Provides a dictionary full of floats, used to alter the line weight between nodes as assigned in GraphSpace.
    #Initialize variables:
    smallestValue = 0
    biggestValue = 0
    maxDifference = 0
    edgeList = []
    valuesToBeWeighted = []
    l = pathwaydict.keys()
    keysToBeWeighted = fnmatch.filter(l, ("*" + measuredVariable))
    #pathwaydict[node1 + node2 + "lineweight"]
    #Variables initialized.

    print("The keys being looked at for line weighting: ", keysToBeWeighted)
    for keys in range(len(keysToBeWeighted)):
        valuesToBeWeighted.append(pathwaydict[keysToBeWeighted[keys]])
    biggestValue = max(valuesToBeWeighted)
    smallestValue = min(valuesToBeWeighted)
    maxDifference = biggestValue - smallestValue

    for edges in range(len(filePathways)):
        node1 = namedict[filePathways[edges][0]]
        node1val = pathwaydict[node1 + measuredVariable]
        node2 = namedict[filePathways[edges][1]]
        node2val = pathwaydict[node2 + measuredVariable]
        edge = node1 + node2 #String of both nodes concatenated
        lineweight = 0.1/(abs(node1val - node2val)/maxDifference)
        pathwaydict[node1 + node2 + "lineweight"] = lineweight
        print("The lineweight for" + node1 + node2 + "will be", lineweight)

    return pathwaydict

if __name__ == '__main__':
    main()
