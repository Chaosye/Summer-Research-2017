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
### Weight mRNAseq data by their z-scores
##(3) "map" the data to the NetPath nodes (NP nodes are in the UniProtID)
### Use the weighted data to produce a color gradient from white to red of their significance, as well as possibly line thickness and shape size.
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
from GSwnt import fileReader, NetPath_pathwayGenerator, NetPath_nodeGenerator, NetPath_GraphSpace
#Pulling in packages
#sys.path.insert(0, '/afs/reed.edu/user/n/i/nicegan/Desktop/Summer-Research-2017/packages')
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

        NetPath_GraphSpace(filePathways, fileNodes, ("(Full Pathway Version) " + pathways[paths] + "Pathway (NetPath with fbget data)"), "The "+ pathways[paths] + " pathway, uploaded to GraphSpace and generated using NetPath data.", "Common", nodedata, namedict)
    #print(nodedata["BTRCweightedavg"]) #DB

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
    valuesToBeWeighted = []
    #Variables initialized.

    for node in range(len(pathwayNodes)):

        print("The node currently being read from fbget is: "+ pathwayNodes[node][2])
        fbget_Output = fbget_Reader(fbget.mrnaseq(gene=pathwayNodes[node][2], cohort="coad"), {3})
        #Nodes in a pathway, and then eventually nodes for an entire type of cancer on the pathways it has been seen to affect.
        value = fbget_PreWeighting(fbget_Output)
        pathwaydict[pathwayNodes[node][2]] = fbget_Output
        pathwaydict[(pathwayNodes[node][2]) + "unchangedvalue"] = value
        # pathwaydict[(pathwayNodes[node][2]) + "weightedvalue"] = weightedvalue
        # pathwaydict[(pathwayNodes[node][2]) + "avg"] = avg
        # pathwaydict[(pathwayNodes[node][2]) + "weightedavg"] = weightedavg
        # pathwaydict[(pathwayNodes[node][2]) + "huevalue"] = huevalue

        #e.g. To access WNT7's weightedavg, use the value "pathwaydict["WNT7weightedavg"]""
    smallBigValues = fbget_PathwayWeighting(pathwaydict)
    for node in range(len(pathwayNodes)):
        print("The node currently being weighted is: "+ pathwayNodes[node][2])
        pathwaydict.update(fbget_NodeWeighting(pathwayNodes[node][2], pathwaydict, smallBigValues))
    pathwaydict.update(hueWeighter(pathwaydict))
    #print(pathwaydict["WNT7Aweightedavg"])  #DB

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

def fbget_PreWeighting(fbget_Output):
    #Input: An array, containing the info of a text file separated by col internally and then line from fbget file.
    #Process: Turns strings into floats.
    #Ouput:

    #Initialize variables:
    valuesToBeWeighted = [] #Single array of values to be weighted in a range
    weightedValues = [] #Array of values after being weighted in a range
    smallestValue = 0.0
    biggestValue = 0.0
    huevalue = 0
    #Variables initialized.

    #The loop below pulls the values out of a row & col format and puts them into a single layer array.
    for row in range(len(fbget_Output)):
        for col in range(len(fbget_Output[0])):
            fbget_Output[row][col] = float(fbget_Output[row][col])
            valuesToBeWeighted.append(fbget_Output[row][col])
            #print("PreWeight's valuesToBeWeighted are: ", valuesToBeWeighted)
    #Done!
    #print(valuesToBeWeighted) #DB


    return valuesToBeWeighted

def fbget_PathwayWeighting(pathwaydict):
    #Input: Every value frome every node as part of a dictionary
    #Process: Takes every single value from each node in an entire pathway by reading in the pathway's nodes from pathwaydict.keys() and finds the min/max of them so that each value can be scaled from the same [0,1] scale once returning to the node.
    #Ouput: Minimum value, maximum value.

    #Initialize variables:
    valuesToBeWeighted = [] #Single array of values to be weighted in a range for all nodes
    smallestValue = 0.0
    biggestValue = 0.0
    huevalue = 0
    currentnode = None
    #Variables initialized.

    l = pathwaydict.keys()
    keysToBeWeighted = fnmatch.filter(l, "*unchangedvalue")
    print("The keys being looked at are: ", keysToBeWeighted)
    for keys in range(len(keysToBeWeighted)):
        for items in range(len(keysToBeWeighted[0])):
            valuesToBeWeighted.append(pathwaydict[keysToBeWeighted[keys]][items])

    #Overall weighting of values below from [0,1]
    smallestValue = min(valuesToBeWeighted)
    biggestValue = max(valuesToBeWeighted)
    #print("big = ", biggestValue, " small = ", smallestValue, " diff = ", difference, "avg = ", avg) #DB

    #Function where (max-min)/(max-min) = 1, (min-min) = 0, and so (x-min)/(max-min)is less than 1.  Therefore normalizing to a range of [0,1]
    print("An example of a value to be weighted is: ", valuesToBeWeighted[1])
    print("Smallest value is: ", smallestValue)
    print("Biggest value is: ", biggestValue)

    smallBigValues = [smallestValue, biggestValue]

    print("Done!")
    #print(weightedValues) #DB

    return smallBigValues

def fbget_NodeWeighting(nodeName, pathwaydict, smallBigValues):
    #Input: An array, containing the info of a text file separated by col internally and then line from fbget file.
    #Process: Turns strings into floats.
    #Ouput:

    #Initialize variables:
    valuesToBeWeighted = [] #Single array of values to be weighted in a range for this specific node
    smallestValue = smallBigValues[0] #Smallest value of all nodes in the pathway
    biggestValue = smallBigValues[1] #Biggest value of all nodes in the pathway
    huevalue = 0
    pathwaydict[nodeName + "weightedvalue"] = []
    pathwaydict[nodeName + "weightedavg"] = 0

    #Variables initialized.

    valuesToBeWeighted = (pathwaydict[nodeName + "unchangedvalue"])
    #print("The valuesToBeWeighted are ", valuesToBeWeighted)
    for entry in range(len(valuesToBeWeighted)):
        pathwaydict[nodeName + "weightedvalue"].append((valuesToBeWeighted[entry]-smallestValue)/(biggestValue - smallestValue))
    print("The weighted values are ", pathwaydict[nodeName + "weightedvalue"])
    pathwaydict[nodeName + "weightedavg"] = sum(pathwaydict[nodeName + "weightedvalue"])/len(pathwaydict[nodeName + "weightedvalue"])
    print("The weighted average is ", pathwaydict[nodeName + "weightedavg"])



    return pathwaydict

def hex_to_rgb(value):
    """Return (red, green, blue) for the color given as #rrggbb."""
    value = value.lstrip('#')
    lv = len(value)
    return tuple(int(value[i:i + lv // 3], 16) for i in range(0, lv, lv // 3))

def rgb_to_hex(red, green, blue):
    """Return color as #rrggbb for the given color values."""
    return '#%02x%02x%02x' % (red, green, blue)

def hueWeighter(pathwaydict):
    huesToBeWeighted = []
    l = pathwaydict.keys()
    keysToBeWeighted = fnmatch.filter(l, "*weightedavg")
    print("The keys being looked at are: ", keysToBeWeighted)
    for keys in range(len(keysToBeWeighted)):
        huesToBeWeighted.append(pathwaydict[keysToBeWeighted[keys]])
    print("The huesToBeWeighted are: ", huesToBeWeighted)

    largestavg = max(huesToBeWeighted)
    smallestavg = min(huesToBeWeighted)

    print("The largest hue avg is: ", largestavg, " and the smallest is: ", smallestavg)

    for entry in range(len(huesToBeWeighted)):
        print("The node being looked at is: ", keysToBeWeighted[entry].replace("weightedavg", "weightedhueavg"))

        print("The huesToBeWeightedp[entry] is: ", huesToBeWeighted[entry])

        pathwaydict[keysToBeWeighted[entry].replace("weightedavg", "weightedhueavg")] = (huesToBeWeighted[entry]-smallestavg)/(largestavg - smallestavg)

        print("The hue pre-processing values are: ", pathwaydict[keysToBeWeighted[entry].replace("weightedavg","weightedhueavg")], "and ints they are: ", int(255-(255*pathwaydict[keysToBeWeighted[entry].replace("weightedavg","weightedhueavg")])))

        pathwaydict[keysToBeWeighted[entry].replace("weightedavg", "huevalue")] = rgb_to_hex(255, int(255-(255*pathwaydict[keysToBeWeighted[entry].replace("weightedavg","weightedhueavg")])), int(255-(255*pathwaydict[keysToBeWeighted[entry].replace("weightedavg","weightedhueavg")])))

    return pathwaydict

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
