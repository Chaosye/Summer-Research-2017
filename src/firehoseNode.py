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
print(sys.path)
import itertools
from collections import defaultdict
from NetPathProcessor import fileReader, NetPath_pathwayGenerator, NetPath_nodeGenerator, NetPath_GraphSpace
from utils import rgb_to_hex, UniProtID_Common_Dictionary
#Pulling in packages
#sys.path.insert(0, '../packages')
print(sys.path)
from scipy import stats
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
        filePathways = NetPath_pathwayGenerator(fileText, "NetPath")

        #print(fileNodes)
        nodedata.update(fbget_PerPathwayInfo(fileNodes, True)) #Produces the dictionary holding unprocessed and processed data for each node.

        NetPath_GraphSpace(filePathways, fileNodes, "(Node Version) "+ (pathways[paths] + "Pathway (NetPath with fbget data and ttest)"), pathways[paths] + " pathway, uploaded to GraphSpace and generated using NetPath data.", "Common", nodedata, namedict, "T-test")
    #print(nodedata["BTRCnormalizedavg"]) #DB

    ##Other options for fbget analysis below
    #fbget.maf(gene="wnt3", cohort="coad") #Collapsible nodes? popup panel that links to other graphs.  Note that exporting graphs of laid out nodes saves location of nodes, and so you can copy graph layouts between each other.
    #fbget.smg(gene="wnt3", cohort="coad")
    #fbget.cn_genes_all()
    #GraphSpace Uploading

    return

def fbget_PerPathwayInfo(pathwayNodes, pairedcases):
    #Input: Nodes from NetPath for a particular pathway
    #Process: Produces fbget data for every single node in the pathway, puts the data into a dictionary for that node that is formatted as "[Node][output]" in a string.
    #Output: An array>array containing the information in the text files separated by column and line.

    #Initialize variables
    pathwaydict = {} #Holds all the fbget information per node in a dictionary.  Processed information is accessible by "[Node][output]" - which is the node, and the type of information in one string.
    #Variables initialized.

    for node in range(len(pathwayNodes)):
        #For determining what kind of data you'd like to pull from firebrowse - for gene expression data, 2 means expression_log2, while 3 means the z-score of that data point.
        print(fbget.mrnaseq(gene=pathwayNodes[node][2], cohort="coad"))

        print("The node currently being analyzed is: "+ pathwayNodes[node][2])
        #Pulls the mRNA sequencing activity
        fbget_Output = fbget_Reader(fbget.mrnaseq(gene=pathwayNodes[node][2], cohort="coad", page_size = 2000), {2})
        #Pulls out patient ID and sample type.
        fbget_Patient = fbget_Reader(fbget.mrnaseq(gene=pathwayNodes[node][2], cohort="coad", page_size = 2000),{0,5})

        #Production of the information dictionary pathwaydict
        pathwaydict[pathwayNodes[node][2]] = fbget_Output
        pathwaydict[(pathwayNodes[node][2]) + "ID"] = fbget_Patient
        if pairedcases == True:
            Ndict, Tdict = fbget_PatientPairer(pathwaydict, ((pathwayNodes[node][2])))

            Narray = list(itertools.chain.from_iterable(Ndict.values()))
            Tarray = list(itertools.chain.from_iterable(Tdict.values()))
            #Performing a t-test, and then other statistical analysis.
	    print("The amount of normal samples is: ", len(Narray))
	    print("The amount of tumor samples is: ", len(Tarray))
            tstat, tsignificance =  stats.ttest_ind(Narray, Tarray, equal_var = True)

            if tsignificance < 0.05:
                print("The difference is significant for ", pathwayNodes[node][2])

            pathwaydict[(pathwayNodes[node][2]) + "tstat"] = tstat
            pathwaydict[(pathwayNodes[node][2]) + "tsignificance"] = tsignificance
            print("The t-test results for ", pathwayNodes[node][2] ,"are ", pathwaydict[(pathwayNodes[node][2]) + "tstat"], " and ", pathwaydict[(pathwayNodes[node][2]) + "tsignificance"])
            print("Normalizing Narray: ")
            fbget_Normalizing(Narray, False)
            print("Normalizing Tarray: ")
            fbget_Normalizing(Tarray, False)

        #Nodes in a pathway, and then eventually nodes for an entire type of cancer on the pathways it has been seen to affect.

        normalizedvalue, avg , normalizedavg, zscores = fbget_Normalizing(fbget_Output, True)

        pathwaydict[(pathwayNodes[node][2]) + "normalizedvalue"] = normalizedvalue
        pathwaydict[(pathwayNodes[node][2]) + "avg"] = avg
        pathwaydict[(pathwayNodes[node][2]) + "normalizedavg"] = normalizedavg
        pathwaydict[(pathwayNodes[node][2]) + "zscores"] = zscores
        if tsignificance < 0.01:
            pathwaydict[(pathwayNodes[node][2]) + "huevalue"] = rgb_to_hex(255, int(255*tsignificance), int(255*tsignificance))
        else:
            pathwaydict[(pathwayNodes[node][2]) + "huevalue"] = rgb_to_hex(255, int(255), int(255))
        #e.g. To access WNT7's normalizedavg, use the value "pathwaydict["WNT7normalizedavg"]""

    #print(pathwaydict["WNT7Anormalizedavg"])  #DB

    return pathwaydict

def fbget_PatientPairer(pathwaydict, nodename):
    #Input: An array with the format [patientID, sample_type]
    #Process: Isolate the patients that reappear within the patient set with different sample types.
    #Output: Return the indices of these patients' data, with each patient as an individual array.

    #Initialize variables:
    patientID = []
    sample_type = []
    IDindices = defaultdict(list)
    locationdict = {}
    Ndict = defaultdict(list)
    Tdict = defaultdict(list)


    #Variables initialized.
    #Ports info from fbget_Patient to two flat arrays
    for row in range(len(pathwaydict[nodename])):
        patientID.append(pathwaydict[nodename + "ID"][row][0])
        sample_type.append(pathwaydict[nodename + "ID"][row][1])

    #Creates dictionary of appearing indices, and removes any entries that don't have duplicates.
    for index,ID in enumerate(patientID):
        #print(IDindices[ID])  #DB
        #print(index)  #DB
        IDindices[ID].append(index)
    for key,indices in IDindices.items():
        if len(indices) < 2:
            del IDindices[key]

    #Deletes any indices that don't have both a N and T component
    tempIDList = list(IDindices)

    for patient in tempIDList:
        print("Scanning for patient ", patient)
        #Keeps track of how many times a patient reappears, with what sampletype.
        Ncounter = 0
        Tcounter = 0
        for entries in IDindices[patient]:
            print("Line entry number is ", entries)  #DB
            patientsample = sample_type[entries]
            if pathwaydict[nodename][entries][0] == "None":
                pathwaydict[nodename][entries][0] = None
                continue
            if patientsample[0] == "N": #If it has a Normal tag in the TCGA
                Ncounter += 1
                print("Ncounter =", Ncounter)  #DB
                Ndict[patient].append(float(pathwaydict[nodename][entries][0])) #0 is to prevent unnecessary nesting
                #print(Ndict)
            elif patientsample[0] == "T": #If it has a Tumor tag in the TCGA
                Tcounter += 1
                print("Tcounter =", Tcounter)  #DB
                Tdict[patient].append(float(pathwaydict[nodename][entries][0]))
                #print(Tdict)
            else:
                print("Patient sample is: ", patientsample)  #Error
                raise NameError("Unconsidered datatype present for patient.")
                quit()
        if Ncounter == 0 or Tcounter == 0: #requires that there be at least 1 N and T sample, otherwise comparison can't be made.
            print("Deleting ", patient)
            del IDindices[patient]
    #print("The IDindices dict is: ", sorted(IDindices)) #DB

    print("The set of values produced by Patient Pairer are:")  #DB
    #print("patientID: ", patientID)  #DB
    #print("sample_type: ", sample_type)
    #print("IDindices: ", IDindices)
    #print("The dictionary of N-tagged values is: ", Ndict)
    #print("The dictionary of T-tagged values is: ", Tdict)
    #Whats best to do here statistically analyze in a separate method.


    return Ndict, Tdict

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

def fbget_Normalizing(dataToBeNormalized, directOutput):
    #Input: An array, containing the info of a text file separated by col internally and then line from fbget file.  directOutput is a bool asking if the data source is from fbget directly, as it may be pre-processed in some other way.
    #Process: Turns strings into floats.
    #Ouput:

    #Initialize variables:
    valuesToBeNormalized = [] #Single array of values to be normalized in a range
    normalizedValues = [] #Array of values after being normalized in a range
    smallestValue = 0.0
    biggestValue = 0.0
    #Variables initialized.

    #The loop below pulls the values out of a row & col format and puts them into a single layer array.
    if directOutput == True:
        for row in range(len(dataToBeNormalized)):
            for col in range(len(dataToBeNormalized[0])):
                if dataToBeNormalized[row][col] == "None" or dataToBeNormalized[row][col] == None:
                    dataToBeNormalized[row][col] = None
                    continue
                dataToBeNormalized[row][col] = float(dataToBeNormalized[row][col])
                valuesToBeNormalized.append(dataToBeNormalized[row][col])
    elif directOutput == False:
        for item in dataToBeNormalized:
            #print("Item is ", item)  #DB
            valuesToBeNormalized.append(item)
    #Overall normalizing of values below from [0,1]
    smallestValue = min(valuesToBeNormalized)
    biggestValue = max(valuesToBeNormalized)
    difference = biggestValue - smallestValue
    avg = sum(valuesToBeNormalized)/len(valuesToBeNormalized)
    zscores = stats.zscore(valuesToBeNormalized)
    print("big = ", biggestValue, " small = ", smallestValue, " diff = ", difference, "avg = ", avg) #DB
    #Function where (max-min)/(max-min) = 1, (min-min) = 0, and so (x-min)/(max-min)is less than 1.  Therefore normalizing to a range of [0,1]
    for x in range(len(valuesToBeNormalized)):
        normalizedValues.append((valuesToBeNormalized[x]-smallestValue)/(biggestValue - smallestValue))
    normalizedavg = sum(normalizedValues)/len(normalizedValues)
    #Done!
    #print(normalizedValues) #DB

    return normalizedValues, avg, normalizedavg, zscores

if __name__ == '__main__':
    main()
