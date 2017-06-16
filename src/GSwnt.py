## Wnt Pathway through Graphspace

## 6/15/17
## Uploads the Wnt pathway to Graphspace

#Citations:

#Notes:
#DB = Debugging

#####################################
## Nick Egan
## Advisor: Anna Ritz
#####################################

# Plan of Attack!  (\\ means complete)
# Figure out how to handle nu Graphspace
## Set up for three shapes (receptor, tf, none)
## Set up for colored/individually directed connections (physical v phosphorylation)
## Set up for multiple possible connections to the same node
# Read in Wnt-edges.txt and Wnt-nodes.txt

##Start all functions with:
#Input:
#Processes:
#Output:
#Initialize variables:
#Variables initialized.

def main():

    fileText = fileReader("Wnt-edges.txt", True)
    filePathways = NetPath_pathwayGenerator(fileText)

    return


def fileReader(fileName, separatecols):
    #Input: A text file, the lines to be read(set), the columns to be read(set), and if columns will be separated into different entries into the final array (bool)
    #Process: Processes the lines and columns into
    #Output: An array containing the information in the text files separated by column and line.

    #Initialize variables:
    fileContents = [] # The final array to contain the file information
    #Variables initialized.

    myFile = open(fileName,'r') ## open the file
    lineString = myFile.read() ## Read the file into one long string
    myFile.close() ## close the file
    lineString = lineString.strip() ## remove extra whitespace
    lineList = lineString.split('\n') ## Split the string by line
    if separatecols == True:
        for pattern in range(len(lineList)):
            lineList[pattern] = lineList[pattern].split('\t')  ## Split the list by col

    #print(lineList) #DB

    return lineList

def NetPath_pathwayGenerator(fileText):
    #Input: An array, containing the info of a text file separated by col internally and then line from a NetPath-edges.txt file.
    #Process: Takes the node and edge information of this text file, in the format of Wnt-edges.txt, so it will take the columns 0, 1, 5
    #Ouput: A smaller array, this time with only useful node and edge information

    #Initialize variables:
    colPathway = [] # Columns of the eventual pathway in a single line
    genPathway = [] # Final generated pathway in format of [input, output, edge type]
    #Variables initialized.
    for row in range(len(fileText)):
        for col in range(len(fileText)):
            if col in [0, 1, 5]:
                colPathway.append(fileText[row][col])
                print("colPathway is:", colPathway)
        genPathway.append(colPathway)
        colPathway = []

    print("The final pathway from this NetPath edges text files is: ")
    for pattern in genPathway:
        print(pattern)

    return genPathway

def NetPath_GraphSpace():
     # Input: Array of arrays containing concatenated from a NetPath-edges.txt file in the format of [Input, Output, Pathway type]
     # As well as a concatenated NetPath-nodes.txt files with additional information for node type.
     # Process: Create nodes with graphical information using NetPath-nodes.txt, and then connect them with graphically detailed connections from the edges information.

     
    return

if __name__ == '__main__':
    main()
