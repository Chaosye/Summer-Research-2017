### hex_to_rgb and rgb_to_hex functions are from https://stackoverflow.com/questions/214359/converting-hex-color-to-rgb-and-vice-versa

def hex_to_rgb(value):
    """Return (red, green, blue) for the color given as #rrggbb."""
    value = value.lstrip('#')
    lv = len(value)
    return tuple(int(value[i:i + lv // 3], 16) for i in range(0, lv, lv // 3))

def rgb_to_hex(red, green, blue):
    """Return color as #rrggbb for the given color values."""
    return '#%02x%02x%02x' % (red, green, blue)

###

### Utilities used for Firehose data processing.

def UniProtID_Common_Dictionary():
    #Input: Reads from "human-gene-map.txt"
    #Process: Processes the text file, pulling the first and last columns, and making entries of each other.
    #Ouput: Dictionary with keys of both UniProtID and Common names for genes, with entries of each other.

    namedict = {}

    namedatabase = fileReader("human-gene-map.txt", True)
    for names in range(len(namedatabase)):
        namedict[namedatabase[names][0]] = namedatabase[names][len(namedatabase[names])-1]
        namedict[namedatabase[names][len(namedatabase[names])-1]] = namedatabase[names][0]
    #print(namedict)
    return namedict

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
    
def fileWriter(fileName, textOutput):
    #Input: Takes a text ouput from something, and then writes it directly into a file.
    target = open(fileName, 'w')
    target.truncate()
    target.write(textOutput)
    target.close()
    return

def normalizeValues(valRange, valArray):
    #Input: Reads in a desired range for values to placed in, and a single layer array.
    #Process: First performs the function where x is an array of values:
    # f(x) = (x - min(x))/(max(x) - min(x))  #Which normalizes to a range of [0,1]
    # If a different range is desired, then the function reads valRange as [x,y]:
    # where range = y - x
    # f(a) = (a * range) + x  #Which normalizes to the range [x,y]
    #Ouput: Dictionary with keys of both UniProtID and Common names for genes, with entries of each other.
    normalizedArray = []
    tempArray = []
    maxDifference = max(valArray) - min(valArray)
    for x in range(len(valArray)):
        normalizedArray.append(1-((valArray[x] - min(valArray))/maxDifference))
    if valRange != [0,1]:
        print("Desired lineweight range is not [0,1]. It is ", valRange)
        tempArray = normalizedArray
        normalizedArray = []
        normalDifference = max(valRange) - min(valRange)

        for a in range(len(valArray)):
            normalizedArray.append((tempArray[a]*normalDifference) + min(valRange))
    print("The final normalizedArray is: ", normalizedArray) #DB
    return normalizedArray

###
