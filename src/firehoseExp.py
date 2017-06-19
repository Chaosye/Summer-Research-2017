## Firehose Experimentation

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
#Write a program to:
##(1) get ONE TYPE of data from ONE TYPE of cancer from the API
##(2) "map" the data to the NetPath nodes (NP nodes are in the UniProtID)
##(3) get the data for every NODE (scaled/normalized to [0,1] (p-value?))
##(4) annotate graph with this info.
##(5) post to GS

##Start all functions with:
#Input:
#Processes:
#Output:
#Initialize variables:
#Variables initialized.

#Pulls in firehose packages
import sys
sys.path.insert(0, '/afs/reed.edu/user/n/i/nicegan/Desktop/Summer-Research-2017/packages')
import firebrowse
from firebrowse import fbget

def main():
    print (fbget.mrnaseq(gene="wnt3", cohort="coad"))
    return



if __name__ == '__main__':
    main()
