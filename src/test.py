import sys
import fnmatch
from utils import rgb_to_hex, UniProtID_Common_Dictionary, normalizeValues
from NetPathProcessor import fileReader, NetPath_pathwayGenerator, NetPath_nodeGenerator, NetPath_GraphSpace
#Pulling in packages from the packages directory.
sys.path.insert(0, '../packages')
import firebrowse
from firebrowse import fbget

print(fbget.mrnaseq(gene="wnt3", cohort="coad"))
