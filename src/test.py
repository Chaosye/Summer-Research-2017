import sys
sys.path.insert(0, '../packages')
import firebrowse
from firebrowse import fbget

print(fbget.mrnaseq(gene="WNT3", cohort="coad"))
