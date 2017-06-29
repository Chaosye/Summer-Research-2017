import json
from graphspace_python.graphs.classes.gsgraph import GSGraph
from graphspace_python.api.client import GraphSpace
from optparse import OptionParser
##
##Explain sif format
##To do:
##-Modify to work with more than two columns
##http://wiki.cytoscape.org/Cytoscape_User_Manual/Network_Formats
##
parser = OptionParser()

parser.add_option('-f', '--File Input', dest='FileInput', action='store', default=None, help='File you want to convert. (Required)')

parser.add_option('-n', '--Graph Name', dest='GraphName', action='store', default='Practice Graph', help='Name for graph to post on Graphspace.')

parser.add_option('-c', '--Graph Description', dest='GraphDescription', action='store', default='No description', help='Description for the graph.')

parser.add_option('-l', '--Graphspace Username', dest='GraphspaceLogin', action='store', default='ramirezg@reed.edu', help='Username for graphspace.')

parser.add_option('-p', '--Graphspace Password', dest='GraphspacePassword', action='store', default='summer2k17', help='Password for graphsapce.')

parser.add_option('-u', '--Upstream Input', dest='upstream', action='store', default=None, help='Upstream Input for TieDIE')

parser.add_option('-d', '--Downstream Input', dest='downstream', action='store', default=None, help='Downstream Input for TieDIE')

parser.add_option('-i', '--Unit ProID', dest='unitproID', action='store', default=None, help='Translate ID to common name for protein')

parser.add_option('-o', '--Other Pathway', dest='OtherPathway', action='store', default=None, help='Pathway you want to compare.')

(opts, args) = parser.parse_args()

##
##Input Validation
##
##Fatal Errors
##

if opts.FileInput is None:
	sys.stderr.write("Warning: Must supply and input file")
	sys.exit(1)
elif opts.GraphName is None:
	sys.stderr.write("Warning: Must supply a Graph name")
	sys.exit(1)


G = GSGraph()

##OutputFile = opts.FileInput

def ConvertToNodes(file): ##Converts the Upstream/Downstream inputs into nodes so they can be compared to the nodes in the pathway to find overlaps.
	nodes = set()

	if file == None:
		return nodes

	with open(file) as t:
		g = [line.split('\t') for line in t.readlines().copy()]
		for x in g:
			nodes.add(x[0])
	return nodes

def FileToDictionary(fileName):

	namedict={}
	if fileName == None:
		return namedict

	with open(fileName) as t:
		g = [line.split('\t') for line in t.readlines().copy()]
		for x in g:
			for y in [0,5]:
				x[y] = x[y].strip('\n')
				
				namedict[x[5]] = x[0]

	return namedict

def get_node_style(node, USnodes, DSnodes,OtherPath):
	if node in USnodes:
		NodeShape = 'triangle'
		NodeColor = 'yellow'

	elif node in DSnodes:
		NodeShape = 'ellipse'
		NodeColor = 'yellow'

	elif node in OtherPath:
		NodeShape = 'star'
		NodeColor = '#F33F6F'

	elif node in OtherPath and node in USnodes:
		NodeShape = 'pentagon'
		NodeColor = '#F33F6F'

	elif node in OtherPath and node in DSnodes:
		NodeShape = 'hexagon'
		NodeColor = '#B200A1'

	else:
		NodeShape = 'rectangle'
		NodeColor = 'red'

	return NodeShape,NodeColor

def Convert(opts):
	nodes=set()
	uni2name=FileToDictionary(opts.unitproID)
	DSnodes=ConvertToNodes(opts.downstream)
	USnodes=ConvertToNodes(opts.upstream)
	OtherPath=ConvertToNodes(opts.OtherPathway)
	with open(opts.FileInput) as t:
		g = [line.split('\t') for line in t.readlines().copy()]
		for x in g:
			for y in [0,2]: ##Only care about columns 1 and 3.
				x[y] = x[y].strip('\n')
				if x[y] not in nodes:
					nodes.add(x[y])
					G.add_node(x[y], label=uni2name.get(x[y], x[y]))
					print('Current Node=', uni2name.get(x[y]))
					NodeShape, NodeColor = get_node_style(x[y], USnodes, DSnodes,OtherPath)
					G.add_node_style(x[y], shape=NodeShape, color=NodeColor, width=90, height=90)

			G.add_edge(x[0], x[2], directed=True)
			G.add_edge_style(x[0], x[2], directed=True, edge_style='dotted')


		G.set_name(opts.GraphName)

		G.set_data(data={'Graph Key':'Yellow Triangle = Upstream Input; Yellow Ellipse = Downstream Input; Red Rectangle = TieDie Result Path'})

		G.set_data(data={'Description':opts.GraphDescription})

		G.set_tags(['No Tags'])
		graphspace = GraphSpace(opts.GraphspaceLogin, opts.GraphspacePassword)
		graphspace.post_graph(G)

Convert(opts)
