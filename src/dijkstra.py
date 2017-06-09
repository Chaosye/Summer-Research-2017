## Dijkstra's Algorithm - an algorithm for finding the shortest paths between nodes in a graph.

## 6/6/17
## For a given source node, find a shortest path between that node and every other node.

# Citations:
# http://www.geeksforgeeks.org/greedy-algorithms-set-6-dijkstras-shortest-path-algorithm/
# https://en.wikipedia.org/wiki/Dijkstra%27s_algorithm

# Notes:
# DB = Debugging lines

################################
## Nick Egan
## Advisor: Anna Ritz
################################

import math

def main():
	
	vertexSet = {}
	dijkstra(staticTestGraph(), 'a')
	
	return

def staticTestGraph(): #Static example graph for Dijkstra's to work through
	
	# Format is ['destination', distance from current node to that destination]
	adjacencyDict = {
	'a': ['b', 1,'c', 5],
	'b': ['f', 10],
	'c': ['d', 1, 'e', 1],
	'd': ['e', 2, 'g', 1],
	'e': ['i', 8],
	'f': ['g', 2],
	'g': ['h', 2]
	,
	'i': []
	}
	
	return adjacencyDict
	
def dijkstra(Graph, source):
	
	sptSet = [] #Shortest path tree set - the possible vertices in the shortest path tree
	spt = set()  #Shortest path tree - the tree being constructed.
	dist = {} #Distance from source to a vertex
	prev = {} #Holder for previous node in the optimal path from source to target
	pathRunner = True
	#Currentnode relevant variables
	currentnode = source
	neighbors = []
	immediatedist = {}
	#

	#Initialization
	sptSet = sorted(Graph.keys())
	print("The nodes in sptSet are: ", sptSet)
	for n in range(len(sptSet)):
		dist[sptSet[n-1]] = math.inf  #Unknown nodes are infinitely far away
		prev[sptSet[n-1]] = None      #Previous node is unknown
	
	dist[source] = 0  #Setting source distance to 0
	
	
	#Creating path length
	destination = {} #Holder for the possible destinations of a path from the current node
	nextup = [] #Holds the next node to come.
	
	while len(sptSet) > 0:
		#Emptying temporary holders
		neighbors = []
		immediatedist = {}
		
		
		sptSet.remove(currentnode)
		print("The nodes in sptSet to be analyzed are: ", sptSet)
		
		#Find the neighbors of currentnode
		for n in range(int(len(Graph[currentnode])/2)):  #For all destinations/distances for a particular node...
			neighbors.append((Graph[currentnode][(n)*2])) #Destination pulls out the nodes that will be traveled to (e.g. for 'a' it is ['b','c']
			immediatedist[neighbors[n]] = (Graph[currentnode][n*2 + 1])
		print("The neighbors of this node are: ", neighbors," and the distances to them are: ", immediatedist)

		#Checking if new distances are shorter than existing ones for nodes.
		for n in range(len(neighbors)):
			if dist[neighbors[n]] > immediatedist[neighbors[n]]:
				dist[neighbors[n]] = immediatedist[neighbors[n]] + dist[currentnode]
			print("The distance of ", neighbors[n], " from ", currentnode, " is ", dist[neighbors[n]])
		
		#Determining next node to be analyzed
		smallestneighbor = min(immediatedist.values())
		print("The shortest distance neighbors from ", currentnode, " is ", smallestneighbor)
		for n in range(len(sptSet)):
			if dist[sptSet[n]] <= smallestneighbor:
				currentnode = sptSet[n]
		#If there are any nodes in the other nodes with lower value neighbors, however, switch to those instead.

		print("The next node to be analyzed is ", currentnode)
		
		
		#Next iteration of loop
	return
	
if __name__ == '__main__':
	main()