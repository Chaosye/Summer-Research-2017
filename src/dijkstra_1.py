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
	dijkstra(staticTestGraph())
	
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
	'g': ['h', 1],
	'h': ['i', 1],
	'i': []
	}
	
	return adjacencyDict
	
def dynamicTestGraph(nodeNum, edgeNum):
			
	
	return

def dijkstra(Graph):
	
	sptSet = {} #Shortest path tree set - the possible vertices in the shortest path tree
	spt = set()  #Shortest path tree - the tree being constructed.
	dist = {} #Distance from source to a vertex
	prev = {} #Holder for previous node in the optimal path from source to target
	source = 'a' #Initializing the source as 'a'
	pathRunner = True
	currentnode = source
	destinationdistHolder = []
	
	#Initialization
	sptSet = sorted(Graph.keys())
	for n in range(len(sptSet)):
		dist[sptSet[n-1]] = math.inf  #Unknown nodes are infinitely far away
		prev[sptSet[n-1]] = None      #Previous node is unknown
	
	dist['a'] = 0  #Setting source distance to 0
	
	
	#Creating path length
	destination = [] #Holder for the possible destinations of a path from the current node
	shortdist = [] #Holder for the lengths of the possible destinations
	sourceUsed = False #Stating for the coming loop that the source has not yet been considered
	nextup = 'a'
	
	for n in range(len(sptSet)): #For every unfilled node in sptSet...
		if sourceUsed == False:  #If source hasn't been considered - must be to start off the graph
			spt.add(source)
			for n in range(int(len(Graph[source])/2)):  #For all destinations/distances for a particular node...
				destination.append((Graph[source][(n-1)*2])) #Destination pulls out the nodes that will be traveled to (e.g. for 'a' it is ['b','c']
				dist[destination[n]] = (Graph[source][(n*2)-1])
				spt.add(destination[n])
				prev[destination[n]] = 'a'
					
			sourceUsed = True #Source node and connected nodes have been added to spt.
		else:  ##The main section of the loop for all nodes beyond the source
			'''
			Starting this section with c = 5 and b = 1, with all three included in spt
			build a set of all possible destinations from both b and c
			pick the smallest distance destination of all of them, remove it from the consider set, and add that node to spt
			then add all of its destinations and repeat
			finish when there are no more destinations left
			'''
			for n in range(len(destination)):
				if dist[destination[n]] < math.inf:
					destinationdistHolder.append(dist[destination[n]])
					print("The distances of our possible destinations ", destination," are ",destinationdistHolder)  #DB
			currentnode = destination[destinationdistHolder.index(min(destinationdistHolder))]
			print("The current node being examined is: ", currentnode)   #DB
			destinationdistHolder.pop(destination.index(currentnode))
			destination.remove(currentnode) #Removes the node from nodes to be considered
			
			for n in range(int(len(Graph[currentnode])/2)):  #For all destinations/distances for a particular node...
				destination.append((Graph[currentnode][(n-1)*2])) #Destination pulls out the nodes that will be traveled to (e.g. for 'a' it is 'b' then 'c')
				tempdist = (Graph[currentnode][(n*2)-1]) + dist[currentnode]
				if tempdist <= dist[destination[n]]:
					dist[destination[n]] = tempdist
				prev[destination[n]] = currentnode
				
		if currentnode == "i":
			break
			
	
	print("dist: ", dist)  #DB
	print("prev: ", prev)  #DB
	
	printVar = "i"
	while pathRunner == True:
		if printVar == "a":
			pathRunner = False
		print(printVar)
		printVar = prev[printVar]
	
	return

if __name__ == '__main__':
	main()