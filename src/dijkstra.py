## Dijkstra's Algorithm - an algorithm for finding the shortest paths between nodes in a graph.

## 6/6/17
## For a given source node, find a shortest path between that node and every other node.

# Citations:
# http://www.geeksforgeeks.org/greedy-algorithms-set-6-dijkstras-shortest-path-algorithm/
# https://en.wikipedia.org/wiki/Dijkstra%27s_algorithm

################################
## Nick Egan
## Estimate time spent: 
################################

def main():
	
	vertexSet = {}
	dijkstra(staticTestGraph())
	
	return

def staticTestGraph(): #Static example graph for Dijkstra's to work through
	
	adjacencyDict = {
	'a': ['b', 1,'c', 5],
	'b': ['f', 10],
	'c': ['e', 1],
	'd': ['e', 2, 'g', 1],
	'e': ['d', 2, 'i', 8],
	'f': ['d', 1, 'g', 2],
	'g': ['h', 1],
	'h': ['i', 1]
	}
	
	return adjacencyDict
	
def dynamicTestGraph(nodeNum, edgeNum):
			
	
	return

def dijkstra(Graph):
	
	sptSet = {} #Shortest path tree set - the set that will hold vertices in the shortest path tree
	dist = {} #Distance from source to a vertex
	prev = {} #Holder for previous node in the optimal path from source to target
	
	
	
	'''
	for all keys in adjacencyDict, give them all infinite value, and the source a value of 0
	'''
	
	return

if __name__ == '__main__':
	main()