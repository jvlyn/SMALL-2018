from copy import deepcopy

class Graph:
    """Graph module. Internal representation consists of a list
    of integers in ascending order (vertices) and an adjacency
    dictionary which keeps track of edges of the graph.
    The adjacency dictionary consists of vertex, vertex list
    pairs where the vertex list can optionally contain tuples
    if a given vertex has multiple edges connecting it.
    """
    def __init__(self, vert=None, edges=None):
        """Constructs a Graph object from a given vertex list vert
        and an adjacency dictionary.
        """
        if (vert == None or edges == None): # empty graph
            self.vertices = []
            self.edges = {}
        else:
            self.vertices = vert
            newDict = {}
            for v in edges:
                newList = []
                for elmt in edges[v]:
                    if isinstance(elmt, tuple):
                        newList += [elmt]
                    else:
                        newList += [(elmt, 1)]
                newDict[v] = newList
            self.edges = newDict

    def getNumVertices(self):
        """Returns number of vertices in the graph.
        """
        return len(self.vertices)

    def getNumEdges(self):
        """Returns number of edges in the graph.
        """
        num = 0
        for v in self.edges:
            num += len(self.edges[v])
        return num / 2 # each edges is counted exactly twice

    def connected(self, v1, v2):
        """Given input vertices v1, v2, determines if v1 is connected
        to v2 by at least one edge.
        """
        return v2 in map(lambda x : x[0], self.edges[v1]) or \
                v1 in map(lambda x : x[0], self.edges[v2])

    def vertexIn(self, v):
        """Determine whether the vertex v is in the graph.
        """
        return v in self.vertices

    def Laplacian(self):
        """Generate the Laplacian of the graph as a matrix object over QQ.
        """
        diagonalMatrix = [[0]*len(self.vertices) for i in range(len(self.vertices))]
        for v in range(len(self.vertices)):
            diagonalMatrix[v][v] = sum(map(lambda x : x[1], self.edges[v+1]))
        adjacencyMatrix = []
        for i in range(len(self.vertices)):
            row = []
            for j in range(len(self.vertices)):
                if self.connected(i+1, j+1):
                    for adjVertex, num in self.edges[i+1]:
                        if adjVertex == j+1:
                            row += [num]
                else:
                    row += [0]
            adjacencyMatrix += [row]
        laplacian = []
        for i in range(len(diagonalMatrix)):
            row = []
            for j in range(len(diagonalMatrix)):
                row += [diagonalMatrix[i][j] - adjacencyMatrix[i][j]]
            laplacian += [row]
        return matrix(QQ, laplacian) # Laplacian over rationals

    def genus(self):
        """Determine the genus of the graph where
        genus = |edges| - |vertices| + 1.
        """
        return self.getNumEdges() - self.getNumVertices() + 1

    def __eq__(self, g):
        """Overloaded method to determine if two graph objects are equal.
        We define two graphs to be equivalent if they have identical
        vertex and edge labels. This does NOT check for equivalence
        up to isomorphism.
        """
        if set(self.vertices) == set(g.vertices):
            for vertex in self.vertices:
                if set(self.edges[vertex]) != set(g.edges[vertex]):
                    return False
            return True
        return False

    def __ne__(self, g):
        """Overloaded method to determine if two graph objects are not equal.
        """
        if set(self.vertices) == set(g.vertices):
            for vertex in self.vertices:
                if set(self.edges[vertex]) != set(g.edges[vertex]):
                    return True
            return False
        return True
