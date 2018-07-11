from itertools import combinations
from itertools import combinations_with_replacement

class Divisor:
    """Divisor class. Internal representation consists of the graph
    that it is built over and a labels dictionary that contains
    vertex, integer pairs.
    """

    def __init__(self, graph, labels):
        """Constructs a Divisor object based upon an input Graph object
        and a dictionary of vertex, integer pairs. Will warn the user if
        dictionary contains vertices not in the graph, but will not warn
        about other possible errors.
        """
        for v in labels:
            if v not in graph.vertices:
                print("Error! Incompatible labelling of vertices!")
                return
        self.graph = graph
        self.labels = labels

    def isEffective(self):
        """Determine if each integer label is nonnegative.
        """
        for n in self.labels.values():
            if n < 0:
                return False
        return True

    def isPrincipal(self):
        """Determine if the divisor is in the image of the graph Laplacian.
        Relies on SAGE's .image() method.
        """
        LP = self.graph.Laplacian()
        divList = []
        for v in sorted(self.labels.keys()):
            divList += [self.labels[v]]
        return vector(QQ, divList) in LP.image()

    def degree(self):
        """Determine the degree (sum of all integer labels) of the divisor.
        """
        total = 0
        for n in self.labels.values():
            total += n
        return total

    def __add__(self, d):
        """Overloaded method to add two divisors by computing the
        pointwise sum of each of their labels.
        """
        if (self.graph != d.graph):
            print("Invalid operands!")
            return
        labels = {}
        for v in self.labels:
            labels[v] = self.labels[v] + d.labels[v]
        return Divisor(self.graph, labels)

    def __sub__(self, d):
        """Overloaded method to subtract two divisors.
        """
        if (self.graph != d.graph):
            print("Invalid operands!")
            return
        labels = []
        for v in self.labels:
            labels += [(v, self.labels[v] - d.labels[v])]
        return Divisor(self.graph, labels)

    def chipFire(self, vertex):
        """Performs the operation of removing val(vertex) chips from the
        given vertex and adding one chip to each of the given vertex's
        adjacent vertices.
        """
        newDict = deepcopy(self.labels)
        d = Divisor(self.graph, newDict)
        for v, n in d.graph.edges[vertex]:
            d.labels[v] += n # number of edges firing
        d.labels[vertex] -= sum(map(lambda x : x[1], self.graph.edges[vertex]))
        return d

    def __checkAllExcept(self, v):
        """Helper method to determine if all vertices except
        the given vertex v have nonnegative values.
        """
        for vertex in self.labels:
            if self.labels[vertex] < 0 and vertex != v:
                return False
        return True

    def allButOne(self, v, num):
        """Helper method to determine if all vertices except
        the given vertex v have values above a certain threshold, num.
        """
        for vertex in self.labels:
            if vertex != v and self.labels[vertex] < 0:
                return False
            elif vertex == v and self.labels[vertex] < num:
                return False
        return True

    def getSemiReducedDivisor(self, v):
        """Helper method to determine the v-semireduced divisor associated
        to the divisor. Uses an algorithm outlined in de Bruyn's thesis.
        """
        tempDict = {}
        dictSum = 0
        for vertex in self.graph.vertices:
            if vertex != v:
                valence = sum(map(lambda x : x[1], self.graph.edges[vertex]))
                value = self.labels[vertex]
                tempDict[vertex] = valence - value
                dictSum += valence - value
        tempDict[v] = -dictSum
        x = Divisor(self.graph, tempDict)
        LP = self.graph.Laplacian()
        tempList = []
        for vertex in sorted(x.labels):
            tempList += [x.labels[vertex]]
        vec = vector(QQ, tempList)
        M = MatrixSpace(QQ, self.graph.getNumVertices(), self.graph.getNumVertices())
        VS = VectorSpace(QQ, self.graph.getNumVertices())
        solution = M(LP).solve_right(VS(vec))
        tempList = []
        for term in solution:
            tempList += [int(term)]
        while tempList[v - 1] != 0:
            if tempList[v - 1] < 0:
                tempList = [val+1 for val in tempList]
            else:
                tempList = [val-1 for val in tempList]
        image = M(LP) * VS(tempList)
        tempDict, i = {}, 1
        for term in image:
            tempDict[i] = term
            i += 1
        newDiv = Divisor(self.graph, tempDict)
        return self + newDiv

    def getReducedDivisor(self, v):
        """Determines the v-reduced divisor associated to the divisor. This
        divisor is shown to exist and be unique in de Bruyn's thesis.
        """
        newDiv = self.getSemiReducedDivisor(v)
        while True:
            verticesOnFire, edgesOnFire = [v], []
            # initially, we want to set just the v vertex on fire
            numOnFire = 1
            for adjVertex, num in newDiv.graph.edges[v]:
                edgesOnFire += [((v, adjVertex), num)] # setting all adjacent edges on fire
            finishedBurning = False
            while not finishedBurning:
                for vertex in newDiv.graph.vertices:
                    numEdges = sum([n for ((v1, v2), n) in edgesOnFire if \
                                v1 == vertex or v2 == vertex])
                    if numEdges > newDiv.labels[vertex] and vertex not in verticesOnFire:
                        # first condition comes from Dhar's burning algorithm
                        # second condition is to avoid adding vertices twice
                        verticesOnFire += [vertex]
                        # if a vertex is added, its adjacent edges also need
                        # to be added
                        for adjVertex, num in newDiv.graph.edges[vertex]:
                            if ((vertex, adjVertex), num) not in edgesOnFire and \
                                ((adjVertex, vertex), num) not in edgesOnFire:
                                edgesOnFire += [((vertex, adjVertex), num)]
                if numOnFire == len(verticesOnFire):
                    # if we haven't set any new vertices on fire compared to
                    # the last time through, then we are done
                    finishedBurning = True
                else:
                    numOnFire = len(verticesOnFire)
            if len(verticesOnFire) == newDiv.graph.getNumVertices():
                # if all vertices have been burned, then we have found the v-red. divisor
                return newDiv
            else:
                # otherwise, we fire from every unburnt vertex
                for vertex in newDiv.graph.vertices:
                    if vertex not in verticesOnFire:
                        newDiv = newDiv.chipFire(vertex)

def gonality(g):
    """Given a graph g, determine its gonality using the methods built
    into the Divisor class.
    """
    k = 1
    while True:
        combos = combinations_with_replacement(g.vertices, k)
        # combos = combinations(g.vertices, k)
        for combo in combos:
            d = {}
            for number in combo:
                d[number] = d.get(number, 0) + 1
            for v in g.vertices:
                if d.get(v, 0) == 0:
                    d[v] = 0
            div = Divisor(g, d)
            works = True
            for v in g.vertices:
                if not div.getReducedDivisor(v).allButOne(v, 1):
                    works = False
            if works:
                return k, div.labels
        k += 1
