author = 'Surabhi S Nath'
import numpy as np
import sys

"""
This module implements the intricacy measure using https://www.geeksforgeeks.org/program-to-count-number-of-connected-components-in-an-undirected-graph/.
Refer to Nath et al., 2023, Relating Objective Complexity, Subjective Complexity and Beauty in Binary Pixel Patterns, PACA, Section 2.2 (5) and AII for details.
"""

class Graph:
    """
    This class defines the class graph
    """
    def __init__(self, V):
 
        self.V = V      # No. of vertices
        self.adj = [[] for i in range(self.V)]      # Pointer to an array containing adjacency lists
 
    # Function to return the number of connected components in an undirected graph
    def NumberOfconnectedComponents(self):
        """
        Calculates the number of components in the pattern

        Returns:
            integer: the number of components
        """
        # Mark all the vertices as not visited
        visited = [False for i in range(self.V)]
        
        # To store the number of connected components
        count = 0
         
        for v in range(self.V):
            if (visited[v] == False):
                self.DFSUtil(v, visited)
                count += 1
                 
        return count
         
    def DFSUtil(self, v, visited):
        """
        Recursive function for DFS

        Args:
            v (integer): vertex id
            visited (array): indicates if the vertex is visited
        """
        # Mark the current node as visited
        visited[v] = True
 
        # Recur for all the vertices adjacent to this vertex
        for i in self.adj[v]:
            if (not visited[i]):
                self.DFSUtil(i, visited)
                 
    # Add an undirected edge
    def addEdge(self, v, w):
        self.adj[v].append(w)
        self.adj[w].append(v)

def get_neigh(grid, rnum, cnum, dirn, grid_size):   # Origin: top left of matrix
    """
    Gets neighbouring cell based on direction

    Args:
        grid (array): CA grid
        rnum (integer): row index
        cnum (integer): column index
        dirn (integer): 0-7 indicating the direction of nieghbour
        grid_size (integer): grid size of CA grid

    Returns:
        _type_: _description_
    """
    if dirn == 0: # Up
        if rnum > 0:
            return grid[rnum-1, cnum], rnum-1, cnum
    if dirn == 1: # Right
        if cnum < grid_size - 1:
            return grid[rnum, cnum+1], rnum, cnum+1
    if dirn == 2: # Down
        if rnum < grid_size - 1:
            return grid[rnum+1, cnum], rnum+1, cnum
    if dirn == 3: # Left
        if cnum > 0:
            return grid[rnum, cnum-1], rnum, cnum-1
    if dirn == 4: # Up, Right
        if rnum > 0 and cnum < grid_size - 1:
            return grid[rnum-1, cnum+1], rnum-1, cnum+1
    if dirn == 5: # Right, Down
        if rnum < grid_size - 1 and cnum < grid_size - 1:
            return grid[rnum+1, cnum+1], rnum+1, cnum+1
    if dirn == 6: # Down Left
        if rnum < grid_size - 1 and cnum > 0:
            return grid[rnum+1, cnum-1], rnum+1, cnum-1
    if dirn == 7: # Left, Up
        if rnum > 0 and cnum > 0:
            return grid[rnum-1, cnum-1], rnum-1, cnum-1
    return -1, -1, -1

def set_graph_edges(graph, grid, grid_size, neighbourhood_size):
    """
    This functions initialises the graph edges based on neighbourhood size

    Returns:
        Graph: updated graph with edges
    """
    for rnum in range(grid_size):
        for cnum in range(grid_size):
            for dirn in range(neighbourhood_size):
                neigh = get_neigh(grid, rnum, cnum, dirn, grid_size)
                if neigh[0] == grid[rnum, cnum]:
                    graph.addEdge(neigh[1] * grid_size + neigh[2], rnum * grid_size + cnum)
    return graph


def calculate_intricacy(grid, grid_size):
    """
    This function calculates the 4-neighbourhood and 8-neighbourhood intricacy of a pattern.
    
    Returns:
        tuple: 4-neighbourhood and 8-netghbourhood intricacy
    """
    grid = grid.reshape(grid_size, grid_size)
    g1 = Graph(grid_size * grid_size)       # 4-neighbourhood
    g2 = Graph(grid_size * grid_size)       # 8-neighbourhood

    # set graph edges based on neighbourhood size
    g1 = set_graph_edges(g1, grid, grid_size, 4)    # 4-neighbourhood
    g2 = set_graph_edges(g2, grid, grid_size, 8)    # 8-neighbourhood

    return g1.NumberOfconnectedComponents(), g2.NumberOfconnectedComponents()

if __name__ == "__main__":
    grid_size = 27

    # calculate density for some example patterns:
    #    1) Full black pattern
    #    2) White grid with one central black CellType
    #    3) checkerboard
    #    4) Random pattern
    
    all_black = np.ones(grid_size * grid_size)
    assert calculate_intricacy(all_black, grid_size) == (1, 1)
    print(calculate_intricacy(all_black, grid_size))

    grid_with_one_centre = np.zeros(grid_size * grid_size)
    grid_with_one_centre[(grid_size * grid_size) //2] = 1

    # assert calculate_intricacy(grid_with_one_centre, grid_size) == (2, 2)
    print(calculate_intricacy(grid_with_one_centre, grid_size))

    checkerboard = np.zeros(grid_size * grid_size)
    checkerboard[1::2] = 1
    assert calculate_intricacy(checkerboard, grid_size) == (729, 2)
    print(calculate_intricacy(checkerboard, grid_size))

    random = np.random.choice([0, 1], size=(grid_size, grid_size))
    print(calculate_intricacy(random, grid_size))