import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

#Union find
class UnionFind:
    def __init__(self, n):
        self.parent = list(range(n))
        self.rank = [0] * n  
    
    #path compression
    def find(self, i):
        if self.parent[i] != i:
            self.parent[i] = self.find(self.parent[i])
        return self.parent[i]

    #weighted union
    def union(self, i, j):
        root_i = self.find(i)
        root_j = self.find(j)

        #tree union by ranks
        if root_i != root_j:           
            if self.rank[root_i] < self.rank[root_j]:
                self.parent[root_i] = root_j
            elif self.rank[root_i] > self.rank[root_j]:
                self.parent[root_j] = root_i
            else:
                self.parent[root_j] = root_i
                self.rank[root_i] += 1
            return True 
        return False

#data load
def load_graph_data(filename):
    edges = []
    nodes = set()
    node_coords = {} 

    #temp map to cord transform
    coord_to_id = {}
    next_id = 0

    with open(filename, 'r') as f:
        for line in f:
            parts = line.split()
            if len(parts) < 5: continue
            
            x1, y1, x2, y2, w = map(float, parts)
           
            #add ID to nodes
            p1 = (x1, y1)
            p2 = (x2, y2)
            
            if p1 not in coord_to_id:
                coord_to_id[p1] = next_id
                node_coords[next_id] = p1
                next_id += 1
            if p2 not in coord_to_id:
                coord_to_id[p2] = next_id
                node_coords[next_id] = p2
                next_id += 1
            
            #edge append    
            u, v = coord_to_id[p1], coord_to_id[p2]
            edges.append((w, u, v, p1, p2))
            
    return edges, next_id, node_coords

#Boruvka/Kruskal min spanning tree
def bk_mst(num_nodes, edges):
    uf = UnionFind(num_nodes)
    mst_edges = []
    total_weight = 0
    
    #edge sort and union
    edges.sort(key=lambda x: x[0])
        
    for w, u, v, p1, p2 in edges:
        if uf.union(u, v):
            mst_edges.append((p1, p2))
            total_weight += w
            
    return mst_edges, total_weight

#main
if __name__ == "__main__":
    file = 'graph.txt'
    
    
    all_edges, num_nodes, coords = load_graph_data(file)
    print(f"Number of nodes:{num_nodes}    Total number of edges {len(all_edges)}")
    
    mst_edges, mst_weight = bk_mst(num_nodes, all_edges)
    
    print(f'Result:')
    print(f'Total weight: {mst_weight:.2f}')
    print(f'Number of edges: {len(mst_edges)}')
     
    #graph
    fig, ax = plt.subplots(figsize=(12, 10))

    #roads
    bg_lines = [[e[3], e[4]] for e in all_edges]
    lc_bg = LineCollection(bg_lines, colors='lightgray', linewidths=0.5, alpha=0.5)
    ax.add_collection(lc_bg)
    
    #tree
    lc_mst = LineCollection(mst_edges, colors='red', linewidths=1.5, label='Min spanning tree')
    ax.add_collection(lc_mst)
    
    ax.autoscale()
    ax.set_aspect('equal')
    plt.title(f"Minimum spanning tree Boruvka/Kruskal\nTotal weight: {mst_weight:.1f}")
    plt.legend()
    plt.show()
   
        
   