import heapq
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from collections import defaultdict

#data
#adjacency list
def graph_adjacency(filename):
    adj = defaultdict(list) 
    node_coords = {}        
    all_edges_geom = []    
        
    coord_to_id = {}
    next_id = 0

    #data load
    with open(filename, 'r') as f:
        for line in f:
            parts = line.split()
            if len(parts) < 5: continue #stop split if less then 5
            
            x1, y1, x2, y2, w = map(float, parts)
            p1 = (x1, y1)
            p2 = (x2, y2)
          
          #assign id
            if p1 not in coord_to_id:
                coord_to_id[p1] = next_id
                node_coords[next_id] = p1
                next_id += 1
            if p2 not in coord_to_id:
                coord_to_id[p2] = next_id
                node_coords[next_id] = p2
                next_id += 1
                
            u, v = coord_to_id[p1], coord_to_id[p2]
           
            #append edges both ways
            adj[u].append((w, v))
            adj[v].append((w, u))
          
            all_edges_geom.append([p1, p2])
            
    return adj, next_id, node_coords, all_edges_geom

#Jarnik/Prime algorithm
def jarnik_prime_mst(num_nodes, graph_adj, start_node=0):
    mst_edges = []      
    total_weight = 0
    
    visited = [False] * num_nodes
    
    #priority que
    pq = [(0, -1, start_node)]
        
    visit_count = 0
    
    while pq and visit_count < num_nodes:
        weight, u, v = heapq.heappop(pq)
      
        if visited[v]:
            continue #skip uf visited
      
        visited[v] = True
        visit_count += 1
       
        if u != -1:
            mst_edges.append((u, v))
            total_weight += weight

        #add neighbors to que
        for w_neighbor, neighbor_id in graph_adj[v]:
            if not visited[neighbor_id]:
                heapq.heappush(pq, (w_neighbor, v, neighbor_id))
                
    return mst_edges, total_weight

#graph
if __name__ == '__main__':
    file = 'graph.txt'
    
    #data load
    adj, num_nodes, coords, bg_lines =  graph_adjacency(file)
    print(f'Number of nodes: {num_nodes}.')
    
    #start
    mst_links, mst_weight = jarnik_prime_mst(num_nodes, adj, start_node=0)
    
    #logs
    print(f'Result:')
    print(f'Weight: {mst_weight:.2f}')
    print(f'Number of edges: {len(mst_links)}')
    
    
    #visualization
    fig, ax = plt.subplots(figsize=(12, 10))

    #road background
    lc_bg = LineCollection(bg_lines, colors='lightgray', linewidths=0.5, alpha=0.5)
    ax.add_collection(lc_bg)
    
    #min tree visualization    
    mst_lines_geom = []
    for u, v in mst_links:
        p1 = coords[u]
        p2 = coords[v]
        mst_lines_geom.append([p1, p2])
        
    lc_mst = LineCollection(mst_lines_geom, colors='red', linewidths=1.5, label='Minimum spanning tree Jarník/Prime')
    ax.add_collection(lc_mst)
    
    #graph settings
    ax.autoscale()
    ax.set_aspect('equal')
    plt.title(f"Minimum spanning tree Jarník/Prim\nTotal weight: {mst_weight:.1f}")
    plt.legend()
    plt.show()
    
   