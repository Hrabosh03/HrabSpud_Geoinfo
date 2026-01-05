from queue import PriorityQueue
from math import inf
import osmnx as ox

#data import
try:
    from lines_to_graph2 import G, PSE
    print("Data was sucessfully imported.")
except ImportError:
    print("Error: file not found")
    exit()

#dijkstra algorithm
def dijkstra(G, start, end):
   #attributes init 
    n = len(PSE)
    
    d = [inf] * n     
    p = [-1] * n      
    
    Q = PriorityQueue()
    #start node    
    d[start] = 0
    Q.put((0, start)) 

    #get closest node
    while not Q.empty():        
        du, u = Q.get()
        
        #break when end node is found     
        if u == end:
            break
                 
        if du > d[u]:
            continue

        
        if u in G:
            for v, wuv in G[u].items():
                if d[v] > d[u] + wuv: 
                    d[v] = d[u] + wuv 
                    p[v] = u          
                    Q.put((d[v], v))  
                    
    return p, d

#backwards path reconstruction
def pathreconstruction(p, start, end):    
    path = []
    curr = end
    
    while curr != start and curr != -1:
        path.append(curr)
        curr = p[curr]
    
    path.append(start)
    
    return path[::-1]

#start and end point allocation
def start_end_point(name, point_list_pse):   
    try:        
        lat, lon = ox.geocode(name)
        print(f"{name} -> {lat:.5f}, {lon:.5f}")
        
        closest_id = -1
        min_dist = float('inf')
        
        for i in range(len(point_list_pse)):            
            x = point_list_pse[i][0]
            y = point_list_pse[i][1]
                       
            if x > 1000: continue
            
            dist = (x - lon)**2 + (y - lat)**2
            
            if dist < min_dist:
                min_dist = dist
                closest_id = i
                
        return closest_id

    except Exception as e:
        print(f"Error: '{name}' not found.")
        return None

#input municipality
start = input('Enter start point (municipality): ')
end = input('Enter end point (municipality): ') 

start_id = start_end_point(start, PSE)
end_id = start_end_point(end, PSE)

#path finding ĺogs
if start_id is not None and end_id is not None:
    print('Finding a path..')
        
    predecessors, distances = dijkstra(G, start_id, end_id)
    finW = distances[end_id]
    
    if finW == inf:
        print('\nPath does not exist')
        path_id = []
    else:
        path_id = pathreconstruction(predecessors, start_id, end_id)
        
        print('\nPath found.')
        
        #min = int(finW // 60)
        #sec = int(finW % 60)
        #print(f'Time: {min} min {sec} s ({finW:.2f} s)')
        print(f'Final weight: {finW}')
        print(f'Number of nodes: {len(path_id)}')
else:
    print('Start or end point not found')
    path_id = []

#graph
if path_id:    
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection
    from lines_to_graph2 import PS, PE
        
    fig, ax = plt.subplots(figsize=(20, 16))
    
    lines = []
    for i in range(len(PS)):
        lines.append([(PS[i][0], PS[i][1]), (PE[i][0], PE[i][1])])
        
    lc = LineCollection(lines, colors='lightgray', linewidths=0.5, alpha=0.7)
    ax.add_collection(lc)
    
    path_x = [PSE[id_node][0] for id_node in path_id]
    path_y = [PSE[id_node][1] for id_node in path_id]
    
    ax.plot(path_x, path_y, color='blue', linewidth=2, label='Path', zorder=2)
    ax.scatter(path_x, path_y, color='blue', s=10, zorder=3)
    
    ax.scatter(path_x[0], path_y[0], color='green', s=100, label=start, edgecolors='black', zorder=4)
    ax.scatter(path_x[-1], path_y[-1], color='red', s=100, label=end, edgecolors='black', zorder=4)

    ax.autoscale() 
    ax.set_aspect('equal') 
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.title(f'Path: {start} -> {end}')
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.3)                
    plt.show()

