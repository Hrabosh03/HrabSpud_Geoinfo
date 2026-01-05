import osmnx as ox
import pandas as pd
from math import radians, cos, sin, asin, sqrt

#calculating orthodromic distance using Haversine formula
def ortho_dis(lon1, lat1, lon2, lat2):
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])
    dlon = lon2 - lon1 
    dlat = lat2 - lat1 
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a)) 
    r = 6371000 
    return c * r

def data():
    place_name = 'Okres Kutná Hora, Czechia'
    output_file = 'graph.txt'
    
    
    #variant = 'time'      
    variant = 'distance'   
    #variant = 'sinuosity'  
    
    print(f'Preparing the data for variant: {variant}')

    #max speed map
    speed_map = {
        'motorway': 130, 'motorway_link': 80,
        'trunk': 110, 'trunk_link': 80,
        'primary': 90, 'primary_link': 80,
        'secondary': 90, 'secondary_link': 80,
        'tertiary': 90, 'tertiary_link': 80,
        'unclassified': 90, 'residential': 50,
        'living_street': 20, 'service': 30, 'track': 30
    }
    default_speed = 50

    #downloading the data
    G = ox.graph_from_place(place_name, network_type='drive')
    nodes, edges = ox.graph_to_gdfs(G, nodes=True, edges=True)

    with open(output_file, 'w') as f:
        for index, row in edges.iterrows():
            u, v = index[0], index[1]
            try:
                x1, y1 = nodes.loc[u].x, nodes.loc[u].y
                x2, y2 = nodes.loc[v].x, nodes.loc[v].y
            except KeyError: continue

            #edges 
            highway = row['highway']
            if isinstance(highway, list): highway = highway[0]
            
            real_length = row['length'] 
            speed_kmh = speed_map.get(highway, default_speed)
            speed_ms = speed_kmh / 3.6
            
            #W calculation based on chosen variant
            if variant == 'distance':                
                weight = real_length
                
            elif variant == 'time':                
                weight = real_length / speed_ms
                
            elif variant == 'sinuosity':               
                air_distance = ortho_dis(x1, y1, x2, y2)                
                
                if air_distance < 0.1: 
                    kappa = 1.0
                else:
                    kappa = real_length / air_distance
                                
                base_time = real_length / speed_ms
                
                weight = kappa * base_time

           
            f.write(f"{x1} {y1} {x2} {y2} {weight:.4f}\n")

    print(f" . The file {output_file} is prepared.")

if __name__ == '__main__':
    data()