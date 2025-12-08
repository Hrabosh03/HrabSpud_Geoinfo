# Prostorová indexace - Geoinformatika
# Hrabal M., Spudilová K. (1.N-GKDPZ)
# Bodové mračno: tree_18.txt
# Metoda naivního hledání

import math
import time
import matplotlib.pyplot as plt 

# Load data
file = "tree_18.txt"

def loadPoints (file):
    
    X, Y, Z = [], [], []
    
    with open (file) as f:
        
        for line in f:
            x, y, z = line.split("\t")
            
            X.append(float(x))
            Y.append(float(y))
            Z.append(float(z))
    
    return X, Y, Z

# Compute spatial density
X, Y, Z = loadPoints(file)
points_sum = len(X)     # number of points

nearest_list = []

time_start = time.time()

for i in range(points_sum):
    min_dist_sq = float("inf")
    
    px = X[i]
    py = Y[i]
    pz = Z[i]
    
    for j in range(points_sum):
        if i == j:
            continue        # skip if i is the same as j
        
        dist_sq = (px - X[j])**2 + (py - Y[j])**2 + (pz - Z[j])**2      # compute distance squared
        if dist_sq < min_dist_sq:
            min_dist_sq = dist_sq       # new min distance squared

    nearest_dist = math.sqrt(min_dist_sq)

    nearest_list.append(nearest_dist)        # append nearest distance to the list

time_end = time.time()
naive_time = time_end - time_start

if len(nearest_list) > 0:
    average_dist = sum(nearest_list) / len(nearest_list)
    if average_dist > 0:
        pc_avg_dens = 1 / (average_dist**3)
    else:
        pc_avg_dens = 0
else:
    pc_avg_dens = 0     # compute average density in the point cloud

print(f"Čas: {naive_time:.2f} s")
print(f"Hustota: {pc_avg_dens}")

# Time visualisation
tests = [100, 500, 1000, 2000, 3000, 5000, 7500, 10000]
measured_times = []

for size in tests:
    sub_X = X[:size]
    sub_Y = Y[:size]
    sub_Z = Z[:size]        # create subset of data

    start = time.time()
    
    # Naive algorithm for the subset
    for i in range(size):
        min_dist_sq = float("inf")
        px, py, pz = sub_X[i], sub_Y[i], sub_Z[i]
        
        for j in range(size):
            if i == j: continue
            
            d_sq = (px - sub_X[j])**2 + (py - sub_Y[j])**2 + (pz - sub_Z[j])**2
            
            if d_sq < min_dist_sq:
                min_dist_sq = d_sq
        math.sqrt(min_dist_sq) 
        
    duration = time.time() - start
    measured_times.append(duration)

# Graph
plt.figure(figsize=(10, 6))

plt.plot(tests[:len(measured_times)], measured_times, 'b-o', label='Naivní metoda (O(N²))')

plt.title('Závislost výpočetního času na počtu bodů')
plt.xlabel('Počet bodů')
plt.ylabel('Čas (s)')
plt.grid(True)
plt.legend()
plt.show()
