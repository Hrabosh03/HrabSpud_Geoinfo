# Prostorová indexace - Geoinformatika
# Hrabal, Spudilová (1.N-GKDPZ)
# Bodové mračno: tree_18.txt
# Metoda kd-tree

import math
import time
import matplotlib.pyplot as plt
import numpy as np

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

X, Y, Z = loadPoints(file)
points_sum = len(X)     # number of points

class KDNode:
    def __init__(self, point, left = None, right = None, axis = 0, index = None):
        self.point = point      # [x, y, z]
        self.left = left
        self.right = right
        self.axis = axis
        self.index = index      # store original index to identify self

def build_kdtree(points_with_index, depth = 0):
    if not points_with_index:
        return None
    
    k = 3
    axis = depth % k
    
    # Sort points by current axis
    points_with_index.sort(key = lambda x: x[axis])
    median = len(points_with_index) // 2
    
    node_data = points_with_index[median]
    
    return KDNode(
        point = node_data[:3],
        left = build_kdtree(points_with_index[:median], depth + 1),
        right = build_kdtree(points_with_index[median + 1:], depth + 1),
        axis = axis,
        index = node_data[3]   # original index
    )

def get_dist_sq(p1, p2):
    return (p1[0] - p2[0])**2 + (p1[1] - p2[1])**2 + (p1[2] - p2[2])**2

def search_knn(node, target, target_idx, k, heap):
    if node is None:
        return

    dist_sq = get_dist_sq(node.point, target)
    
    if node.index != target_idx:
        heap.append((dist_sq, node.point))
        heap.sort(key = lambda x: x[0])
        if len(heap) > k:
            heap.pop()
            
    axis = node.axis
    diff = target[axis] - node.point[axis]
    
    if diff < 0:
        near, far = node.left, node.right
    else:
        near, far = node.right, node.left
        
    search_knn(near, target, target_idx, k, heap)
    
    if len(heap) < k or (diff**2 < heap[-1][0]):
        search_knn(far, target, target_idx, k, heap)

def compute_curvature_pca(neighbors_coords):
    if len(neighbors_coords) < 3: return 0
    
    data = np.array(neighbors_coords)
    cov_matrix = np.cov(data, rowvar=False)
    eig, _ = np.linalg.eig(cov_matrix)
    eig = np.sort(eig)
    
    l1, l2, l3 = eig[0], eig[1], eig[2]
    s = l1 + l2 + l3
    
    return l1 / s if s != 0 else 0

# Prepare data
points_with_ids = [[X[i], Y[i], Z[i], i] for i in range(points_sum)]

nearest_list = []
curvatures = []
knn = 30    # k-nearest neighbors for curvature

time_build = time.time()
tree_root = build_kdtree(list(points_with_ids)) 
time_start = time.time()

for i in range(points_sum):
    target = [X[i], Y[i], Z[i]]
    neighbors = []
    
    search_knn(tree_root, target, i, knn, neighbors)
    
    # density
    if len(neighbors) > 0:
        d = math.sqrt(neighbors[0][0])
        nearest_list.append(d)
    else:
        nearest_list.append(0)
        
    # curvature
    if len(neighbors) >= 3:
        coords = [n[1] for n in neighbors]
        curvatures.append(compute_curvature_pca(coords))
    else:
        curvatures.append(0)

kdtree_time = time.time() - time_start

# Compute results
if nearest_list:
    d_aver = sum(nearest_list) / len(nearest_list)
    if d_aver > 0:
        pc_avg_dens = 1 / (d_aver**3)
    else:
        pc_avg_dens = 0
else: pc_avg_dens = 0

print(f"Čas: {kdtree_time:.4f} s")
print(f"Hustota: {pc_avg_dens}")
if curvatures:
    print(f"Křivost: {sum(curvatures)/len(curvatures):.4f}")

tests = [1000, 5000, 10000, 15000, 20000, 25000]
measured_times = []
tests = [n for n in tests if n <= points_sum]

for size in tests:
    sub_points = points_with_ids[:size]
    
    start = time.time()
    
    root = build_kdtree(list(sub_points))
    
    for i in range(size):
        heap = []
        target = sub_points[i][:3]
        target_idx = sub_points[i][3]
        search_knn(root, target, target_idx, 1, heap)
        
    duration = time.time() - start
    measured_times.append(duration)
    
plt.figure(figsize=(10, 6))
plt.plot(tests, measured_times, 'r-o', label='Kd-tree')
plt.title('Závislost výpočetního času na počtu bodů')
plt.xlabel('Počet bodů')
plt.ylabel('Čas (s)')
plt.grid(True)
plt.legend()
plt.show()