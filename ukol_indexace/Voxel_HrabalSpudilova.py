# Prostorová indexace - Geoinformatika
# Hrabal, Spudilová (1.N-GKDPZ)
# Bodové mračno: tree_18.txt
# Metoda voxelizace

import math
import time
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict

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

def get_grid_params(X, Y, Z, n):
    # Bounding box limits
    x_min, x_max = min(X), max(X)
    y_min, y_max = min(Y), max(Y)
    z_min, z_max = min(Z), max(Z)
    
    # Point cloud dimensions
    dx = x_max - x_min
    dy = y_max - y_min
    dz = z_max - z_min
    
    # Number of cells per edge
    nr = int(n**(1/3))
    if nr < 1: nr = 1
    
    return x_min, x_max, y_min, y_max, z_min, z_max, dx, dy, dz, nr

def get_voxel_indices(px, py, pz, x_min, x_max, y_min, y_max, z_min, z_max, nr):
    # Coordinate normalization
    denom_x = (x_max - x_min) if (x_max - x_min) != 0 else 1
    denom_y = (y_max - y_min) if (y_max - y_min) != 0 else 1
    denom_z = (z_max - z_min) if (z_max - z_min) != 0 else 1

    x_norm = (px - x_min) / denom_x
    y_norm = (py - y_min) / denom_y
    z_norm = (pz - z_min) / denom_z
    
    # 3D Index calculation
    c = 0.99
    jx = int(x_norm * nr * c)
    jy = int(y_norm * nr * c)
    jz = int(z_norm * nr * c)
    
    return jx, jy, jz

def get_hash(jx, jy, jz, nr):
    # Hashing function
    return jx + (jy * nr) + (jz * (nr**2))

def compute_curvature_pca(neighbors_coords):
    if len(neighbors_coords) < 3: return 0
    
    data = np.array(neighbors_coords)
    cov_matrix = np.cov(data, rowvar=False)
    eig, _ = np.linalg.eig(cov_matrix)
    eig = np.sort(eig)
    
    l1, l2, l3 = eig[0], eig[1], eig[2]
    s = l1 + l2 + l3
    
    return l1 / s if s != 0 else 0

# Compute spatial density and curvature

nearest_list = []
curvatures = []
knn = 30

time_start = time.time()

x_min, x_max, y_min, y_max, z_min, z_max, dx, dy, dz, nr = get_grid_params(X, Y, Z, points_sum)

H = defaultdict(list)
for i in range(points_sum):
    jx, jy, jz = get_voxel_indices(X[i], Y[i], Z[i], x_min, x_max, y_min, y_max, z_min, z_max, nr)
    h_idx = get_hash(jx, jy, jz, nr)
    H[h_idx].append(i)

for i in range(points_sum):
    jx, jy, jz = get_voxel_indices(X[i], Y[i], Z[i], x_min, x_max, y_min, y_max, z_min, z_max, nr)
    
    candidates = []
    
    for kx in [-1, 0, 1]:
        for ky in [-1, 0, 1]:
            for kz in [-1, 0, 1]:
                nh = get_hash(jx+kx, jy+ky, jz+kz, nr)
                if nh in H:
                    for idx in H[nh]:
                        if i == idx: continue
                        d_sq = (X[i]-X[idx])**2 + (Y[i]-Y[idx])**2 + (Z[i]-Z[idx])**2
                        candidates.append((d_sq, idx))
    
    candidates.sort(key=lambda x: x[0])
    
    # density
    if candidates:
        nearest_list.append(math.sqrt(candidates[0][0]))
    else:
        nearest_list.append(0)
        
    # curvature
    k_subset = candidates[:knn]
    if len(k_subset) >= 3:
        coords = [[X[c[1]], Y[c[1]], Z[c[1]]] for c in k_subset]
        curvatures.append(compute_curvature_pca(coords))
    else:
        curvatures.append(0)

voxel_time = time.time() - time_start

if nearest_list:
    valid = [d for d in nearest_list if d > 0]
    if valid:
        d_aver = sum(valid) / len(valid)
        pc_avg_dens = 1 / (d_aver**3)
    else: pc_avg_dens = 0
else: pc_avg_dens = 0

print(f"Čas: {voxel_time:.4f} s")
print(f"Hustota: {pc_avg_dens}")
if curvatures:
    print(f"Křivost: {sum(curvatures)/len(curvatures):.4f}")

test_nrs = [10, 20, 30, 40, 50, 60] 
measured_times = []
voxel_sizes = [] 
sample_size = 1000

for test_nr in test_nrs:
    bx = dx / test_nr
    voxel_sizes.append(bx)
    start = time.time()
    
    temp_H = defaultdict(list)
    for i in range(points_sum):
        jx, jy, jz = get_voxel_indices(X[i], Y[i], Z[i], x_min, x_max, y_min, y_max, z_min, z_max, test_nr)
        h = get_hash(jx, jy, jz, test_nr)
        temp_H[h].append(i)
        
    for i in range(sample_size):
        jx, jy, jz = get_voxel_indices(X[i], Y[i], Z[i], x_min, x_max, y_min, y_max, z_min, z_max, test_nr)
        min_sq = float("inf")
        for kx in [-1, 0, 1]:
            for ky in [-1, 0, 1]:
                for kz in [-1, 0, 1]:
                    nh = get_hash(jx+kx, jy+ky, jz+kz, test_nr)
                    if nh in temp_H:
                        for idx in temp_H[nh]:
                            if i == idx: continue
                            ds = (X[i]-X[idx])**2 + (Y[i]-Y[idx])**2 + (Z[i]-Z[idx])**2
                            if ds < min_sq: min_sq = ds
        if min_sq < float("inf"): math.sqrt(min_sq)
    
    measured_times.append(time.time() - start)

plt.figure(1, figsize=(10, 6))
plt.plot(voxel_sizes, measured_times, 'g-o', label='Voxelizace')
plt.title('Závislost času na velikosti voxelu')
plt.xlabel('Velikost voxelu (m)')
plt.ylabel('Čas (s)')
plt.grid(True)
plt.legend()

test_Ns = [1000, 5000, 10000, 15000, 20000, 25000]
measured_times_N = []
test_Ns = [n for n in test_Ns if n <= points_sum]

for n in test_Ns:
    sub_X = X[:n]
    sub_Y = Y[:n]
    sub_Z = Z[:n]
    
    start_N = time.time()
    
    sub_nr = int(n**(1/3))
    sx_min, sx_max, sy_min, sy_max, sz_min, sz_max, s_dx, s_dy, s_dz, _ = get_grid_params(sub_X, sub_Y, sub_Z, n)
    
    sub_H = defaultdict(list)
    for i in range(n):
        jx, jy, jz = get_voxel_indices(sub_X[i], sub_Y[i], sub_Z[i], sx_min, sx_max, sy_min, sy_max, sz_min, sz_max, sub_nr)
        h = get_hash(jx, jy, jz, sub_nr)
        sub_H[h].append(i)
        
    for i in range(n):
        jx, jy, jz = get_voxel_indices(sub_X[i], sub_Y[i], sub_Z[i], sx_min, sx_max, sy_min, sy_max, sz_min, sz_max, sub_nr)
        min_sq = float("inf")
        for kx in [-1, 0, 1]:
            for ky in [-1, 0, 1]:
                for kz in [-1, 0, 1]:
                    nh = get_hash(jx+kx, jy+ky, jz+kz, sub_nr)
                    if nh in sub_H:
                        for idx in sub_H[nh]:
                            if i == idx: continue
                            ds = (sub_X[i]-sub_X[idx])**2 + (sub_Y[i]-sub_Y[idx])**2 + (sub_Z[i]-sub_Z[idx])**2
                            if ds < min_sq: min_sq = ds
        if min_sq < float("inf"): math.sqrt(min_sq)

    measured_times_N.append(time.time() - start_N)

plt.figure(2, figsize=(10, 6))
plt.plot(test_Ns, measured_times_N, 'b-o', label='Voxelizace')
plt.title('Závislost výpočetního času na počtu bodů')
plt.xlabel('Počet bodů')
plt.ylabel('Čas (s)')
plt.grid(True)
plt.legend()
plt.show()