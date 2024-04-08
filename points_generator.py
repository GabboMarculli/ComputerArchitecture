import numpy as np
import random
import matplotlib.pyplot as plt

def generate_points(k, d, min, max):
    points = set()
    while len(points) < k:
        points.add(tuple(np.random.uniform(min, max, d)))
    points_list = list(points)
    random.shuffle(points_list)
    return points_list

def dataset_generator(mu, sigma, k, d, n, min, max):
    n_points_per_cluster = n / k
    reminder = n % k
    points = set()

    for i in range(0, k):
        initial_len = len(points)
        if ((reminder) != 0) and (i == (k - 1)): # last cluster will have (n % k) points more than others
            n_points_per_cluster += (reminder)

        while len(points) < (initial_len + n_points_per_cluster):
            points.add(tuple(np.random.normal(mu[i], sigma[i], d)))

    points_list = list(points)
    random.shuffle(points_list)
    return points_list  

def scrivi_su_file(punti, nome_file):
    with open(nome_file, 'w') as file:
        for punto in punti:
            file.write(','.join(str(coordinate) for coordinate in punto) + '\n')

# Richiesta delle dimensioni dei punti
#d = int(input("Insert number of dimensions of points: "))

# Richiesta del numero di centroidi
#k = int(input("Insert number of centroids: "))

# Richiesta del numero di punti nel dataset
#n = int(input("Insert number of points in the synthetic dataset: "))

d_arr = [2]
n_arr = [4000000]
k_arr = [3]

for n in n_arr:
    for d in d_arr:
        for k in k_arr:          
            for i in range(1):
                # Genera i punti casuali
                MIN_VALUE = -100
                MAX_VALUE = 100
                mu = generate_points(k, d, MIN_VALUE, MAX_VALUE)
                sigma = generate_points(k, d, (n / 1000*k), (n / 300*k))
                dataset = dataset_generator(mu, sigma, k, d, n, MIN_VALUE, MAX_VALUE)
                dataset = np.array(dataset)
                plt.scatter(dataset[:,0], dataset[:,1])
                plt.show()
                # Scrive i punti su un file esterno
                nome_file = "points" + str(n) + '_' + str(d) +'_' + str(k) + '.txt'
                scrivi_su_file(dataset, nome_file)
                print(f"I punti sono stati scritti su {nome_file}.")
