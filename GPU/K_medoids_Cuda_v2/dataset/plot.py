import matplotlib.pyplot as plt

# Function to load points from a file
def load_points_from_file(filename):
    points = []
    with open(filename, 'r') as file:
        for line in file:
            x, y = map(float, line.strip().split(','))
            points.append((x, y))
    return points

# Load points from dataset file
filename = 'd:/Università/AIDE/Computer Architecture/K_medoids_3/dataset/points2000_2_3.txt'
points = load_points_from_file(filename)

# Plotting the original points
plt.figure(figsize=(8, 6))
plt.scatter([p[0] for p in points], [p[1] for p in points], color='blue', label='Original Points')

# Generate new set of points (example)
new_points = [(-97.7222,278.147), (490.94,93.563), (-39.0594,-113.349)]

# Plotting the new points
plt.scatter([p[0] for p in new_points], [p[1] for p in new_points], color='red', label='New Points')

# Customize plot
plt.title('Original and New Points')
plt.xlabel('X')
plt.ylabel('Y')
plt.legend()

# Show plot
plt.grid(True)
plt.show()
