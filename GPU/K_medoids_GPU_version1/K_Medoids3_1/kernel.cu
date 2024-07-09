#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <sstream>
#include <fstream>
#include <chrono>
#include <cstdlib>
#include <cfloat>
#include <random>

using namespace std;

#define MAX_ITERATIONS 10
#define NUM_EXECUTIONS 1
#define NUM_ARGUMENTS 3


struct Point {
    double x;
    double y;
};

class KMedoidsData {
private:
    int numPoints; // Number of points in the dataset
    int numDimensions = 2; // Number of dimensions for each point
    vector<Point> points; // Vector of points

public:
    // Constructor that takes the file path and loads the points
    KMedoidsData(const string& filename) {
        loadPoints(filename);
    }

    // Method to load points from a file
    void loadPoints(const string& filename) {
        ifstream file(filename);
        if (!file.is_open()) {
            cerr << "Error opening file" << endl;
            return;
        }

        string line;
        numPoints = 0;
        numDimensions = 2; // We know the struct has 2 dimensions (x, y)

        while (getline(file, line)) {
            stringstream ss(line);
            string coord;
            double coords[2];
            int dimIdx = 0;

            // Parse each coordinate of the point
            while (getline(ss, coord, ',')) {
                if (dimIdx < 2) {
                    coords[dimIdx++] = stod(coord);
                }
                else {
                    cerr << "Error: More than two dimensions in the point" << endl;
                    return;
                }
            }

            if (dimIdx != 2) {
                cerr << "Error: Point does not have two dimensions" << endl;
                return;
            }

            points.push_back({ coords[0], coords[1] });
            numPoints++;
        }
        file.close();
    }

    // Method to get the points
    const vector<Point>& getPoints() const {
        return points;
    }

    // Method to get the number of points
    int getNumPoints() const {
        return numPoints;
    }

    // Method to get the number of dimensions
    int getNumDimensions() const {
        return numDimensions;
    }
};

void printCudaMemoryInfo(const char* stage) {
    size_t free_mem, total_mem;
    cudaMemGetInfo(&free_mem, &total_mem);
    cout << stage << " - Free memory: " << free_mem / (1024 * 1024) << " MB, Total memory: " << total_mem / (1024 * 1024) << " MB" << endl;
}

// Kernel function to calculate the Euclidean distance between two points
__device__ double computeDistance(const Point p1, const Point p2) {
    return sqrt((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y));

}

// Kernel function to calculate the nearest medoid for each point and update the cluster assignment
__global__ void calculateClusters(Point* d_p, int np, Point* d_m, int nc, int* d_c) {
    for (int i = blockIdx.x * blockDim.x + threadIdx.x; i < np; i += blockDim.x * gridDim.x) {
        int min_index;
        double minDistance = DBL_MAX;

        for (int j = 0; j < nc; j++) {
            double distance = computeDistance(d_p[i], d_m[j]);
            if (distance < minDistance) {
                minDistance = distance;
                min_index = j;
            }
        }
        d_c[i] = min_index;
    }
}


__global__ void updateMedoids(Point* d_p, int np, Point* d_m, int nc, int* d_c, double* d_distances, int k) {

    for (int i = blockIdx.x * blockDim.x + threadIdx.x; i < np; i += blockDim.x * gridDim.x) {
        int clusterIdx = d_c[i];
        Point p = d_p[i];
        int start = k;
        int end = k + gridDim.x*1000;
        double distance = 0.0;
        if (start > 0)
            distance = d_distances[i];
        // Calculate distances and find best medoids for each cluster
        for (int j = start; j < end && j < np; ++j) {
            if (d_c[j] == clusterIdx) {
                distance += computeDistance(p, d_p[j]);
            }
        }
        d_distances[i] = distance;
    }
}

void retrieveBestMedoids(Point* d_p, Point* d_m, int nc, int np, double* d_distances, vector<int> c) {
    // Host arrays to store intermediate results
    vector<double> SharedDistances(np);

    // Copy from device to host
    cudaMemcpy(SharedDistances.data(), d_distances, np * sizeof(double), cudaMemcpyDeviceToHost);

    // Final reduction on host
    vector<double> finalMinDistances(nc, DBL_MAX);
    vector<int> finalBestMedoids(nc, -1);

    for (int i = 0; i < np; ++i) {
        //cout << SharedDistances[i] << endl;
        if (SharedDistances[i] < finalMinDistances[c[i]]) {
            finalMinDistances[c[i]] = SharedDistances[i];
            finalBestMedoids[c[i]] = i;
            //cout << i << ": " << finalMinDistances[c[i]] << ", " << finalBestMedoids[c[i]] << endl;
        }
    }

    // Update medoids
    for (int i = 0; i < nc; ++i) {
        if (finalBestMedoids[i] != -1) {
            cudaMemcpy(&d_m[i], &d_p[finalBestMedoids[i]], sizeof(Point), cudaMemcpyDeviceToDevice);
        }
    }
}

vector<Point> initializeMedoids(vector<Point> p, int np, int nc) {

    vector<Point> medoids;

    // Initialize medoids at random positions
    int* medoidsIndices = new int[np];
    for (int i = 0; i < np; ++i)
        medoidsIndices[i] = i;

    random_shuffle(medoidsIndices, medoidsIndices + np); // shuffle indices randomly

    // Assign points as medoids
    for (int i = 0; i < nc; ++i)
        medoids.push_back(p[i]);

    return medoids;
}

long long kMedoids(vector<Point> p, int np, int nc, int nt, int nb) {

    printCudaMemoryInfo("Start execution");

    int size_p = np * sizeof(Point);
    Point* d_p;
    cudaMalloc((void**)&d_p, size_p);
    cudaMemcpy(d_p, p.data(), size_p, cudaMemcpyHostToDevice);

    // Launch the kernel with enough threads to cover all points
    int threadsPerBlock = nt;
    int blocksPerGrid = nb;
    //cout << threadsPerBlock << ", " << blocksPerGrid << endl;

    vector<Point> medoids = initializeMedoids(p, np, nc);
    cout << "Old Medoids" << endl;
    for (int i = 0; i < medoids.size(); i++) {
        cout << medoids[i].x << "," << medoids[i].y << endl;
    }

    int size_m = nc * sizeof(Point);
    Point* d_m;
    cudaMalloc((void**)&d_m, size_m);
    cudaMemcpy(d_m, medoids.data(), size_m, cudaMemcpyHostToDevice);

    int size_c = np * sizeof(int);
    int* d_c;
    vector<int> cluster(np);
    cudaMalloc((void**)&d_c, size_c);
    cudaMemset(d_c, -1, size_c);

    // Allocate memory for intermediate results
    double* d_distances;
    cudaMalloc(&d_distances, np * sizeof(double));
    cudaMemset(d_distances, -1, np * sizeof(double));

    printCudaMemoryInfo("After Creation");

    // Execute the K-Medoids algorithm on GPU
    auto start = chrono::steady_clock::now();

    int m = nb*1000 * ((np / (nb*1000)) + 1);

    for (int i = 0; i < MAX_ITERATIONS; i++) {
        calculateClusters << <blocksPerGrid, threadsPerBlock >> > (d_p, np, d_m, nc, d_c);
        cudaDeviceSynchronize();

        cudaMemcpy(cluster.data(), d_c, size_c, cudaMemcpyDeviceToHost);
        /*for (int j = 0; j < np; j++) {
            if (cluster[j] == -1) {
                cout << j << " NO CLUSTER" << endl;
            }
        }*/

        for (int k = 0; k < m; k += nb*1000) {
            updateMedoids << <blocksPerGrid, threadsPerBlock >> > (d_p, np, d_m, nc, d_c, d_distances, k);
            cudaDeviceSynchronize();
        }

        retrieveBestMedoids(d_p, d_m, nc, np, d_distances, cluster);
        cudaMemcpy(medoids.data(), d_m, size_m, cudaMemcpyDeviceToHost);


    }

    // Retrieve updated medoids from device to host
    cudaMemcpy(medoids.data(), d_m, size_m, cudaMemcpyDeviceToHost);

    auto end = chrono::steady_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(end - start).count();

    cout << "Final Medoids" << endl;
    for (int i = 0; i < nc; i++) {
        cout << medoids[i].x << "," << medoids[i].y << endl;
    }

    cudaFree(d_p);
    cudaFree(d_m);
    cudaFree(d_c);
    cudaFree(d_distances);

    printCudaMemoryInfo("After deletion");

    return duration;
}

// Function to calculate the standard deviation of an array of values
// Parameters:
// - values: Array of values
// - mean: Mean value of the values array
// Returns:
// - Standard deviation of the values
long long standardDeviation(long long values[], long long mean) {
    long long std_deviation = 0;
    // Calculate the sum of squared differences from the mean
    for (int i = 0; i < NUM_EXECUTIONS; i++)
        std_deviation += pow(values[i] - mean, 2);

    std_deviation /= NUM_EXECUTIONS; // Divide by the number of values to get the mean of squared differences
    return sqrt(std_deviation); // Calculate the square root to get the standard deviation
}

// Function to find the minimum value in an array of values
// Parameters:
// - values: Array of values
// Returns:
// - Minimum value in the array
long long min_value(long long values[]) {
    long long minimum = values[0];
    // Iterate through the array to find the minimum value
    for (int i = 1; i < NUM_EXECUTIONS; i++)
        if (values[i] < minimum) // If the current value is less than the current minimum
            minimum = values[i]; // Update the minimum value

    return minimum; // Return the minimum value found
}

// Function to find the maximum value in an array of values
// Parameters:
// - values: Array of values
// Returns:
// - Maximum value in the array
long long max_value(long long values[]) {
    long long maximum = values[0];
    // Iterate through the array to find the maximum value
    for (int i = 1; i < NUM_EXECUTIONS; i++)
        if (values[i] > maximum) // If the current value is greater than the current maximum
            maximum = values[i]; // Update the maximum value

    return maximum; // Return the maximum value found
}

int main(int argc, char* argv[])
{
    /*if (argc != NUM_ARGUMENTS) {
        cout << "Wrong number of parameters" << endl;
        return -1;
    }*/

    int arr_np[] = { 1000 };
    int n_t[] = { 32 };
    int n_b[] = { 1, 2, 4, 8, 16, 32, 64, 128, 512, 1024 };
    // Parse command line arguments
    for (int i_np = 0; i_np < 10; i_np++) {
        for (int i_t = 0; i_t < 1; i_t++) {
                int np = arr_np[i_t]; // std::atoi(argv[1]); // Number of points
                int nc = 5;// std::atoi(argv[2]); // Number of clusters (K)

                // Specify the path to the file containing the points
                string filename = "../dataset/points" + to_string(np) + "_" + to_string(2) + "_" + to_string(nc) + ".txt";
                // Create an instance of KMedoidsData to load the points from the file
                KMedoidsData kmedoidsData(filename);
                long long total = 0.0;
                long long values[NUM_EXECUTIONS];

                const vector<Point>& p = kmedoidsData.getPoints();

                for (int i = 0; i < NUM_EXECUTIONS; i++) {
                    long long execution_time = kMedoids(p, np, nc, n_t[i_t], n_b[i_np]);
                    total += execution_time;
                    values[i] = execution_time;
                }
                // Calculate the arithmetic mean of execution times
                long long mean = total / NUM_EXECUTIONS;
                cout << np << "," << nc << "," << n_t[i_t] << "," << n_b[i_np] << "," << total << "," << mean << "," << standardDeviation(values, mean) << ',' << min_value(values) << ',' << max_value(values) << endl;

                // Prepare results in CSV format
                string output_filename = "results.csv";
                ofstream outfile(output_filename, ios::app); // Open file for appending
                if (!outfile.is_open()) {
                    cerr << "Error opening output file: " << output_filename << endl;
                    return -1;
                }

                // Write results to CSV file
                outfile << np << "," << nc << "," << n_t[i_t] << "," << n_b[i_np] << "," << total << "," << mean << ","
                    << standardDeviation(values, mean) << ',' << min_value(values) << ',' << max_value(values) << endl;

                // Close the file
                outfile.close();
        }
    }

    return 0;
}
