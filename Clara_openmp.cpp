#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <limits>
#include <cmath>
#include <random>
#include <chrono>
#include <sstream>
#include <omp.h>

using namespace std;

// Define a structure for a point in the dataset
struct Point {
    vector<double> coordinates;

    // Constructor to initialize the coordinates
    Point(const vector<double>& coords) : coordinates(coords) {}
};

// Function to calculate the Euclidean distance between two points
double calculateDistance(const Point& point1, const Point& point2, const int size) {
    double distance = 0.0;
    for (size_t i = 0; i < size; ++i)
        distance += (point1.coordinates[i] - point2.coordinates[i])*(point1.coordinates[i] - point2.coordinates[i]);
    
    return distance; 
}

// Function to calculate the total cost of a set of medoids
double calculateTotalCost(const vector<Point>& dataset, const vector<Point>& medoids, const int dimensione, const int numClusters, const int npunti) {
    double totalCost = 0.0;

    #pragma omp parallel for shared(npunti, dataset,medoids,dimensione)
    for (int i =0; i < npunti; i++) {
        double minDistance = numeric_limits<double>::max();
        for (int j = 0; j < numClusters; j++) {
            double distance = calculateDistance(dataset[i], medoids[j], dimensione);
            if (distance < minDistance) {
                minDistance = distance;
            }
        }
        totalCost += minDistance;
    }
    return totalCost;
}

// Function to generate a random sample of the dataset
vector<Point> generateRandomSample(const vector<Point>& dataset, size_t sampleSize, const int npunti) {
    vector<Point> sample;
    vector<size_t> indices(npunti+sampleSize);

    // Initialize indices from 0 to dataset.size() - 1
    #pragma omp parallel for shared(indices)
    for (size_t i = 0; i < npunti + sampleSize; ++i) {
        indices[i] = i;
    }

    // Shuffle the indices
    #pragma omp parallel for shared(indices)
    for (size_t i = 0; i < npunti + sampleSize; ++i) {
        size_t j = rand() % (i + 1);
        swap(indices[i], indices[j]);
    }

    // Select the first 'sampleSize' indices to form the sample
    #pragma omp parallel for shared(sample, dataset, indices)
    for (size_t i = 0; i < sampleSize; ++i) {
        sample.push_back(dataset[indices[i]]);
    }

    return sample;
}

// Function to calculate the total distance of a point to all other points in the same cluster
double calculateTotalDistance(const Point& point, const vector<Point>& dataset, const vector<size_t>& clusterAssignment, size_t pointIndex, const int dimensione, const int numClusters) {
    double totalDistance = 0.0;
    for (size_t i = 0; i < numClusters; ++i)
        if (clusterAssignment[i] == clusterAssignment[pointIndex])
            totalDistance += calculateDistance(point, dataset[i], dimensione);

    return totalDistance;
}

// Function to perform Partitioning Around Medoids (PAM) clustering
vector<Point> pam(const vector<Point>& dataset, const int numClusters, const int dimensione) {
    vector<Point> medoids;
    
    // Initialize medoids by selecting the first 'numClusters' points as initial medoids
    for (size_t i = 0; i < numClusters; ++i)
        medoids.push_back(dataset[i]);

    bool converged = false;
    vector<size_t> clusterAssignment(numClusters);

    while (!converged) {
        converged = true;

        for (size_t i = 0; i < numClusters; ++i) {
            double minDistance = numeric_limits<double>::max();
            size_t closestMedoidIndex = 0;

            for (size_t j = 0; j < numClusters; ++j) {
                double distance = calculateDistance(dataset[i], medoids[j], dimensione);
                if (distance < minDistance) {
                    minDistance = distance;
                    closestMedoidIndex = j;
                }
            }

            if (clusterAssignment[i] != closestMedoidIndex) {
                clusterAssignment[i] = closestMedoidIndex;
                converged = false;
            }
        }

        // Update medoids by calculating the medoid of each cluster
        for (size_t i = 0; i < numClusters; ++i) {
            Point newMedoid = medoids[i];

            for (size_t j = 0; j < numClusters; ++j) {
                if (clusterAssignment[j] == i) {
                    double totalDistance = calculateTotalDistance(dataset[j], dataset, clusterAssignment, j, dimensione, numClusters);
                    if (totalDistance < calculateTotalDistance(newMedoid, dataset, clusterAssignment, j, dimensione, numClusters)) {
                        newMedoid = dataset[j];
                    }
                }
            }

            medoids[i] = newMedoid;
        }
    }

    return medoids;
}

// Function to perform CLARA clustering
vector<Point> clara(const vector<Point>& dataset, const int numClusters, const int npunti, const int dimensione) {
    vector<Point> bestMedoids;
    double minCost = numeric_limits<double>::max();

    #pragma omp parallel for shared(bestMedoids, minCost, dataset, numClusters, dimensione, npunti) default(none)
    for (size_t i = 0; i < npunti; ++i) {
        // Generate a random sample of the dataset
        vector<Point> sample = generateRandomSample(dataset, numClusters, npunti);

        // Perform PAM (Partitioning Around Medoids) clustering on the sample
        vector<Point> medoids = pam(sample, numClusters, dimensione);

        // Calculate the total cost of the clustering
        double cost = calculateTotalCost(dataset, medoids, dimensione, numClusters, npunti);
        
        // Update the best medoids if the current clustering has lower cost
        #pragma omp critical
        {
            if (cost < minCost) {
                minCost = cost;
                bestMedoids = medoids;
            }
        }
    }

    return bestMedoids;
}

// Function to load points from a file
vector<Point> loadPoints(const string& filename) {
    vector<Point> points;
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error opening file" << endl;
        return points;
    }

    string line;
    while (getline(file, line)) {
        stringstream ss(line);
        vector<double> coordinates;
        string coord;
        while (getline(ss, coord, ',')) {
            coordinates.push_back(stod(coord));
        }
        points.emplace_back(coordinates);
    }
    file.close();

    return points;
}

int main(int argc, char* argv[]) {
    if (argc != 5) {
        cout << "Wrong number of parameters" << endl;
        return -1;
    }

    // Parse command line arguments
    int npunti = std::atoi(argv[1]); // Number of points
    int dimensione = std::atoi(argv[2]); // Dimension of each point
    int nClusters = std::atoi(argv[3]); // Number of clusters (K)
    int nthread = std::atoi(argv[4]); // Number of threads for parallelization

    // Specify the path to the file containing the points
    string filename = "points" + to_string(npunti) + "_" + to_string(dimensione) + "_" + to_string(nClusters) + ".txt";
    // Load points from file
    vector<Point> dataset = loadPoints(filename);

    //for(int i = 0; i < npunti; i++)
    //    cout << dataset[i].coordinates[0] << ' ';

    // Perform CLARA clustering
    auto total = 0;
    for (int iter = 0; iter < 36; iter++) {
        auto start = chrono::steady_clock::now(); // Start measuring time
        vector<Point> medoids = clara(dataset, nClusters, npunti, dimensione);
        auto end = chrono::steady_clock::now(); // Stop measuring time

        for (int i = 0; i < nClusters; ++i) {
            cout << "Centroid " << i << ": ";
            cout << medoids[i].coordinates[0] << ","<< medoids[i].coordinates[1] << ' ';
            cout << endl;
        }

        auto duration =  std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count(); // Calculate duration in milliseconds
        cout << "Tempo: " << duration << endl;
        total+=duration;
    }

    cout << "Tempo totale: " << total ;

    return 0;
}
