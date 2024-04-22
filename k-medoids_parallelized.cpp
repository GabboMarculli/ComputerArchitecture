#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <sstream>
#include <fstream>
#include <chrono>
#include <stdlib.h>
#include <cfloat>
using namespace std;

#define MAX_ITERATIONS 10
#define NUM_ARGUMENTS 5
#define NUM_EXECUTIONS 36

// Definition of the points array to represent a point in the dataset
double** points = nullptr;

class KMedoidsData {
private:
    int numPoints; // Number of points in the dataset
    int numDimensions; // Number of dimensions for each point

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
        numDimensions = 0;
        while (getline(file, line)) {
            stringstream ss(line);
            numDimensions = 0;
            string coord;
            // Count the number of dimensions in the point
            while (getline(ss, coord, ','))
                numDimensions++;

            if (numDimensions == 0) {
                cerr << "Error: Number of dimensions is inconsistent among points" << endl;
                return;
            }

            numPoints++;
        }
        file.close();

        // Allocate memory for points
        points = new double*[numPoints];
        for (int i = 0; i < numPoints; ++i)
            points[i] = new double[numDimensions];

        // Reload the file to get the actual points
        file.open(filename);
        int idx = 0;
        while (getline(file, line)) {
            stringstream ss(line);
            int dimIdx = 0;
            string coord;

            // Parse each coordinate of the point
            while (getline(ss, coord, ','))
                points[idx][dimIdx++] = stod(coord);

            idx++;
        }
        file.close();
    }

    // Method to get the points
    double** getPoints() const {
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

    // Destructor
    ~KMedoidsData() {
        if (points != nullptr) {
            for (int i = 0; i < numPoints; ++i) 
                delete[] points[i];
            delete[] points;
        }
    }
};


// Function to calculate the Euclidean distance between two points
// Parameters:
// - point1: Pointer to the coordinates of the first point
// - point2: Pointer to the coordinates of the second point
// - numDimensions: Number of dimensions in each point
double computeDistance(const double* point1, const double* point2, const int numDimensions) {
    double distance = 0.0;
    for (int i = 0; i < numDimensions; ++i) 
        distance += (point1[i] - point2[i]) * (point1[i] - point2[i]);

    return sqrt(distance); // va tolto???
}

// Function to initialize medoids for the K-Medoids algorithm
// Parameters:
// - nClusters: Number of clusters (K)
// - points: Array of pointers to the dataset points
// - numPoints: Number of points in the dataset
// - numDimensions: Number of dimensions in each point
// - medoids: Array of pointers to store the medoids
// - dimensione: Dimension of each medoid (same as numDimensions)
void initializeMedoids(const int nClusters, double** points, const int numPoints, const int numDimensions, double** medoids, const int dimensione) {
    // Allocate memory for medoids
    for (int i = 0; i < nClusters; ++i) 
        medoids[i] = new double[dimensione];

    // Initialize medoids at random positions
    int* medoidsIndices = new int[numPoints];
    for (int i = 0; i < numPoints; ++i)
        medoidsIndices[i] = i;

    random_shuffle(medoidsIndices, medoidsIndices + numPoints); // shuffle indices randomly

    // Assign points as medoids
    for (int i = 0; i < nClusters; ++i)
        for (int j = 0; j < numDimensions; ++j) 
            medoids[i][j] = points[medoidsIndices[i]][j];
}

// Function to print the coordinates of each medoid
// Parameters:
// - medoids: Array of pointers to the medoids
// - nClusters: Number of clusters
// - dimensione: Dimension of each medoid (same as numDimensions)
void printMedoids(double** medoids, const int nClusters, const int dimensione){
    for (int i = 0; i < nClusters; ++i) {
        cout << "centroid " << i << ": ";

        // Print the coordinates of the medoid
        for (int j = 0; j < dimensione; ++j) 
            cout << medoids[i][j] << ",";

        cout << endl;
    }
}

// Function to assign each point to its nearest medoid and update the clusters accordingly
// Parameters:
// - points: Array of pointers to the dataset points
// - npunti: Number of points in the dataset
// - dimensione: Dimension of each point
// - medoids: Array of pointers to the medoids
// - nClusters: Number of clusters
// - clusters: Array to store the assigned cluster index for each point
// - nthread: Number of threads for parallelization
void calculateClusters(double** points, const int npunti, const int dimensione, double** medoids, const int nClusters, int* clusters, const int nthread) {
    int min_index;
    double minDistance;

    // OpenMP parallelization for loop
    #pragma omp parallel for default(shared) private(min_index, minDistance) num_threads(nthread)
    for (int i = 0; i < npunti; i++) {
        minDistance = numeric_limits<double>::max(); // Initialize minDistance to maximum possible value
        for (int j = 0; j < nClusters; j++) {
            double distance = computeDistance(points[i], medoids[j], dimensione); // Calculate distance between point and medoid
            if (distance < minDistance) { // If distance is smaller than the current minimum distance
                minDistance = distance; // Update minDistance
                min_index = j; // Update index of the nearest medoid
            }
        }
        clusters[i] = min_index; // Assign the point to the nearest cluster
    }
}


// Function to update the medoids based on the assigned clusters
// Parameters:
// - points: Array of pointers to the dataset points
// - npunti: Number of points in the dataset
// - dimensione: Dimension of each point
// - medoids: Array of pointers to the medoids
// - nClusters: Number of clusters
// - clusters: Array storing the assigned cluster index for each point
// - nthread: Number of threads for parallelization
void updateMedoids(double** points, const int npunti, const int dimensione, double** medoids, const int nClusters, int* clusters, const int nthread) {
    double totalDistanceArray[nClusters];
    int medoidIndexArray[nClusters];
    // OpenMP parallelization for loop
    #pragma omp parallel default(none) shared(points, clusters, totalDistanceArray, medoidIndexArray, npunti, nClusters, dimensione) num_threads(nthread)
    {
        double localTotalDistanceArray[nClusters];
        int localMedoidIndexArray[nClusters];
        // Initialize local arrays
        for (int i = 0; i < nClusters; i++) {
            localTotalDistanceArray[i] = DBL_MAX;
            localMedoidIndexArray[i] = -1;
        }

        // Each thread calculates distances and updates its local arrays
        #pragma omp for
        for (int j = 0; j < npunti; j++) {
            for (int i = 0; i < nClusters; i++) {
                if (clusters[j] == i) { // Check if the point belongs to the current cluster
                    double distanceToCluster = 0.0;

                    // Calculate the total distance of the current point to all other points in the same cluster
                    for (int k = 0; k < npunti; k++) {
                        if (clusters[k] == i)
                            distanceToCluster += computeDistance(points[j], points[k], dimensione);
                    }

                    // If the total distance is smaller than the current total distance or it's the first iteration,
                    // update the total distance and the index of the new medoid
                    if (distanceToCluster < localTotalDistanceArray[i]) {
                        localTotalDistanceArray[i] = distanceToCluster;
                        localMedoidIndexArray[i] = j;
                    }
                }
            }
        }

        // Combine local arrays into global arrays using a critical section
        #pragma omp critical
        {
            for (int i = 0; i < nClusters; i++) {
                if (localTotalDistanceArray[i] < totalDistanceArray[i]) {
                    totalDistanceArray[i] = localTotalDistanceArray[i];
                    medoidIndexArray[i] = localMedoidIndexArray[i];
                }
            }
        }
    }
}



// Function to perform K-Medoids clustering algorithm
// Parameters:
// - points: Array of pointers to the dataset points
// - npunti: Number of points in the dataset
// - dimensione: Dimension of each point
// - nClusters: Number of clusters (K)
// - nthread: Number of threads for parallelization
// Returns:
// - Time taken for the algorithm to run in milliseconds
long long kMedoids(double** points, const int npunti, const int dimensione, const int nClusters, const int nthread) {
    double** medoids = new double*[nClusters];
    // Initialize medoids
    initializeMedoids(nClusters, points, npunti, dimensione, medoids, dimensione);

    int iter, clusters[npunti];
    auto start = chrono::steady_clock::now(); // Start measuring time

    // Run K-Medoids algorithm for a maximum number of iterations
    for (iter = 0; iter < MAX_ITERATIONS; iter++) {
        calculateClusters(points, npunti, dimensione, medoids, nClusters, clusters, nthread);
        updateMedoids(points, npunti, dimensione, medoids, nClusters, clusters, nthread);
    }

    auto end = chrono::steady_clock::now(); // Stop measuring time
    auto duration =  std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count(); // Calculate duration in milliseconds

    // Print information about the execution
    //printf("nPunti: %d, Time: %lld ms, nThread: %d\n", npunti, duration, nthread);

    // Print the coordinates of the medoids
    //printMedoids(medoids, nClusters, dimensione);

    // Deallocate memory for medoids
    for (int i = 0; i < nClusters; ++i)
        delete[] medoids[i];
    delete[] medoids;

    return duration; // Return the time taken for the algorithm to run
}

// Function to calculate the standard deviation of an array of values
// Parameters:
// - values: Array of values
// - mean: Mean value of the values array
// Returns:
// - Standard deviation of the values
long long standardDeviation(long long values[], long long mean){
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
long long min_value(long long values[]){
    long long minimum = values[0];
    // Iterate through the array to find the minimum value
    for (int i = 1; i < NUM_EXECUTIONS; i++) 
        if(values[i] < minimum) // If the current value is less than the current minimum
            minimum = values[i]; // Update the minimum value

    return minimum; // Return the minimum value found
}

// Function to find the maximum value in an array of values
// Parameters:
// - values: Array of values
// Returns:
// - Maximum value in the array
long long max_value(long long values[]){
    long long maximum = values[0];
    // Iterate through the array to find the maximum value
    for (int i = 1; i < NUM_EXECUTIONS; i++) 
        if(values[i] > maximum) // If the current value is greater than the current maximum
            maximum = values[i]; // Update the maximum value

    return maximum; // Return the maximum value found
}

int main(int argc, char* argv[]) {
    if (argc != NUM_ARGUMENTS) {
        cout << "Wrong number of parameters" << endl;
        return -1;
    }

    // Parse command line arguments
    int npunti = std::atoi(argv[1]); // Number of points
    int dimensione = std::atoi(argv[2]); // Dimension of each point
    int nClusters = std::atoi(argv[3]); // Number of clusters (K)
    int nthread = std::atoi(argv[4]); // Number of threads for parallelization

    // Specify the path to the file containing the points
    string filename = "dataset/points" + to_string(npunti) + "_" + to_string(dimensione) + "_" + to_string(nClusters) + ".txt";
    // Create an instance of KMedoidsData to load the points from the file
    KMedoidsData kmedoidsData(filename);

    // Execute the K-Medoids algorithm multiple times to calculate statistics
    long long total = 0;
    long long values[NUM_EXECUTIONS];
    for(int i = 0; i < NUM_EXECUTIONS; i++){
        // Run K-Medoids algorithm and record the execution time
        long long execution_time = kMedoids(kmedoidsData.getPoints(), npunti, dimensione, nClusters, nthread);
        total += execution_time;
        values[i] = execution_time; // Store execution time for each iteration
    }
    //cout << "Total execution time: " << total << endl;

    // Calculate the arithmetic mean of execution times
    long long mean = total / NUM_EXECUTIONS;
    //cout << "Arithmetic mean: " << mean << endl;

    // Calculate and print standard deviation of execution times
    //cout << "Standard deviation: " << standardDeviation(values, mean) << endl;
    // Find and print minimum execution time
    //cout << "Minimum value: " << min_value(values) << endl;
    // Find and print maximum execution time
    //cout << "Maximum value: " << max_value(values) << endl;

    cout<< npunti << "," << nClusters << "," << nthread << "," << total << "," << mean << "," << standardDeviation(values, mean) << ',' << min_value(values) << ',' << max_value(values) << endl;

    return 0;
}


