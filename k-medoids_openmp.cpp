#include <iostream>
#include <cmath>
#include <algorithm>
#include <sstream>
#include <fstream>
#include <omp.h>
#include <stdlib.h>
using namespace std;

#define MAX_ITERATIONS 100
#define NUM_ARGUMENTS 5

// Definizione del vettore di punti per rappresentare un punto nel dataset
double** points = nullptr;

class KMedoidsData {
private:
    int numPoints;
    int numDimensions;

public:
    // Costruttore che prende il percorso del file e carica i punti
    KMedoidsData(const string& filename) {
        loadPoints(filename);
    }

    // Metodo per caricare i punti da un file
    void loadPoints(const string& filename) {
        ifstream file(filename);
        if (!file.is_open()) {
            cerr << "Errore nell'apertura del file" << endl;
            return;
        }

        string line;
        numPoints = 0;
        numDimensions = 0;
        while (getline(file, line)) {
            stringstream ss(line);
            numDimensions = 0;
            string coord;
            while (getline(ss, coord, ','))
                numDimensions++;

            if (numDimensions == 0) {
                cerr << "Errore: Il numero di dimensioni non è consistente tra i punti" << endl;
                return;
            }

            numPoints++;
        }
        file.close();

        // Alloca memoria per i punti
        points = new double*[numPoints];
        for (int i = 0; i < numPoints; ++i)
            points[i] = new double[numDimensions];

        // Ricarica il file per ottenere i punti effettivi
        file.open(filename);
        int idx = 0;
        while (getline(file, line)) {
            stringstream ss(line);
            int dimIdx = 0;
            string coord;

            while (getline(ss, coord, ','))
                points[idx][dimIdx++] = stod(coord);

            idx++;
        }
        file.close();
    }

    // Metodo per ottenere i punti
    double** getPoints() const {
        return points;
    }

    // Metodo per ottenere il numero di punti
    int getNumPoints() const {
        return numPoints;
    }

    // Metodo per ottenere il numero di dimensioni
    int getNumDimensions() const {
        return numDimensions;
    }

    // Distruttore per liberare la memoria
    ~KMedoidsData() {
        if (points != nullptr) {
            for (int i = 0; i < numPoints; ++i) 
                delete[] points[i];
            delete[] points;
        }
    }
};

// Funzione per calcolare la distanza euclidea tra due punti
double computeDistance(const double* point1, const double* point2, const int numDimensions) {
    double distance = 0.0;
    for (int i = 0; i < numDimensions; ++i) 
        distance += (point1[i] - point2[i]) * (point1[i] - point2[i]);

    return sqrt(distance); // togliere??
}

// Funzione per inizializzare i medoidi
void initializeMedoids(const int nClusters, double** points, const int numPoints, const int numDimensions, double** medoids, const int dimensione) {
    for (int i = 0; i < nClusters; ++i) 
        medoids[i] = new double[dimensione];

    // Inizializza i medoidi in posizioni casuali
    int medoidsIndices[numPoints];
    for (int i = 0; i < numPoints; ++i)
        medoidsIndices[i] = i;

    random_shuffle(medoidsIndices, medoidsIndices + numPoints);

    for (int i = 0; i < nClusters; ++i)
        for (int j = 0; j < numDimensions; ++j) 
            medoids[i][j] = points[medoidsIndices[i]][j];
}

// Stampa i centroidi
void printMedoids(double** medoids, const int nClusters, const int dimensione){
    for (int i = 0; i < nClusters; ++i) {
        cout << "centroide " << i << ": ";

        for (int j = 0; j < dimensione; ++j) 
            cout << medoids[i][j] << ",";

        cout << endl;
    }
}

void calculateClusters(double** points, const int npunti, const int dimensione, double** medoids, const int nClusters, int* clusters, const int nthread) {
    int min_index;
    double minDistance;

    #pragma omp parallel for default(shared) private(min_index, minDistance) num_threads(nthread)
    for (int i = 0; i < npunti; i++) {
        minDistance = numeric_limits<double>::max();
        for (int j = 0; j < nClusters; j++) {
            double distance = computeDistance(points[i], medoids[j], dimensione);
            if (distance < minDistance) {
                minDistance = distance;
                min_index = j;
            }
        }
        clusters[i] = min_index;
    }
}

void updateMedoids(double** points, const int npunti, const int dimensione, double** medoids, const int nClusters, int* clusters, const int nthread) {
    double totalDistance;
    int newMedoidIndex;

    #pragma omp parallel for default(shared) private(totalDistance, newMedoidIndex) num_threads(nthread)
    for (int i = 0; i < nClusters; i++) {
        totalDistance = 0.0;
        newMedoidIndex = -1;
        for (int j = 0; j < npunti; j++) {
            if (clusters[j] == i) {
                double distanceToCluster = 0.0;

                for (int k = 0; k < npunti; k++)
                    if (clusters[k] == i) 
                        distanceToCluster += computeDistance(points[j], points[k], dimensione);

                if (distanceToCluster < totalDistance || newMedoidIndex == -1) {
                    totalDistance = distanceToCluster;
                    newMedoidIndex = j;
                }
            }
        }
        for (int k = 0; k < dimensione; ++k) 
            medoids[i][k] = points[newMedoidIndex][k];
    }
}

void kMedoids(double** points, const int npunti, const int dimensione, const int nClusters, const int nthread) {
    double** medoids = new double*[nClusters];
    initializeMedoids(nClusters, points, npunti, dimensione, medoids, dimensione);

    int iter;
    int clusters[npunti];

    double start = omp_get_wtime();

    for (iter = 0; iter < MAX_ITERATIONS; iter++) {
        calculateClusters(points, npunti, dimensione, medoids, nClusters, clusters, nthread);
        updateMedoids(points, npunti, dimensione, medoids, nClusters, clusters, nthread);
    }

    double end = omp_get_wtime();
    printf("nPunti: %d, Time: %f, nThread: %d\n", npunti, end - start, nthread);

    printMedoids(medoids, nClusters, dimensione);

    // Liberare la memoria allocata per i medoidi
    for (int i = 0; i < nClusters; ++i)
        delete[] medoids[i];
    delete[] medoids;
}

int main(int argc, char* argv[]) {
    if (argc != NUM_ARGUMENTS) {
        cout << "Numero di parametri sbagliato" << endl;
        return -1;
    }

    int npunti = std::atoi(argv[1]);
    int dimensione = std::atoi(argv[2]);
    int nClusters = std::atoi(argv[3]);
    int nthread = std::atoi(argv[4]);

    // Specifica il percorso del file contenente i punti
    string filename = "points" + to_string(npunti) + "_" + to_string(dimensione) + "_" + to_string(nClusters) + ".txt";
    KMedoidsData kmedoidsData(filename);

    // Esegui l'algoritmo K-Medoids
    kMedoids(kmedoidsData.getPoints(), npunti, dimensione, nClusters, nthread);

    return 0;
}
