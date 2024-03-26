#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <sstream>
#include <fstream>
#include <omp.h>

using namespace std;

// Definizione della struttura Point per rappresentare un punto nel dataset
struct Point {
    vector<double> coordinates;
    Point(const vector<double>& coords) : coordinates(coords) {}
};

class KMeansData {
private:
    vector<Point> points;
    int numDimensions;
public:
    // Costruttore che prende il percorso del file e carica i punti
    KMeansData(const string& filename) {
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
        while (getline(file, line)) {
            stringstream ss(line);
            vector<double> coords;
            numDimensions = 0;
            string coord;
            while (getline(ss, coord, ',')) {
                coords.push_back(stod(coord));
            }
            if (numDimensions == 0) {
                numDimensions = coords.size();
            } else if (coords.size() != numDimensions) {
                cerr << "Errore: Il numero di dimensioni non è consistente tra i punti" << endl;
                return;
            }
            points.push_back(Point(coords));
        }
        file.close();
    }
    // Metodo per ottenere i punti
    vector<Point> getPoints() const {
        return points;
    }
    // Metodo per ottenere il numero di dimensioni
    int getNumDimensions() const {
        return numDimensions;
    }
};
// Calcola la distanza euclidea tra due punti
double distance(const Point& p1, const Point& p2, int numDimensions) {
    double somma = 0;
    for(int i = 0; i < numDimensions; i++){
        somma += pow(p1.coordinates[i] - p2.coordinates[i], 2);
    }
    return sqrt(somma);
}

// Trova il centroide più vicino per un dato punto
int findClosestCentroid(const Point& point, const vector<Point>& centroids, int numDimensions) {
    double minDist = numeric_limits<double>::max();
    int closestCentroid = -1;
    for (int i = 0; i < centroids.size(); ++i) {
        double dist = distance(point, centroids[i], numDimensions);
        if (dist < minDist) {
            minDist = dist;
            closestCentroid = i;
        }
    }
    return closestCentroid;
}

// Aggiorna la posizione dei centroidi
void updateCentroids(vector<Point>& centroids, const vector<vector<Point>>& clusters, int numDimensions) {
    for (int i = 0; i < centroids.size(); ++i) {
        vector<double> sum(numDimensions, 0);
        int clusterSize = clusters[i].size();
        if (clusterSize == 0) continue;
        for (const auto& point : clusters[i]) {
            for(int j = 0; j < numDimensions; j++){
                sum[j] += point.coordinates[j];
            }
        }
        for(int j = 0; j < numDimensions; j++)
        centroids[i].coordinates[j] = sum[j]/clusterSize;
    }
}
// Implementazione di K-Means
vector<vector<Point>> kMeans(const vector<Point>& points, int k, int maxIterations, int numDimensions) {
    // Inizializza i centroidi in posizioni casuali
    vector<Point> centroids;
    //random_shuffle(points.begin(), points.end());
    for (int i = 0; i < k; ++i) {
        centroids.push_back(points[i]);
        cout<<"centroide " <<i<< ": "<<centroids[i].coordinates[0]<<","<<centroids[i].coordinates[1]<<endl;
    }
    vector<vector<Point>> clusters(k);
    // Esegui le iterazioni di K-Means
    for (int iter = 0; iter < maxIterations; iter++) {
        // Assegna ogni punto al cluster più vicino
        for (const auto& point : points) {
            int closestCentroid = findClosestCentroid(point, centroids, numDimensions);
            clusters[closestCentroid].push_back(point);
        }
        // Aggiorna la posizione dei centroidi
        updateCentroids(centroids, clusters, numDimensions);
        
        // Svuota i cluster per la prossima iterazione
        if(!(iter == maxIterations-1))
            for (auto& cluster : clusters) {
                cluster.clear();
            }
    }
    cout << endl;
    return clusters;
}
int main(int argc, char* argv[]) {
    if(argc != 5){
        cout<<"Numero di parametri sbagliato"<<endl;
        return -1;
    }   
    int npunti = std::atoi(argv[1]);
    int dimensione = std::atoi(argv[2]);
    int cluster = std::atoi(argv[3]);
    int nthread = std::atoi(argv[4]);
    // Definisce il dataset
    // Specifica il percorso del file contenente i punti
    string filename = "points" + to_string(npunti) + "_" + to_string(dimensione) + "_" + to_string(cluster) + ".txt"; 
    KMeansData kmeansData(filename);
    // Ottieni i punti e il numero di dimensioni
    vector<Point> points = kmeansData.getPoints();

    int numDimensions = kmeansData.getNumDimensions();

    cout << "Numero di dimensioni: " << numDimensions << endl;
    // Ora puoi applicare l'algoritmo K-Means ai punti caricati
    // Specifica il numero di cluster
    int k = 3;
    // Specifica il numero massimo di iterazioni
    int maxIterations = 5;
    // Esegui l'algoritmo K-Means
    double start = omp_get_wtime();
    vector<vector<Point>> clusters = kMeans(points, k, maxIterations, numDimensions);
    double end = omp_get_wtime();
    // Stampa i centroidi e i punti appartenenti a ciascun cluster
    for (int i = 0; i < k; i++) {
        cout << "Cluster " << i + 1 << " centroid: (" << clusters[i][0].coordinates[0] << ", " << clusters[i][0].coordinates[1] << ")" << endl;
    }
    printf("Time: %f", end-start);
    return 0;
}