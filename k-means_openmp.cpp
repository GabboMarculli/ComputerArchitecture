#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <sstream>
#include <fstream>
#include <omp.h>
#include <stdlib.h>

using namespace std;

#define maxIterations 100
#define dimensione 2
#define nClusters 3

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

int main(int argc, char* argv[]) {
    if(argc != 3){
        cout<<"Numero di parametri sbagliato"<<endl;
        return -1;
    }   
    int npunti = std::atoi(argv[1]);
    int nthread = std::atoi(argv[2]);

    // Definisce il dataset
    // Specifica il percorso del file contenente i punti
    string filename = "points" + to_string(npunti) + "_" + to_string(dimensione) + "_" + to_string(nClusters) + ".txt"; 
    KMeansData kmeansData(filename);
    // Ottieni i punti e il numero di dimensioni
    vector<Point> points = kmeansData.getPoints();

    // Inizializza i centroidi in posizioni casuali
    vector<Point> centroids;
    for (int i = 0; i < nClusters; ++i) {
        centroids.push_back(points[i]);
        //cout<<"centroide " <<i<< ": "<<centroids[i].coordinates[0]<<","<<centroids[i].coordinates[1]<<endl;
    }
    
    vector<Point> oldCentroids;
    int iter, i, j, min_index, kl;
    double distance, minDistance;
    double clusterSum[nClusters][dimensione+1];
    int* clusters = new int[npunti];    
    double start = omp_get_wtime();

    omp_set_num_threads(nthread);
    // Inizializzazione dei cluster all'inizio di ogni iterazione
    for (iter = 0; iter < maxIterations; iter++) {

        // Svuota i cluster all'inizio di ogni iterazione
        #pragma omp parallel for
            for (i = 0; i < npunti; ++i) {
                clusters[i] = 0;
            }

        // Calcolo dei nuovi cluster
        #pragma omp parallel for default(shared) private(j, min_index, distance, minDistance)
        for (i = 0; i < npunti; i++) {
            minDistance = numeric_limits<double>::max();
            for (j = 0; j < nClusters; j++){
                distance = ((points[i].coordinates[0]-centroids[j].coordinates[0])*(points[i].coordinates[0]-centroids[j].coordinates[0]) + (points[i].coordinates[1]-centroids[j].coordinates[1])*(points[i].coordinates[1]-centroids[j].coordinates[1]));
                if(distance < minDistance){
                    minDistance = distance;
                    min_index = j;
                }
            }
            clusters[i] = min_index;
        }

        // Resetto ClusterSum
        #pragma omp master
        for (i = 0; i < nClusters; i++){
            for(j = 0; j<dimensione+1; j++){
                clusterSum[i][j] = 0;
            }
        }

        // Aggiornamento dei centroidi
        #pragma omp parallel for default(shared) private(j) reduction(+:clusterSum)
        for (i = 0; i < npunti; i++) {
            for (j = 0; j < dimensione; j++) {
                clusterSum[clusters[i]][j] += points[i].coordinates[j];
            }
            clusterSum[clusters[i]][dimensione]++;
        }

        // Normalize clusterSum
        #pragma omp single
        for (i = 0; i < nClusters; i++) {
            for (j = 0; j < dimensione; j++) {
                centroids[i].coordinates[j] = clusterSum[i][j]/clusterSum[i][dimensione];
            }
        }
    }

    double end = omp_get_wtime();
    printf("%d, %d, %f\n", npunti, nthread, end-start);
    // Stampa i centroidi e i punti appartenenti a ciascun cluster

    for (int i = 0; i < nClusters; ++i) {
        cout<<"centroide " <<i<< ": "<<centroids[i].coordinates[0]<<","<<centroids[i].coordinates[1]<<endl;
    }
    
    return 0;
}