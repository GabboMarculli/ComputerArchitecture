#include <iostream>
#include <cmath>
#include <algorithm>
#include <sstream>
#include <fstream>
#include <omp.h>
#include <stdlib.h>

using namespace std;

#define maxIterations 100
#define nClusters 3
#define dimensione 2

void loadPoints(const string& filename, double* points, int npunti) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Errore nell'apertura del file" << endl;
        return;
    }

    int index = 0;
    string line;
    while (getline(file, line)) {
        stringstream ss(line);
        string coord;
        int dim = 0;
        while (getline(ss, coord, ',')) {
            if (index >= npunti * dimensione) {
                cerr << "Troppi punti nel file rispetto alla dimensione specificata" << endl;
                return;
            }
            points[index++] = stod(coord);
            dim++;
        }
        if (dim != dimensione) {
            cerr << "Errore: Il numero di dimensioni non è consistente tra i punti" << endl;
            return;
        }
    }
    file.close();
}

int main(int argc, char* argv[]) {
    if(argc != 3){
        cout<<"Numero di parametri sbagliato"<<endl;
        return -1;
    }   
    int npunti = std::atoi(argv[1]);
    int nthread = std::atoi(argv[2]);

    double* points = new double[npunti*dimensione];
    string filename = "points" + to_string(npunti) + "_" + to_string(dimensione) + "_" + to_string(nClusters) + ".txt"; 

    loadPoints(filename, points, npunti);

    for(int i = 0; i < npunti; i++){
        for(int j = 0; j < dimensione; j++){
        }
    }

    double* centroids = new double[nClusters*dimensione];
    for (int i = 0; i < nClusters; ++i) {
        for(int j = 0; j < dimensione; j++){
            centroids[i*dimensione+j] = points[i*dimensione+j];
        }
    }

    int iter, i, j, min_index, kl;
    double distance, minDistance;
    double clusterSum[nClusters][dimensione+1];
    int* clusters = new int[npunti];    
    double start = omp_get_wtime();
    // Inizializzazione dei cluster all'inizio di ogni iterazione
    for (iter = 0; iter < maxIterations; iter++) {

        // Svuota i cluster all'inizio di ogni iterazione
        #pragma omp parallel for num_threads(nthread)
            for (i = 0; i < npunti; ++i) {
                clusters[i] = 0;
            }

        // Calcolo dei nuovi cluster
        #pragma omp parallel for default(shared) private(j, min_index, distance, minDistance) num_threads(nthread)
        for (i = 0; i < npunti; i++) {
            minDistance = numeric_limits<double>::max();
            for (j = 0; j < nClusters; j++){
                distance = (sqrt((points[i*dimensione] - centroids[j*dimensione])*(points[i*dimensione] - centroids[j*dimensione]) + (points[i*dimensione+1]- centroids[j*dimensione+1])*(points[i*dimensione+1]- centroids[j*dimensione+1])));
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
        #pragma omp parallel for default(shared) private(j) reduction(+:clusterSum) num_threads(nthread)
        for (i = 0; i < npunti; i++) {
            for (j = 0; j < dimensione; j++) {
                clusterSum[clusters[i]][j] += points[i*dimensione+j];
            }
            clusterSum[clusters[i]][dimensione]++;
        }

        // Normalize clusterSum
        #pragma omp single
        for (i = 0; i < nClusters; i++) {
            for (j = 0; j < dimensione; j++) {
                centroids[i*dimensione+j] = clusterSum[i][j]/clusterSum[i][dimensione];
            }
        }
    }

    double end = omp_get_wtime();
    printf("%d, %d, %f\n", npunti, nthread, end-start);
    // Stampa i centroidi e i punti appartenenti a ciascun cluster
    /*for(int i = 0; i < nClusters; i++){
        cout<<centroids[i*dimensione]<<","<<centroids[i*dimensione+1]<<endl;
    }*/
    return 0;
}