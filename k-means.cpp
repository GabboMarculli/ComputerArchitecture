#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <sstream>
#include <fstream>
using namespace std;
// Definizione della struttura Point per rappresentare un punto nel dataset
/*struct Point {
    double x, y;
    Point(double x_, double y_) : x(x_), y(y_) {}
};*/
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
    }
    vector<vector<Point>> clusters(k);
    // Esegui le iterazioni di K-Means
    for (int iter = 0; iter < maxIterations; ++iter) {
        // Assegna ogni punto al cluster più vicino
        for (const auto& point : points) {
            int closestCentroid = findClosestCentroid(point, centroids, numDimensions);
            clusters[closestCentroid].push_back(point);
        }
        // Aggiorna la posizione dei centroidi
        updateCentroids(centroids, clusters, numDimensions);
        // Svuota i cluster per la prossima iterazione
        for (auto& cluster : clusters) {
            cluster.clear();
        }
    }
    return clusters;
}
int main() {
    // Definisci il dataset di esempio
    //vector<Point> points = { {1, 2}, {1.5, 1.8}, {5, 8}, {8, 8}, {1, 0.6}, {9, 11}, {8, 2}, {10, 2}, {9, 3} };
    // Specifica il percorso del file contenente i punti
    string filename = "points20_2_3.txt";
    // Crea un'istanza di KMeansData e carica i punti dal file
    KMeansData kmeansData(filename);
    // Ottieni i punti e il numero di dimensioni
    vector<Point> points = kmeansData.getPoints();
    int numDimensions = kmeansData.getNumDimensions();
    // Stampa i punti caricati
    cout << "Punti caricati:" << endl;
    for (const auto& point : points) {
        cout << "Punto: ";
        for (const auto& coord : point.coordinates) {
            cout << coord << " ";
        }
        cout << endl;
    }
    cout << "Numero di dimensioni: " << numDimensions << endl;
    // Ora puoi applicare l'algoritmo K-Means ai punti caricati
    // Specifica il numero di cluster
    int k = 3;
    // Specifica il numero massimo di iterazioni
    int maxIterations = 100;
    // Esegui l'algoritmo K-Means
    vector<vector<Point>> clusters = kMeans(points, k, maxIterations, numDimensions);
    // Stampa i centroidi e i punti appartenenti a ciascun cluster
    for (int i = 0; i < clusters.size(); ++i) {
        cout << "Cluster " << i + 1 << " centroid: (" << clusters[i][0].coordinates[0] << ", " << clusters[i][0].coordinates[1] << ")" << endl;
        cout << "Points in cluster " << i + 1 << ":" << endl;
        for (const auto& point : clusters[i]) {
            cout << "(" << point.coordinates[0] << ", " << point.coordinates[1] << ")" << endl;
        }
        cout << endl;
    }
    return 0;
}