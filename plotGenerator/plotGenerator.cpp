#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

using namespace std;

typedef long int ulint;

bool pointComp (pair < pair <double, double>, double > lhs, pair < pair <double, double>, double > rhs) {
    return lhs.second < rhs.second;
}

int main (int argc, char * argv[]) {
    string path, N, D, K, T;

    if (argc == 6) {
        path = string(argv[1]);
        N = string(argv[2]);
        D = string(argv[3]);
        K = string(argv[4]);
        T = string(argv[5]);
    } else {
        cin >> path >> N >> D >> K >> T;
    }

    ifstream inputFile ("../input/instanceN" + N + "D" + D + "K" + K + "T" + T + ".in");
    ifstream resultFile ("../" + path + "/output/N" + N + "D" + D + "K" + K + "T" + T + "/result.out");

    if (inputFile.is_open() && resultFile.is_open()) {
        ulint n, mComplete, m, k;
        inputFile >> n >> mComplete >> m >> k;
        vector < pair < pair <double, double>, double > > vertices (n, make_pair(make_pair(0, 0), 0));
        cout << "Vertices: " << endl;
        // reading vertices' coordinates and penalty
        for (ulint v = 0; v < n; v++) {
            inputFile >> vertices[v].first.first;
            inputFile >> vertices[v].first.second;
            inputFile >> vertices[v].second;
            cout << v << " : (" << vertices[v].first.first << ", " << vertices[v].first.second << ")" << endl;
        }
        // limiting penalities to the interval [0.25, 1]
        double minPenalty = vertices[0].second;
        double maxPenalty = vertices[0].second;
        for (int v = 1; v < n; v++) {
            if (minPenalty > vertices[v].second) {
                minPenalty = vertices[v].second;
            }
            if (maxPenalty < vertices[v].second) {
                maxPenalty = vertices[v].second;
            }
        }
        for (int v = 0; v < n; v++) {
            vertices[v].second -= minPenalty;
        }
        minPenalty -= minPenalty;
        maxPenalty -= minPenalty;
        for (int v = 0; v < n; v++) {
            vertices[v].second *= (1-0.25)/maxPenalty;
        }
        minPenalty *= (1-0.25)/maxPenalty;
        maxPenalty *= (1-0.25)/maxPenalty;
        for (int v = 0; v < n; v++) {
            vertices[v].second += 0.25;
        }

        // printing vertices' coordinates and penalty
        ofstream pointsFile ("../" + path + "/output/N" + N + "D" + D + "K" + K + "T" + T + "/vertices.txt");
        for (ulint v = 0; v < n; v++) {
            pointsFile << vertices[v].first.first << ' ';
            pointsFile << vertices[v].first.second << ' ';
            pointsFile << vertices[v].second << endl;
        }
        pointsFile.close();
        
        // reading (and ignoring) the complete graph's edges
        for (int j = 0; j < mComplete; j++) {
            int u, v, w;
            inputFile >> u >> v >> w;
        }

        /*
        // reading the graph's edges and printing its coordinates
        ofstream edgesFile ("../" + path + "/output/N" + N + "D" + D + "K" + K + "T" + T + "/edges.txt");
        for (int e = 0; e < m; e++) {
            int u, v, w;
            inputFile >> u >> v >> w;
            edgesFile << vertices[u].first.first << ' ';
            edgesFile << vertices[u].first.second << endl;
            edgesFile << vertices[v].first.first << ' ';
            edgesFile << vertices[v].first.second << endl << endl;
        }
        edgesFile.close();
        */

        ulint nSolution, mSolution;
        double costSolution;
        resultFile >> nSolution >> mSolution >> costSolution;
        // reading the solution's vertices and printing its coordinates and neighborhood radio
        ofstream solutionVerticesFile ("../" + path + "/output/N" + N + "D" + D + "K" + K + "T" + T + "/solutionVertices.txt");
        for (ulint i = 0; i < nSolution; i++) {
            int v;
            resultFile >> v;
            solutionVerticesFile << vertices[i].first.first << ' ';
            solutionVerticesFile << vertices[i].first.second << ' ';
            solutionVerticesFile << ((double) k)/100.0 << endl;
        }
        solutionVerticesFile.close();

        cout << "Edges:" << endl;
        // reading the solution's edges and printing its coordinates
        ofstream solutionEdgesFile ("../" + path + "/output/N" + N + "D" + D + "K" + K + "T" + T + "/solutionEdges.txt");
        for (ulint e = 0; e < mSolution; e++) {
            ulint u, v;
            resultFile >> u >> v;
            cout << u << " <-> " << v << " : (" << vertices[u].first.first << ", " << vertices[u].first.second << ") <-> (" << vertices[v].first.first << ", " << vertices[v].first.second << ")" << endl;
            solutionEdgesFile << vertices[u].first.first << ' ';
            solutionEdgesFile << vertices[u].first.second << endl;
            solutionEdgesFile << vertices[v].first.first << ' ';
            solutionEdgesFile << vertices[v].first.second << endl << endl;
        }
        solutionEdgesFile.close();
    }
    return 0;
}
