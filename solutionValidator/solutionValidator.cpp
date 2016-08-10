#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

using namespace std;

typedef long int ulint;

int main (int argc, char * argv[]) {
    string path, N, D, K, T, R;

    cout << argc << endl;

    if (argc == 7) {
        path = string(argv[1]);
        N = string(argv[2]);
        D = string(argv[3]);
        K = string(argv[4]);
        T = string(argv[5]);
        R = string(argv[6]);
    } else {
        cin >> path >> N >> D >> K >> T >> R;
    }

    ifstream inputFile ("../input/instanceN" + N + "D" + D + "K" + K + "T" + T + "R" + R + ".in");
    ifstream resultFile ("../" + path + "/output/N" + N + "D" + D + "K" + K + "T" + T + "R" + R + "/result.out");

    if (inputFile.is_open() && resultFile.is_open()) {
        ulint n, mComplete, m, k, r;
        inputFile >> n >> mComplete >> m >> k >> r;

        if (N == n) {
            cout << "N == n : OK" << endl;
        } else {
            cout << "N == n : ERROR" << endl;
        }

        if (K == k) {
            cout << "K == k : OK" << endl;
        } else {
            cout << "K == k : ERROR" << endl;
        }

        if (R == r) {
            cout << "R == r : OK" << endl;
        } else {
            cout << "R == r : ERROR" << endl;
        }

        vector < pair < pair <double, double>, double > > vertices (n, make_pair(make_pair(0, 0), 0));
        cout << "Vertices: " << endl;
        // reading vertices' coordinates and penalty
        for (ulint v = 0; v < n; v++) {
            inputFile >> vertices[v].first.first;
            inputFile >> vertices[v].first.second;
            inputFile >> vertices[v].second;
        }

        // reading (and ignoring) the complete graph's edges
        for (int j = 0; j < mComplete; j++) {
            int u, v, w;
            inputFile >> u >> v >> w;
        }

        ulint nSolution, mSolution;
        double costSolution;
        resultFile >> nSolution >> mSolution >> costSolution;
        // reading the solution's vertices and printing its coordinates and neighborhood radio
        ofstream solutionVerticesFile ("../" + path + "/output/N" + N + "D" + D + "K" + K + "T" + T + "R" + R +"/solutionVertices.txt");
        for (ulint i = 0; i < nSolution; i++) {
            int v;
            resultFile >> v;
            solutionVerticesFile << vertices[v].first.first << ' ';
            solutionVerticesFile << vertices[v].first.second << ' ';
            solutionVerticesFile << ((double) k)/100.0 << endl;
        }
        solutionVerticesFile.close();

        cout << "Edges:" << endl;
        // reading the solution's edges and printing its coordinates
        ofstream solutionEdgesFile ("../" + path + "/output/N" + N + "D" + D + "K" + K + "T" + T + "R" + R + "/solutionEdges.txt");
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
