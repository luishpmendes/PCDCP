#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

using namespace std;

typedef long int ulint;

int main (int argc, char * argv[]) {
    string path, N, D, K, T, R, P;

    cout << argc << endl;

    if (argc == 8) {
        path = string(argv[1]);
        N = string(argv[2]);
        D = string(argv[3]);
        K = string(argv[4]);
        T = string(argv[5]);
        R = string(argv[6]);
        P = string(argv[7]);
    } else {
        cin >> path >> N >> D >> K >> T >> R;
    }

    ifstream inputFile ("../input/instanceN" + N + "D" + D + "K" + K + "T" + T + "R" + R + "P" + P + ".in");
    ifstream resultFile ("../" + path + "/output/N" + N + "D" + D + "K" + K + "T" + T + "R" + R + "P" + P + "/result.out");

    if (inputFile.is_open() && resultFile.is_open()) {
        ulint n, mComplete, m, k, t, r, root;
        double d, p;
        inputFile >> n >> d >> k >> t >> r >> p >> mComplete >> m >> root;

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

        vector < pair <double, double> > verticesCoordinates (n, make_pair(0, 0));
        vector <double> penalties (n, 0);
        // reading vertices' coordinates and penalty
        for (ulint v = 0; v < n; v++) {
            inputFile >> vertices[v].first;
            inputFile >> vertices[v].second;
            inputFile >> penalties[v];
        }

        // reading (and ignoring) the complete graph's edges
        for (ulint j = 0; j < mComplete; j++) {
            ulint u, v, w;
            inputFile >> u >> v >> w;
        }

        vector < pair <double, double> > edges (m, make_pair(0, 0));
        vector <double> weights (m, 0);
        // reading the graph's edges and printing its coordinates
        for (ulint e = 0; e < m; e++) {
            int u, v, w;
            inputFile >> edges[e].first;
            inputFile >> edges[e].second;
            inputFile >> weights[e];
        }

        ulint nSolution, mSolution;
        double costSolution;
        resultFile >> nSolution >> mSolution >> costSolution;

        if (nSolution <= n) {
            cout << "nSolution <= n : OK" << endl;
        } else {
            cout << "nSolution <= n : ERROR" << endl;
        }

        if (mSolution <= m) {
            cout << "mSolution <= m : OK" << endl;
        } else {
            cout << "mSolution <= m : ERROR" << endl;
        }

        vector <ulint> solutionVertices (nSolution, 0);
        // reading the solution's vertices
        for (ulint v = 0; v < nSolution; v++) {
            resultFile >> solutionVertices[v];
        }

        vector < pair <ulint, ulint> > solutionEdges (mSolution, make_pair(0, 0));
        // reading the solution's edges
        for (ulint e = 0; e < mSolution; e++) {
            resultFile >> solutionEdges[e].first;
            resultFile >> solutionEdges[e].second;
        }
    }
    return 0;
}
