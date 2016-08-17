#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <utility>
#include <vector>
#include <sstream>
#include <set>
#include <map>

using namespace std;

typedef long int ulint;

string itos(ulint i) {
    stringstream s;
    s << i;
    return s.str();
}

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

        if (N == itos(n)) {
            cout << "N == n : OK" << endl;
        } else {
            cout << "N == n : ERROR" << endl;
        }

        if (K == itos(k)) {
            cout << "K == k : OK" << endl;
        } else {
            cout << "K == k : ERROR" << endl;
        }

        if (R == itos(r)) {
            cout << "R == r : OK" << endl;
        } else {
            cout << "R == r : ERROR" << endl;
        }

        vector < pair <double, double> > verticesCoordinates (n, make_pair(0, 0));
        vector <ulint> penalties (n, 0);
        // reading vertices' coordinates and penalty
        for (ulint v = 0; v < n; v++) {
            inputFile >> verticesCoordinates[v].first;
            inputFile >> verticesCoordinates[v].second;
            inputFile >> penalties[v];
        }

        // reading (and ignoring) the complete graph's edges
        for (ulint j = 0; j < mComplete; j++) {
            ulint u, v, w;
            inputFile >> u >> v >> w;
        }

        vector < pair <ulint, ulint> > edges (m, make_pair(0, 0));
        map < pair <ulint, ulint>, ulint > weights;
        // reading the graph's edges and printing its coordinates
        for (ulint e = 0; e < m; e++) {
            inputFile >> edges[e].first;
            inputFile >> edges[e].second;
            ulint w;
            inputFile >> w;
            weights[edges[e]] = w;
        }

        ulint nSolution, mSolution, costSolution;
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
        vector <ulint> solutionVerticesInvalid;
        // reading the solution's vertices
        for (ulint v = 0; v < nSolution; v++) {
            resultFile >> solutionVertices[v];
            if (solutionVertices[v] < 0 || solutionVertices[v] >= n) {
                solutionVerticesInvalid.push_back(solutionVertices[v]);
            }
        }
        if (solutionVerticesInvalid.size() <= 0) {
            cout << "Solution Vertices : OK" << endl;
        } else {
            cout << "Solution Vertices : ERROR - ";
            for (ulint i = 0; i < (ulint) solutionVerticesInvalid.size(); i++) {
                cout << solutionVerticesInvalid[i] << " ";
            }
            cout << endl;
        }

        set <ulint> setSolutionVertices (solutionVertices.begin(), solutionVertices.end());
        vector < pair <ulint, ulint> > solutionEdges (mSolution, make_pair(0, 0));
        vector < pair <ulint, ulint> > solutionEdgesInvalid;
        // reading the solution's edges
        for (ulint e = 0; e < mSolution; e++) {
            resultFile >> solutionEdges[e].first;
            resultFile >> solutionEdges[e].second;
            if (setSolutionVertices.find(solutionEdges[e].first) == setSolutionVertices.end() 
            || setSolutionVertices.find(solutionEdges[e].second) == setSolutionVertices.end()) {
                solutionEdgesInvalid.push_back(solutionEdges[e]);
            }
        }
        if (solutionEdgesInvalid.size() <= 0) {
            cout << "Solution Edges : OK" << endl;            
        } else {
            cout << "Solution Edges : ERROR - ";
            for (ulint i = 0; i < (ulint) solutionEdgesInvalid.size(); i++) {
                cout << "(" << solutionEdgesInvalid[i].first << ", " << solutionEdgesInvalid[i].second << ") ";
            }
            cout << endl;
        }

        ulint cost = 0;
        for (ulint e = 0; e < mSolution; e++) {
            cost += weights[solutionEdges[e]];
        }

        ulint allPenalties = 0;
        for (ulint v = 0; v < n; v++) {
            allPenalties += penalties[v];
        }

        ulint chosenPenalties = 0;
        for (ulint v = 0; v < nSolution; v++) {
            chosenPenalties += penalties[solutionVertices[v]];
        }

        cost += (allPenalties - chosenPenalties);

        if (costSolution == cost) {
            cout << "Solution Cost : OK" << endl;
        } else {
            cout << "Solution Cost : ERROR" << endl;
        }
    }
    return 0;
}
