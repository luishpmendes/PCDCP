#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <utility>
#include <vector>
#include <sstream>
#include <set>
#include <map>
#include <list>
#include <queue>

#ifndef INFINITE
#define INFINITE 15 << 25
#endif

using namespace std;

typedef long int ulint;
typedef vector < vector <ulint> > matrix;

string itos(ulint i) {
    stringstream s;
    s << i;
    return s.str();
}

vector < set <ulint> > neighbourhoods (matrix W, ulint k) {
    vector < set <ulint> > result (W.size());
/*
    matrix D = matrix (W.size(), vector <ulint> (W.size(), INFINITE));

    for (ulint i = 0; i < W.size(); i++) {
        for (ulint j = 0; j < W.size(); j++) {
            D[i][j] = W[i][j];
        }
    }

    for (ulint l = 0; l < W.size(); l++) {
        for (ulint i = 0; i < W.size(); i++) {
            for (ulint j = 0; j < W.size(); j++) {
                if (D[i][j] > D[i][l] + D[l][j]) {
                    D[i][j] = D[i][l] + D[l][j];
                }
            }
        }
    }

    for (ulint u = 0; u < W.size(); u++) {
        for (ulint v = 0; v < W[u].size(); v++) {
            if (D[u][v] <= k) {
                result[u].insert(v);
            }
        }
    }
*/
    for (ulint u = 0; u < (ulint) W.size(); u++) {
        for (ulint v = 0; v < (ulint) W[u].size(); v++) {
            if (W[u][v] <= k) {
                result[u].insert(v);
            }
        }
    }
    return result;
}

int main (int argc, char * argv[]) {
    string path, N, D, K, T, P;

    cout << argc << endl;

    if (argc == 7) {
        path = string(argv[1]);
        N = string(argv[2]);
        D = string(argv[3]);
        K = string(argv[4]);
        T = string(argv[5]);
        P = string(argv[6]);
    } else {
        cin >> path >> N >> D >> K >> T >> P;
    }

    ifstream inputFile ("../input/instanceN" + N + "D" + D + "K" + K + "T" + T + "P" + P + ".in");
    ifstream resultFile ("../" + path + "/output/N" + N + "D" + D + "K" + K + "T" + T + "P" + P + "/result.out");

    if (inputFile.is_open() && resultFile.is_open()) {
        ulint n, mComplete, m, k, t, root;
        double d, p;
        int errorFlag = 0;
        inputFile >> n >> d >> k >> t >> p >> mComplete >> m >> root;

        if (N == itos(n)) {
            cout << "N == n : OK" << endl;
        } else {
            cout << "N == n : ERROR" << endl;
            errorFlag = 1;
        }

        if (K == itos(k)) {
            cout << "K == k : OK" << endl;
        } else {
            cout << "K == k : ERROR" << endl;
            errorFlag = 1;
        }

        vector < pair < pair <double, double>, ulint > > vertices (n, make_pair(make_pair(0, 0), 0));
        // reading vertices' coordinates and penalty
        for (ulint v = 0; v < n; v++) {
            inputFile >> vertices[v].first.first;
            inputFile >> vertices[v].first.second;
            inputFile >> vertices[v].second;
        }

        matrix W (n, vector <ulint> (n, INFINITE)); // adjacency matrix for the complete graph
        for (ulint i = 0; i < n; i++) {
            W[i][i] = 0;
        }

        // reading (and ignoring) the complete graph's edges
        for (ulint j = 0; j < mComplete; j++) {
            ulint u, v, w;
            inputFile >> u >> v >> w;
            W[u][v] = w;
            W[v][u] = w;
        }

        vector < pair <ulint, ulint> > edges (m, make_pair(0, 0));
        map < pair <ulint, ulint>, ulint > weights;
        // reading the graph's edges and its weights
        for (ulint e = 0; e < m; e++) {
            ulint u, v, w;
            inputFile >> u >> v >> w;
            edges[e].first = u;
            edges[e].second = v;
            weights[edges[e]] = w;
        }

        ulint nSolution, mSolution, costSolution;
        resultFile >> nSolution >> mSolution >> costSolution;

        if (nSolution <= n) {
            cout << "nSolution <= n : OK" << endl;
        } else {
            cout << "nSolution <= n : ERROR" << endl;
            errorFlag = 1;
        }

        if (mSolution <= m) {
            cout << "mSolution <= m : OK" << endl;
        } else {
            cout << "mSolution <= m : ERROR" << endl;
            errorFlag = 1;
        }

        vector <ulint> solutionVertices (nSolution, 0);
        vector <ulint> solutionVerticesInvalid;
        // reading the solution's vertices
        for (ulint v = 0; v < nSolution; v++) {
            resultFile >> solutionVertices[v];
            // check if the solution vertice is a valid vertice
            if (solutionVertices[v] < 0 || solutionVertices[v] >= n) {
                solutionVerticesInvalid.push_back(solutionVertices[v]);
            }
        }
        if (solutionVerticesInvalid.size() <= 0) {
            cout << "Solution Vertices : OK" << endl;
        } else {
            cout << "Solution Vertices : ERROR - ";
            errorFlag = 1;
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
            // check if the ends of the solution edge are solution vertices
            if (setSolutionVertices.find(solutionEdges[e].first) == setSolutionVertices.end() 
            || setSolutionVertices.find(solutionEdges[e].second) == setSolutionVertices.end()) {
                solutionEdgesInvalid.push_back(solutionEdges[e]);
            }
        }
        if (solutionEdgesInvalid.size() <= 0) {
            cout << "Solution Edges : OK" << endl;            
        } else {
            cout << "Solution Edges : ERROR - ";
            errorFlag = 1;
            for (ulint i = 0; i < (ulint) solutionEdgesInvalid.size(); i++) {
                cout << "(" << solutionEdgesInvalid[i].first << ", " << solutionEdgesInvalid[i].second << ") ";
            }
            cout << endl;
        }

        // compute solution cost
        ulint cost = 0;
        for (ulint e = 0; e < mSolution; e++) {
            cost += weights[solutionEdges[e]];
        }

        ulint allPenalties = 0;
        for (ulint v = 0; v < n; v++) {
            allPenalties += vertices[v].second;
        }

        ulint chosenPenalties = 0;
        for (ulint v = 0; v < nSolution; v++) {
            chosenPenalties += vertices[solutionVertices[v]].second;
        }

        cost += (allPenalties - chosenPenalties);

        // check if solution cost is correct
        if (costSolution == cost) {
            cout << "Solution Cost : OK" << endl;
        } else {
            cout << "Solution Cost : ERROR" << endl;
            errorFlag = 1;
        }

        vector <ulint> verticesInSolution (n, 0);

        for (ulint i = 0; i < nSolution; i++) {
            ulint v = solutionVertices[i];
            verticesInSolution[v] = 1;
        }

        vector < list < pair <ulint, ulint> > > adj (n);

        for (ulint e = 0; e < mSolution; e++) {
            ulint u, v, w;
            u = solutionEdges[e].first;
            v = solutionEdges[e].second;
            w = weights[solutionEdges[e]];
            adj[u].push_back(make_pair(v, w));
            adj[v].push_back(make_pair(u, w));
        }

        // check if all vertices are covered
        vector < set <ulint> > Ns = neighbourhoods (W, k);

        set <ulint> uncoveredVertices;
        for (ulint u = 0; u < n; u++) {
            if (verticesInSolution[u] == 0) {
                int flag = 0;
                for (ulint i = 0; i < nSolution && flag == 0; i++) {
                    ulint v = solutionVertices[i];
                    if (Ns[v].find(u) != Ns[v].end()) {
                        flag = 1;
                    }
                }
                if (flag == 0) {
                    uncoveredVertices.insert(u);
                }
            }
        }

        if (uncoveredVertices.size() <= 0) {
            cout << "Cover : OK" << endl;
        } else {
            cout << "Cover : Error - Uncovered vertices: ";
            errorFlag = 1;
            for (set <ulint> :: iterator it = uncoveredVertices.begin(); it != uncoveredVertices.end(); it++) {
                ulint u = *it;
                cout << u << " ";
            }
            cout << endl;
        }

        vector <ulint> verticesInMainCycle (n, 0);

        queue <ulint> Q;

        Q.push(root);

        while (!Q.empty()) {
            ulint u = Q.front();
            Q.pop();

            for (list < pair <ulint, ulint> > :: iterator it = adj[u].begin(); it != adj[u].end(); it++) {
                ulint v = it->first;
                if (verticesInMainCycle[v] == 0) {
                    verticesInMainCycle[v] = 1;
                    Q.push(v);
                }
            }
        }

        vector <ulint> solutionVerticesNotInMainCycle;
        vector <ulint> nonSolutionVerticesInMainCycle;

        for (int v = 0; v < n; v++) {
            if (verticesInSolution[v] == 1) {
                if (verticesInMainCycle[v] == 0) {
                    solutionVerticesNotInMainCycle.push_back(v);
                }
            } else {
                if (verticesInMainCycle[v] == 1) {
                    nonSolutionVerticesInMainCycle.push_back(v);
                }
            }
        }

        if (solutionVerticesNotInMainCycle.size() <= 0) {
            cout << "Solution Vertice Not in Main Cycle: OK" << endl;
        } else {
            cout << "Solution Vertice Not in Main Cycle: Error - ";
            errorFlag = 1;
            for (int i = 0; i < (int) solutionVerticesNotInMainCycle.size() - 1; i++) {
                cout << solutionVerticesNotInMainCycle[i] << " ";
            }
            cout << solutionVerticesNotInMainCycle[solutionVerticesNotInMainCycle.size() - 1] << endl;
        }

        if (nonSolutionVerticesInMainCycle.size() <= 0) {
            cout << "Non Solution Vertice in Main Cycle: OK" << endl;
        } else {
            cout << "Non Solution Vertice in Main Cycle: Error - ";
            errorFlag = 1;
            for (int i = 0; i < (int) nonSolutionVerticesInMainCycle.size() - 1; i++) {
                cout << nonSolutionVerticesInMainCycle[i] << " ";
            }
            cout << nonSolutionVerticesInMainCycle[nonSolutionVerticesInMainCycle.size() - 1] << endl;
        }

        ofstream errorFlagFile ("../" + path + "/output/N" + N + "D" + D + "K" + K + "T" + T + "P" + P + "/errorFlag.txt", , ofstream::out);
        errorFlagFile << errorFlag;
        errorFlagFile.close();
    }
    return 0;
}
