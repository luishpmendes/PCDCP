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

using namespace std;

typedef long int lint;
typedef unsigned long int ulint;
typedef vector < vector <lint> > matrix;

string itos(ulint i) {
    stringstream s;
    s << i;
    return s.str();
}

vector < set <ulint> > neighbourhoods (matrix W, ulint k) {
    vector < set <ulint> > result (W.size());
    for (ulint u = 0; u < W.size(); u++) {
        for (ulint v = 0; v < W[u].size(); v++) {
            if (W[u][v] >= 0 && (ulint) W[u][v] <= k) {
                result[u].insert(v);
            }
        }
    }
    return result;
}

int main (int argc, char * argv[]) {
    string path, N, D, K, T, I, A, PS, MR;

    if (argc == 7) {
        path = string(argv[1]);
        N = string(argv[2]);
        D = string(argv[3]);
        K = string(argv[4]);
        T = string(argv[5]);
        I = string(argv[6]);
    } else if (argc == 8) {
        path = string(argv[1]);
        N = string(argv[2]);
        D = string(argv[3]);
        K = string(argv[4]);
        T = string(argv[5]);
        I = string(argv[6]);
        A = string(argv[7]);
    } else if (argc == 9) {
        path = string(argv[1]);
        N = string(argv[2]);
        D = string(argv[3]);
        K = string(argv[4]);
        T = string(argv[5]);
        I = string(argv[6]);
        PS = string(argv[7]);
        MR = string(argv[8]);
    } else {
        cin >> path >> N >> D >> K >> T >> I >> A >> PS >> MR;
    }

    ifstream inputFile ("../input/instanceN" + N + "D" + D + "K" + K + "T" + T + "I" + I + ".in");
    ifstream resultFile;

    if (path.compare("grasp") == 0) {
        resultFile = ifstream("../" + path + "/output/N" + N + "D" + D + "K" + K + "T" + T + "I" + I + "A" + A + "/result.out");
    } else if (path.compare("geneticAlgorithm") == 0) {
        resultFile = ifstream("../" + path + "/output/N" + N + "D" + D + "K" + K + "T" + T + "I" + I + "PS" + PS + "MR" + MR + "/result.out");
    } else {
        resultFile = ifstream("../" + path + "/output/N" + N + "D" + D + "K" + K + "T" + T + "I" + I + "/result.out");
    }

    if (inputFile.is_open() && resultFile.is_open()) {
        ulint n, mComplete, m, k, t, root;
        double d;
        int nFlag, kFlag, nSolutionFlag, mSolutionFlag, solutionVerticesFlag, solutionEdgesFlag, solutionCostFlag, coverFlag, solutionVerticesNotInMainCycleFlag, nonSolutionVerticesInMainCycleFlag, errorFlag;
        nFlag = kFlag = nSolutionFlag = mSolutionFlag = solutionVerticesFlag = solutionEdgesFlag = solutionCostFlag = coverFlag = solutionVerticesNotInMainCycleFlag = nonSolutionVerticesInMainCycleFlag = errorFlag = 0;
        inputFile >> n >> d >> k >> t >> mComplete >> m >> root;

        if (N == itos(n)) {
            cout << "N == n : OK" << endl;
        } else {
            cout << "N == n : ERROR" << endl;
            nFlag = errorFlag = 1;
        }

        if (K == itos(k)) {
            cout << "K == k : OK" << endl;
        } else {
            cout << "K == k : ERROR" << endl;
            kFlag = errorFlag = 1;
        }

        vector < pair < pair <double, double>, ulint > > vertices (n, make_pair(make_pair(0, 0), 0));
        // reading vertices' coordinates and penalty
        for (ulint v = 0; v < n; v++) {
            inputFile >> vertices[v].first.first;
            inputFile >> vertices[v].first.second;
            inputFile >> vertices[v].second;
        }

        matrix W (n, vector <lint> (n, -1)); // adjacency matrix for the complete graph
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
            nSolutionFlag = errorFlag = 1;
        }

        if (mSolution <= m) {
            cout << "mSolution <= m : OK" << endl;
        } else {
            cout << "mSolution <= m : ERROR" << endl;
            mSolutionFlag = errorFlag = 1;
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
            solutionVerticesFlag = errorFlag = 1;
            for (ulint i = 0; i < solutionVerticesInvalid.size(); i++) {
                cout << solutionVerticesInvalid[i] << " ";
            }
            cout << endl;
        }

        set <ulint> setSolutionVertices (solutionVertices.begin(), solutionVertices.end());
        vector < pair <ulint, ulint> > solutionEdges (mSolution, make_pair(0, 0));
        vector < pair <ulint, ulint> > solutionEdgesInvalid;
        // reading the solution's edges
        for (ulint e = 0; e < mSolution; e++) {
            ulint u, v;
            resultFile >> u >> v;
            solutionEdges[e] = minmax(u, v);
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
            solutionEdgesFlag = errorFlag = 1;
            for (ulint i = 0; i < solutionEdgesInvalid.size(); i++) {
                cout << "(" << solutionEdgesInvalid[i].first << ", " << solutionEdgesInvalid[i].second << ") ";
            }
            cout << endl;
        }

        // compute solution cost
        ulint cost = 0;
        for (ulint e = 0; e < mSolution; e++) {
            cost += weights[solutionEdges[e]];
        }

        vector <int> isInSolution(n, 0);

        for (ulint i = 0; i < solutionVertices.size(); i++) {
            isInSolution[solutionVertices[i]] = 1;
        }

        for (ulint v = 0; v < n; v++) {
            if (isInSolution[v] == 0) {
                cost += vertices[v].second;
            }
        }

        // check if solution cost is correct
        if (costSolution == cost) {
            cout << "Solution Cost : OK" << endl;
        } else {
            cout << "Solution Cost : ERROR - " << costSolution << " != " << cost << endl;
            solutionCostFlag = errorFlag = 1;
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
            coverFlag = errorFlag = 1;
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

        for (ulint v = 0; v < n; v++) {
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
            solutionVerticesNotInMainCycleFlag = errorFlag = 1;
            for (int i = 0; i < (int) solutionVerticesNotInMainCycle.size() - 1; i++) {
                cout << solutionVerticesNotInMainCycle[i] << " ";
            }
            cout << solutionVerticesNotInMainCycle[solutionVerticesNotInMainCycle.size() - 1] << endl;
        }

        if (nonSolutionVerticesInMainCycle.size() <= 0) {
            cout << "Non Solution Vertice in Main Cycle: OK" << endl;
        } else {
            cout << "Non Solution Vertice in Main Cycle: Error - ";
            nonSolutionVerticesInMainCycleFlag = errorFlag = 1;
            for (int i = 0; i < (int) nonSolutionVerticesInMainCycle.size() - 1; i++) {
                cout << nonSolutionVerticesInMainCycle[i] << " ";
            }
            cout << nonSolutionVerticesInMainCycle[nonSolutionVerticesInMainCycle.size() - 1] << endl;
        }

        ofstream nFlagFile;
        if (path.compare("grasp") == 0) {
            nFlagFile = ofstream("../" + path + "/output/N" + N + "D" + D + "K" + K + "T" + T + "I" + I + "A" + A + "/nFlag.txt", ofstream::out);
        } else if (path.compare("geneticAlgorithm") == 0) {
            nFlagFile = ofstream("../" + path + "/output/N" + N + "D" + D + "K" + K + "T" + T + "I" + I + "PS" + PS + "MR" + MR + "/nFlag.txt", ofstream::out);
        } else {
            nFlagFile = ofstream("../" + path + "/output/N" + N + "D" + D + "K" + K + "T" + T + "I" + I + "/nFlag.txt", ofstream::out);
        }
        nFlagFile << nFlag;
        nFlagFile.close();

        ofstream kFlagFile;
        if (path.compare("grasp") == 0) {
            kFlagFile = ofstream("../" + path + "/output/N" + N + "D" + D + "K" + K + "T" + T + "I" + I + "A" + A + "/kFlag.txt", ofstream::out);
        } else if (path.compare("geneticAlgorithm") == 0) {
            kFlagFile = ofstream("../" + path + "/output/N" + N + "D" + D + "K" + K + "T" + T + "I" + I + "PS" + PS + "MR" + MR + "/kFlag.txt", ofstream::out);
        } else {
            kFlagFile = ofstream("../" + path + "/output/N" + N + "D" + D + "K" + K + "T" + T + "I" + I + "/kFlag.txt", ofstream::out);
        }
        kFlagFile << kFlag;
        kFlagFile.close();

        ofstream nSolutionFlagFile;
        if (path.compare("grasp") == 0) {
            nSolutionFlagFile = ofstream("../" + path + "/output/N" + N + "D" + D + "K" + K + "T" + T + "I" + I + "A" + A + "/nSolutionFlag.txt", ofstream::out);
        } else if (path.compare("geneticAlgorithm") == 0) {
            nSolutionFlagFile = ofstream("../" + path + "/output/N" + N + "D" + D + "K" + K + "T" + T + "I" + I + "PS" + PS + "MR" + MR + "/nSolutionFlag.txt", ofstream::out);
        } else {
            nSolutionFlagFile = ofstream("../" + path + "/output/N" + N + "D" + D + "K" + K + "T" + T + "I" + I + "/nSolutionFlag.txt", ofstream::out);
        }
        nSolutionFlagFile << nSolutionFlag;
        nSolutionFlagFile.close();

        ofstream mSolutionFlagFile;
        if (path.compare("grasp") == 0) {
            mSolutionFlagFile = ofstream("../" + path + "/output/N" + N + "D" + D + "K" + K + "T" + T + "I" + I + "A" + A + "/mSolutionFlag.txt", ofstream::out);
        } else if (path.compare("geneticAlgorithm") == 0) {
            mSolutionFlagFile = ofstream("../" + path + "/output/N" + N + "D" + D + "K" + K + "T" + T + "I" + I + "PS" + PS + "MR" + MR + "/mSolutionFlag.txt", ofstream::out);
        } else {
            mSolutionFlagFile = ofstream("../" + path + "/output/N" + N + "D" + D + "K" + K + "T" + T + "I" + I + "/mSolutionFlag.txt", ofstream::out);
        }
        mSolutionFlagFile << mSolutionFlag;
        mSolutionFlagFile.close();

        ofstream solutionVerticesFlagFile;
        if (path.compare("grasp") == 0) {
            solutionVerticesFlagFile = ofstream("../" + path + "/output/N" + N + "D" + D + "K" + K + "T" + T + "I" + I + "A" + A + "/solutionVerticesFlag.txt", ofstream::out);
        } else if (path.compare("geneticAlgorithm") == 0) {
            solutionVerticesFlagFile = ofstream("../" + path + "/output/N" + N + "D" + D + "K" + K + "T" + T + "I" + I + "PS" + PS + "MR" + MR + "/solutionVerticesFlag.txt", ofstream::out);
        } else {
            solutionVerticesFlagFile = ofstream("../" + path + "/output/N" + N + "D" + D + "K" + K + "T" + T + "I" + I + "/solutionVerticesFlag.txt", ofstream::out);
        }
        solutionVerticesFlagFile << solutionVerticesFlag;
        solutionVerticesFlagFile.close();

        ofstream solutionEdgesFlagFile;
        if (path.compare("grasp") == 0) {
            solutionEdgesFlagFile = ofstream("../" + path + "/output/N" + N + "D" + D + "K" + K + "T" + T + "I" + I + "A" + A + "/solutionEdgesFlag.txt", ofstream::out);
        } else if (path.compare("geneticAlgorithm") == 0) {
            solutionEdgesFlagFile = ofstream("../" + path + "/output/N" + N + "D" + D + "K" + K + "T" + T + "I" + I + "PS" + PS + "MR" + MR + "/solutionEdgesFlag.txt", ofstream::out);
        } else {
            solutionEdgesFlagFile = ofstream("../" + path + "/output/N" + N + "D" + D + "K" + K + "T" + T + "I" + I + "/solutionEdgesFlag.txt", ofstream::out);
        }
        solutionEdgesFlagFile << solutionEdgesFlag;
        solutionEdgesFlagFile.close();

        ofstream solutionCostFlagFile;
        if (path.compare("grasp") == 0) {
            solutionCostFlagFile = ofstream("../" + path + "/output/N" + N + "D" + D + "K" + K + "T" + T + "I" + I + "A" + A + "/solutionCostFlag.txt", ofstream::out);
        } else if (path.compare("geneticAlgorithm") == 0) {
            solutionCostFlagFile = ofstream("../" + path + "/output/N" + N + "D" + D + "K" + K + "T" + T + "I" + I + "PS" + PS + "MR" + MR + "/solutionCostFlag.txt", ofstream::out);
        } else {
            solutionCostFlagFile = ofstream("../" + path + "/output/N" + N + "D" + D + "K" + K + "T" + T + "I" + I + "/solutionCostFlag.txt", ofstream::out);
        }
        solutionCostFlagFile << solutionCostFlag;
        solutionCostFlagFile.close();

        ofstream coverFlagFile;
        if (path.compare("grasp") == 0) {
            coverFlagFile = ofstream("../" + path + "/output/N" + N + "D" + D + "K" + K + "T" + T + "I" + I + "A" + A + "/coverFlag.txt", ofstream::out);
        } else if (path.compare("geneticAlgorithm") == 0) {
            coverFlagFile = ofstream("../" + path + "/output/N" + N + "D" + D + "K" + K + "T" + T + "I" + I + "PS" + PS + "MR" + MR + "/coverFlag.txt", ofstream::out);
        } else {
            coverFlagFile = ofstream("../" + path + "/output/N" + N + "D" + D + "K" + K + "T" + T + "I" + I + "/coverFlag.txt", ofstream::out);
        }
        coverFlagFile << coverFlag;
        coverFlagFile.close();

        ofstream solutionVerticesNotInMainCycleFlagFile;
        if (path.compare("grasp") == 0) {
            solutionVerticesNotInMainCycleFlagFile = ofstream("../" + path + "/output/N" + N + "D" + D + "K" + K + "T" + T + "I" + I + "A" + A + "/solutionVerticesNotInMainCycleFlag.txt", ofstream::out);
        } else if (path.compare("geneticAlgorithm") == 0) {
            solutionVerticesNotInMainCycleFlagFile = ofstream("../" + path + "/output/N" + N + "D" + D + "K" + K + "T" + T + "I" + I + "PS" + PS + "MR" + MR + "/solutionVerticesNotInMainCycleFlag.txt", ofstream::out);
        } else {
            solutionVerticesNotInMainCycleFlagFile = ofstream("../" + path + "/output/N" + N + "D" + D + "K" + K + "T" + T + "I" + I + "/solutionVerticesNotInMainCycleFlag.txt", ofstream::out);
        }
        solutionVerticesNotInMainCycleFlagFile << solutionVerticesNotInMainCycleFlag;
        solutionVerticesNotInMainCycleFlagFile.close();

        ofstream nonSolutionVerticesInMainCycleFlagFile;
        if (path.compare("grasp") == 0) {
            nonSolutionVerticesInMainCycleFlagFile = ofstream("../" + path + "/output/N" + N + "D" + D + "K" + K + "T" + T + "I" + I + "A" + A + "/nonSolutionVerticesInMainCycleFlag.txt", ofstream::out);
        } else if (path.compare("geneticAlgorithm") == 0) {
            nonSolutionVerticesInMainCycleFlagFile = ofstream("../" + path + "/output/N" + N + "D" + D + "K" + K + "T" + T + "I" + I + "PS" + PS + "MR" + MR + "/nonSolutionVerticesInMainCycleFlag.txt", ofstream::out);
        } else {
            nonSolutionVerticesInMainCycleFlagFile = ofstream("../" + path + "/output/N" + N + "D" + D + "K" + K + "T" + T + "I" + I + "/nonSolutionVerticesInMainCycleFlag.txt", ofstream::out);
        }
        nonSolutionVerticesInMainCycleFlagFile << nonSolutionVerticesInMainCycleFlag;
        nonSolutionVerticesInMainCycleFlagFile.close();

        ofstream errorFlagFile;
        if (path.compare("grasp") == 0) {
            errorFlagFile = ofstream("../" + path + "/output/N" + N + "D" + D + "K" + K + "T" + T + "I" + I + "A" + A + "/errorFlag.txt", ofstream::out);
        } else if (path.compare("geneticAlgorithm") == 0) {
            errorFlagFile = ofstream("../" + path + "/output/N" + N + "D" + D + "K" + K + "T" + T + "I" + I + "PS" + PS + "MR" + MR + "/errorFlag.txt", ofstream::out);
        } else {
            errorFlagFile = ofstream("../" + path + "/output/N" + N + "D" + D + "K" + K + "T" + T + "I" + I + "/errorFlag.txt", ofstream::out);
        }
        errorFlagFile << errorFlag;
        errorFlagFile.close();
    }
    return 0;
}
