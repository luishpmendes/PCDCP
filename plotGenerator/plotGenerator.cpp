#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

using namespace std;

typedef unsigned long int ulint;

bool pointComp (pair < pair <double, double>, double > lhs, pair < pair <double, double>, double > rhs) {
    return lhs.second < rhs.second;
}

int main (int argc, char * argv[]) {
    string path, N, D, K, T, P, I, A;

    if (argc == 8) {
        path = string(argv[1]);
        N = string(argv[2]);
        D = string(argv[3]);
        K = string(argv[4]);
        T = string(argv[5]);
        P = string(argv[6]);
        I = string(argv[7]);
    } else if (argc == 9) {
        path = string(argv[1]);
        N = string(argv[2]);
        D = string(argv[3]);
        K = string(argv[4]);
        T = string(argv[5]);
        P = string(argv[6]);
        I = string(argv[7]);
        A = string(argv[8]);
    } else {
        cin >> path >> N >> D >> K >> T >> P >> I >> A;
    }

    ifstream inputFile ("../input/instanceN" + N + "D" + D + "K" + K + "T" + T + "P" + P + "I" + I + ".in");
    ifstream resultFile;

    if (path.compare("grasp") != 0) {
        resultFile = ifstream("../" + path + "/output/N" + N + "D" + D + "K" + K + "T" + T + "P" + P + "I" + I + "/result.out");
    } else {
        resultFile = ifstream("../" + path + "/output/N" + N + "D" + D + "K" + K + "T" + T + "P" + P + "I" + I + "A" + A + "/result.out");
    }

    if (inputFile.is_open() && resultFile.is_open()) {
        ulint n, mComplete, m, k, t, root;
        double d, p;
        inputFile >> n >> d >> k >> t >> p >> mComplete >> m >> root;
        vector < pair < pair <double, double>, double > > vertices (n, make_pair(make_pair(0, 0), 0));
        // reading vertices' coordinates and penalty
        for (ulint v = 0; v < n; v++) {
            inputFile >> vertices[v].first.first;
            inputFile >> vertices[v].first.second;
            inputFile >> vertices[v].second;
        }
        // limiting penalities to the interval [0.25, 1]
        double minPenalty = vertices[0].second;
        double maxPenalty = vertices[0].second;
        for (ulint v = 1; v < n; v++) {
            if (minPenalty > vertices[v].second) {
                minPenalty = vertices[v].second;
            }
            if (maxPenalty < vertices[v].second) {
                maxPenalty = vertices[v].second;
            }
        }
        for (ulint v = 0; v < n; v++) {
            vertices[v].second -= minPenalty;
        }
        minPenalty -= minPenalty;
        maxPenalty -= minPenalty;
        for (ulint v = 0; v < n; v++) {
            vertices[v].second *= (1-0.25)/maxPenalty;
        }
        minPenalty *= (1-0.25)/maxPenalty;
        maxPenalty *= (1-0.25)/maxPenalty;
        for (ulint v = 0; v < n; v++) {
            vertices[v].second += 0.25;
        }

        // printing vertices' coordinates and penalty
        ofstream verticesFile;
        if (path.compare("grasp") != 0) {
            verticesFile = ofstream("../" + path + "/output/N" + N + "D" + D + "K" + K + "T" + T + "P" + P + "I" + I + "/vertices.txt");
        } else {
            verticesFile = ofstream("../" + path + "/output/N" + N + "D" + D + "K" + K + "T" + T + "P" + P + "I" + I + "A" + A + "/vertices.txt");
        }
        for (ulint v = 0; v < n; v++) {
            verticesFile << vertices[v].first.first << ' ';
            verticesFile << vertices[v].first.second << ' ';
            verticesFile << vertices[v].second << endl;
        }
        verticesFile.close();
        
        // reading (and ignoring) the complete graph's edges
        for (ulint j = 0; j < mComplete; j++) {
            ulint u, v, w;
            inputFile >> u >> v >> w;
        }

        // reading the graph's edges and printing its coordinates
        ofstream edgesFile;
        if (path.compare("grasp") != 0) {
            edgesFile = ofstream("../" + path + "/output/N" + N + "D" + D + "K" + K + "T" + T + "P" + P + "I" + I + "/edges.txt");
        } else {
            edgesFile = ofstream("../" + path + "/output/N" + N + "D" + D + "K" + K + "T" + T + "P" + P + "I" + I + "A" + A + "/edges.txt");
        }
        for (ulint e = 0; e < m; e++) {
            ulint u, v, w;
            inputFile >> u >> v >> w;
            edgesFile << vertices[u].first.first << ' ';
            edgesFile << vertices[u].first.second << endl;
            edgesFile << vertices[v].first.first << ' ';
            edgesFile << vertices[v].first.second << endl << endl;
        }
        edgesFile.close();

        ulint nSolution, mSolution;
        double costSolution;
        resultFile >> nSolution >> mSolution >> costSolution;
        // reading the solution's vertices and printing its coordinates and neighborhood radio
        ofstream solutionVerticesFile;
        if (path.compare("grasp") != 0) {
            solutionVerticesFile = ofstream("../" + path + "/output/N" + N + "D" + D + "K" + K + "T" + T + "P" + P + "I" + I + "/solutionVertices.txt");
        } else {
            solutionVerticesFile = ofstream("../" + path + "/output/N" + N + "D" + D + "K" + K + "T" + T + "P" + P + "I" + I + "A" + A + "/solutionVertices.txt");
        }
        for (ulint i = 0; i < nSolution; i++) {
            ulint v;
            resultFile >> v;
            solutionVerticesFile << vertices[v].first.first << ' ';
            solutionVerticesFile << vertices[v].first.second << ' ';
            solutionVerticesFile << ((double) k)/100.0 << endl;
        }
        solutionVerticesFile.close();

        // reading the solution's edges and printing its coordinates
        ofstream solutionEdgesFile;
        if (path.compare("grasp") != 0) {
            solutionEdgesFile = ofstream("../" + path + "/output/N" + N + "D" + D + "K" + K + "T" + T + "P" + P + "I" + I + "/solutionEdges.txt");
        } else {
            solutionEdgesFile = ofstream("../" + path + "/output/N" + N + "D" + D + "K" + K + "T" + T + "P" + P + "I" + I + "A" + A + "/solutionEdges.txt");
        }
        for (ulint e = 0; e < mSolution; e++) {
            ulint u, v;
            resultFile >> u >> v;
            solutionEdgesFile << vertices[u].first.first << ' ';
            solutionEdgesFile << vertices[u].first.second << endl;
            solutionEdgesFile << vertices[v].first.first << ' ';
            solutionEdgesFile << vertices[v].first.second << endl << endl;
        }
        solutionEdgesFile.close();
    }
    return 0;
}
