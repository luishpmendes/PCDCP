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
        vector < pair < pair <double, double>, double > > points (n, make_pair(make_pair(0, 0), 0));
        cout << "Vertices: " << endl;
        for (ulint v = 0; v < n; v++) {
            inputFile >> points[v].first.first;
            inputFile >> points[v].first.second;
            inputFile >> points[v].second;
            cout << v << " : (" << points[v].first.first << ", " << points[v].first.second << ")" << endl;
        }
        vector < pair < pair <double, double>, double > > sortedPoints (points.begin(), points.end());
        sort(sortedPoints.begin(), sortedPoints.end(), pointComp);
        double dw = (1.0 - 0.25) / ((double) (n - 1));
        sortedPoints[0].second = 0.25;
        for (int v = 1; v < n; v++) {
            sortedPoints[v].second = sortedPoints[v - 1].second + dw;
        }
        ofstream pointsFile ("../" + path + "/output/N" + N + "D" + D + "K" + K + "T" + T + "/points.txt");
        for (ulint v = 0; v < n; v++) {
            pointsFile << sortedPoints[v].first.first << ' ';
            pointsFile << sortedPoints[v].first.second << ' ';
            pointsFile << sortedPoints[v].second << endl;
        }
        pointsFile.close();
        ulint nSolution, mSolution;
        double costSolution;
        resultFile >> nSolution >> mSolution >> costSolution;
        ofstream domSetFile ("../" + path + "/output/N" + N + "D" + D + "K" + K + "T" + T + "/dom_set.txt");
        for (ulint i = 0; i < nSolution; i++) {
            int v;
            resultFile >> v;
            domSetFile << points[i].first.first << ' ';
            domSetFile << points[i].first.second << ' ';
            domSetFile << ((double) k)/100.0 << endl;
        }
        vector < pair < pair <double, double>, pair <double, double> > > edgesCoordinates (mSolution, make_pair(make_pair(0.0, 0.0), make_pair(0.0, 0.0)));
        cout << "Edges:" << endl;
        for (ulint e = 0; e < mSolution; e++) {
            ulint u, v;
            resultFile >> u >> v;
            edgesCoordinates[e].first.first = points[u].first.first;
            edgesCoordinates[e].first.second = points[u].first.second;
            edgesCoordinates[e].second.first = points[v].first.first;
            edgesCoordinates[e].second.second = points[v].first.second;
            cout << u << " <-> " << v << " : (" << points[u].first.first << ", " << points[u].first.second << ") <-> (" << points[v].first.first << ", " << points[v].first.second << ")" << endl;
        }
        ofstream edgesFile ("../" + path + "/output/N" + N + "D" + D + "K" + K + "T" + T + "/edges.txt");
        for (ulint e = 0; e < mSolution; e++) {
            edgesFile << edgesCoordinates[e].first.first << ' ';
            edgesFile << edgesCoordinates[e].first.second << endl;
            edgesFile << edgesCoordinates[e].second.first << ' ';
            edgesFile << edgesCoordinates[e].second.second << endl << endl;
        }
        edgesFile.close();
    }
    return 0;
}
