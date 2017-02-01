#include <iostream>
#include <vector>
#include <sstream>
#include <set>
#include <list>
#include <chrono>
#include <algorithm>
#include <iomanip>
#include <fstream>

using namespace std;

typedef long int lint;
typedef unsigned long int ulint;
typedef vector < vector <lint> > matrix;

string itos (ulint i) {
    stringstream s;
    s << i;
    return s.str();
}

vector < set <ulint> > neighbourhoods (matrix W, lint k) {
    vector < set <ulint> > result (W.size());
    for (ulint u = 0; u < (ulint) W.size(); u++) {
        for (ulint v = 0; v < (ulint) W[u].size(); v++) {
            if (W[u][v] <= k) {
                result[u].insert(v);
            }
        }
    }
    return result;
}

void floydWarshall (matrix W, matrix * D, matrix * PI) {
    *D = matrix (W.size(), vector <lint> (W.size(), -1));
    *PI = matrix (W.size(), vector <lint> (W.size(), -1));
    for (ulint i = 0; i < W.size(); i++) {
        for (ulint j = 0; j < W.size(); j++) {
            (*D)[i][j] = W[i][j];
            if (i != j && W[i][j] >= 0) {
                (*PI)[i][j] = i;
            }
        }
    }
    for (ulint k = 0; k < W.size(); k++) {
        for (ulint i = 0; i < W.size(); i++) {
            for (ulint j = 0; j < W.size(); j++) {
                if ((*D)[i][k] >= 0 && (*D)[k][j] >= 0) {
                    if ((*D)[i][j] < 0 || (*D)[i][j] > (*D)[i][k] + (*D)[k][j]) {
                        (*D)[i][j] = (*D)[i][k] + (*D)[k][j];
                        (*PI)[i][j] = (*PI)[k][j];
                    }
                }
            }
        }
    }
}

void geneticAlgorithm (matrix W, vector < list < pair <ulint, ulint> > > adj, matrix Dist, matrix PI, vector <ulint> penalty, ulint root, vector < set <ulint> > Ns, ulint maxIterations, double alpha, ulint seed, vector <ulint> * solution, ulint * solutionCost) {
}

int main (int argc, char * argv[]) {
    chrono :: steady_clock :: time_point tBegin = chrono :: steady_clock :: now();
    string I ("0");
    ulint maxIterations = 100;
    double alpha = 0.3;

    if (argc >= 2) {
        I = string (argv[1]);
    }

    if (argc >= 3) {
        maxIterations = atoi(argv[2]);
    }

    if (argc >= 4) {
        alpha = atof(argv[3]);
    }

    if (alpha < 0.0) {
        alpha = 0.0;
    } else if (alpha > 1.0) {
        alpha = 1.0;
    }

    ulint n, mComplete, m, k, t, root;
    double d, p;

    cin >> n >> d >> k >> t >> p >> mComplete >> m >> root;

    vector <ulint> penalty (n); // vector with de penalties of each vectex
    matrix WComplete (n, vector <lint> (n, -1)); // adjacency matrix for the complete graph
    matrix W (n, vector <lint> (n, -1)); // adjacency matrix for the graph
    vector < list < pair <ulint, ulint> > > adj (n); // adjacency lists for the graph
    matrix Dist, PI; // matrices for the distance and precedence on the graph

    for (ulint v = 0; v < n; v++) {
        WComplete[v][v] = 0;
        W[v][v] = 0;
    }

    for (ulint v = 0; v < n; v++) {
        double x, y;
        cin >> x >> y >> penalty[v];
    }

    // reading complete graph
    for (ulint e = 0; e < mComplete; e++) {
        ulint u, v, w;
        cin >> u >> v >> w;
        WComplete[u][v] = w;
        WComplete[v][u] = w;
    }

    // reading graph
    for (ulint e = 0; e < m; e++) {
        ulint u, v, w;
        cin >> u >> v >> w;
        W[u][v] = w;
        W[v][u] = w;
        adj[u].push_back(make_pair(v, w));
        adj[v].push_back(make_pair(u, w));
    }

    vector < set <ulint> > Ns = neighbourhoods (WComplete, k);

    floydWarshall (W, &Dist, &PI);

    ulint seed = chrono::system_clock::now().time_since_epoch().count();

    vector <ulint> solution;
    ulint solutionCost = 0;

    geneticAlgorithm(W, adj, Dist, PI, penalty, root, Ns, maxIterations, alpha, seed, &solution, &solutionCost);

    cout << solution.size() << ' ' << solution.size() << ' ' << solutionCost << endl;

    for (vector <ulint> :: iterator it = solution.begin(); it != solution.end(); it++) {
        ulint v = *it;
        cout << v << endl;
    }

    for (ulint i = 0; i < solution.size() - 1; i++) {
        pair <ulint, ulint> e = minmax(solution[i], solution[i + 1]);
        cout << e.first << ' ' << e.second << endl;
    }
    pair <ulint, ulint> e = minmax(solution[solution.size() - 1], solution[0]);
    cout << e.first << ' ' << e.second << endl;

    string N = itos(n);
    stringstream ssD;
    ssD << fixed << setprecision(1) << d;
    string D = ssD.str();
    D.erase(remove(D.begin(), D.end(), '.'), D.end());
    string K = itos(k);
    string T = itos(t);
    stringstream ssP;
    ssP << fixed << setprecision(1) << p;
    string P = ssP.str();
    P.erase(remove(P.begin(), P.end(), '.'), P.end());
    stringstream ssA;
    ssA << fixed << setprecision(1) << alpha;
    string A = ssA.str();
    A.erase(remove(A.begin(), A.end(), '.'), A.end());

    ofstream objValFile ("./output/N" + N + "D" + D + "K" + K + "T" + T + "P" + P + "I" + I + "A" + A + "/objVal.txt", ofstream :: out);
    objValFile << solutionCost;
    objValFile.close();

    chrono :: steady_clock :: time_point tEnd = chrono :: steady_clock :: now();
    chrono :: nanoseconds elapsedTime = chrono :: duration_cast <chrono :: nanoseconds> (tEnd - tBegin);
    ofstream elapsedTimeFile ("./output/N" + N + "D" + D + "K" + K + "T" + T + "P" + P + "I" + I + "A" + A + "/elapsedTime.txt", ofstream :: out);
    elapsedTimeFile << elapsedTime.count();
    elapsedTimeFile.close();

    return 0;
}