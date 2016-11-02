#include <iostream>
#include <vector>
#include <list>
#include <queue>
#include <set>
#include <sstream>
#include <map>
#include <string>
#include <algorithm>
#include <iomanip>
#include <cmath>

#ifndef INFINITE
#define INFINITE 15 << 25
#endif

#ifndef NIL
#define NIL - (15 << 25)
#endif

#ifndef WHITE
#define WHITE 0
#endif

#ifndef GRAY
#define GRAY  1
#endif

#ifndef BLACK
#define BLACK 2
#endif

using namespace std;

typedef long int ulint;
typedef vector < vector <ulint> > matrix;

string itos (ulint i) {
    stringstream s;
    s << i;
    return s.str();
}

string ftos (double i) {
    stringstream s;
    s << i;
    return s.str();
}
vector < set <ulint> > neighbourhoods (matrix W, ulint k) {
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

int main (int argc, char * argv[]) {
    ulint maxIterations;
    double alpha;

    if (argc == 3) {
        maxIterations = atoi(argv[1]);
        alpha = atof(argv[2]);
    } else {
        maxIterations = 100;
        alpha = 0.5;
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
    matrix W (n, vector <ulint> (n, INFINITE)); // adjacency matrix for the complete graph
    vector < list < pair <ulint, ulint> > > adj (n); // adjacency lists for the graph

    for (ulint i = 0; i < n; i++) {
        W[i][i] = 0;
    }

    for (ulint i = 0; i < n; i++) {
        double x, y;
        cin >> x >> y >> penalty[i];
    }

    vector < pair < pair <ulint, ulint> , ulint> > E (m); // vector of edges with the format ((u, v), w)
    map < pair <ulint, ulint>, ulint> mE; // map an edge to its ID

    // reading complete graph
    for (ulint e = 0; e < mComplete; e++) {
        ulint u, v, w;
        cin >> u >> v >> w;
        W[u][v] = w;
        W[v][u] = w;
    }

    // reading graph
    for (ulint e = 0; e < m; e++) {
        ulint u, v, w;
        cin >> u >> v >> w;
        adj[u].push_back(make_pair(v, w));
        adj[v].push_back(make_pair(u, w));
        E[e] = make_pair(make_pair(u, v), w);
        mE[make_pair(u, v)] = e;
        mE[make_pair(v, u)] = e;
    }

    vector < set <ulint> > Ns = neighbourhoods (W, k);

    return 0;
}
