#include "gurobi_c++.h"
#include <iostream>
#include <vector>
#include <list>
#include <queue>
#include <set>

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

class myCompare {
    public:
        bool operator ()(pair <int, int> a, pair <int, int> b) {
            return (a.second > b.second || (a.second == b.second && a.first > b.first));
        }
};

set <int> neighbourhood (vector < list < pair <int, int> > > adj, int s, int k) {
    set <int> result;
    vector <int> d (adj.size(), INFINITE);
    vector <bool> Q (adj.size(), true);

    d[s] = 0;

    for (int i = 0; i < (int) adj.size(); i++) {
        int u = -1;
        for (int j = 0; j < (int) adj.size(); j++) {
            if (Q[j] && (u < 0 || d[u] > d[j])) {
                u = j;
            }
        }

        Q[u] = false;

        for (list < pair <int, int> >::iterator it = adj[u].begin(); it != adj[u].end(); ++it) {
            int v = (*it).first;
            int w = (*it).second;
            if (d[v] > d[u] + w) {
                d[v] = d[u] + w;
            }
        }
    }

    for (int u = 0; u < (int) adj.size(); u++) {
        if (d[u] <= k) {
            result.insert(u);
        }
    }
    return result;
}

int main () {
    int n, m, k;
    cin >> n >> m >> k;

    vector <int> penalty (n);
    vector < list < pair <int, int> > > adj (n);

    for (int i = 0; i < n; i++) {
        cin >> penalty[i];
    }

    vector < pair < pair <int, int> , int> > E (m);

    for (int i = 0; i < m; i++) {
        int u, v, w;
        cin >> u >> v >> w;
        adj[u].push_back(make_pair(v, w));
        adj[v].push_back(make_pair(u, w));
        E[i] = make_pair(make_pair(u, v), w);
    }

    for (int u = 0; u < n; u++) {
        try {
            GRBEnv env = GRBEnv();

            GRBModel model = GRBModel(env);

            vector <GRBVar> y (n);

            for (int v = 0; v < n; v++) {
                y[v] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "y_" + itos(v));
            }

            vector <GRBVar> x (m);

            for (int e = 0; e < m; e++) {
                int u, v;
                u = E[e].first.first;
                v = E[e].first.second;
                x[e] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "x_" + itos(u) + "_" + itos(v));
            }

            model.update();

            GRBLinExpr obj = 0.0;

            for (int e = 0; e < m; e++) {
                int w;
                w = E[e].second;
                obj += w * x[e];
            }

            for (int v = 0; v < n; v++) {
                obj += p[v] * (1.0 - y[v]);
            }

            model.setObjective(obj, GRB_MINIMIZE);

            model.addConstr(y[u] == 1, "c_0");


            // lazy constraint
        } catch (GRBException e) {
            cout << "Error code = " << e.getErrorCode() << endl;
            cout << e.getMessage() << endl;
        } catch (...) {
            cout << "Exception during opstimisation" << endl;
        }
    }
}
