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

using namespace std;

class myCompare {
    public:
        bool operator ()(pair <int, int> a, pair <int, int> b) {
            return (a.second > b.second || (a.second == b.second && a.first > b.first));
        }
};

set <int> neighbourhood (vector < list < pair <int, int> > > adj, int s, int k) {
    set <int> result;
    vector <int> color (adj.size(), WHITE);
    priority_queue < pair <int, int>, vector < pair <int, int> >,  myCompare> Q;

    vector <int> d (adj.size(), INFINITE);

    for (int u = 0; u < (int) adj.size(); u++) {
        color[u] = WHITE;
        d[u] = INFINITE;
    }

    color[s] = GRAY;
    d[s] = 0;

    Q.push(s);

    while (!Q.empty()) {
        int u = Q.top();
        Q.pop();

        for (list < pair <int, int> >::iterator it = adj[u].begin(); it != adj[u].end(); ++it) {
            int v = (*it).first;
            int w = (*it).second;
            if (color[v] == WHITE) {
                d[v] = d[u] + w;
                if (d[v] <= k) {
                    color[v] = GRAY;
                    Q.push(v);
                } else {
                    color[v] = BLACK;
                }
            }
        }
        color[u] = BLACK;
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
