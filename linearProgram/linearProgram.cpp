#include "gurobi_c++.h"
#include <iostream>
#include <vector>
#include <list>
#include <queue>
#include <set>
#include <sstream>
#include <map>
#include <climits>

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

typedef vector < vector <int> > matrix;

string itos(int i) {
    stringstream s;
    s << i;
    return s.str();
}

class myCompare {
    public:
        bool operator ()(pair <int, int> a, pair <int, int> b) {
            return (a.second > b.second || (a.second == b.second && a.first > b.first));
        }
};

vector < set <int> > neighbourhoods (matrix W, int k) {
    vector < set <int> > result (W.size());

    matrix D = matrix (W.size(), vector <int> (W.size(), INFINITE));

    for (int i = 0; i < (int) W.size(); i++) {
        for (int j = 0; j < (int) W.size(); j++) {
            D[i][j] = W[i][j];
        }
    }

    for (int l = 0; l < (int) W.size(); l++) {
        for (int i = 0; i < (int) W.size(); i++) {
            for (int j = 0; j < (int) W.size(); j++) {
                if (D[i][j] > D[i][l] + D[l][j]) {
                    D[i][j] = D[i][l] + D[l][j];
                }
            }
        }
    }

    for (int u = 0; u < (int) W.size(); u++) {
        for (int v = 0; v < (int) W[u].size(); v++) {
            if (D[u][v] <= k) {
                result[u].insert(v);
            }
        }
    }

    return result;
}

void bfs (vector < list < pair <int, int> > > adj, map < pair <int, int>, int> mE, int n, int m, int s, vector <int> * vertex, vector <int> * edge, int * cycleID) {
    vector <int> color (adj.size(), WHITE);
    queue <int> Q;

    *vertex = vector <int> (n, 0);
    *edge = vector <int> (m, 0);

    color[s] = GRAY;

    Q.push(s);

    while (!Q.empty()) {
        int u = Q.front();
        Q.pop();

        for (list < pair <int, int> >::iterator it = adj[u].begin(); it != adj[u].end(); ++it) {
            int v = (*it).first;
            if (color[v] == WHITE) {
                color[v] = GRAY;
                int e = mE[make_pair(u, v)];
                (*edge)[e] = * cycleID;
                Q.push(v);
            }
        }
        color[u] = BLACK;
    }
    for (int v = 0; v < n; v++) {
        if (color[v] == BLACK) {
            (*vertex)[v] = * cycleID;
        }
    }
    (*cycleID) = (*cycleID) + 1;
}

class subtourelim: public GRBCallback {
    public:
        vector <GRBVar> y;
        vector <GRBVar> x;
        int n;
        int m;
        vector < pair < pair <int, int> , int> > E;
        map < pair <int, int>, int> mE;
        int u;

        subtourelim (vector <GRBVar> _y, vector <GRBVar> _x, int _n, int _m, vector < pair < pair <int, int> , int> > _E, map < pair <int, int>, int> _mE, int _u) {
            y = _y;
            x = _x;
            n = _n;
            m = _m;
            E = _E;
            mE = _mE;
            u = _u;
        }
    protected:
        // sera C o ciclo que contém u, passando por n1 vértices
        // seja S = n - n1
        // a soma das arestas que não estão em C tem que ser menor que S-1
        void callback () {
            try {
                if (where == GRB_CB_MIPSOL) {
                    // find subtour containing vertex u, and add constraint to edges not in the tour
                    vector < list < pair <int, int> > > adj (n);

                    for (int e = 0; e < m; e++) {
                        if (getSolution(x[e]) > 0) {
                            int a = E[e].first.first;
                            int b = E[e].first.second;
                            int c = E[e].second;
                            adj[a].push_back(make_pair(b, c));
                            adj[b].push_back(make_pair(a, c));
                        }
                    }

                    // bfs from u, recording which vertex and edges are visited
                    vector <int> vertexCycle (n, -1);
                    vector <int> edgeCycle (m, -1);

                    int numCycles = 0;
                    while (true) {
                        int v = 0;
                        for (; v < n; v++) {
                            if (getSolution(y[v]) > 0) { // v is in solution
                                if (vertexCycle[v] < 0) { // cycle containing v not discovered yet
                                    break;
                                }
                            }
                        }
                        if (v >= n) {
                            break;
                        }
                        bfs (adj, mE, n, m, u, &vertexCycle, &edgeCycle, &numCycles);
                    }

                    vector < set <int> > cyclesVertices (numCycles);
                    vector < set <int> > cyclesEdges (numCycles);

                    for (int v = 0; v < n; v++) {
                        if (vertexCycle[v] >= 0) {
                            cyclesVertices[vertexCycle[v]].insert(v);
                        }
                    }

                    for (int e = 0; e < m; e++) {
                        if (edgeCycle[e] >= 0) {
                            cyclesEdges[edgeCycle[e]].insert(e);
                        }
                    }

                    for (int c = 0; c < numCycles; c++) {
                        GRBLinExpr expr = 0;
                        for (set<int>::iterator it = cyclesEdges[c].begin(); it != cyclesEdges[c].end(); it++) {
                            int e = *it;
                            expr += x[e];//???ok?
                        }
                        addLazy(expr <= cyclesVertices[c].size() - 1);
                    }
                }
            } catch (GRBException e) {
                cout << "Error number: " << e.getErrorCode() << endl;
                cout << e.getMessage() << endl;
            } catch (...) {
                cout << "Error during callback" << endl;
            }
        }
};

int main () {
    int n, m, k;
    cin >> n >> m >> k;

    vector <int> penalty (n);
    matrix W (n, vector <int> (n, INFINITE));
    vector < list < pair <int, int> > > adj (n);

    for (int i = 0; i < n; i++) {
        W[i][i] = 0;
    }

    for (int i = 0; i < n; i++) {
        cin >> penalty[i];
    }

    vector < pair < pair <int, int> , int> > E (m); // vector of edges with the format ((u, v), w)
    map < pair <int, int>, int> mE; // map an edge to its ID

    for (int e = 0; e < m; e++) {
        int u, v, w;
        cin >> u >> v >> w;
        W[u][v] = w;
        W[v][u] = w;
        adj[u].push_back(make_pair(v, w));
        adj[v].push_back(make_pair(u, w));
        E[e] = make_pair(make_pair(u, v), w);
        cout << "E[" << e << "] = ((" << u << ", " << v << "), " << w << ")" << endl;
        mE[make_pair(u, v)] = e;
        mE[make_pair(v, u)] = e;
    }

    vector < set <int> > N = neighbourhoods (W, k);

    for (int v = 0; v < n; v++) {
        cout << "N[" << v << "] = {";
        for (set <int>::iterator it = N[v].begin(); it != N[v].end(); it++) {
            int u = *it;
            cout << u << ", ";
        }
        cout << "\b\b}" << endl;
    }

    bool flag = false;
    int bestSolutionCost;
    set <int> bestSolutionVectices;
    set <int> bestSolutionEdges;

    for (int u = 0; u < n; u++) {
        try {
            GRBEnv env = GRBEnv();

            GRBModel model = GRBModel(env);

            vector <GRBVar> y (n);

            // ∀ v ∈ V
            for (int v = 0; v < n; v++) {
                // y_v ∈ {0.0, 1.0}
                y[v] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "y_" + itos(v));
            }

            vector <GRBVar> x (m);

            // ∀ (u, v) ∈ E
            for (int e = 0; e < m; e++) {
                int u, v;
                u = E[e].first.first;
                v = E[e].first.second;
                // y_u_v ∈ {0.0, 1.0}
                x[e] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "x_" + itos(u) + "_" + itos(v));
            }

            model.update();

            GRBLinExpr obj = 0.0;

            // obj = ∑ ce * ye
            for (int e = 0; e < m; e++) {
                int w;
                w = E[e].second;
                obj += w * x[e];
            }

            // obj += ∑ πv * xv
            for (int v = 0; v < n; v++) {
                obj += penalty[v] * (1.0 - y[v]);
            }

            model.setObjective(obj, GRB_MINIMIZE);

            // yu == 1
            model.addConstr(y[u] == 1.0, "c_0");

            // dominance
            // ∀ v ∈ V
            for (int v = 0; v < n; v++) {
                // ∑ yw >= 1 , w ∈ Nk(v)
                GRBLinExpr constr = 0.0;
                for (set <int> :: iterator it = N[v].begin(); it != N[v].end(); it++) {
                    int w = *it;
                    constr += y[w];
                }
                model.addConstr(constr >= 1.0, "c_1_" + itos(v));
            }

            // each vertex must have exactly two edges adjacent to itself
            // ∀ v ∈ V
            for (int v = 0; v < n; v++) {
                // ∑ xe == 2 * yv , e ∈ δ({v})
                GRBLinExpr constr = 0.0;
                for (list < pair <int, int> >::iterator it = adj[v].begin(); it != adj[v].end(); it++) {
                    int w = (*it).first; // destino
                    int e = mE[make_pair(v, w)];
                    constr += x[e];
                }
                model.addConstr(constr == 2.0 * y[v], "c_2_" + itos(v));
            }

            subtourelim cb = subtourelim(y, x, n, m, E, mE, u);
            model.setCallback(&cb);

            model.optimize();

            if (model.get(GRB_IntAttr_SolCount) > 0) {
                if (!flag || model.get(GRB_DoubleAttr_ObjVal) < bestSolutionCost) {
                    flag = true;
                    bestSolutionCost = model.get(GRB_DoubleAttr_ObjVal);
                    bestSolutionVectices.clear();
                    for (int v = 0; v < n; v++) {
                        if (y[v].get(GRB_DoubleAttr_X) == 1) {
                            bestSolutionVectices.insert(v);
                        }
                    }
                    bestSolutionEdges.clear();
                    for (int e = 0; e < m; e++) {
                        if (x[e].get(GRB_DoubleAttr_X) == 1) {
                            bestSolutionEdges.insert(e);
                        }
                    }
                }


            } else {
                cout << "Solution not found." << endl;
            }

        } catch (GRBException e) {
            cout << "Error code = " << e.getErrorCode() << endl;
            cout << e.getMessage() << endl;
        } catch (...) {
            cout << "Exception during opstimisation" << endl;
        }
    }
    cout << "Cost: " << bestSolutionCost << endl;
    cout << "vertex: " << endl;
    for (set <int>::iterator it = bestSolutionVectices.begin(); it != bestSolutionVectices.end(); it++) {
        int v = *it;
        cout << v << endl;
    }
    cout << "edges: " << endl;
    for (set <int>::iterator it = bestSolutionEdges.begin(); it != bestSolutionEdges.end(); it++) {
        int e = *it;
        int a = E[e].first.first;
        int b = E[e].first.second;
        cout << a << " " << b << endl;
    }
}
