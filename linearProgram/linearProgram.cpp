#include "gurobi_c++.h"
#include <iostream>
#include <vector>
#include <list>
#include <queue>
#include <set>
#include <sstream>
#include <map>
#include <climits>
#include <string>
#include <algorithm>
#include <iomanip>

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

string itos(ulint i) {
    stringstream s;
    s << i;
    return s.str();
}

string ftos(double i) {
    stringstream s;
    s << i;
    return s.str();
}
/*
class myCompare {
    public:
        bool operator ()(pair <int, int> a, pair <int, int> b) {
            return (a.second > b.second || (a.second == b.second && a.first > b.first));
        }
};
*/
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

void bfs (vector < list < pair <ulint, ulint> > > adj, map < pair <ulint, ulint>, ulint> mE, ulint n, ulint m, ulint s, vector <ulint> * vertex, vector <ulint> * edge, ulint * cycleID) {
    vector <int> color (adj.size(), WHITE);
    queue <ulint> Q;

    *vertex = vector <ulint> (n, 0);
    *edge = vector <ulint> (m, 0);

    color[s] = GRAY;

    Q.push(s);

    while (!Q.empty()) {
        ulint u = Q.front();
        Q.pop();

        for (list < pair <ulint, ulint> >::iterator it = adj[u].begin(); it != adj[u].end(); ++it) {
            ulint v = (*it).first;
            if (color[v] == WHITE) {
                color[v] = GRAY;
                ulint e = mE[make_pair(u, v)];
                (*edge)[e] = * cycleID;
                Q.push(v);
            }
        }
        color[u] = BLACK;
    }
    for (ulint v = 0; v < n; v++) {
        if (color[v] == BLACK) {
            (*vertex)[v] = * cycleID;
        }
    }
    (*cycleID) = (*cycleID) + 1;
}

class subtourelim: public GRBCallback {
    public :
        vector <GRBVar> y;
        vector <GRBVar> x;
        ulint n;
        ulint m;
        vector < pair < pair <ulint, ulint> , ulint> > E;
        map < pair <ulint, ulint>, ulint> mE;
        ulint u;

        subtourelim (vector <GRBVar> _y, vector <GRBVar> _x, ulint _n, ulint _m, vector < pair < pair <ulint, ulint>, ulint> > _E, map < pair <ulint, ulint>, ulint> _mE, ulint _u) {
            y = _y;
            x = _x;
            n = _n;
            m = _m;
            E = _E;
            mE = _mE;
            u = _u;
        }
    protected :
        // let C be the cycle that contains u, with n1 vertices
        // let S = n - n1
        // the sum of the edges that are not in C must be less than S-1
        void callback () {
            try {
                if (where == GRB_CB_MIPSOL) {
                    // find subtour containing vertex u, and add constraint to edges not in the tour
                    vector < list < pair <ulint, ulint> > > adj (n);

                    for (ulint e = 0; e < m; e++) {
                        if (getSolution(x[e]) > 0) {
                            ulint a = E[e].first.first;
                            ulint b = E[e].first.second;
                            ulint c = E[e].second;
                            adj[a].push_back(make_pair(b, c));
                            adj[b].push_back(make_pair(a, c));
                        }
                    }

                    // bfs from u, recording which vertex and edges are visited
                    vector <ulint> vertexCycle (n, -1);
                    vector <ulint> edgeCycle (m, -1);

                    ulint numCycles = 0;
                    while (true) {
                        ulint v = 0;
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

                    vector < set <ulint> > cyclesVertices (numCycles);
                    vector < set <ulint> > cyclesEdges (numCycles);

                    for (ulint v = 0; v < n; v++) {
                        if (vertexCycle[v] >= 0) {
                            cyclesVertices[vertexCycle[v]].insert(v);
                        }
                    }

                    for (ulint e = 0; e < m; e++) {
                        if (edgeCycle[e] >= 0) {
                            cyclesEdges[edgeCycle[e]].insert(e);
                        }
                    }
                    for (ulint c = 1; c < numCycles; c++) {
                        GRBLinExpr expr = 0;
                        for (set <ulint> ::iterator it = cyclesEdges[c].begin(); it != cyclesEdges[c].end(); it++) {
                            ulint e = *it;
                            expr += x[e];
                        }
                        addLazy(expr <= cyclesVertices[c].size() - 1);
                    }
                }
            } catch (GRBException e) {
                cout << "Error number: " << e.getErrorCode() << endl;
                cout << "Error message: " << e.getMessage() << endl;
            } catch (...) {
                cout << "Error during callback" << endl;
            }
        }
};

int main () {
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

    ulint solutionCost = 0;
    set <ulint> solutionVectices;
    set <ulint> solutionEdges;

    try {
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

        GRBEnv env = GRBEnv();

        env.set(GRB_IntParam_LogToConsole, 0);
        env.set(GRB_StringParam_LogFile, "./output/N" + N + "D" + D + "K" + K + "T" + T + "P" + P + "/log.txt");

        GRBModel model = GRBModel(env);

        vector <GRBVar> y (n);

        // ∀ v ∈ V
        for (ulint v = 0; v < n; v++) {
            // y_v ∈ {0.0, 1.0}
            y[v] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "y_" + itos(v));
        }

        vector <GRBVar> x (m);

        // ∀ (u, v) ∈ E
        for (ulint e = 0; e < m; e++) {
            ulint u, v;
            u = E[e].first.first;
            v = E[e].first.second;
            // y_u_v ∈ {0.0, 1.0}
            x[e] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "x_" + itos(u) + "_" + itos(v));
        }

        model.update();

        GRBLinExpr obj = 0.0;

        // obj = ∑ ce * ye
        for (ulint e = 0; e < m; e++) {
            ulint w;
            w = E[e].second;
            obj += w * x[e];
        }

        // obj += ∑ πv * xv
        for (ulint v = 0; v < n; v++) {
            obj += penalty[v] * (1.0 - y[v]);
        }

        model.setObjective(obj, GRB_MINIMIZE);

        // yu == 1
        model.addConstr(y[root] == 1.0, "c_0");

        // dominance
        // ∀ v ∈ V
        for (ulint v = 0; v < n; v++) {
            // ∑ yw >= 1 , w ∈ Nk(v)
            GRBLinExpr constr = 0.0;
            for (set <ulint> :: iterator it = Ns[v].begin(); it != Ns[v].end(); it++) {
                ulint w = *it;
                constr += y[w];
            }
            model.addConstr(constr >= 1.0, "c_1_" + itos(v));
        }

        // each vertex must have exactly two edges adjacent to itself
        // ∀ v ∈ V
        for (ulint v = 0; v < n; v++) {
            // ∑ xe == 2 * yv , e ∈ δ({v})
            GRBLinExpr constr = 0.0;
            for (list < pair <ulint, ulint> >::iterator it = adj[v].begin(); it != adj[v].end(); it++) {
                ulint w = (*it).first; // destination
                ulint e = mE[make_pair(v, w)];
                constr += x[e];
            }
            model.addConstr(constr == 2.0 * y[v], "c_2_" + itos(v));
        }

        subtourelim cb = subtourelim(y, x, n, m, E, mE, root);
        model.setCallback(&cb);

        model.optimize();

        if (model.get(GRB_IntAttr_SolCount) > 0) {
            solutionCost = model.get(GRB_DoubleAttr_ObjVal);
            solutionVectices.clear();
            for (ulint v = 0; v < n; v++) {
                if (y[v].get(GRB_DoubleAttr_X) == 1) {
                    solutionVectices.insert(v);
                }
            }
            solutionEdges.clear();
            for (ulint e = 0; e < m; e++) {
                if (x[e].get(GRB_DoubleAttr_X) == 1) {
                    solutionEdges.insert(e);
                }
            }
        } else {
            cout << "Solution not found." << endl;
        }

        // exporting model
        model.write("./output/N" + N + "D" + D + "K" + K + "T" + T + "P" + P + "/model.lp");
    } catch (GRBException e) {
        cout << "Error code = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    } catch (...) {
        cout << "Exception during opstimisation" << endl;
    }

    cout << solutionVectices.size() << ' ' << solutionEdges.size() << ' ' << solutionCost << endl;
    for (set <ulint> :: iterator it = solutionVectices.begin(); it != solutionVectices.end(); it++) {
        ulint v = *it;
        cout << v << endl;
    }
    for (set <ulint> :: iterator it = solutionEdges.begin(); it != solutionEdges.end(); it++) {
        ulint e = *it;
        ulint a = E[e].first.first;
        ulint b = E[e].first.second;
        cout << a << " " << b << endl;
    }
    return 0;
}
