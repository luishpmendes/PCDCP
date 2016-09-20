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
#include <cmath>

#ifndef INFINITE
#define INFINITE 15 << 25
#endif

using namespace std;

typedef long int ulint;
typedef vector < vector <ulint> > matrix;

string itos (ulint i) {
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

int main () {
    ulint n, mComplete, m, k, t, root;
    double d, p;

    cin >> n >> d >> k >> t >> p >> mComplete >> m >> root;

    vector <ulint> penalty (n); // vector with the penalties of each vertex
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

        env.set(GRB_IntParam_LazyConstraints, 1);
        env.set(GRB_IntParam_LogToConsole, 0);
        env.set(GRB_StringParam_LogFile, "./output/N" + N + "D" + D + "K" + K + "T" + T + "P" + P + "/log1.txt");
        env.set(GRB_DoubleParam_TimeLimit, 10);

        GRBModel model = GRBModel(env);

        model.getEnv().set(GRB_IntParam_LazyConstraints, 1);
        model.getEnv().set(GRB_IntParam_LogToConsole, 0);
        model.getEnv().set(GRB_StringParam_LogFile, "./output/N" + N + "D" + D + "K" + K + "T" + T + "P" + P + "/log1.txt");
        model.getEnv().set(GRB_DoubleParam_TimeLimit, 10);

        vector <GRBVar> y (n);

        // ∀ v ∈ V
        for (ulint v = 0; v < n; v++) {
            // y_v ∈ {0.0, 1.0}
            y[v] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "y_" + itos(v));
        }

        model.update();

        GRBLinExpr obj = 0.0;

        // obj = ∑ πv * xv
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

        model.optimize();

        if (model.get(GRB_IntAttr_SolCount) > 0) {
            vector <ulint> solutionV;
            for (ulint v = 0; v < n; v++) {
                if (y[v].get(GRB_DoubleAttr_X) > 0.5) {
                    solutionV.push_back(v);
                }
            }
            vector < pair < pair <ulint, ulint> , ulint> > solutionE;
            for (ulint e = 0; e < m; e++) {
                if (y[E[e].first.first].get(GRB_DoubleAttr_X) > 0.5 && y[E[e].first.second].get(GRB_DoubleAttr_X) > 0.5) {
                    solutionE.push_back(E[e]);
                }
            }

            cout << n << ' ' << d << ' ' << k << ' ' << t << ' ' << p << ' ' << solutionV.size() << ' ' << solutionE.size() << ' ' << root;

            for (ulint i = 0; i < n; i++) {
                cout << penalty[i] << endl;
            }

            for (vector <ulint> :: iterator it = solutionV.begin(); it != solutionV.end(); it++) {
                int v = *it;
                cout << v << endl;
            }

            for (vector < pair < pair <ulint, ulint> , ulint> > :: iterator it = solutionE.begin(); it != solutionE.end(); it++) {
                pair < pair <ulint, ulint> , ulint> e = *it;
                ulint u, v, w;
                u = e.first.first;
                v = e.first.second;
                w = e.second;
                cout << u << ' ' << v << ' ' << w << endl;
            }
        }

        // exporting model
        model.write("./output/N" + N + "D" + D + "K" + K + "T" + T + "P" + P + "/model1.lp");
    } catch (GRBException e) {
        cout << "Error code = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    } catch (...) {
        cout << "Exception during opstimisation" << endl;
    }
    return 0;
}
