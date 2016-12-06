#include "gurobi_c++.h"
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
#include <chrono>
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

class subtourelim: public GRBCallback {
    public :
        vector <GRBVar> y;
        vector <GRBVar> x;
        ulint n;
        ulint m;
        vector < pair < pair <ulint, ulint> , ulint> > E;
        map < pair <ulint, ulint>, ulint> mE;
        ulint root;

        subtourelim (vector <GRBVar> _y, vector <GRBVar> _x, ulint _n, ulint _m, vector < pair < pair <ulint, ulint>, ulint> > _E, map < pair <ulint, ulint>, ulint> _mE, ulint _root) {
            y = _y;
            x = _x;
            n = _n;
            m = _m;
            E = _E;
            mE = _mE;
            root = _root;
        }
    protected :
        // let C be the cycle that contains root, with n1 vertices
        // let S = n - n1
        // the sum of the edges that are not in C must be less than S-1
        void callback () {
            try {
                if (where == GRB_CB_MIPSOL) {
                    // find subtour containing vertex root, and add constraint to edges not in the tour
                    vector < list < pair <ulint, ulint> > > adj (n);

                    // build graph from solution vertices
                    for (ulint e = 0; e < m; e++) {
                        if (getSolution(x[e]) >= 0.5) {
                            ulint u = E[e].first.first;
                            ulint v = E[e].first.second;
                            ulint w = E[e].second;
                            adj[u].push_back(make_pair(v, w));
                            adj[v].push_back(make_pair(u, w));
                        }
                    }

                    vector <ulint> visitedVertices (n, 0);
                    vector <ulint> visitedEdges (m, 0);

                    queue <ulint> Q;

                    Q.push(root);

                    while (!Q.empty()) {
                        ulint u = Q.front();
                        Q.pop();

                        for (list < pair <ulint, ulint> > :: iterator it = adj[u].begin(); it != adj[u].end(); it++) {
                            ulint v = it->first;
                            if (visitedVertices[v] == 0 && getSolution(y[v]) >= 0.5) {
                                visitedVertices[v] = 1;
                                ulint e = mE[make_pair(u, v)];
                                visitedEdges[e] = 1;
                                Q.push(v);
                            }
                        }
                    }

                    ulint nUnvisitedVertices = 0;
                    for (ulint v = 0; v < n; v++) {
                        if (visitedVertices[v] == 0 && getSolution(y[v]) >= 0.5) {
                            nUnvisitedVertices++;
                        }
                    }

                    if (nUnvisitedVertices > 0) {
                        GRBLinExpr expr = 0;
                        for (ulint e = 0; e < m; e++) {
                            if (visitedEdges[e] == 0 && getSolution(x[e]) >= 0.5) {
                                if (visitedVertices[E[e].first.first] == 0 && getSolution(y[E[e].first.first]) >= 0.5
                                && visitedVertices[E[e].first.second] == 0 && getSolution(y[E[e].first.second]) >= 0.5)
                                expr += x[e];
                            }
                        }
                        addLazy(expr <= nUnvisitedVertices - 1);
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

int main (int argc, char * argv[]) {
    chrono :: steady_clock :: time_point tBegin = chrono :: steady_clock :: now();
    double timeLimit;

    if (argc == 2) {
        timeLimit = atof(argv[1]);
    } else {
        timeLimit = 10.0;
    }

    ulint n, mComplete, m, k, t, root;
    double d, p;

    cin >> n >> d >> k >> t >> p >> mComplete >> m >> root;

    vector <ulint> penalty (n); // vector with de penalties of each vectex
    matrix W (n, vector <lint> (n, -1)); // adjacency matrix for the complete graph
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
        env.set(GRB_StringParam_LogFile, "./output/N" + N + "D" + D + "K" + K + "T" + T + "P" + P + "/log.txt");
        env.set(GRB_DoubleParam_TimeLimit, timeLimit);

        GRBModel model = GRBModel(env);

        model.getEnv().set(GRB_IntParam_LazyConstraints, 1);
        model.getEnv().set(GRB_IntParam_LogToConsole, 0);
        model.getEnv().set(GRB_StringParam_LogFile, "./output/N" + N + "D" + D + "K" + K + "T" + T + "P" + P + "/log.txt");
        model.getEnv().set(GRB_DoubleParam_TimeLimit, timeLimit);

        vector <GRBVar> y (n);

        // ∀ v ∈ V
        for (ulint v = 0; v < n; v++) {
            // y_v ∈ {0.0, 1.0}
            y[v] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "y_" + itos(v));
        }

        vector <GRBVar> x (m);

        // ∀ e ∈ E
        for (ulint e = 0; e < m; e++) {
            ulint u, v;
            u = E[e].first.first;
            v = E[e].first.second;
            // x_e ∈ {0.0, 1.0}
            x[e] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "x_" + itos(u) + "_" + itos(v));
        }

        model.update();

        GRBLinExpr obj = 0.0;

        // obj = ∑ ce * xe
        for (ulint e = 0; e < m; e++) {
            ulint w;
            w = E[e].second;
            obj += w * x[e];
        }

        // obj += ∑ πv * (1 - yv)
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
            for (list < pair <ulint, ulint> > :: iterator it = adj[v].begin(); it != adj[v].end(); it++) {
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
            ulint solutionCost = 0;
            set <ulint> solutionVectices;
            set <ulint> solutionEdges;
            solutionCost = round(model.get(GRB_DoubleAttr_ObjVal));
            for (ulint v = 0; v < n; v++) {
                if (y[v].get(GRB_DoubleAttr_X) >= 0.5) {
                    solutionVectices.insert(v);
                }
            }
            for (ulint e = 0; e < m; e++) {
                if (x[e].get(GRB_DoubleAttr_X) >= 0.5) {
                    solutionEdges.insert(e);
                }
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
        } else {
            cout << "0 0 0" << endl;
        }

        // exporting model
        model.write("./output/N" + N + "D" + D + "K" + K + "T" + T + "P" + P + "/model.lp");

        chrono :: steady_clock :: time_point tEnd = chrono :: steady_clock :: now();
        chrono :: nanoseconds elapsedTime = chrono :: duration_cast <chrono :: nanoseconds> (tEnd - tBegin);
        ofstream elapsedTimeFile ("./output/N" + N + "D" + D + "K" + K + "T" + T + "P" + P + "/elapsedTime.txt", ofstream :: out);
        elapsedTimeFile << elapsedTime.count();
        elapsedTimeFile.close();

        ofstream gapFile ("./output/N" + N + "D" + D + "K" + K + "T" + T + "P" + P + "/gap.txt", ofstream :: out);
        gapFile << model.get(GRB_DoubleAttr_MIPGap);
        gapFile.close();

        ofstream objValFile ("./output/N" + N + "D" + D + "K" + K + "T" + T + "P" + P + "/objVal.txt", ofstream :: out);
        objVal << model.get(GRB_DoubleAttr_ObjVal);
        objValFile.close();
    } catch (GRBException e) {
        cout << "Error code = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    } catch (...) {
        cout << "Exception during opstimisation" << endl;
    }
    return 0;
}
