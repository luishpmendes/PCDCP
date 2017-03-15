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

string itos (ulint i) {
    stringstream s;
    s << i;
    return s.str();
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
        // let C be the cycle that contains u, with n1 vertices
        // let S = n - n1
        // the sum of the edges that are not in C must be less than S-1
        void callback () {
            try {
                if (where == GRB_CB_MIPSOL) {
                    // find subtour containing vertex u, and add constraint to edges not in the tour
                    vector < list < pair <ulint, ulint> > > adj (n);

                    // build graph from solution vertices
                    for (ulint e = 0; e < m; e++) {
                        if (getSolution(x[e]) >= 0.5) {
                            ulint a = E[e].first.first;
                            ulint b = E[e].first.second;
                            ulint c = E[e].second;
                            adj[a].push_back(make_pair(b, c));
                            adj[b].push_back(make_pair(a, c));
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
    string I ("0");
    ulint timeLimit = 10;

    if (argc >= 2) {
        I = string (argv[1]);
    }

    if (argc >= 3) {
        timeLimit = atoi(argv[2]);
    }

    ulint nComplete, k, t, n, m, root;
    double d;

    cin >> nComplete >> d >> k >> t >> n >> m >> root;

    vector <ulint> penalty (nComplete); // vector with de penalties of each vectex
    vector < list < pair <ulint, ulint> > > adj (nComplete); // adjacency lists for the graph

    for (ulint v = 0; v < nComplete; v++) {
        cin >> penalty[v];
    }

    vector <ulint> solutionV (nComplete, 0);

    // reading solution vertices
    for (ulint i = 0; i < n; i++) {
        ulint v;
        cin >> v;
        solutionV[v] = 1;
    }

    vector < pair < pair <ulint, ulint> , ulint> > E (m); // vector of edges with the format ((u, v), w)
    map < pair <ulint, ulint>, ulint> mE; // map an edge to its ID
    vector < vector <ulint> > paths (m);

    // reading graph
    for (ulint e = 0; e < m; e++) {
        ulint u, v, w, pathSize;
        cin >> u >> v >> w >> pathSize;
        adj[u].push_back(make_pair(v, w));
        adj[v].push_back(make_pair(u, w));
        E[e] = make_pair(make_pair(u, v), w);
        mE[make_pair(u, v)] = e;
        mE[make_pair(v, u)] = e;
        paths[e] = vector <ulint> (pathSize);
        for (ulint i = 0; i < pathSize; i++) {
            cin >> paths[e][i];
        }
    }

    try {
        string N = itos(nComplete);
        stringstream ssD;
        ssD << fixed << setprecision(1) << d;
        string D = ssD.str();
        D.erase(remove(D.begin(), D.end(), '.'), D.end());
        string K = itos(k);
        string T = itos(t);

        ifstream remainingTimeFile ("./output/N" + N + "D" + D + "K" + K + "T" + T + "I" + I + "/remainingTime.txt");
        lint remainingTime = 0;
        if (remainingTimeFile.is_open()) {
            remainingTimeFile >> remainingTime;
        }
        if (remainingTime > 0) {
            timeLimit += remainingTime;
        }

        GRBEnv env = GRBEnv();

        env.set(GRB_IntParam_LazyConstraints, 1);
        env.set(GRB_IntParam_LogToConsole, 0);
        env.set(GRB_StringParam_LogFile, "./output/N" + N + "D" + D + "K" + K + "T" + T + "I" + I + "/log2.txt");
        env.set(GRB_DoubleParam_TimeLimit, ((double) timeLimit));

        GRBModel model = GRBModel(env);

        model.getEnv().set(GRB_IntParam_LazyConstraints, 1);
        model.getEnv().set(GRB_IntParam_LogToConsole, 0);
        model.getEnv().set(GRB_StringParam_LogFile, "./output/N" + N + "D" + D + "K" + K + "T" + T + "I" + I + "/log2.txt");
        model.getEnv().set(GRB_DoubleParam_TimeLimit, ((double) timeLimit));

        vector <GRBVar> y (nComplete);

        // ∀ v ∈ V
        for (ulint v = 0; v < nComplete; v++) {
            // y_v ∈ {0.0, 1.0}
            y[v] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "y_" + itos(v));
        }

        vector <GRBVar> x (m);

        // ∀ e ∈ E
        for (ulint e = 0; e < m; e++) {
            ulint u, v;
            u = E[e].first.first;
            v = E[e].first.second;
            // y_e ∈ {0.0, 1.0}
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
        for (ulint v = 0; v < nComplete; v++) {
            obj += penalty[v] * (1.0 - y[v]);
        }

        model.setObjective(obj, GRB_MINIMIZE);

        // yu == 1
        model.addConstr(y[root] == 1.0, "c_0");

        // dominance
        // ∀ v ∈ V
        for (ulint v = 0; v < nComplete; v++) {
            if (solutionV[v] == 1) {
                GRBLinExpr constr = 0.0;
                constr += y[v];
                model.addConstr(constr == 1, "c_1_" + itos(v));
            }
        }

        // each vertex must have exactly two edges adjacent to itself
        // ∀ v ∈ V
        for (ulint v = 0; v < nComplete; v++) {
            // ∑ xe == 2 * yv , e ∈ δ({v})
            GRBLinExpr constr = 0.0;
            for (list < pair <ulint, ulint> > :: iterator it = adj[v].begin(); it != adj[v].end(); it++) {
                ulint w = (*it).first; // destination
                ulint e = mE[make_pair(v, w)];
                constr += x[e];
            }
            model.addConstr(constr == 2.0 * y[v], "c_2_" + itos(v));
        }

        subtourelim cb = subtourelim(y, x, nComplete, m, E, mE, root);
        model.setCallback(&cb);

        model.optimize();

        if (model.get(GRB_IntAttr_SolCount) > 0) {
            ulint solutionCost = 0;
            set <ulint> solutionVectices;
            vector < pair <ulint, ulint> > solutionEdges;
            solutionCost = round(model.get(GRB_DoubleAttr_ObjVal));
            for (ulint v = 0; v < nComplete; v++) {
                if (y[v].get(GRB_DoubleAttr_X) >= 0.5) {
                    solutionVectices.insert(v);
                }
            }
            for (ulint e = 0; e < m; e++) {
                if (x[e].get(GRB_DoubleAttr_X) >= 0.5) {
                    for (ulint i = 0; i < paths[e].size() - 1; i++) {
                        pair <ulint, ulint> edge;
                        if (paths[e][i] < paths[e][i + 1]) {
                            edge.first = paths[e][i];
                            edge.second = paths[e][i + 1];
                        } else {
                            edge.first = paths[e][i + 1];
                            edge.second = paths[e][i];
                        }
                        solutionEdges.push_back(edge);
                    }
                }
            }
            cout << solutionVectices.size() << ' ' << solutionEdges.size() << ' ' << solutionCost << endl;
            for (set <ulint> :: iterator it = solutionVectices.begin(); it != solutionVectices.end(); it++) {
                ulint v = *it;
                cout << v << endl;
            }
            for (vector < pair <ulint, ulint> > :: iterator it = solutionEdges.begin(); it != solutionEdges.end(); it++) {
                pair <ulint, ulint> e = *it;
                cout << e.first << " " << e.second << endl;
            }
        } else {
            cout << "0 0 0" << endl;
        }

        // exporting model
        model.write("./output/N" + N + "D" + D + "K" + K + "T" + T + "I" + I + "/model2.lp");

        ofstream objValFile ("./output/N" + N + "D" + D + "K" + K + "T" + T + "I" + I + "/objVal2.txt", ofstream :: out);
        objValFile << model.get(GRB_DoubleAttr_ObjVal);
        objValFile.close();

        ofstream gapFile ("./output/N" + N + "D" + D + "K" + K + "T" + T + "I" + I + "/gap2.txt", ofstream :: out);
        gapFile << model.get(GRB_DoubleAttr_MIPGap);
        gapFile.close();

        chrono :: steady_clock :: time_point tEnd = chrono :: steady_clock :: now();
        chrono :: nanoseconds elapsedTime = chrono :: duration_cast <chrono :: nanoseconds> (tEnd - tBegin);
        ofstream elapsedTimeFile ("./output/N" + N + "D" + D + "K" + K + "T" + T + "I" + I + "/elapsedTime2.txt", ofstream :: out);
        elapsedTimeFile << elapsedTime.count();
        elapsedTimeFile.close();
    } catch (GRBException e) {
        cout << "Error code = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    } catch (...) {
        cout << "Exception during opstimisation" << endl;
    }
    return 0;
}
