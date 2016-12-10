#include "gurobi_c++.h"
#include <iostream>
#include <vector>
#include <set>
#include <sstream>
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

int main (int argc, char * argv[]) {
    chrono :: steady_clock :: time_point tBegin = chrono :: steady_clock :: now();
    string I ("0");
    double timeLimit = 10.0;

    if (argc >= 2) {
        I = string (argv[1]);
    }

    if (argc >= 3) {
        timeLimit = atof(argv[2]);
    }

    ulint n, mComplete, m, k, t, root;
    double d, p;

    cin >> n >> d >> k >> t >> p >> mComplete >> m >> root;

    vector <ulint> penalty (n); // vector with the penalties of each vertex
    matrix Wcomplete (n, vector <lint> (n, -1)); // adjacency matrix for the complete graph
    matrix W (n, vector <lint> (n, -1)); // adjacency matrix for the graph
    vector < pair < pair <ulint, ulint> , ulint> > E (m); // vector of edges with the format ((u, v), w)
    matrix Dist, PI; // matrices for the distance and precedence on the graph

    for (ulint i = 0; i < n; i++) {
        Wcomplete[i][i] = 0;
        W[i][i] = 0;
    }

    for (ulint i = 0; i < n; i++) {
        double x, y;
        cin >> x >> y >> penalty[i];
    }

    // reading complete graph
    for (ulint e = 0; e < mComplete; e++) {
        ulint u, v, w;
        cin >> u >> v >> w;
        Wcomplete[u][v] = w;
        Wcomplete[v][u] = w;
    }

    // reading graph
    for (ulint e = 0; e < m; e++) {
        ulint u, v, w;
        cin >> u >> v >> w;
        E[e] = make_pair(make_pair(u, v), w);
        W[u][v] = w;
        W[v][u] = w;
    }

    floydWarshall (W, &Dist, &PI);

    vector < set <ulint> > Ns = neighbourhoods (Wcomplete, k);

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
        env.set(GRB_StringParam_LogFile, "./output/N" + N + "D" + D + "K" + K + "T" + T + "P" + P + "I" + I + "/log1.txt");
        env.set(GRB_DoubleParam_TimeLimit, timeLimit);

        GRBModel model = GRBModel(env);

        model.getEnv().set(GRB_IntParam_LazyConstraints, 1);
        model.getEnv().set(GRB_IntParam_LogToConsole, 0);
        model.getEnv().set(GRB_StringParam_LogFile, "./output/N" + N + "D" + D + "K" + K + "T" + T + "P" + P + "I" + I + "/log1.txt");
        model.getEnv().set(GRB_DoubleParam_TimeLimit, timeLimit);

        vector <GRBVar> y (n);

        // ∀ v ∈ V
        for (ulint v = 0; v < n; v++) {
            // y_v ∈ {0.0, 1.0}
            y[v] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "y_" + itos(v));
        }

        model.update();

        GRBLinExpr obj = 0.0;

        // obj = ∑ πv * (1 - yv)
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
                if (y[v].get(GRB_DoubleAttr_X) >= 0.5) {
                    solutionV.push_back(v);
                }
            }
            vector < pair < pair <ulint, ulint> , ulint> > solutionE;
            for (ulint u = 0; u < n; u++) {
                for (ulint v = u + 1; v < n; v++) {
                    if (y[u].get(GRB_DoubleAttr_X) >= 0.5 && y[v].get(GRB_DoubleAttr_X) >= 0.5) {
                        pair < pair <ulint, ulint> , ulint> edge;
                        edge.first.first = u;
                        edge.first.second = v;
                        edge.second = Dist[u][v];
                        solutionE.push_back(edge);
                    }
                }
            }

            cout << n << ' ' << d << ' ' << k << ' ' << t << ' ' << p << ' ' << solutionV.size() << ' ' << solutionE.size() << ' ' << root << endl;

            for (ulint i = 0; i < n; i++) {
                cout << penalty[i] << endl;
            }

            for (vector <ulint> :: iterator it = solutionV.begin(); it != solutionV.end(); it++) {
                ulint v = *it;
                cout << v << endl;
            }

            for (vector < pair < pair <ulint, ulint> , ulint> > :: iterator it = solutionE.begin(); it != solutionE.end(); it++) {
                pair < pair <ulint, ulint> , ulint> e = *it;
                ulint u, v, w;
                u = e.first.first;
                v = e.first.second;
                w = e.second;
                cout << u << ' ' << v << ' ' << w;
                vector <ulint> path;
                while (v != u) {
                    path.push_back(v);
                    v = PI[u][v];
                }
                path.push_back(v);
                cout << ' ' << path.size();
                for (vector <ulint> :: reverse_iterator it = path.rbegin(); it != path.rend(); it++) {
                    cout << ' ' << *it;
                }
                cout << endl;
            }
        }

        // exporting model
        model.write("./output/N" + N + "D" + D + "K" + K + "T" + T + "P" + P + "I" + I + "/model1.lp");

        ofstream objValFile ("./output/N" + N + "D" + D + "K" + K + "T" + T + "P" + P + "I" + I + "/objVal1.txt", ofstream :: out);
        objValFile << model.get(GRB_DoubleAttr_ObjVal);
        objValFile.close();

        ofstream gapFile ("./output/N" + N + "D" + D + "K" + K + "T" + T + "P" + P + "I" + I + "/gap1.txt", ofstream :: out);
        gapFile << model.get(GRB_DoubleAttr_MIPGap);
        gapFile.close();

        chrono :: steady_clock :: time_point tEnd = chrono :: steady_clock :: now();
        chrono :: nanoseconds elapsedTime = chrono :: duration_cast <chrono :: nanoseconds> (tEnd - tBegin);
        ofstream elapsedTimeFile ("./output/N" + N + "D" + D + "K" + K + "T" + T + "P" + P + "I" + I + "/elapsedTime1.txt", ofstream :: out);
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
