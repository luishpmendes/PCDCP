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
#include <numeric>
#include <iterator>
#include <chrono>
#include <fstream>

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

void greedyRandomizedConstruction (matrix W, vector < list < pair <ulint, ulint> > > adj, vector <ulint> penalty, vector < set <ulint> > Ns, double alpha, ulint seed, vector <ulint> * solution, ulint * solutionCost) {
    cout << "entrou: greedyRandomizedConstruction" << endl;
    list <ulint> lSolution;
    ulint newSolutionCost = 0;
    for (ulint u = 0; u < W.size(); u++) {
        newSolutionCost += penalty[u];
    }
    vector <ulint> isInSolution(W.size(), 0);
    vector <ulint> isDominated(W.size(), 0);
    // while Solution is not a complete solution do
    while ((ulint) accumulate(isDominated.begin(), isDominated.end(), 0) < isDominated.size()) {
        // populate canditate list with its respectives incremental costs and position
        vector < pair <ulint, pair <lint, list <ulint> :: iterator> > > cantidateList;
        // try to populate with only not dominated vertices
        for (ulint u = 0; u < W.size(); u++) {
            if (isDominated[u] == 0) {
                lint minCostVW = -1;
                if (lSolution.size() <= 0) {
                    minCostVW = -penalty[u];
                } else if (W[u][*(lSolution.begin())] >= 0 && W[u][*(prev(lSolution.end()))] >= 0) {
                    minCostVW = -penalty[u];
                    minCostVW += W[u][*(lSolution.begin())];
                    minCostVW += W[u][*(prev(lSolution.end()))];
                    minCostVW -= W[*(lSolution.begin())][*(prev(lSolution.end()))];                    
                }
                list <ulint> :: iterator minWIterator = lSolution.begin();
                list <ulint> :: iterator vIterator = lSolution.begin();
                while (vIterator != lSolution.end()) {
                    list <ulint> :: iterator wIterator = next(vIterator);
                    ulint v = *vIterator;
                    ulint w = *wIterator;
                    lint costVW = INFINITE;
                    if (W[u][v] >= 0 && W[u][w] >= 0) {
                        costVW = -penalty[u];
                        costVW += W[u][v];
                        costVW += W[u][w];
                        costVW -= W[v][w];
                    }
                    if (minCostVW < 0 || minCostVW > costVW) {
                        minCostVW = costVW;
                        minWIterator = wIterator;
                    }
                    vIterator = next(vIterator);
                }
                if (minCostVW >= INFINITE) {
                    cantidateList.push_back(make_pair(u, make_pair(minCostVW, minWIterator)));
                }
            }
        }
        // if not possible, try to populate with vertices not in solution, but already dominated
        if (cantidateList.size() <= 0) {
            for (ulint u = 0; u < W.size(); u++) {
                if (isInSolution[u] == 0) {
                    lint minCostVW = -1;
                    if (lSolution.size() <= 0) {
                        minCostVW = -penalty[u];
                    } else if (W[u][*(lSolution.begin())] >= 0 && W[u][*(prev(lSolution.end()))] >= 0) {
                        minCostVW = -penalty[u];
                        minCostVW += W[u][*(lSolution.begin())];
                        minCostVW += W[u][*(prev(lSolution.end()))];
                        minCostVW -= W[*(lSolution.begin())][*(prev(lSolution.end()))];
                    }
                    list <ulint> :: iterator minWIterator = lSolution.begin();
                    list <ulint> :: iterator vIterator = lSolution.begin();
                    while (vIterator != lSolution.end()) {
                        list <ulint> :: iterator wIterator = next(vIterator);
                        ulint v = *vIterator;
                        ulint w = *wIterator;
                        lint costVW = INFINITE;
                        if (W[u][v] >= 0 && W[u][w] >= 0) {
                            costVW = -penalty[u];
                            costVW += W[u][v];
                            costVW += W[u][w];
                            costVW -= W[v][w];
                        }
                        if (minCostVW > costVW) {
                            minCostVW = costVW;
                            minWIterator = wIterator;
                        }
                        vIterator = next(vIterator);
                    }
                    if (minCostVW < INFINITE) {
                        cantidateList.push_back(make_pair(u, make_pair(minCostVW, minWIterator)));
                    }
                }
            }
        }
        lint minCost, maxCost;
        minCost = maxCost = 0;
        if (cantidateList.size() > 0) {
            minCost = maxCost = cantidateList[0].second.first;
        }
        for (ulint i = 1; i < cantidateList.size(); i++) {
            if (minCost > cantidateList[i].second.first) {
                minCost = cantidateList[i].second.first;
            }
            if (maxCost < cantidateList[i].second.first) {
                maxCost = cantidateList[i].second.first;
            }
        }
        // build the restricted canditate list (RCL)
        double aux = ((double) (maxCost - minCost));
        aux *= alpha;
        aux += ((double) minCost);
        vector < pair <ulint, pair <ulint, list <ulint> :: iterator> > > restrictedCantidateList;
        for (ulint i = 0; i < cantidateList.size(); i++) {
            if (cantidateList[i].second.first <= aux) {
                restrictedCantidateList.push_back(cantidateList[i]);
            }
        }
        // Select an element s from the RCL at random;
        default_random_engine generator (seed);
        uniform_int_distribution <ulint> distribution (0, restrictedCantidateList.size() - 1);
        ulint s = distribution(generator);
        lSolution.insert(restrictedCantidateList[s].second.second, restrictedCantidateList[s].first);
        newSolutionCost += restrictedCantidateList[s].second.first;
        isInSolution[restrictedCantidateList[s].first] = 1;
        isDominated[restrictedCantidateList[s].first] = 1;
        for (set <ulint> :: iterator it = Ns[restrictedCantidateList[s].first].begin(); it != Ns[restrictedCantidateList[s].first].end(); it++) {
            isDominated[*it] = 1;
        }
    }
    (*solution) = vector <ulint> (lSolution.begin(), lSolution.end());
    (*solutionCost) = newSolutionCost;
}

// maneira mais generica: contador de quantos vértices da rota dominam cada vértice do grafo
// vejo os vertices v dominados por um vértice u da rota, se o contador[v] >= 2, posso eliminar u
// exclusao
bool mergeDominantVertices (matrix W, vector < list < pair <ulint, ulint> > > adj, vector <ulint> penalty, vector < set <ulint> > Ns, vector <ulint> * solution, ulint * solutionCost) {
    cout << "entrou: mergeDominantVertices" << endl;
    bool result = false;
    vector <ulint> dominatorsCount (W.size(), 0);
    for (ulint i = 0; i < (*solution).size(); i++) {
        ulint u = (*solution)[i];
        for (set <ulint> :: iterator it = Ns[u].begin(); it != Ns[u].end(); it++) {
            ulint v = *it;
            dominatorsCount[v]++;
        }
    }
    ulint flag = 0;
    while (flag == 0) {
        flag = 1;
        for (ulint i = 0; i < (*solution).size(); i++) {
            ulint u = (*solution)[i];
            ulint flag2 = 0;
            for (set <ulint> :: iterator it = Ns[u].begin(); it != Ns[u].end() && flag == 0; it++) {
                ulint v = *it;
                if (dominatorsCount[v] < 2) {
                    flag2 = 1;
                }
            }
            // all vertices dominated by u are dominated by at least one other vertex in solution
            if (flag2 == 0) {
                ulint prevU, nextU;
                prevU = (*solution)[(*solution).size() - 1];
                if (i > 0) {
                    prevU = (*solution)[i - 1];
                }
                nextU = (*solution)[0];
                if (i < (*solution).size() - 1) {
                    nextU = (*solution)[i + 1];
                }
                // if there is an edge linking prevU with nextU
                if (W[prevU][nextU] > 0) {
                    lint deltaCost = W[prevU][nextU] - W[prevU][u] - W[u][nextU] + penalty[u];
                    if (deltaCost < 0) {
                        result = true;
                        flag = 0;
                        (*solution).erase((*solution).begin() + i);
                        (*solutionCost) += deltaCost;
                        for (set <ulint> :: iterator it = Ns[u].begin(); it != Ns[u].end(); it++) {
                            ulint v = *it;
                            dominatorsCount[v]--;
                        }
                    }
                }
            }
        }
    }
    return result;
}

// troca
bool swapDominantVertices (matrix W, vector < list < pair <ulint, ulint> > > adj, vector <ulint> penalty, vector < set <ulint> > Ns, vector <ulint> * solution, ulint * solutionCost) {
    cout << "entrou: swapDominantVertices" << endl;
    bool result = false;
    ulint flag = 0;
    while (flag == 0) {
        flag = 1;
        for (ulint i = 0; i < (*solution).size() && flag == 1; i++) {
            ulint u, prevU, nextU;
            u = (*solution)[i];
            prevU = (*solution)[(*solution).size() - 1];
            if (i > 0) {
                prevU = (*solution)[i - 1];
            }
            nextU = (*solution)[0];
            if (i < (*solution).size() - 1) {
                nextU = (*solution)[i + 1];
            }
            for (set <ulint> :: iterator it = Ns[u].begin(); it != Ns[u].end(); it++) {
                ulint v = *it;
                vector <ulint> setDiff;
                set_difference(Ns[u].begin(), Ns[u].end(), Ns[v].begin(), Ns[v].end(), setDiff.begin());
                // if the neighborhood of v contains the neighborhood of u
                if (setDiff.size() <= 0) {
                    // if there is edges linking v with the 'neighbors' of u
                    if (W[prevU][v] > 0 && W[v][nextU] > 0) {
                        lint deltaCost = -penalty[v] + W[prevU][v] + W[v][nextU] + penalty[u] - W[prevU][u] - W[u][nextU];
                        if (deltaCost < 0) {
                            result = true;
                            flag = 0;
                            (*solution)[i] = v;
                            (*solutionCost) += deltaCost;
                        }
                    }
                }
            }
        }
    }
    return result;
}

bool twoOpt (matrix W, vector < list < pair <ulint, ulint> > > adj, vector <ulint> * solution, ulint * solutionCost) {
    cout << "entrou: twoOpt" << endl;
    bool result = false;
    ulint flag = 0;
    while (flag == 0) {
        flag = 1;
        for (ulint i = 0; i < (*solution).size() - 1 && flag == 1; i++) {
            ulint u, v;
            if (i > 0) {
                u = (*solution)[i - 1];
            } else {
                u = (*solution)[(*solution).size() - 1];
            }
            v = (*solution)[i];
            for (ulint j = i + 1; j < (*solution).size() && flag == 1; j++) {
                ulint x, y;
                x = (*solution)[j];
                if (j < (*solution).size()) {
                    y = (*solution)[j + 1];
                } else {
                    y = (*solution)[j];
                }
                if (W[u][x] > 0 && W[v][y] > 0) {
                    if (W[u][v] + W[x][y] > W[u][x] + W[v][y]) {
                        result = true;
                        flag = 0;
                        for (ulint k = 0; k <= j - i; k++) {
                            (*solution)[i + k] = (*solution)[j - k];
                        }
                        (*solutionCost) += W[u][x] + W[v][y] - W[u][v] - W[x][y];
                    }
                }
            }
        }
    }
    return result;
}

void localSearch (matrix W, vector < list < pair <ulint, ulint> > > adj, vector <ulint> penalty, vector < set <ulint> > Ns, vector <ulint> * solution, ulint * solutionCost) {
    ulint flag = 0;
    while (flag == 0) {
        flag = 1;
        if (mergeDominantVertices (W, adj, penalty, Ns, solution, solutionCost)) {
            flag = 0;
        }
        if (swapDominantVertices (W, adj, penalty, Ns, solution, solutionCost)) {
            flag = 0;
        }
        if (twoOpt (W, adj, solution, solutionCost)) {
            flag = 0;
        }
    }
}

vector <ulint> grasp (matrix W, vector < list < pair <ulint, ulint> > > adj, vector <ulint> penalty, vector < set <ulint> > Ns, ulint maxIterations, double alpha, ulint seed) {

    vector <ulint> bestSolution;
    ulint bestSolutionCost = 0;

    for (ulint k = 0; k < maxIterations; k++) {
        vector <ulint> solution;
        ulint solutionCost = 0;

        greedyRandomizedConstruction(W, adj, penalty, Ns, alpha, seed, &solution, &solutionCost);
        localSearch(W, adj, penalty, Ns, &solution, &solutionCost);

        if (k == 0 || bestSolutionCost > solutionCost) {
            bestSolution = solution;
            bestSolutionCost = solutionCost;
        }
    }

    return bestSolution;
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

    ulint seed = chrono::system_clock::now().time_since_epoch().count();

    vector <ulint> solution = grasp(W, adj, penalty, Ns, maxIterations, alpha, seed);

    ulint solutionCost = 0;

    for (vector <ulint> :: iterator it = solution.begin(); it != solution.end(); it++) {
        ulint v = *it;
        solutionCost += penalty[v];
    }

    for (ulint i = 0; i < solution.size() - 1; i++) {
        solutionCost += W[solution[i]][solution[i + 1]];
    }
    solutionCost += W[solution[solution.size() - 1]][solution[0]];

    cout << solution.size() << ' ' << solution.size() << ' ' << solutionCost << endl;

    for (vector <ulint> :: iterator it = solution.begin(); it != solution.end(); it++) {
        ulint v = *it;
        cout << v << endl;
    }

    for (ulint i = 0; i < solution.size() - 1; i++) {
        cout << solution[i] << ' ' << solution[i + 1] << endl;
    }
    cout << solution[solution.size() - 1] << ' ' << solution[0] << endl;

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

    ofstream objValFile ("./output/N" + N + "D" + D + "K" + K + "T" + T + "P" + P + "I" + I + "/objVal.txt", ofstream :: out);
    objValFile << solutionCost;
    objValFile.close();

    chrono :: steady_clock :: time_point tEnd = chrono :: steady_clock :: now();
    chrono :: nanoseconds elapsedTime = chrono :: duration_cast <chrono :: nanoseconds> (tEnd - tBegin);
    ofstream elapsedTimeFile ("./output/N" + N + "D" + D + "K" + K + "T" + T + "P" + P + "I" + I + "/elapsedTime.txt", ofstream :: out);
    elapsedTimeFile << elapsedTime.count();
    elapsedTimeFile.close();

    return 0;
}
