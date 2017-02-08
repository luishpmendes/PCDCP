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

using namespace std;

typedef long int lint;
typedef unsigned long int ulint;
typedef vector < vector <lint> > matrix;
typedef pair < vector <ulint>, lint > tSolution;

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

tSolution greedyRandomizedConstruction (matrix W, matrix Dist, matrix PI, vector <ulint> penalty, ulint root, vector < set <ulint> > Ns, double alpha, ulint seed) {
    tSolution result;
    list <ulint> lSolution;
    result.second = 0;
    for (ulint u = 0; u < W.size(); u++) {
        result.second += penalty[u];
    }

    vector <ulint> isInSolution(W.size(), 0);
    vector <ulint> isDominated(W.size(), 0);

    lSolution.insert(lSolution.begin(), root);
    result.second -= penalty[root];

    isInSolution[root] = 1;
    for (set <ulint> :: iterator it = Ns[root].begin(); it != Ns[root].end(); it++) {
        isDominated[*it] = 1;
    }

    // while Solution is not a complete solution do
    while ((ulint) accumulate(isDominated.begin(), isDominated.end(), 0) < isDominated.size()) {
        // populate canditate list with its respectives incremental costs and position
        vector < pair <list <ulint>, pair <lint, list <ulint> :: iterator> > > cantidateList;
        // try to populate the list
        for (ulint u = 0; u < W.size(); u++) {
            // check if the vertice u is not in the solution
            if (isInSolution[u] == 0) {
                // check if the vertice u dominates a not already dominated vertice
                int flag = 0;
                for (set <ulint> :: iterator it = Ns[u].begin(); it != Ns[u].end() && flag == 0; it++) {
                    if (isDominated[*it] == 0) {
                        flag = 1;
                    }
                }
                if (flag == 1) {
                    list <ulint> minU;
                    minU.push_back(u);
                    list <ulint> :: iterator minWIterator = lSolution.begin();
                    lint minCostU = -penalty[u];
                    ulint v = *(prev(lSolution.end()));
                    ulint w = *(minWIterator);
                    minCostU += Dist[v][u];
                    minCostU += Dist[u][w];
                    minCostU -= W[v][w];
                    set <ulint> verticesToBeAdded;
                    ulint x = PI[v][u]; // x = predecessor of u on the path from v to u
                    while (x >= 0 && x != v) {
                        verticesToBeAdded.insert(x);
                        minU.push_front(x);
                        x = PI[v][x];
                    }
                    x = PI[w][u]; // x = predecessor of u on the path from w to u
                    while (x >= 0 && x != w) {
                        verticesToBeAdded.insert(x);
                        minU.push_back(x);
                        x = PI[w][x];
                    }
                    for (set <ulint> :: iterator it = verticesToBeAdded.begin(); it != verticesToBeAdded.end(); it++) {
                        minCostU -= (1 - isInSolution[*it]) * penalty[*it];
                    }
                    list <ulint> :: iterator vIterator = lSolution.begin();
                    while (vIterator != lSolution.end()) {
                        list <ulint> :: iterator wIterator = next(vIterator);
                        if (wIterator != lSolution.end()) {
                            ulint v = *vIterator;
                            ulint w = *wIterator;
                            list <ulint> lU;
                            lU.push_back(u);
                            lint costU = -penalty[u];
                            costU += Dist[v][u];
                            costU += Dist[u][w];
                            costU -= W[v][w];
                            ulint x = PI[v][u];
                            set <ulint> verticesToBeAdded;
                            while (x >= 0 && x != v) {
                                verticesToBeAdded.insert(x);
                                lU.push_front(x);
                                x = PI[v][x];
                            }
                            x = PI[w][u];
                            while (x >= 0 && x != w) {
                                verticesToBeAdded.insert(x);
                                lU.push_back(x);
                                x = PI[w][x];
                            }
                            for (set <ulint> :: iterator it = verticesToBeAdded.begin(); it != verticesToBeAdded.end(); it++) {
                                costU -= (1 - isInSolution[*it]) * penalty[*it];
                            }
                            if (minCostU > costU) {
                                minCostU = costU;
                                minU = list <ulint> (lU.begin(), lU.end());
                                minWIterator = wIterator;
                            }
                        }
                        vIterator = next(vIterator);
                    }
                    cantidateList.push_back(make_pair(minU, make_pair(minCostU, minWIterator)));
                }
            }
        }
        if (cantidateList.size() > 0) {
            lint minCost, maxCost;
            minCost = maxCost = cantidateList[0].second.first;
            for (ulint i = 1; i < cantidateList.size(); i++) {
                if (minCost > cantidateList[i].second.first) {
                    minCost = cantidateList[i].second.first;
                }
                if (maxCost < cantidateList[i].second.first) {
                    maxCost = cantidateList[i].second.first;
                }
            }
            // build the restricted canditate list (RCL)
            double restriction = ((double) (maxCost - minCost));
            restriction *= alpha;
            restriction += ((double) minCost);
            vector < pair <list <ulint>, pair <lint, list <ulint> :: iterator> > > restrictedCantidateList;
            for (ulint i = 0; i < cantidateList.size(); i++) {
                if (cantidateList[i].second.first <= restriction) {
                    restrictedCantidateList.push_back(cantidateList[i]);
                }
            }

            // Select an element s from the RCL at random;
            default_random_engine generator (seed);
            uniform_int_distribution <ulint> distribution (0, restrictedCantidateList.size() - 1);
            ulint s = distribution(generator);
            lSolution.insert(restrictedCantidateList[s].second.second, restrictedCantidateList[s].first.begin(), restrictedCantidateList[s].first.end());
            result.second += restrictedCantidateList[s].second.first;
            for (list <ulint> :: iterator it = restrictedCantidateList[s].first.begin(); it != restrictedCantidateList[s].first.end(); it++) {
                isInSolution[*it] = 1;
                for (set <ulint> :: iterator it2 = Ns[*it].begin(); it2 != Ns[*it].end(); it2++) {
                    isDominated[*it2] = 1;
                }
            }
        }
    }
    result.first = vector <ulint> (lSolution.begin(), lSolution.end());
    return result;
}

// maneira mais generica: contador de quantos vértices da rota dominam cada vértice do grafo
// vejo os vertices v dominados por um vértice u da rota, se o contador[v] >= 2, posso eliminar u
// exclusao
bool mergeDominantVertices (matrix W, vector <ulint> penalty, ulint root, vector < set <ulint> > Ns, tSolution * solution) {
    bool result = false;
    vector <ulint> occurrencesCount (W.size(), 0);
    vector <ulint> dominatorsCount (W.size(), 0);
    for (ulint i = 0; i < (*solution).first.size(); i++) {
        ulint u = (*solution).first[i];
        occurrencesCount[u]++;
        for (set <ulint> :: iterator it = Ns[u].begin(); it != Ns[u].end(); it++) {
            ulint v = *it;
            dominatorsCount[v]++;
        }
    }
    ulint flag = 0;
    while (flag == 0) {
        flag = 1;
        for (ulint i = 0; i < (*solution).first.size(); i++) {
            ulint u = (*solution).first[i];
            // check to avoid removing the only occurrence of the root vertice 
            if (u != root || occurrencesCount[root] > 1) {
                ulint flag2 = 0;
                for (set <ulint> :: iterator it = Ns[u].begin(); it != Ns[u].end() && flag2 == 0; it++) {
                    ulint v = *it;
                    if (dominatorsCount[v] < 2) {
                        flag2 = 1;
                    }
                }
                // all vertices dominated by u are dominated by at least one other vertex in solution
                if (flag2 == 0) {
                    ulint prevU, nextU;
                    prevU = (*solution).first[(*solution).first.size() - 1];
                    if (i > 0) {
                        prevU = (*solution).first[i - 1];
                    }
                    nextU = (*solution).first[0];
                    if (i < (*solution).first.size() - 1) {
                        nextU = (*solution).first[i + 1];
                    }
                    // if there is an edge linking prevU with nextU
                    if (W[prevU][nextU] > 0) {
                        lint deltaCost = W[prevU][nextU] - W[prevU][u] - W[u][nextU];
                        if (occurrencesCount[u] == 1) {
                            deltaCost += penalty[u];
                        }
                        if (deltaCost < 0) {
                            result = true;
                            flag = 0;
                            (*solution).first.erase((*solution).first.begin() + i);
                            (*solution).second += deltaCost;
                            occurrencesCount[u]--;
                            for (set <ulint> :: iterator it = Ns[u].begin(); it != Ns[u].end(); it++) {
                                ulint v = *it;
                                dominatorsCount[v]--;
                            }
                        }
                    }
                }
            }
        }
    }
    return result;
}

// troca
bool swapDominantVertices (matrix W, vector <ulint> penalty, ulint root, vector < set <ulint> > Ns, tSolution * solution) {
    bool result = false;
    vector <ulint> occurrencesCount (W.size(), 0);
    for (ulint i = 0; i < (*solution).first.size(); i++) {
        ulint u = (*solution).first[i];
        occurrencesCount[u]++;
    }
    ulint flag = 0;
    while (flag == 0) {
        flag = 1;
        for (ulint i = 0; i < (*solution).first.size() && flag == 1; i++) {
            ulint u, prevU, nextU;
            u = (*solution).first[i];
            // check to avoid removing the only occurrence of the root vertice 
            if (u != root || occurrencesCount[root] > 1) {
                prevU = (*solution).first[(*solution).first.size() - 1];
                if (i > 0) {
                    prevU = (*solution).first[i - 1];
                }
                nextU = (*solution).first[0];
                if (i < (*solution).first.size() - 1) {
                    nextU = (*solution).first[i + 1];
                }
                for (set <ulint> :: iterator it = Ns[u].begin(); it != Ns[u].end(); it++) {
                    ulint v = *it;
                    vector <ulint> setDiff;
                    set_difference(Ns[u].begin(), Ns[u].end(), Ns[v].begin(), Ns[v].end(), inserter(setDiff, setDiff.begin()));
                    // if the neighborhood of v contains the neighborhood of u
                    if (setDiff.size() <= 0) {
                        // if there is edges linking v with the 'neighbors' of u
                        if (W[prevU][v] > 0 && W[v][nextU] > 0) {
                            lint deltaCost = W[prevU][v] + W[v][nextU] - W[prevU][u] - W[u][nextU];
                            if (occurrencesCount[u] == 1) {
                                deltaCost += penalty[u];
                            }
                            if (occurrencesCount[v] == 0) {
                                deltaCost -= penalty[v];
                            }
                            if (deltaCost < 0) {
                                result = true;
                                flag = 0;
                                (*solution).first[i] = v;
                                (*solution).second += deltaCost;
                                occurrencesCount[u]--;
                                occurrencesCount[v]++;
                            }
                        }
                    }
                }
            }
        }
    }
    return result;
}

bool twoOpt (matrix W, tSolution * solution) {
    bool result = false;
    ulint flag = 0;
    while (flag == 0) {
        flag = 1;
        for (ulint i = 0; i < (*solution).first.size() - 1 && flag == 1; i++) {
            ulint u, v;
            if (i > 0) {
                u = (*solution).first[i - 1];
            } else {
                u = (*solution).first[(*solution).first.size() - 1];
            }
            v = (*solution).first[i];
            for (ulint j = i + 1; j < (*solution).first.size() && flag == 1; j++) {
                ulint x, y;
                x = (*solution).first[j];
                if (j < (*solution).first.size() - 1) {
                    y = (*solution).first[j + 1];
                } else {
                    y = (*solution).first[0];
                }
                if (W[u][x] > 0 && W[v][y] > 0) {
                    if (W[u][v] + W[x][y] > W[u][x] + W[v][y]) {
                        result = true;
                        flag = 0;
                        reverse((*solution).first.begin() + i, (*solution).first.begin() + j + 1);
                        (*solution).second += W[u][x] + W[v][y] - W[u][v] - W[x][y];
                    }
                }
            }
        }
    }
    return result;
}

void localSearch (matrix W, vector <ulint> penalty, ulint root, vector < set <ulint> > Ns, tSolution * solution) {
    ulint flag = 0;
    while (flag == 0) {
        flag = 1;
        if (mergeDominantVertices (W, penalty, root, Ns, solution)) {
            flag = 0;
        }
        if (swapDominantVertices (W, penalty, root, Ns, solution)) {
            flag = 0;
        }
        if (twoOpt (W, solution)) {
            flag = 0;
        }
    }
}

tSolution grasp (matrix W, matrix Dist, matrix PI, vector <ulint> penalty, ulint root, vector < set <ulint> > Ns, ulint maxIterations, double alpha, ulint seed, vector <ulint> * solution, ulint * solutionCost) {
    tSolution result;
    for (ulint k = 0; k < maxIterations; k++) {
        tSolution solution = greedyRandomizedConstruction(W, Dist, PI, penalty, root, Ns, alpha, seed);
        localSearch(W, penalty, root, Ns, &solution);

        if (k == 0 || result.second > solution.second) {
            result.first = vector <ulint> (solution.first.begin(), solution.first.end());
            result.second = solution.second;
        }
    }
    return result;
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
    matrix Dist, PI; // matrices for the distance and precedence on the graph

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
    }

    vector < set <ulint> > Ns = neighbourhoods (WComplete, k);

    floydWarshall (W, &Dist, &PI);

    ulint seed = chrono::system_clock::now().time_since_epoch().count();

    vector <ulint> solution;
    ulint solutionCost = 0;

    grasp(W, Dist, PI, penalty, root, Ns, maxIterations, alpha, seed, &solution, &solutionCost);

    cout << solution.size() << ' ' << solution.size() << ' ' << solutionCost << endl;

    for (vector <ulint> :: iterator it = solution.begin(); it != solution.end(); it++) {
        ulint v = *it;
        cout << v << endl;
    }

    for (ulint i = 0; i < solution.size() - 1; i++) {
        pair <ulint, ulint> e = minmax(solution[i], solution[i + 1]);
        cout << e.first << ' ' << e.second << endl;
    }
    pair <ulint, ulint> e = minmax(solution[solution.size() - 1], solution[0]);
    cout << e.first << ' ' << e.second << endl;

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
    stringstream ssA;
    ssA << fixed << setprecision(1) << alpha;
    string A = ssA.str();
    A.erase(remove(A.begin(), A.end(), '.'), A.end());

    ofstream objValFile ("./output/N" + N + "D" + D + "K" + K + "T" + T + "P" + P + "I" + I + "A" + A + "/objVal.txt", ofstream :: out);
    objValFile << solutionCost;
    objValFile.close();

    chrono :: steady_clock :: time_point tEnd = chrono :: steady_clock :: now();
    chrono :: nanoseconds elapsedTime = chrono :: duration_cast <chrono :: nanoseconds> (tEnd - tBegin);
    ofstream elapsedTimeFile ("./output/N" + N + "D" + D + "K" + K + "T" + T + "P" + P + "I" + I + "A" + A + "/elapsedTime.txt", ofstream :: out);
    elapsedTimeFile << elapsedTime.count();
    elapsedTimeFile.close();

    return 0;
}
