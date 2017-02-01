#include <iostream>
#include <vector>
#include <sstream>
#include <set>
#include <list>
#include <chrono>
#include <algorithm>
#include <iomanip>
#include <fstream>

using namespace std;

typedef long int lint;
typedef unsigned long int ulint;
typedef vector < vector <lint> > matrix;
typedef pair < vector <ulint>, lint > tIndividual;

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

set < tIndividual > initialPopulation (ulint populationSize) {
    set < tIndividual > result;
    return result;
}

void crossOver (matrix W, vector < list < pair <ulint, ulint> > > adj, matrix Dist, matrix PI, vector <ulint> penalty, ulint root, vector < set <ulint> > Ns, tIndividual parent1, tIndividual parent2, tIndividual * offspring1, tIndividual * offspring2) {
}

void mutation (matrix W, vector < list < pair <ulint, ulint> > > adj, matrix Dist, matrix PI, vector <ulint> penalty, ulint root, vector < set <ulint> > Ns, double mutationRate, tIndividual * individual) {
}

tIndividual selection (set < tIndividual > population) {
    tIndividual result;
    ulint at = 0;
    for (set < tIndividual > :: iterator it = population.begin(); it != population.end(); it++) {
        tIndividual individual = *it;
        at += individual.second;
    }
    ulint rouletteWheelSelectionSeed = chrono :: system_clock :: now().time_since_epoch().count();
    default_random_engine rouletteWheelSelectionGenerator (rouletteWheelSelectionSeed);
    uniform_int_distribution <ulint> rouletteWheelSelectionDistribution (0, at);
    ulint r = rouletteWheelSelectionDistribution(rouletteWheelSelectionGenerator);
    ulint sum = 0;
    for (set < tIndividual > :: iterator it = population.begin(); it != population.end(); it++) {
        tIndividual individual = *it;
        sum += individual.second;
        if (sum >= r) {
            result = individual;
            break;
        }
    }
    return result;
}

bool indValueComp (tIndividual lhs, tIndividual rhs) {
    return lhs.second < rhs.second;
}

// elitism
set < tIndividual > populationSubstitution (matrix W, vector < list < pair <ulint, ulint> > > adj, matrix Dist, matrix PI, vector <ulint> penalty, ulint root, vector < set <ulint> > Ns, ulint populationSize, double mutationRate, set < tIndividual > population) {
    set < tIndividual > result;
    vector < tIndividual > sortedPopulation (population.begin(), population.end());
    sort(sortedPopulation.begin(), sortedPopulation.end(), indValueComp);
    result.insert(sortedPopulation[0]);
    while (result.size() < populationSize) {
        tIndividual parent1 = selection(population);
        tIndividual parent2 = selection(population);
        tIndividual offspring1;
        tIndividual offspring2;
        crossOver (W, adj, Dist, PI, penalty, root, Ns, parent1, parent2, &offspring1, &offspring2);
        mutation(W, adj, Dist, PI, penalty, root, Ns, mutationRate, &offspring1);
        mutation(W, adj, Dist, PI, penalty, root, Ns, mutationRate, &offspring2);
        result.insert(offspring1);
        if (result.size() < population.size()) {
            result.insert(offspring2);
        }
    }
    return result;
}

bool termination (chrono :: steady_clock :: time_point tBegin, double timeLimit) {
    chrono :: steady_clock :: time_point tCurrent = chrono :: steady_clock :: now();
    chrono :: seconds elapsedTime = chrono :: duration_cast <chrono :: seconds> (tCurrent - tBegin);
    if (elapsedTime.count() >= timeLimit) {
        return true;
    }
    return false;
}

tIndividual geneticAlgorithm (matrix W, vector < list < pair <ulint, ulint> > > adj, matrix Dist, matrix PI, vector <ulint> penalty, ulint root, vector < set <ulint> > Ns, chrono :: steady_clock :: time_point tBegin, double timeLimit, ulint populationSize, double mutationRate) {
    tIndividual result;
    bool flag = true;
    set < tIndividual > oldPopulation;
    set < tIndividual > newPopulation = initialPopulation (populationSize);
    while (termination (tBegin, timeLimit) != true) {
        oldPopulation = set < tIndividual > (newPopulation.begin(), newPopulation.end());
        newPopulation = populationSubstitution(W, adj, Dist, PI, penalty, root, Ns, populationSize, mutationRate, oldPopulation);
    }
    for (set < tIndividual > :: iterator it = newPopulation.begin(); it != newPopulation.end(); it++) {
        tIndividual individual = *it;
        if (flag || result.second > individual.second) {
            flag = false;
            result = individual;
        }
    }
    return result;
}

int main (int argc, char * argv[]) {
    chrono :: steady_clock :: time_point tBegin = chrono :: steady_clock :: now();
    string I ("0");
    double timeLimit = 100;
    ulint populationSize = 100;
    double mutationRate = 0.1;

    if (argc >= 2) {
        I = string (argv[1]);
    }

    if (argc >= 3) {
        timeLimit = atof(argv[2]);
    }

    if (argc >= 4) {
        populationSize = atoi(argv[3]);
    }

    if (argc >= 5) {
        mutationRate = atof(argv[4]);
    }

    ulint n, mComplete, m, k, t, root;
    double d, p;

    cin >> n >> d >> k >> t >> p >> mComplete >> m >> root;

    vector <ulint> penalty (n); // vector with de penalties of each vectex
    matrix WComplete (n, vector <lint> (n, -1)); // adjacency matrix for the complete graph
    matrix W (n, vector <lint> (n, -1)); // adjacency matrix for the graph
    vector < list < pair <ulint, ulint> > > adj (n); // adjacency lists for the graph
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
        adj[u].push_back(make_pair(v, w));
        adj[v].push_back(make_pair(u, w));
    }

    vector < set <ulint> > Ns = neighbourhoods (WComplete, k);

    floydWarshall (W, &Dist, &PI);

    tIndividual solution = geneticAlgorithm(W, adj, Dist, PI, penalty, root, Ns, tBegin, timeLimit, populationSize, mutationRate);

    cout << solution.first.size() << ' ' << solution.first.size() << ' ' << solution.second << endl;

    for (vector <ulint> :: iterator it = solution.first.begin(); it != solution.first.end(); it++) {
        ulint v = *it;
        cout << v << endl;
    }

    for (ulint i = 0; i < solution.first.size() - 1; i++) {
        pair <ulint, ulint> e = minmax(solution.first[i], solution.first[i + 1]);
        cout << e.first << ' ' << e.second << endl;
    }
    pair <ulint, ulint> e = minmax(solution.first[solution.first.size() - 1], solution.first[0]);
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

    ofstream objValFile ("./output/N" + N + "D" + D + "K" + K + "T" + T + "P" + P + "I" + I + "/objVal.txt", ofstream :: out);
    objValFile << solution.second;
    objValFile.close();

    chrono :: steady_clock :: time_point tEnd = chrono :: steady_clock :: now();
    chrono :: nanoseconds elapsedTime = chrono :: duration_cast <chrono :: nanoseconds> (tEnd - tBegin);
    ofstream elapsedTimeFile ("./output/N" + N + "D" + D + "K" + K + "T" + T + "P" + P + "I" + I + "/elapsedTime.txt", ofstream :: out);
    elapsedTimeFile << elapsedTime.count();
    elapsedTimeFile.close();

    return 0;
}