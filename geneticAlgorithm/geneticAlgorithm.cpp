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
typedef pair < vector <ulint>, lint > tGenotype; // permutation
typedef pair < vector < pair <ulint, bool> >, lint > tPhenotype; // circular list

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

// from circular list
lint fitnessFunction (matrix W, vector <ulint> penalty, vector < pair <ulint, bool> > solution) {
    lint result = 0;
    result = 0;
    vector <int> isInSolution (W.size(), 0);
    for (ulint i = 0; i < solution.size() - 1; i++) {
        isInSolution[solution[i].first] = 1;
        result += W[solution[i].first][solution[i + 1].first];
    }
    isInSolution[solution[solution.size() - 1].first] = 1;
    result += W[solution[solution.size() - 1].first][solution[0].first];
    for (ulint u = 0; u < isInSolution.size() && u < penalty.size(); u++) {
        result += (1 - isInSolution[u]) * penalty[u];
    }
    return result;
}
/*
vector <ulint> circularList2Permutation (ulint n, vector < set <ulint> > Ns, vector <ulint> solution) {
    vector <ulint> result (n, 0);
    vector <ulint> isDominated(n, 0);
    vector <int> isInResult (n, 0);
    ulint i = 0;
    for (ulint j = 0; j < solution.size() && i < n && (ulint) accumulate(isDominated.begin(), isDominated.end(), 0) < isDominated.size(); j++) {
        ulint u = solution[j];
        if (isInResult[u] == 0) {
            result[i++] = u;
            isInResult[u] = 1;
            for (set <ulint> :: iterator it = Ns[u].begin(); it != Ns[u].end(); it++) {
                isDominated[*it] = 1;
            }
        }
    }
    for (ulint u = 0; u < n && i < n; u++) {
        if (isInResult[u] == 0) {
            result[i++] = u;
        }
    }
    return result;
}
*/

/*
vector <ulint> circularList2Permutation (ulint n, matrix PI, vector < set <ulint> > Ns, vector <ulint> solution) {
    vector <ulint> result;

    vector <ulint> aux;
    
    ulint i = 0;

    while (i < solution.size()) {
        ulint j = i;

        while (j < solution.size() - 1 && PI[solution[i]][solution[j + 1]] == solution[j]) {
            j++;
        }

        if (j < solution.size() - 1) {
            aux.push_back(i);
            i = j;
        } else {
            if (PI[solution[i]][solution[0]] != solution[j]) {
                aux.push_back(j);
            }
            break;
        }
    }

    return result;
}
*/

vector <ulint> circularList2Permutation (ulint n, vector < pair <ulint, bool> > solution) {
    vector <ulint> result;
    vector <bool> isInResult (n, false);
    for (ulint i = 0; i < solution.size(); i++) {
        if (solution[i].second) {
            if (!isInResult[solution[i].first]) { // if necessário??
                result.push_back(solution[i].first);
                isInResult[solution[i].first] = true;
            }
        }
    }
    for (ulint u = 0; u < n; u++) {
        if (!isInResult[u]) {
            result.push_back(u);
            isInResult[u] = true;
        }
    }
    return result;
}

// from circular list to permutation
tGenotype encode (matrix W, vector <ulint> penalty, tPhenotype solution) {
    tGenotype result;
    result.first = circularList2Permutation (W.size(), solution.first);
    result.second = fitnessFunction (W, penalty, solution.first);
    return result;
}

vector < pair <ulint, bool> > permutation2circularList (ulint n, matrix PI, vector < set <ulint> > Ns, vector <ulint> solution) {
    vector < pair <ulint, bool> > result;
    vector <ulint> aux;
    vector <ulint> isDominated(n, 0);
    for (ulint i = 0; i < solution.size() && (ulint) accumulate(isDominated.begin(), isDominated.end(), 0) < isDominated.size(); i++) {
        ulint u = solution[i];
        aux.push_back(u);
        for (set <ulint> :: iterator it = Ns[u].begin(); it != Ns[u].end(); it++) {
            isDominated[*it] = 1;
        }
    }
    for (ulint i = 0; i < aux.size() - 1; i++) {
        ulint u, v;
        u = aux[i];
        v = aux[i + 1];
        while (u != v) {
            if (u == aux[i]) {
                result.push_back(make_pair(u, true));
            } else {
                result.push_back(make_pair(u, false));
            }
            u = PI[v][u];
        }
    }
    ulint u, v;
    u = aux[aux.size() - 1];
    v = aux[0];
    while (u != v) {
        if (u == aux[aux.size() - 1]) {
            result.push_back(make_pair(u, true));
        } else {
            result.push_back(make_pair(u, false));
        }
        u = PI[v][u];
    }
    return result;
}

// from permutation to circular list
tPhenotype decode (matrix W, matrix PI, vector <ulint> penalty, vector < set <ulint> > Ns, tGenotype solution) {
    vector < pair <ulint, bool> > circularList = permutation2circularList (W.size(), PI, Ns, solution.first);
    lint solutionCost = fitnessFunction (W, penalty, circularList);
    return make_pair(circularList, solutionCost);
}

void fixCircularList (ulint root, vector < pair <ulint, bool> > * circularList) {
    vector < pair <ulint, bool> > :: iterator it = (*circularList).begin();
    while (it != (*circularList).end()) {
        if ((*it).first == root && (*it).second) {
            break;
        }
        it++;
    }
    if (it != (*circularList).end()) {
        rotate((*circularList).begin(), it, (*circularList).end());
    }
}

// candidata a sair
bool mergeDominantVertices (matrix W, vector <ulint> penalty, ulint root, vector < set <ulint> > Ns, tPhenotype * solution) {
    bool result = false;
    vector <ulint> occurrencesCount (W.size(), 0);
    vector <ulint> dominatorsCount (W.size(), 0);
    for (ulint i = 0; i < (*solution).first.size(); i++) {
        ulint u = (*solution).first[i].first;
        occurrencesCount[u]++;
        if ((*solution).first[i].second) {
            // if u é representante
            for (set <ulint> :: iterator it = Ns[u].begin(); it != Ns[u].end(); it++) {
                ulint v = *it;
                dominatorsCount[v]++;
            }
        }
    }
    ulint flag = 0;
    while (flag == 0) {
        flag = 1;
        for (ulint i = 0; i < (*solution).first.size(); i++) {
            ulint u = (*solution).first[i].first;
            // check to avoid removing the representant occurrence of the root vertice 
            if (u != root || (*solution).first[i].second == false) {
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
                    prevU = (*solution).first[(*solution).first.size() - 1].first;
                    if (i > 0) {
                        prevU = (*solution).first[i - 1].first;
                    }
                    nextU = (*solution).first[0].first;
                    if (i < (*solution).first.size() - 1) {
                        nextU = (*solution).first[i + 1].first;
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
                            fixCircularList (root, &((*solution).first));
                        }
                    }
                }
            }
        }
    }
    return result;
}

// se u for representante, v (q vai substitui-lo) será tambem
// eliminar esta busca
bool swapDominantVertices (matrix W, vector <ulint> penalty, ulint root, vector < set <ulint> > Ns, tPhenotype * solution) {
    bool result = false;
    vector <ulint> occurrencesCount (W.size(), 0);
    for (ulint i = 0; i < (*solution).first.size(); i++) {
        ulint u = (*solution).first[i].first;
        occurrencesCount[u]++;
    }
    ulint flag = 0;
    while (flag == 0) {
        flag = 1;
        for (ulint i = 0; i < (*solution).first.size() && flag == 1; i++) {
            ulint u, prevU, nextU;
            u = (*solution).first[i].first;
            // check to avoid removing the representant occurrence of the root vertice 
            if (u != root || (*solution).first[i].second == false) {
                prevU = (*solution).first[(*solution).first.size() - 1].first;
                if (i > 0) {
                    prevU = (*solution).first[i - 1].first;
                }
                nextU = (*solution).first[0].first;
                if (i < (*solution).first.size() - 1) {
                    nextU = (*solution).first[i + 1].first;
                }
                for (set <ulint> :: iterator it = Ns[u].begin(); it != Ns[u].end() && flag == 1; it++) {
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
                                (*solution).first[i].first = v;
                                (*solution).second += deltaCost;
                                occurrencesCount[u]--;
                                occurrencesCount[v]++;
                                fixCircularList (root, &((*solution).first));
                            }
                        }
                    }
                }
            }
        }
    }
    return result;
}

// cada nó da lista circular tem uma flag dizendo se ele é representante
// 2opt n muda os representantes
bool twoOpt (matrix W, ulint root, tPhenotype * solution) {
    bool result = false;
    ulint flag = 0;
    while (flag == 0) {
        flag = 1;
        for (ulint i = 0; i < (*solution).first.size() - 1 && flag == 1; i++) {
            ulint u, v;
            if (i > 0) {
                u = (*solution).first[i - 1].first;
            } else {
                u = (*solution).first[(*solution).first.size() - 1].first;
            }
            v = (*solution).first[i].first;
            for (ulint j = i + 1; j < (*solution).first.size() && flag == 1; j++) {
                ulint x, y;
                x = (*solution).first[j].first;
                if (j < (*solution).first.size() - 1) {
                    y = (*solution).first[j + 1].first;
                } else {
                    y = (*solution).first[0].first;
                }
                if (W[u][x] > 0 && W[v][y] > 0) {
                    if (W[u][v] + W[x][y] > W[u][x] + W[v][y]) {
                        result = true;
                        flag = 0;
                        reverse((*solution).first.begin() + i, (*solution).first.begin() + j + 1);
                        (*solution).second += W[u][x] + W[v][y] - W[u][v] - W[x][y];
                        fixCircularList (root, &((*solution).first));
                    }
                }
            }
        }
    }
    return result;
}

void localSearch (matrix W, vector <ulint> penalty, ulint root, vector < set <ulint> > Ns, tPhenotype * solution) {
    ulint flag = 0;
    while (flag == 0) {
        flag = 1;
        if (mergeDominantVertices (W, penalty, root, Ns, solution)) {
            flag = 0;
        }
        if (swapDominantVertices (W, penalty, root, Ns, solution)) {
            flag = 0;
        }
        if (twoOpt (W, root, solution)) {
            flag = 0;
        }
    }
}

void fixPermutation (ulint n, ulint root, vector <ulint> * permutation) {
    vector <ulint> geneCounter (n, 0);
    for (vector <ulint> :: iterator it = (*permutation).begin(); it != (*permutation).end(); it++) {
        geneCounter[*it]++;
    }
    while ((*permutation).size() > n) {
        (*permutation).pop_back();
    }
    while ((*permutation).size() < n) {
        bool flag = false;
        for (ulint i = 0; i < n; i++) {
            if (geneCounter[i] < 1) {
                (*permutation).push_back(i);
                geneCounter[i]++;
                flag = true;
            }
        }
        if (!flag) {
            break;
        }
    }
    for (ulint i = 0; i < (*permutation).size(); i++) {
        ulint u = (*permutation)[i];
        for (ulint v = 0; v < n && geneCounter[u] > 1; v++) {
            if (geneCounter[v] < 1) {
                (*permutation)[i] = v;
                geneCounter[u]--;
                geneCounter[v]++;
            }
        }
    }
    if ((*permutation)[0] != root) {
        ulint i = 0;
        while (i < (*permutation).size() && (*permutation)[i] != root) {
            i++;
        }
        if (i < (*permutation).size()) {
            (*permutation)[i] = (*permutation)[0];
            (*permutation)[0] = root;
        }
    }
}

set < tGenotype > initialPopulation (matrix W, matrix PI, vector <ulint> penalty, ulint root, vector < set <ulint> > Ns, ulint populationSize) {
    set < tGenotype > result;
    vector <ulint> sequence (W.size(), 0);
    for (ulint i = 0; i < sequence.size(); i++) {
        sequence[i] = i;
    }
    ulint seed = chrono :: system_clock :: now().time_since_epoch().count();
    auto engine = default_random_engine(seed);
    while (result.size() < populationSize) {
        vector <ulint> permutation (sequence.begin(), sequence.end());
        shuffle(permutation.begin(), permutation.end(), engine);
        if (permutation[0] != root) {
            ulint i = 0;
            while (i < permutation.size() && permutation[i] != root) {
                i++;
            }
            if (i < permutation.size()) {
                permutation[i] = permutation[0];
                permutation[0] = root;
            }
        }
        vector < pair <ulint, bool> > circularList = permutation2circularList (W.size(), PI, Ns, permutation);
        lint solutionCost = fitnessFunction (W, penalty, circularList);
        tGenotype chromossome = make_pair(permutation, solutionCost);
        tPhenotype individual = decode (W, PI, penalty, Ns, chromossome);
        localSearch (W, penalty, root, Ns, &individual);
        chromossome = encode (W, penalty, individual);
        individual = decode (W, PI, penalty, Ns, chromossome);
        chromossome = encode (W, penalty, individual);
        result.insert(chromossome);
    }
    return result;
}

void crossOver (matrix W, matrix PI, vector <ulint> penalty, ulint root, vector < set <ulint> > Ns, tGenotype parent1, tGenotype parent2, tGenotype * offspring1, tGenotype * offspring2) {
    // chooses crossover point at random
    ulint mutationSeed = chrono :: system_clock :: now().time_since_epoch().count();
    default_random_engine mutationGenerator (mutationSeed);
    uniform_int_distribution <ulint> mutationDistribution (0, W.size() - 1);
    ulint x = mutationDistribution(mutationGenerator);

    vector <int> isInOffspring1 (W.size(), 0);
    vector <int> isInOffspring2 (W.size(), 0);
    // copy genes before crossover point from parent to offspring
    for (ulint i = 0; i <= x; i++) {
        (*offspring1).first.push_back(parent1.first[i]);
        isInOffspring1[parent1.first[i]] = 1;
        (*offspring2).first.push_back(parent2.first[i]);
        isInOffspring2[parent2.first[i]] = 1;
    }
    // fill chromosome with genes from the other parent
    for (ulint i = 0; i <= parent1.first.size() && i < parent2.first.size(); i++) {
        if (isInOffspring1[parent2.first[i]] == 0) {
            (*offspring1).first.push_back(parent2.first[i]);
            isInOffspring1[parent2.first[i]] = 1;
        }
        if (isInOffspring2[parent1.first[i]] == 0) {
            (*offspring2).first.push_back(parent1.first[i]);
            isInOffspring2[parent1.first[i]] = 1;
        }
    }
    // fixing offsprings
    fixPermutation(W.size(), root, &((*offspring1).first));
    fixPermutation(W.size(), root, &((*offspring2).first));
    // computing its costs
    (*offspring1).second = fitnessFunction (W, penalty, permutation2circularList (W.size(), PI, Ns, (*offspring1).first));
    (*offspring2).second = fitnessFunction (W, penalty, permutation2circularList (W.size(), PI, Ns, (*offspring2).first));
}

void mutation (matrix W, matrix PI, vector <ulint> penalty, ulint root, vector < set <ulint> > Ns, double mutationRate, tGenotype * individual) {
    ulint mutationSeed1 = chrono :: system_clock :: now().time_since_epoch().count() + 1;
    default_random_engine mutationGenerator1 (mutationSeed1);
    uniform_real_distribution <double> mutationDistribution1 (0.0, 1.0);
    if (mutationDistribution1(mutationGenerator1) <= mutationRate) {
        ulint mutationSeed2 = chrono :: system_clock :: now().time_since_epoch().count() + 2;
        default_random_engine mutationGenerator2 (mutationSeed2);
        uniform_int_distribution <ulint> mutationDistribution2 (1, W.size() - 1);        
        ulint i = mutationDistribution2(mutationGenerator2);
        ulint j = mutationDistribution2(mutationGenerator2);
        ulint aux = (*individual).first[i];
        (*individual).first[i] = (*individual).first[j];
        (*individual).first[j] = aux;
        fixPermutation(W.size(), root, &((*individual).first));
        (*individual).second = fitnessFunction (W, penalty, permutation2circularList (W.size(), PI, Ns, (*individual).first));
    }
}

tGenotype selection (set < tGenotype > population) {
    tGenotype result;
    ulint at = 0;
    for (set < tGenotype > :: iterator it = population.begin(); it != population.end(); it++) {
        tGenotype individual = *it;
        at += individual.second;
    }
    ulint rouletteWheelSelectionSeed = chrono :: system_clock :: now().time_since_epoch().count();
    default_random_engine rouletteWheelSelectionGenerator (rouletteWheelSelectionSeed);
    uniform_int_distribution <ulint> rouletteWheelSelectionDistribution (0, at);
    ulint r = rouletteWheelSelectionDistribution(rouletteWheelSelectionGenerator);
    ulint sum = 0;
    for (set < tGenotype > :: iterator it = population.begin(); it != population.end(); it++) {
        tGenotype individual = *it;
        sum += individual.second;
        if (sum >= r) {
            result = individual;
            break;
        }
    }
    return result;
}

bool indValueComp (tGenotype lhs, tGenotype rhs) {
    return lhs.second < rhs.second;
}

// elitism
set < tGenotype > populationSubstitution (matrix W, matrix PI, vector <ulint> penalty, ulint root, vector < set <ulint> > Ns, ulint populationSize, double mutationRate, set < tGenotype > population) {
    set < tGenotype > result;
    tGenotype bestChromossome;
    bool flag = true;
    for (set <tGenotype> :: iterator it = population.begin(); it != population.end(); it++) {
        if (flag || bestChromossome.second > (*it).second) {
            flag = false;
            bestChromossome = *it;
        }
    }
    result.insert(bestChromossome);
    while (result.size() < populationSize) {
        tGenotype parent1 = selection(population);
        tGenotype parent2 = selection(population);
        tGenotype offspring1;
        tGenotype offspring2;

        crossOver (W, PI, penalty, root, Ns, parent1, parent2, &offspring1, &offspring2);

        mutation(W, PI, penalty, root, Ns, mutationRate, &offspring1);
        mutation(W, PI, penalty, root, Ns, mutationRate, &offspring2);

        tPhenotype individual1 = decode (W, PI, penalty, Ns, offspring1);
        localSearch (W, penalty, root, Ns, &individual1);
        offspring1 = encode (W, penalty, individual1);
        individual1 = decode (W, PI, penalty, Ns, offspring1);
        offspring1 = encode (W, penalty, individual1);
        result.insert(offspring1);

        if (result.size() < populationSize) {
            tPhenotype individual2 = decode (W, PI, penalty, Ns, offspring2);
            localSearch (W, penalty, root, Ns, &individual2);
            offspring2 = encode (W, penalty, individual2);
            individual2 = decode (W, PI, penalty, Ns, offspring2);
            offspring2 = encode (W, penalty, individual2);
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

tGenotype geneticAlgorithm (matrix W, matrix PI, vector <ulint> penalty, ulint root, vector < set <ulint> > Ns, chrono :: steady_clock :: time_point tBegin, double timeLimit, ulint populationSize, double mutationRate) {
    tGenotype result;
    bool flag = true;
    set < tGenotype > oldPopulation;
    set < tGenotype > newPopulation = initialPopulation (W, PI, penalty, root, Ns, populationSize);
    if (termination (tBegin, timeLimit) != true) {
        oldPopulation = set < tGenotype > (newPopulation.begin(), newPopulation.end());
        newPopulation = populationSubstitution(W, PI, penalty, root, Ns, populationSize, mutationRate, oldPopulation);
    }
    for (set < tGenotype > :: iterator it = newPopulation.begin(); it != newPopulation.end(); it++) {
        tGenotype individual = *it;
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

    tGenotype chromossome = geneticAlgorithm(W, PI, penalty, root, Ns, tBegin, timeLimit, populationSize, mutationRate);

    tPhenotype individual = decode (W, PI, penalty, Ns, chromossome);

    cout << individual.first.size() << ' ' << individual.first.size() << ' ' << individual.second << endl;

    for (vector < pair <ulint, bool> > :: iterator it = individual.first.begin(); it != individual.first.end(); it++) {
        ulint v = (*it).first;
        cout << v << endl;
    }

    for (ulint i = 0; i < individual.first.size() - 1; i++) {
        pair <ulint, ulint> e = minmax(individual.first[i].first, individual.first[i + 1].first);
        cout << e.first << ' ' << e.second << endl;
    }
    pair <ulint, ulint> e = minmax(individual.first[individual.first.size() - 1].first, individual.first[0].first);
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
    string PS = itos(populationSize);
    stringstream ssMR;
    ssMR << fixed << setprecision(1) << d;
    string MR = ssMR.str();
    MR.erase(remove(MR.begin(), MR.end(), '.'), MR.end());

    ofstream objValFile ("./output/N" + N + "D" + D + "K" + K + "T" + T + "P" + P + "I" + I + "PS" + PS + "MR" + MR + "/objVal.txt", ofstream :: out);
    objValFile << individual.second;
    objValFile.close();

    chrono :: steady_clock :: time_point tEnd = chrono :: steady_clock :: now();
    chrono :: nanoseconds elapsedTime = chrono :: duration_cast <chrono :: nanoseconds> (tEnd - tBegin);
    ofstream elapsedTimeFile ("./output/N" + N + "D" + D + "K" + K + "T" + T + "P" + P + "I" + I + "PS" + PS + "MR" + MR + "/elapsedTime.txt", ofstream :: out);
    elapsedTimeFile << elapsedTime.count();
    elapsedTimeFile.close();

    return 0;
}
