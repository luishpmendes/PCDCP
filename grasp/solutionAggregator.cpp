#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <map>
#include <set>

using namespace std;

typedef long int lint;
typedef unsigned long int ulint;
typedef vector < vector <lint> > matrix;

double average (vector <double> v) {
    if (v.size() <= 0) {
        return 1E100;
    }
    double result = ((double) accumulate(v.begin(), v.end(), 0.0));
    result /= ((double) v.size());
    return result;
}

double stdev (vector <double> v) {
    if (v.size() <= 0) {
        return 1E100;
    }
    vector <double> aux (v.size(), 0);
    double mean = average(v);
    for (ulint i = 0; i < v.size(); i++) {
        aux[i] = v[i];
        aux[i] -= mean;
        aux[i] *= aux[i];
    }
    double result = ((double) accumulate(aux.begin(), aux.end(), 0.0));
    result /= ((double) v.size());
    result = sqrt(result);
    return result;
}

struct keyComp {
    bool operator() (const pair < double, ulint > & lhs, const pair < double, ulint > & rhs) const {
        if (lhs.first == rhs.first) {
            return lhs.second < rhs.second;
        }
        return lhs.first < rhs.first;
    }
};

int main () {
    vector <ulint> vN;
    vN.push_back(50);
    vN.push_back(100);
    vN.push_back(200);
    vector <double> vD;
    vD.push_back(0.3);
    vD.push_back(0.5);
    vD.push_back(0.7);
    vector <ulint> vK;
    vK.push_back(0);
    vK.push_back(10);
    vK.push_back(20);
    vector <ulint> vT;
    vT.push_back(0);
    vT.push_back(1);
    vector <ulint> vI;
    vI.push_back(0);
    vI.push_back(1);
    vI.push_back(2);
    vector <double> vA;
    vA.push_back(0.3);
    vA.push_back(0.5);
    vA.push_back(0.7);

    set <ulint> sM;

    map <pair < double, ulint >, vector <double>, keyComp > vMN;
    map <pair < double, ulint >, vector <double>, keyComp > vObjValN;
    map <pair < double, ulint >, vector <double>, keyComp > vElapsedTimeN;
    map <pair < double, ulint >, vector <double>, keyComp > vNSolutionN;
    map <pair < double, ulint >, vector <double>, keyComp > vSumPenaltyN;
    map <pair < double, ulint >, vector <double>, keyComp > vSumEdgeCostN;
    map <pair < double, ulint >, vector <double>, keyComp > vObjValTimesElapsedTimeN;
    map <pair < double, ulint >, vector <double>, keyComp > vNSolutionDividedByNN;
    map <pair < double, ulint >, vector <double>, keyComp > vSumPenaltyDividedByObjValN;
    map <pair < double, ulint >, vector <double>, keyComp > vSumEdgeCostDividedByObjValN;
    for (vector <double> :: iterator itA = vA.begin(); itA != vA.end(); itA++) {
        double a = *itA;
        for (vector <ulint> :: iterator itN = vN.begin(); itN != vN.end(); itN++) {
            ulint n = *itN;
            vMN[make_pair(a, n)] = vector <double> ();
            vObjValN[make_pair(a, n)] = vector <double> ();
            vElapsedTimeN[make_pair(a, n)] = vector <double> ();
            vNSolutionN[make_pair(a, n)] = vector <double> ();
            vSumPenaltyN[make_pair(a, n)] = vector <double> ();
            vSumEdgeCostN[make_pair(a, n)] = vector <double> ();
            vObjValTimesElapsedTimeN[make_pair(a, n)] = vector <double> ();
            vNSolutionDividedByNN[make_pair(a, n)] = vector <double> ();
            vSumPenaltyDividedByObjValN[make_pair(a, n)] = vector <double> ();
            vSumEdgeCostDividedByObjValN[make_pair(a, n)] = vector <double> ();
        }
    }

    map <pair < double, ulint >, vector <double>, keyComp > vMD;
    map <pair < double, ulint >, vector <double>, keyComp > vObjValD;
    map <pair < double, ulint >, vector <double>, keyComp > vElapsedTimeD;
    map <pair < double, ulint >, vector <double>, keyComp > vNSolutionD;
    map <pair < double, ulint >, vector <double>, keyComp > vSumPenaltyD;
    map <pair < double, ulint >, vector <double>, keyComp > vSumEdgeCostD;
    map <pair < double, ulint >, vector <double>, keyComp > vObjValTimesElapsedTimeD;
    map <pair < double, ulint >, vector <double>, keyComp > vNSolutionDividedByND;
    map <pair < double, ulint >, vector <double>, keyComp > vSumPenaltyDividedByObjValD;
    map <pair < double, ulint >, vector <double>, keyComp > vSumEdgeCostDividedByObjValD;
    for (vector <double> :: iterator itA = vA.begin(); itA != vA.end(); itA++) {
        double a = *itA;
        for (vector <double> :: iterator itD = vD.begin(); itD != vD.end(); itD++) {
            double d = *itD;
            vMD[make_pair(a, d)] = vector <double> ();
            vObjValD[make_pair(a, d)] = vector <double> ();
            vElapsedTimeD[make_pair(a, d)] = vector <double> ();
            vNSolutionD[make_pair(a, d)] = vector <double> ();
            vSumPenaltyD[make_pair(a, d)] = vector <double> ();
            vSumEdgeCostD[make_pair(a, d)] = vector <double> ();
            vObjValTimesElapsedTimeD[make_pair(a, d)] = vector <double> ();
            vNSolutionDividedByND[make_pair(a, d)] = vector <double> ();
            vSumPenaltyDividedByObjValD[make_pair(a, d)] = vector <double> ();
            vSumEdgeCostDividedByObjValD[make_pair(a, d)] = vector <double> ();
        }
    }

    map <pair < double, ulint >, vector <double>, keyComp > vMK;
    map <pair < double, ulint >, vector <double>, keyComp > vObjValK;
    map <pair < double, ulint >, vector <double>, keyComp > vElapsedTimeK;
    map <pair < double, ulint >, vector <double>, keyComp > vNSolutionK;
    map <pair < double, ulint >, vector <double>, keyComp > vSumPenaltyK;
    map <pair < double, ulint >, vector <double>, keyComp > vSumEdgeCostK;
    map <pair < double, ulint >, vector <double>, keyComp > vObjValTimesElapsedTimeK;
    map <pair < double, ulint >, vector <double>, keyComp > vNSolutionDividedByNK;
    map <pair < double, ulint >, vector <double>, keyComp > vSumPenaltyDividedByObjValK;
    map <pair < double, ulint >, vector <double>, keyComp > vSumEdgeCostDividedByObjValK;
    for (vector <double> :: iterator itA = vA.begin(); itA != vA.end(); itA++) {
        double a = *itA;
        for (vector <ulint> :: iterator itK = vK.begin(); itK != vK.end(); itK++) {
            ulint k = *itK;
            vMK[make_pair(a, k)] = vector <double> ();
            vObjValK[make_pair(a, k)] = vector <double> ();
            vElapsedTimeK[make_pair(a, k)] = vector <double> ();
            vNSolutionK[make_pair(a, k)] = vector <double> ();
            vSumPenaltyK[make_pair(a, k)] = vector <double> ();
            vSumEdgeCostK[make_pair(a, k)] = vector <double> ();
            vObjValTimesElapsedTimeK[make_pair(a, k)] = vector <double> ();
            vNSolutionDividedByNK[make_pair(a, k)] = vector <double> ();
            vSumPenaltyDividedByObjValK[make_pair(a, k)] = vector <double> ();
            vSumEdgeCostDividedByObjValK[make_pair(a, k)] = vector <double> ();
        }
    }

    map <pair < double, ulint >, vector <double>, keyComp > vMT;
    map <pair < double, ulint >, vector <double>, keyComp > vObjValT;
    map <pair < double, ulint >, vector <double>, keyComp > vElapsedTimeT;
    map <pair < double, ulint >, vector <double>, keyComp > vNSolutionT;
    map <pair < double, ulint >, vector <double>, keyComp > vSumPenaltyT;
    map <pair < double, ulint >, vector <double>, keyComp > vSumEdgeCostT;
    map <pair < double, ulint >, vector <double>, keyComp > vObjValTimesElapsedTimeT;
    map <pair < double, ulint >, vector <double>, keyComp > vNSolutionDividedByNT;
    map <pair < double, ulint >, vector <double>, keyComp > vSumPenaltyDividedByObjValT;
    map <pair < double, ulint >, vector <double>, keyComp > vSumEdgeCostDividedByObjValT;
    for (vector <double> :: iterator itA = vA.begin(); itA != vA.end(); itA++) {
        double a = *itA;
        for (vector <ulint> :: iterator itT = vT.begin(); itT != vT.end(); itT++) {
            ulint t = *itT;
            vMT[make_pair(a, t)] = vector <double> ();
            vObjValT[make_pair(a, t)] = vector <double> ();
            vElapsedTimeT[make_pair(a, t)] = vector <double> ();
            vNSolutionT[make_pair(a, t)] = vector <double> ();
            vSumPenaltyT[make_pair(a, t)] = vector <double> ();
            vSumEdgeCostT[make_pair(a, t)] = vector <double> ();
            vObjValTimesElapsedTimeT[make_pair(a, t)] = vector <double> ();
            vNSolutionDividedByNT[make_pair(a, t)] = vector <double> ();
            vSumPenaltyDividedByObjValT[make_pair(a, t)] = vector <double> ();
            vSumEdgeCostDividedByObjValT[make_pair(a, t)] = vector <double> ();
        }
    }

    map <pair < double, ulint >, vector <double>, keyComp > vObjValM;
    map <pair < double, ulint >, vector <double>, keyComp > vElapsedTimeM;
    map <pair < double, ulint >, vector <double>, keyComp > vNSolutionM;
    map <pair < double, ulint >, vector <double>, keyComp > vSumPenaltyM;
    map <pair < double, ulint >, vector <double>, keyComp > vSumEdgeCostM;
    map <pair < double, ulint >, vector <double>, keyComp > vObjValTimesElapsedTimeM;
    map <pair < double, ulint >, vector <double>, keyComp > vNSolutionDividedByNM;
    map <pair < double, ulint >, vector <double>, keyComp > vSumPenaltyDividedByObjValM;
    map <pair < double, ulint >, vector <double>, keyComp > vSumEdgeCostDividedByObjValM;

    map <double, vector <double> > vMA;
    map <double, vector <double> > vObjValA;
    map <double, vector <double> > vElapsedTimeA;
    map <double, vector <double> > vNSolutionA;
    map <double, vector <double> > vSumPenaltyA;
    map <double, vector <double> > vSumEdgeCostA;
    map <double, vector <double> > vObjValTimesElapsedTimeA;
    map <double, vector <double> > vNSolutionDividedByNA;
    map <double, vector <double> > vSumPenaltyDividedByObjValA;
    map <double, vector <double> > vSumEdgeCostDividedByObjValA;
    for (vector <double> :: iterator itA = vA.begin(); itA != vA.end(); itA++) {
        double a = *itA;
        vMA[a] = vector <double> ();
        vObjValA[a] = vector <double> ();
        vElapsedTimeA[a] = vector <double> ();
        vNSolutionA[a] = vector <double> ();
        vSumPenaltyA[a] = vector <double> ();
        vSumEdgeCostA[a] = vector <double> ();
        vObjValTimesElapsedTimeA[a] = vector <double> ();
        vNSolutionDividedByNA[a] = vector <double> ();
        vSumPenaltyDividedByObjValA[a] = vector <double> ();
        vSumEdgeCostDividedByObjValA[a] = vector <double> ();
    }

    cout << "n,d,k,t,i,m,a,objVal,elapsedTime,nSolution,sumPenalty,sumEdgeCost,objVal*elapsedTime,nSolution/n,sumPenalty/objVal,sumEdgeCost/objVal" << endl;

    for (vector <ulint> :: iterator itN = vN.begin(); itN != vN.end(); itN++) {
        ulint n = *itN;
        stringstream ssN;
        ssN << n;
        string N = ssN.str();
        for (vector <double> :: iterator itD = vD.begin(); itD != vD.end(); itD++) {
            double d = *itD;
            stringstream ssD;
            ssD << fixed << setprecision(1) << d;
            string D = ssD.str();
            D.erase(remove(D.begin(), D.end(), '.'), D.end());
            for (vector <ulint> :: iterator itK = vK.begin(); itK != vK.end(); itK++) {
                ulint k = *itK;
                stringstream ssK;
                ssK << k;
                string K = ssK.str();
                for (vector <ulint> :: iterator itT = vT.begin(); itT != vT.end(); itT++) {
                    ulint t = *itT;
                    stringstream ssT;
                    ssT << t;
                    string T = ssT.str();
                    for (vector <ulint> :: iterator itI = vI.begin(); itI != vI.end(); itI++) {
                        ulint i = *itI;
                        stringstream ssI;
                        ssI << i;
                        string I = ssI.str();
                        for (vector <double> :: iterator itA = vA.begin(); itA != vA.end(); itA++) {
                            double a = *itA;
                            stringstream ssA;
                            ssA << fixed << setprecision(1) << a;
                            string A = ssA.str();
                            A.erase(remove(A.begin(), A.end(), '.'), A.end());

                            ulint m = 0;

                            double objVal = 0.0;
                            ifstream objValFile ("./output/N" + N + "D" + D + "K" + K + "T" + T + "I" + I + "A" + A + "/objVal.txt");
                            if (objValFile.is_open()) {
                                objValFile >> objVal;
                            }

                            ulint elapsedTime = 0.0;
                            ifstream elapsedTimeFile ("./output/N" + N + "D" + D + "K" + K + "T" + T + "I" + I + "A" + A + "/elapsedTime.txt");
                            if (elapsedTimeFile.is_open()) {
                                elapsedTimeFile >> elapsedTime;
                            }

                            double nSolution = 0.0, sumPenalty = 0.0, sumEdgeCost = 0.0;

                            ifstream inputFile ("../input/instanceN" + N + "D" + D + "K" + K + "T" + T + "I" + I + ".in");
                            ifstream resultFile ("./output/N" + N + "D" + D + "K" + K + "T" + T + "I" + I + "A" + A + "/result.out");

                            if (inputFile.is_open() && resultFile.is_open()) {
                                ulint n, mComplete, k, t, root;
                                double d;
                                inputFile >> n >> d >> k >> t >> mComplete >> m >> root;

                                vector < pair < pair <double, double>, ulint > > vertices (n, make_pair(make_pair(0, 0), 0));
                                // reading vertices' coordinates and penalty
                                for (ulint v = 0; v < n; v++) {
                                    inputFile >> vertices[v].first.first;
                                    inputFile >> vertices[v].first.second;
                                    inputFile >> vertices[v].second;
                                }

                                matrix W (n, vector <lint> (n, -1)); // adjacency matrix for the complete graph
                                for (ulint i = 0; i < n; i++) {
                                    W[i][i] = 0;
                                }

                                // reading (and ignoring) the complete graph's edges
                                for (ulint j = 0; j < mComplete; j++) {
                                    ulint u, v, w;
                                    inputFile >> u >> v >> w;
                                    W[u][v] = w;
                                    W[v][u] = w;
                                }

                                vector < pair <ulint, ulint> > edges (m, make_pair(0, 0));
                                map < pair <ulint, ulint>, ulint > weights;
                                // reading the graph's edges and its weights
                                for (ulint e = 0; e < m; e++) {
                                    ulint u, v, w;
                                    inputFile >> u >> v >> w;
                                    edges[e].first = u;
                                    edges[e].second = v;
                                    weights[edges[e]] = w;
                                }

                                ulint mSolution, costSolution;
                                resultFile >> nSolution >> mSolution >> costSolution;

                                vector <ulint> solutionVertices (nSolution, 0);
                                // reading the solution's vertices
                                for (ulint v = 0; v < nSolution; v++) {
                                    resultFile >> solutionVertices[v];
                                }

                                vector < pair <ulint, ulint> > solutionEdges (mSolution, make_pair(0, 0));
                                // reading the solution's edges
                                for (ulint e = 0; e < mSolution; e++) {
                                    resultFile >> solutionEdges[e].first;
                                    resultFile >> solutionEdges[e].second;
                                }

                                ulint allPenalties = 0;
                                for (ulint v = 0; v < n; v++) {
                                    allPenalties += vertices[v].second;
                                }

                                ulint chosenPenalties = 0;
                                for (ulint v = 0; v < nSolution; v++) {
                                    chosenPenalties += vertices[solutionVertices[v]].second;
                                }

                                sumPenalty = (allPenalties - chosenPenalties);

                                for (ulint e = 0; e < mSolution; e++) {
                                    sumEdgeCost += weights[solutionEdges[e]];
                                }
                            }

                            double objValTimesElapsedTime = (double) (objVal * ((double) elapsedTime));
                            double nSolutionDividedByN = (double) (((double) nSolution)/ ((double) n));
                            double sumPenaltyDividedByObjVal = sumPenalty / objVal;
                            double sumEdgeCostDividedByObjVal = sumEdgeCost / objVal;

                            sM.insert(m);

                            vMN[make_pair(a, n)].push_back((double) m);
                            vObjValN[make_pair(a, n)].push_back((double) objVal);
                            vElapsedTimeN[make_pair(a, n)].push_back((double) elapsedTime);
                            vNSolutionN[make_pair(a, n)].push_back((double) nSolution);
                            vSumPenaltyN[make_pair(a, n)].push_back((double) sumPenalty);
                            vSumEdgeCostN[make_pair(a, n)].push_back((double) sumEdgeCost);
                            vObjValTimesElapsedTimeN[make_pair(a, n)].push_back((double) objValTimesElapsedTime);
                            vNSolutionDividedByNN[make_pair(a, n)].push_back((double) nSolutionDividedByN);
                            vSumPenaltyDividedByObjValN[make_pair(a, n)].push_back((double) sumPenaltyDividedByObjVal);
                            vSumEdgeCostDividedByObjValN[make_pair(a, n)].push_back((double) sumEdgeCostDividedByObjVal);

                            vMD[make_pair(a, d)].push_back((double) m);
                            vObjValD[make_pair(a, d)].push_back((double) objVal);
                            vElapsedTimeD[make_pair(a, d)].push_back((double) elapsedTime);
                            vNSolutionD[make_pair(a, d)].push_back((double) nSolution);
                            vSumPenaltyD[make_pair(a, d)].push_back((double) sumPenalty);
                            vSumEdgeCostD[make_pair(a, d)].push_back((double) sumEdgeCost);
                            vObjValTimesElapsedTimeD[make_pair(a, d)].push_back((double) objValTimesElapsedTime);
                            vNSolutionDividedByND[make_pair(a, d)].push_back((double) nSolutionDividedByN);
                            vSumPenaltyDividedByObjValD[make_pair(a, d)].push_back((double) sumPenaltyDividedByObjVal);
                            vSumEdgeCostDividedByObjValD[make_pair(a, d)].push_back((double) sumEdgeCostDividedByObjVal);

                            vMK[make_pair(a, k)].push_back((double) m);
                            vObjValK[make_pair(a, k)].push_back((double) objVal);
                            vElapsedTimeK[make_pair(a, k)].push_back((double) elapsedTime);
                            vNSolutionK[make_pair(a, k)].push_back((double) nSolution);
                            vSumPenaltyK[make_pair(a, k)].push_back((double) sumPenalty);
                            vSumEdgeCostK[make_pair(a, k)].push_back((double) sumEdgeCost);
                            vObjValTimesElapsedTimeK[make_pair(a, k)].push_back((double) objValTimesElapsedTime);
                            vNSolutionDividedByNK[make_pair(a, k)].push_back((double) nSolutionDividedByN);
                            vSumPenaltyDividedByObjValK[make_pair(a, k)].push_back((double) sumPenaltyDividedByObjVal);
                            vSumEdgeCostDividedByObjValK[make_pair(a, k)].push_back((double) sumEdgeCostDividedByObjVal);

                            vMT[make_pair(a, t)].push_back((double) m);
                            vObjValT[make_pair(a, t)].push_back((double) objVal);
                            vElapsedTimeT[make_pair(a, t)].push_back((double) elapsedTime);
                            vNSolutionT[make_pair(a, t)].push_back((double) nSolution);
                            vSumPenaltyT[make_pair(a, t)].push_back((double) sumPenalty);
                            vSumEdgeCostT[make_pair(a, t)].push_back((double) sumEdgeCost);
                            vObjValTimesElapsedTimeT[make_pair(a, t)].push_back((double) objValTimesElapsedTime);
                            vNSolutionDividedByNT[make_pair(a, t)].push_back((double) nSolutionDividedByN);
                            vSumPenaltyDividedByObjValT[make_pair(a, t)].push_back((double) sumPenaltyDividedByObjVal);
                            vSumEdgeCostDividedByObjValT[make_pair(a, t)].push_back((double) sumEdgeCostDividedByObjVal);

                            if (vObjValM.find(make_pair(a, m)) == vObjValM.end()) {
                                vObjValM[make_pair(a, m)] = vector <double> ();
                            }
                            if (vElapsedTimeM.find(make_pair(a, m)) == vElapsedTimeM.end()) {
                                vElapsedTimeM[make_pair(a, m)] = vector <double> ();
                            }
                            if (vNSolutionM.find(make_pair(a, m)) == vNSolutionM.end()) {
                                vNSolutionM[make_pair(a, m)] = vector <double> ();
                            }
                            if (vSumPenaltyM.find(make_pair(a, m)) == vSumPenaltyM.end()) {
                                vSumPenaltyM[make_pair(a, m)] = vector <double> ();
                            }
                            if (vSumEdgeCostM.find(make_pair(a, m)) == vSumEdgeCostM.end()) {
                                vSumEdgeCostM[make_pair(a, m)] = vector <double> ();
                            }
                            if (vObjValTimesElapsedTimeM.find(make_pair(a, m)) == vObjValTimesElapsedTimeM.end()) {
                                vObjValTimesElapsedTimeM[make_pair(a, m)] = vector <double> ();
                            }
                            if (vNSolutionDividedByNM.find(make_pair(a, m)) == vNSolutionDividedByNM.end()) {
                                vNSolutionDividedByNM[make_pair(a, m)] = vector <double> ();
                            }
                            if (vSumPenaltyDividedByObjValM.find(make_pair(a, m)) == vSumPenaltyDividedByObjValM.end()) {
                                vSumPenaltyDividedByObjValM[make_pair(a, m)] = vector <double> ();
                            }
                            if (vSumEdgeCostDividedByObjValM.find(make_pair(a, m)) == vSumEdgeCostDividedByObjValM.end()) {
                                vSumEdgeCostDividedByObjValM[make_pair(a, m)] = vector <double> ();
                            }

                            vObjValM[make_pair(a, m)].push_back((double) objVal);
                            vElapsedTimeM[make_pair(a, m)].push_back((double) elapsedTime);
                            vNSolutionM[make_pair(a, m)].push_back((double) nSolution);
                            vSumPenaltyM[make_pair(a, m)].push_back((double) sumPenalty);
                            vSumEdgeCostM[make_pair(a, m)].push_back((double) sumEdgeCost);
                            vObjValTimesElapsedTimeM[make_pair(a, m)].push_back((double) objValTimesElapsedTime);
                            vNSolutionDividedByNM[make_pair(a, m)].push_back((double) nSolutionDividedByN);
                            vSumPenaltyDividedByObjValM[make_pair(a, m)].push_back((double) sumPenaltyDividedByObjVal);
                            vSumEdgeCostDividedByObjValM[make_pair(a, m)].push_back((double) sumEdgeCostDividedByObjVal);

                            vMA[a].push_back((double) m);
                            vObjValA[a].push_back((double) objVal);
                            vElapsedTimeA[a].push_back((double) elapsedTime);
                            vNSolutionA[a].push_back((double) nSolution);
                            vSumPenaltyA[a].push_back((double) sumPenalty);
                            vSumEdgeCostA[a].push_back((double) sumEdgeCost);
                            vObjValTimesElapsedTimeA[a].push_back((double) objValTimesElapsedTime);
                            vNSolutionDividedByNA[a].push_back((double) nSolutionDividedByN);
                            vSumPenaltyDividedByObjValA[a].push_back((double) sumPenaltyDividedByObjVal);
                            vSumEdgeCostDividedByObjValA[a].push_back((double) sumEdgeCostDividedByObjVal);

                            cout << n << ',' << d << ',' << k << ',' << t << ',' << i << ',' << m << ',' << a << ',' << objVal << ',' << elapsedTime << ',' << nSolution << ',' << sumPenalty << ',' << sumEdgeCost << ',' << objValTimesElapsedTime << ',' << nSolutionDividedByN << ',' << sumPenaltyDividedByObjVal << ',' << sumEdgeCostDividedByObjVal << endl;
                        }
                    }
                }
            }
        }
    }

    ofstream solutionNFile ("./output/solutionN.csv", ofstream :: out);
    solutionNFile << "a,n,m,deltam,objVal,deltaObjVal,elapsedTime,deltaElapsedTime,nSolution,deltaNSolution,sumPenalty,deltaSumPenalty,sumEdgeCost,deltaSumEdgeCost,objValTimesElapsedTime,deltaObjValTimesElapsedTime,nSolutionDividedByN,deltaNSolutionDividedByN,sumPenaltyDividedByObjVal,deltaSumPenaltyDividedByObjVal,sumEdgeCostDividedByObjVal,deltaSumEdgeCostDividedByObjVal" << endl;
    for (vector <double> :: iterator itA = vA.begin(); itA != vA.end(); itA++) {
        double a = *itA;
        for (vector <ulint> :: iterator itN = vN.begin(); itN != vN.end(); itN++) {
            ulint n = *itN;
            solutionNFile << a << ',' << n << ',' << average(vMN[make_pair(a, n)]) << ',' << stdev(vMN[make_pair(a, n)]) << ',' << average(vObjValN[make_pair(a, n)]) << ',' << stdev(vObjValN[make_pair(a, n)]) << ',' << average(vElapsedTimeN[make_pair(a, n)]) << ',' << stdev(vElapsedTimeN[make_pair(a, n)]) << ',' << average(vNSolutionN[make_pair(a, n)]) << ',' << stdev(vNSolutionN[make_pair(a, n)]) << ',' << average(vSumPenaltyN[make_pair(a, n)]) << ',' << stdev(vSumPenaltyN[make_pair(a, n)]) << ',' << average(vSumEdgeCostN[make_pair(a, n)]) << ',' << stdev(vSumEdgeCostN[make_pair(a, n)]) << ',' << average(vObjValTimesElapsedTimeN[make_pair(a, n)]) << ',' << stdev(vObjValTimesElapsedTimeN[make_pair(a, n)]) << ',' << average(vNSolutionDividedByNN[make_pair(a, n)]) << ',' << stdev(vNSolutionDividedByNN[make_pair(a, n)]) << ',' << average(vSumPenaltyDividedByObjValN[make_pair(a, n)]) << ',' << stdev(vSumPenaltyDividedByObjValN[make_pair(a, n)]) << ',' << average(vSumEdgeCostDividedByObjValN[make_pair(a, n)]) << ',' << stdev(vSumEdgeCostDividedByObjValN[make_pair(a, n)]) << endl;
        }
    }
    solutionNFile.close();

    ofstream solutionDFile ("./output/solutionD.csv", ofstream :: out);
    solutionDFile << "a,n,m,deltam,objVal,deltaObjVal,elapsedTime,deltaElapsedTime,nSolution,deltaNSolution,sumPenalty,deltaSumPenalty,sumEdgeCost,deltaSumEdgeCost,objValTimesElapsedTime,deltaObjValTimesElapsedTime,nSolutionDividedByN,deltaNSolutionDividedByN,sumPenaltyDividedByObjVal,deltaSumPenaltyDividedByObjVal,sumEdgeCostDividedByObjVal,deltaSumEdgeCostDividedByObjVal" << endl;
    for (vector <double> :: iterator itA = vA.begin(); itA != vA.end(); itA++) {
        double a = *itA;
        for (vector <double> :: iterator itD = vD.begin(); itD != vD.end(); itD++) {
            double d = *itD;
            solutionDFile << a << ',' << d << ',' << average(vMD[make_pair(a, d)]) << ',' << stdev(vMD[make_pair(a, d)]) << ',' << average(vObjValD[make_pair(a, d)]) << ',' << stdev(vObjValD[make_pair(a, d)]) << ',' << average(vElapsedTimeD[make_pair(a, d)]) << ',' << stdev(vElapsedTimeD[make_pair(a, d)]) << ',' << average(vNSolutionD[make_pair(a, d)]) << ',' << stdev(vNSolutionD[make_pair(a, d)]) << ',' << average(vSumPenaltyD[make_pair(a, d)]) << ',' << stdev(vSumPenaltyD[make_pair(a, d)]) << ',' << average(vSumEdgeCostD[make_pair(a, d)]) << ',' << stdev(vSumEdgeCostD[make_pair(a, d)]) << ',' << average(vObjValTimesElapsedTimeD[make_pair(a, d)]) << ',' << stdev(vObjValTimesElapsedTimeD[make_pair(a, d)]) << ',' << average(vNSolutionDividedByND[make_pair(a, d)]) << ',' << stdev(vNSolutionDividedByND[make_pair(a, d)]) << ',' << average(vSumPenaltyDividedByObjValD[make_pair(a, d)]) << ',' << stdev(vSumPenaltyDividedByObjValD[make_pair(a, d)]) << ',' << average(vSumEdgeCostDividedByObjValD[make_pair(a, d)]) << ',' << stdev(vSumEdgeCostDividedByObjValD[make_pair(a, d)]) << endl;
        }
    }
    solutionDFile.close();

    ofstream solutionKFile ("./output/solutionK.csv", ofstream :: out);
    solutionKFile << "a,k,m,deltam,objVal,deltaObjVal,elapsedTime,deltaElapsedTime,nSolution,deltaNSolution,sumPenalty,deltaSumPenalty,sumEdgeCost,deltaSumEdgeCost,objValTimesElapsedTime,deltaObjValTimesElapsedTime,nSolutionDividedByN,deltaNSolutionDividedByN,sumPenaltyDividedByObjVal,deltaSumPenaltyDividedByObjVal,sumEdgeCostDividedByObjVal,deltaSumEdgeCostDividedByObjVal" << endl;
    for (vector <double> :: iterator itA = vA.begin(); itA != vA.end(); itA++) {
        double a = *itA;
        for (vector <ulint> :: iterator itK = vK.begin(); itK != vK.end(); itK++) {
            ulint k = *itK;
            solutionKFile << a << ',' << k << ',' << average(vMK[make_pair(a, k)]) << ',' << stdev(vMK[make_pair(a, k)]) << ',' << average(vObjValK[make_pair(a, k)]) << ',' << stdev(vObjValK[make_pair(a, k)]) << ',' << average(vElapsedTimeK[make_pair(a, k)]) << ',' << stdev(vElapsedTimeK[make_pair(a, k)]) << ',' << average(vNSolutionK[make_pair(a, k)]) << ',' << stdev(vNSolutionK[make_pair(a, k)]) << ',' << average(vSumPenaltyK[make_pair(a, k)]) << ',' << stdev(vSumPenaltyK[make_pair(a, k)]) << ',' << average(vSumEdgeCostK[make_pair(a, k)]) << ',' << stdev(vSumEdgeCostK[make_pair(a, k)]) << ',' << average(vObjValTimesElapsedTimeK[make_pair(a, k)]) << ',' << stdev(vObjValTimesElapsedTimeK[make_pair(a, k)]) << ',' << average(vNSolutionDividedByNK[make_pair(a, k)]) << ',' << stdev(vNSolutionDividedByNK[make_pair(a, k)]) << ',' << average(vSumPenaltyDividedByObjValK[make_pair(a, k)]) << ',' << stdev(vSumPenaltyDividedByObjValK[make_pair(a, k)]) << ',' << average(vSumEdgeCostDividedByObjValK[make_pair(a, k)]) << ',' << stdev(vSumEdgeCostDividedByObjValK[make_pair(a, k)]) << endl;
        }
    }
    solutionKFile.close();

    ofstream solutionTFile ("./output/solutionT.csv", ofstream :: out);
    solutionTFile << "a,t,m,deltam,objVal,deltaObjVal,elapsedTime,deltaElapsedTime,nSolution,deltaNSolution,sumPenalty,deltaSumPenalty,sumEdgeCost,deltaSumEdgeCost,objValTimesElapsedTime,deltaObjValTimesElapsedTime,nSolutionDividedByN,deltaNSolutionDividedByN,sumPenaltyDividedByObjVal,deltaSumPenaltyDividedByObjVal,sumEdgeCostDividedByObjVal,deltaSumEdgeCostDividedByObjVal" << endl;
    for (vector <double> :: iterator itA = vA.begin(); itA != vA.end(); itA++) {
        double a = *itA;
        for (vector <ulint> :: iterator itT = vT.begin(); itT != vT.end(); itT++) {
            ulint t = *itT;
            solutionTFile << a << ',' << t << ',' << average(vMT[make_pair(a, t)]) << ',' << stdev(vMT[make_pair(a, t)]) << ',' << average(vObjValT[make_pair(a, t)]) << ',' << stdev(vObjValT[make_pair(a, t)]) << ',' << average(vElapsedTimeT[make_pair(a, t)]) << ',' << stdev(vElapsedTimeT[make_pair(a, t)]) << ',' << average(vNSolutionT[make_pair(a, t)]) << ',' << stdev(vNSolutionT[make_pair(a, t)]) << ',' << average(vSumPenaltyT[make_pair(a, t)]) << ',' << stdev(vSumPenaltyT[make_pair(a, t)]) << ',' << average(vSumEdgeCostT[make_pair(a, t)]) << ',' << stdev(vSumEdgeCostT[make_pair(a, t)]) << ',' << average(vObjValTimesElapsedTimeT[make_pair(a, t)]) << ',' << stdev(vObjValTimesElapsedTimeT[make_pair(a, t)]) << ',' << average(vNSolutionDividedByNT[make_pair(a, t)]) << ',' << stdev(vNSolutionDividedByNT[make_pair(a, t)]) << ',' << average(vSumPenaltyDividedByObjValT[make_pair(a, t)]) << ',' << stdev(vSumPenaltyDividedByObjValT[make_pair(a, t)]) << ',' << average(vSumEdgeCostDividedByObjValT[make_pair(a, t)]) << ',' << stdev(vSumEdgeCostDividedByObjValT[make_pair(a, t)]) << endl;
        }
    }
    solutionTFile.close();

    ofstream solutionMFile ("./output/solutionM.csv", ofstream :: out);
    solutionMFile << "a,m,objVal,deltaObjVal,elapsedTime,deltaElapsedTime,nSolution,deltaNSolution,sumPenalty,deltaSumPenalty,sumEdgeCost,deltaSumEdgeCost,objValTimesElapsedTime,deltaObjValTimesElapsedTime,nSolutionDividedByN,deltaNSolutionDividedByN,sumPenaltyDividedByObjVal,deltaSumPenaltyDividedByObjVal,sumEdgeCostDividedByObjVal,deltaSumEdgeCostDividedByObjVal" << endl;
    for (vector <double> :: iterator itA = vA.begin(); itA != vA.end(); itA++) {
        double a = *itA;
        for (set <ulint> :: iterator itM = sM.begin(); itM != sM.end(); itM++) {
            ulint m = *itM;
            solutionMFile << a << ',' << m << ',' << average(vObjValM[make_pair(a, m)]) << ',' << stdev(vObjValM[make_pair(a, m)]) << ',' << average(vElapsedTimeM[make_pair(a, m)]) << ',' << stdev(vElapsedTimeM[make_pair(a, m)]) << ',' << average(vNSolutionM[make_pair(a, m)]) << ',' << stdev(vNSolutionM[make_pair(a, m)]) << ',' << average(vSumPenaltyM[make_pair(a, m)]) << ',' << stdev(vSumPenaltyM[make_pair(a, m)]) << ',' << average(vSumEdgeCostM[make_pair(a, m)]) << ',' << stdev(vSumEdgeCostM[make_pair(a, m)]) << ',' << average(vObjValTimesElapsedTimeM[make_pair(a, m)]) << ',' << stdev(vObjValTimesElapsedTimeM[make_pair(a, m)]) << ',' << average(vNSolutionDividedByNM[make_pair(a, m)]) << ',' << stdev(vNSolutionDividedByNM[make_pair(a, m)]) << ',' << average(vSumPenaltyDividedByObjValM[make_pair(a, m)]) << ',' << stdev(vSumPenaltyDividedByObjValM[make_pair(a, m)]) << ',' << average(vSumEdgeCostDividedByObjValM[make_pair(a, m)]) << ',' << stdev(vSumEdgeCostDividedByObjValM[make_pair(a, m)]) << endl;
        }
    }
    solutionMFile.close();

    ofstream solutionAFile ("./output/solutionA.csv", ofstream :: out);
    solutionAFile << "a,m,deltam,objVal,deltaObjVal,elapsedTime,deltaElapsedTime,nSolution,deltaNSolution,sumPenalty,deltaSumPenalty,sumEdgeCost,deltaSumEdgeCost,objValTimesElapsedTime,deltaObjValTimesElapsedTime,nSolutionDividedByN,deltaNSolutionDividedByN,sumPenaltyDividedByObjVal,deltaSumPenaltyDividedByObjVal,sumEdgeCostDividedByObjVal,deltaSumEdgeCostDividedByObjVal" << endl;
    for (vector <double> :: iterator itA = vA.begin(); itA != vA.end(); itA++) {
        double a = *itA;
        solutionAFile << a << ',' << average(vMA[a]) << ',' << stdev(vMA[a]) << ',' << average(vObjValA[a]) << ',' << stdev(vObjValA[a]) << ',' << average(vElapsedTimeA[a]) << ',' << stdev(vElapsedTimeA[a]) << ',' << average(vNSolutionA[a]) << ',' << stdev(vNSolutionA[a]) << ',' << average(vSumPenaltyA[a]) << ',' << stdev(vSumPenaltyA[a]) << ',' << average(vSumEdgeCostA[a]) << ',' << stdev(vSumEdgeCostA[a]) << ',' << average(vObjValTimesElapsedTimeA[a]) << ',' << stdev(vObjValTimesElapsedTimeA[a]) << ',' << average(vNSolutionDividedByNA[a]) << ',' << stdev(vNSolutionDividedByNA[a]) << ',' << average(vSumPenaltyDividedByObjValA[a]) << ',' << stdev(vSumPenaltyDividedByObjValA[a]) << ',' << average(vSumEdgeCostDividedByObjValA[a]) << ',' << stdev(vSumEdgeCostDividedByObjValA[a]) << endl;
    }
    solutionAFile.close();

    return 0;
}
