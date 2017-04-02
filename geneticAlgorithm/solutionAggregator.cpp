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
    bool operator() (const pair < pair <ulint, double>, ulint > & lhs, const pair < pair <ulint, double>, ulint > & rhs) const {
        if (lhs.first.first == rhs.first.first) {
            if (lhs.first.second == rhs.first.second) {
                return lhs.second < rhs.second;
            }
            return lhs.first.second < rhs.first.second;
        }
        return lhs.first.first < rhs.first.first;
    }
};

struct keyComp2 {
    bool operator() (const pair <ulint, double> & lhs, const pair <ulint, double> & rhs) const {
        if (lhs.first == rhs.first) {
            return lhs.second < rhs.second;
        }
        return lhs.first < rhs.first;
    }
};

int main () {
    vector <ulint> vPS;
    vPS.push_back(10);
    vPS.push_back(50);
    vPS.push_back(100);
    vector <double> vMR;
    vMR.push_back(0.1);
    vMR.push_back(0.2);
    vMR.push_back(0.3);
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

    set <ulint> sM;

    map <pair < pair <ulint, double>, ulint >, vector <double>, keyComp > vMN;
    map <pair < pair <ulint, double>, ulint >, vector <double>, keyComp > vObjValN;
    map <pair < pair <ulint, double>, ulint >, vector <double>, keyComp > vElapsedTimeN;
    map <pair < pair <ulint, double>, ulint >, vector <double>, keyComp > vNSolutionN;
    map <pair < pair <ulint, double>, ulint >, vector <double>, keyComp > vSumPenaltyN;
    map <pair < pair <ulint, double>, ulint >, vector <double>, keyComp > vSumEdgeCostN;
    map <pair < pair <ulint, double>, ulint >, vector <double>, keyComp > vObjValTimesElapsedTimeN;
    map <pair < pair <ulint, double>, ulint >, vector <double>, keyComp > vNSolutionDividedByNN;
    map <pair < pair <ulint, double>, ulint >, vector <double>, keyComp > vSumPenaltyDividedByObjValN;
    map <pair < pair <ulint, double>, ulint >, vector <double>, keyComp > vSumEdgeCostDividedByObjValN;
    for (vector <ulint> :: iterator itPS = vPS.begin(); itPS != vPS.end(); itPS++) {
        ulint ps = *itPS;
        for (vector <double> :: iterator itMR = vMR.begin(); itMR != vMR.end(); itMR++) {
            double mr = *itMR;
            for (vector <ulint> :: iterator itN = vN.begin(); itN != vN.end(); itN++) {
                ulint n = *itN;
                vMN[make_pair(make_pair(ps, mr), n)] = vector <double> ();
                vObjValN[make_pair(make_pair(ps, mr), n)] = vector <double> ();
                vElapsedTimeN[make_pair(make_pair(ps, mr), n)] = vector <double> ();
                vNSolutionN[make_pair(make_pair(ps, mr), n)] = vector <double> ();
                vSumPenaltyN[make_pair(make_pair(ps, mr), n)] = vector <double> ();
                vSumEdgeCostN[make_pair(make_pair(ps, mr), n)] = vector <double> ();
                vObjValTimesElapsedTimeN[make_pair(make_pair(ps, mr), n)] = vector <double> ();
                vNSolutionDividedByNN[make_pair(make_pair(ps, mr), n)] = vector <double> ();
                vSumPenaltyDividedByObjValN[make_pair(make_pair(ps, mr), n)] = vector <double> ();
                vSumEdgeCostDividedByObjValN[make_pair(make_pair(ps, mr), n)] = vector <double> ();
            }
        }
    }

    map <pair < pair <ulint, double>, double >, vector <double>, keyComp > vMD;
    map <pair < pair <ulint, double>, double >, vector <double>, keyComp > vObjValD;
    map <pair < pair <ulint, double>, double >, vector <double>, keyComp > vElapsedTimeD;
    map <pair < pair <ulint, double>, double >, vector <double>, keyComp > vNSolutionD;
    map <pair < pair <ulint, double>, double >, vector <double>, keyComp > vSumPenaltyD;
    map <pair < pair <ulint, double>, double >, vector <double>, keyComp > vSumEdgeCostD;
    map <pair < pair <ulint, double>, double >, vector <double>, keyComp > vObjValTimesElapsedTimeD;
    map <pair < pair <ulint, double>, double >, vector <double>, keyComp > vNSolutionDividedByND;
    map <pair < pair <ulint, double>, double >, vector <double>, keyComp > vSumPenaltyDividedByObjValD;
    map <pair < pair <ulint, double>, double >, vector <double>, keyComp > vSumEdgeCostDividedByObjValD;
    for (vector <ulint> :: iterator itPS = vPS.begin(); itPS != vPS.end(); itPS++) {
        ulint ps = *itPS;
        for (vector <double> :: iterator itMR = vMR.begin(); itMR != vMR.end(); itMR++) {
            double mr = *itMR;
            for (vector <double> :: iterator itD = vD.begin(); itD != vD.end(); itD++) {
                double d = *itD;
                vMD[make_pair(make_pair(ps, mr), d)] = vector <double> ();
                vObjValD[make_pair(make_pair(ps, mr), d)] = vector <double> ();
                vElapsedTimeD[make_pair(make_pair(ps, mr), d)] = vector <double> ();
                vNSolutionD[make_pair(make_pair(ps, mr), d)] = vector <double> ();
                vSumPenaltyD[make_pair(make_pair(ps, mr), d)] = vector <double> ();
                vSumEdgeCostD[make_pair(make_pair(ps, mr), d)] = vector <double> ();
                vObjValTimesElapsedTimeD[make_pair(make_pair(ps, mr), d)] = vector <double> ();
                vNSolutionDividedByND[make_pair(make_pair(ps, mr), d)] = vector <double> ();
                vSumPenaltyDividedByObjValD[make_pair(make_pair(ps, mr), d)] = vector <double> ();
                vSumEdgeCostDividedByObjValD[make_pair(make_pair(ps, mr), d)] = vector <double> ();
            }
        }
    }

    map <pair < pair <ulint, double>, ulint >, vector <double>, keyComp > vMK;
    map <pair < pair <ulint, double>, ulint >, vector <double>, keyComp > vObjValK;
    map <pair < pair <ulint, double>, ulint >, vector <double>, keyComp > vElapsedTimeK;
    map <pair < pair <ulint, double>, ulint >, vector <double>, keyComp > vNSolutionK;
    map <pair < pair <ulint, double>, ulint >, vector <double>, keyComp > vSumPenaltyK;
    map <pair < pair <ulint, double>, ulint >, vector <double>, keyComp > vSumEdgeCostK;
    map <pair < pair <ulint, double>, ulint >, vector <double>, keyComp > vObjValTimesElapsedTimeK;
    map <pair < pair <ulint, double>, ulint >, vector <double>, keyComp > vNSolutionDividedByNK;
    map <pair < pair <ulint, double>, ulint >, vector <double>, keyComp > vSumPenaltyDividedByObjValK;
    map <pair < pair <ulint, double>, ulint >, vector <double>, keyComp > vSumEdgeCostDividedByObjValK;
    for (vector <ulint> :: iterator itPS = vPS.begin(); itPS != vPS.end(); itPS++) {
        ulint ps = *itPS;
        for (vector <double> :: iterator itMR = vMR.begin(); itMR != vMR.end(); itMR++) {
            double mr = *itMR;
            for (vector <ulint> :: iterator itK = vK.begin(); itK != vK.end(); itK++) {
                ulint k = *itK;
                vMK[make_pair(make_pair(ps, mr), k)] = vector <double> ();
                vObjValK[make_pair(make_pair(ps, mr), k)] = vector <double> ();
                vElapsedTimeK[make_pair(make_pair(ps, mr), k)] = vector <double> ();
                vNSolutionK[make_pair(make_pair(ps, mr), k)] = vector <double> ();
                vSumPenaltyK[make_pair(make_pair(ps, mr), k)] = vector <double> ();
                vSumEdgeCostK[make_pair(make_pair(ps, mr), k)] = vector <double> ();
                vObjValTimesElapsedTimeK[make_pair(make_pair(ps, mr), k)] = vector <double> ();
                vNSolutionDividedByNK[make_pair(make_pair(ps, mr), k)] = vector <double> ();
                vSumPenaltyDividedByObjValK[make_pair(make_pair(ps, mr), k)] = vector <double> ();
                vSumEdgeCostDividedByObjValK[make_pair(make_pair(ps, mr), k)] = vector <double> ();
            }
        }
    }

    map <pair < pair <ulint, double>, ulint >, vector <double>, keyComp > vMT;
    map <pair < pair <ulint, double>, ulint >, vector <double>, keyComp > vObjValT;
    map <pair < pair <ulint, double>, ulint >, vector <double>, keyComp > vElapsedTimeT;
    map <pair < pair <ulint, double>, ulint >, vector <double>, keyComp > vNSolutionT;
    map <pair < pair <ulint, double>, ulint >, vector <double>, keyComp > vSumPenaltyT;
    map <pair < pair <ulint, double>, ulint >, vector <double>, keyComp > vSumEdgeCostT;
    map <pair < pair <ulint, double>, ulint >, vector <double>, keyComp > vObjValTimesElapsedTimeT;
    map <pair < pair <ulint, double>, ulint >, vector <double>, keyComp > vNSolutionDividedByNT;
    map <pair < pair <ulint, double>, ulint >, vector <double>, keyComp > vSumPenaltyDividedByObjValT;
    map <pair < pair <ulint, double>, ulint >, vector <double>, keyComp > vSumEdgeCostDividedByObjValT;
    for (vector <ulint> :: iterator itPS = vPS.begin(); itPS != vPS.end(); itPS++) {
        ulint ps = *itPS;
        for (vector <double> :: iterator itMR = vMR.begin(); itMR != vMR.end(); itMR++) {
            double mr = *itMR;
            for (vector <ulint> :: iterator itT = vT.begin(); itT != vT.end(); itT++) {
                ulint t = *itT;
                vMT[make_pair(make_pair(ps, mr), t)] = vector <double> ();
                vObjValT[make_pair(make_pair(ps, mr), t)] = vector <double> ();
                vElapsedTimeT[make_pair(make_pair(ps, mr), t)] = vector <double> ();
                vNSolutionT[make_pair(make_pair(ps, mr), t)] = vector <double> ();
                vSumPenaltyT[make_pair(make_pair(ps, mr), t)] = vector <double> ();
                vSumEdgeCostT[make_pair(make_pair(ps, mr), t)] = vector <double> ();
                vObjValTimesElapsedTimeT[make_pair(make_pair(ps, mr), t)] = vector <double> ();
                vNSolutionDividedByNT[make_pair(make_pair(ps, mr), t)] = vector <double> ();
                vSumPenaltyDividedByObjValT[make_pair(make_pair(ps, mr), t)] = vector <double> ();
                vSumEdgeCostDividedByObjValT[make_pair(make_pair(ps, mr), t)] = vector <double> ();
            }
        }
    }

    map <pair < pair <ulint, double>, ulint >, vector <double>, keyComp > vObjValM;
    map <pair < pair <ulint, double>, ulint >, vector <double>, keyComp > vElapsedTimeM;
    map <pair < pair <ulint, double>, ulint >, vector <double>, keyComp > vNSolutionM;
    map <pair < pair <ulint, double>, ulint >, vector <double>, keyComp > vSumPenaltyM;
    map <pair < pair <ulint, double>, ulint >, vector <double>, keyComp > vSumEdgeCostM;
    map <pair < pair <ulint, double>, ulint >, vector <double>, keyComp > vObjValTimesElapsedTimeM;
    map <pair < pair <ulint, double>, ulint >, vector <double>, keyComp > vNSolutionDividedByNM;
    map <pair < pair <ulint, double>, ulint >, vector <double>, keyComp > vSumPenaltyDividedByObjValM;
    map <pair < pair <ulint, double>, ulint >, vector <double>, keyComp > vSumEdgeCostDividedByObjValM;

    map <pair <ulint, double>, vector <double>, keyComp2> vMPSMR;
    map <pair <ulint, double>, vector <double>, keyComp2> vObjValPSMR;
    map <pair <ulint, double>, vector <double>, keyComp2> vElapsedTimePSMR;
    map <pair <ulint, double>, vector <double>, keyComp2> vNSolutionPSMR;
    map <pair <ulint, double>, vector <double>, keyComp2> vSumPenaltyPSMR;
    map <pair <ulint, double>, vector <double>, keyComp2> vSumEdgeCostPSMR;
    map <pair <ulint, double>, vector <double>, keyComp2> vObjValTimesElapsedTimePSMR;
    map <pair <ulint, double>, vector <double>, keyComp2> vNSolutionDividedByNPSMR;
    map <pair <ulint, double>, vector <double>, keyComp2> vSumPenaltyDividedByObjValPSMR;
    map <pair <ulint, double>, vector <double>, keyComp2> vSumEdgeCostDividedByObjValPSMR;
    for (vector <ulint> :: iterator itPS = vPS.begin(); itPS != vPS.end(); itPS++) {
        ulint ps = *itPS;
        for (vector <double> :: iterator itMR = vMR.begin(); itMR != vMR.end(); itMR++) {
            double mr = *itMR;
            vMPSMR[make_pair(ps, mr)] = vector <double> ();
            vObjValPSMR[make_pair(ps, mr)] = vector <double> ();
            vElapsedTimePSMR[make_pair(ps, mr)] = vector <double> ();
            vNSolutionPSMR[make_pair(ps, mr)] = vector <double> ();
            vSumPenaltyPSMR[make_pair(ps, mr)] = vector <double> ();
            vSumEdgeCostPSMR[make_pair(ps, mr)] = vector <double> ();
            vObjValTimesElapsedTimePSMR[make_pair(ps, mr)] = vector <double> ();
            vNSolutionDividedByNPSMR[make_pair(ps, mr)] = vector <double> ();
            vSumPenaltyDividedByObjValPSMR[make_pair(ps, mr)] = vector <double> ();
            vSumEdgeCostDividedByObjValPSMR[make_pair(ps, mr)] = vector <double> ();
        }
    }

    cout << "ps,mr,n,d,k,t,i,m,objVal,elapsedTime,nSolution,sumPenalty,sumEdgeCost,objVal*elapsedTime,nSolution/n,sumPenalty/objVal,sumEdgeCost/objVal" << endl;

    for (vector <ulint> :: iterator itPS = vPS.begin(); itPS != vPS.end(); itPS++) {
        ulint ps = *itPS;
        stringstream ssPS;
        ssPS << ps;
        string PS = ssPS.str();
        for (vector <double> :: iterator itMR = vMR.begin(); itMR != vMR.end(); itMR++) {
            double mr = *itMR;
            stringstream ssMR;
            ssMR << fixed << setprecision(1) << mr;
            string MR = ssMR.str();
            MR.erase(remove(MR.begin(), MR.end(), '.'), MR.end());
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

                                ulint m = 0;

                                double objVal = 0.0;
                                ifstream objValFile ("./output/N" + N + "D" + D + "K" + K + "T" + T + "I" + I + "PS" + PS + "MR" + MR + "/objVal.txt");
                                if (objValFile.is_open()) {
                                    objValFile >> objVal;
                                }

                                ulint elapsedTime = 0.0;
                                ifstream elapsedTimeFile ("./output/N" + N + "D" + D + "K" + K + "T" + T + "I" + I + "PS" + PS + "MR" + MR + "/elapsedTime.txt");
                                if (elapsedTimeFile.is_open()) {
                                    elapsedTimeFile >> elapsedTime;
                                }

                                double nSolution = 0.0, sumPenalty = 0.0, sumEdgeCost = 0.0;

                                ifstream inputFile ("../input/instanceN" + N + "D" + D + "K" + K + "T" + T + "I" + I + ".in");
                                ifstream resultFile ("./output/N" + N + "D" + D + "K" + K + "T" + T + "I" + I + "PS" + PS + "MR" + MR + "/result.out");

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

                                    vector <int> isInSolution (n, 0);

                                    vector <ulint> solutionVertices (nSolution, 0);
                                    // reading the solution's vertices
                                    for (ulint v = 0; v < nSolution; v++) {
                                        resultFile >> solutionVertices[v];
                                        isInSolution[solutionVertices[v]] = 1;
                                    }

                                    vector < pair <ulint, ulint> > solutionEdges (mSolution, make_pair(0, 0));
                                    // reading the solution's edges
                                    for (ulint e = 0; e < mSolution; e++) {
                                        resultFile >> solutionEdges[e].first;
                                        resultFile >> solutionEdges[e].second;
                                    }

                                    for (ulint v = 0; v < n; v++) {
                                        sumPenalty += (1 - isInSolution[v]) * vertices[v].second;
                                    }

                                    for (ulint e = 0; e < mSolution; e++) {
                                        sumEdgeCost += weights[solutionEdges[e]];
                                    }
                                }

                                double objValTimesElapsedTime = (double) (objVal * ((double) elapsedTime));
                                double nSolutionDividedByN = (double) (((double) nSolution)/ ((double) n));
                                double sumPenaltyDividedByObjVal = sumPenalty / objVal;
                                double sumEdgeCostDividedByObjVal = sumEdgeCost / objVal;

                                sM.insert(m);

                                vMN[make_pair(make_pair(ps, mr), n)].push_back((double) m);
                                vObjValN[make_pair(make_pair(ps, mr), n)].push_back((double) objVal);
                                vElapsedTimeN[make_pair(make_pair(ps, mr), n)].push_back((double) elapsedTime);
                                vNSolutionN[make_pair(make_pair(ps, mr), n)].push_back((double) nSolution);
                                vSumPenaltyN[make_pair(make_pair(ps, mr), n)].push_back((double) sumPenalty);
                                vSumEdgeCostN[make_pair(make_pair(ps, mr), n)].push_back((double) sumEdgeCost);
                                vObjValTimesElapsedTimeN[make_pair(make_pair(ps, mr), n)].push_back((double) objValTimesElapsedTime);
                                vNSolutionDividedByNN[make_pair(make_pair(ps, mr), n)].push_back((double) nSolutionDividedByN);
                                vSumPenaltyDividedByObjValN[make_pair(make_pair(ps, mr), n)].push_back((double) sumPenaltyDividedByObjVal);
                                vSumEdgeCostDividedByObjValN[make_pair(make_pair(ps, mr), n)].push_back((double) sumEdgeCostDividedByObjVal);

                                vMD[make_pair(make_pair(ps, mr), d)].push_back((double) m);
                                vObjValD[make_pair(make_pair(ps, mr), d)].push_back((double) objVal);
                                vElapsedTimeD[make_pair(make_pair(ps, mr), d)].push_back((double) elapsedTime);
                                vNSolutionD[make_pair(make_pair(ps, mr), d)].push_back((double) nSolution);
                                vSumPenaltyD[make_pair(make_pair(ps, mr), d)].push_back((double) sumPenalty);
                                vSumEdgeCostD[make_pair(make_pair(ps, mr), d)].push_back((double) sumEdgeCost);
                                vObjValTimesElapsedTimeD[make_pair(make_pair(ps, mr), d)].push_back((double) objValTimesElapsedTime);
                                vNSolutionDividedByND[make_pair(make_pair(ps, mr), d)].push_back((double) nSolutionDividedByN);
                                vSumPenaltyDividedByObjValD[make_pair(make_pair(ps, mr), d)].push_back((double) sumPenaltyDividedByObjVal);
                                vSumEdgeCostDividedByObjValD[make_pair(make_pair(ps, mr), d)].push_back((double) sumEdgeCostDividedByObjVal);

                                vMK[make_pair(make_pair(ps, mr), k)].push_back((double) m);
                                vObjValK[make_pair(make_pair(ps, mr), k)].push_back((double) objVal);
                                vElapsedTimeK[make_pair(make_pair(ps, mr), k)].push_back((double) elapsedTime);
                                vNSolutionK[make_pair(make_pair(ps, mr), k)].push_back((double) nSolution);
                                vSumPenaltyK[make_pair(make_pair(ps, mr), k)].push_back((double) sumPenalty);
                                vSumEdgeCostK[make_pair(make_pair(ps, mr), k)].push_back((double) sumEdgeCost);
                                vObjValTimesElapsedTimeK[make_pair(make_pair(ps, mr), k)].push_back((double) objValTimesElapsedTime);
                                vNSolutionDividedByNK[make_pair(make_pair(ps, mr), k)].push_back((double) nSolutionDividedByN);
                                vSumPenaltyDividedByObjValK[make_pair(make_pair(ps, mr), k)].push_back((double) sumPenaltyDividedByObjVal);
                                vSumEdgeCostDividedByObjValK[make_pair(make_pair(ps, mr), k)].push_back((double) sumEdgeCostDividedByObjVal);

                                vMT[make_pair(make_pair(ps, mr), t)].push_back((double) m);
                                vObjValT[make_pair(make_pair(ps, mr), t)].push_back((double) objVal);
                                vElapsedTimeT[make_pair(make_pair(ps, mr), t)].push_back((double) elapsedTime);
                                vNSolutionT[make_pair(make_pair(ps, mr), t)].push_back((double) nSolution);
                                vSumPenaltyT[make_pair(make_pair(ps, mr), t)].push_back((double) sumPenalty);
                                vSumEdgeCostT[make_pair(make_pair(ps, mr), t)].push_back((double) sumEdgeCost);
                                vObjValTimesElapsedTimeT[make_pair(make_pair(ps, mr), t)].push_back((double) objValTimesElapsedTime);
                                vNSolutionDividedByNT[make_pair(make_pair(ps, mr), t)].push_back((double) nSolutionDividedByN);
                                vSumPenaltyDividedByObjValT[make_pair(make_pair(ps, mr), t)].push_back((double) sumPenaltyDividedByObjVal);
                                vSumEdgeCostDividedByObjValT[make_pair(make_pair(ps, mr), t)].push_back((double) sumEdgeCostDividedByObjVal);

                                if (vObjValM.find(make_pair(make_pair(ps, mr), m)) == vObjValM.end()) {
                                    vObjValM[make_pair(make_pair(ps, mr), m)] = vector <double> ();
                                }
                                if (vElapsedTimeM.find(make_pair(make_pair(ps, mr), m)) == vElapsedTimeM.end()) {
                                    vElapsedTimeM[make_pair(make_pair(ps, mr), m)] = vector <double> ();
                                }
                                if (vNSolutionM.find(make_pair(make_pair(ps, mr), m)) == vNSolutionM.end()) {
                                    vNSolutionM[make_pair(make_pair(ps, mr), m)] = vector <double> ();
                                }
                                if (vSumPenaltyM.find(make_pair(make_pair(ps, mr), m)) == vSumPenaltyM.end()) {
                                    vSumPenaltyM[make_pair(make_pair(ps, mr), m)] = vector <double> ();
                                }
                                if (vSumEdgeCostM.find(make_pair(make_pair(ps, mr), m)) == vSumEdgeCostM.end()) {
                                    vSumEdgeCostM[make_pair(make_pair(ps, mr), m)] = vector <double> ();
                                }
                                if (vObjValTimesElapsedTimeM.find(make_pair(make_pair(ps, mr), m)) == vObjValTimesElapsedTimeM.end()) {
                                    vObjValTimesElapsedTimeM[make_pair(make_pair(ps, mr), m)] = vector <double> ();
                                }
                                if (vNSolutionDividedByNM.find(make_pair(make_pair(ps, mr), m)) == vNSolutionDividedByNM.end()) {
                                    vNSolutionDividedByNM[make_pair(make_pair(ps, mr), m)] = vector <double> ();
                                }
                                if (vSumPenaltyDividedByObjValM.find(make_pair(make_pair(ps, mr), m)) == vSumPenaltyDividedByObjValM.end()) {
                                    vSumPenaltyDividedByObjValM[make_pair(make_pair(ps, mr), m)] = vector <double> ();
                                }
                                if (vSumEdgeCostDividedByObjValM.find(make_pair(make_pair(ps, mr), m)) == vSumEdgeCostDividedByObjValM.end()) {
                                    vSumEdgeCostDividedByObjValM[make_pair(make_pair(ps, mr), m)] = vector <double> ();
                                }

                                vObjValM[make_pair(make_pair(ps, mr), m)].push_back((double) objVal);
                                vElapsedTimeM[make_pair(make_pair(ps, mr), m)].push_back((double) elapsedTime);
                                vNSolutionM[make_pair(make_pair(ps, mr), m)].push_back((double) nSolution);
                                vSumPenaltyM[make_pair(make_pair(ps, mr), m)].push_back((double) sumPenalty);
                                vSumEdgeCostM[make_pair(make_pair(ps, mr), m)].push_back((double) sumEdgeCost);
                                vObjValTimesElapsedTimeM[make_pair(make_pair(ps, mr), m)].push_back((double) objValTimesElapsedTime);
                                vNSolutionDividedByNM[make_pair(make_pair(ps, mr), m)].push_back((double) nSolutionDividedByN);
                                vSumPenaltyDividedByObjValM[make_pair(make_pair(ps, mr), m)].push_back((double) sumPenaltyDividedByObjVal);
                                vSumEdgeCostDividedByObjValM[make_pair(make_pair(ps, mr), m)].push_back((double) sumEdgeCostDividedByObjVal);

                                vMPSMR[make_pair(ps, mr)].push_back((double) m);
                                vObjValPSMR[make_pair(ps, mr)].push_back((double) objVal);
                                vElapsedTimePSMR[make_pair(ps, mr)].push_back((double) elapsedTime);
                                vNSolutionPSMR[make_pair(ps, mr)].push_back((double) nSolution);
                                vSumPenaltyPSMR[make_pair(ps, mr)].push_back((double) sumPenalty);
                                vSumEdgeCostPSMR[make_pair(ps, mr)].push_back((double) sumEdgeCost);
                                vObjValTimesElapsedTimePSMR[make_pair(ps, mr)].push_back((double) objValTimesElapsedTime);
                                vNSolutionDividedByNPSMR[make_pair(ps, mr)].push_back((double) nSolutionDividedByN);
                                vSumPenaltyDividedByObjValPSMR[make_pair(ps, mr)].push_back((double) sumPenaltyDividedByObjVal);
                                vSumEdgeCostDividedByObjValPSMR[make_pair(ps, mr)].push_back((double) sumEdgeCostDividedByObjVal);

                                cout << ps << ',' << mr << ',' << n << ',' << d << ',' << k << ',' << t << ',' << i << ',' << m << ',' << objVal << ',' << elapsedTime << ',' << nSolution << ',' << sumPenalty << ',' << sumEdgeCost << ',' << objValTimesElapsedTime << ',' << nSolutionDividedByN << ',' << sumPenaltyDividedByObjVal << ',' << sumEdgeCostDividedByObjVal << endl;
                            }
                        }
                    }
                }
            }
        }
    }

    ofstream solutionNFile ("./output/solutionN.csv", ofstream :: out);
    solutionNFile << "ps,mr,n,m,deltam,objVal,deltaObjVal,elapsedTime,deltaElapsedTime,nSolution,deltaNSolution,sumPenalty,deltaSumPenalty,sumEdgeCost,deltaSumEdgeCost,objValTimesElapsedTime,deltaObjValTimesElapsedTime,nSolutionDividedByN,deltaNSolutionDividedByN,sumPenaltyDividedByObjVal,deltaSumPenaltyDividedByObjVal,sumEdgeCostDividedByObjVal,deltaSumEdgeCostDividedByObjVal" << endl;
    for (vector <ulint> :: iterator itPS = vPS.begin(); itPS != vPS.end(); itPS++) {
        ulint ps = *itPS;
        for (vector <double> :: iterator itMR = vMR.begin(); itMR != vMR.end(); itMR++) {
            double mr = *itMR;
            for (vector <ulint> :: iterator itN = vN.begin(); itN != vN.end(); itN++) {
                ulint n = *itN;
                solutionNFile << ps << ',' << mr << ',' << n << ',' << average(vMN[make_pair(make_pair(ps, mr), n)]) << ',' << stdev(vMN[make_pair(make_pair(ps, mr), n)]) << ',' << average(vObjValN[make_pair(make_pair(ps, mr), n)]) << ',' << stdev(vObjValN[make_pair(make_pair(ps, mr), n)]) << ',' << average(vElapsedTimeN[make_pair(make_pair(ps, mr), n)]) << ',' << stdev(vElapsedTimeN[make_pair(make_pair(ps, mr), n)]) << ',' << average(vNSolutionN[make_pair(make_pair(ps, mr), n)]) << ',' << stdev(vNSolutionN[make_pair(make_pair(ps, mr), n)]) << ',' << average(vSumPenaltyN[make_pair(make_pair(ps, mr), n)]) << ',' << stdev(vSumPenaltyN[make_pair(make_pair(ps, mr), n)]) << ',' << average(vSumEdgeCostN[make_pair(make_pair(ps, mr), n)]) << ',' << stdev(vSumEdgeCostN[make_pair(make_pair(ps, mr), n)]) << ',' << average(vObjValTimesElapsedTimeN[make_pair(make_pair(ps, mr), n)]) << ',' << stdev(vObjValTimesElapsedTimeN[make_pair(make_pair(ps, mr), n)]) << ',' << average(vNSolutionDividedByNN[make_pair(make_pair(ps, mr), n)]) << ',' << stdev(vNSolutionDividedByNN[make_pair(make_pair(ps, mr), n)]) << ',' << average(vSumPenaltyDividedByObjValN[make_pair(make_pair(ps, mr), n)]) << ',' << stdev(vSumPenaltyDividedByObjValN[make_pair(make_pair(ps, mr), n)]) << ',' << average(vSumEdgeCostDividedByObjValN[make_pair(make_pair(ps, mr), n)]) << ',' << stdev(vSumEdgeCostDividedByObjValN[make_pair(make_pair(ps, mr), n)]) << endl;
            }
        }
    }
    solutionNFile.close();

    ofstream solutionDFile ("./output/solutionD.csv", ofstream :: out);
    solutionDFile << "ps,mr,d,m,deltam,objVal,deltaObjVal,elapsedTime,deltaElapsedTime,nSolution,deltaNSolution,sumPenalty,deltaSumPenalty,sumEdgeCost,deltaSumEdgeCost,objValTimesElapsedTime,deltaObjValTimesElapsedTime,nSolutionDividedByN,deltaNSolutionDividedByN,sumPenaltyDividedByObjVal,deltaSumPenaltyDividedByObjVal,sumEdgeCostDividedByObjVal,deltaSumEdgeCostDividedByObjVal" << endl;
    for (vector <ulint> :: iterator itPS = vPS.begin(); itPS != vPS.end(); itPS++) {
        ulint ps = *itPS;
        for (vector <double> :: iterator itMR = vMR.begin(); itMR != vMR.end(); itMR++) {
            double mr = *itMR;
            for (vector <double> :: iterator itD = vD.begin(); itD != vD.end(); itD++) {
                double d = *itD;
                solutionDFile << ps << ',' << mr << ',' << d << ',' << average(vMD[make_pair(make_pair(ps, mr), d)]) << ',' << stdev(vMD[make_pair(make_pair(ps, mr), d)]) << ',' << average(vObjValD[make_pair(make_pair(ps, mr), d)]) << ',' << stdev(vObjValD[make_pair(make_pair(ps, mr), d)]) << ',' << average(vElapsedTimeD[make_pair(make_pair(ps, mr), d)]) << ',' << stdev(vElapsedTimeD[make_pair(make_pair(ps, mr), d)]) << ',' << average(vNSolutionD[make_pair(make_pair(ps, mr), d)]) << ',' << stdev(vNSolutionD[make_pair(make_pair(ps, mr), d)]) << ',' << average(vSumPenaltyD[make_pair(make_pair(ps, mr), d)]) << ',' << stdev(vSumPenaltyD[make_pair(make_pair(ps, mr), d)]) << ',' << average(vSumEdgeCostD[make_pair(make_pair(ps, mr), d)]) << ',' << stdev(vSumEdgeCostD[make_pair(make_pair(ps, mr), d)]) << ',' << average(vObjValTimesElapsedTimeD[make_pair(make_pair(ps, mr), d)]) << ',' << stdev(vObjValTimesElapsedTimeD[make_pair(make_pair(ps, mr), d)]) << ',' << average(vNSolutionDividedByND[make_pair(make_pair(ps, mr), d)]) << ',' << stdev(vNSolutionDividedByND[make_pair(make_pair(ps, mr), d)]) << ',' << average(vSumPenaltyDividedByObjValD[make_pair(make_pair(ps, mr), d)]) << ',' << stdev(vSumPenaltyDividedByObjValD[make_pair(make_pair(ps, mr), d)]) << ',' << average(vSumEdgeCostDividedByObjValD[make_pair(make_pair(ps, mr), d)]) << ',' << stdev(vSumEdgeCostDividedByObjValD[make_pair(make_pair(ps, mr), d)]) << endl;
            }
        }
    }
    solutionDFile.close();

    ofstream solutionKFile ("./output/solutionK.csv", ofstream :: out);
    solutionKFile << "ps,mr,k,m,deltam,objVal,deltaObjVal,elapsedTime,deltaElapsedTime,nSolution,deltaNSolution,sumPenalty,deltaSumPenalty,sumEdgeCost,deltaSumEdgeCost,objValTimesElapsedTime,deltaObjValTimesElapsedTime,nSolutionDividedByN,deltaNSolutionDividedByN,sumPenaltyDividedByObjVal,deltaSumPenaltyDividedByObjVal,sumEdgeCostDividedByObjVal,deltaSumEdgeCostDividedByObjVal" << endl;
    for (vector <ulint> :: iterator itPS = vPS.begin(); itPS != vPS.end(); itPS++) {
        ulint ps = *itPS;
        for (vector <double> :: iterator itMR = vMR.begin(); itMR != vMR.end(); itMR++) {
            double mr = *itMR;
            for (vector <ulint> :: iterator itK = vK.begin(); itK != vK.end(); itK++) {
                ulint k = *itK;
                solutionKFile << ps << ',' << mr << ',' << k << ',' << average(vMK[make_pair(make_pair(ps, mr), k)]) << ',' << stdev(vMK[make_pair(make_pair(ps, mr), k)]) << ',' << average(vObjValK[make_pair(make_pair(ps, mr), k)]) << ',' << stdev(vObjValK[make_pair(make_pair(ps, mr), k)]) << ',' << average(vElapsedTimeK[make_pair(make_pair(ps, mr), k)]) << ',' << stdev(vElapsedTimeK[make_pair(make_pair(ps, mr), k)]) << ',' << average(vNSolutionK[make_pair(make_pair(ps, mr), k)]) << ',' << stdev(vNSolutionK[make_pair(make_pair(ps, mr), k)]) << ',' << average(vSumPenaltyK[make_pair(make_pair(ps, mr), k)]) << ',' << stdev(vSumPenaltyK[make_pair(make_pair(ps, mr), k)]) << ',' << average(vSumEdgeCostK[make_pair(make_pair(ps, mr), k)]) << ',' << stdev(vSumEdgeCostK[make_pair(make_pair(ps, mr), k)]) << ',' << average(vObjValTimesElapsedTimeK[make_pair(make_pair(ps, mr), k)]) << ',' << stdev(vObjValTimesElapsedTimeK[make_pair(make_pair(ps, mr), k)]) << ',' << average(vNSolutionDividedByNK[make_pair(make_pair(ps, mr), k)]) << ',' << stdev(vNSolutionDividedByNK[make_pair(make_pair(ps, mr), k)]) << ',' << average(vSumPenaltyDividedByObjValK[make_pair(make_pair(ps, mr), k)]) << ',' << stdev(vSumPenaltyDividedByObjValK[make_pair(make_pair(ps, mr), k)]) << ',' << average(vSumEdgeCostDividedByObjValK[make_pair(make_pair(ps, mr), k)]) << ',' << stdev(vSumEdgeCostDividedByObjValK[make_pair(make_pair(ps, mr), k)]) << endl;
            }
        }
    }
    solutionKFile.close();

    ofstream solutionTFile ("./output/solutionT.csv", ofstream :: out);
    solutionTFile << "ps,mr,t,m,deltam,objVal,deltaObjVal,elapsedTime,deltaElapsedTime,nSolution,deltaNSolution,sumPenalty,deltaSumPenalty,sumEdgeCost,deltaSumEdgeCost,objValTimesElapsedTime,deltaObjValTimesElapsedTime,nSolutionDividedByN,deltaNSolutionDividedByN,sumPenaltyDividedByObjVal,deltaSumPenaltyDividedByObjVal,sumEdgeCostDividedByObjVal,deltaSumEdgeCostDividedByObjVal" << endl;
    for (vector <ulint> :: iterator itPS = vPS.begin(); itPS != vPS.end(); itPS++) {
        ulint ps = *itPS;
        for (vector <double> :: iterator itMR = vMR.begin(); itMR != vMR.end(); itMR++) {
            double mr = *itMR;
            for (vector <ulint> :: iterator itT = vT.begin(); itT != vT.end(); itT++) {
                ulint t = *itT;
                solutionTFile << ps << ',' << mr << ',' << t << ',' << average(vMT[make_pair(make_pair(ps, mr), t)]) << ',' << stdev(vMT[make_pair(make_pair(ps, mr), t)]) << ',' << average(vObjValT[make_pair(make_pair(ps, mr), t)]) << ',' << stdev(vObjValT[make_pair(make_pair(ps, mr), t)]) << ',' << average(vElapsedTimeT[make_pair(make_pair(ps, mr), t)]) << ',' << stdev(vElapsedTimeT[make_pair(make_pair(ps, mr), t)]) << ',' << average(vNSolutionT[make_pair(make_pair(ps, mr), t)]) << ',' << stdev(vNSolutionT[make_pair(make_pair(ps, mr), t)]) << ',' << average(vSumPenaltyT[make_pair(make_pair(ps, mr), t)]) << ',' << stdev(vSumPenaltyT[make_pair(make_pair(ps, mr), t)]) << ',' << average(vSumEdgeCostT[make_pair(make_pair(ps, mr), t)]) << ',' << stdev(vSumEdgeCostT[make_pair(make_pair(ps, mr), t)]) << ',' << average(vObjValTimesElapsedTimeT[make_pair(make_pair(ps, mr), t)]) << ',' << stdev(vObjValTimesElapsedTimeT[make_pair(make_pair(ps, mr), t)]) << ',' << average(vNSolutionDividedByNT[make_pair(make_pair(ps, mr), t)]) << ',' << stdev(vNSolutionDividedByNT[make_pair(make_pair(ps, mr), t)]) << ',' << average(vSumPenaltyDividedByObjValT[make_pair(make_pair(ps, mr), t)]) << ',' << stdev(vSumPenaltyDividedByObjValT[make_pair(make_pair(ps, mr), t)]) << ',' << average(vSumEdgeCostDividedByObjValT[make_pair(make_pair(ps, mr), t)]) << ',' << stdev(vSumEdgeCostDividedByObjValT[make_pair(make_pair(ps, mr), t)]) << endl;
            }
        }
    }
    solutionTFile.close();

    ofstream solutionMFile ("./output/solutionM.csv", ofstream :: out);
    solutionMFile << "ps,mr,m,objVal,deltaObjVal,elapsedTime,deltaElapsedTime,nSolution,deltaNSolution,sumPenalty,deltaSumPenalty,sumEdgeCost,deltaSumEdgeCost,objValTimesElapsedTime,deltaObjValTimesElapsedTime,nSolutionDividedByN,deltaNSolutionDividedByN,sumPenaltyDividedByObjVal,deltaSumPenaltyDividedByObjVal,sumEdgeCostDividedByObjVal,deltaSumEdgeCostDividedByObjVal" << endl;
    for (vector <ulint> :: iterator itPS = vPS.begin(); itPS != vPS.end(); itPS++) {
        ulint ps = *itPS;
        for (vector <double> :: iterator itMR = vMR.begin(); itMR != vMR.end(); itMR++) {
            double mr = *itMR;
            for (set <ulint> :: iterator itM = sM.begin(); itM != sM.end(); itM++) {
                ulint m = *itM;
                solutionMFile << ps << ',' << mr << ',' << m << ',' << average(vObjValM[make_pair(make_pair(ps, mr), m)]) << ',' << stdev(vObjValM[make_pair(make_pair(ps, mr), m)]) << ',' << average(vElapsedTimeM[make_pair(make_pair(ps, mr), m)]) << ',' << stdev(vElapsedTimeM[make_pair(make_pair(ps, mr), m)]) << ',' << average(vNSolutionM[make_pair(make_pair(ps, mr), m)]) << ',' << stdev(vNSolutionM[make_pair(make_pair(ps, mr), m)]) << ',' << average(vSumPenaltyM[make_pair(make_pair(ps, mr), m)]) << ',' << stdev(vSumPenaltyM[make_pair(make_pair(ps, mr), m)]) << ',' << average(vSumEdgeCostM[make_pair(make_pair(ps, mr), m)]) << ',' << stdev(vSumEdgeCostM[make_pair(make_pair(ps, mr), m)]) << ',' << average(vObjValTimesElapsedTimeM[make_pair(make_pair(ps, mr), m)]) << ',' << stdev(vObjValTimesElapsedTimeM[make_pair(make_pair(ps, mr), m)]) << ',' << average(vNSolutionDividedByNM[make_pair(make_pair(ps, mr), m)]) << ',' << stdev(vNSolutionDividedByNM[make_pair(make_pair(ps, mr), m)]) << ',' << average(vSumPenaltyDividedByObjValM[make_pair(make_pair(ps, mr), m)]) << ',' << stdev(vSumPenaltyDividedByObjValM[make_pair(make_pair(ps, mr), m)]) << ',' << average(vSumEdgeCostDividedByObjValM[make_pair(make_pair(ps, mr), m)]) << ',' << stdev(vSumEdgeCostDividedByObjValM[make_pair(make_pair(ps, mr), m)]) << endl;
            }
        }
    }
    solutionMFile.close();

    ofstream solutionPSMRFile ("./output/solutionPSMR.csv", ofstream :: out);
    solutionPSMRFile << "ps,mr,m,deltam,objVal,deltaObjVal,elapsedTime,deltaElapsedTime,nSolution,deltaNSolution,sumPenalty,deltaSumPenalty,sumEdgeCost,deltaSumEdgeCost,objValTimesElapsedTime,deltaObjValTimesElapsedTime,nSolutionDividedByN,deltaNSolutionDividedByN,sumPenaltyDividedByObjVal,deltaSumPenaltyDividedByObjVal,sumEdgeCostDividedByObjVal,deltaSumEdgeCostDividedByObjVal" << endl;
    for (vector <ulint> :: iterator itPS = vPS.begin(); itPS != vPS.end(); itPS++) {
        ulint ps = *itPS;
        for (vector <double> :: iterator itMR = vMR.begin(); itMR != vMR.end(); itMR++) {
            double mr = *itMR;
            solutionPSMRFile << ps << ',' << mr << ',' << average(vMPSMR[make_pair(ps, mr)]) << ',' << stdev(vMPSMR[make_pair(ps, mr)]) << ',' << average(vObjValPSMR[make_pair(ps, mr)]) << ',' << stdev(vObjValPSMR[make_pair(ps, mr)]) << ',' << average(vElapsedTimePSMR[make_pair(ps, mr)]) << ',' << stdev(vElapsedTimePSMR[make_pair(ps, mr)]) << ',' << average(vNSolutionPSMR[make_pair(ps, mr)]) << ',' << stdev(vNSolutionPSMR[make_pair(ps, mr)]) << ',' << average(vSumPenaltyPSMR[make_pair(ps, mr)]) << ',' << stdev(vSumPenaltyPSMR[make_pair(ps, mr)]) << ',' << average(vSumEdgeCostPSMR[make_pair(ps, mr)]) << ',' << stdev(vSumEdgeCostPSMR[make_pair(ps, mr)]) << ',' << average(vObjValTimesElapsedTimePSMR[make_pair(ps, mr)]) << ',' << stdev(vObjValTimesElapsedTimePSMR[make_pair(ps, mr)]) << ',' << average(vNSolutionDividedByNPSMR[make_pair(ps, mr)]) << ',' << stdev(vNSolutionDividedByNPSMR[make_pair(ps, mr)]) << ',' << average(vSumPenaltyDividedByObjValPSMR[make_pair(ps, mr)]) << ',' << stdev(vSumPenaltyDividedByObjValPSMR[make_pair(ps, mr)]) << ',' << average(vSumEdgeCostDividedByObjValPSMR[make_pair(ps, mr)]) << ',' << stdev(vSumEdgeCostDividedByObjValPSMR[make_pair(ps, mr)]) << endl;
        }
    }
    solutionPSMRFile.close();

    return 0;
}
