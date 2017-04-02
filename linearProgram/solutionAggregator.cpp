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

    set <ulint> sM;

    map <ulint, vector <double> > vMN;
    map <ulint, vector <double> > vFlagN;
    map <ulint, vector <double> > vObjValN;
    map <ulint, vector <double> > vGapN;
    map <ulint, vector <double> > vElapsedTimeN;
    map <ulint, vector <double> > vNSolutionN;
    map <ulint, vector <double> > vSumPenaltyN;
    map <ulint, vector <double> > vSumEdgeCostN;
    map <ulint, vector <double> > vObjValTimesElapsedTimeN;
    map <ulint, vector <double> > vNSolutionDividedByNN;
    map <ulint, vector <double> > vSumPenaltyDividedByObjValN;
    map <ulint, vector <double> > vSumEdgeCostDividedByObjValN;
    for (vector <ulint> :: iterator itN = vN.begin(); itN != vN.end(); itN++) {
        ulint n = *itN;
        vMN[n] = vector <double> ();
        vFlagN[n] = vector <double> ();
        vObjValN[n] = vector <double> ();
        vGapN[n] = vector <double> ();
        vElapsedTimeN[n] = vector <double> ();
        vNSolutionN[n] = vector <double> ();
        vSumPenaltyN[n] = vector <double> ();
        vSumEdgeCostN[n] = vector <double> ();
        vObjValTimesElapsedTimeN[n] = vector <double> ();
        vNSolutionDividedByNN[n] = vector <double> ();
        vSumPenaltyDividedByObjValN[n] = vector <double> ();
        vSumEdgeCostDividedByObjValN[n] = vector <double> ();
    }

    map <double, vector <double> > vMD;
    map <double, vector <double> > vFlagD;
    map <double, vector <double> > vObjValD;
    map <double, vector <double> > vGapD;
    map <double, vector <double> > vElapsedTimeD;
    map <double, vector <double> > vNSolutionD;
    map <double, vector <double> > vSumPenaltyD;
    map <double, vector <double> > vSumEdgeCostD;
    map <double, vector <double> > vObjValTimesElapsedTimeD;
    map <double, vector <double> > vNSolutionDividedByND;
    map <double, vector <double> > vSumPenaltyDividedByObjValD;
    map <double, vector <double> > vSumEdgeCostDividedByObjValD;
    for (vector <double> :: iterator itD = vD.begin(); itD != vD.end(); itD++) {
        double d = *itD;
        vMD[d] = vector <double> ();
        vFlagD[d] = vector <double> ();
        vObjValD[d] = vector <double> ();
        vGapD[d] = vector <double> ();
        vElapsedTimeD[d] = vector <double> ();
        vNSolutionD[d] = vector <double> ();
        vSumPenaltyD[d] = vector <double> ();
        vSumEdgeCostD[d] = vector <double> ();
        vObjValTimesElapsedTimeD[d] = vector <double> ();
        vNSolutionDividedByND[d] = vector <double> ();
        vSumPenaltyDividedByObjValD[d] = vector <double> ();
        vSumEdgeCostDividedByObjValD[d] = vector <double> ();
    }

    map <ulint, vector <double> > vMK;
    map <ulint, vector <double> > vFlagK;
    map <ulint, vector <double> > vObjValK;
    map <ulint, vector <double> > vGapK;
    map <ulint, vector <double> > vElapsedTimeK;
    map <ulint, vector <double> > vNSolutionK;
    map <ulint, vector <double> > vSumPenaltyK;
    map <ulint, vector <double> > vSumEdgeCostK;
    map <ulint, vector <double> > vObjValTimesElapsedTimeK;
    map <ulint, vector <double> > vNSolutionDividedByNK;
    map <ulint, vector <double> > vSumPenaltyDividedByObjValK;
    map <ulint, vector <double> > vSumEdgeCostDividedByObjValK;
    for (vector <ulint> :: iterator itK = vK.begin(); itK != vK.end(); itK++) {
        ulint k = *itK;
        vMK[k] = vector <double> ();
        vFlagK[k] = vector <double> ();
        vGapK[k] = vector <double> ();
        vObjValK[k] = vector <double> ();
        vElapsedTimeK[k] = vector <double> ();
        vNSolutionK[k] = vector <double> ();
        vSumPenaltyK[k] = vector <double> ();
        vSumEdgeCostK[k] = vector <double> ();
        vObjValTimesElapsedTimeK[k] = vector <double> ();
        vNSolutionDividedByNK[k] = vector <double> ();
        vSumPenaltyDividedByObjValK[k] = vector <double> ();
        vSumEdgeCostDividedByObjValK[k] = vector <double> ();
    }

    map <ulint, vector <double> > vMT;
    map <ulint, vector <double> > vFlagT;
    map <ulint, vector <double> > vObjValT;
    map <ulint, vector <double> > vGapT;
    map <ulint, vector <double> > vElapsedTimeT;
    map <ulint, vector <double> > vNSolutionT;
    map <ulint, vector <double> > vSumPenaltyT;
    map <ulint, vector <double> > vSumEdgeCostT;
    map <ulint, vector <double> > vObjValTimesElapsedTimeT;
    map <ulint, vector <double> > vNSolutionDividedByNT;
    map <ulint, vector <double> > vSumPenaltyDividedByObjValT;
    map <ulint, vector <double> > vSumEdgeCostDividedByObjValT;
    for (vector <ulint> :: iterator itT = vT.begin(); itT != vT.end(); itT++) {
        ulint t = *itT;
        vMT[t] = vector <double> ();
        vFlagT[t] = vector <double> ();
        vObjValT[t] = vector <double> ();
        vGapT[t] = vector <double> ();
        vElapsedTimeT[t] = vector <double> ();
        vNSolutionT[t] = vector <double> ();
        vSumPenaltyT[t] = vector <double> ();
        vSumEdgeCostT[t] = vector <double> ();
        vObjValTimesElapsedTimeT[t] = vector <double> ();
        vNSolutionDividedByNT[t] = vector <double> ();
        vSumPenaltyDividedByObjValT[t] = vector <double> ();
        vSumEdgeCostDividedByObjValT[t] = vector <double> ();
    }

    map <ulint, vector <double> > vFlagM;
    map <ulint, vector <double> > vObjValM;
    map <ulint, vector <double> > vGapM;
    map <ulint, vector <double> > vElapsedTimeM;
    map <ulint, vector <double> > vNSolutionM;
    map <ulint, vector <double> > vSumPenaltyM;
    map <ulint, vector <double> > vSumEdgeCostM;
    map <ulint, vector <double> > vObjValTimesElapsedTimeM;
    map <ulint, vector <double> > vNSolutionDividedByNM;
    map <ulint, vector <double> > vSumPenaltyDividedByObjValM;
    map <ulint, vector <double> > vSumEdgeCostDividedByObjValM;

    cout << "n,d,k,t,i,m,flag,objVal,gap,elapsedTime,nSolution,sumPenalty,sumEdgeCost,objVal*elapsedTime,nSolution/n,sumPenalty/objVal,sumEdgeCost/objVal" << endl;

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
                        ifstream objValFile ("./output/N" + N + "D" + D + "K" + K + "T" + T + "I" + I + "/objVal.txt");
                        if (objValFile.is_open()) {
                            objValFile >> objVal;
                        }

                        double gap = 0.0;
                        ifstream gapFile ("./output/N" + N + "D" + D + "K" + K + "T" + T + "I" + I + "/gap.txt");
                        if (gapFile.is_open()) {
                            gapFile >> gap;
                        }

                        ulint elapsedTime = 0.0;
                        ifstream elapsedTimeFile ("./output/N" + N + "D" + D + "K" + K + "T" + T + "I" + I + "/elapsedTime.txt");
                        if (elapsedTimeFile.is_open()) {
                            elapsedTimeFile >> elapsedTime;
                        }

                        double nSolution = 0.0, sumPenalty = 0.0, sumEdgeCost = 0.0;

                        ifstream inputFile ("../input/instanceN" + N + "D" + D + "K" + K + "T" + T + "I" + I + ".in");
                        ifstream resultFile ("./output/N" + N + "D" + D + "K" + K + "T" + T + "I" + I + "/result.out");

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

                        ulint flag = 0;

                        if (objVal < 1E100) {
                            flag = 1;
                        }

                        vFlagN[n].push_back((double) flag);
                        vFlagD[d].push_back((double) flag);
                        vFlagK[k].push_back((double) flag);
                        vFlagT[t].push_back((double) flag);
                        if (vFlagM.find(m) == vFlagM.end()) {
                            vFlagM[m] = vector <double> ();
                        }
                        vFlagM[m].push_back((double) flag);

                        if (flag == 1) {
                            sM.insert(m);

                            vMN[n].push_back((double) m);
                            vObjValN[n].push_back((double) objVal);
                            vGapN[n].push_back((double) gap);
                            vElapsedTimeN[n].push_back((double) elapsedTime);
                            vNSolutionN[n].push_back((double) nSolution);
                            vSumPenaltyN[n].push_back((double) sumPenalty);
                            vSumEdgeCostN[n].push_back((double) sumEdgeCost);
                            vObjValTimesElapsedTimeN[n].push_back((double) objValTimesElapsedTime);
                            vNSolutionDividedByNN[n].push_back((double) nSolutionDividedByN);
                            vSumPenaltyDividedByObjValN[n].push_back((double) sumPenaltyDividedByObjVal);
                            vSumEdgeCostDividedByObjValN[n].push_back((double) sumEdgeCostDividedByObjVal);

                            vMD[d].push_back((double) m);
                            vObjValD[d].push_back((double) objVal);
                            vGapD[d].push_back((double) gap);
                            vElapsedTimeD[d].push_back((double) elapsedTime);
                            vNSolutionD[d].push_back((double) nSolution);
                            vSumPenaltyD[d].push_back((double) sumPenalty);
                            vSumEdgeCostD[d].push_back((double) sumEdgeCost);
                            vObjValTimesElapsedTimeD[d].push_back((double) objValTimesElapsedTime);
                            vNSolutionDividedByND[d].push_back((double) nSolutionDividedByN);
                            vSumPenaltyDividedByObjValD[d].push_back((double) sumPenaltyDividedByObjVal);
                            vSumEdgeCostDividedByObjValD[d].push_back((double) sumEdgeCostDividedByObjVal);

                            vMK[k].push_back((double) m);
                            vObjValK[k].push_back((double) objVal);
                            vGapK[k].push_back((double) gap);
                            vElapsedTimeK[k].push_back((double) elapsedTime);
                            vNSolutionK[k].push_back((double) nSolution);
                            vSumPenaltyK[k].push_back((double) sumPenalty);
                            vSumEdgeCostK[k].push_back((double) sumEdgeCost);
                            vObjValTimesElapsedTimeK[k].push_back((double) objValTimesElapsedTime);
                            vNSolutionDividedByNK[k].push_back((double) nSolutionDividedByN);
                            vSumPenaltyDividedByObjValK[k].push_back((double) sumPenaltyDividedByObjVal);
                            vSumEdgeCostDividedByObjValK[k].push_back((double) sumEdgeCostDividedByObjVal);

                            vMT[t].push_back((double) m);
                            vObjValT[t].push_back((double) objVal);
                            vGapT[t].push_back((double) gap);
                            vElapsedTimeT[t].push_back((double) elapsedTime);
                            vNSolutionT[t].push_back((double) nSolution);
                            vSumPenaltyT[t].push_back((double) sumPenalty);
                            vSumEdgeCostT[t].push_back((double) sumEdgeCost);
                            vObjValTimesElapsedTimeT[t].push_back((double) objValTimesElapsedTime);
                            vNSolutionDividedByNT[t].push_back((double) nSolutionDividedByN);
                            vSumPenaltyDividedByObjValT[t].push_back((double) sumPenaltyDividedByObjVal);
                            vSumEdgeCostDividedByObjValT[t].push_back((double) sumEdgeCostDividedByObjVal);

                            if (vObjValM.find(m) == vObjValM.end()) {
                                vObjValM[m] = vector <double> ();
                            }
                            if (vGapM.find(m) == vGapM.end()) {
                                vGapM[m] = vector <double> ();
                            }
                            if (vElapsedTimeM.find(m) == vElapsedTimeM.end()) {
                                vElapsedTimeM[m] = vector <double> ();
                            }
                            if (vNSolutionM.find(m) == vNSolutionM.end()) {
                                vNSolutionM[m] = vector <double> ();
                            }
                            if (vSumPenaltyM.find(m) == vSumPenaltyM.end()) {
                                vSumPenaltyM[m] = vector <double> ();
                            }
                            if (vSumEdgeCostM.find(m) == vSumEdgeCostM.end()) {
                                vSumEdgeCostM[m] = vector <double> ();
                            }
                            if (vObjValTimesElapsedTimeM.find(m) == vObjValTimesElapsedTimeM.end()) {
                                vObjValTimesElapsedTimeM[m] = vector <double> ();
                            }
                            if (vNSolutionDividedByNM.find(m) == vNSolutionDividedByNM.end()) {
                                vNSolutionDividedByNM[m] = vector <double> ();
                            }
                            if (vSumPenaltyDividedByObjValM.find(m) == vSumPenaltyDividedByObjValM.end()) {
                                vSumPenaltyDividedByObjValM[m] = vector <double> ();
                            }
                            if (vSumEdgeCostDividedByObjValM.find(m) == vSumEdgeCostDividedByObjValM.end()) {
                                vSumEdgeCostDividedByObjValM[m] = vector <double> ();
                            }

                            vObjValM[m].push_back((double) objVal);
                            vGapM[m].push_back((double) gap);
                            vElapsedTimeM[m].push_back((double) elapsedTime);
                            vNSolutionM[m].push_back((double) nSolution);
                            vSumPenaltyM[m].push_back((double) sumPenalty);
                            vSumEdgeCostM[m].push_back((double) sumEdgeCost);
                            vObjValTimesElapsedTimeM[m].push_back((double) objValTimesElapsedTime);
                            vNSolutionDividedByNM[m].push_back((double) nSolutionDividedByN);
                            vSumPenaltyDividedByObjValM[m].push_back((double) sumPenaltyDividedByObjVal);
                            vSumEdgeCostDividedByObjValM[m].push_back((double) sumEdgeCostDividedByObjVal);
                        }

                        cout << n << ',' << d << ',' << k << ',' << t << ',' << i << ',' << m << ',' << flag << ',' << objVal << ',' << gap << ',' << elapsedTime << ',' << nSolution << ',' << sumPenalty << ',' << sumEdgeCost << ',' << objValTimesElapsedTime << ',' << nSolutionDividedByN << ',' << sumPenaltyDividedByObjVal << ',' << sumEdgeCostDividedByObjVal << endl;
                    }
                }
            }
        }
    }

    ofstream solutionNFile ("./output/solutionN.csv", ofstream :: out);
    solutionNFile << "n,m,deltaM,flag,deltaFlag,objVal,deltaObjVal,gap,deltaGap,elapsedTime,deltaElapsedTime,nSolution,deltaNSolution,sumPenalty,deltaSumPenalty,sumEdgeCost,deltaSumEdgeCost,objValTimesElapsedTime,deltaObjValTimesElapsedTime,nSolutionDividedByN,deltaNSolutionDividedByN,sumPenaltyDividedByObjVal,deltaSumPenaltyDividedByObjVal,sumEdgeCostDividedByObjVal,deltaSumEdgeCostDividedByObjVal" << endl;
    for (vector <ulint> :: iterator itN = vN.begin(); itN != vN.end(); itN++) {
        ulint n = *itN;
        solutionNFile << n << ',' << average(vMN[n]) << ',' << stdev(vMN[n]) << ',' << average(vFlagN[n]) << ',' << stdev(vFlagN[n]) << ',' << average(vObjValN[n]) << ',' << stdev(vObjValN[n]) << ',' << average(vGapN[n]) << ',' << stdev(vGapN[n]) << ',' << average(vElapsedTimeN[n]) << ',' << stdev(vElapsedTimeN[n]) << ',' << average(vNSolutionN[n]) << ',' << stdev(vNSolutionN[n]) << ',' << average(vSumPenaltyN[n]) << ',' << stdev(vSumPenaltyN[n]) << ',' << average(vSumEdgeCostN[n]) << ',' << stdev(vSumEdgeCostN[n]) << ',' << average(vObjValTimesElapsedTimeN[n]) << ',' << stdev(vObjValTimesElapsedTimeN[n]) << ',' << average(vNSolutionDividedByNN[n]) << ',' << stdev(vNSolutionDividedByNN[n]) << ',' << average(vSumPenaltyDividedByObjValN[n]) << ',' << stdev(vSumPenaltyDividedByObjValN[n]) << ',' << average(vSumEdgeCostDividedByObjValN[n]) << ',' << stdev(vSumEdgeCostDividedByObjValN[n]) << endl;
    }
    solutionNFile.close();

    ofstream solutionDFile ("./output/solutionD.csv", ofstream :: out);
    solutionDFile << "d,m,deltaM,flag,deltaFlag,objVal,deltaObjVal,gap,deltaGap,elapsedTime,deltaElapsedTime,nSolution,deltaNSolution,sumPenalty,deltaSumPenalty,sumEdgeCost,deltaSumEdgeCost,objValTimesElapsedTime,deltaObjValTimesElapsedTime,nSolutionDividedByN,deltaNSolutionDividedByN,sumPenaltyDividedByObjVal,deltaSumPenaltyDividedByObjVal,sumEdgeCostDividedByObjVal,deltaSumEdgeCostDividedByObjVal" << endl;
    for (vector <double> :: iterator itD = vD.begin(); itD != vD.end(); itD++) {
        double d = *itD;
        solutionDFile << d << ',' << average(vMD[d]) << ',' << stdev(vMD[d]) << ',' << average(vFlagD[d]) << ',' << stdev(vFlagD[d]) << ',' << average(vObjValD[d]) << ',' << stdev(vObjValD[d]) << ',' << average(vGapD[d]) << ',' << stdev(vGapD[d]) << ',' << average(vElapsedTimeD[d]) << ',' << stdev(vElapsedTimeD[d]) << ',' << average(vNSolutionD[d]) << ',' << stdev(vNSolutionD[d]) << ',' << average(vSumPenaltyD[d]) << ',' << stdev(vSumPenaltyD[d]) << ',' << average(vSumEdgeCostD[d]) << ',' << stdev(vSumEdgeCostD[d]) << ',' << average(vObjValTimesElapsedTimeD[d]) << ',' << stdev(vObjValTimesElapsedTimeD[d]) << ',' << average(vNSolutionDividedByND[d]) << ',' << stdev(vNSolutionDividedByND[d]) << ',' << average(vSumPenaltyDividedByObjValD[d]) << ',' << stdev(vSumPenaltyDividedByObjValD[d]) << ',' << average(vSumEdgeCostDividedByObjValD[d]) << ',' << stdev(vSumEdgeCostDividedByObjValD[d]) << endl;
    }
    solutionDFile.close();

    ofstream solutionKFile ("./output/solutionK.csv", ofstream :: out);
    solutionKFile << "k,m,deltaM,flag,deltaFlag,objVal,deltaObjVal,gap,deltaGap,elapsedTime,deltaElapsedTime,nSolution,deltaNSolution,sumPenalty,deltaSumPenalty,sumEdgeCost,deltaSumEdgeCost,objValTimesElapsedTime,deltaObjValTimesElapsedTime,nSolutionDividedByN,deltaNSolutionDividedByN,sumPenaltyDividedByObjVal,deltaSumPenaltyDividedByObjVal,sumEdgeCostDividedByObjVal,deltaSumEdgeCostDividedByObjVal" << endl;
    for (vector <ulint> :: iterator itK = vK.begin(); itK != vK.end(); itK++) {
        ulint k = *itK;
        solutionKFile << k << ',' << average(vMK[k]) << ',' << stdev(vMK[k]) << ',' << average(vFlagK[k]) << ',' << stdev(vFlagK[k]) << ',' << average(vObjValK[k]) << ',' << stdev(vObjValK[k]) << ',' << average(vGapK[k]) << ',' << stdev(vGapK[k]) << ',' << average(vElapsedTimeK[k]) << ',' << stdev(vElapsedTimeK[k]) << ',' << average(vNSolutionK[k]) << ',' << stdev(vNSolutionK[k]) << ',' << average(vSumPenaltyK[k]) << ',' << stdev(vSumPenaltyK[k]) << ',' << average(vSumEdgeCostK[k]) << ',' << stdev(vSumEdgeCostK[k]) << ',' << average(vObjValTimesElapsedTimeK[k]) << ',' << stdev(vObjValTimesElapsedTimeK[k]) << ',' << average(vNSolutionDividedByNK[k]) << ',' << stdev(vNSolutionDividedByNK[k]) << ',' << average(vSumPenaltyDividedByObjValK[k]) << ',' << stdev(vSumPenaltyDividedByObjValK[k]) << ',' << average(vSumEdgeCostDividedByObjValK[k]) << ',' << stdev(vSumEdgeCostDividedByObjValK[k]) << endl;
    }
    solutionKFile.close();

    ofstream solutionTFile ("./output/solutionT.csv", ofstream :: out);
    solutionTFile << "t,m,deltaM,flag,deltaFlag,objVal,deltaObjVal,gap,deltaGap,elapsedTime,deltaElapsedTime,nSolution,deltaNSolution,sumPenalty,deltaSumPenalty,sumEdgeCost,deltaSumEdgeCost,objValTimesElapsedTime,deltaObjValTimesElapsedTime,nSolutionDividedByN,deltaNSolutionDividedByN,sumPenaltyDividedByObjVal,deltaSumPenaltyDividedByObjVal,sumEdgeCostDividedByObjVal,deltaSumEdgeCostDividedByObjVal" << endl;
    for (vector <ulint> :: iterator itT = vT.begin(); itT != vT.end(); itT++) {
        ulint t = *itT;
        solutionTFile << t << ',' << average(vMT[t]) << ',' << stdev(vMT[t]) << ',' << average(vFlagT[t]) << ',' << stdev(vFlagT[t]) << ',' << average(vObjValT[t]) << ',' << stdev(vObjValT[t]) << ',' << average(vGapT[t]) << ',' << stdev(vGapT[t]) << ',' << average(vElapsedTimeT[t]) << ',' << stdev(vElapsedTimeT[t]) << ',' << average(vNSolutionT[t]) << ',' << stdev(vNSolutionT[t]) << ',' << average(vSumPenaltyT[t]) << ',' << stdev(vSumPenaltyT[t]) << ',' << average(vSumEdgeCostT[t]) << ',' << stdev(vSumEdgeCostT[t]) << ',' << average(vObjValTimesElapsedTimeT[t]) << ',' << stdev(vObjValTimesElapsedTimeT[t]) << ',' << average(vNSolutionDividedByNT[t]) << ',' << stdev(vNSolutionDividedByNT[t]) << ',' << average(vSumPenaltyDividedByObjValT[t]) << ',' << stdev(vSumPenaltyDividedByObjValT[t]) << ',' << average(vSumEdgeCostDividedByObjValT[t]) << ',' << stdev(vSumEdgeCostDividedByObjValT[t]) << endl;
    }
    solutionTFile.close();

    ofstream solutionMFile ("./output/solutionM.csv", ofstream :: out);
    solutionMFile << "m,flag,deltaFlag,objVal,deltaObjVal,gap,deltaGap,elapsedTime,deltaElapsedTime,nSolution,deltaNSolution,sumPenalty,deltaSumPenalty,sumEdgeCost,deltaSumEdgeCost,objValTimesElapsedTime,deltaObjValTimesElapsedTime,nSolutionDividedByN,deltaNSolutionDividedByN,sumPenaltyDividedByObjVal,deltaSumPenaltyDividedByObjVal,sumEdgeCostDividedByObjVal,deltaSumEdgeCostDividedByObjVal" << endl;
    for (set <ulint> :: iterator itM = sM.begin(); itM != sM.end(); itM++) {
        ulint m = *itM;
        solutionMFile << m << ',' << average(vFlagM[m]) << ',' << stdev(vFlagM[m]) << ',' << average(vObjValM[m]) << ',' << stdev(vObjValM[m]) << ',' << average(vGapM[m]) << ',' << stdev(vGapM[m]) << ',' << average(vElapsedTimeM[m]) << ',' << stdev(vElapsedTimeM[m]) << ',' << average(vNSolutionM[m]) << ',' << stdev(vNSolutionM[m]) << ',' << average(vSumPenaltyM[m]) << ',' << stdev(vSumPenaltyM[m]) << ',' << average(vSumEdgeCostM[m]) << ',' << stdev(vSumEdgeCostM[m]) << ',' << average(vObjValTimesElapsedTimeM[m]) << ',' << stdev(vObjValTimesElapsedTimeM[m]) << ',' << average(vNSolutionDividedByNM[m]) << ',' << stdev(vNSolutionDividedByNM[m]) << ',' << average(vSumPenaltyDividedByObjValM[m]) << ',' << stdev(vSumPenaltyDividedByObjValM[m]) << ',' << average(vSumEdgeCostDividedByObjValM[m]) << ',' << stdev(vSumEdgeCostDividedByObjValM[m]) << endl;
    }
    solutionMFile.close();

    return 0;
}
