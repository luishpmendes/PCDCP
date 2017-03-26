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
    map <ulint, vector <double> > vObjVal1N;
    map <ulint, vector <double> > vGap1N;
    map <ulint, vector <double> > vElapsedTime1N;
    map <ulint, vector <double> > vObjVal2N;
    map <ulint, vector <double> > vGap2N;
    map <ulint, vector <double> > vElapsedTime2N;
    map <ulint, vector <double> > vTotalElapsedTimeN;
    map <ulint, vector <double> > vNSolutionN;
    map <ulint, vector <double> > vSumPenaltyN;
    map <ulint, vector <double> > vSumEdgeCostN;
    map <ulint, vector <double> > vObjVal2TimesTotalElapsedTimeN;
    map <ulint, vector <double> > vNSolutionDividedByNN;
    map <ulint, vector <double> > vSumPenaltyDividedByObjVal2N;
    map <ulint, vector <double> > vSumEdgeCostDividedByObjVal2N;
    for (vector <ulint> :: iterator itN = vN.begin(); itN != vN.end(); itN++) {
        ulint n = *itN;
        vMN[n] = vector <double> ();
        vFlagN[n] = vector <double> ();
        vObjVal1N[n] = vector <double> ();
        vGap1N[n] = vector <double> ();
        vElapsedTime1N[n] = vector <double> ();
        vObjVal2N[n] = vector <double> ();
        vGap2N[n] = vector <double> ();
        vElapsedTime2N[n] = vector <double> ();
        vTotalElapsedTimeN[n] = vector <double> ();
        vNSolutionN[n] = vector <double> ();
        vSumPenaltyN[n] = vector <double> ();
        vSumEdgeCostN[n] = vector <double> ();
        vObjVal2TimesTotalElapsedTimeN[n] = vector <double> ();
        vNSolutionDividedByNN[n] = vector <double> ();
        vSumPenaltyDividedByObjVal2N[n] = vector <double> ();
        vSumEdgeCostDividedByObjVal2N[n] = vector <double> ();
    }

    map <double, vector <double> > vMD;
    map <double, vector <double> > vFlagD;
    map <double, vector <double> > vObjVal1D;
    map <double, vector <double> > vGap1D;
    map <double, vector <double> > vElapsedTime1D;
    map <double, vector <double> > vObjVal2D;
    map <double, vector <double> > vGap2D;
    map <double, vector <double> > vElapsedTime2D;
    map <double, vector <double> > vTotalElapsedTimeD;
    map <double, vector <double> > vNSolutionD;
    map <double, vector <double> > vSumPenaltyD;
    map <double, vector <double> > vSumEdgeCostD;
    map <double, vector <double> > vObjVal2TimesTotalElapsedTimeD;
    map <double, vector <double> > vNSolutionDividedByND;
    map <double, vector <double> > vSumPenaltyDividedByObjVal2D;
    map <double, vector <double> > vSumEdgeCostDividedByObjVal2D;
    for (vector <double> :: iterator itD = vD.begin(); itD != vD.end(); itD++) {
        double d = *itD;
        vMD[d] = vector <double> ();
        vFlagD[d] = vector <double> ();
        vObjVal1D[d] = vector <double> ();
        vGap1D[d] = vector <double> ();
        vElapsedTime1D[d] = vector <double> ();
        vObjVal2D[d] = vector <double> ();
        vGap2D[d] = vector <double> ();
        vElapsedTime2D[d] = vector <double> ();
        vTotalElapsedTimeD[d] = vector <double> ();
        vNSolutionD[d] = vector <double> ();
        vSumPenaltyD[d] = vector <double> ();
        vSumEdgeCostD[d] = vector <double> ();
        vObjVal2TimesTotalElapsedTimeD[d] = vector <double> ();
        vNSolutionDividedByND[d] = vector <double> ();
        vSumPenaltyDividedByObjVal2D[d] = vector <double> ();
        vSumEdgeCostDividedByObjVal2D[d] = vector <double> ();
    }

    map <ulint, vector <double> > vMK;
    map <ulint, vector <double> > vFlagK;
    map <ulint, vector <double> > vObjVal1K;
    map <ulint, vector <double> > vGap1K;
    map <ulint, vector <double> > vElapsedTime1K;
    map <ulint, vector <double> > vObjVal2K;
    map <ulint, vector <double> > vGap2K;
    map <ulint, vector <double> > vElapsedTime2K;
    map <ulint, vector <double> > vTotalElapsedTimeK;
    map <ulint, vector <double> > vNSolutionK;
    map <ulint, vector <double> > vSumPenaltyK;
    map <ulint, vector <double> > vSumEdgeCostK;
    map <ulint, vector <double> > vObjVal2TimesTotalElapsedTimeK;
    map <ulint, vector <double> > vNSolutionDividedByNK;
    map <ulint, vector <double> > vSumPenaltyDividedByObjVal2K;
    map <ulint, vector <double> > vSumEdgeCostDividedByObjVal2K;
    for (vector <ulint> :: iterator itK = vK.begin(); itK != vK.end(); itK++) {
        ulint k = *itK;
        vMK[k] = vector <double> ();
        vFlagK[k] = vector <double> ();
        vObjVal1K[k] = vector <double> ();
        vGap1K[k] = vector <double> ();
        vElapsedTime1K[k] = vector <double> ();
        vObjVal2K[k] = vector <double> ();
        vGap2K[k] = vector <double> ();
        vElapsedTime2K[k] = vector <double> ();
        vTotalElapsedTimeK[k] = vector <double> ();
        vNSolutionK[k] = vector <double> ();
        vSumPenaltyK[k] = vector <double> ();
        vSumEdgeCostK[k] = vector <double> ();
        vObjVal2TimesTotalElapsedTimeK[k] = vector <double> ();
        vNSolutionDividedByNK[k] = vector <double> ();
        vSumPenaltyDividedByObjVal2K[k] = vector <double> ();
        vSumEdgeCostDividedByObjVal2K[k] = vector <double> ();
    }

    map <ulint, vector <double> > vMT;
    map <ulint, vector <double> > vFlagT;
    map <ulint, vector <double> > vObjVal1T;
    map <ulint, vector <double> > vGap1T;
    map <ulint, vector <double> > vElapsedTime1T;
    map <ulint, vector <double> > vObjVal2T;
    map <ulint, vector <double> > vGap2T;
    map <ulint, vector <double> > vElapsedTime2T;
    map <ulint, vector <double> > vTotalElapsedTimeT;
    map <ulint, vector <double> > vNSolutionT;
    map <ulint, vector <double> > vSumPenaltyT;
    map <ulint, vector <double> > vSumEdgeCostT;
    map <ulint, vector <double> > vObjVal2TimesTotalElapsedTimeT;
    map <ulint, vector <double> > vNSolutionDividedByNT;
    map <ulint, vector <double> > vSumPenaltyDividedByObjVal2T;
    map <ulint, vector <double> > vSumEdgeCostDividedByObjVal2T;
    for (vector <ulint> :: iterator itT = vT.begin(); itT != vT.end(); itT++) {
        ulint t = *itT;
        vMT[t] = vector <double> ();
        vFlagT[t] = vector <double> ();
        vObjVal1T[t] = vector <double> ();
        vGap1T[t] = vector <double> ();
        vElapsedTime1T[t] = vector <double> ();
        vObjVal2T[t] = vector <double> ();
        vGap2T[t] = vector <double> ();
        vElapsedTime2T[t] = vector <double> ();
        vTotalElapsedTimeT[t] = vector <double> ();
        vNSolutionT[t] = vector <double> ();
        vSumPenaltyT[t] = vector <double> ();
        vSumEdgeCostT[t] = vector <double> ();
        vObjVal2TimesTotalElapsedTimeT[t] = vector <double> ();
        vNSolutionDividedByNT[t] = vector <double> ();
        vSumPenaltyDividedByObjVal2T[t] = vector <double> ();
        vSumEdgeCostDividedByObjVal2T[t] = vector <double> ();
    }

    map <ulint, vector <double> > vFlagM;
    map <ulint, vector <double> > vObjVal1M;
    map <ulint, vector <double> > vGap1M;
    map <ulint, vector <double> > vElapsedTime1M;
    map <ulint, vector <double> > vObjVal2M;
    map <ulint, vector <double> > vGap2M;
    map <ulint, vector <double> > vElapsedTime2M;
    map <ulint, vector <double> > vTotalElapsedTimeM;
    map <ulint, vector <double> > vNSolutionM;
    map <ulint, vector <double> > vSumPenaltyM;
    map <ulint, vector <double> > vSumEdgeCostM;
    map <ulint, vector <double> > vObjVal2TimesTotalElapsedTimeM;
    map <ulint, vector <double> > vNSolutionDividedByNM;
    map <ulint, vector <double> > vSumPenaltyDividedByObjVal2M;
    map <ulint, vector <double> > vSumEdgeCostDividedByObjVal2M;

    cout << "n,d,k,t,i,m,flag,objVal1,gap1,elapsedTime1,objVal2,gap2,elapsedTime2,totalElapsedTime,nSolution,sumPenalty,sumEdgeCost,objVal2*totalElapsedTime,nSolution/n,sumPenalty/objVal2,sumEdgeCost/objVal2" << endl;

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

                        double objVal1 = 0.0;
                        ifstream objVal1File ("./output/N" + N + "D" + D + "K" + K + "T" + T + "I" + I + "/objVal1.txt");
                        if (objVal1File.is_open()) {
                            objVal1File >> objVal1;
                        }

                        double gap1 = 0.0;
                        ifstream gap1File ("./output/N" + N + "D" + D + "K" + K + "T" + T + "I" + I + "/gap1.txt");
                        if (gap1File.is_open()) {
                            gap1File >> gap1;
                        }

                        ulint elapsedTime1 = 0.0;
                        ifstream elapsedTime1File ("./output/N" + N + "D" + D + "K" + K + "T" + T + "I" + I + "/elapsedTime1.txt");
                        if (elapsedTime1File.is_open()) {
                            elapsedTime1File >> elapsedTime1;
                        }

                        double objVal2 = 0.0;
                        ifstream objVal2File ("./output/N" + N + "D" + D + "K" + K + "T" + T + "I" + I + "/objVal2.txt");
                        if (objVal2File.is_open()) {
                            objVal2File >> objVal2;
                        }

                        double gap2 = 0.0;
                        ifstream gap2File ("./output/N" + N + "D" + D + "K" + K + "T" + T + "I" + I + "/gap2.txt");
                        if (gap2File.is_open()) {
                            gap2File >> gap2;
                        }

                        ulint elapsedTime2 = 0.0;
                        ifstream elapsedTime2File ("./output/N" + N + "D" + D + "K" + K + "T" + T + "I" + I + "/elapsedTime2.txt");
                        if (elapsedTime2File.is_open()) {
                            elapsedTime2File >> elapsedTime2;
                        }

                        ulint totalElapsedTime = elapsedTime1 + elapsedTime2;

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

                        double objVal2TimesTotalElapsedTime = (double) (objVal2 * ((double) totalElapsedTime));
                        double nSolutionDividedByN = (double) (((double) nSolution)/ ((double) n));
                        double sumPenaltyDividedByObjVal2 = sumPenalty / objVal2;
                        double sumEdgeCostDividedByObjVal2 = sumEdgeCost / objVal2;

                        ulint flag = 0;

                        if (objVal2 < 1E100) {
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
                            vObjVal1N[n].push_back((double) objVal1);
                            vGap1N[n].push_back((double) gap1);
                            vElapsedTime1N[n].push_back((double) elapsedTime1);
                            vObjVal2N[n].push_back((double) objVal2);
                            vGap2N[n].push_back((double) gap2);
                            vElapsedTime2N[n].push_back((double) elapsedTime2);
                            vTotalElapsedTimeN[n].push_back((double) totalElapsedTime);
                            vNSolutionN[n].push_back((double) nSolution);
                            vSumPenaltyN[n].push_back((double) sumPenalty);
                            vSumEdgeCostN[n].push_back((double) sumEdgeCost);
                            vObjVal2TimesTotalElapsedTimeN[n].push_back((double) objVal2TimesTotalElapsedTime);
                            vNSolutionDividedByNN[n].push_back((double) nSolutionDividedByN);
                            vSumPenaltyDividedByObjVal2N[n].push_back((double) sumPenaltyDividedByObjVal2);
                            vSumEdgeCostDividedByObjVal2N[n].push_back((double) sumEdgeCostDividedByObjVal2);

                            vMD[d].push_back((double) m);
                            vObjVal1D[d].push_back((double) objVal1);
                            vGap1D[d].push_back((double) gap1);
                            vElapsedTime1D[d].push_back((double) elapsedTime1);
                            vObjVal2D[d].push_back((double) objVal2);
                            vGap2D[d].push_back((double) gap2);
                            vElapsedTime2D[d].push_back((double) elapsedTime2);
                            vTotalElapsedTimeD[d].push_back((double) totalElapsedTime);
                            vNSolutionD[d].push_back((double) nSolution);
                            vSumPenaltyD[d].push_back((double) sumPenalty);
                            vSumEdgeCostD[d].push_back((double) sumEdgeCost);
                            vObjVal2TimesTotalElapsedTimeD[d].push_back((double) objVal2TimesTotalElapsedTime);
                            vNSolutionDividedByND[d].push_back((double) nSolutionDividedByN);
                            vSumPenaltyDividedByObjVal2D[d].push_back((double) sumPenaltyDividedByObjVal2);
                            vSumEdgeCostDividedByObjVal2D[d].push_back((double) sumEdgeCostDividedByObjVal2);

                            vMK[k].push_back((double) m);
                            vObjVal1K[k].push_back((double) objVal1);
                            vGap1K[k].push_back((double) gap1);
                            vElapsedTime1K[k].push_back((double) elapsedTime1);
                            vObjVal2K[k].push_back((double) objVal2);
                            vGap2K[k].push_back((double) gap2);
                            vElapsedTime2K[k].push_back((double) elapsedTime2);
                            vTotalElapsedTimeK[k].push_back((double) totalElapsedTime);
                            vNSolutionK[k].push_back((double) nSolution);
                            vSumPenaltyK[k].push_back((double) sumPenalty);
                            vSumEdgeCostK[k].push_back((double) sumEdgeCost);
                            vObjVal2TimesTotalElapsedTimeK[k].push_back((double) objVal2TimesTotalElapsedTime);
                            vNSolutionDividedByNK[k].push_back((double) nSolutionDividedByN);
                            vSumPenaltyDividedByObjVal2K[k].push_back((double) sumPenaltyDividedByObjVal2);
                            vSumEdgeCostDividedByObjVal2K[k].push_back((double) sumEdgeCostDividedByObjVal2);

                            vMT[t].push_back((double) m);
                            vObjVal1T[t].push_back((double) objVal1);
                            vGap1T[t].push_back((double) gap1);
                            vElapsedTime1T[t].push_back((double) elapsedTime1);
                            vObjVal2T[t].push_back((double) objVal2);
                            vGap2T[t].push_back((double) gap2);
                            vElapsedTime2T[t].push_back((double) elapsedTime2);
                            vTotalElapsedTimeT[t].push_back((double) totalElapsedTime);
                            vNSolutionT[t].push_back((double) nSolution);
                            vSumPenaltyT[t].push_back((double) sumPenalty);
                            vSumEdgeCostT[t].push_back((double) sumEdgeCost);
                            vObjVal2TimesTotalElapsedTimeT[t].push_back((double) objVal2TimesTotalElapsedTime);
                            vNSolutionDividedByNT[t].push_back((double) nSolutionDividedByN);
                            vSumPenaltyDividedByObjVal2T[t].push_back((double) sumPenaltyDividedByObjVal2);
                            vSumEdgeCostDividedByObjVal2T[t].push_back((double) sumEdgeCostDividedByObjVal2);

                            if (vObjVal1M.find(m) == vObjVal1M.end()) {
                                vObjVal1M[m] = vector <double> ();
                            }
                            if (vGap1M.find(m) == vGap1M.end()) {
                                vGap1M[m] = vector <double> ();
                            }
                            if (vElapsedTime1M.find(m) == vElapsedTime1M.end()) {
                                vElapsedTime1M[m] = vector <double> ();
                            }
                            if (vObjVal2M.find(m) == vObjVal2M.end()) {
                                vObjVal2M[m] = vector <double> ();
                            }
                            if (vGap2M.find(m) == vGap2M.end()) {
                                vGap2M[m] = vector <double> ();
                            }
                            if (vElapsedTime2M.find(m) == vElapsedTime2M.end()) {
                                vElapsedTime2M[m] = vector <double> ();
                            }
                            if (vTotalElapsedTimeM.find(m) == vTotalElapsedTimeM.end()) {
                                vTotalElapsedTimeM[m] = vector <double> ();
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
                            if (vObjVal2TimesTotalElapsedTimeM.find(m) == vObjVal2TimesTotalElapsedTimeM.end()) {
                                vObjVal2TimesTotalElapsedTimeM[m] = vector <double> ();
                            }
                            if (vNSolutionDividedByNM.find(m) == vNSolutionDividedByNM.end()) {
                                vNSolutionDividedByNM[m] = vector <double> ();
                            }
                            if (vSumPenaltyDividedByObjVal2M.find(m) == vSumPenaltyDividedByObjVal2M.end()) {
                                vSumPenaltyDividedByObjVal2M[m] = vector <double> ();
                            }
                            if (vSumEdgeCostDividedByObjVal2M.find(m) == vSumEdgeCostDividedByObjVal2M.end()) {
                                vSumEdgeCostDividedByObjVal2M[m] = vector <double> ();
                            }

                            vObjVal1M[m].push_back((double) objVal1);
                            vGap1M[m].push_back((double) gap1);
                            vElapsedTime1M[m].push_back((double) elapsedTime1);
                            vObjVal2M[m].push_back((double) objVal2);
                            vGap2M[m].push_back((double) gap2);
                            vElapsedTime2M[m].push_back((double) elapsedTime2);
                            vTotalElapsedTimeM[m].push_back((double) totalElapsedTime);
                            vNSolutionM[m].push_back((double) nSolution);
                            vSumPenaltyM[m].push_back((double) sumPenalty);
                            vSumEdgeCostM[m].push_back((double) sumEdgeCost);
                            vObjVal2TimesTotalElapsedTimeM[m].push_back((double) objVal2TimesTotalElapsedTime);
                            vNSolutionDividedByNM[m].push_back((double) nSolutionDividedByN);
                            vSumPenaltyDividedByObjVal2M[m].push_back((double) sumPenaltyDividedByObjVal2);
                            vSumEdgeCostDividedByObjVal2M[m].push_back((double) sumEdgeCostDividedByObjVal2);
                        }

                        cout << n << ',' << d << ',' << k << ',' << t << ',' << i << ',' << m << ',' << flag << ',' << objVal1 << ',' << gap1 << ',' << elapsedTime1 << ',' << objVal2 << ',' << gap2 << ',' << elapsedTime2 << ',' << totalElapsedTime << ',' << nSolution << ',' << sumPenalty << ',' << sumEdgeCost << ',' << objVal2TimesTotalElapsedTime << ',' << nSolutionDividedByN << ',' << sumPenaltyDividedByObjVal2 << ',' << sumEdgeCostDividedByObjVal2 << endl;
                    }
                }
            }
        }
    }

    ofstream solutionNFile ("./output/solutionN.csv", ofstream :: out);
    solutionNFile << "n,m,deltaM,flag,deltaFlag,objVal1,deltaObjVal1,gap1,deltaGap1,elapsedTime1,deltaElapsedTime1,objVal2,deltaObjVal2,gap2,deltaGap2,elapsedTime2,deltaElapsedTime2,totalElapsedTime,deltaTotalElapsedTime,nSolution,deltaNSolution,sumPenalty,deltaSumPenalty,sumEdgeCost,deltaSumEdgeCost,objVal2TimesTotalElapsedTime,deltaObjVal2TimesTotalElapsedTime,nSolutionDividedByN,deltaNSolutionDividedByN,sumPenaltyDividedByObjVal2,deltaSumPenaltyDividedByObjVal2,sumEdgeCostDividedByObjVal2,deltaSumEdgeCostDividedByObjVal2" << endl;
    for (vector <ulint> :: iterator itN = vN.begin(); itN != vN.end(); itN++) {
        ulint n = *itN;
        solutionNFile << n << ',' << average(vMN[n]) << ',' << stdev(vMN[n]) << ',' << average(vFlagN[n]) << ',' << stdev(vFlagN[n]) << ',' << average(vObjVal1N[n]) << ',' << stdev(vObjVal1N[n]) << ',' << average(vGap1N[n]) << ',' << stdev(vGap1N[n]) << ',' << average(vElapsedTime1N[n]) << ',' << stdev(vElapsedTime1N[n]) << ',' << average(vObjVal2N[n]) << ',' << stdev(vObjVal2N[n]) << ',' << average(vGap2N[n]) << ',' << stdev(vGap2N[n]) << ',' << average(vElapsedTime2N[n]) << ',' << stdev(vElapsedTime2N[n]) << ',' << average(vTotalElapsedTimeN[n]) << ',' << stdev(vTotalElapsedTimeN[n]) << ',' << average(vNSolutionN[n]) << ',' << stdev(vNSolutionN[n]) << ',' << average(vSumPenaltyN[n]) << ',' << stdev(vSumPenaltyN[n]) << ',' << average(vSumEdgeCostN[n]) << ',' << stdev(vSumEdgeCostN[n]) << ',' << average(vObjVal2TimesTotalElapsedTimeN[n]) << ',' << stdev(vObjVal2TimesTotalElapsedTimeN[n]) << ',' << average(vNSolutionDividedByNN[n]) << ',' << stdev(vNSolutionDividedByNN[n]) << ',' << average(vSumPenaltyDividedByObjVal2N[n]) << ',' << stdev(vSumPenaltyDividedByObjVal2N[n]) << ',' << average(vSumEdgeCostDividedByObjVal2N[n]) << ',' << stdev(vSumEdgeCostDividedByObjVal2N[n]) << endl;
    }
    solutionNFile.close();

    ofstream solutionDFile ("./output/solutionD.csv", ofstream :: out);
    solutionDFile << "d,m,deltaM,flag,deltaFlag,objVal1,deltaObjVal1,gap1,deltaGap1,elapsedTime1,deltaElapsedTime1,objVal2,deltaObjVal2,gap2,deltaGap2,elapsedTime2,deltaElapsedTime2,totalElapsedTime,deltaTotalElapsedTime,nSolution,deltaNSolution,sumPenalty,deltaSumPenalty,sumEdgeCost,deltaSumEdgeCost,objVal2TimesTotalElapsedTime,deltaObjVal2TimesTotalElapsedTime,nSolutionDividedByN,deltaNSolutionDividedByN,sumPenaltyDividedByObjVal2,deltaSumPenaltyDividedByObjVal2,sumEdgeCostDividedByObjVal2,deltaSumEdgeCostDividedByObjVal2" << endl;
    for (vector <double> :: iterator itD = vD.begin(); itD != vD.end(); itD++) {
        double d = *itD;
        solutionDFile << d << ',' << average(vMD[d]) << ',' << stdev(vMD[d]) << ',' << average(vFlagD[d]) << ',' << stdev(vFlagD[d]) << ',' << average(vObjVal1D[d]) << ',' << stdev(vObjVal1D[d]) << ',' << average(vGap1D[d]) << ',' << stdev(vGap1D[d]) << ',' << average(vElapsedTime1D[d]) << ',' << stdev(vElapsedTime1D[d]) << ',' << average(vObjVal2D[d]) << ',' << stdev(vObjVal2D[d]) << ',' << average(vGap2D[d]) << ',' << stdev(vGap2D[d]) << ',' << average(vElapsedTime2D[d]) << ',' << stdev(vElapsedTime2D[d]) << ',' << average(vTotalElapsedTimeD[d]) << ',' << stdev(vTotalElapsedTimeD[d]) << ',' << average(vNSolutionD[d]) << ',' << stdev(vNSolutionD[d]) << ',' << average(vSumPenaltyD[d]) << ',' << stdev(vSumPenaltyD[d]) << ',' << average(vSumEdgeCostD[d]) << ',' << stdev(vSumEdgeCostD[d]) << ',' << average(vObjVal2TimesTotalElapsedTimeD[d]) << ',' << stdev(vObjVal2TimesTotalElapsedTimeD[d]) << ',' << average(vNSolutionDividedByND[d]) << ',' << stdev(vNSolutionDividedByND[d]) << ',' << average(vSumPenaltyDividedByObjVal2D[d]) << ',' << stdev(vSumPenaltyDividedByObjVal2D[d]) << ',' << average(vSumEdgeCostDividedByObjVal2D[d]) << ',' << stdev(vSumEdgeCostDividedByObjVal2D[d]) << endl;
    }
    solutionDFile.close();

    ofstream solutionKFile ("./output/solutionK.csv", ofstream :: out);
    solutionKFile << "k,m,deltaM,flag,deltaFlag,objVal1,deltaObjVal1,gap1,deltaGap1,elapsedTime1,deltaElapsedTime1,objVal2,deltaObjVal2,gap2,deltaGap2,elapsedTime2,deltaElapsedTime2,totalElapsedTime,deltaTotalElapsedTime,nSolution,deltaNSolution,sumPenalty,deltaSumPenalty,sumEdgeCost,deltaSumEdgeCost,objVal2TimesTotalElapsedTime,deltaObjVal2TimesTotalElapsedTime,nSolutionDividedByN,deltaNSolutionDividedByN,sumPenaltyDividedByObjVal2,deltaSumPenaltyDividedByObjVal2,sumEdgeCostDividedByObjVal2,deltaSumEdgeCostDividedByObjVal2" << endl;
    for (vector <ulint> :: iterator itK = vK.begin(); itK != vK.end(); itK++) {
        ulint k = *itK;
        solutionKFile << k << ',' << average(vMK[k]) << ',' << stdev(vMK[k]) << ',' << average(vFlagK[k]) << ',' << stdev(vFlagK[k]) << ',' << average(vObjVal1K[k]) << ',' << stdev(vObjVal1K[k]) << ',' << average(vGap1K[k]) << ',' << stdev(vGap1K[k]) << ',' << average(vElapsedTime1K[k]) << ',' << stdev(vElapsedTime1K[k]) << ',' << average(vObjVal2K[k]) << ',' << stdev(vObjVal2K[k]) << ',' << average(vGap2K[k]) << ',' << stdev(vGap2K[k]) << ',' << average(vElapsedTime2K[k]) << ',' << stdev(vElapsedTime2K[k]) << ',' << average(vTotalElapsedTimeK[k]) << ',' << stdev(vTotalElapsedTimeK[k]) << ',' << average(vNSolutionK[k]) << ',' << stdev(vNSolutionK[k]) << ',' << average(vSumPenaltyK[k]) << ',' << stdev(vSumPenaltyK[k]) << ',' << average(vSumEdgeCostK[k]) << ',' << stdev(vSumEdgeCostK[k]) << ',' << average(vObjVal2TimesTotalElapsedTimeK[k]) << ',' << stdev(vObjVal2TimesTotalElapsedTimeK[k]) << ',' << average(vNSolutionDividedByNK[k]) << ',' << stdev(vNSolutionDividedByNK[k]) << ',' << average(vSumPenaltyDividedByObjVal2K[k]) << ',' << stdev(vSumPenaltyDividedByObjVal2K[k]) << ',' << average(vSumEdgeCostDividedByObjVal2K[k]) << ',' << stdev(vSumEdgeCostDividedByObjVal2K[k]) << endl;
    }
    solutionKFile.close();

    ofstream solutionTFile ("./output/solutionT.csv", ofstream :: out);
    solutionTFile << "t,m,deltaM,flag,deltaFlag,objVal1,deltaObjVal1,gap1,deltaGap1,elapsedTime1,deltaElapsedTime1,objVal2,deltaObjVal2,gap2,deltaGap2,elapsedTime2,deltaElapsedTime2,totalElapsedTime,deltaTotalElapsedTime,nSolution,deltaNSolution,sumPenalty,deltaSumPenalty,sumEdgeCost,deltaSumEdgeCost,objVal2TimesTotalElapsedTime,deltaObjVal2TimesTotalElapsedTime,nSolutionDividedByN,deltaNSolutionDividedByN,sumPenaltyDividedByObjVal2,deltaSumPenaltyDividedByObjVal2,sumEdgeCostDividedByObjVal2,deltaSumEdgeCostDividedByObjVal2" << endl;
    for (vector <ulint> :: iterator itT = vT.begin(); itT != vT.end(); itT++) {
        ulint t = *itT;
        solutionTFile << t << ',' << average(vMT[t]) << ',' << stdev(vMT[t]) << ',' << average(vFlagT[t]) << ',' << stdev(vFlagT[t]) << ',' << average(vObjVal1T[t]) << ',' << stdev(vObjVal1T[t]) << ',' << average(vGap1T[t]) << ',' << stdev(vGap1T[t]) << ',' << average(vElapsedTime1T[t]) << ',' << stdev(vElapsedTime1T[t]) << ',' << average(vObjVal2T[t]) << ',' << stdev(vObjVal2T[t]) << ',' << average(vGap2T[t]) << ',' << stdev(vGap2T[t]) << ',' << average(vElapsedTime2T[t]) << ',' << stdev(vElapsedTime2T[t]) << ',' << average(vTotalElapsedTimeT[t]) << ',' << stdev(vTotalElapsedTimeT[t]) << ',' << average(vNSolutionT[t]) << ',' << stdev(vNSolutionT[t]) << ',' << average(vSumPenaltyT[t]) << ',' << stdev(vSumPenaltyT[t]) << ',' << average(vSumEdgeCostT[t]) << ',' << stdev(vSumEdgeCostT[t]) << ',' << average(vObjVal2TimesTotalElapsedTimeT[t]) << ',' << stdev(vObjVal2TimesTotalElapsedTimeT[t]) << ',' << average(vNSolutionDividedByNT[t]) << ',' << stdev(vNSolutionDividedByNT[t]) << ',' << average(vSumPenaltyDividedByObjVal2T[t]) << ',' << stdev(vSumPenaltyDividedByObjVal2T[t]) << ',' << average(vSumEdgeCostDividedByObjVal2T[t]) << ',' << stdev(vSumEdgeCostDividedByObjVal2T[t]) << endl;
    }
    solutionTFile.close();

    ofstream solutionMFile ("./output/solutionM.csv", ofstream :: out);
    solutionMFile << "m,flag,deltaFlag,objVal1,deltaObjVal1,gap1,deltaGap1,elapsedTime1,deltaElapsedTime1,objVal2,deltaObjVal2,gap2,deltaGap2,elapsedTime2,deltaElapsedTime2,totalElapsedTime,deltaTotalElapsedTime,nSolution,deltaNSolution,sumPenalty,deltaSumPenalty,sumEdgeCost,deltaSumEdgeCost,objVal2TimesTotalElapsedTime,deltaObjVal2TimesTotalElapsedTime,nSolutionDividedByN,deltaNSolutionDividedByN,sumPenaltyDividedByObjVal2,deltaSumPenaltyDividedByObjVal2,sumEdgeCostDividedByObjVal2,deltaSumEdgeCostDividedByObjVal2" << endl;
    for (set <ulint> :: iterator itM = sM.begin(); itM != sM.end(); itM++) {
        ulint m = *itM;
        solutionMFile << m << ',' << average(vFlagM[m]) << ',' << stdev(vFlagM[m]) << ',' << average(vObjVal1M[m]) << ',' << stdev(vObjVal1M[m]) << ',' << average(vGap1M[m]) << ',' << stdev(vGap1M[m]) << ',' << average(vElapsedTime1M[m]) << ',' << stdev(vElapsedTime1M[m]) << ',' << average(vObjVal2M[m]) << ',' << stdev(vObjVal2M[m]) << ',' << average(vGap2M[m]) << ',' << stdev(vGap2M[m]) << ',' << average(vElapsedTime2M[m]) << ',' << stdev(vElapsedTime2M[m]) << ',' << average(vTotalElapsedTimeM[m]) << ',' << stdev(vTotalElapsedTimeM[m]) << ',' << average(vNSolutionM[m]) << ',' << stdev(vNSolutionM[m]) << ',' << average(vSumPenaltyM[m]) << ',' << stdev(vSumPenaltyM[m]) << ',' << average(vSumEdgeCostM[m]) << ',' << stdev(vSumEdgeCostM[m]) << ',' << average(vObjVal2TimesTotalElapsedTimeM[m]) << ',' << stdev(vObjVal2TimesTotalElapsedTimeM[m]) << ',' << average(vNSolutionDividedByNM[m]) << ',' << stdev(vNSolutionDividedByNM[m]) << ',' << average(vSumPenaltyDividedByObjVal2M[m]) << ',' << stdev(vSumPenaltyDividedByObjVal2M[m]) << ',' << average(vSumEdgeCostDividedByObjVal2M[m]) << ',' << stdev(vSumEdgeCostDividedByObjVal2M[m]) << endl;
    }
    solutionMFile.close();

    return 0;
}
