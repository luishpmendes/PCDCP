#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <map>

using namespace std;

typedef long int lint;
typedef unsigned long int ulint;
typedef vector < vector <lint> > matrix;

int main () {
    vector <int> vN;
    vN.push_back(10);
    vN.push_back(20);
    vN.push_back(50);
    vN.push_back(100);
    vector <double> vD;
    vD.push_back(0.3);
    vD.push_back(0.5);
    vD.push_back(0.7);
    vector <int> vK;
    vK.push_back(0);
    vK.push_back(10);
    vK.push_back(20);
    vector <int> vT;
    vT.push_back(0);
    vT.push_back(1);
    vector <double> vP;
    vP.push_back(0.2);
    vector <int> vI;
    vI.push_back(0);

    cout << "n,d,k,t,p,i,objVal,gap,elapsedTime,nSolution,sumPenalty,sumEdgeCost" << endl;

    for (vector <int> :: iterator itN = vN.begin(); itN != vN.end(); itN++) {
        int n = *itN;
        stringstream ssN;
        ssN << n;
        string N = ssN.str();
        for (vector <double> :: iterator itD = vD.begin(); itD != vD.end(); itD++) {
            double d = *itD;
            stringstream ssD;
            ssD << fixed << setprecision(1) << d;
            string D = ssD.str();
            D.erase(remove(D.begin(), D.end(), '.'), D.end());
            for (vector <int> :: iterator itK = vK.begin(); itK != vK.end(); itK++) {
                int k = *itK;
                stringstream ssK;
                ssK << k;
                string K = ssK.str();
                for (vector <int> :: iterator itT = vT.begin(); itT != vT.end(); itT++) {
                    int t = *itT;
                    stringstream ssT;
                    ssT << t;
                    string T = ssT.str();
                    for (vector <double> :: iterator itP = vP.begin(); itP != vP.end(); itP++) {
                        double p = *itP;
                        stringstream ssP;
                        ssP << fixed << setprecision(1) << p;
                        string P = ssP.str();
                        P.erase(remove(P.begin(), P.end(), '.'), P.end());
                        for (vector <int> :: iterator itI = vI.begin(); itI != vI.end(); itI++) {
                            int i = *itI;
                            stringstream ssI;
                            ssI << i;
                            string I = ssI.str();

                            double objVal = 0.0;
                            ifstream objValFile ("./output/N" + N + "D" + D + "K" + K + "T" + T + "P" + P + "I" + I + "/objVal.txt");
                            if (objValFile.is_open()) {
                                objValFile >> objVal;
                            }

                            double gap = 0.0;
                            ifstream gapFile ("./output/N" + N + "D" + D + "K" + K + "T" + T + "P" + P + "I" + I + "/gap.txt");
                            if (gapFile.is_open()) {
                                gapFile >> gap;
                            }

                            ulint elapsedTime = 0.0;
                            ifstream elapsedTimeFile ("./output/N" + N + "D" + D + "K" + K + "T" + T + "P" + P + "I" + I + "/elapsedTime.txt");
                            if (elapsedTimeFile.is_open()) {
                                elapsedTimeFile >> elapsedTime;
                            }

                            double nSolution = 0.0, sumPenalty = 0.0, sumEdgeCost = 0.0;

                            ifstream inputFile ("../input/instanceN" + N + "D" + D + "K" + K + "T" + T + "P" + P + "I" + I + ".in");
                            ifstream resultFile ("./output/N" + N + "D" + D + "K" + K + "T" + T + "P" + P + "I" + I + "/result.out");

                            if (inputFile.is_open() && resultFile.is_open()) {
                                ulint n, mComplete, m, k, t, root;
                                double d, p;
                                inputFile >> n >> d >> k >> t >> p >> mComplete >> m >> root;

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

                            cout << n << ',' << d << ',' << k << ',' << t << ',' << p << ',' << i << ',' << objVal << ',' << gap << ',' << elapsedTime << ',' << nSolution << ',' << sumPenalty << ',' << sumEdgeCost << endl;
                        }
                    }
                }
            }
        }
    }
    return 0;
}
