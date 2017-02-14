#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>

using namespace std;

int main (int argc, char * argv[]) {
    string path;

    if (argc >= 2) {
        path = string(argv[1]);
    } else {
        cin >> path;
    }

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

    if (path.compare("grasp") == 0) {
        cout << "n,d,k,t,p,i,a,nFlag,kFlag,nSolutionFlag,mSolutionFlag, solutionVerticesFlag,solutionEdgesFlag,solutionCostFlag,coverFlag,solutionVerticesNotInMainCycleFlag,nonSolutionVerticesInMainCycleFlag,errorFlag" << endl;
    } else if (path.compare("geneticAlgorithm") == 0) {
        cout << "n,d,k,t,p,i,ps,mr,nFlag,kFlag,nSolutionFlag,mSolutionFlag, solutionVerticesFlag,solutionEdgesFlag,solutionCostFlag,coverFlag,solutionVerticesNotInMainCycleFlag,nonSolutionVerticesInMainCycleFlag,errorFlag" << endl;
    } else {
        cout << "n,d,k,t,p,i,nFlag,kFlag,nSolutionFlag,mSolutionFlag, solutionVerticesFlag,solutionEdgesFlag,solutionCostFlag,coverFlag,solutionVerticesNotInMainCycleFlag,nonSolutionVerticesInMainCycleFlag,errorFlag" << endl;
    }
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

                            if (path.compare("grasp") == 0) {
                                vector <double> vA;
                                vA.push_back(0.3);
                                vA.push_back(0.5);
                                vA.push_back(0.7);
                                for (vector <double> :: iterator itA = vA.begin(); itA != vA.end(); itA++) {
                                    double a = *itA;
                                    stringstream ssA;
                                    ssA << fixed << setprecision(1) << a;
                                    string A = ssA.str();
                                    A.erase(remove(A.begin(), A.end(), '.'), A.end());

                                    int nFlag, kFlag, nSolutionFlag, mSolutionFlag, solutionVerticesFlag, solutionEdgesFlag, solutionCostFlag, coverFlag, solutionVerticesNotInMainCycleFlag, nonSolutionVerticesInMainCycleFlag, errorFlag;
                                    nFlag = kFlag = nSolutionFlag = mSolutionFlag = solutionVerticesFlag = solutionEdgesFlag = solutionCostFlag = coverFlag = solutionVerticesNotInMainCycleFlag = nonSolutionVerticesInMainCycleFlag = errorFlag = 1;

                                    ifstream nFlagFile = ifstream("../" + path + "/output/N" + N + "D" + D + "K" + K + "T" + T + "P" + P + "I" + I + "A" + A + "/nFlag.txt");
                                    if (nFlagFile.is_open()) {
                                        nFlagFile >> nFlag;
                                    }

                                    ifstream kFlagFile = ifstream("../" + path + "/output/N" + N + "D" + D + "K" + K + "T" + T + "P" + P + "I" + I + "A" + A + "/kFlag.txt");
                                    if (kFlagFile.is_open()) {
                                        kFlagFile >> kFlag;
                                    }

                                    ifstream nSolutionFlagFile = ifstream("../" + path + "/output/N" + N + "D" + D + "K" + K + "T" + T + "P" + P + "I" + I + "A" + A + "/nSolutionFlag.txt");
                                    if (nSolutionFlagFile.is_open()) {
                                        nSolutionFlagFile >> nSolutionFlag;
                                    }

                                    ifstream mSolutionFlagFile = ifstream("../" + path + "/output/N" + N + "D" + D + "K" + K + "T" + T + "P" + P + "I" + I + "A" + A + "/mSolutionFlag.txt");
                                    if (mSolutionFlagFile.is_open()) {
                                        mSolutionFlagFile >> mSolutionFlag;
                                    }

                                    ifstream solutionVerticesFlagFile = ifstream("../" + path + "/output/N" + N + "D" + D + "K" + K + "T" + T + "P" + P + "I" + I + "A" + A + "/solutionVerticesFlag.txt");
                                    if (solutionVerticesFlagFile.is_open()) {
                                        solutionVerticesFlagFile >> solutionVerticesFlag;
                                    }

                                    ifstream solutionEdgesFlagFile = ifstream("../" + path + "/output/N" + N + "D" + D + "K" + K + "T" + T + "P" + P + "I" + I + "A" + A + "/solutionEdgesFlag.txt");
                                    if (solutionEdgesFlagFile.is_open()) {
                                        solutionEdgesFlagFile >> solutionEdgesFlag;
                                    }

                                    ifstream solutionCostFlagFile = ifstream("../" + path + "/output/N" + N + "D" + D + "K" + K + "T" + T + "P" + P + "I" + I + "A" + A + "/solutionCostFlag.txt");
                                    if (solutionCostFlagFile.is_open()) {
                                        solutionCostFlagFile >> solutionCostFlag;
                                    }

                                    ifstream coverFlagFile = ifstream("../" + path + "/output/N" + N + "D" + D + "K" + K + "T" + T + "P" + P + "I" + I + "A" + A + "/coverFlag.txt");
                                    if (coverFlagFile.is_open()) {
                                        coverFlagFile >> coverFlag;
                                    }

                                    ifstream solutionVerticesNotInMainCycleFlagFile = ifstream("../" + path + "/output/N" + N + "D" + D + "K" + K + "T" + T + "P" + P + "I" + I + "A" + A + "/solutionVerticesNotInMainCycleFlag.txt");
                                    if (solutionVerticesNotInMainCycleFlagFile.is_open()) {
                                        solutionVerticesNotInMainCycleFlagFile >> solutionVerticesNotInMainCycleFlag;
                                    }

                                    ifstream nonSolutionVerticesInMainCycleFlagFile = ifstream("../" + path + "/output/N" + N + "D" + D + "K" + K + "T" + T + "P" + P + "I" + I + "A" + A + "/nonSolutionVerticesInMainCycleFlag.txt");
                                    if (nonSolutionVerticesInMainCycleFlagFile.is_open()) {
                                        nonSolutionVerticesInMainCycleFlagFile >> nonSolutionVerticesInMainCycleFlag;
                                    }

                                    ifstream errorFlagFile = ifstream("../" + path + "/output/N" + N + "D" + D + "K" + K + "T" + T + "P" + P + "I" + I + "A" + A + "/errorFlag.txt");
                                    if (errorFlagFile.is_open()) {
                                        errorFlagFile >> errorFlag;
                                    }

                                    cout << n << ',' << d << ',' << k << ',' << t << ',' << p << ',' << i << ',' << a << ',' << nFlag << ',' << kFlag << ',' << nSolutionFlag << ',' << mSolutionFlag << ',' << solutionVerticesFlag << ',' << solutionEdgesFlag << ',' << solutionCostFlag << ',' << coverFlag << ',' << solutionVerticesNotInMainCycleFlag << ',' << nonSolutionVerticesInMainCycleFlag << ',' << errorFlag << endl;
                                }
                            } else if (path.compare("geneticAlgorithm") == 0) {
                                vector <int> vPS;
                                vPS.push_back(10);
                                vPS.push_back(50);
                                vPS.push_back(100);
                                for (vector <int> :: iterator itPS = vPS.begin(); itPS != vPS.end(); itPS++) {
                                    int ps = *itPS;
                                    stringstream ssPS;
                                    ssPS << ps;
                                    string PS = ssPS.str();

                                    vector <double> vMR;
                                    vMR.push_back(0.3);
                                    vMR.push_back(0.5);
                                    vMR.push_back(0.7);
                                    for (vector <double> :: iterator itMR = vMR.begin(); itMR != vMR.end(); itMR++) {
                                        double mr = *itMR;
                                        stringstream ssMR;
                                        ssMR << fixed << setprecision(1) << mr;
                                        string MR = ssMR.str();
                                        MR.erase(remove(MR.begin(), MR.end(), '.'), MR.end());

                                        int nFlag, kFlag, nSolutionFlag, mSolutionFlag, solutionVerticesFlag, solutionEdgesFlag, solutionCostFlag, coverFlag, solutionVerticesNotInMainCycleFlag, nonSolutionVerticesInMainCycleFlag, errorFlag;
                                        nFlag = kFlag = nSolutionFlag = mSolutionFlag = solutionVerticesFlag = solutionEdgesFlag = solutionCostFlag = coverFlag = solutionVerticesNotInMainCycleFlag = nonSolutionVerticesInMainCycleFlag = errorFlag = 1;
                                    
                                        ifstream nFlagFile = ifstream("../" + path + "/output/N" + N + "D" + D + "K" + K + "T" + T + "P" + P + "I" + I + "PS" + PS + "MR" + MR + "/nFlag.txt");
                                        if (nFlagFile.is_open()) {
                                            nFlagFile >> nFlag;
                                        }

                                        ifstream kFlagFile = ifstream("../" + path + "/output/N" + N + "D" + D + "K" + K + "T" + T + "P" + P + "I" + I + "PS" + PS + "MR" + MR + "/kFlag.txt");
                                        if (kFlagFile.is_open()) {
                                            kFlagFile >> kFlag;
                                        }

                                        ifstream nSolutionFlagFile = ifstream("../" + path + "/output/N" + N + "D" + D + "K" + K + "T" + T + "P" + P + "I" + I + "PS" + PS + "MR" + MR + "/nSolutionFlag.txt");
                                        if (nSolutionFlagFile.is_open()) {
                                            nSolutionFlagFile >> nSolutionFlag;
                                        }

                                        ifstream mSolutionFlagFile = ifstream("../" + path + "/output/N" + N + "D" + D + "K" + K + "T" + T + "P" + P + "I" + I + "PS" + PS + "MR" + MR + "/mSolutionFlag.txt");
                                        if (mSolutionFlagFile.is_open()) {
                                            mSolutionFlagFile >> mSolutionFlag;
                                        }

                                        ifstream solutionVerticesFlagFile = ifstream("../" + path + "/output/N" + N + "D" + D + "K" + K + "T" + T + "P" + P + "I" + I + "PS" + PS + "MR" + MR + "/solutionVerticesFlag.txt");
                                        if (solutionVerticesFlagFile.is_open()) {
                                            solutionVerticesFlagFile >> solutionVerticesFlag;
                                        }

                                        ifstream solutionEdgesFlagFile = ifstream("../" + path + "/output/N" + N + "D" + D + "K" + K + "T" + T + "P" + P + "I" + I + "PS" + PS + "MR" + MR + "/solutionEdgesFlag.txt");
                                        if (solutionEdgesFlagFile.is_open()) {
                                            solutionEdgesFlagFile >> solutionEdgesFlag;
                                        }

                                        ifstream solutionCostFlagFile = ifstream("../" + path + "/output/N" + N + "D" + D + "K" + K + "T" + T + "P" + P + "I" + I + "PS" + PS + "MR" + MR + "/solutionCostFlag.txt");
                                        if (solutionCostFlagFile.is_open()) {
                                            solutionCostFlagFile >> solutionCostFlag;
                                        }

                                        ifstream coverFlagFile = ifstream("../" + path + "/output/N" + N + "D" + D + "K" + K + "T" + T + "P" + P + "I" + I + "PS" + PS + "MR" + MR + "/coverFlag.txt");
                                        if (coverFlagFile.is_open()) {
                                            coverFlagFile >> coverFlag;
                                        }

                                        ifstream solutionVerticesNotInMainCycleFlagFile = ifstream("../" + path + "/output/N" + N + "D" + D + "K" + K + "T" + T + "P" + P + "I" + I + "PS" + PS + "MR" + MR + "/solutionVerticesNotInMainCycleFlag.txt");
                                        if (solutionVerticesNotInMainCycleFlagFile.is_open()) {
                                            solutionVerticesNotInMainCycleFlagFile >> solutionVerticesNotInMainCycleFlag;
                                        }

                                        ifstream nonSolutionVerticesInMainCycleFlagFile = ifstream("../" + path + "/output/N" + N + "D" + D + "K" + K + "T" + T + "P" + P + "I" + I + "PS" + PS + "MR" + MR + "/nonSolutionVerticesInMainCycleFlag.txt");
                                        if (nonSolutionVerticesInMainCycleFlagFile.is_open()) {
                                            nonSolutionVerticesInMainCycleFlagFile >> nonSolutionVerticesInMainCycleFlag;
                                        }

                                        ifstream errorFlagFile = ifstream("../" + path + "/output/N" + N + "D" + D + "K" + K + "T" + T + "P" + P + "I" + I + "PS" + PS + "MR" + MR + "/errorFlag.txt");
                                        if (errorFlagFile.is_open()) {
                                            errorFlagFile >> errorFlag;
                                        }

                                        cout << n << ',' << d << ',' << k << ',' << t << ',' << p << ',' << i << ',' << ps << ',' << mr << ',' << nFlag << ',' << kFlag << ',' << nSolutionFlag << ',' << mSolutionFlag << ',' << solutionVerticesFlag << ',' << solutionEdgesFlag << ',' << solutionCostFlag << ',' << coverFlag << ',' << solutionVerticesNotInMainCycleFlag << ',' << nonSolutionVerticesInMainCycleFlag << ',' << errorFlag << endl;
                                    }
                                }
                            } else {
                                int nFlag, kFlag, nSolutionFlag, mSolutionFlag, solutionVerticesFlag, solutionEdgesFlag, solutionCostFlag, coverFlag, solutionVerticesNotInMainCycleFlag, nonSolutionVerticesInMainCycleFlag, errorFlag;
                                nFlag = kFlag = nSolutionFlag = mSolutionFlag = solutionVerticesFlag = solutionEdgesFlag = solutionCostFlag = coverFlag = solutionVerticesNotInMainCycleFlag = nonSolutionVerticesInMainCycleFlag = errorFlag = 1;

                                ifstream nFlagFile = ifstream("../" + path + "/output/N" + N + "D" + D + "K" + K + "T" + T + "P" + P + "I" + I + "/nFlag.txt");
                                if (nFlagFile.is_open()) {
                                    nFlagFile >> nFlag;
                                }

                                ifstream kFlagFile = ifstream("../" + path + "/output/N" + N + "D" + D + "K" + K + "T" + T + "P" + P + "I" + I + "/kFlag.txt");
                                if (kFlagFile.is_open()) {
                                    kFlagFile >> kFlag;
                                }

                                ifstream nSolutionFlagFile = ifstream("../" + path + "/output/N" + N + "D" + D + "K" + K + "T" + T + "P" + P + "I" + I + "/nSolutionFlag.txt");
                                if (nSolutionFlagFile.is_open()) {
                                    nSolutionFlagFile >> nSolutionFlag;
                                }

                                ifstream mSolutionFlagFile = ifstream("../" + path + "/output/N" + N + "D" + D + "K" + K + "T" + T + "P" + P + "I" + I + "/mSolutionFlag.txt");
                                if (mSolutionFlagFile.is_open()) {
                                    mSolutionFlagFile >> mSolutionFlag;
                                }

                                ifstream solutionVerticesFlagFile = ifstream("../" + path + "/output/N" + N + "D" + D + "K" + K + "T" + T + "P" + P + "I" + I + "/solutionVerticesFlag.txt");
                                if (solutionVerticesFlagFile.is_open()) {
                                    solutionVerticesFlagFile >> solutionVerticesFlag;
                                }

                                ifstream solutionEdgesFlagFile = ifstream("../" + path + "/output/N" + N + "D" + D + "K" + K + "T" + T + "P" + P + "I" + I + "/solutionEdgesFlag.txt");
                                if (solutionEdgesFlagFile.is_open()) {
                                    solutionEdgesFlagFile >> solutionEdgesFlag;
                                }

                                ifstream solutionCostFlagFile = ifstream("../" + path + "/output/N" + N + "D" + D + "K" + K + "T" + T + "P" + P + "I" + I + "/solutionCostFlag.txt");
                                if (solutionCostFlagFile.is_open()) {
                                    solutionCostFlagFile >> solutionCostFlag;
                                }

                                ifstream coverFlagFile = ifstream("../" + path + "/output/N" + N + "D" + D + "K" + K + "T" + T + "P" + P + "I" + I + "/coverFlag.txt");
                                if (coverFlagFile.is_open()) {
                                    coverFlagFile >> coverFlag;
                                }

                                ifstream solutionVerticesNotInMainCycleFlagFile = ifstream("../" + path + "/output/N" + N + "D" + D + "K" + K + "T" + T + "P" + P + "I" + I + "/solutionVerticesNotInMainCycleFlag.txt");
                                if (solutionVerticesNotInMainCycleFlagFile.is_open()) {
                                    solutionVerticesNotInMainCycleFlagFile >> solutionVerticesNotInMainCycleFlag;
                                }

                                ifstream nonSolutionVerticesInMainCycleFlagFile = ifstream("../" + path + "/output/N" + N + "D" + D + "K" + K + "T" + T + "P" + P + "I" + I + "/nonSolutionVerticesInMainCycleFlag.txt");
                                if (nonSolutionVerticesInMainCycleFlagFile.is_open()) {
                                    nonSolutionVerticesInMainCycleFlagFile >> nonSolutionVerticesInMainCycleFlag;
                                }

                                ifstream errorFlagFile = ifstream("../" + path + "/output/N" + N + "D" + D + "K" + K + "T" + T + "P" + P + "I" + I + "/errorFlag.txt");
                                if (errorFlagFile.is_open()) {
                                    errorFlagFile >> errorFlag;
                                }

                                cout << n << ',' << d << ',' << k << ',' << t << ',' << p << ',' << i << ',' << nFlag << ',' << kFlag << ',' << nSolutionFlag << ',' << mSolutionFlag << ',' << solutionVerticesFlag << ',' << solutionEdgesFlag << ',' << solutionCostFlag << ',' << coverFlag << ',' << solutionVerticesNotInMainCycleFlag << ',' << nonSolutionVerticesInMainCycleFlag << ',' << errorFlag << endl;
                            }

                        }
                    }
                }
            }
        }
    }
    return 0;
}
