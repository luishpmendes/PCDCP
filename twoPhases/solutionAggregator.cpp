#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>

using namespace std;

typedef unsigned long int ulint;

int main () {
    vector <int> vN;
    vN.push_back(10);
    vN.push_back(20);
    vN.push_back(50);
    vN.push_back(100);
    vN.push_back(250);
    vector <double> vD;
    vD.push_back(0.3);
    vD.push_back(0.5);
    vD.push_back(0.7);
    vector <int> vK;
    vK.push_back(0);
    vK.push_back(5);
    vK.push_back(10);
    vK.push_back(20);
    vector <int> vT;
    vT.push_back(0);
    vT.push_back(1);
    vector <double> vP;
    vP.push_back(0.1);
    vP.push_back(0.5);
    vP.push_back(1.0);
    vector <int> vI;
    vI.push_back(0);
    vI.push_back(1);
    vI.push_back(2);
    vI.push_back(3);
    vI.push_back(4);
    vI.push_back(5);
    vI.push_back(6);
    vI.push_back(7);
    vI.push_back(8);
    vI.push_back(9);

    cout << "n,d,k,t,p,i,objVal1,gap1,elapsedTime1,objVal2,gap2,elapsedTime2" << endl;

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
                    ssK << t;
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

                            double objVal1 = 0.0;
                            ifstream objVal1File ("./output/N" + N + "D" + D + "K" + K + "T" + T + "P" + P + "I" + I + "/objVal1.txt");
                            if (objVal1File.is_open()) {
                                objVal1File >> objVal1;
                            }

                            double gap1 = 0.0;
                            ifstream gap1File ("./output/N" + N + "D" + D + "K" + K + "T" + T + "P" + P + "I" + I + "/gap1.txt");
                            if (gap1File.is_open()) {
                                gap1File >> gap1;
                            }

                            ulint elapsedTime1 = 0.0;
                            ifstream elapsedTime1File ("./output/N" + N + "D" + D + "K" + K + "T" + T + "P" + P + "I" + I + "/elapsedTime1.txt");
                            if (elapsedTime1File.is_open()) {
                                elapsedTime1File >> elapsedTime1;
                            }

                            double objVal2 = 0.0;
                            ifstream objVal2File ("./output/N" + N + "D" + D + "K" + K + "T" + T + "P" + P + "I" + I + "/objVal2.txt");
                            if (objVal2File.is_open()) {
                                objVal2File >> objVal2;
                            }

                            double gap2 = 0.0;
                            ifstream gap2File ("./output/N" + N + "D" + D + "K" + K + "T" + T + "P" + P + "I" + I + "/gap2.txt");
                            if (gap2File.is_open()) {
                                gap2File >> gap2;
                            }

                            ulint elapsedTime2 = 0.0;
                            ifstream elapsedTime2File ("./output/N" + N + "D" + D + "K" + K + "T" + T + "P" + P + "I" + I + "/elapsedTime2.txt");
                            if (elapsedTime2File.is_open()) {
                                elapsedTime2File >> elapsedTime2;
                            }

                            cout << n << ',' << d << ',' << k << ',' << t << ',' << p << ',' << i << ',' << objVal1 << ',' << gap1 << ',' << elapsedTime1 << ',' << objVal2 << ',' << gap2 << ',' << elapsedTime2 << endl;
                        }
                    }
                }
            }
        }
    }
    return 0;
}
