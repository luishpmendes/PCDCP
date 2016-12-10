#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>

using namespace std;

typedef unsigned long int ulint;

int main (int argc, char * argv[]) {
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

    cout << "n,d,k,t,p,i,objVal,gap,elapsedTime" << endl;

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
                        for (vector <int> :: iterator itI = vI.begin(); itI != vI.end(); itI++) {
                            int i = *itI;
                            stringstream ssI;
                            ssK << i;
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

                            cout << n << ',' << d << ',' << k << ',' << t << ',' << p << ',' << i << ',' << objVal << ',' << gap << ',' << elapsedTime << endl;
                        }
                    }
                }
            }
        }
    }
    return 0;
}
