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
    vector <int> vI;
    vI.push_back(0);
    vector <double> vA;
    vA.push_back(0.3);
    vA.push_back(0.5);
    vA.push_back(0.7);

    cout << "n,d,k,t,i,a,objVal,elapsedTime" << endl;

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
                    for (vector <int> :: iterator itI = vI.begin(); itI != vI.end(); itI++) {
                        int i = *itI;
                        stringstream ssI;
                        ssI << i;
                        string I = ssI.str();
                        for (vector <double> :: iterator itA = vA.begin(); itA != vA.end(); itA++) {
                            double a = *itA;
                            stringstream ssA;
                            ssA << fixed << setprecision(1) << a;
                            string A = ssA.str();
                            A.erase(remove(A.begin(), A.end(), '.'), A.end());

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

                            cout << n << ',' << d << ',' << k << ',' << t << ',' << i << ',' << a << ',' << objVal << ',' << elapsedTime << endl;
                        }
                    }
                }
            }
        }
    }
    return 0;
}
