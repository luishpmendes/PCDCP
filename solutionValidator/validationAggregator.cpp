#include <iostream>
#include <vector>
#include <string>
#include <fstream>

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
    vP.push_back(0);
    vP.push_back(1);
    vP.push_back(2);
    vP.push_back(3);
    vP.push_back(4);
    vP.push_back(5);
    vP.push_back(6);
    vP.push_back(7);
    vP.push_back(8);
    vP.push_back(9);

    cout << "n,d,k,t,p,i,error" << endl;

    for (vector <int> :: iterator itN = vN.begin(); itN != vN.end(); itN++) {
        int n = *itN;
        for (vector <double> :: iterator itD = vD.begin(); itD != vD.end(); itD++) {
            double d = *itD;
            stringstream ssD;
            ssD << fixed << setprecision(1) << d;
            string D = ssD.str();
            D.erase(remove(D.begin(), D.end(), '.'), D.end());
            for (vector <int> :: iterator itK = vK.begin(); itK != vK.end(); itK++) {
                int k = *itK;
                for (vector <int> :: iterator itT = vT.begin(); itT != vT.end(); itT++) {
                    int t = *itT;
                    for (vector <double> :: iterator itP = vP.begin(); itP != vP.end(); itP++) {
                        double p = *itP;
                        stringstream ssP;
                        ssP << fixed << setprecision(1) << p;
                        string P = ssP.str();
                        P.erase(remove(P.begin(), P.end(), '.'), P.end());
                        for (vector <int> :: iterator itI = vI.begin(); itI != vI.end(); itI++) {
                            int i = *itI;

                            int errorFlag = 1;
                            ifstream errorFlagFile ("../" + path + "/output/N" + n + "D" + D + "K" + k + "T" + t + "P" + P + "I" + i + "/errorFlag.txt");
                            if (errorFlagFile.is_open()) {
                                errorFlagFile >> errorFlag;
                            }

                            out << n << ',' << d << ',' << k << ',' << t << ',' << p << ',' << i << ',' << errorFlag << endl;
                        }
                    }
                }
            }
        }
    }
    return 0;
}
