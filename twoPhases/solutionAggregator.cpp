#include <iostream>
#include <vector>
#include <string>
#include <fstream>

using namespace std;

typedef unsigned long int ulint;

int main (int argc, char * argv[]) {
    string path;

    if (argc == 2) {
        path = string(argv[1]);
    } else {
        cin >> path;
    }

    vector <string> N;
    N.push_back("10");
    N.push_back("20");
    N.push_back("50");
    N.push_back("100");
    vector <string> D;
    D.push_back("03");
    D.push_back("05");
    D.push_back("07");
    vector <string> K;
    K.push_back("0");
    K.push_back("5");
    K.push_back("10");
    K.push_back("20");
    vector <string> T;
    T.push_back("0");
    T.push_back("1");
    vector <string> P;
    P.push_back("01");
    P.push_back("05");
    P.push_back("10");
    vector <string> I;
    I.push_back("0");
    I.push_back("1");
    I.push_back("2");
    I.push_back("3");
    I.push_back("4");
    I.push_back("5");
    I.push_back("6");
    I.push_back("7");
    I.push_back("8");
    I.push_back("9");

    cout << "n,d,k,t,p,i,objVal1,gap1,elapsedTime1,objVal2,gap2,elapsedTime2" << endl;

    for (vector <string> :: iterator itN = N.begin(); itN != N.end(); itN++) {
        string n = *itN;
        for (vector <string> :: iterator itD = D.begin(); itD != D.end(); itD++) {
            string d = *itD;
            for (vector <string> :: iterator itK = K.begin(); itK != K.end(); itK++) {
                string k = *itK;
                for (vector <string> :: iterator itT = T.begin(); itT != T.end(); itT++) {
                    string t = *itT;
                    for (vector <string> :: iterator itP = P.begin(); itP != P.end(); itP++) {
                        string p = *itP;
                        for (vector <string> :: iterator itI = I.begin(); itI != I.end(); itI++) {
                            string i = *itI;

                            double objVal1 = 0.0;
                            ifstream objVal1File ("./output/N" + n + "D" + d + "K" + k + "T" + t + "P" + p + "I" + i + "/objVal1.txt");
                            if (objVal1File.is_open()) {
                                objVal1File >> objVal1;
                            }

                            double gap1 = 0.0;
                            ifstream gap1File ("./output/N" + n + "D" + d + "K" + k + "T" + t + "P" + p + "I" + i + "/gap1.txt");
                            if (gap1File.is_open()) {
                                gap1File >> gap1;
                            }

                            ulint elapsedTime1 = 0.0;
                            ifstream elapsedTime1File ("./output/N" + n + "D" + d + "K" + k + "T" + t + "P" + p + "I" + i + "/elapsedTime1.txt");
                            if (elapsedTime1File.is_open()) {
                                elapsedTime1File >> elapsedTime1;
                            }

                            double objVal2 = 0.0;
                            ifstream objVal2File ("./output/N" + n + "D" + d + "K" + k + "T" + t + "P" + p + "I" + i + "/objVal2.txt");
                            if (objVal2File.is_open()) {
                                objVal2File >> objVal2;
                            }

                            double gap2 = 0.0;
                            ifstream gap2File ("./output/N" + n + "D" + d + "K" + k + "T" + t + "P" + p + "I" + i + "/gap2.txt");
                            if (gap2File.is_open()) {
                                gap2File >> gap2;
                            }

                            ulint elapsedTime2 = 0.0;
                            ifstream elapsedTime2File ("./output/N" + n + "D" + d + "K" + k + "T" + t + "P" + p + "I" + i + "/elapsedTime2.txt");
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
