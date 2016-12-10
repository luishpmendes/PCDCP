#include <iostream>
#include <vector>
#include <string>
#include <fstream>

using namespace std;

typedef unsigned long int ulint;

int main (int argc, char * argv[]) {
    vector <string> N;
    N.push_back("10");
    N.push_back("20");
    N.push_back("50");
    N.push_back("100");
    N.push_back("250");
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
    vector <string> A;
    A.push_back("03");
    A.push_back("05");
    A.push_back("07");

    cout << "n,d,k,t,p,i,a,objVal,elapsedTime" << endl;

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
                            for (vector <string> :: iterator itA = A.begin(); itA != A.end(); itA++) {
                                string a = *itA;

                                double objVal = 0.0;
                                ifstream objValFile ("./output/N" + n + "D" + d + "K" + k + "T" + t + "P" + p + "I" + i + "A" + a + "/objVal.txt");
                                if (objValFile.is_open()) {
                                    objValFile >> objVal;
                                }

                                ulint elapsedTime = 0.0;
                                ifstream elapsedTimeFile ("./output/N" + n + "D" + d + "K" + k + "T" + t + "P" + p + "I" + i + "A" + a + "/elapsedTime.txt");
                                if (elapsedTimeFile.is_open()) {
                                    elapsedTimeFile >> elapsedTime;
                                }

                                cout << n << ',' << d << ',' << k << ',' << t << ',' << p << ',' << i << ',' << a << ',' << objVal << ',' << elapsedTime << endl;
                            }
                        }
                    }
                }
            }
        }
    }
    return 0;
}
