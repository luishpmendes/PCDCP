#include <iostream>
#include <vector>
#include <string>
#include <fstream>

using namespace std;

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
						int errorFlag = 1;
						ifstream errorFlagFile ("../" + path + "/output/N" + n + "D" + d + "K" + k + "T" + t + "P" + p + "/errorFlag.txt");
						if (errorFlagFile.is_open()) {
							errorFlagFile >> errorFlag;
						}
						cout << "N" << n << "D" << d << "K" << k << "T" << t << "P" << p + ": ";
						if (errorFlag == 0) {
							cout << "OK" << endl;
						} else {
							cout << "ERROR" << endl;
						}
					}
				}
			}
		}
	}
	return 0;
}
