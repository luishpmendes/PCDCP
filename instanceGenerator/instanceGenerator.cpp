#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <random>
#include <set>
#include <vector>

using namespace std;

// para garantir conexidade, criar arvore escolhendo pontos n visitados aleatoriamente
// demais arestas, dependem da densidade do grafo

int main () {
    int n, k;
    double d, m;
    cin >> n >> d >> k;

    double dMin = 0.0;
    if (n > 1) {
        dMin = 2.0 / (((double) n) - 1.0);
    }
    double dMax = 1.0;

    if (d < dMin) {
        d = dMin;
    }
    if (d > dMax) {
        d = dMax;
    }

    m = (double) n;
    m *= ((double) n) - 1.0;
    m /= 2.0;
    m *= d;
    m = round(m);

    unsigned xSeed = chrono::system_clock::now().time_since_epoch().count();
    default_random_engine xGenerator (xSeed);
    uniform_real_distribution <double> xDistribution (0.0, 1.0);

    unsigned ySeed = chrono::system_clock::now().time_since_epoch().count();
    default_random_engine yGenerator (ySeed);
    uniform_real_distribution <double> yDistribution (0.0, 1.0);

    vector < pair <double, double> > points;
    set < pair <int, int> > allEdges;
    set < pair <int, int> > chosenEdges;
    vector <int> vAux;

    cout << n << ' ' << m << ' ' << k << endl;

    for (int i = 0; i < n; i++) {
        double x = xDistribution(xGenerator);
        double y = yDistribution(yGenerator);
        points.push_back (make_pair(x, y));
        double d = sqrt((x - 0.5) * (x - 0.5) + (y - 0.5) * (y - 0.5));
        cout << round(100 * d) << endl;
        vAux.push_back(i);
    }

    shuffle(vAux.begin(), vAux.end(), default_random_engine(chrono::system_clock::now().time_since_epoch().count()));

    for (int i = 0; i < (int) vAux.size() - 1; i++) {
        int u = vAux[i];
        int v = vAux[i + 1];
        if (u > v) {
            int aux = u;
            u = v;
            v = aux;
        }
        chosenEdges.insert(make_pair(u, v));
    }
    chosenEdges.insert(make_pair(vAux[vAux.size() - 1], vAux[0]));

    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            if (chosenEdges.find(make_pair(i, j)) == chosenEdges.end()) {
                allEdges.insert(make_pair(i, j));
            }
        }
    }

    vector < pair <int, int> > vAllEdges (allEdges.begin(), allEdges.end());

    shuffle(vAllEdges.begin(), vAllEdges.end(), default_random_engine(chrono::system_clock::now().time_since_epoch().count()));

    for (int i = 0; (int) chosenEdges.size() < m; i++) {
        chosenEdges.insert(vAllEdges[i]);
    }

    vector < pair <int, int> > vChosenEdges (chosenEdges.begin(), chosenEdges.end());

    for (int i = 0; i < m; i++) {
        int u = vChosenEdges[i].first;
        int v = vChosenEdges[i].second;
        double dx = points[u].first - points[v].first;
        double dy = points[u].second - points[v].second;
        double d = sqrt(dx*dx+dy*dy);
        cout << u << ' ' << v << ' ' << round(100 * d) << endl;
    }
}
