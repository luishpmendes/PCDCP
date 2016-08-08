#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <random>
#include <set>
#include <vector>

using namespace std;

typedef unsigned long int ulint;

int main (int argc, char * argv[]) {
    int n, k, t, r;
    double density, m;

    if (argc == 6) {
        n = atoi(argv[1]);
        density = atof(argv[2]);
        k = atoi(argv[3]);
        t = atoi(argv[4]);
        r = atoi(argv[5]);
    } else {
        cin >> n >> density >> k >> t >> r;
    }

    double minDensity = 0.0;
    if (n > 1) {
        minDensity = 2.0 / (((double) n) - 1.0);
    }
    double maxDensity = 1.0;

    // validate density value
    if (density < minDensity) {
        density = minDensity;
    }
    if (density > maxDensity) {
        density = maxDensity;
    }

    // m = density*((n*(n-1))/2)
    m = (double) n;
    m *= ((double) n) - 1.0;
    m /= 2.0;
    m *= density;
    m = round(m);

    unsigned xSeed = chrono::system_clock::now().time_since_epoch().count();
    default_random_engine xGenerator (xSeed);
    uniform_real_distribution <double> xDistribution (0.0, 1.0);

    unsigned ySeed = chrono::system_clock::now().time_since_epoch().count();
    default_random_engine yGenerator (ySeed);
    uniform_real_distribution <double> yDistribution (0.0, 1.0);

    vector < pair <double, double> > points (n);
    vector <int> vAux (n);

    int closest = 0;
    int distClosest = 200;
    int farthest = 0;
    int distFarthest = 0;
    int root;

    // creating points array
    // calculating the vertices that are closest and farthest from the middle (0.5, 0.5)
    for (int i = 0; i < n; i++) {
        double x = xDistribution(xGenerator);
        double y = yDistribution(yGenerator);
        points[i] = (make_pair(x, y));
        double dx = x - 0.5;
        double dy = y - 0.5;
        int dist = round(100 * sqrt(dx * dx + dy * dy));
        if (distClosest > dist) {
            closest = i;
            distClosest = dist;
        }
        if (distFarthest < dist) {
            farthest = i;
            distFarthest = dist;
        }
        vAux[i] = i;
    }

    // setting reference point for vertices penalty
    int referencePoint = closest;
    if (t == 1) {
        referencePoint = farthest;
    }

    // setting root vertex
    if (r == 0) {
        root = closest;
    } else if (r == 0) {
        root = farthest;
    } else {
        unsigned rootSeed = chrono::system_clock::now().time_since_epoch().count();
        default_random_engine rootGenerator (rootSeed);
        uniform_int_distribution <double> rootDistribution (0, n-1);
        root = rootDistribution(rootGenerator);
    }

    // number of vertices, number of edges in the complete graph, number of 'chosen' edges and neighborhood radio
    cout << n << ' ' << (n * (n - 1)) / 2 << ' ' << m << ' ' << k << ' ' << root << endl;

    // printing vertices' coordinates and its penalty
    for (int i = 0; i < n; i++) {
        double x = points[i].first;
        double y = points[i].second;
        double x0 = points[referencePoint].first;
        double y0 = points[referencePoint].second;
        double dx = x - x0;
        double dy = y - y0;
        int dist = round(100 * sqrt(dx * dx + dy * dy));
        cout << x << ' ' << y << ' ' << dist << endl;
    }

    // generating and printing the edges of the complete graph and its weights
    set < pair <int, int> > allEdges;
    for (int u = 0; u < n; u++) {
        for (int v = u + 1; v < n; v++) {
            double dx = points[u].first - points[v].first;
            double dy = points[u].second - points[v].second;
            int dist = round(100 * sqrt(dx * dx + dy * dy));
            cout << u << ' ' << v << ' ' << dist << endl;
            allEdges.insert(make_pair(u, v));
        }
    }

    // choosing the edges to form the basic euclidean tour in the graph
    set < pair <int, int> > chosenEdges;
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

    // put all edges in a array and shuffle its
    vector < pair <int, int> > vAllEdges (allEdges.begin(), allEdges.end());
    shuffle(vAllEdges.begin(), vAllEdges.end(), default_random_engine(chrono::system_clock::now().time_since_epoch().count()));

    // inserting random edges in the graph until it reaches the desired density
    for (int i = 0; (int) chosenEdges.size() < m; i++) {
        chosenEdges.insert(vAllEdges[i]);
    }

    vector < pair <int, int> > vChosenEdges (chosenEdges.begin(), chosenEdges.end());

    // printing the graph edges and its weigths
    for (int i = 0; i < m; i++) {
        int u = vChosenEdges[i].first;
        int v = vChosenEdges[i].second;
        double dx = points[u].first - points[v].first;
        double dy = points[u].second - points[v].second;
        int dist = round(100 * sqrt(dx * dx + dy * dy));
        cout << u << ' ' << v << ' ' << dist << endl;
    }

    return 0;
}
