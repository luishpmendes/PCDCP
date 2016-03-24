#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <chrono>

using namespace std;

int main () {
    int n;
    cin >> n;

    unsigned xSeed = chrono::system_clock::now().time_since_epoch().count();
    default_random_engine xGenerator (xSeed);
    uniform_real_distribution <double> xDistribution (0.0, 1.0);

    unsigned ySeed = chrono::system_clock::now().time_since_epoch().count();
    default_random_engine yGenerator (ySeed);
    uniform_real_distribution <double> yDistribution (0.0, 1.0);

    vector < pair <double, double> > points;

    cout << n << ' ' << (n * (n - 1)) / 2 << endl;

    for (int i = 0; i < n; i++) {
        double x = xDistribution(xGenerator);
        double y = yDistribution(yGenerator);
        points.push_back (make_pair(x, y));
        double d = sqrt((x - 0.5) * (x - 0.5) + (y - 0.5) * (y - 0.5));
        cout << round(100 * d) << endl;
    }

    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            double d = sqrt((points[i].first - points[j].first) * (points[i].first - points[j].first) + (points[i].second - points[j].second) * (points[i].second - points[j].second));
            cout << i << ' ' << j << ' ' << round(100 * d) << endl;
        }
    }
}
