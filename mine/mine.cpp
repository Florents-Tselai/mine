#include <iostream>
#include "assert.h"
#include <vector>
#include <algorithm>
#include <iterator>
#include <map>
#include "Point.h"

using namespace std;


template <typename T>
void
PRINT_ELEMENTS(const vector <T> &v)
{
    for (auto elem : v)
        cout << elem;
    cout << endl;
}

map<Point,int>
equipartition_y_axis(const vector<Point> &d, int y) {
    int n = d.size();
    map<Point,int> q;
    PRINT_ELEMENTS(d);
    vector<const Point*> dy(n);
    //TODO REFACTOR: use transform with lambda and back_inserter
    //TODO OPTIMIZATION: exploit y ordering
    for (int i = 0; i < n; i++)
        dy[i] = &(d[i]);
    sort(dy.begin(), dy.end(), less_y);

    double desiredRowSize = double(n) / y;
    int i = 0;
    int sharp = 0;
    int currRow = 0;

    while (i < n) {
        int s=0;
        //TODO REFACTOR: use count_if
        for (int j=0; j<n; j++)
            if (dy[j]->y == dy[i]->y)
                s++;
        double lhs = abs(double(sharp) + double(s) - desiredRowSize);
        double rhs = abs(double(sharp) - desiredRowSize);

        if (sharp != 0 and lhs >= rhs) {
            sharp = 0;
            currRow += 1;
            double temp1 = double(n) - double(i);
            double temp2 = double(y) - double(currRow);
            desiredRowSize = temp1 / temp2;
        }

        for (int j=0; j<s; j++) {
            Point p = *dy[i+j];
            cout << p << ", " << currRow << endl;
            //TODO BUG: insert in q

        }

        i += s;
        sharp += s;

    }

    return q;
}

void test_equipartition_y_axis() {
    Point p[] = {{1,1}, {1,2}, {9,2}, {9,1}, {1,3}, {1,4}, {2,3}, {2,4}, {8,3}, {3,5}, {4,6}, {5,6}, {6,6}, {7,5} };
    vector <Point> points(p, p + 14);
    equipartition_y_axis(points, 3);

}

void run_tests() {
    test_equipartition_y_axis();
}

int main(void) {
    run_tests();
}
