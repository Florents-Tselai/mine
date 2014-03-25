#include <iostream>
#include "assert.h"
#include <vector>
#include <algorithm>
#include <iterator>
#include <map>
#include <set>
#include "Point.h"
#include "mine.h"

using namespace std;


Partition_map
equipartition_y_axis(const vector<Point> &d, int y) {
    int n = d.size();
    vector<const Point*> dy(n);
    //TODO REFACTOR: use transform with lambda and back_inserter
    //TODO OPTIMIZATION: exploit y ordering
    for (int i = 0; i < n; i++)
        dy[i] = &(d[i]);
    sort(dy.begin(), dy.end(), less_y);

    map<const Point*,int> q;

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

        for (int j=0; j<s; j++)
            q.insert(make_pair(dy[i+j], currRow));

        i += s;
        sharp += s;

    }

    return q;
}


Partition_map
get_clumps_partition(const vector<Point> &d, const Partition_map &q) {
    int n = d.size();
    Partition_map q_tilde;
    q_tilde.insert(q.begin(), q.end());

    vector<const Point*> dx(n);
    //TODO REFACTOR: use transform with lambda and back_inserter
    //TODO OPTIMIZATION: exploit y ordering
    for (int i = 0; i < n; i++)
        dx[i] = &(d[i]);
    sort(dx.begin(), dx.end(), less_x);

    int i = 0;
    int c = -1;

    while (i<n) {
        int s = 1;
        bool flag = false;
        for (int j=i+1; j<n; j++) {
            if (dx[i]->x == dx[j]->x) {
                s++;
                if (q_tilde[dx[i]] != q_tilde[dx[j]]) {
                    flag = true;
                }
            }
            else {
                break;
            }
        }
        if ((s>1) & flag) {
            for (int j=0; j<s; j++)
                q_tilde.insert(make_pair(dx[i+j], c));
            c--;
        }
        i += s;
    }

    i = 0;
    Partition_map p;
    p.insert(make_pair(dx[0], i));

    for(int j=1; j<n; j++) {
        if (q_tilde[dx[j]] !=q_tilde[dx[j-1]]) {
            i++;
        }
        p.insert(make_pair(dx[j], i));
    }
    return p;
}

Partition_map
get_superclumps_partition(const vector<Point> &dx, Partition_map q, int k_hat) {
    Partition_map p_tilde = get_clumps_partition(dx, q);

    set<int> temp;
    for(auto e : q) {
        temp.insert(e.second);
    }
    int k = temp.size();

    if (k>k_hat) {
        vector<Point> d_p_tilde(dx.size());
        for(auto point : dx) {
            Point p(0, p_tilde[&point]);
        }
        Partition_map p_hat = equipartition_y_axis(d_p_tilde, k_hat);
        Partition_map p;
        for(auto point : dx) {
            int i = p_tilde[&point];
            Point temp = Point(0, i);
            int j = p_hat[&temp];
            p.insert(make_pair(&point, j));
        }
        return p;
    }

    else {
        return p_tilde;
    }
}

void test_equipartition_y_axis() {
    Point p[] = {{1,1}, {1,2}, {9,2}, {9,1}, {1,3}, {1,4}, {2,3}, {2,4}, {8,3}, {3,5}, {4,6}, {5,6}, {6,6}, {7,5} };
    vector <Point> points(p, p + 14);
    map<const Point*, int> q = equipartition_y_axis(points, 3);
}

void test_get_clumps_partition() {
    Point points[] = {{1,1}, {1,2}, {9,2}, {9,1}, {1,3}, {1,4}, {2,3}, {2,4}, {8,3}, {3,5}, {4,6}, {5,6}, {6,6}, {7,5} };
    vector <Point> data(points, points + 14);
    Partition_map q = equipartition_y_axis(data, 3);
    Partition_map p = get_clumps_partition(data, q);
    cout << p;
}

void run_tests() {
    test_equipartition_y_axis();
    test_get_clumps_partition();
}

int main(void) {
    run_tests();
}
