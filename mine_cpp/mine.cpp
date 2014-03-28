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


SortedData
sort_by_x(const Data &d) {
    SortedData sorted_d(d.size());
    for (int i = 0; i < d.size(); i++)
        sorted_d[i] = &(d[i]);
    sort(sorted_d.begin(), sorted_d.end(), less_x);
    return sorted_d;
}

SortedData
sort_by_y(const Data &d) {
    SortedData sorted_d(d.size());
    for (int i = 0; i < d.size(); i++)
        sorted_d[i] = &(d[i]);
    sort(sorted_d.begin(), sorted_d.end(), less_y);
    return sorted_d;
}


int
number_of_clumps(const Partition_map &p) {
    set<int> temp;
    for(auto e : p) {
        temp.insert(e.second);
    }
    return temp.size();
}

Partition_map
equipartition_y_axis(const Data &points, int y) {
    int n = points.size();
    SortedData dy = sort_by_y(points);


    Partition_map q;

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
get_clumps_partition(const SortedData &dx, const Partition_map &q) {
    int n = dx.size();
    Partition_map q_tilde;
    q_tilde.insert(q.cbegin(), q.cend());

    int i = 0;
    int c = -1;

    while (i<n) {
        int s = 1;
        bool flag = false;
        for (int j=i+1; j<n; j++) {
            if (dx[i]->x == dx[j]->x) {
                s += 1;
                if (q_tilde[dx[i]] != q_tilde[dx[j]]) {
                    flag = true;
                }
            }
            else {
                break;
            }
        }
        if (s>1 and flag) {
            for (int j=0; j<s; j++)
                q_tilde.insert(make_pair(dx[i+j], c));
            c -= 1;
        }
        i += s;
    }

    i = 0;
    Partition_map p;
    p.insert(make_pair(dx[0], i));

    for(int j=1; j<n; j++) {
        if (q_tilde[dx[j]] !=q_tilde[dx[j-1]]) {
            i += 1;
        }
        p.insert(make_pair(dx[j], i));
    }
    return p;
}


Partition_map
get_superclumps_partition(const SortedData &dx, Partition_map q, int k_hat) {
    Partition_map p_tilde = get_clumps_partition(dx, q);

    int k = number_of_clumps(p_tilde);
    if (k > k_hat) {
        vector<Point> d_p_tilde(dx.size());
        for(auto point : dx) {
            Point p(0, p_tilde[point]);
        }
        Partition_map p_hat = equipartition_y_axis(d_p_tilde, k_hat);
        Partition_map p;
        for(auto point : dx) {
            int i = p_tilde[point];
            Point temp = Point(0, i);
            int j = p_hat[&temp];
            p.insert(make_pair(point, j));
        }
        return p;
    }

    else {
        return p_tilde;
    }
}



void
test_equipartition_y_axis() {
    cout << "EQUIPARTITION Y AXIS" << endl;
    Data d({{1,1}, {1,2}, {9,2}, {9,1}, {1,3}, {1,4}, {2,3}, {2,4}, {8,3}, {3,5}, {4,6}, {5,6}, {6,6}, {7,5} });

    Partition_map q = equipartition_y_axis(d, 3);
    cout << q;
}

void
test_get_clumps_partition() {
    cout << "GET CLUMPS PARTITION" << endl;
    Data d({{1,1}, {1,2}, {9,2}, {9,1}, {1,3}, {1,4}, {2,3}, {2,4}, {8,3}, {3,5}, {4,6}, {5,6}, {6,6}, {7,5} });
    Partition_map q = equipartition_y_axis(d, 3);
    SortedData dx = sort_by_x(d);
    Partition_map p = get_clumps_partition(dx, q);
    cout << p << endl;
    cout << number_of_clumps(p);

}

void
test_get_superclumps_partition() {
    cout << "GET SUPERCLUMPS PARTITION" << endl;
    Data d({{1,1}, {1,2}, {9,2}, {9,1}, {1,3}, {1,4}, {2,3}, {2,4}, {8,3}, {3,5}, {4,6}, {5,6}, {6,6}, {7,5} });
    Partition_map q = equipartition_y_axis(d, 3);
    SortedData dx = sort_by_x(d);
    Partition_map p = get_superclumps_partition(dx, q, 4);
    cout << p << endl;
    cout << number_of_clumps(p);
}

void
test_sort_by_x() {
    Data d({{1,1}, {1,2}, {9,2}, {9,1}, {1,3}, {1,4}, {2,3}, {2,4}, {8,3}, {3,5}, {4,6}, {5,6}, {6,6}, {7,5} });
    PRINT_ELEMENTS(sort_by_x(d));
}

void
test_sort_by_y() {
    Data d({{1,1}, {1,2}, {9,2}, {9,1}, {1,3}, {1,4}, {2,3}, {2,4}, {8,3}, {3,5}, {4,6}, {5,6}, {6,6}, {7,5} });
    PRINT_ELEMENTS(sort_by_y(d));
}

void
run_tests(){
    test_sort_by_x();
    test_sort_by_y();
    test_equipartition_y_axis();
    test_get_clumps_partition();
    test_get_superclumps_partition();
}

int main(void) {
    run_tests();
}
