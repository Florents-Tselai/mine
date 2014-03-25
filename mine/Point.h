#include <iostream>


using namespace std;


struct Point {
double x, y;
Point() {}
Point(double px, double py) : x(px), y(py) {}

};

ostream&
operator<<(ostream& o, const Point &p) {
    o << "("<< p.x << ',' << p.y << ") ";
    return o;
}

ostream&
operator<<(ostream& o, const Point *p)
{
o << *p;
return o;
}

bool less_x(const Point *a, const Point *b) { return a->x < b->x; }

// Functor for sorting points on y
bool less_y(const Point *a, const Point *b) { return a->y < b->y; }
