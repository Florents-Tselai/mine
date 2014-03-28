typedef map<const Point *,int> Partition_map;

typedef vector<Point> Data;

typedef vector<const Point *> SortedData;

ostream&
operator<<(ostream& o, const Partition_map p)
{
    o << "{ " << endl;
    for ( auto p_assignment : p ) {
        cout << p_assignment.first << p_assignment.second << endl;
    }
    o << "}" << endl;
    return o;
}

ostream&
operator<<(ostream& o, const Data d)
{
    o << "{ " << endl;
    for ( auto point : d ) {
        cout << point << ", " ;
    }
    o << "}" << endl;
    return o;
}

template <typename T>
void
PRINT_ELEMENTS(const vector <T> &v)
{
    for (auto elem : v)
        cout << elem;
    cout << endl;
}

void
PRINT_ELEMENTS(const Data &d)
{
    for (auto elem : d)
        cout << elem;
    cout << endl;
}

void
PRINT_ELEMENTS(const SortedData &d)
{
    for (auto elem : d)
        cout << elem;
    cout << endl;
}
