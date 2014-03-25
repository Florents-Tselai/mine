typedef map<const Point *,int> Partition_map;

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

template <typename T>
void
PRINT_ELEMENTS(const vector <T> &v)
{
    for (auto elem : v)
        cout << elem;
    cout << endl;
}
