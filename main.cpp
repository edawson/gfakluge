#include "gfakluge.hpp"

using namespace std;
using namespace gfak;
int main(){
    GFAKluge gg = GFAKluge();
    gg.parse_gfa_file("test.gfa");
    cerr << "GFA file parsed successfully" << endl;
    cerr << "Testing to_string" << endl;
    cout << gg.to_string() << endl;
    cerr << "If it looks like GFA, it worked!" << endl;
    cerr << "Testing ostream (<<) operator..." << endl;
    cout << gg << endl;
    cerr << "If it looks like GFA, it worked!" << endl;
    return 0;
}
