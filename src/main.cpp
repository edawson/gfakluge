#include "gfakluge.hpp"

using namespace std;
using namespace gfak;
int main(){
    GFAKluge gg = GFAKluge();
    gg.parse_gfa_file("../data/q_redundant.gfa");
    cerr << "GFA file parsed successfully" << endl;
    cerr << "Testing to_string" << endl;
    cerr << gg.to_string();

    ofstream ff;
    ff.open("../data/q_test.gfa");
    ff << gg;
    ff.close();
    
    ff.open("../data/test_test.gfa");

    GFAKluge gt;
    gt.parse_gfa_file("../data/test.gfa");
    ff << gt;
    ff.close();

    //GFAKluge gc;
    //gc.parse_gfa_file("comment.gfa");


    return 0;
}
