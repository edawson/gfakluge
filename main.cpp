#include "gfakluge.hpp"

using namespace std;
using namespace gfak;
int main(){
    GFAKluge gg = GFAKluge();
    gg.parse_gfa_file("reads.gfa");
    cout << gg.to_string() << endl;
    return 0;
}
