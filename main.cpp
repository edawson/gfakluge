#include "gfakluge.hpp"

int main(){
    GFAKluge gg = GFAKluge();
    gg.parse_gfa_file("reads.gfa");
    return 0;
}
