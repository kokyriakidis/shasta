#include "ShortBaseSequence.hpp"
#include "ShortBaseSequenceEdit.hpp"
#include "algorithm.hpp"
#include "SHASTA_ASSERT.hpp"
using namespace shasta;

#include <iomanip>



void shasta::testShortBaseSequence()
{
    const string s = "CGATCTTGAT";
    const uint64_t k = s.size();
    ShortBaseSequence64 x;
    for(uint64_t i=0; i<k; i++) {
        x.set(i, Base::fromCharacter(s[i]));
    }

    vector<ShortBaseSequence64> v;
    applySingleEdit(x, k, v);

    cout << "x:" << endl;
    x.write(cout, k);
    cout << endl;
    cout << "y:" << endl;


    for(const auto& y: v) {
        y.write(cout, k);
        cout << endl;
    }
}
