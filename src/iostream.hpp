#ifndef SHASTA_IOSTREAM_HPP
#define SHASTA_IOSTREAM_HPP

#include <iostream>

namespace shasta {
    using std::cin;
    using std::cout;
    using std::dec;
    using std::endl;
    using std::flush;
    using std::hex;
    using std::istream;
    using std::ostream;

    // In Shasta we don't use cerr. All log output is to cout.
    // using std::cerr;

    // Output 128-bit integer.
    inline ostream& operator<<(ostream& s, __uint128_t)
    {
        s << "(__uint128_t)";
        return s;
    }
}

#endif
