#include <wibble/tests.h>
#include <vector>
#include "utils.h"

using namespace wibble::tests;
using namespace elaboradar;
using namespace std;

namespace tut {

struct utils_shar {
    Logging logging;
};
TESTGRP(utils);

template<> template<>
void to::test<1>()
{
    char buf[] = ",,, 1, 2, 3  ";
    vector<string> split;
    str_split(buf, ", ", [&](const char* tok) { split.push_back(tok); });
    wassert(actual(split.size()) == 3);
    wassert(actual(split[0]) == "1");
    wassert(actual(split[1]) == "2");
    wassert(actual(split[2]) == "3");
}

}
