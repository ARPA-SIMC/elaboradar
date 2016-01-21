#include <radarelab/utils/tests.h>
#include <radarelab/utils.h>
#include <vector>

using namespace radarelab::utils::tests;
using namespace radarelab;
using namespace std;

namespace {

class Tests : public TestCase
{
    using TestCase::TestCase;

    void register_tests() override;
} test("utils");

void Tests::register_tests() {

add_method("str_split", []() {
    char buf[] = ",,, 1, 2, 3  ";
    vector<string> split;
    str_split(buf, ", ", [&](const char* tok) { split.push_back(tok); });
    wassert(actual(split.size()) == 3);
    wassert(actual(split[0]) == "1");
    wassert(actual(split[1]) == "2");
    wassert(actual(split[2]) == "3");
});

}

}
