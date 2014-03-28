#include <wibble/tests.h>
#include <stdexcept>
#include <cstdlib>
#include "assets.h"
#include "utils.h"
#include "logging.h"

using namespace wibble::tests;
using namespace std;
using namespace cumbac;

namespace tut {

struct assets_shar {
    Logging logging;
};
TESTGRP(assets);

namespace {

float fscanf_float_and_close(FILE* in)
{
    float res;
    if (fscanf(in, "%f ", &res) != 1)
        throw runtime_error("failed to read one float");
    if (fclose(in) != 0)
        throw runtime_error("failed to close the file");
    return res;
}

}

template<> template<>
void to::test<1>()
{
    Assets assets;
    assets.configure("GAT", 1389108600);
    wassert(actual(fscanf_float_and_close(assets.open_file_dem())) == 34);

    assets.configure("SPC", 1389108600);
    wassert(actual(fscanf_float_and_close(assets.open_file_dem())) == 9);
}

template<> template<>
void to::test<2>()
{
    setenv("FIRST_LEVEL_FILE", "../dati/FIRST_LEVEL_corto_GAT_2006_INV", 1);
    Assets assets;
    assets.configure("GAT", 1389108600);
    FILE* in = fopen_checked(assets.fname_first_level().c_str(), "rb", "first level");
    wassert(actual(in != 0).istrue());
    fclose(in);
}

template<> template<>
void to::test<3>()
{
    setenv("DIR_OUT_PP_BLOC", "testdata", 1);
    Assets assets;
    assets.configure("GAT", 1389108600);
    FILE* in = fopen_checked(assets.fname_first_level_bb_el().c_str(), "rb", "elev BB");
    wassert(actual(in != 0).istrue());
    fclose(in);
}

template<> template<>
void to::test<4>()
{
    setenv("DIR_OUT_PP_BLOC", "testdata", 1);
    Assets assets;
    assets.configure("GAT", 1389108600);
    FILE* in = fopen_checked(assets.fname_first_level_bb_bloc().c_str(), "rb", "elev BB");
    wassert(actual(in != 0).istrue());
    fclose(in);
}

template<> template<>
void to::test<5>()
{
    Assets assets;
    assets.configure("GAT", 1389108600);
    wassert(actual((int)(fscanf_float_and_close(assets.open_file_hray()) * 100)) == 54406);
}

template<> template<>
void to::test<6>()
{
    Assets assets;
    assets.configure("GAT", 1389108600);
    wassert(actual((int)(fscanf_float_and_close(assets.open_file_hray_inf()) * 100)) == 54406);
}

template<> template<>
void to::test<7>()
{
    setenv("FILE_T", "testdata/temperature.txt", 1);
    Assets assets;
    assets.configure("GAT", 1389108600);
    wassert(actual((int)(assets.read_t_ground() * 100)) == 1010);
}

template<> template<>
void to::test<8>()
{
    setenv("LAST_VPR", "testdata/last_vpr.tmp", 1);
    Assets assets;
    assets.configure("GAT", 1389108600);
    assets.write_last_vpr();
    wassert(actual(assets.read_profile_gap()) == 0);
}

template<> template<>
void to::test<9>()
{
    setenv("FIRST_LEVEL_FILE", "../dati/FIRST_LEVEL_corto+medio_GAT_PRI-EST_2011", 1);
    Assets assets;
    assets.configure("GAT", 1389108600);
    FILE* in = fopen_checked(assets.fname_first_level().c_str(), "rb", "mappa statica");
    wassert(actual(in != 0).istrue());
    fclose(in);
}
}
