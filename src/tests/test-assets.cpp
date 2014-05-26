#include <wibble/tests.h>
#include <stdexcept>
#include <cstdlib>
#include <unistd.h>
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
    Matrix2D<float> dem(400, 512);
    assets.load_dem(dem);
    wassert(actual(dem(0, 0)) == 36);

    assets.configure("SPC", 1389108600);
    assets.load_dem(dem);
    wassert(actual(dem(0, 0)) == 9);
}

template<> template<>
void to::test<2>()
{
    setenv("FIRST_LEVEL_FILE", "../dati/FIRST_LEVEL_corto_GAT_2006_INV", 1);
    Assets assets;
    assets.configure("GAT", 1389108600);
    Matrix2D<unsigned char> m(400, 512);
    assets.load_first_level(m);
}

template<> template<>
void to::test<3>()
{
    setenv("DIR_OUT_PP_BLOC", "testdata", 1);
    Assets assets;
    assets.configure("GAT", 1389108600);
    Matrix2D<unsigned char> m(400, 1024);
    assets.load_first_level_bb_el(m);
}

template<> template<>
void to::test<4>()
{
    setenv("DIR_OUT_PP_BLOC", "testdata", 1);
    Assets assets;
    assets.configure("GAT", 1389108600);
    Matrix2D<unsigned char> m(400, 1024);
    assets.load_first_level_bb_bloc(m);
}

template<> template<>
void to::test<5>()
{
    Assets assets;
    assets.configure("GAT", 1389108600);
    //wassert(actual((int)(fscanf_float_and_close(assets.open_file_hray()) * 100)) == 54406);
}

template<> template<>
void to::test<6>()
{
    Assets assets;
    assets.configure("GAT", 1389108600);
    //wassert(actual((int)(fscanf_float_and_close(assets.open_file_hray_inf()) * 100)) == 54406);
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
    Matrix2D<unsigned char> m(400, 1024);
    assets.load_first_level(m);
}

template<> template<>
void to::test<10>()
{
    static const char* fname = "testdata/test_last_file";
    Assets assets;
    assets.configure("GAT", 1389108600);

    // True because $LAST_FILE is not set
    unlink(fname);
    wassert(actual(assets.save_acq_time(1)).istrue());

    // True because $LAST_FILE does not exist
    unlink(fname);
    setenv("LAST_FILE", fname, 1);
    wassert(actual(assets.save_acq_time(1)).istrue());

    // False because $LAST_FILE exists and has the same date as the current one
    wassert(actual(assets.save_acq_time(1)).isfalse());

    // True because 1 is higher than the value in last_file
    wassert(actual(assets.save_acq_time(2)).istrue());

    // False because $LAST_FILE exists and has the same date as the current one
    wassert(actual(assets.save_acq_time(2)).isfalse());
}

template<> template<>
void to::test<11>()
{
    Assets assets;
    assets.configure("GAT", 1389108600);

    // Result is 0 when the env var is not set
    wassert(actual(assets.read_vpr_heating()) == 0);

    // Unfortunately, our test file also contains 0
    setenv("VPR_HEATING", "testdata/vpr_heat_GAT", 1);
    wassert(actual(assets.read_vpr_heating()) == 0);
}

}
