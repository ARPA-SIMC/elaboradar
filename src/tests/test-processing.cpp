#include <wibble/tests.h>
#include "cum_bac.h"
#include "logging.h"
#include <stdio.h>
#include <vector>

using namespace wibble::tests;
using namespace cumbac;

namespace tut {

struct process_shar {
    Logging logging;
};
TESTGRP(process);

namespace {

template<typename T>
struct ArrayStats
{
    bool first;
    T min;
    T max;
    double avg;

    ArrayStats() : first(true) {}

    template<int A, int B>
    void fill2(const T (&arr)[A][B])
    {
        for (int i = 0; i < A; ++i)
            for (int j = 0; j < B; ++j)
            {
                if (first)
                {
                    min = arr[i][j];
                    max = arr[i][j];
                    avg = min / (A * B);
                    first = false;
                }
                else
                {
                    if (arr[i][j] < min)
                        min = arr[i][j];
                    if (arr[i][j] > max)
                        max = arr[i][j];
                    avg += (double)arr[i][j] / (A * B);
                }
            }
    }

    template<int A, int B, int C>
    void fill3(const T (&arr)[A][B][C])
    {
        for (int i = 0; i < A; ++i)
            for (int j = 0; j < B; ++j)
                for (int k = 0; k < C; ++k)
                {
                    if (first)
                    {
                        min = arr[i][j][k];
                        max = arr[i][j][k];
                        avg = min / (A * B * C);
                        first = false;
                    }
                    else
                    {
                        if (arr[i][j][k] < min)
                            min = arr[i][j][k];
                        if (arr[i][j][k] > max)
                            max = arr[i][j][k];
                        avg += (double)arr[i][j][k] / (A * B * C);
                    }
                }
    }
};

}

template<> template<>
void to::test<3>()
{
    // Test elabora_dato, con tutti i do_* a false
    static const char* fname = "testdata/DBP2_070120141530_GATTATICO";

    setenv("FIRST_LEVEL_DIM_FILE", "../dati/FL_2006.DIM", 1);
    setenv("FIRST_LEVEL_FILE", "../dati/FIRST_LEVEL_corto_GAT_2006_INV", 1);
    setenv("DIR_OUT_PP_BLOC", "testdata", 1);
    setenv("FILE_T", "testdata/temperature.txt", 1);

    CUM_BAC* cb = new CUM_BAC;
    cb->read_sp20_volume(fname, "GAT", 0);
    cb->setup_elaborazione(fname, "GAT");

    wassert(actual(cb->t_ground) == NODATAVPR);

    int ier = cb->elabora_dato();
    wassert(actual(ier) == 0);

    // Check results
    VolumeStats stats;
    cb->volume.compute_stats(stats);
    wassert(actual(stats.count_zeros[0]) == 0);
    wassert(actual(stats.count_zeros[1]) == 0);
    wassert(actual(stats.count_zeros[2]) == 0);
    wassert(actual(stats.count_zeros[3]) == 0);
    wassert(actual(stats.count_zeros[4]) == 0);
    wassert(actual(stats.count_ones[0]) == 192288);
    wassert(actual(stats.count_ones[1]) == 193052);
    wassert(actual(stats.count_ones[2]) == 194525);
    wassert(actual(stats.count_ones[3]) == 196576);
    wassert(actual(stats.count_ones[4]) == 196160);
    wassert(actual(stats.count_others[0]) == 5312);
    wassert(actual(stats.count_others[1]) == 4548);
    wassert(actual(stats.count_others[2]) == 3075);
    wassert(actual(stats.count_others[3]) == 1024);
    wassert(actual(stats.count_others[4]) == 1440);
    wassert(actual(stats.sum_others[0]) == 493518);
    wassert(actual(stats.sum_others[1]) == 358463);
    wassert(actual(stats.sum_others[2]) == 220768);
    wassert(actual(stats.sum_others[3]) ==  41677);
    wassert(actual(stats.sum_others[4]) ==  78321);

    cb->caratterizzo_volume();

    ArrayStats<unsigned char> qual_stats;
    qual_stats.fill3(cb->qual);
    wassert(actual((unsigned)qual_stats.first).isfalse());
    wassert(actual((unsigned)qual_stats.min) == 0);
    wassert(actual((unsigned)qual_stats.max) == 99);
    wassert(actual((unsigned)(qual_stats.avg * 100)) == 4234);

    ArrayStats<unsigned char> vpr_stats;
    vpr_stats.fill3(cb->flag_vpr);
    wassert(actual((unsigned)vpr_stats.first).isfalse());
    wassert(actual((unsigned)vpr_stats.min) == 0);
    wassert(actual((unsigned)vpr_stats.max) == 0);
    wassert(actual((unsigned)(vpr_stats.avg * 100)) == 0);

    ArrayStats<unsigned char> top_stats;
    top_stats.fill2(cb->top);
    wassert(actual((unsigned)top_stats.first).isfalse());
    wassert(actual((unsigned)top_stats.min) == 0);
    wassert(actual((unsigned)top_stats.max) == 35);
    wassert(actual((unsigned)(top_stats.avg * 100)) == 5);

    delete cb;
}

template<> template<>
void to::test<4>()
{
    // Test elabora_dato, con tutti i do_* a true
    static const char* fname = "testdata/DBP2_070120141530_GATTATICO";

    setenv("FIRST_LEVEL_DIM_FILE", "../dati/FL_2006.DIM", 1);
    setenv("FIRST_LEVEL_FILE", "../dati/FIRST_LEVEL_corto_GAT_2006_INV", 1);
    setenv("DIR_OUT_PP_BLOC", "testdata", 1);
    setenv("FILE_T", "testdata/temperature.txt", 1);

    CUM_BAC* cb = new CUM_BAC;
    cb->do_quality = true;
    cb->do_beamblocking = true;
    cb->do_declutter = true;
    cb->do_bloccorr = true;
    cb->do_vpr = true;
    cb->do_class = true;
    cb->read_sp20_volume(fname, "GAT", 0);
    cb->setup_elaborazione(fname, "GAT");

    wassert(actual((int)(cb->t_ground * 100)) == 1010);

    int ier = cb->elabora_dato();
    wassert(actual(ier) == 0);

    // Check results
    VolumeStats stats;
    cb->volume.compute_stats(stats);
    wassert(actual(stats.count_zeros[0]) == 0);
    wassert(actual(stats.count_zeros[1]) == 0);
    wassert(actual(stats.count_zeros[2]) == 0);
    wassert(actual(stats.count_zeros[3]) == 0);
    wassert(actual(stats.count_zeros[4]) == 0);
    wassert(actual(stats.count_ones[0]) == 176003);
    wassert(actual(stats.count_ones[1]) == 190293);
    wassert(actual(stats.count_ones[2]) == 193588);
    wassert(actual(stats.count_ones[3]) == 196292);
    wassert(actual(stats.count_ones[4]) == 196160);
    wassert(actual(stats.count_others[0]) == 21597);
    wassert(actual(stats.count_others[1]) ==  7307);
    wassert(actual(stats.count_others[2]) ==  4012);
    wassert(actual(stats.count_others[3]) ==  1308);
    wassert(actual(stats.count_others[4]) ==  1440);
    wassert(actual(stats.sum_others[0]) == 1538560);
    wassert(actual(stats.sum_others[1]) ==  478421);
    wassert(actual(stats.sum_others[2]) ==  246658);
    wassert(actual(stats.sum_others[3]) ==   45968);
    wassert(actual(stats.sum_others[4]) ==   78321);

    cb->caratterizzo_volume();

    ArrayStats<unsigned char> qual_stats;
    qual_stats.fill3(cb->qual);
    wassert(actual((unsigned)qual_stats.first).isfalse());
    wassert(actual((unsigned)qual_stats.min) == 0);
    wassert(actual((unsigned)qual_stats.max) == 99);
    wassert(actual((unsigned)(qual_stats.avg * 100)) == 5955);

    ArrayStats<unsigned char> vpr_stats;
    vpr_stats.fill3(cb->flag_vpr);
    wassert(actual((unsigned)vpr_stats.first).isfalse());
    wassert(actual((unsigned)vpr_stats.min) == 0);
    wassert(actual((unsigned)vpr_stats.max) == 1);
    wassert(actual((unsigned)(vpr_stats.avg * 100)) == 92);

    ArrayStats<unsigned char> top_stats;
    top_stats.fill2(cb->top);
    wassert(actual((unsigned)top_stats.first).isfalse());
    wassert(actual((unsigned)top_stats.min) == 0);
    wassert(actual((unsigned)top_stats.max) == 35);
    wassert(actual((unsigned)(top_stats.avg * 100)) == 10);

    cb->classifica_rain();

    ArrayStats<unsigned char> stratiform_stats;
    stratiform_stats.fill2(cb->stratiform);
    wassert(actual((unsigned)stratiform_stats.first).isfalse());
    wassert(actual((unsigned)stratiform_stats.min) == 0);
    wassert(actual((unsigned)stratiform_stats.max) == 0);
    wassert(actual((unsigned)(stratiform_stats.avg * 100)) == 0);

    delete cb;
}

template<> template<>
void to::test<5>()
{
    // Test elabora_dato, con tutti i do_* a true
    static const char* fname = "testdata/DBP2_070120141530_GATTATICO";

    setenv("FIRST_LEVEL_DIM_FILE", "../dati/FL_2006.DIM", 1);
    setenv("FIRST_LEVEL_FILE", "../dati/FIRST_LEVEL_corto_GAT_2006_INV", 1);
    setenv("DIR_OUT_PP_BLOC", "testdata", 1);

    CUM_BAC* cb = new CUM_BAC;
    cb->do_quality = true;
    cb->do_beamblocking = true;
    cb->do_declutter = false;
    cb->do_bloccorr = true;
    cb->do_vpr = true;
    cb->read_sp20_volume(fname, "GAT", 0);
    cb->setup_elaborazione(fname, "GAT");

    int ier = cb->elabora_dato();
    wassert(actual(ier) == 0);

    // Check results
    VolumeStats stats;
    cb->volume.compute_stats(stats);
    //stats.print(stdout);

    wassert(actual(stats.count_zeros[0]) == 0);
    wassert(actual(stats.count_zeros[1]) == 0);
    wassert(actual(stats.count_zeros[2]) == 0);
    wassert(actual(stats.count_zeros[3]) == 0);
    wassert(actual(stats.count_zeros[4]) == 0);
    wassert(actual(stats.count_ones[0]) == 193091);
    wassert(actual(stats.count_ones[1]) == 193929);
    wassert(actual(stats.count_ones[2]) == 194526);
    wassert(actual(stats.count_ones[3]) == 196576);
    wassert(actual(stats.count_ones[4]) == 196160);
    wassert(actual(stats.count_others[0]) ==  4509);
    wassert(actual(stats.count_others[1]) ==  3671);
    wassert(actual(stats.count_others[2]) ==  3074);
    wassert(actual(stats.count_others[3]) ==  1024);
    wassert(actual(stats.count_others[4]) ==  1440);
    wassert(actual(stats.sum_others[0]) ==  387010);
    wassert(actual(stats.sum_others[1]) ==  276931);
    wassert(actual(stats.sum_others[2]) ==  220729);
    wassert(actual(stats.sum_others[3]) ==   41677);
    wassert(actual(stats.sum_others[4]) ==   78321);

    delete cb;
}

template<> template<>
void to::test<6>()
{
    // versione BB che corrisponde al parametro algo_corto
    static const char* fname = "testdata/DBP2_070120141530_GATTATICO";
    setenv("FIRST_LEVEL_DIM_FILE", "../dati/FL_2006.DIM", 1);
    setenv("FIRST_LEVEL_FILE", "../dati/FIRST_LEVEL_corto_GAT_2006_INV", 1);
    setenv("DIR_OUT_PP_BLOC", "testdata", 1);

    CUM_BAC* cb = new CUM_BAC;
    cb->do_quality = true;
    cb->do_beamblocking = true;
    cb->do_declutter = true;
    cb->do_class = false;
    cb->do_bloccorr = false;
    cb->do_vpr = false;
    cb->read_sp20_volume(fname, "GAT", 0);
    cb->setup_elaborazione(fname, "GAT");

    int ier = cb->elabora_dato();
    wassert(actual(ier) == 0);

    cb->caratterizzo_volume();
    cb->creo_cart();
    cb->creo_cart_z_lowris();

    // TODO: scrivo_out_file_bin

    delete cb;
}

template<> template<>
void to::test<7>()
{
    // versione BB_VPR che corrisponde al parametro algo_corto_dev
    static const char* fname = "testdata/DBP2_070120141530_GATTATICO";
    setenv("FIRST_LEVEL_DIM_FILE", "../dati/FL_2006.DIM", 1);
    setenv("FIRST_LEVEL_FILE", "../dati/FIRST_LEVEL_corto_GAT_2006_INV", 1);
    setenv("DIR_OUT_PP_BLOC", "testdata", 1);
    setenv("VPR_HEATING", "testdata/vpr_heat_GAT", 1);
    unlink("testdata/vpr_heat_GAT");
    setenv("FILE_T", "testdata/temperature.txt", 1);

    CUM_BAC* cb = new CUM_BAC;
    cb->do_quality = true;
    cb->do_beamblocking = true;
    cb->do_declutter = true;
    cb->do_class = false;
    cb->do_bloccorr = false;
    cb->do_vpr = true;
    cb->read_sp20_volume(fname, "GAT", 0);
    cb->setup_elaborazione(fname, "GAT");

    int ier = cb->elabora_dato();
    wassert(actual(ier) == 0);

    cb->caratterizzo_volume();

    // la combina_profili restituisce 1 se non riesce a costruire un profilo
    // perchè non piove o piove poco
    ier = cb->combina_profili("GAT");
    wassert(actual(ier) == 1);

    cb->heating = cb->profile_heating();
    wassert(actual(cb->heating) == 0);

    ier = cb->corr_vpr("GAT");
    wassert(actual(ier) == 1);

    // TODO: cb->stampa_vpr()

    cb->creo_cart();
    cb->creo_cart_z_lowris();

    // TODO: scrivo_out_file_bin

    delete cb;
}

template<> template<>
void to::test<8>()
{
    // versione BB_VPR_CLASS che corrisponde al parametro algo_corto_dev
    static const char* fname = "testdata/DBP2_070120141530_GATTATICO";
    setenv("FIRST_LEVEL_DIM_FILE", "../dati/FL_2006.DIM", 1);
    setenv("FIRST_LEVEL_FILE", "../dati/FIRST_LEVEL_corto_GAT_2006_INV", 1);
    setenv("DIR_OUT_PP_BLOC", "testdata", 1);
    setenv("VPR_HEATING", "testdata/vpr_heat_GAT", 1);
    unlink("testdata/vpr_heat_GAT");
    setenv("FILE_T", "testdata/temperature.txt", 1);

    CUM_BAC* cb = new CUM_BAC;
    cb->do_quality = true;
    cb->do_beamblocking = true;
    cb->do_declutter = true;
    cb->do_class = true;
    cb->do_bloccorr = false;
    cb->do_vpr = true;
    cb->read_sp20_volume(fname, "GAT", 0);
    cb->setup_elaborazione(fname, "GAT");

    int ier = cb->elabora_dato();
    wassert(actual(ier) == 0);

    cb->caratterizzo_volume();

    cb->classifica_rain();

    // la combina_profili restituisce 1 se non riesce a costruire un profilo
    // perchè non piove o piove poco
    ier = cb->combina_profili("GAT");
    wassert(actual(ier) == 1);

    cb->heating = cb->profile_heating();
    wassert(actual(cb->heating) == 0);

    ier = cb->corr_vpr("GAT");
    wassert(actual(ier) == 1); // TODO: cosa deve dare?

    // TODO: cb->stampa_vpr()

    cb->class_conv_fixme_find_a_name();

    cb->creo_cart();
    cb->creo_cart_z_lowris();

    // TODO: scrivo_out_file_bin

    delete cb;
}

template<> template<>
void to::test<9>()
{
    // versione ' m' che corrisponde al parametro algo_medio FIXME: serve?
    static const char* fname = "testdata/DBP2_070120141530_GATTATICO";
    setenv("FIRST_LEVEL_DIM_FILE", "../dati/FL_2006.DIM", 1);
    setenv("FIRST_LEVEL_FILE", "../dati/FIRST_LEVEL_corto_GAT_2006_INV", 1);
    setenv("DIR_OUT_PP_BLOC", "testdata", 1);

    CUM_BAC* cb = new CUM_BAC;
    cb->do_quality = false;
    cb->do_beamblocking = false;
    cb->do_declutter = false;
    cb->do_class = true;
    cb->do_bloccorr = false;
    cb->do_vpr = true;
    cb->read_sp20_volume(fname, "GAT", 0);
    cb->setup_elaborazione(fname, "GAT");

    int ier = cb->elabora_dato();
    wassert(actual(ier) == 0);

    // FIXME: anche con quality = false?
    //cb->caratterizzo_volume();

    cb->creo_cart();
    cb->creo_cart_z_lowris();

    // TODO: scrivo_out_file_bin

    delete cb;
}

template<> template<>
void to::test<10>()
{
    // versione BB che corrisponde al parametro algo_corto
    static const char* fname = "testdata/DBP2_060220140140_GATTATICO";
    setenv("FIRST_LEVEL_DIM_FILE", "../dati/FL_2006.DIM", 1);
    setenv("FIRST_LEVEL_FILE", "../dati/FIRST_LEVEL_corto_GAT_2006_INV", 1);
    setenv("DIR_OUT_PP_BLOC", "testdata", 1);

    CUM_BAC* cb = new CUM_BAC;
    cb->do_quality = true;
    cb->do_beamblocking = true;
    cb->do_declutter = true;
    cb->do_class = false;
    cb->do_bloccorr = false;
    cb->do_vpr = false;
    cb->read_sp20_volume(fname, "GAT", 0);
    cb->setup_elaborazione(fname, "GAT");

    int ier = cb->elabora_dato();
    wassert(actual(ier) == 0);

    cb->caratterizzo_volume();
    cb->creo_cart();
    cb->creo_cart_z_lowris();

    // TODO: scrivo_out_file_bin

    delete cb;
}

template<> template<>
void to::test<11>()
{
    // versione BB_VPR che corrisponde al parametro algo_corto_dev
    static const char* fname = "testdata/DBP2_060220140140_GATTATICO";
    setenv("FIRST_LEVEL_DIM_FILE", "../dati/FL_2006.DIM", 1);
    setenv("FIRST_LEVEL_FILE", "../dati/FIRST_LEVEL_corto_GAT_2006_INV", 1);
    setenv("DIR_OUT_PP_BLOC", "testdata", 1);
    setenv("VPR_HEATING", "testdata/vpr_heat_GAT", 1);
    unlink("testdata/vpr_heat_GAT");
    setenv("FILE_T", "testdata/temperature.txt", 1);

    CUM_BAC* cb = new CUM_BAC;
    cb->do_quality = true;
    cb->do_beamblocking = true;
    cb->do_declutter = true;
    cb->do_class = false;
    cb->do_bloccorr = false;
    cb->do_vpr = true;
    cb->read_sp20_volume(fname, "GAT", 0);
    cb->setup_elaborazione(fname, "GAT");

    int ier = cb->elabora_dato();
    wassert(actual(ier) == 0);

    cb->caratterizzo_volume();

    // la combina_profili restituisce 1 se non riesce a costruire un profilo
    // perchè non piove o piove poco
    ier = cb->combina_profili("GAT");
    wassert(actual(ier) == 1);

    cb->heating = cb->profile_heating();
    wassert(actual(cb->heating) == 0);

    ier = cb->corr_vpr("GAT");
    wassert(actual(ier) == 1);

    // TODO: cb->stampa_vpr()

    cb->creo_cart();
    cb->creo_cart_z_lowris();

    // TODO: scrivo_out_file_bin

    delete cb;
}

template<> template<>
void to::test<12>()
{
    // versione BB_VPR_CLASS che corrisponde al parametro algo_corto_dev
    static const char* fname = "testdata/DBP2_060220140140_GATTATICO";
    setenv("FIRST_LEVEL_DIM_FILE", "../dati/FL_2006.DIM", 1);
    setenv("FIRST_LEVEL_FILE", "../dati/FIRST_LEVEL_corto_GAT_2006_INV", 1);
    setenv("DIR_OUT_PP_BLOC", "testdata", 1);
    setenv("VPR_HEATING", "testdata/vpr_heat_GAT", 1);
    unlink("testdata/vpr_heat_GAT");
    setenv("FILE_T", "testdata/temperature.txt", 1);

    CUM_BAC* cb = new CUM_BAC;
    cb->do_quality = true;
    cb->do_beamblocking = true;
    cb->do_declutter = true;
    cb->do_class = true;
    cb->do_bloccorr = false;
    cb->do_vpr = true;
    cb->read_sp20_volume(fname, "GAT", 0);
    cb->setup_elaborazione(fname, "GAT");

    int ier = cb->elabora_dato();
    wassert(actual(ier) == 0);

    cb->caratterizzo_volume();

    cb->classifica_rain();

    // la combina_profili restituisce 1 se non riesce a costruire un profilo
    // perchè non piove o piove poco
    ier = cb->combina_profili("GAT");
    wassert(actual(ier) == 1);

    cb->heating = cb->profile_heating();
    wassert(actual(cb->heating) == 0);

    ier = cb->corr_vpr("GAT");
    wassert(actual(ier) == 1); // TODO: cosa deve dare?

    // TODO: cb->stampa_vpr()

    cb->class_conv_fixme_find_a_name();

    cb->creo_cart();
    cb->creo_cart_z_lowris();

    // TODO: scrivo_out_file_bin

    delete cb;
}

#if 0
template<> template<>
void to::test<6>()
{
    // Test di tutto quello che può essere chiamato
    static const char* fname = "testdata/DBP2_070120141530_GATTATICO";

    setenv("FIRST_LEVEL_DIM_FILE", "../dati/FL_2006.DIM", 1);
    setenv("FIRST_LEVEL_FILE", "../dati/FIRST_LEVEL_corto_GAT_2006_INV", 1);
    setenv("DIR_OUT_PP_BLOC", "testdata", 1);
    setenv("LAST_VPR", "testdata/last_vpr_GAT", 1);
    setenv("VPR0_FILE", "testdata/vpr_SPC", 1);
    setenv("VPR_HEATING", "vpr_heat_GAT", 1);

    CUM_BAC* cb = new CUM_BAC;
    cb->do_quality = true;
    cb->do_beamblocking = true;
    cb->do_declutter = false;
    cb->do_bloccorr = true;
    cb->do_vpr = true;
    cb->read_sp20_volume(fname, "GAT", 0);
    cb->setup_elaborazione(fname, "GAT");

    int ier = cb->elabora_dato();
    wassert(actual(ier) == 0);

    cb->caratterizzo_volume();

    ArrayStats<unsigned char> qual_stats;
    qual_stats.fill3(cb->qual);
    wassert(actual((unsigned)qual_stats.min) == 0);
    wassert(actual((unsigned)qual_stats.max) == 99);
    wassert(actual((unsigned)(qual_stats.avg * 100)) == 5874);

    ArrayStats<unsigned char> vpr_stats;
    vpr_stats.fill3(cb->flag_vpr);
    wassert(actual((unsigned)vpr_stats.min) == 0);
    wassert(actual((unsigned)vpr_stats.max) == 1);
    wassert(actual((unsigned)(vpr_stats.avg * 100)) == 91);

    ArrayStats<unsigned char> top_stats;
    top_stats.fill2(cb->top);
    wassert(actual((unsigned)top_stats.min) == 0);
    wassert(actual((unsigned)top_stats.max) == 35);
    wassert(actual((unsigned)(top_stats.avg * 100)) == 3);

    cb->classifica_rain();

    ier = cb->combina_profili("GAT");
    wassert(actual(ier) == 1);

    cb->heating = cb->profile_heating();
    wassert(actual(cb->heating) == 0);

    delete cb;
}
#endif

}
