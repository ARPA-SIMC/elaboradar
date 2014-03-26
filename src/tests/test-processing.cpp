#include <wibble/tests.h>
#include "cum_bac.h"
#include "logging.h"
#include <cstdio>
#include <cstdlib>
#include <vector>
#include "setwork.h"
#include <unistd.h>

using namespace wibble::tests;
using namespace cumbac;

namespace tut {

struct process_shar {
    Logging logging;
};
TESTGRP(process);

namespace {

template<typename T>
struct ArrayStats : public cumbac::ArrayStats<T>
{
    void fill(const PolarMap<T>& arr)
    {
        for (int i = 0; i < arr.beam_count * arr.beam_size; ++i)
            this->count_sample(arr.data[i], arr.beam_count * arr.beam_size);
    }

    template<int A, int B>
    void fill2(const T (&arr)[A][B])
    {
        for (int i = 0; i < A; ++i)
            for (int j = 0; j < B; ++j)
                this->count_sample(arr[i][j], A * B);
    }

    template<int A, int B, int C>
    void fill3(const T (&arr)[A][B][C])
    {
        for (int i = 0; i < A; ++i)
            for (int j = 0; j < B; ++j)
                for (int k = 0; k < C; ++k)
                    this->count_sample(arr[i][j][k], A * B * C);
    }

    void print()
    {
        fprintf(stderr, "min %f max %f avg %f, zeros: %u, ones: %u\n",
                (double)this->min, (double)this->max, this->avg,
                this->count_zeros, this->count_ones);
    }
};

}

template<> template<>
void to::test<1>()
{
    // Test elabora_dato, con tutti i do_* a false
    //printwork();
    unsetwork();
    static const char* fname = "testdata/DBP2_070120141530_GATTATICO";


    setenv("FIRST_LEVEL_FILE", "../dati/FIRST_LEVEL_corto_GAT_2006_INV", 1);
    setenv("DIR_OUT_PP_BLOC", "testdata", 1);
    setenv("FILE_T", "testdata/temperature.txt", 1);

    CUM_BAC* cb = new CUM_BAC("GAT");
    cb->do_clean= true;
    cb->do_readStaticMap=true;
    cb->read_sp20_volume(fname, 0);
    cb->setup_elaborazione(fname);

    wassert(actual(cb->calcolo_vpr->t_ground) == NODATAVPR);

    cb->elabora_dato();

    // Check results
    VolumeStats stats;
    cb->volume.compute_stats(stats);
    wassert(actual(stats.count_zeros[0]) == 0);
    wassert(actual(stats.count_zeros[1]) == 0);
    wassert(actual(stats.count_zeros[2]) == 0);
    wassert(actual(stats.count_zeros[3]) == 0);
    wassert(actual(stats.count_zeros[4]) == 0);
    wassert(actual(stats.count_ones[0]) == 197154);
    wassert(actual(stats.count_ones[1]) == 197495);
    wassert(actual(stats.count_ones[2]) == 197420);
    wassert(actual(stats.count_ones[3]) == 197522);
    wassert(actual(stats.count_ones[4]) == 197474);
    wassert(actual(stats.count_others[0]) == 446);
    wassert(actual(stats.count_others[1]) == 105);
    wassert(actual(stats.count_others[2]) == 180);
    wassert(actual(stats.count_others[3]) ==  78);
    wassert(actual(stats.count_others[4]) == 126);
    wassert(actual(stats.sum_others[0]) == 36449);
    wassert(actual(stats.sum_others[1]) ==  8452);
    wassert(actual(stats.sum_others[2]) == 16605);
    wassert(actual(stats.sum_others[3]) ==  4291);
    wassert(actual(stats.sum_others[4]) ==  7079);

    cb->caratterizzo_volume();

    ArrayStats<unsigned char> qual_stats;
    cb->qual->fill_array_stats(qual_stats);
    wassert(actual((unsigned)qual_stats.first).isfalse());
    wassert(actual((unsigned)qual_stats.min) == 1);
    wassert(actual((unsigned)qual_stats.max) == 99);
    wassert(actual((unsigned)(qual_stats.avg * 100)) == 1172);


    ArrayStats<unsigned char> vpr_stats;
    cb->calcolo_vpr->flag_vpr->fill_array_stats(vpr_stats);
    wassert(actual((unsigned)vpr_stats.first).isfalse());
    wassert(actual((unsigned)vpr_stats.min) == 0);
    wassert(actual((unsigned)vpr_stats.max) == 0);
    wassert(actual((unsigned)(vpr_stats.avg * 100)) == 0);

    ArrayStats<unsigned char> top_stats;
    top_stats.fill(cb->top);
    wassert(actual((unsigned)top_stats.first).isfalse());
    wassert(actual((unsigned)top_stats.min) == 0);
    wassert(actual((unsigned)top_stats.max) == 15);
    wassert(actual((unsigned)(top_stats.avg * 100)) == 0);

    delete cb;
}

template<> template<>
void to::test<2>()
{
    // Test elabora_dato, con tutti i do_* a true
    static const char* fname = "testdata/DBP2_070120141530_GATTATICO";

    unsetwork();
    setenv("FIRST_LEVEL_FILE", "../dati/FIRST_LEVEL_corto_GAT_2006_INV", 1);
    setenv("DIR_OUT_PP_BLOC", "testdata", 1);
    setenv("FILE_T", "testdata/temperature.txt", 1);
    setenv("LAST_VPR","testdata/last_vpr",1);
    setenv("FILE_ZERO_TERMICO","testdata/zero_termico.txt",1);

    CUM_BAC* cb = new CUM_BAC("GAT");
    cb->do_quality = true;
    cb->do_beamblocking = true;
    cb->do_declutter = false;
    cb->do_bloccorr = true;
    cb->do_vpr = true;
    cb->do_class = true;
    cb->do_clean= true;
    cb->do_readStaticMap=true;
    cb->read_sp20_volume(fname, 0);
    cb->setup_elaborazione(fname);

    wassert(actual((int)(cb->calcolo_vpr->t_ground * 100)) == 1010);

    cb->elabora_dato();

    // Check results
    VolumeStats stats;
    cb->volume.compute_stats(stats);
    wassert(actual(stats.count_zeros[0]) == 0);
    wassert(actual(stats.count_zeros[1]) == 0);
    wassert(actual(stats.count_zeros[2]) == 0);
    wassert(actual(stats.count_zeros[3]) == 0);
    wassert(actual(stats.count_zeros[4]) == 0);
    wassert(actual(stats.count_zeros[5]) == 0);
    wassert(actual(stats.count_ones[0]) == 197194);
    wassert(actual(stats.count_ones[1]) == 197495);
    wassert(actual(stats.count_ones[2]) == 197420);
    wassert(actual(stats.count_ones[3]) == 197522);
    wassert(actual(stats.count_ones[4]) == 197474);
    wassert(actual(stats.count_ones[5]) == 196872);
    wassert(actual(stats.count_others[0]) == 406);
    wassert(actual(stats.count_others[1]) == 105);
    wassert(actual(stats.count_others[2]) == 180);
    wassert(actual(stats.count_others[3]) ==  78);
    wassert(actual(stats.count_others[4]) == 126);
    wassert(actual(stats.count_others[5]) == 728);
    wassert(actual(stats.sum_others[0]) == 33757);
    wassert(actual(stats.sum_others[1]) ==  8452);
    wassert(actual(stats.sum_others[2]) == 16605);
    wassert(actual(stats.sum_others[3]) ==  4291);
    wassert(actual(stats.sum_others[4]) ==  7079);
    wassert(actual(stats.sum_others[5]) == 55505);

    cb->caratterizzo_volume();

    ArrayStats<unsigned char> qual_stats;
    cb->qual->fill_array_stats(qual_stats);
    wassert(actual((unsigned)qual_stats.first).isfalse());
    wassert(actual((unsigned)qual_stats.min) == 1);
    wassert(actual((unsigned)qual_stats.max) == 99);
    wassert(actual((unsigned)(qual_stats.avg * 100)) == 5506);

    ArrayStats<unsigned char> vpr_stats;
    cb->calcolo_vpr->flag_vpr->fill_array_stats(vpr_stats);
    wassert(actual((unsigned)vpr_stats.first).isfalse());
    wassert(actual((unsigned)vpr_stats.min) == 0);
    wassert(actual((unsigned)vpr_stats.max) == 1);
    wassert(actual((unsigned)(vpr_stats.avg * 100)) == 87);

    ArrayStats<unsigned char> top_stats;
    top_stats.fill(cb->top);
    wassert(actual((unsigned)top_stats.first).isfalse());
    wassert(actual((unsigned)top_stats.min) == 0);
    wassert(actual((unsigned)top_stats.max) == 15);
    wassert(actual((unsigned)(top_stats.avg * 100)) == 0);

    cb->calcolo_vpr->classifica_rain();

    ArrayStats<unsigned char> stratiform_stats;
    stratiform_stats.fill2(cb->calcolo_vpr->stratiform);
    wassert(actual((unsigned)stratiform_stats.first).isfalse());
    wassert(actual((unsigned)stratiform_stats.min) == 0);
    wassert(actual((unsigned)stratiform_stats.max) == 1);
    wassert(actual((unsigned)(stratiform_stats.avg * 100)) == 0);

    // calcolo_vpr->esegui_tutto();
    // conversione_convettiva();
    // creo_cart();
    // creo_cart_z_lowris();

    delete cb;
}

template<> template<>
void to::test<3>()
{
    // Test elabora_dato, con tutti i do_* a true
    static const char* fname = "testdata/DBP2_070120141530_GATTATICO";

    unsetwork();
    setenv("FIRST_LEVEL_FILE", "../dati/FIRST_LEVEL_corto_GAT_2006_INV", 1);
    setenv("DIR_OUT_PP_BLOC", "testdata", 1);

    CUM_BAC* cb = new CUM_BAC("GAT");
    cb->do_quality = true;
    cb->do_beamblocking = true;
    cb->do_declutter = false;
    cb->do_bloccorr = true;
    cb->do_vpr = true;
    cb->do_clean= true;
    cb->do_readStaticMap=true;
    cb->read_sp20_volume(fname, 0);
    cb->setup_elaborazione(fname);

    cb->elabora_dato();

    // Check results
    VolumeStats stats;
    cb->volume.compute_stats(stats);
    //stats.print(stdout);

    wassert(actual(stats.count_zeros[0]) == 0);
    wassert(actual(stats.count_zeros[1]) == 0);
    wassert(actual(stats.count_zeros[2]) == 0);
    wassert(actual(stats.count_zeros[3]) == 0);
    wassert(actual(stats.count_zeros[4]) == 0);
    wassert(actual(stats.count_ones[0]) == 197194);
    wassert(actual(stats.count_ones[1]) == 197495);
    wassert(actual(stats.count_ones[2]) == 197420);
    wassert(actual(stats.count_ones[3]) == 197522);
    wassert(actual(stats.count_ones[4]) == 197474);
    wassert(actual(stats.count_others[0]) == 406);
    wassert(actual(stats.count_others[1]) == 105);
    wassert(actual(stats.count_others[2]) == 180);
    wassert(actual(stats.count_others[3]) ==  78);
    wassert(actual(stats.count_others[4]) == 126);
    wassert(actual(stats.sum_others[0]) == 33757);
    wassert(actual(stats.sum_others[1]) ==  8452);
    wassert(actual(stats.sum_others[2]) == 16605);
    wassert(actual(stats.sum_others[3]) ==  4291);
    wassert(actual(stats.sum_others[4]) ==  7079);

    delete cb;
}

template<> template<>
void to::test<4>()
{
    // versione BB che corrisponde al parametro algo_corto
    static const char* fname = "testdata/DBP2_070120141530_GATTATICO";
    unsetwork();
    setenv("FIRST_LEVEL_FILE", "../dati/FIRST_LEVEL_corto_GAT_2006_INV", 1);
    setenv("DIR_OUT_PP_BLOC", "testdata", 1);
    setenv("VPR_HEATING", "testdata/vpr_heat_GAT", 1);
    unlink("testdata/vpr_heat_GAT");
    setenv("FILE_T", "testdata/temperature.txt", 1);

    CUM_BAC* cb = new CUM_BAC("GAT");
    cb->do_quality = true;
    cb->do_beamblocking = true;
    cb->do_declutter = true;
    cb->do_class = false;
    cb->do_bloccorr = false;
    cb->do_vpr = false;
    cb->do_clean= true;
    cb->do_readStaticMap=true;
    cb->read_sp20_volume(fname, 0);
    cb->setup_elaborazione(fname);

    cb->elabora_dato();

    cb->caratterizzo_volume();
    // in hray_inf
    // out hray
    // in dato_corrotto
    // in beam_blocking
    // in elev_fin
    // out qual
    // out flag_vpr
    // out top
    ArrayStats<unsigned char> stats_qual;
    cb->qual->fill_array_stats(stats_qual);
    wassert(actual((unsigned)stats_qual.first).isfalse());
    //wassert(actual((unsigned)stats_qual.count_zeros) == 1886400);
    wassert(actual((unsigned)stats_qual.count_zeros) == 0);
    wassert(actual((unsigned)stats_qual.count_ones) == 141553);
    wassert(actual((unsigned)stats_qual.min) == 1);
    wassert(actual((unsigned)stats_qual.max) == 99);
    wassert(actual((unsigned)(stats_qual.avg * 100)) == 5508);
    ArrayStats<unsigned char> stats_flag_vpr;
    cb->calcolo_vpr->flag_vpr->fill_array_stats(stats_flag_vpr);
    wassert(actual((unsigned)stats_flag_vpr.first).isfalse());
    wassert(actual((unsigned)stats_flag_vpr.count_zeros) == 1185600);
    wassert(actual((unsigned)stats_flag_vpr.count_ones) == 0);
    wassert(actual((unsigned)(stats_flag_vpr.avg * 100)) == 0);
    ArrayStats<unsigned char> stats_top;
    stats_top.fill(cb->top);
    wassert(actual((unsigned)stats_top.first).isfalse());
    wassert(actual((unsigned)stats_top.count_zeros) == 204764);
    wassert(actual((unsigned)stats_top.count_ones) == 0);
    wassert(actual((unsigned)stats_top.min) == 0);
    wassert(actual((unsigned)stats_top.max) == 15);
    wassert(actual((unsigned)(stats_top.avg * 100)) == 0);

    cb->creo_cart();
    wassert(actual((unsigned)cb->cart.min()) == 0);
    wassert(actual((unsigned)cb->cart.avg()) == 1);
    wassert(actual((unsigned)cb->cart.max()) == 227);
    wassert(actual(cb->cartm.min()) == 0);
    wassert(actual(cb->cartm.max()) == 0);
    wassert(actual((unsigned)cb->topxy.min()) == 0);
    wassert(actual((unsigned)cb->topxy.avg()) == 0);
    wassert(actual((unsigned)cb->topxy.max()) == 15);
    wassert(actual((unsigned)cb->qual_Z_cart.min()) == 0);
    wassert(actual((unsigned)cb->qual_Z_cart.avg()) == 30);
    wassert(actual((unsigned)cb->qual_Z_cart.max()) == 98);
    wassert(actual((unsigned)cb->quota_cart.min()) == 0);
    wassert(actual((unsigned)cb->quota_cart.max()) == 0);
    wassert(actual((unsigned)cb->dato_corr_xy.min()) == 0);
    wassert(actual((unsigned)cb->dato_corr_xy.max()) == 0);
    wassert(actual((unsigned)cb->beam_blocking_xy.min()) == 0);
    wassert(actual((unsigned)cb->beam_blocking_xy.avg()) == 9);
    wassert(actual((unsigned)cb->beam_blocking_xy.max()) == 51);
    wassert(actual((unsigned)cb->elev_fin_xy.min()) == 0);
    wassert(actual((unsigned)cb->elev_fin_xy.avg()) == 0);
    wassert(actual((unsigned)cb->elev_fin_xy.max()) == 3);
    wassert(actual((unsigned)cb->neve_cart.min()) == 0);
    wassert(actual((unsigned)cb->neve_cart.max()) == 0);
    wassert(actual((unsigned)cb->corr_cart.min()) == 0);
    wassert(actual((unsigned)cb->corr_cart.max()) == 0);
    wassert(actual((unsigned)cb->conv_cart.min()) == 0);
    wassert(actual((unsigned)cb->conv_cart.max()) == 0);

    cb->creo_cart_z_lowris();
    wassert(actual((unsigned)cb->z_out.min()) == 0);
    wassert(actual((unsigned)cb->z_out.avg()) == 1);
    wassert(actual((unsigned)cb->z_out.max()) == 227);
    wassert(actual((unsigned)cb->qual_Z_1x1.min()) == 0);
    wassert(actual((unsigned)cb->qual_Z_1x1.avg()) == 30);
    wassert(actual((unsigned)cb->qual_Z_1x1.max()) == 98);
    wassert(actual((unsigned)cb->quota_1x1.min()) == 128);
    wassert(actual((unsigned)cb->quota_1x1.max()) == 128);
    wassert(actual((unsigned)cb->dato_corr_1x1.min()) == 0);
    wassert(actual((unsigned)cb->dato_corr_1x1.max()) == 0);
    wassert(actual((unsigned)cb->elev_fin_1x1.min()) == 0);
    wassert(actual((unsigned)cb->elev_fin_1x1.avg()) == 0);
    wassert(actual((unsigned)cb->elev_fin_1x1.max()) == 3);
    wassert(actual((unsigned)cb->beam_blocking_1x1.min()) == 0);
    wassert(actual((unsigned)cb->beam_blocking_1x1.avg()) == 8);
    wassert(actual((unsigned)cb->beam_blocking_1x1.max()) == 51);
    wassert(actual((unsigned)cb->top_1x1.min()) == 0);
    wassert(actual((unsigned)cb->top_1x1.avg()) == 0);
    wassert(actual((unsigned)cb->top_1x1.max()) == 15);
    wassert(actual((unsigned)cb->neve_1x1.min()) == 0);
    wassert(actual((unsigned)cb->neve_1x1.max()) == 0);
    wassert(actual((unsigned)cb->corr_1x1.min()) == 0);
    wassert(actual((unsigned)cb->corr_1x1.max()) == 0);
    wassert(actual((unsigned)cb->conv_1x1.min()) == 0);
    wassert(actual((unsigned)cb->conv_1x1.max()) == 0);

    delete cb;
}

template<> template<>
void to::test<5>()
{
    // versione BB_VPR che corrisponde al parametro algo_corto_dev
	LOG_CATEGORY("Test");
	LOG_INFO ("Start test 7");

    static const char* fname = "testdata/DBP2_070120141530_GATTATICO";
    unsetwork();
    setenv("FIRST_LEVEL_FILE", "../dati/FIRST_LEVEL_corto_GAT_2006_INV", 1);
    setenv("DIR_OUT_PP_BLOC", "testdata", 1);
    setenv("VPR_HEATING", "testdata/vpr_heat_GAT", 1);
    setenv("VPR0_FILE", "testdata/ultimo_vpr", 1);
    unlink("testdata/vpr_heat_GAT");
    setenv("FILE_T", "testdata/temperature.txt", 1);
    setenv("LAST_VPR","testdata/last_vpr",1);
    printwork();
    CUM_BAC* cb = new CUM_BAC("GAT");
    cb->do_quality = true;
    cb->do_beamblocking = true;
    cb->do_declutter = true;
    cb->do_class = false;
    cb->do_bloccorr = false;
    cb->do_vpr = true;
    cb->do_readStaticMap=true;
    cb->read_sp20_volume(fname, 0);
    cb->setup_elaborazione(fname);
LOG_INFO ("Chiamo elabora_dato");
    cb->elabora_dato();
LOG_INFO("Chiamo caratterizzo volumi");

    cb->caratterizzo_volume();
    ArrayStats<unsigned char> stats_qual;
    cb->qual->fill_array_stats(stats_qual);
    wassert(actual((unsigned)stats_qual.first).isfalse());
    wassert(actual((unsigned)stats_qual.count_zeros) == 0);
    wassert(actual((unsigned)stats_qual.count_ones) == 141553);
    wassert(actual((unsigned)stats_qual.min) == 1);
    wassert(actual((unsigned)stats_qual.max) == 99);
    wassert(actual((unsigned)(stats_qual.avg * 100)) == 5508);
    ArrayStats<unsigned char> stats_flag_vpr;
    cb->calcolo_vpr->flag_vpr->fill_array_stats(stats_flag_vpr);
    wassert(actual((unsigned)stats_flag_vpr.first).isfalse());
    wassert(actual((unsigned)stats_flag_vpr.count_zeros) == 143320);
    wassert(actual((unsigned)stats_flag_vpr.count_ones) == 1042280);
    wassert(actual((unsigned)(stats_flag_vpr.avg * 100)) == 87);
    ArrayStats<unsigned char> stats_top;
    stats_top.fill(cb->top);
    wassert(actual((unsigned)stats_top.first).isfalse());
    wassert(actual((unsigned)stats_top.count_zeros) == 203211);
    wassert(actual((unsigned)stats_top.count_ones) == 4);
    wassert(actual((unsigned)stats_top.min) == 0);
    wassert(actual((unsigned)stats_top.max) == 36);
    wassert(actual((unsigned)(stats_top.avg * 100)) == 10);

    // la combina_profili restituisce 1 se non riesce a costruire un profilo
    // perchè non piove o piove poco
    int ier = cb->calcolo_vpr->combina_profili();
    wassert(actual(ier) == 1);

    cb->calcolo_vpr->heating = cb->calcolo_vpr->profile_heating();
    wassert(actual(cb->calcolo_vpr->heating) == 0);

    //ier = cb->calcolo_vpr->corr_vpr();
    //wassert(actual(ier) == 1);

    // TODO: cb->stampa_vpr()

    cb->creo_cart();
    wassert(actual((unsigned)cb->cart.min()) == 0);
    wassert(actual((unsigned)cb->cart.avg()) == 10);
    wassert(actual((unsigned)cb->cart.max()) == 227);
    wassert(actual(cb->cartm.min()) == 0);
    wassert(actual(cb->cartm.max()) == 0);
    wassert(actual((unsigned)cb->topxy.min()) == 0);
    wassert(actual((unsigned)cb->topxy.avg()) == 0);
    wassert(actual((unsigned)cb->topxy.max()) == 36);
    wassert(actual((unsigned)cb->qual_Z_cart.min()) == 0);
    wassert(actual((unsigned)cb->qual_Z_cart.avg()) == 30);
    wassert(actual((unsigned)cb->qual_Z_cart.max()) == 98);
    wassert(actual((unsigned)cb->quota_cart.min()) == 0);
    wassert(actual((unsigned)cb->quota_cart.max()) == 0);
    wassert(actual((unsigned)cb->dato_corr_xy.min()) == 0);
    wassert(actual((unsigned)cb->dato_corr_xy.max()) == 0);
    wassert(actual((unsigned)cb->beam_blocking_xy.min()) == 0);
    wassert(actual((unsigned)cb->beam_blocking_xy.avg()) == 9);
    wassert(actual((unsigned)cb->beam_blocking_xy.max()) == 51);
    wassert(actual((unsigned)cb->elev_fin_xy.min()) == 0);
    wassert(actual((unsigned)cb->elev_fin_xy.avg()) == 0);
    wassert(actual((unsigned)cb->elev_fin_xy.max()) == 3);
    wassert(actual((unsigned)cb->neve_cart.min()) == 0);
    wassert(actual((unsigned)cb->neve_cart.max()) == 1);
    wassert(actual((unsigned)cb->corr_cart.min()) == 0);
    wassert(actual((unsigned)cb->corr_cart.max()) == 0);
    wassert(actual((unsigned)cb->conv_cart.min()) == 0);
    wassert(actual((unsigned)cb->conv_cart.max()) == 0);

    cb->creo_cart_z_lowris();
    wassert(actual((unsigned)cb->z_out.min()) == 0);
    wassert(actual((unsigned)cb->z_out.avg()) == 16);
    wassert(actual((unsigned)cb->z_out.max()) == 227);
    wassert(actual((unsigned)cb->qual_Z_1x1.min()) == 0);
    wassert(actual((unsigned)cb->qual_Z_1x1.avg()) == 29);
    wassert(actual((unsigned)cb->qual_Z_1x1.max()) == 98);
    wassert(actual((unsigned)cb->quota_1x1.min()) == 128);
    wassert(actual((unsigned)cb->quota_1x1.max()) == 128);
    wassert(actual((unsigned)cb->dato_corr_1x1.min()) == 0);
    wassert(actual((unsigned)cb->dato_corr_1x1.max()) == 0);
    wassert(actual((unsigned)cb->elev_fin_1x1.min()) == 0);
    wassert(actual((unsigned)cb->elev_fin_1x1.avg()) == 0);
    wassert(actual((unsigned)cb->elev_fin_1x1.max()) == 3);
    wassert(actual((unsigned)cb->beam_blocking_1x1.min()) == 0);
    wassert(actual((unsigned)cb->beam_blocking_1x1.avg()) == 9);
    wassert(actual((unsigned)cb->beam_blocking_1x1.max()) == 51);
    wassert(actual((unsigned)cb->top_1x1.min()) == 0);
    wassert(actual((unsigned)cb->top_1x1.max()) == 36);
    wassert(actual((unsigned)cb->top_1x1.avg()) == 0);
    wassert(actual((unsigned)cb->neve_1x1.min()) == 0);
    wassert(actual((unsigned)cb->neve_1x1.max()) == 1);
    wassert(actual((unsigned)cb->neve_1x1.avg()) == 1);
    wassert(actual((unsigned)cb->corr_1x1.min()) == 0);
    wassert(actual((unsigned)cb->corr_1x1.max()) == 0);
    wassert(actual((unsigned)cb->conv_1x1.min()) == 0);
    wassert(actual((unsigned)cb->conv_1x1.max()) == 0);

    delete cb;
    LOG_INFO("End test 5");
}

template<> template<>
void to::test<6>()
{
    // versione BB_VPR_CLASS che corrisponde al parametro algo_corto_dev
	LOG_CATEGORY("Test");
	LOG_INFO ("Start test 8");

    static const char* fname = "testdata/DBP2_070120141530_GATTATICO";
    unsetwork();
    setenv("FIRST_LEVEL_FILE", "../dati/FIRST_LEVEL_corto_GAT_2006_INV", 1);
    setenv("DIR_OUT_PP_BLOC", "testdata", 1);
    setenv("VPR_HEATING", "testdata/vpr_heat_GAT", 1);
    unlink("testdata/vpr_heat_GAT");
    setenv("FILE_T", "testdata/temperature.txt", 1);
    setenv("LAST_VPR","testdata/last_vpr",1);
unlink("LAST_VPR");
    setenv("FILE_ZERO_TERMICO","testdata/zero_termico.txt",1);
    printwork();
    CUM_BAC* cb = new CUM_BAC("GAT");
    cb->do_quality = true;
    cb->do_beamblocking = true;
    cb->do_declutter = true;
    cb->do_class = true;
    cb->do_bloccorr = false;
    cb->do_vpr = true;
    cb->do_clean= true;
    cb->do_readStaticMap=true;
    cb->read_sp20_volume(fname, 0);
    cb->setup_elaborazione(fname);

    cb->elabora_dato();

    cb->caratterizzo_volume();
    ArrayStats<unsigned char> stats_qual;
    cb->qual->fill_array_stats(stats_qual);
    wassert(actual((unsigned)stats_qual.first).isfalse());
    wassert(actual((unsigned)stats_qual.count_zeros) == 0);
    wassert(actual((unsigned)stats_qual.count_ones) == 141553);
    wassert(actual((unsigned)stats_qual.min) == 1);
    wassert(actual((unsigned)stats_qual.max) == 99);
    wassert(actual((unsigned)(stats_qual.avg * 100)) == 5508);
    ArrayStats<unsigned char> stats_flag_vpr;
    cb->calcolo_vpr->flag_vpr->fill_array_stats(stats_flag_vpr);
    wassert(actual((unsigned)stats_flag_vpr.first).isfalse());
    wassert(actual((unsigned)stats_flag_vpr.count_zeros) == 143320);
    wassert(actual((unsigned)stats_flag_vpr.count_ones) == 1042280);
    wassert(actual((unsigned)(stats_flag_vpr.avg * 100)) == 87);
    ArrayStats<unsigned char> stats_top;
    stats_top.fill(cb->top);
    wassert(actual((unsigned)stats_top.first).isfalse());
    wassert(actual((unsigned)stats_top.count_zeros) == 204764);
    wassert(actual((unsigned)stats_top.count_ones) == 0);
    wassert(actual((unsigned)stats_top.min) == 0);
    wassert(actual((unsigned)stats_top.max) == 15);
    wassert(actual((unsigned)(stats_top.avg * 100)) == 0);

    cb->calcolo_vpr->classifica_rain();

    // la combina_profili restituisce 1 se non riesce a costruire un profilo
    // perchè non piove o piove poco
    int ier = cb->calcolo_vpr->combina_profili();
    wassert(actual(ier) == 1);

    cb->calcolo_vpr->heating = cb->calcolo_vpr->profile_heating();
    wassert(actual(cb->calcolo_vpr->heating) == 0);

    ier = cb->calcolo_vpr->corr_vpr();
    wassert(actual(ier) == 1); // TODO: cosa deve dare?

    // TODO: cb->stampa_vpr()

    cb->conversione_convettiva();

    cb->creo_cart();

LOG_INFO("cart         min avg max :%d %d %d",(unsigned)cb->cart.min(),(unsigned)cb->cart.avg(),(unsigned)cb->cart.max());
LOG_INFO("cartm        min avg max :%d %d %d",(unsigned)cb->cartm.min(),(unsigned)cb->cartm.avg(),(unsigned)cb->cartm.max());
LOG_INFO("topxy        min avg max :%d %d %d",(unsigned)cb->topxy.min(),(unsigned)cb->topxy.avg(),(unsigned)cb->topxy.max());
LOG_INFO("qual_Z_cart  min avg max :%d %d %d",(unsigned)cb->qual_Z_cart.min(),(unsigned)cb->qual_Z_cart.avg(),(unsigned)cb->qual_Z_cart.max());
LOG_INFO("quota_cart   min avg max :%d %d %d",(unsigned)cb->quota_cart.min(),(unsigned)cb->quota_cart.avg(),(unsigned)cb->quota_cart.max());
LOG_INFO("dato_corr_xy min avg max :%d %d %d",(unsigned)cb->dato_corr_xy.min(),(unsigned)cb->dato_corr_xy.avg(),(unsigned)cb->dato_corr_xy.max());
LOG_INFO("beam_blocking_xy min avg max :%d %d %d",(unsigned)cb->beam_blocking_xy.min(),(unsigned)cb->beam_blocking_xy.avg(),(unsigned)cb->beam_blocking_xy.max());
LOG_INFO("elev_fin_xy  min avg max :%d %d %d",(unsigned)cb->elev_fin_xy.min(),(unsigned)cb->elev_fin_xy.avg(),(unsigned)cb->elev_fin_xy.max());
LOG_INFO("neve_cart    min avg max :%d %d %d",(unsigned)cb->neve_cart.min(),(unsigned)cb->neve_cart.avg(),(unsigned)cb->neve_cart.max());
LOG_INFO("corr_cart    min avg max :%d %d %d",(unsigned)cb->corr_cart.min(),(unsigned)cb->corr_cart.avg(),(unsigned)cb->corr_cart.max());
LOG_INFO("conv_cart    min avg max :%d %d %d",(unsigned)cb->conv_cart.min(),(unsigned)cb->conv_cart.avg(),(unsigned)cb->conv_cart.max());
LOG_INFO("------------------------------");
    wassert(actual((unsigned)cb->cart.min()) == 0);
    wassert(actual((unsigned)cb->cart.avg()) == 1);
    wassert(actual((unsigned)cb->cart.max()) == 227);
    wassert(actual(cb->cartm.min()) == 0);
    wassert(actual(cb->cartm.max()) == 0);
    wassert(actual((unsigned)cb->topxy.min()) == 0);
    wassert(actual((unsigned)cb->topxy.avg()) == 0);
    wassert(actual((unsigned)cb->topxy.max()) == 15);
    wassert(actual((unsigned)cb->qual_Z_cart.min()) == 0);
    wassert(actual((unsigned)cb->qual_Z_cart.avg()) == 30);
    wassert(actual((unsigned)cb->qual_Z_cart.max()) == 98);
    wassert(actual((unsigned)cb->quota_cart.min()) == 0);
    wassert(actual((unsigned)cb->quota_cart.max()) == 0);
    wassert(actual((unsigned)cb->dato_corr_xy.min()) == 0);
    wassert(actual((unsigned)cb->dato_corr_xy.max()) == 0);
    wassert(actual((unsigned)cb->beam_blocking_xy.min()) == 0);
    wassert(actual((unsigned)cb->beam_blocking_xy.avg()) == 9);
    wassert(actual((unsigned)cb->beam_blocking_xy.max()) == 51);
    wassert(actual((unsigned)cb->elev_fin_xy.min()) == 0);
    wassert(actual((unsigned)cb->elev_fin_xy.avg()) == 0);
    wassert(actual((unsigned)cb->elev_fin_xy.max()) == 3);
    wassert(actual((unsigned)cb->neve_cart.min()) == 0);
    wassert(actual((unsigned)cb->neve_cart.max()) == 1);
    wassert(actual((unsigned)cb->corr_cart.min()) == 0);
    wassert(actual((unsigned)cb->corr_cart.max()) == 0);
    wassert(actual((unsigned)cb->conv_cart.min()) == 0);
    wassert(actual((unsigned)cb->conv_cart.max()) == 0);

    cb->creo_cart_z_lowris();
LOG_INFO("z_out         min avg max :%d %d %d",(unsigned)cb->z_out.min(),(unsigned)cb->z_out.avg(),(unsigned)cb->z_out.max());
LOG_INFO("qual_Z_1x1    min avg max :%d %d %d",(unsigned)cb->qual_Z_1x1.min(),(unsigned)cb->qual_Z_1x1.avg(),(unsigned)cb->qual_Z_1x1.max());
LOG_INFO("quota_1x1     min avg max :%d %d %d",(unsigned)cb->quota_1x1.min(),(unsigned)cb->quota_1x1.avg(),(unsigned)cb->quota_1x1.max());
LOG_INFO("dato_corr_1x1 min avg max :%d %d %d",(unsigned)cb->dato_corr_1x1.min(),(unsigned)cb->dato_corr_1x1.avg(),(unsigned)cb->dato_corr_1x1.max());
LOG_INFO("elev_fin_1x1  min avg max :%d %d %d",(unsigned)cb->elev_fin_1x1.min(),(unsigned)cb->elev_fin_1x1.avg(),(unsigned)cb->elev_fin_1x1.max());
LOG_INFO("beam_blocking_1x1 min avg max :%d %d %d",(unsigned)cb->beam_blocking_1x1.min(),(unsigned)cb->beam_blocking_1x1.avg(),(unsigned)cb->beam_blocking_1x1.max());
LOG_INFO("top_1x1       min avg max :%d %d %d",(unsigned)cb->top_1x1.min(),(unsigned)cb->top_1x1.avg(),(unsigned)cb->top_1x1.max());
LOG_INFO("neve_1x1      min avg max :%d %d %d",(unsigned)cb->neve_1x1.min(),(unsigned)cb->neve_1x1.avg(),(unsigned)cb->neve_1x1.max());
LOG_INFO("corr_1x1      min avg max :%d %d %d",(unsigned)cb->corr_1x1.min(),(unsigned)cb->corr_1x1.avg(),(unsigned)cb->corr_1x1.max());
LOG_INFO("conv_1x1      min avg max :%d %d %d",(unsigned)cb->conv_1x1.min(),(unsigned)cb->conv_1x1.avg(),(unsigned)cb->conv_1x1.max());
    wassert(actual((unsigned)cb->z_out.min()) == 0);
    wassert(actual((unsigned)cb->z_out.avg()) == 1);
    wassert(actual((unsigned)cb->z_out.max()) == 227);
    wassert(actual((unsigned)cb->qual_Z_1x1.min()) == 0);
    wassert(actual((unsigned)cb->qual_Z_1x1.avg()) == 30);
    wassert(actual((unsigned)cb->qual_Z_1x1.max()) == 98);
    wassert(actual((unsigned)cb->quota_1x1.min()) == 128);
    wassert(actual((unsigned)cb->quota_1x1.max()) == 128);
    wassert(actual((unsigned)cb->dato_corr_1x1.min()) == 0);
    wassert(actual((unsigned)cb->dato_corr_1x1.max()) == 0);
    wassert(actual((unsigned)cb->elev_fin_1x1.min()) == 0);
    wassert(actual((unsigned)cb->elev_fin_1x1.avg()) == 0);
    wassert(actual((unsigned)cb->elev_fin_1x1.max()) == 3);
    wassert(actual((unsigned)cb->beam_blocking_1x1.min()) == 0);
    wassert(actual((unsigned)cb->beam_blocking_1x1.avg()) == 8);
    wassert(actual((unsigned)cb->beam_blocking_1x1.max()) == 51);
    wassert(actual((unsigned)cb->top_1x1.min()) == 0);
    wassert(actual((unsigned)cb->top_1x1.avg()) == 0);
    wassert(actual((unsigned)cb->top_1x1.max()) == 15);
    wassert(actual((unsigned)cb->neve_1x1.min()) == 0);
    wassert(actual((unsigned)cb->neve_1x1.max()) == 1);
    wassert(actual((unsigned)cb->neve_1x1.avg()) == 1);
    wassert(actual((unsigned)cb->corr_1x1.min()) == 0);
    wassert(actual((unsigned)cb->corr_1x1.max()) == 0);
    wassert(actual((unsigned)cb->conv_1x1.min()) == 0);
    wassert(actual((unsigned)cb->conv_1x1.max()) == 0);

    // TODO: scrivo_out_file_bin

    delete cb;
}

template<> template<>
void to::test<7>()
{
	LOG_CATEGORY("Test");
	LOG_INFO ("Start test 9");

    // versione BB che corrisponde al parametro algo_corto
    static const char* fname = "testdata/DBP2_060220140140_GATTATICO";
    unsetwork();
    setenv("FIRST_LEVEL_FILE", "../dati/FIRST_LEVEL_corto_GAT_2006_INV", 1);
    setenv("DIR_OUT_PP_BLOC", "testdata", 1);
    setenv("VPR_HEATING", "testdata/vpr_heat_GAT", 1);
    setenv("FILE_ZERO_TERMICO", "testdata/20140206/0termico.prev", 1);
    unlink("testdata/vpr_heat_GAT");
    setenv("VPR0_FILE", "testdata/ultimo_vpr", 1);
    unlink("testdata/ultimo_vpr");
    setenv("FILE_T", "testdata/temperature.txt", 1);
	printwork();

    CUM_BAC* cb = new CUM_BAC("GAT");
    cb->do_quality = true;
    cb->do_beamblocking = true;
    cb->do_declutter = true;
    cb->do_class = false;
    cb->do_bloccorr = false;
    cb->do_vpr = false;
    cb->do_clean= true;
    cb->do_readStaticMap=true;
    cb->read_sp20_volume(fname, 0);
    cb->setup_elaborazione(fname);
    wassert(actual(cb->calcolo_vpr) != (void*)0);

    cb->elabora_dato();

    VolumeStats stats;
    cb->volume.compute_stats(stats);
    wassert(actual(stats.count_zeros[0]) == 0);
    wassert(actual(stats.count_zeros[1]) == 0);
    wassert(actual(stats.count_zeros[2]) == 0);
    wassert(actual(stats.count_zeros[3]) == 0);
    wassert(actual(stats.count_zeros[4]) == 0);
    wassert(actual(stats.count_ones[0]) ==  61406);
    wassert(actual(stats.count_ones[1]) ==  72916);
    wassert(actual(stats.count_ones[2]) ==  96948);
    wassert(actual(stats.count_ones[3]) == 110168);
    wassert(actual(stats.count_ones[4]) == 119120);
    wassert(actual(stats.count_others[0]) == 136194);
    wassert(actual(stats.count_others[1]) == 124684);
    wassert(actual(stats.count_others[2]) == 100652);
    wassert(actual(stats.count_others[3]) ==  87432);
    wassert(actual(stats.count_others[4]) ==  78480);
    wassert(actual(stats.sum_others[0]) == 18104214);
    wassert(actual(stats.sum_others[1]) == 16004410);
    wassert(actual(stats.sum_others[2]) == 12562044);
    wassert(actual(stats.sum_others[3]) == 10553980);
    wassert(actual(stats.sum_others[4]) ==  9141377);

    ArrayStats<unsigned char> beam_blocking_stats;
    beam_blocking_stats.fill(cb->beam_blocking);
    wassert(actual((unsigned)beam_blocking_stats.first).isfalse());
    wassert(actual((unsigned)beam_blocking_stats.min) == 0);
    wassert(actual((unsigned)beam_blocking_stats.max) == 51);
    wassert(actual((unsigned)(beam_blocking_stats.avg * 100)) == 1364);

    ArrayStats<int> stat_anap_stats;
    stat_anap_stats.fill2(cb->stat_anap);
    wassert(actual((unsigned)stat_anap_stats.first).isfalse());
    wassert(actual((unsigned)stat_anap_stats.min) == 0);
    wassert(actual((unsigned)stat_anap_stats.max) == 0);

    ArrayStats<int> stat_anap_tot_stats;
    stat_anap_tot_stats.fill2(cb->stat_anap_tot);
    wassert(actual((unsigned)stat_anap_tot_stats.first).isfalse());
    wassert(actual((unsigned)stat_anap_tot_stats.min) == 0);
    wassert(actual((unsigned)stat_anap_tot_stats.max) == 0);

    ArrayStats<long int> stat_bloc_stats;
    stat_bloc_stats.fill2(cb->stat_bloc);
    wassert(actual((unsigned)stat_bloc_stats.first).isfalse());
    wassert(actual((unsigned)stat_bloc_stats.min) == 0);
    wassert(actual((unsigned)stat_bloc_stats.max) == 0);

    ArrayStats<int> stat_elev_stats;
    stat_elev_stats.fill2(cb->stat_elev);
    wassert(actual((unsigned)stat_elev_stats.first).isfalse());
    wassert(actual((unsigned)stat_elev_stats.min) == 0);
    wassert(actual((unsigned)stat_elev_stats.max) == 0);

    ArrayStats<unsigned char> dato_corrotto_stats;
    dato_corrotto_stats.fill(cb->dato_corrotto);
    wassert(actual((unsigned)dato_corrotto_stats.first).isfalse());
    wassert(actual((unsigned)dato_corrotto_stats.min) == 0);
    wassert(actual((unsigned)dato_corrotto_stats.max) == 0);


    cb->caratterizzo_volume();
    wassert(actual(cb->calcolo_vpr) != (void*)0);

    ArrayStats<unsigned char> stats_qual;
    cb->qual->fill_array_stats(stats_qual);
    wassert(actual((unsigned)stats_qual.first).isfalse());
    wassert(actual((unsigned)stats_qual.count_zeros) == 0);
    wassert(actual((unsigned)stats_qual.count_ones) == 159263);
    wassert(actual((unsigned)stats_qual.min) == 1);
    wassert(actual((unsigned)stats_qual.max) == 99);
    wassert(actual((unsigned)(stats_qual.avg * 100)) == 5762);

    ArrayStats<unsigned char> stats_flag_vpr;
    cb->calcolo_vpr->flag_vpr->fill_array_stats(stats_flag_vpr);
    wassert(actual((unsigned)stats_flag_vpr.first).isfalse());
    wassert(actual((unsigned)stats_flag_vpr.count_zeros) == 2173600);
    wassert(actual((unsigned)stats_flag_vpr.count_ones) == 0);
    wassert(actual((unsigned)(stats_flag_vpr.avg * 100)) == 0);

    ArrayStats<unsigned char> stats_top;
    stats_top.fill(cb->top);
    wassert(actual((unsigned)stats_top.first).isfalse());
    wassert(actual((unsigned)stats_top.count_zeros) == 112982);
    wassert(actual((unsigned)stats_top.count_ones) == 159);
    wassert(actual((unsigned)stats_top.min) == 0);
    wassert(actual((unsigned)stats_top.max) == 76);
    wassert(actual((unsigned)(stats_top.avg * 100)) == 552);

    ArrayStats<float> stats_hray;
    stats_hray.fill2(cb->hray);
    wassert(actual((unsigned)stats_hray.first).isfalse());
    wassert(actual((unsigned)stats_hray.count_zeros) == 2138);
    wassert(actual((unsigned)stats_hray.min) == 0);
    wassert(actual((unsigned)stats_hray.max) == 52796);
    wassert(actual((unsigned)(stats_hray.avg * 100)) == 637556);


    cb->creo_cart();

LOG_INFO("cart         min avg max :%d %d %d",(unsigned)cb->cart.min(),(unsigned)cb->cart.avg(),(unsigned)cb->cart.max());
LOG_INFO("cartm        min avg max :%d %d %d",(unsigned)cb->cartm.min(),(unsigned)cb->cartm.avg(),(unsigned)cb->cartm.max());
LOG_INFO("topxy        min avg max :%d %d %d",(unsigned)cb->topxy.min(),(unsigned)cb->topxy.avg(),(unsigned)cb->topxy.max());
LOG_INFO("qual_Z_cart  min avg max :%d %d %d",(unsigned)cb->qual_Z_cart.min(),(unsigned)cb->qual_Z_cart.avg(),(unsigned)cb->qual_Z_cart.max());
LOG_INFO("quota_cart   min avg max :%d %d %d",(unsigned)cb->quota_cart.min(),(unsigned)cb->quota_cart.avg(),(unsigned)cb->quota_cart.max());
LOG_INFO("dato_corr_xy min avg max :%d %d %d",(unsigned)cb->dato_corr_xy.min(),(unsigned)cb->dato_corr_xy.avg(),(unsigned)cb->dato_corr_xy.max());
LOG_INFO("beam_blocking_xy min avg max :%d %d %d",(unsigned)cb->beam_blocking_xy.min(),(unsigned)cb->beam_blocking_xy.avg(),(unsigned)cb->beam_blocking_xy.max());
LOG_INFO("elev_fin_xy  min avg max :%d %d %d",(unsigned)cb->elev_fin_xy.min(),(unsigned)cb->elev_fin_xy.avg(),(unsigned)cb->elev_fin_xy.max());
LOG_INFO("neve_cart    min avg max :%d %d %d",(unsigned)cb->neve_cart.min(),(unsigned)cb->neve_cart.avg(),(unsigned)cb->neve_cart.max());
LOG_INFO("corr_cart    min avg max :%d %d %d",(unsigned)cb->corr_cart.min(),(unsigned)cb->corr_cart.avg(),(unsigned)cb->corr_cart.max());
LOG_INFO("conv_cart    min avg max :%d %d %d",(unsigned)cb->conv_cart.min(),(unsigned)cb->conv_cart.avg(),(unsigned)cb->conv_cart.max());
    wassert(actual((unsigned)(unsigned)cb->cart.min()) == 0); 
    wassert(actual((unsigned)cb->cart.avg()) == 53);
    wassert(actual((unsigned)cb->cart.max()) == 255);
    wassert(actual(cb->cartm.min()) == 0);
    wassert(actual(cb->cartm.max()) == 0);
    wassert(actual((unsigned)cb->topxy.min()) == 0);
    wassert(actual((unsigned)cb->topxy.avg()) == 3);
    wassert(actual((unsigned)cb->topxy.max()) == 76);
    wassert(actual((unsigned)cb->qual_Z_cart.min()) == 0);
    wassert(actual((unsigned)cb->qual_Z_cart.avg()) == 25);
    wassert(actual((unsigned)cb->qual_Z_cart.max()) == 98);
    wassert(actual((unsigned)cb->quota_cart.min()) == 0);
    wassert(actual((unsigned)cb->quota_cart.max()) == 0);
    wassert(actual((unsigned)cb->dato_corr_xy.min()) == 0);
    wassert(actual((unsigned)cb->dato_corr_xy.max()) == 0);
    wassert(actual((unsigned)cb->beam_blocking_xy.min()) == 0);
    wassert(actual((unsigned)cb->beam_blocking_xy.avg()) == 14);
    wassert(actual((unsigned)cb->beam_blocking_xy.max()) == 51);
    wassert(actual((unsigned)cb->elev_fin_xy.min()) == 0);
    wassert(actual((unsigned)cb->elev_fin_xy.avg()) == 0);
    wassert(actual((unsigned)cb->elev_fin_xy.max()) == 3);
    wassert(actual((unsigned)cb->neve_cart.min()) == 0);
    wassert(actual((unsigned)cb->neve_cart.max()) == 0);
    wassert(actual((unsigned)cb->corr_cart.min()) == 0);
    wassert(actual((unsigned)cb->corr_cart.max()) == 0);
    wassert(actual((unsigned)cb->conv_cart.min()) == 0);
    wassert(actual((unsigned)cb->conv_cart.max()) == 0);

    cb->creo_cart_z_lowris();
LOG_INFO("z_out         min avg max :%d %d %d",(unsigned)cb->z_out.min(),(unsigned)cb->z_out.avg(),(unsigned)cb->z_out.max());
LOG_INFO("qual_Z_1x1    min avg max :%d %d %d",(unsigned)cb->qual_Z_1x1.min(),(unsigned)cb->qual_Z_1x1.avg(),(unsigned)cb->qual_Z_1x1.max());
LOG_INFO("quota_1x1     min avg max :%d %d %d",(unsigned)cb->quota_1x1.min(),(unsigned)cb->quota_1x1.avg(),(unsigned)cb->quota_1x1.max());
LOG_INFO("dato_corr_1x1 min avg max :%d %d %d",(unsigned)cb->dato_corr_1x1.min(),(unsigned)cb->dato_corr_1x1.avg(),(unsigned)cb->dato_corr_1x1.max());
LOG_INFO("elev_fin_1x1  min avg max :%d %d %d",(unsigned)cb->elev_fin_1x1.min(),(unsigned)cb->elev_fin_1x1.avg(),(unsigned)cb->elev_fin_1x1.max());
LOG_INFO("beam_blocking_1x1 min avg max :%d %d %d",(unsigned)cb->beam_blocking_1x1.min(),(unsigned)cb->beam_blocking_1x1.avg(),(unsigned)cb->beam_blocking_1x1.max());
LOG_INFO("top_1x1       min avg max :%d %d %d",(unsigned)cb->top_1x1.min(),(unsigned)cb->top_1x1.avg(),(unsigned)cb->top_1x1.max());
LOG_INFO("neve_1x1      min avg max :%d %d %d",(unsigned)cb->neve_1x1.min(),(unsigned)cb->neve_1x1.avg(),(unsigned)cb->neve_1x1.max());
LOG_INFO("corr_1x1      min avg max :%d %d %d",(unsigned)cb->corr_1x1.min(),(unsigned)cb->corr_1x1.avg(),(unsigned)cb->corr_1x1.max());
LOG_INFO("conv_1x1      min avg max :%d %d %d",(unsigned)cb->conv_1x1.min(),(unsigned)cb->conv_1x1.avg(),(unsigned)cb->conv_1x1.max());
    wassert(actual((unsigned)cb->z_out.min()) == 0);
    wassert(actual((unsigned)cb->z_out.avg()) == 63);
    wassert(actual((unsigned)cb->z_out.max()) == 255);
    wassert(actual((unsigned)cb->qual_Z_1x1.min()) == 0);
    wassert(actual((unsigned)cb->qual_Z_1x1.avg()) == 24);
    wassert(actual((unsigned)cb->qual_Z_1x1.max()) == 97);
    wassert(actual((unsigned)cb->quota_1x1.min()) == 128);
    wassert(actual((unsigned)cb->quota_1x1.max()) == 128);
    wassert(actual((unsigned)cb->dato_corr_1x1.min()) == 0);
    wassert(actual((unsigned)cb->dato_corr_1x1.max()) == 0);
    wassert(actual((unsigned)cb->elev_fin_1x1.min()) == 0);
    wassert(actual((unsigned)cb->elev_fin_1x1.avg()) == 0);
    wassert(actual((unsigned)cb->elev_fin_1x1.max()) == 3);
    wassert(actual((unsigned)cb->beam_blocking_1x1.min()) == 0);
    wassert(actual((unsigned)cb->beam_blocking_1x1.avg()) == 13);
    wassert(actual((unsigned)cb->beam_blocking_1x1.max()) == 51);
    wassert(actual((unsigned)cb->top_1x1.min()) == 0);
    wassert(actual((unsigned)cb->top_1x1.avg()) == 3);
    wassert(actual((unsigned)cb->top_1x1.max()) == 42);
    wassert(actual((unsigned)cb->neve_1x1.min()) == 0);
    wassert(actual((unsigned)cb->neve_1x1.max()) == 0);
    wassert(actual((unsigned)cb->corr_1x1.min()) == 0);
    wassert(actual((unsigned)cb->corr_1x1.max()) == 0);
    wassert(actual((unsigned)cb->conv_1x1.min()) == 0);
    wassert(actual((unsigned)cb->conv_1x1.max()) == 0);

    // TODO: scrivo_out_file_bin

    delete cb;
}

template<> template<>
void to::test<8>()
{
	LOG_CATEGORY("Test");
	LOG_INFO ("Start test 10");

    // versione BB_VPR che corrisponde al parametro algo_corto_dev
    static const char* fname = "testdata/DBP2_060220140140_GATTATICO";
    unsetwork();
    setenv("FIRST_LEVEL_FILE", "../dati/FIRST_LEVEL_corto_GAT_2006_INV", 1);
    setenv("DIR_OUT_PP_BLOC", "testdata", 1);
    setenv("FILE_T", "testdata/temperature.txt", 1);
    setenv("FILE_ZERO_TERMICO", "testdata/20140206/0termico.prev", 1);
    setenv("VPR_HEATING", "testdata/vpr_heat_GAT", 1);
    unlink("testdata/vpr_heat_GAT");
    setenv("VPR0_FILE", "testdata/ultimo_vpr", 1);
    unlink("testdata/ultimo_vpr");
    setenv("LAST_VPR", "testdata/last_vpr", 1);
    unlink("testdata/last_vpr");
    setenv("VPR_HMAX", "testdata/vpr_hmax", 1);
    setenv("TEST_VPR", "testdata/test_vpr", 1);
	printwork();

    CUM_BAC* cb = new CUM_BAC("GAT");
    cb->do_quality = true;
    cb->do_beamblocking = true;
    cb->do_declutter = true;
    cb->do_class = false;
    cb->do_bloccorr = false;
    cb->do_vpr = true;
    cb->do_clean= true;
    cb->do_readStaticMap=true;
    std::cout<<"INIZIO - "<<fname<<std::endl;
    cb->read_sp20_volume(fname, 0);
    cb->setup_elaborazione(fname);

    cb->elabora_dato();

    cb->caratterizzo_volume();

    // la combina_profili restituisce 1 se non riesce a costruire un profilo
    // perchè non piove o piove poco
   cb->calcolo_vpr->test_vpr=fopen("testdata/test_vpr","a+");

    int ier = cb->calcolo_vpr->combina_profili();
    wassert(actual(ier) == 0);

    cb->calcolo_vpr->heating = cb->calcolo_vpr->profile_heating();
    wassert(actual(cb->calcolo_vpr->heating) == 0);

    ier = cb->calcolo_vpr->corr_vpr();
    wassert(actual(ier) == 0);

    // TODO: cb->stampa_vpr()

    cb->creo_cart();
    cb->creo_cart_z_lowris();

    // TODO: scrivo_out_file_bin

    delete cb;
}

template<> template<>
void to::test<9>()
{
    // versione BB_VPR_CLASS che corrisponde al parametro algo_corto_dev
    static const char* fname = "testdata/DBP2_060220140140_GATTATICO";
    unsetwork();
    setenv("FIRST_LEVEL_FILE", "../dati/FIRST_LEVEL_corto_GAT_2006_INV", 1);
    setenv("DIR_OUT_PP_BLOC", "testdata", 1);
    setenv("FILE_T", "testdata/temperature.txt", 1);
    setenv("FILE_ZERO_TERMICO", "testdata/20140206/0termico.prev", 1);
    setenv("VPR_HEATING", "testdata/vpr_heat_GAT", 1);
    unlink("testdata/vpr_heat_GAT");
    setenv("VPR0_FILE", "testdata/ultimo_vpr", 1);
    unlink("testdata/ultimo_vpr");
    setenv("LAST_VPR", "testdata/last_vpr", 1);
    unlink("testdata/last_vpr");
    setenv("VPR_HMAX", "testdata/vpr_hmax", 1);

    CUM_BAC* cb = new CUM_BAC("GAT");
    cb->do_quality = true;
    cb->do_beamblocking = true;
    cb->do_declutter = true;
    cb->do_class = true;
    cb->do_bloccorr = false;
    cb->do_vpr = true;
    cb->do_clean= true;
    cb->do_readStaticMap=true;
    cb->read_sp20_volume(fname, 0);
    cb->setup_elaborazione(fname);

    cb->elabora_dato();

    cb->caratterizzo_volume();

    cb->calcolo_vpr->classifica_rain();

    // la combina_profili restituisce 1 se non riesce a costruire un profilo
    // perchè non piove o piove poco
    cb->calcolo_vpr->test_vpr=fopen("testdata/test_vpr","a+");
    int ier = cb->calcolo_vpr->combina_profili();
    wassert(actual(ier) == 0);

    cb->calcolo_vpr->heating = cb->calcolo_vpr->profile_heating();
    wassert(actual(cb->calcolo_vpr->heating) == 0);

    ier = cb->calcolo_vpr->corr_vpr();
    wassert(actual(ier) == 0);

    // TODO: cb->stampa_vpr()

    cb->conversione_convettiva();

    cb->creo_cart();
    cb->creo_cart_z_lowris();

    // TODO: scrivo_out_file_bin

    delete cb;
}

template<> template<>
void to::test<10>()
{
	LOG_CATEGORY("Test");
	LOG_INFO ("Start test 10");

    // versione BB che corrisponde al parametro algo_corto
    static const char* fname = "testdata/DBP2_060220140140_GATTATICO";
    unsetwork();
    setenv("FIRST_LEVEL_FILE", "../dati/FIRST_LEVEL_corto+medio_GAT_INV_2011", 1);
    printwork();

    CUM_BAC* cb = new CUM_BAC("GAT",true,1024);
    cb->do_medium= true;
    cb->do_declutter = true;
    cb->do_readStaticMap=true;
    cb->read_sp20_volume(fname, 0);
    cb->setup_elaborazione(fname);

    cb->elabora_dato();

    VolumeStats stats;
    cb->volume.compute_stats(stats);
    wassert(actual(stats.count_zeros[0]) == 0);
    wassert(actual(stats.count_zeros[1]) == 0);
    wassert(actual(stats.count_zeros[2]) == 0);
    wassert(actual(stats.count_zeros[3]) == 0);
    wassert(actual(stats.count_zeros[4]) == 0);
    wassert(actual(stats.count_ones[0]) ==  61406);
    wassert(actual(stats.count_ones[1]) ==  72916);
    wassert(actual(stats.count_ones[2]) ==  96948);
    wassert(actual(stats.count_ones[3]) == 110168);
    wassert(actual(stats.count_ones[4]) == 119120);
    wassert(actual(stats.count_others[0]) == 136194);
    wassert(actual(stats.count_others[1]) == 124684);
    wassert(actual(stats.count_others[2]) == 100652);
    wassert(actual(stats.count_others[3]) ==  87432);
    wassert(actual(stats.count_others[4]) ==  78480);
    wassert(actual(stats.sum_others[0]) == 18104214);
    wassert(actual(stats.sum_others[1]) == 16004410);
    wassert(actual(stats.sum_others[2]) == 12562044);
    wassert(actual(stats.sum_others[3]) == 10553980);
    wassert(actual(stats.sum_others[4]) ==  9141377);

    cb->creo_cart();

LOG_INFO("cart         min avg max :%d %d %d",(unsigned)cb->cart.min(),(unsigned)cb->cart.avg(),(unsigned)cb->cart.max());
LOG_INFO("cartm        min avg max :%d %d %d",(unsigned)cb->cartm.min(),(unsigned)cb->cartm.avg(),(unsigned)cb->cartm.max());
LOG_INFO("topxy        min avg max :%d %d %d",(unsigned)cb->topxy.min(),(unsigned)cb->topxy.avg(),(unsigned)cb->topxy.max());
LOG_INFO("qual_Z_cart  min avg max :%d %d %d",(unsigned)cb->qual_Z_cart.min(),(unsigned)cb->qual_Z_cart.avg(),(unsigned)cb->qual_Z_cart.max());
LOG_INFO("quota_cart   min avg max :%d %d %d",(unsigned)cb->quota_cart.min(),(unsigned)cb->quota_cart.avg(),(unsigned)cb->quota_cart.max());
LOG_INFO("dato_corr_xy min avg max :%d %d %d",(unsigned)cb->dato_corr_xy.min(),(unsigned)cb->dato_corr_xy.avg(),(unsigned)cb->dato_corr_xy.max());
LOG_INFO("beam_blocking_xy min avg max :%d %d %d",(unsigned)cb->beam_blocking_xy.min(),(unsigned)cb->beam_blocking_xy.avg(),(unsigned)cb->beam_blocking_xy.max());
LOG_INFO("elev_fin_xy  min avg max :%d %d %d",(unsigned)cb->elev_fin_xy.min(),(unsigned)cb->elev_fin_xy.avg(),(unsigned)cb->elev_fin_xy.max());
LOG_INFO("neve_cart    min avg max :%d %d %d",(unsigned)cb->neve_cart.min(),(unsigned)cb->neve_cart.avg(),(unsigned)cb->neve_cart.max());
LOG_INFO("corr_cart    min avg max :%d %d %d",(unsigned)cb->corr_cart.min(),(unsigned)cb->corr_cart.avg(),(unsigned)cb->corr_cart.max());
LOG_INFO("conv_cart    min avg max :%d %d %d",(unsigned)cb->conv_cart.min(),(unsigned)cb->conv_cart.avg(),(unsigned)cb->conv_cart.max());
    wassert(actual((unsigned)(unsigned)cb->cart.min()) == 0); 
    wassert(actual((unsigned)cb->cart.avg()) == 53);
    wassert(actual((unsigned)cb->cart.max()) == 255);
    wassert(actual(cb->cartm.min()) == 0);
    wassert(actual(cb->cartm.max()) == 0);
    wassert(actual((unsigned)cb->topxy.min()) == 0);
    wassert(actual((unsigned)cb->topxy.avg()) == 3);
    wassert(actual((unsigned)cb->topxy.max()) == 76);
    wassert(actual((unsigned)cb->qual_Z_cart.min()) == 0);
    wassert(actual((unsigned)cb->qual_Z_cart.avg()) == 25);
    wassert(actual((unsigned)cb->qual_Z_cart.max()) == 98);
    wassert(actual((unsigned)cb->quota_cart.min()) == 0);
    wassert(actual((unsigned)cb->quota_cart.max()) == 0);
    wassert(actual((unsigned)cb->dato_corr_xy.min()) == 0);
    wassert(actual((unsigned)cb->dato_corr_xy.max()) == 0);
    wassert(actual((unsigned)cb->beam_blocking_xy.min()) == 0);
    wassert(actual((unsigned)cb->beam_blocking_xy.avg()) == 14);
    wassert(actual((unsigned)cb->beam_blocking_xy.max()) == 51);
    wassert(actual((unsigned)cb->elev_fin_xy.min()) == 0);
    wassert(actual((unsigned)cb->elev_fin_xy.avg()) == 0);
    wassert(actual((unsigned)cb->elev_fin_xy.max()) == 3);
    wassert(actual((unsigned)cb->neve_cart.min()) == 0);
    wassert(actual((unsigned)cb->neve_cart.max()) == 0);
    wassert(actual((unsigned)cb->corr_cart.min()) == 0);
    wassert(actual((unsigned)cb->corr_cart.max()) == 0);
    wassert(actual((unsigned)cb->conv_cart.min()) == 0);
    wassert(actual((unsigned)cb->conv_cart.max()) == 0);

    cb->creo_cart_z_lowris();
LOG_INFO("z_out         min avg max :%d %d %d",(unsigned)cb->z_out.min(),(unsigned)cb->z_out.avg(),(unsigned)cb->z_out.max());
LOG_INFO("qual_Z_1x1    min avg max :%d %d %d",(unsigned)cb->qual_Z_1x1.min(),(unsigned)cb->qual_Z_1x1.avg(),(unsigned)cb->qual_Z_1x1.max());
LOG_INFO("quota_1x1     min avg max :%d %d %d",(unsigned)cb->quota_1x1.min(),(unsigned)cb->quota_1x1.avg(),(unsigned)cb->quota_1x1.max());
LOG_INFO("dato_corr_1x1 min avg max :%d %d %d",(unsigned)cb->dato_corr_1x1.min(),(unsigned)cb->dato_corr_1x1.avg(),(unsigned)cb->dato_corr_1x1.max());
LOG_INFO("elev_fin_1x1  min avg max :%d %d %d",(unsigned)cb->elev_fin_1x1.min(),(unsigned)cb->elev_fin_1x1.avg(),(unsigned)cb->elev_fin_1x1.max());
LOG_INFO("beam_blocking_1x1 min avg max :%d %d %d",(unsigned)cb->beam_blocking_1x1.min(),(unsigned)cb->beam_blocking_1x1.avg(),(unsigned)cb->beam_blocking_1x1.max());
LOG_INFO("top_1x1       min avg max :%d %d %d",(unsigned)cb->top_1x1.min(),(unsigned)cb->top_1x1.avg(),(unsigned)cb->top_1x1.max());
LOG_INFO("neve_1x1      min avg max :%d %d %d",(unsigned)cb->neve_1x1.min(),(unsigned)cb->neve_1x1.avg(),(unsigned)cb->neve_1x1.max());
LOG_INFO("corr_1x1      min avg max :%d %d %d",(unsigned)cb->corr_1x1.min(),(unsigned)cb->corr_1x1.avg(),(unsigned)cb->corr_1x1.max());
LOG_INFO("conv_1x1      min avg max :%d %d %d",(unsigned)cb->conv_1x1.min(),(unsigned)cb->conv_1x1.avg(),(unsigned)cb->conv_1x1.max());
    wassert(actual((unsigned)cb->z_out.min()) == 0);
    wassert(actual((unsigned)cb->z_out.avg()) == 63);
    wassert(actual((unsigned)cb->z_out.max()) == 255);
    wassert(actual((unsigned)cb->qual_Z_1x1.min()) == 0);
    wassert(actual((unsigned)cb->qual_Z_1x1.avg()) == 24);
    wassert(actual((unsigned)cb->qual_Z_1x1.max()) == 97);
    wassert(actual((unsigned)cb->quota_1x1.min()) == 128);
    wassert(actual((unsigned)cb->quota_1x1.max()) == 128);
    wassert(actual((unsigned)cb->dato_corr_1x1.min()) == 0);
    wassert(actual((unsigned)cb->dato_corr_1x1.max()) == 0);
    wassert(actual((unsigned)cb->elev_fin_1x1.min()) == 0);
    wassert(actual((unsigned)cb->elev_fin_1x1.avg()) == 0);
    wassert(actual((unsigned)cb->elev_fin_1x1.max()) == 3);
    wassert(actual((unsigned)cb->beam_blocking_1x1.min()) == 0);
    wassert(actual((unsigned)cb->beam_blocking_1x1.avg()) == 13);
    wassert(actual((unsigned)cb->beam_blocking_1x1.max()) == 51);
    wassert(actual((unsigned)cb->top_1x1.min()) == 0);
    wassert(actual((unsigned)cb->top_1x1.avg()) == 3);
    wassert(actual((unsigned)cb->top_1x1.max()) == 42);
    wassert(actual((unsigned)cb->neve_1x1.min()) == 0);
    wassert(actual((unsigned)cb->neve_1x1.max()) == 0);
    wassert(actual((unsigned)cb->corr_1x1.min()) == 0);
    wassert(actual((unsigned)cb->corr_1x1.max()) == 0);
    wassert(actual((unsigned)cb->conv_1x1.min()) == 0);
    wassert(actual((unsigned)cb->conv_1x1.max()) == 0);

    // TODO: scrivo_out_file_bin

    delete cb;
}
#if 0
template<> template<>
void to::test<6>()
{
    // Test di tutto quello che può essere chiamato
    static const char* fname = "testdata/DBP2_070120141530_GATTATICO";

    setenv("FIRST_LEVEL_FILE", "../dati/FIRST_LEVEL_corto_GAT_2006_INV", 1);
    setenv("DIR_OUT_PP_BLOC", "testdata", 1);
    setenv("LAST_VPR", "testdata/last_vpr_GAT", 1);
    setenv("VPR0_FILE", "testdata/vpr_SPC", 1);
    setenv("VPR_HEATING", "vpr_heat_GAT", 1);

    CUM_BAC* cb = new CUM_BAC("GAT");
    cb->do_quality = true;
    cb->do_beamblocking = true;
    cb->do_declutter = false;
    cb->do_bloccorr = true;
    cb->do_vpr = true;
    cb->read_sp20_volume(fname, 0);
    cb->setup_elaborazione(fname);

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

    ier = cb->combina_profili();
    wassert(actual(ier) == 1);

    cb->heating = cb->profile_heating();
    wassert(actual(cb->heating) == 0);

    delete cb;
}
#endif

}
