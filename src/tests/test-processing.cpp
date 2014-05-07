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
    void fill(const T* arr, unsigned size)
    {
        for (unsigned i = 0; i < size; ++i)
            this->count_sample(arr[i], size);
    }

    void fill(const Matrix2D<T>& arr)
    {
        for (int i = 0; i < arr.rows() * arr.cols(); ++i)
            this->count_sample(arr.data()[i], arr.rows() * arr.cols());
    }

    void fill(const PolarScan<T>& arr)
    {
        for (int i = 0; i < arr.rows() * arr.cols(); ++i)
            this->count_sample(arr.data()[i], arr.rows() * arr.cols());
    }

    void fill(const Volume<T>& vol)
    {
        for (unsigned i = 0; i < vol.size(); ++i)
            fill(vol.scan(i));
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

template<typename T>
int avg(const Matrix2D<T>& m)
{
    unsigned size = m.rows() * m.cols();
    double mean = 0;
    for (unsigned i = 0; i < size; ++i)
        mean += (double)m.data()[i] / size;
    return round(mean);
}

}

template<> template<>
void to::test<1>()
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
    qual_stats.fill(*cb->qual);
    wassert(actual((unsigned)qual_stats.first).isfalse());
    wassert(actual((unsigned)qual_stats.min) == 1);
    wassert(actual((unsigned)qual_stats.max) == 99);
    wassert(actual((unsigned)(qual_stats.avg * 100)) == 33038);

    ArrayStats<unsigned char> vpr_stats;
    vpr_stats.fill(*cb->calcolo_vpr->flag_vpr);
    wassert(actual((unsigned)vpr_stats.first).isfalse());
    wassert(actual((unsigned)vpr_stats.min) == 0);
    wassert(actual((unsigned)vpr_stats.max) == 1);
    wassert(actual((unsigned)(vpr_stats.avg * 100)) == 527);

    ArrayStats<unsigned char> top_stats;
    top_stats.fill(cb->top);
    wassert(actual((unsigned)top_stats.first).isfalse());
    wassert(actual((unsigned)top_stats.min) == 0);
    wassert(actual((unsigned)top_stats.max) == 15);
    wassert(actual((unsigned)(top_stats.avg * 100)) == 0);

    cb->calcolo_vpr->classifica_rain();

    /*
    ArrayStats<unsigned char> stratiform_stats;
    stratiform_stats.fill(cb->calcolo_vpr->stratiform);
    wassert(actual((unsigned)stratiform_stats.first).isfalse());
    wassert(actual((unsigned)stratiform_stats.min) == 0);
    wassert(actual((unsigned)stratiform_stats.max) == 1);
    wassert(actual((unsigned)(stratiform_stats.avg * 100)) == 0);
    */

    // calcolo_vpr->esegui_tutto();
    // conversione_convettiva();
    // creo_cart();
    // creo_cart_z_lowris();

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
void to::test<3>()
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
    cb->do_declutter =false ;
    cb->do_class = false;
    cb->do_bloccorr = false;
    cb->do_vpr = true;
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
    stats_qual.fill(*cb->qual);
    wassert(actual((unsigned)stats_qual.first).isfalse());
    //wassert(actual((unsigned)stats_qual.count_zeros) == 1886400);
    wassert(actual((unsigned)stats_qual.count_zeros) == 0);
    wassert(actual((unsigned)stats_qual.count_ones) == 141695);
    wassert(actual((unsigned)stats_qual.min) == 1);
    wassert(actual((unsigned)stats_qual.max) == 99);
    wassert(actual((unsigned)(stats_qual.avg * 100)) == 33042);
    ArrayStats<unsigned char> stats_flag_vpr;
    stats_flag_vpr.fill(*cb->calcolo_vpr->flag_vpr);
    wassert(actual((unsigned)stats_flag_vpr.first).isfalse());
    wassert(actual((unsigned)stats_flag_vpr.count_zeros) == 143635);
    wassert(actual((unsigned)stats_flag_vpr.count_ones) == 1041965);
    wassert(actual((unsigned)(stats_flag_vpr.avg * 100)) == 527);
    ArrayStats<unsigned char> stats_top;
    stats_top.fill(cb->top);
    wassert(actual((unsigned)stats_top.first).isfalse());
    wassert(actual((unsigned)stats_top.count_zeros) == 204765);
    wassert(actual((unsigned)stats_top.count_ones) == 0);
    wassert(actual((unsigned)stats_top.min) == 0);
    wassert(actual((unsigned)stats_top.max) == 15);
    wassert(actual((unsigned)(stats_top.avg * 100)) == 0);

    cb->creo_cart();
    wassert(actual((unsigned)cb->cart.minCoeff()) == 0);
    wassert(actual(avg(cb->cart)) == 1);
    wassert(actual((unsigned)cb->cart.maxCoeff()) == 227);
    wassert(actual(cb->cartm.minCoeff()) == 0);
    wassert(actual(cb->cartm.maxCoeff()) == 0);
    wassert(actual((unsigned)cb->topxy.minCoeff()) == 0);
    wassert(actual(avg(cb->topxy)) == 0);
    wassert(actual((unsigned)cb->topxy.maxCoeff()) == 15);
    wassert(actual((unsigned)cb->qual_Z_cart.minCoeff()) == 0);
    wassert(actual(avg(cb->qual_Z_cart)) == 30);
    wassert(actual((unsigned)cb->qual_Z_cart.maxCoeff()) == 98);
    wassert(actual((unsigned)cb->quota_cart.minCoeff()) == 0);
    wassert(actual((unsigned)cb->quota_cart.maxCoeff()) == 5732);
    wassert(actual((unsigned)cb->dato_corr_xy.minCoeff()) == 0);
    wassert(actual((unsigned)cb->dato_corr_xy.maxCoeff()) == 1);
    wassert(actual((unsigned)cb->beam_blocking_xy.minCoeff()) == 0);
    wassert(actual(avg(cb->beam_blocking_xy)) == 9);
    wassert(actual((unsigned)cb->beam_blocking_xy.maxCoeff()) == 51);
    wassert(actual((unsigned)cb->elev_fin_xy.minCoeff()) == 0);
    wassert(actual(avg(cb->elev_fin_xy)) == 0);
    wassert(actual((unsigned)cb->elev_fin_xy.maxCoeff()) == 4);
    wassert(actual((unsigned)cb->neve_cart.minCoeff()) == 0);
    wassert(actual((unsigned)cb->neve_cart.maxCoeff()) == 1);
    wassert(actual((unsigned)cb->corr_cart.minCoeff()) == 0);
    wassert(actual((unsigned)cb->corr_cart.maxCoeff()) == 0);
    wassert(actual((unsigned)cb->conv_cart.minCoeff()) == 0);
    wassert(actual((unsigned)cb->conv_cart.maxCoeff()) == 0);

    cb->creo_cart_z_lowris();
    wassert(actual((unsigned)cb->z_out.minCoeff()) == 0);
    wassert(actual(avg(cb->z_out)) == 1);
    wassert(actual((unsigned)cb->z_out.maxCoeff()) == 227);
    wassert(actual((unsigned)cb->qual_Z_1x1.minCoeff()) == 0);
    wassert(actual(avg(cb->qual_Z_1x1)) == 30);
    wassert(actual((unsigned)cb->qual_Z_1x1.maxCoeff()) == 98);
    wassert(actual((unsigned)cb->quota_1x1.minCoeff()) == 128);
    wassert(actual((unsigned)cb->quota_1x1.maxCoeff()) == 185);
    wassert(actual((unsigned)cb->dato_corr_1x1.minCoeff()) == 0);
    wassert(actual((unsigned)cb->dato_corr_1x1.maxCoeff()) == 1);
    wassert(actual((unsigned)cb->elev_fin_1x1.minCoeff()) == 0);
    wassert(actual(avg(cb->elev_fin_1x1)) == 0);
    wassert(actual((unsigned)cb->elev_fin_1x1.maxCoeff()) == 3);
    wassert(actual((unsigned)cb->beam_blocking_1x1.minCoeff()) == 0);
    wassert(actual(avg(cb->beam_blocking_1x1)) == 8);
    wassert(actual((unsigned)cb->beam_blocking_1x1.maxCoeff()) == 51);
    wassert(actual((unsigned)cb->top_1x1.minCoeff()) == 0);
    wassert(actual(avg(cb->top_1x1)) == 0);
    wassert(actual((unsigned)cb->top_1x1.maxCoeff()) == 15);
    wassert(actual((unsigned)cb->neve_1x1.minCoeff()) == 0);
    wassert(actual((unsigned)cb->neve_1x1.maxCoeff()) == 1);
    wassert(actual((unsigned)cb->corr_1x1.minCoeff()) == 0);
    wassert(actual((unsigned)cb->corr_1x1.maxCoeff()) == 0);
    wassert(actual((unsigned)cb->conv_1x1.minCoeff()) == 0);
    wassert(actual((unsigned)cb->conv_1x1.maxCoeff()) == 0);

    delete cb;
}

template<> template<>
void to::test<4>()
{
    // versione BB_VPR che corrisponde al parametro algo_corto_dev
	LOG_CATEGORY("Test");
	LOG_INFO ("Start test 4");

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
    stats_qual.fill(*cb->qual);
    wassert(actual((unsigned)stats_qual.first).isfalse());
    wassert(actual((unsigned)stats_qual.count_zeros) == 0);
    wassert(actual((unsigned)stats_qual.count_ones) == 141553);
    wassert(actual((unsigned)stats_qual.min) == 1);
    wassert(actual((unsigned)stats_qual.max) == 99);
    wassert(actual((unsigned)(stats_qual.avg * 100)) == 33051);
    ArrayStats<unsigned char> stats_flag_vpr;
    stats_flag_vpr.fill(*cb->calcolo_vpr->flag_vpr);
    wassert(actual((unsigned)stats_flag_vpr.first).isfalse());
    wassert(actual((unsigned)stats_flag_vpr.count_zeros) == 143320);
    wassert(actual((unsigned)stats_flag_vpr.count_ones) == 1042280);
    wassert(actual((unsigned)(stats_flag_vpr.avg * 100)) == 527);
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
    wassert(actual((unsigned)cb->cart.minCoeff()) == 0);
    wassert(actual(avg(cb->cart)) == 10);
    wassert(actual((unsigned)cb->cart.maxCoeff()) == 227);
    wassert(actual(cb->cartm.minCoeff()) == 0);
    wassert(actual(cb->cartm.maxCoeff()) == 0);
    wassert(actual((unsigned)cb->topxy.minCoeff()) == 0);
    wassert(actual(avg(cb->topxy)) == 0);
    wassert(actual((unsigned)cb->topxy.maxCoeff()) == 36);
    wassert(actual((unsigned)cb->qual_Z_cart.minCoeff()) == 0);
    wassert(actual(avg(cb->qual_Z_cart)) == 30);
    wassert(actual((unsigned)cb->qual_Z_cart.maxCoeff()) == 98);
    wassert(actual((unsigned)cb->quota_cart.minCoeff()) == 0);
    wassert(actual((unsigned)cb->quota_cart.maxCoeff()) == 0);
    wassert(actual((unsigned)cb->dato_corr_xy.minCoeff()) == 0);
    wassert(actual((unsigned)cb->dato_corr_xy.maxCoeff()) == 0);
    wassert(actual((unsigned)cb->beam_blocking_xy.minCoeff()) == 0);
    wassert(actual(avg(cb->beam_blocking_xy)) == 9);
    wassert(actual((unsigned)cb->beam_blocking_xy.maxCoeff()) == 51);
    wassert(actual((unsigned)cb->elev_fin_xy.minCoeff()) == 0);
    wassert(actual(avg(cb->elev_fin_xy)) == 0);
    wassert(actual((unsigned)cb->elev_fin_xy.maxCoeff()) == 3);
    wassert(actual((unsigned)cb->neve_cart.minCoeff()) == 0);
    wassert(actual((unsigned)cb->neve_cart.maxCoeff()) == 1);
    wassert(actual((unsigned)cb->corr_cart.minCoeff()) == 0);
    wassert(actual((unsigned)cb->corr_cart.maxCoeff()) == 0);
    wassert(actual((unsigned)cb->conv_cart.minCoeff()) == 0);
    wassert(actual((unsigned)cb->conv_cart.maxCoeff()) == 0);

    cb->creo_cart_z_lowris();
    wassert(actual((unsigned)cb->z_out.minCoeff()) == 0);
    wassert(actual(avg(cb->z_out)) == 16);
    wassert(actual((unsigned)cb->z_out.maxCoeff()) == 227);
    wassert(actual((unsigned)cb->qual_Z_1x1.minCoeff()) == 0);
    wassert(actual(avg(cb->qual_Z_1x1)) == 29);
    wassert(actual((unsigned)cb->qual_Z_1x1.maxCoeff()) == 98);
    wassert(actual((unsigned)cb->quota_1x1.minCoeff()) == 128);
    wassert(actual((unsigned)cb->quota_1x1.maxCoeff()) == 128);
    wassert(actual((unsigned)cb->dato_corr_1x1.minCoeff()) == 0);
    wassert(actual((unsigned)cb->dato_corr_1x1.maxCoeff()) == 0);
    wassert(actual((unsigned)cb->elev_fin_1x1.minCoeff()) == 0);
    wassert(actual(avg(cb->elev_fin_1x1)) == 0);
    wassert(actual((unsigned)cb->elev_fin_1x1.maxCoeff()) == 3);
    wassert(actual((unsigned)cb->beam_blocking_1x1.minCoeff()) == 0);
    wassert(actual(avg(cb->beam_blocking_1x1)) == 9);
    wassert(actual((unsigned)cb->beam_blocking_1x1.maxCoeff()) == 51);
    wassert(actual((unsigned)cb->top_1x1.minCoeff()) == 0);
    wassert(actual((unsigned)cb->top_1x1.maxCoeff()) == 35);
    wassert(actual(avg(cb->top_1x1)) == 0);
    wassert(actual((unsigned)cb->neve_1x1.minCoeff()) == 0);
    wassert(actual((unsigned)cb->neve_1x1.maxCoeff()) == 1);
    wassert(actual(avg(cb->neve_1x1)) == 1);
    wassert(actual((unsigned)cb->corr_1x1.minCoeff()) == 0);
    wassert(actual((unsigned)cb->corr_1x1.maxCoeff()) == 0);
    wassert(actual((unsigned)cb->conv_1x1.minCoeff()) == 0);
    wassert(actual((unsigned)cb->conv_1x1.maxCoeff()) == 0);

    delete cb;
    LOG_INFO("End test 5");
}

template<> template<>
void to::test<5>()
{
    // versione BB_VPR_CLASS che corrisponde al parametro algo_corto_dev
	LOG_CATEGORY("Test");
	LOG_INFO ("Start test 5");

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
    stats_qual.fill(*cb->qual);
    wassert(actual((unsigned)stats_qual.first).isfalse());
    wassert(actual((unsigned)stats_qual.count_zeros) == 0);
    wassert(actual((unsigned)stats_qual.count_ones) == 141553);
    wassert(actual((unsigned)stats_qual.min) == 1);
    wassert(actual((unsigned)stats_qual.max) == 99);
    wassert(actual((unsigned)(stats_qual.avg * 100)) == 33052);
    ArrayStats<unsigned char> stats_flag_vpr;
    stats_flag_vpr.fill(*cb->calcolo_vpr->flag_vpr);
    wassert(actual((unsigned)stats_flag_vpr.first).isfalse());
    wassert(actual((unsigned)stats_flag_vpr.count_zeros) == 143320);
    wassert(actual((unsigned)stats_flag_vpr.count_ones) == 1042280);
    wassert(actual((unsigned)(stats_flag_vpr.avg * 100)) == 527);
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

LOG_INFO("cart         min avg max :%d %d %d",(unsigned)cb->cart.minCoeff(),avg(cb->cart),(unsigned)cb->cart.maxCoeff());
LOG_INFO("cartm        min avg max :%d %d %d",(unsigned)cb->cartm.minCoeff(),avg(cb->cartm),(unsigned)cb->cartm.maxCoeff());
LOG_INFO("topxy        min avg max :%d %d %d",(unsigned)cb->topxy.minCoeff(),avg(cb->topxy),(unsigned)cb->topxy.maxCoeff());
LOG_INFO("qual_Z_cart  min avg max :%d %d %d",(unsigned)cb->qual_Z_cart.minCoeff(),avg(cb->qual_Z_cart),(unsigned)cb->qual_Z_cart.maxCoeff());
LOG_INFO("quota_cart   min avg max :%d %d %d",(unsigned)cb->quota_cart.minCoeff(),avg(cb->quota_cart),(unsigned)cb->quota_cart.maxCoeff());
LOG_INFO("dato_corr_xy min avg max :%d %d %d",(unsigned)cb->dato_corr_xy.minCoeff(),avg(cb->dato_corr_xy),(unsigned)cb->dato_corr_xy.maxCoeff());
LOG_INFO("beam_blocking_xy min avg max :%d %d %d",(unsigned)cb->beam_blocking_xy.minCoeff(),avg(cb->beam_blocking_xy),(unsigned)cb->beam_blocking_xy.maxCoeff());
LOG_INFO("elev_fin_xy  min avg max :%d %d %d",(unsigned)cb->elev_fin_xy.minCoeff(),avg(cb->elev_fin_xy),(unsigned)cb->elev_fin_xy.maxCoeff());
LOG_INFO("neve_cart    min avg max :%d %d %d",(unsigned)cb->neve_cart.minCoeff(),avg(cb->neve_cart),(unsigned)cb->neve_cart.maxCoeff());
LOG_INFO("corr_cart    min avg max :%d %d %d",(unsigned)cb->corr_cart.minCoeff(),avg(cb->corr_cart),(unsigned)cb->corr_cart.maxCoeff());
LOG_INFO("conv_cart    min avg max :%d %d %d",(unsigned)cb->conv_cart.minCoeff(),avg(cb->conv_cart),(unsigned)cb->conv_cart.maxCoeff());
LOG_INFO("------------------------------");
    wassert(actual((unsigned)cb->cart.minCoeff()) == 0);
    wassert(actual(avg(cb->cart)) == 1);
    wassert(actual((unsigned)cb->cart.maxCoeff()) == 227);
    wassert(actual(cb->cartm.minCoeff()) == 0);
    wassert(actual(cb->cartm.maxCoeff()) == 0);
    wassert(actual((unsigned)cb->topxy.minCoeff()) == 0);
    wassert(actual(avg(cb->topxy)) == 0);
    wassert(actual((unsigned)cb->topxy.maxCoeff()) == 15);
    wassert(actual((unsigned)cb->qual_Z_cart.minCoeff()) == 0);
    wassert(actual(avg(cb->qual_Z_cart)) == 30);
    wassert(actual((unsigned)cb->qual_Z_cart.maxCoeff()) == 98);
    wassert(actual((unsigned)cb->quota_cart.minCoeff()) == 0);
    wassert(actual((unsigned)cb->quota_cart.maxCoeff()) == 0);
    wassert(actual((unsigned)cb->dato_corr_xy.minCoeff()) == 0);
    wassert(actual((unsigned)cb->dato_corr_xy.maxCoeff()) == 0);
    wassert(actual((unsigned)cb->beam_blocking_xy.minCoeff()) == 0);
    wassert(actual(avg(cb->beam_blocking_xy)) == 9);
    wassert(actual((unsigned)cb->beam_blocking_xy.maxCoeff()) == 51);
    wassert(actual((unsigned)cb->elev_fin_xy.minCoeff()) == 0);
    wassert(actual(avg(cb->elev_fin_xy)) == 0);
    wassert(actual((unsigned)cb->elev_fin_xy.maxCoeff()) == 3);
    wassert(actual((unsigned)cb->neve_cart.minCoeff()) == 0);
    wassert(actual((unsigned)cb->neve_cart.maxCoeff()) == 1);
    wassert(actual((unsigned)cb->corr_cart.minCoeff()) == 0);
    wassert(actual((unsigned)cb->corr_cart.maxCoeff()) == 0);
    wassert(actual((unsigned)cb->conv_cart.minCoeff()) == 0);
    wassert(actual((unsigned)cb->conv_cart.maxCoeff()) == 0);

    cb->creo_cart_z_lowris();
LOG_INFO("z_out         min avg max :%d %d %d",(unsigned)cb->z_out.minCoeff(),avg(cb->z_out),(unsigned)cb->z_out.maxCoeff());
LOG_INFO("qual_Z_1x1    min avg max :%d %d %d",(unsigned)cb->qual_Z_1x1.minCoeff(),avg(cb->qual_Z_1x1),(unsigned)cb->qual_Z_1x1.maxCoeff());
LOG_INFO("quota_1x1     min avg max :%d %d %d",(unsigned)cb->quota_1x1.minCoeff(),avg(cb->quota_1x1),(unsigned)cb->quota_1x1.maxCoeff());
LOG_INFO("dato_corr_1x1 min avg max :%d %d %d",(unsigned)cb->dato_corr_1x1.minCoeff(),avg(cb->dato_corr_1x1),(unsigned)cb->dato_corr_1x1.maxCoeff());
LOG_INFO("elev_fin_1x1  min avg max :%d %d %d",(unsigned)cb->elev_fin_1x1.minCoeff(),avg(cb->elev_fin_1x1),(unsigned)cb->elev_fin_1x1.maxCoeff());
LOG_INFO("beam_blocking_1x1 min avg max :%d %d %d",(unsigned)cb->beam_blocking_1x1.minCoeff(),avg(cb->beam_blocking_1x1),(unsigned)cb->beam_blocking_1x1.maxCoeff());
LOG_INFO("top_1x1       min avg max :%d %d %d",(unsigned)cb->top_1x1.minCoeff(),avg(cb->top_1x1),(unsigned)cb->top_1x1.maxCoeff());
LOG_INFO("neve_1x1      min avg max :%d %d %d",(unsigned)cb->neve_1x1.minCoeff(),avg(cb->neve_1x1),(unsigned)cb->neve_1x1.maxCoeff());
LOG_INFO("corr_1x1      min avg max :%d %d %d",(unsigned)cb->corr_1x1.minCoeff(),avg(cb->corr_1x1),(unsigned)cb->corr_1x1.maxCoeff());
LOG_INFO("conv_1x1      min avg max :%d %d %d",(unsigned)cb->conv_1x1.minCoeff(),avg(cb->conv_1x1),(unsigned)cb->conv_1x1.maxCoeff());
    wassert(actual((unsigned)cb->z_out.minCoeff()) == 0);
    wassert(actual(avg(cb->z_out)) == 1);
    wassert(actual((unsigned)cb->z_out.maxCoeff()) == 227);
    wassert(actual((unsigned)cb->qual_Z_1x1.minCoeff()) == 0);
    wassert(actual(avg(cb->qual_Z_1x1)) == 30);
    wassert(actual((unsigned)cb->qual_Z_1x1.maxCoeff()) == 98);
    wassert(actual((unsigned)cb->quota_1x1.minCoeff()) == 128);
    wassert(actual((unsigned)cb->quota_1x1.maxCoeff()) == 128);
    wassert(actual((unsigned)cb->dato_corr_1x1.minCoeff()) == 0);
    wassert(actual((unsigned)cb->dato_corr_1x1.maxCoeff()) == 0);
    wassert(actual((unsigned)cb->elev_fin_1x1.minCoeff()) == 0);
    wassert(actual(avg(cb->elev_fin_1x1)) == 0);
    wassert(actual((unsigned)cb->elev_fin_1x1.maxCoeff()) == 3);
    wassert(actual((unsigned)cb->beam_blocking_1x1.minCoeff()) == 0);
    wassert(actual(avg(cb->beam_blocking_1x1)) == 8);
    wassert(actual((unsigned)cb->beam_blocking_1x1.maxCoeff()) == 51);
    wassert(actual((unsigned)cb->top_1x1.minCoeff()) == 0);
    wassert(actual(avg(cb->top_1x1)) == 0);
    wassert(actual((unsigned)cb->top_1x1.maxCoeff()) == 15);
    wassert(actual((unsigned)cb->neve_1x1.minCoeff()) == 0);
    wassert(actual((unsigned)cb->neve_1x1.maxCoeff()) == 1);
    wassert(actual(avg(cb->neve_1x1)) == 1);
    wassert(actual((unsigned)cb->corr_1x1.minCoeff()) == 0);
    wassert(actual((unsigned)cb->corr_1x1.maxCoeff()) == 0);
    wassert(actual((unsigned)cb->conv_1x1.minCoeff()) == 0);
    wassert(actual((unsigned)cb->conv_1x1.maxCoeff()) == 0);

    // TODO: scrivo_out_file_bin

    delete cb;
}

template<> template<>
void to::test<6>()
{
	LOG_CATEGORY("Test");
	LOG_INFO ("Start test 6");

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
    cb->do_declutter = false;
    cb->do_class = false;
    cb->do_bloccorr = false;
    cb->do_vpr = true;
    cb->do_clean= true;
    cb->do_readStaticMap=true;
    cb->read_sp20_volume(fname, 0);
    cb->setup_elaborazione(fname);
    wassert(actual(cb->calcolo_vpr) != (void*)0);

    cb->elabora_dato();

    VolumeStats stats;
    cb->volume.compute_stats(stats);
    //stats.print(stdout);
    wassert(actual(stats.count_zeros[0]) == 0);
    wassert(actual(stats.count_zeros[1]) == 0);
    wassert(actual(stats.count_zeros[2]) == 0);
    wassert(actual(stats.count_zeros[3]) == 0);
    wassert(actual(stats.count_zeros[4]) == 0);
    wassert(actual(stats.count_ones[0]) ==  66946);
    wassert(actual(stats.count_ones[1]) ==  79821);
    wassert(actual(stats.count_ones[2]) == 104195);
    wassert(actual(stats.count_ones[3]) == 110179);
    wassert(actual(stats.count_ones[4]) == 119120);
    wassert(actual(stats.count_others[0]) == 137854);
    wassert(actual(stats.count_others[1]) == 124979);
    wassert(actual(stats.count_others[2]) == 100605);
    wassert(actual(stats.count_others[3]) ==  87421);
    wassert(actual(stats.count_others[4]) ==  78480);
    wassert(actual(stats.sum_others[0]) == 18301988);
    wassert(actual(stats.sum_others[1]) == 16033156);
    wassert(actual(stats.sum_others[2]) == 12558154);
    wassert(actual(stats.sum_others[3]) == 10553141);
    wassert(actual(stats.sum_others[4]) ==  9141377);

    ArrayStats<unsigned char> beam_blocking_stats;
    beam_blocking_stats.fill(cb->beam_blocking);
//    beam_blocking_stats.print();
    wassert(actual((unsigned)beam_blocking_stats.first).isfalse());
    wassert(actual((unsigned)beam_blocking_stats.min) == 0);
    wassert(actual((unsigned)beam_blocking_stats.max) == 51);
    wassert(actual((unsigned)(beam_blocking_stats.avg * 100)) == 1349);

    unsigned stats_size = cb->grid_stats.size_az * cb->grid_stats.size_beam;

    ArrayStats<unsigned> stat_anap_stats;
    stat_anap_stats.fill(cb->grid_stats.stat_anap, stats_size);
//    stat_anap_stats.print();
    wassert(actual((unsigned)stat_anap_stats.first).isfalse());
    wassert(actual((unsigned)stat_anap_stats.min) == 0);
    wassert(actual((unsigned)stat_anap_stats.max) == 75);

    ArrayStats<unsigned> stat_anap_tot_stats;
    stat_anap_tot_stats.fill(cb->grid_stats.stat_tot, stats_size);
//    stat_anap_tot_stats.print();
    wassert(actual((unsigned)stat_anap_tot_stats.first).isfalse());
    wassert(actual((unsigned)stat_anap_tot_stats.min) ==  0);
    wassert(actual((unsigned)stat_anap_tot_stats.max) == 1000);

    ArrayStats<unsigned> stat_bloc_stats;
    stat_bloc_stats.fill(cb->grid_stats.stat_bloc, stats_size);
//    stat_bloc_stats.print();
    wassert(actual((unsigned)stat_bloc_stats.first).isfalse());
    wassert(actual((unsigned)stat_bloc_stats.min) == 0);
    wassert(actual((unsigned)stat_bloc_stats.max) == 0);

    ArrayStats<unsigned> stat_elev_stats;
    stat_elev_stats.fill(cb->grid_stats.stat_elev, stats_size);
//    stat_elev_stats.print();
    wassert(actual((unsigned)stat_elev_stats.first).isfalse());
    wassert(actual((unsigned)stat_elev_stats.min) == 0);
    wassert(actual((unsigned)stat_elev_stats.max) == 75);

    ArrayStats<unsigned char> dato_corrotto_stats;
    dato_corrotto_stats.fill(cb->dato_corrotto);
//    dato_corrotto_stats.print();
    wassert(actual((unsigned)dato_corrotto_stats.first).isfalse());
    wassert(actual((unsigned)dato_corrotto_stats.min) == 0);
    wassert(actual((unsigned)dato_corrotto_stats.max) == 1);


    cb->caratterizzo_volume();
    wassert(actual(cb->calcolo_vpr) != (void*)0);

    ArrayStats<unsigned char> stats_qual;
    stats_qual.fill(*cb->qual);
//    stats_qual.print();
    wassert(actual((unsigned)stats_qual.first).isfalse());
    wassert(actual((unsigned)stats_qual.count_zeros) == 0);
    wassert(actual((unsigned)stats_qual.count_ones) == 162599);
    wassert(actual((unsigned)stats_qual.min) == 1);
    wassert(actual((unsigned)stats_qual.max) == 99);
    wassert(actual((unsigned)(stats_qual.avg * 100)) == 68792);

    ArrayStats<unsigned char> stats_flag_vpr;
    stats_flag_vpr.fill(*cb->calcolo_vpr->flag_vpr);
//    stats_flag_vpr.print();
    wassert(actual((unsigned)stats_flag_vpr.first).isfalse());
    wassert(actual((unsigned)stats_flag_vpr.count_zeros) == 167079);
    wassert(actual((unsigned)stats_flag_vpr.count_ones) == 1650521);
    wassert(actual((unsigned)(stats_flag_vpr.avg * 100)) == 1018);

    ArrayStats<unsigned char> stats_top;
    stats_top.fill(cb->top);
//    stats_top.print();
    wassert(actual((unsigned)stats_top.first).isfalse());
    wassert(actual((unsigned)stats_top.count_zeros) == 112778);
    wassert(actual((unsigned)stats_top.count_ones) == 159);
    wassert(actual((unsigned)stats_top.min) == 0);
    wassert(actual((unsigned)stats_top.max) == 76);
    wassert(actual((unsigned)(stats_top.avg * 100)) == 554);

    cb->creo_cart();

LOG_INFO("cart         min avg max :%d %d %d",(unsigned)cb->cart.minCoeff(),avg(cb->cart),(unsigned)cb->cart.maxCoeff());
LOG_INFO("cartm        min avg max :%d %d %d",(unsigned)cb->cartm.minCoeff(),avg(cb->cartm),(unsigned)cb->cartm.maxCoeff());
LOG_INFO("topxy        min avg max :%d %d %d",(unsigned)cb->topxy.minCoeff(),avg(cb->topxy),(unsigned)cb->topxy.maxCoeff());
LOG_INFO("qual_Z_cart  min avg max :%d %d %d",(unsigned)cb->qual_Z_cart.minCoeff(),avg(cb->qual_Z_cart),(unsigned)cb->qual_Z_cart.maxCoeff());
LOG_INFO("quota_cart   min avg max :%d %d %d",(unsigned)cb->quota_cart.minCoeff(),avg(cb->quota_cart),(unsigned)cb->quota_cart.maxCoeff());
LOG_INFO("dato_corr_xy min avg max :%d %d %d",(unsigned)cb->dato_corr_xy.minCoeff(),avg(cb->dato_corr_xy),(unsigned)cb->dato_corr_xy.maxCoeff());
LOG_INFO("beam_blocking_xy min avg max :%d %d %d",(unsigned)cb->beam_blocking_xy.minCoeff(),avg(cb->beam_blocking_xy),(unsigned)cb->beam_blocking_xy.maxCoeff());
LOG_INFO("elev_fin_xy  min avg max :%d %d %d",(unsigned)cb->elev_fin_xy.minCoeff(),avg(cb->elev_fin_xy),(unsigned)cb->elev_fin_xy.maxCoeff());
LOG_INFO("neve_cart    min avg max :%d %d %d",(unsigned)cb->neve_cart.minCoeff(),avg(cb->neve_cart),(unsigned)cb->neve_cart.maxCoeff());
LOG_INFO("corr_cart    min avg max :%d %d %d",(unsigned)cb->corr_cart.minCoeff(),avg(cb->corr_cart),(unsigned)cb->corr_cart.maxCoeff());
LOG_INFO("conv_cart    min avg max :%d %d %d",(unsigned)cb->conv_cart.minCoeff(),avg(cb->conv_cart),(unsigned)cb->conv_cart.maxCoeff());
    wassert(actual((unsigned)(unsigned)cb->cart.minCoeff()) == 0); 
    wassert(actual(avg(cb->cart)) == 55);
    wassert(actual((unsigned)cb->cart.maxCoeff()) == 255);
    wassert(actual(cb->cartm.minCoeff()) == 0);
    wassert(actual(cb->cartm.maxCoeff()) == 0);
    wassert(actual((unsigned)cb->topxy.minCoeff()) == 0);
    wassert(actual(avg(cb->topxy)) == 3);
    wassert(actual((unsigned)cb->topxy.maxCoeff()) == 76);
    wassert(actual((unsigned)cb->qual_Z_cart.minCoeff()) == 0);
    wassert(actual(avg(cb->qual_Z_cart)) == 26);
    wassert(actual((unsigned)cb->qual_Z_cart.maxCoeff()) == 98);
    wassert(actual((unsigned)cb->quota_cart.minCoeff()) == 0);
    wassert(actual((unsigned)cb->quota_cart.maxCoeff()) == 5825);
    wassert(actual((unsigned)cb->dato_corr_xy.minCoeff()) == 0);
    wassert(actual((unsigned)cb->dato_corr_xy.maxCoeff()) == 1);
    wassert(actual((unsigned)cb->beam_blocking_xy.minCoeff()) == 0);
    wassert(actual(avg(cb->beam_blocking_xy)) == 14);
    wassert(actual((unsigned)cb->beam_blocking_xy.maxCoeff()) == 51);
    wassert(actual((unsigned)cb->elev_fin_xy.minCoeff()) == 0);
    wassert(actual(avg(cb->elev_fin_xy)) == 0);
    wassert(actual((unsigned)cb->elev_fin_xy.maxCoeff()) == 3);
    wassert(actual((unsigned)cb->neve_cart.minCoeff()) == 0);
    wassert(actual((unsigned)cb->neve_cart.maxCoeff()) == 1);
    wassert(actual((unsigned)cb->corr_cart.minCoeff()) == 0);
    wassert(actual((unsigned)cb->corr_cart.maxCoeff()) == 0);
    wassert(actual((unsigned)cb->conv_cart.minCoeff()) == 0);
    wassert(actual((unsigned)cb->conv_cart.maxCoeff()) == 0);

    cb->creo_cart_z_lowris();
LOG_INFO("z_out         min avg max :%d %d %d",(unsigned)cb->z_out.minCoeff(),avg(cb->z_out),(unsigned)cb->z_out.maxCoeff());
LOG_INFO("qual_Z_1x1    min avg max :%d %d %d",(unsigned)cb->qual_Z_1x1.minCoeff(),avg(cb->qual_Z_1x1),(unsigned)cb->qual_Z_1x1.maxCoeff());
LOG_INFO("quota_1x1     min avg max :%d %d %d",(unsigned)cb->quota_1x1.minCoeff(),avg(cb->quota_1x1),(unsigned)cb->quota_1x1.maxCoeff());
LOG_INFO("dato_corr_1x1 min avg max :%d %d %d",(unsigned)cb->dato_corr_1x1.minCoeff(),avg(cb->dato_corr_1x1),(unsigned)cb->dato_corr_1x1.maxCoeff());
LOG_INFO("elev_fin_1x1  min avg max :%d %d %d",(unsigned)cb->elev_fin_1x1.minCoeff(),avg(cb->elev_fin_1x1),(unsigned)cb->elev_fin_1x1.maxCoeff());
LOG_INFO("beam_blocking_1x1 min avg max :%d %d %d",(unsigned)cb->beam_blocking_1x1.minCoeff(),avg(cb->beam_blocking_1x1),(unsigned)cb->beam_blocking_1x1.maxCoeff());
LOG_INFO("top_1x1       min avg max :%d %d %d",(unsigned)cb->top_1x1.minCoeff(),avg(cb->top_1x1),(unsigned)cb->top_1x1.maxCoeff());
LOG_INFO("neve_1x1      min avg max :%d %d %d",(unsigned)cb->neve_1x1.minCoeff(),avg(cb->neve_1x1),(unsigned)cb->neve_1x1.maxCoeff());
LOG_INFO("corr_1x1      min avg max :%d %d %d",(unsigned)cb->corr_1x1.minCoeff(),avg(cb->corr_1x1),(unsigned)cb->corr_1x1.maxCoeff());
LOG_INFO("conv_1x1      min avg max :%d %d %d",(unsigned)cb->conv_1x1.minCoeff(),avg(cb->conv_1x1),(unsigned)cb->conv_1x1.maxCoeff());
    wassert(actual((unsigned)cb->z_out.minCoeff()) == 0);
    wassert(actual(avg(cb->z_out)) == 66);
    wassert(actual((unsigned)cb->z_out.maxCoeff()) == 255);
    wassert(actual((unsigned)cb->qual_Z_1x1.minCoeff()) == 0);
    wassert(actual(avg(cb->qual_Z_1x1)) == 25);
    wassert(actual((unsigned)cb->qual_Z_1x1.maxCoeff()) == 97);
    wassert(actual((unsigned)cb->quota_1x1.minCoeff()) == 128);
    wassert(actual((unsigned)cb->quota_1x1.maxCoeff()) == 186);
    wassert(actual((unsigned)cb->dato_corr_1x1.minCoeff()) == 0);
    wassert(actual((unsigned)cb->dato_corr_1x1.maxCoeff()) == 1);
    wassert(actual((unsigned)cb->elev_fin_1x1.minCoeff()) == 0);
    wassert(actual(avg(cb->elev_fin_1x1)) == 0);
    wassert(actual((unsigned)cb->elev_fin_1x1.maxCoeff()) == 3);
    wassert(actual((unsigned)cb->beam_blocking_1x1.minCoeff()) == 0);
    wassert(actual(avg(cb->beam_blocking_1x1)) == 15);
    wassert(actual((unsigned)cb->beam_blocking_1x1.maxCoeff()) == 51);
    wassert(actual((unsigned)cb->top_1x1.minCoeff()) == 0);
    wassert(actual(avg(cb->top_1x1)) == 3);
    wassert(actual((unsigned)cb->top_1x1.maxCoeff()) == 42);
    wassert(actual((unsigned)cb->neve_1x1.minCoeff()) == 0);
    wassert(actual((unsigned)cb->neve_1x1.maxCoeff()) == 1);
    wassert(actual((unsigned)cb->corr_1x1.minCoeff()) == 0);
    wassert(actual((unsigned)cb->corr_1x1.maxCoeff()) == 0);
    wassert(actual((unsigned)cb->conv_1x1.minCoeff()) == 0);
    wassert(actual((unsigned)cb->conv_1x1.maxCoeff()) == 0);

    // TODO: scrivo_out_file_bin

    delete cb;
}

template<> template<>
void to::test<7>()
{
	LOG_CATEGORY("Test");
	LOG_INFO ("Start test 7");

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
void to::test<8>()
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
void to::test<9>()
{
	LOG_CATEGORY("Test");
	LOG_INFO ("Start test 9");

    // versione BB che corrisponde al parametro algo_corto
    static const char* fname = "testdata/DBP2_060220140140_GATTATICO";
    unsetwork();
    setenv("FIRST_LEVEL_FILE", "../dati/FIRST_LEVEL_corto+medio_GAT_INV_2011", 1);
    printwork();

    CUM_BAC* cb = new CUM_BAC("GAT",true,1024);
    cb->do_medium= true;
    cb->do_readStaticMap=true;
//    cb->do_zlr_media=true; 
    cb->read_sp20_volume(fname, 0);
    cb->setup_elaborazione(fname);
    cb->assets.write_vpr_heating(0);

    cb->elabora_dato();
LOG_INFO(" Valuto statistica");
    VolumeStats stats;
    cb->volume.compute_stats(stats);
    wassert(actual(stats.count_zeros[0]) == 0);
    wassert(actual(stats.count_zeros[1]) == 0);
    wassert(actual(stats.count_zeros[2]) == 0);
    wassert(actual(stats.count_zeros[3]) == 0);
    wassert(actual(stats.count_zeros[4]) == 0);
    wassert(actual(stats.count_ones[0]) == 186919);
    wassert(actual(stats.count_ones[1]) == 215538);
    wassert(actual(stats.count_ones[2]) == 163730);
    wassert(actual(stats.count_ones[3]) == 109585);
    wassert(actual(stats.count_ones[4]) == 118617);
    wassert(actual(stats.count_others[0]) == 152681);
    wassert(actual(stats.count_others[1]) == 124062);
    wassert(actual(stats.count_others[2]) == 100270);
    wassert(actual(stats.count_others[3]) ==  88015);
    wassert(actual(stats.count_others[4]) ==  78983);
    wassert(actual(stats.sum_others[0]) == 19978778);
    wassert(actual(stats.sum_others[1]) == 15866483);
    wassert(actual(stats.sum_others[2]) == 12519289);
    wassert(actual(stats.sum_others[3]) == 10597453);
    wassert(actual(stats.sum_others[4]) ==  9178563);

    cb->creo_cart();

LOG_INFO("cart         min avg max :%d %d %d",(unsigned)cb->cart.minCoeff(),avg(cb->cart),(unsigned)cb->cart.maxCoeff());
LOG_INFO("cartm        min avg max :%d %d %d",(unsigned)cb->cartm.minCoeff(),avg(cb->cartm),(unsigned)cb->cartm.maxCoeff());
LOG_INFO("topxy        min avg max :%d %d %d",(unsigned)cb->topxy.minCoeff(),avg(cb->topxy),(unsigned)cb->topxy.maxCoeff());
LOG_INFO("qual_Z_cart  min avg max :%d %d %d",(unsigned)cb->qual_Z_cart.minCoeff(),avg(cb->qual_Z_cart),(unsigned)cb->qual_Z_cart.maxCoeff());
LOG_INFO("quota_cart   min avg max :%d %d %d",(unsigned)cb->quota_cart.minCoeff(),avg(cb->quota_cart),(unsigned)cb->quota_cart.maxCoeff());
LOG_INFO("dato_corr_xy min avg max :%d %d %d",(unsigned)cb->dato_corr_xy.minCoeff(),avg(cb->dato_corr_xy),(unsigned)cb->dato_corr_xy.maxCoeff());
LOG_INFO("beam_blocking_xy min avg max :%d %d %d",(unsigned)cb->beam_blocking_xy.minCoeff(),avg(cb->beam_blocking_xy),(unsigned)cb->beam_blocking_xy.maxCoeff());
LOG_INFO("elev_fin_xy  min avg max :%d %d %d",(unsigned)cb->elev_fin_xy.minCoeff(),avg(cb->elev_fin_xy),(unsigned)cb->elev_fin_xy.maxCoeff());
LOG_INFO("neve_cart    min avg max :%d %d %d",(unsigned)cb->neve_cart.minCoeff(),avg(cb->neve_cart),(unsigned)cb->neve_cart.maxCoeff());
LOG_INFO("corr_cart    min avg max :%d %d %d",(unsigned)cb->corr_cart.minCoeff(),avg(cb->corr_cart),(unsigned)cb->corr_cart.maxCoeff());
LOG_INFO("conv_cart    min avg max :%d %d %d",(unsigned)cb->conv_cart.minCoeff(),avg(cb->conv_cart),(unsigned)cb->conv_cart.maxCoeff());
    wassert(actual((unsigned)(unsigned)cb->cart.minCoeff()) == 0); 
    wassert(actual(avg(cb->cart)) == 19);
    wassert(actual((unsigned)cb->cart.maxCoeff()) == 255);
    wassert(actual(cb->cartm.minCoeff()) == 0);
    wassert(actual(cb->cartm.maxCoeff()) == 0);
    wassert(actual((unsigned)cb->topxy.minCoeff()) == 0);
    wassert(actual(avg(cb->topxy)) == 0);
    wassert(actual((unsigned)cb->topxy.maxCoeff()) == 0);
    wassert(actual((unsigned)cb->qual_Z_cart.minCoeff()) == 0);
    wassert(actual(avg(cb->qual_Z_cart)) == 0);
    wassert(actual((unsigned)cb->qual_Z_cart.maxCoeff()) == 0);
    wassert(actual((unsigned)cb->quota_cart.minCoeff()) == 0);
    wassert(actual((unsigned)cb->quota_cart.maxCoeff()) == 0);
    wassert(actual((unsigned)cb->dato_corr_xy.minCoeff()) == 0);
    wassert(actual((unsigned)cb->dato_corr_xy.maxCoeff()) == 0);
    wassert(actual((unsigned)cb->beam_blocking_xy.minCoeff()) == 0);
    wassert(actual(avg(cb->beam_blocking_xy)) == 0);
    wassert(actual((unsigned)cb->beam_blocking_xy.maxCoeff()) == 0);
    wassert(actual((unsigned)cb->elev_fin_xy.minCoeff()) == 0);
    wassert(actual(avg(cb->elev_fin_xy)) == 0);
    wassert(actual((unsigned)cb->elev_fin_xy.maxCoeff()) == 0);
    wassert(actual((unsigned)cb->neve_cart.minCoeff()) == 0);
    wassert(actual((unsigned)cb->neve_cart.maxCoeff()) == 0);
    wassert(actual((unsigned)cb->corr_cart.minCoeff()) == 0);
    wassert(actual((unsigned)cb->corr_cart.maxCoeff()) == 0);
    wassert(actual((unsigned)cb->conv_cart.minCoeff()) == 0);
    wassert(actual((unsigned)cb->conv_cart.maxCoeff()) == 0);

    cb->creo_cart_z_lowris();
LOG_INFO("z_out         min avg max :%d %d %d",(unsigned)cb->z_out.minCoeff(),avg(cb->z_out),(unsigned)cb->z_out.maxCoeff());
LOG_INFO("qual_Z_1x1    min avg max :%d %d %d",(unsigned)cb->qual_Z_1x1.minCoeff(),avg(cb->qual_Z_1x1),(unsigned)cb->qual_Z_1x1.maxCoeff());
LOG_INFO("quota_1x1     min avg max :%d %d %d",(unsigned)cb->quota_1x1.minCoeff(),avg(cb->quota_1x1),(unsigned)cb->quota_1x1.maxCoeff());
LOG_INFO("dato_corr_1x1 min avg max :%d %d %d",(unsigned)cb->dato_corr_1x1.minCoeff(),avg(cb->dato_corr_1x1),(unsigned)cb->dato_corr_1x1.maxCoeff());
LOG_INFO("elev_fin_1x1  min avg max :%d %d %d",(unsigned)cb->elev_fin_1x1.minCoeff(),avg(cb->elev_fin_1x1),(unsigned)cb->elev_fin_1x1.maxCoeff());
LOG_INFO("beam_blocking_1x1 min avg max :%d %d %d",(unsigned)cb->beam_blocking_1x1.minCoeff(),avg(cb->beam_blocking_1x1),(unsigned)cb->beam_blocking_1x1.maxCoeff());
LOG_INFO("top_1x1       min avg max :%d %d %d",(unsigned)cb->top_1x1.minCoeff(),avg(cb->top_1x1),(unsigned)cb->top_1x1.maxCoeff());
LOG_INFO("neve_1x1      min avg max :%d %d %d",(unsigned)cb->neve_1x1.minCoeff(),avg(cb->neve_1x1),(unsigned)cb->neve_1x1.maxCoeff());
LOG_INFO("corr_1x1      min avg max :%d %d %d",(unsigned)cb->corr_1x1.minCoeff(),avg(cb->corr_1x1),(unsigned)cb->corr_1x1.maxCoeff());
LOG_INFO("conv_1x1      min avg max :%d %d %d",(unsigned)cb->conv_1x1.minCoeff(),avg(cb->conv_1x1),(unsigned)cb->conv_1x1.maxCoeff());
    wassert(actual((unsigned)cb->z_out.minCoeff()) == 0);
    wassert(actual(avg(cb->z_out)) == 24);
    wassert(actual((unsigned)cb->z_out.maxCoeff()) == 255);
    wassert(actual((unsigned)cb->qual_Z_1x1.minCoeff()) == 0);
    wassert(actual(avg(cb->qual_Z_1x1)) == 0);
    wassert(actual((unsigned)cb->qual_Z_1x1.maxCoeff()) == 0);
    wassert(actual((unsigned)cb->quota_1x1.minCoeff()) == 0);
    wassert(actual((unsigned)cb->quota_1x1.maxCoeff()) == 0);
    wassert(actual((unsigned)cb->dato_corr_1x1.minCoeff()) == 0);
    wassert(actual((unsigned)cb->dato_corr_1x1.maxCoeff()) == 0);
    wassert(actual((unsigned)cb->elev_fin_1x1.minCoeff()) == 0);
    wassert(actual(avg(cb->elev_fin_1x1)) == 0);
    wassert(actual((unsigned)cb->elev_fin_1x1.maxCoeff()) == 0);
    wassert(actual((unsigned)cb->beam_blocking_1x1.minCoeff()) == 0);
    wassert(actual(avg(cb->beam_blocking_1x1)) == 0);
    wassert(actual((unsigned)cb->beam_blocking_1x1.maxCoeff()) == 0);
    wassert(actual((unsigned)cb->top_1x1.minCoeff()) == 0);
    wassert(actual(avg(cb->top_1x1)) == 0);
    wassert(actual((unsigned)cb->top_1x1.maxCoeff()) == 0);
    wassert(actual((unsigned)cb->neve_1x1.minCoeff()) == 0);
    wassert(actual((unsigned)cb->neve_1x1.maxCoeff()) == 0);
    wassert(actual((unsigned)cb->corr_1x1.minCoeff()) == 0);
    wassert(actual((unsigned)cb->corr_1x1.maxCoeff()) == 0);
    wassert(actual((unsigned)cb->conv_1x1.minCoeff()) == 0);
    wassert(actual((unsigned)cb->conv_1x1.maxCoeff()) == 0);

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
