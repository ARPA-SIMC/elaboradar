#include <wibble/tests.h>
#include "cum_bac.h"
#include "logging.h"
#include <cstdio>
#include <cstdlib>
#include <vector>
#include "setwork.h"
#include "test-utils.h"
#include <unistd.h>

using namespace wibble::tests;
using namespace cumbac;
using namespace testradar;
using namespace std;

namespace tut {

struct process_shar {
    Logging logging;
};
TESTGRP(process);

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
   // stats.print(stdout);
    wassert(actual(stats.count_zeros[0]) == 0);
    wassert(actual(stats.count_zeros[1]) == 0);
    wassert(actual(stats.count_zeros[2]) == 0);
    wassert(actual(stats.count_zeros[3]) == 0);
    wassert(actual(stats.count_zeros[4]) == 0);
    wassert(actual(stats.count_zeros[5]) == 0);
    wassert(actual(stats.count_ones[0]) == 197399);
    wassert(actual(stats.count_ones[1]) == 197481);
    wassert(actual(stats.count_ones[2]) == 197389);
    wassert(actual(stats.count_ones[3]) == 197502);
    wassert(actual(stats.count_ones[4]) == 197474);
    wassert(actual(stats.count_ones[5]) == 196872);
    wassert(actual(stats.count_others[0]) == 201);
    wassert(actual(stats.count_others[1]) == 119);
    wassert(actual(stats.count_others[2]) == 211);
    wassert(actual(stats.count_others[3]) ==  98);
    wassert(actual(stats.count_others[4]) == 126);
    wassert(actual(stats.count_others[5]) == 728);
    wassert(actual(stats.sum_others[0]) == 12206);
    wassert(actual(stats.sum_others[1]) ==  6649);
    wassert(actual(stats.sum_others[2]) == 17759);
    wassert(actual(stats.sum_others[3]) ==  4731);
    wassert(actual(stats.sum_others[4]) ==  7079);
    wassert(actual(stats.sum_others[5]) == 55505);
    cb->caratterizzo_volume();

//    print_stats("cb", *cb, cerr);
    wassert(actual(*cb->qual).statsEqual(1, 52.2, 99));
    wassert(actual(cb->top).statsEqual(0, 196947, 1, 3.32, 15));
    wassert(actual(cb->first_level).statsEqual(0, 0.86, 3));
    wassert(actual(cb->first_level_static).statsEqual(0, 0.54, 3));
    wassert(actual(cb->bb_first_level).statsEqual(0, 0.49, 2));
    wassert(actual(cb->beam_blocking).statsEqual(0, 4.72, 49));
    wassert(actual(cb->dem).statsEqual(-9999, 47.59, 4070.26));
    wassert(actual(cb->quota).statsEqual(1, 1543.48, 5813));
    wassert(actual(cb->dato_corrotto).statsEqual(0, 0, 0));
    wassert(actual(*cb->calcolo_vpr->flag_vpr).statsEqual(0, 0.85, 1));
    wassert(actual(cb->calcolo_vpr->corr_polar).statsEqual(0, 0, 0));
    wassert(actual(cb->calcolo_vpr->neve).statsEqual(0, 0, 0));

    cb->calcolo_vpr->classifica_rain();

    /*
    ArrayStats<unsigned char> stratiform_stats;
    stratiform_stats.fill(cb->calcolo_vpr->stratiform);
    wassert(actual((unsigned)stratiform_stats.all_missing).isfalse());
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
    wassert(actual(stats.count_zeros[5]) == 0);
    wassert(actual(stats.count_ones[0]) == 197399);
    wassert(actual(stats.count_ones[1]) == 197481);
    wassert(actual(stats.count_ones[2]) == 197389);
    wassert(actual(stats.count_ones[3]) == 197502);
    wassert(actual(stats.count_ones[4]) == 197474);
    wassert(actual(stats.count_ones[5]) == 196872);
    wassert(actual(stats.count_others[0]) == 201);
    wassert(actual(stats.count_others[1]) == 119);
    wassert(actual(stats.count_others[2]) == 211);
    wassert(actual(stats.count_others[3]) ==  98);
    wassert(actual(stats.count_others[4]) == 126);
    wassert(actual(stats.count_others[5]) == 728);
    wassert(actual(stats.sum_others[0]) == 12206);
    wassert(actual(stats.sum_others[1]) ==  6649);
    wassert(actual(stats.sum_others[2]) == 17759);
    wassert(actual(stats.sum_others[3]) ==  4731);
    wassert(actual(stats.sum_others[4]) ==  7079);
    wassert(actual(stats.sum_others[5]) == 55505);

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
//     print_stats("cb", *cb, cerr);
    wassert(actual(*cb->qual).statsEqual(1, 53.07, 99));
    wassert(actual(cb->top).statsEqual(0, 196947, 1, 3.32, 15));
    wassert(actual(cb->first_level).statsEqual(0, 0.54, 3));
    wassert(actual(cb->first_level_static).statsEqual(0, 0.54, 3));
    wassert(actual(cb->bb_first_level).statsEqual(0, 0.49, 2));
    wassert(actual(cb->beam_blocking).statsEqual(0, 10.44, 51));
    wassert(actual(cb->dem).statsEqual(-9999, 47.59, 4070.26));
    wassert(actual(cb->quota).statsEqual(1, 1164.04, 5732));
    wassert(actual(cb->dato_corrotto).statsEqual(0, 0, 0));
    wassert(actual(*cb->calcolo_vpr->flag_vpr).statsEqual(0, 0.86, 1));
    wassert(actual(cb->calcolo_vpr->corr_polar).statsEqual(0, 0, 0));
    wassert(actual(cb->calcolo_vpr->neve).statsEqual(0, 0, 0));

    Cart cart(cb->volume.max_beam_size());
    wassert(actual(cart.max_bin) == 494);
    cart.creo_cart(*cb);
    //setenv("DIR_DEBUG", "testdata/", 1);
    //cart.write_out(*cb, cb->assets);
//     print_stats("cart", cart, cout);
    wassert(actual(cart.cart).statsEqual(0, 213127, 1, 1.23, 227));
    wassert(actual(cart.cartm).statsEqual(0, 976144, 0, 0, 0));
    wassert(actual(cart.topxy).statsEqual(0, 974793, 1, 4.43, 15));
    wassert(actual(cart.qual_Z_cart).statsEqual(0, 213127, 1, 36.6, 98));
    wassert(actual(cart.quota_cart).statsEqual(0, 213127, 3, 1519.57, 5732));
    wassert(actual(cart.dato_corr_xy).statsEqual(0, 976144, 0, 0, 0));
    wassert(actual(cart.beam_blocking_xy).statsEqual(0, 691190, 1, 47.08, 51));
    wassert(actual(cart.elev_fin_xy).statsEqual(0, 761145, 1, 1.15, 3));
    wassert(actual(cart.neve_cart).statsEqual(0, 213127, 1, 1, 1));
    wassert(actual(cart.corr_cart).statsEqual(0, 976144, 0, 0, 0));
    wassert(actual(cart.conv_cart).statsEqual(0, 976144, 0, 0, 0));

    CartLowris clow(cb->do_medium ? 512: 256, *cb, cart);
    wassert(actual(clow.CART_DIM_ZLR) == 256);
    wassert(actual(clow.ZLR_N_ELEMENTARY_PIXEL) == 4);
    wassert(actual(clow.ZLR_OFFSET) == -18);
    clow.creo_cart_z_lowris();
//     print_stats("clow", clow, cerr);
    wassert(actual(clow.z_out).statsEqual(0, 17489, 1, 1.78, 227));
    wassert(actual(clow.qual_Z_1x1).statsEqual(0, 17489, 1, 36.23, 98));
    wassert(actual(clow.quota_1x1).statsEqual(0, 0, 128, 138.79, 185));
    wassert(actual(clow.dato_corr_1x1).statsEqual(0, 65536, 0, 0, 0));
    wassert(actual(clow.elev_fin_1x1).statsEqual(0, 52275, 1, 1.15, 3));
    wassert(actual(clow.beam_blocking_1x1).statsEqual(0, 47339, 1, 47.16, 51));
    wassert(actual(clow.top_1x1).statsEqual(0, 65429, 1, 4.68, 15));
    wassert(actual(clow.neve_1x1).statsEqual(0, 17489, 1, 1, 1));
    wassert(actual(clow.corr_1x1).statsEqual(0, 65536, 0, 0, 0));
    wassert(actual(clow.conv_1x1).statsEqual(0, 65536, 0, 0, 0));

    delete cb;
}

template<> template<>
void to::test<4>()
{
    // versione BB_VPR che corrisponde al parametro algo_corto_dev
	LOG_CATEGORY("Test");

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

    // la combina_profili restituisce 1 se non riesce a costruire un profilo
    // perchè non piove o piove poco
    int ier = cb->calcolo_vpr->combina_profili();
    wassert(actual(ier) == 1);

    cb->calcolo_vpr->heating = cb->calcolo_vpr->profile_heating();
    wassert(actual(cb->calcolo_vpr->heating) == 0);

//    print_stats("cb", *cb, cerr);
    wassert(actual(*cb->qual).statsEqual(1, 53.07, 99));
    wassert(actual(cb->top).statsEqual(0, 189348, 1, 7.12, 36));
    wassert(actual(cb->first_level).statsEqual(0, 0.54, 3));
    wassert(actual(cb->first_level_static).statsEqual(0, 0.54, 3));
    wassert(actual(cb->bb_first_level).statsEqual(0, 0.49, 2));
    wassert(actual(cb->beam_blocking).statsEqual(0, 10.44, 51));
    wassert(actual(cb->dem).statsEqual(-9999, 47.59, 4070.26));
    wassert(actual(cb->quota).statsEqual(0, 0, 0));
    wassert(actual(cb->dato_corrotto).statsEqual(0, 0, 0));
    wassert(actual(*cb->calcolo_vpr->flag_vpr).statsEqual(0, 0.7, 1));
    wassert(actual(cb->calcolo_vpr->corr_polar).statsEqual(0, 0, 0));
    wassert(actual(cb->calcolo_vpr->neve).statsEqual(0, 0, 0));

    //ier = cb->calcolo_vpr->corr_vpr();
    //wassert(actual(ier) == 1);

    // TODO: cb->stampa_vpr()

    Cart cart(cb->volume.max_beam_size());
    wassert(actual(cart.max_bin) == 494);
    cart.creo_cart(*cb);
//    print_stats("cart", cart, cerr);
    wassert(actual(cart.cart).statsEqual(0, 213127, 1, 13.66, 227));
    wassert(actual(cart.cartm).statsEqual(0, 976144, 0, 0, 0));
    wassert(actual(cart.topxy).statsEqual(0, 945974, 1, 11.08, 36));
    wassert(actual(cart.qual_Z_cart).statsEqual(0, 213127, 1, 36.52, 98));
    wassert(actual(cart.quota_cart).statsEqual(0, 976144, 0, 0, 0));
    wassert(actual(cart.dato_corr_xy).statsEqual(0, 976144, 0, 0, 0));
    wassert(actual(cart.beam_blocking_xy).statsEqual(0, 689781, 1, 47.01, 51));
    wassert(actual(cart.elev_fin_xy).statsEqual(0, 764612, 1, 1.15, 3));
    wassert(actual(cart.neve_cart).statsEqual(0, 213127, 1, 1, 1));
    wassert(actual(cart.corr_cart).statsEqual(0, 976144, 0, 0, 0));
    wassert(actual(cart.conv_cart).statsEqual(0, 976144, 0, 0, 0));

    CartLowris clow(cb->do_medium ? 512: 256, *cb, cart);
    wassert(actual(clow.ZLR_N_ELEMENTARY_PIXEL) == 4);
    wassert(actual(clow.ZLR_OFFSET) == -18);
    clow.creo_cart_z_lowris();
//    print_stats("clow", clow, cerr);
    wassert(actual(clow.z_out).statsEqual(0, 17489, 1, 22.08, 227));
    wassert(actual(clow.qual_Z_1x1).statsEqual(0, 17489, 1, 35.63, 98));
    wassert(actual(clow.quota_1x1).statsEqual(0, 0, 128, 128, 128));
    wassert(actual(clow.dato_corr_1x1).statsEqual(0, 65536, 0, 0, 0));
    wassert(actual(clow.elev_fin_1x1).statsEqual(0, 53781, 1, 1.13, 3));
    wassert(actual(clow.beam_blocking_1x1).statsEqual(0, 46773, 1, 47.02, 51));
    wassert(actual(clow.top_1x1).statsEqual(0, 62875, 1, 11.59, 35));
    wassert(actual(clow.neve_1x1).statsEqual(0, 17489, 1, 1, 1));
    wassert(actual(clow.corr_1x1).statsEqual(0, 65536, 0, 0, 0));
    wassert(actual(clow.conv_1x1).statsEqual(0, 65536, 0, 0, 0));

    delete cb;
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
    cb->calcolo_vpr->classifica_rain();

    // la combina_profili restituisce 1 se non riesce a costruire un profilo
    // perchè non piove o piove poco
    int ier = cb->calcolo_vpr->combina_profili();
    wassert(actual(ier) == 1);

    cb->calcolo_vpr->heating = cb->calcolo_vpr->profile_heating();
    wassert(actual(cb->calcolo_vpr->heating) == 0);

    ier = cb->calcolo_vpr->corr_vpr();
    wassert(actual(ier) == 1); // TODO: cosa deve dare?

    //     print_stats("cb", *cb, cerr);
    wassert(actual(*cb->qual).statsEqual(1, 53.07, 99));
    wassert(actual(cb->top).statsEqual(0, 196947, 1, 3.32, 15));
    wassert(actual(cb->first_level).statsEqual(0, 0.54, 3));
    wassert(actual(cb->first_level_static).statsEqual(0, 0.54, 3));
    wassert(actual(cb->bb_first_level).statsEqual(0, 0.49, 2));
    wassert(actual(cb->beam_blocking).statsEqual(0, 10.44, 51));
    wassert(actual(cb->dem).statsEqual(-9999, 47.59, 4070.26));
    wassert(actual(cb->quota).statsEqual(0, 0, 0));
    wassert(actual(cb->dato_corrotto).statsEqual(0, 0, 0));
    wassert(actual(*cb->calcolo_vpr->flag_vpr).statsEqual(0, 0.7, 1));
    wassert(actual(cb->calcolo_vpr->corr_polar).statsEqual(0, 0, 0));
    wassert(actual(cb->calcolo_vpr->neve).statsEqual(0, 0, 0));

    // TODO: cb->stampa_vpr()

    cb->conversione_convettiva();

    Cart cart(cb->volume.max_beam_size());
    wassert(actual(cart.max_bin) == 494);
    cart.creo_cart(*cb);
//     print_stats("cart", cart, cerr);
    wassert(actual(cart.cart).statsEqual(0, 213127, 1, 1.23, 227));
    wassert(actual(cart.cartm).statsEqual(0, 976144, 0, 0, 0));
    wassert(actual(cart.topxy).statsEqual(0, 974793, 1, 4.43, 15));
    wassert(actual(cart.qual_Z_cart).statsEqual(0, 213127, 1, 36.6, 98));
    wassert(actual(cart.quota_cart).statsEqual(0, 976144, 0, 0, 0));
    wassert(actual(cart.dato_corr_xy).statsEqual(0, 976144, 0, 0, 0));
    wassert(actual(cart.beam_blocking_xy).statsEqual(0, 691190, 1, 47.08, 51));
    wassert(actual(cart.elev_fin_xy).statsEqual(0, 761145, 1, 1.15, 3));
    wassert(actual(cart.neve_cart).statsEqual(0, 213127, 1, 1, 1));
    wassert(actual(cart.corr_cart).statsEqual(0, 976144, 0, 0, 0));
    wassert(actual(cart.conv_cart).statsEqual(0, 976144, 0, 0, 0));

    CartLowris clow(cb->do_medium ? 512: 256, *cb, cart);
    wassert(actual(clow.ZLR_N_ELEMENTARY_PIXEL) == 4);
    wassert(actual(clow.ZLR_OFFSET) == -18);
    clow.creo_cart_z_lowris();
//     print_stats("clow", clow, cerr);
    wassert(actual(clow.z_out).statsEqual(0, 17489, 1, 1.78, 227));
    wassert(actual(clow.qual_Z_1x1).statsEqual(0, 17489, 1, 36.23, 98));
    wassert(actual(clow.quota_1x1).statsEqual(0, 0, 128, 128, 128));
    wassert(actual(clow.dato_corr_1x1).statsEqual(0, 65536, 0, 0, 0));
    wassert(actual(clow.elev_fin_1x1).statsEqual(0, 52275, 1, 1.15, 3));
    wassert(actual(clow.beam_blocking_1x1).statsEqual(0, 47339, 1, 47.16, 51));
    wassert(actual(clow.top_1x1).statsEqual(0, 65429, 1, 4.68, 15));
    wassert(actual(clow.neve_1x1).statsEqual(0, 17489, 1, 1, 1));
    wassert(actual(clow.corr_1x1).statsEqual(0, 65536, 0, 0, 0));
    wassert(actual(clow.conv_1x1).statsEqual(0, 65536, 0, 0, 0));
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
//    stats.print(stdout);
    wassert(actual(stats.count_zeros[0]) == 0);
    wassert(actual(stats.count_zeros[1]) == 0);
    wassert(actual(stats.count_zeros[2]) == 0);
    wassert(actual(stats.count_zeros[3]) == 0);
    wassert(actual(stats.count_zeros[4]) == 0);
    wassert(actual(stats.count_ones[0]) ==  66267);
    wassert(actual(stats.count_ones[1]) ==  79736);
    wassert(actual(stats.count_ones[2]) == 104124);
    wassert(actual(stats.count_ones[3]) == 110168);
    wassert(actual(stats.count_ones[4]) == 119120);
    wassert(actual(stats.count_others[0]) == 138533);
    wassert(actual(stats.count_others[1]) == 125064);
    wassert(actual(stats.count_others[2]) == 100676);
    wassert(actual(stats.count_others[3]) ==  87432);
    wassert(actual(stats.count_others[4]) ==  78480);
    wassert(actual(stats.sum_others[0]) == 18376671);
    wassert(actual(stats.sum_others[1]) == 16042047);
    wassert(actual(stats.sum_others[2]) == 12564340);
    wassert(actual(stats.sum_others[3]) == 10553980);
    wassert(actual(stats.sum_others[4]) ==  9141377);

    wassert(actual(cb->beam_blocking).statsEqual(0, 10.36, 51));

    unsigned stats_size = cb->grid_stats.size_az * cb->grid_stats.size_beam;

    ArrayStats<unsigned> stat_anap_stats;
    stat_anap_stats.fill(cb->grid_stats.stat_anap, stats_size);
    //stat_anap_stats.print();
    wassert(actual((unsigned)stat_anap_stats.all_missing).isfalse());
    wassert(actual((unsigned)stat_anap_stats.min) == 0);
    wassert(actual((unsigned)stat_anap_stats.max) == 0);

    ArrayStats<unsigned> stat_anap_tot_stats;
    stat_anap_tot_stats.fill(cb->grid_stats.stat_tot, stats_size);
    //stat_anap_tot_stats.print();
    wassert(actual((unsigned)stat_anap_tot_stats.all_missing).isfalse());
    wassert(actual((unsigned)stat_anap_tot_stats.min) ==  0);
    wassert(actual((unsigned)stat_anap_tot_stats.max) == 1000);

    ArrayStats<unsigned> stat_bloc_stats;
    stat_bloc_stats.fill(cb->grid_stats.stat_bloc, stats_size);
    //stat_bloc_stats.print();
    wassert(actual((unsigned)stat_bloc_stats.all_missing).isfalse());
    wassert(actual((unsigned)stat_bloc_stats.min) == 0);
    wassert(actual((unsigned)stat_bloc_stats.max) == 0);

    ArrayStats<unsigned> stat_elev_stats;
    stat_elev_stats.fill(cb->grid_stats.stat_elev, stats_size);
    //stat_elev_stats.print();
    wassert(actual((unsigned)stat_elev_stats.all_missing).isfalse());
    wassert(actual((unsigned)stat_elev_stats.min) == 0);
    wassert(actual((unsigned)stat_elev_stats.max) == 0);

    ArrayStats<unsigned char> dato_corrotto_stats;
    dato_corrotto_stats.fill(cb->dato_corrotto);
    //dato_corrotto_stats.print();
    wassert(actual((unsigned)dato_corrotto_stats.all_missing).isfalse());
    wassert(actual((unsigned)dato_corrotto_stats.min) == 0);
    wassert(actual((unsigned)dato_corrotto_stats.max) == 0);


    cb->caratterizzo_volume();
    wassert(actual(cb->calcolo_vpr) != (void*)0);
//    print_stats("cb", *cb, cerr);
    wassert(actual(*cb->qual).statsEqual(1, 58.81, 99));
    wassert(actual(cb->top).statsEqual(0, 97281, 1, 11.99, 76));
    wassert(actual(cb->first_level).statsEqual(0, 0.54, 3));
    wassert(actual(cb->first_level_static).statsEqual(0, 0.54, 3));
    wassert(actual(cb->bb_first_level).statsEqual(0, 0.48, 2));
    wassert(actual(cb->beam_blocking).statsEqual(0, 10.36, 51));
    wassert(actual(cb->dem).statsEqual(-9999, 47.59, 4070.26));
    wassert(actual(cb->quota).statsEqual(1, 1195.8, 5825));
    wassert(actual(cb->dato_corrotto).statsEqual(0, 0, 0));
    wassert(actual(*cb->calcolo_vpr->flag_vpr).statsEqual(0, 0.91, 1));
    wassert(actual(cb->calcolo_vpr->corr_polar).statsEqual(0, 0, 0));
    wassert(actual(cb->calcolo_vpr->neve).statsEqual(0, 0, 0));

    Cart cart(cb->volume.max_beam_size());
    wassert(actual(cart.max_bin) == 512);
    cart.creo_cart(*cb);
//    print_stats("cart", cart, cerr);
    wassert(actual(cart.cart).statsEqual(0, 228715, 1, 70.46, 255));
    wassert(actual(cart.cartm).statsEqual(0, 1048576, 0, 0, 0));
    wassert(actual(cart.topxy).statsEqual(0, 770367, 2, 12.18, 76));
    wassert(actual(cart.qual_Z_cart).statsEqual(0, 228715, 1, 32.79, 98));
    wassert(actual(cart.quota_cart).statsEqual(0, 228715, 3, 1533.9, 5825));
    wassert(actual(cart.dato_corr_xy).statsEqual(0, 1048576, 0, 0, 0));
    wassert(actual(cart.beam_blocking_xy).statsEqual(0, 729618, 1, 46.64, 51));
    wassert(actual(cart.elev_fin_xy).statsEqual(0, 847848, 1, 1.15, 3));
    wassert(actual(cart.neve_cart).statsEqual(0, 228715, 1, 1, 1));
    wassert(actual(cart.corr_cart).statsEqual(0, 1048576, 0, 0, 0));
    wassert(actual(cart.conv_cart).statsEqual(0, 1048576, 0, 0, 0));

    CartLowris clow(cb->do_medium ? 512: 256, *cb, cart);
    wassert(actual(clow.ZLR_N_ELEMENTARY_PIXEL) == 4);
    wassert(actual(clow.ZLR_OFFSET) == 0);
    clow.creo_cart_z_lowris();
//     print_stats("clow", clow, cerr);
    wassert(actual(clow.z_out).statsEqual(0, 13904, 1, 84.01, 255));
    wassert(actual(clow.qual_Z_1x1).statsEqual(0, 13904, 1, 32.05, 97));
    wassert(actual(clow.quota_1x1).statsEqual(0, 0, 128, 139.18, 186));
    wassert(actual(clow.dato_corr_1x1).statsEqual(0, 65536, 0, 0, 0));
    wassert(actual(clow.elev_fin_1x1).statsEqual(0, 55671, 1, 1.16, 3));
    wassert(actual(clow.beam_blocking_1x1).statsEqual(0, 44590, 1, 46.32, 51));
    wassert(actual(clow.top_1x1).statsEqual(0, 45451, 2, 12.33, 42));
    wassert(actual(clow.neve_1x1).statsEqual(0, 13904, 1, 1, 1));
    wassert(actual(clow.corr_1x1).statsEqual(0, 65536, 0, 0, 0));
    wassert(actual(clow.conv_1x1).statsEqual(0, 65536, 0, 0, 0));

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

    //print_stats("cb", *cb, cerr);
    wassert(actual(*cb->qual).statsEqual(1, 58.81, 99));
    wassert(actual(cb->top).statsEqual(0, 97281, 1, 11.99, 76));
    wassert(actual(cb->first_level).statsEqual(0, 0.54, 3));
    wassert(actual(cb->first_level_static).statsEqual(0, 0.54, 3));
    wassert(actual(cb->bb_first_level).statsEqual(0, 0.48, 2));
    wassert(actual(cb->beam_blocking).statsEqual(0, 10.36, 51));
    wassert(actual(cb->dem).statsEqual(-9999, 47.59, 4070.26));
    wassert(actual(cb->quota).statsEqual(0, 0, 0));
    wassert(actual(cb->dato_corrotto).statsEqual(0, 0, 0));
    wassert(actual(*cb->calcolo_vpr->flag_vpr).statsEqual(0, 0.78, 1));
    wassert(actual(cb->calcolo_vpr->corr_polar).statsEqual(0, 0, 0));
    wassert(actual(cb->calcolo_vpr->neve).statsEqual(0, 0, 0));
    // TODO: cb->stampa_vpr()

    Cart cart(cb->volume.max_beam_size());
    wassert(actual(cart.max_bin) == 512);
    cart.creo_cart(*cb);
//     print_stats("cart", cart, cerr);
    wassert(actual(cart.cart).statsEqual(0, 228715, 1, 70.46, 255));
    wassert(actual(cart.cartm).statsEqual(0, 1048576, 0, 0, 0));
    wassert(actual(cart.topxy).statsEqual(0, 770367, 2, 12.18, 76));
    wassert(actual(cart.qual_Z_cart).statsEqual(0, 228715, 1, 32.79, 98));
    wassert(actual(cart.quota_cart).statsEqual(0, 1048576, 0, 0, 0));
    wassert(actual(cart.dato_corr_xy).statsEqual(0, 1048576, 0, 0, 0));
    wassert(actual(cart.beam_blocking_xy).statsEqual(0, 729618, 1, 46.64, 51));
    wassert(actual(cart.elev_fin_xy).statsEqual(0, 847848, 1, 1.15, 3));
    wassert(actual(cart.neve_cart).statsEqual(0, 228715, 1, 1, 1));
    wassert(actual(cart.corr_cart).statsEqual(0, 1048576, 0, 0, 0));
    wassert(actual(cart.conv_cart).statsEqual(0, 1048576, 0, 0, 0));

    CartLowris clow(cb->do_medium ? 512: 256, *cb, cart);
    wassert(actual(clow.ZLR_N_ELEMENTARY_PIXEL) == 4);
    wassert(actual(clow.ZLR_OFFSET) == 0);
    clow.creo_cart_z_lowris();
//     print_stats("clow", clow, cerr);
    wassert(actual(clow.z_out).statsEqual(0, 13904, 1, 84.01, 255));
    wassert(actual(clow.qual_Z_1x1).statsEqual(0, 13904, 1, 32.05, 97));
    wassert(actual(clow.quota_1x1).statsEqual(0, 0, 128, 128, 128));
    wassert(actual(clow.dato_corr_1x1).statsEqual(0, 65536, 0, 0, 0));
    wassert(actual(clow.elev_fin_1x1).statsEqual(0, 55671, 1, 1.16, 3));
    wassert(actual(clow.beam_blocking_1x1).statsEqual(0, 44590, 1, 46.32, 51));
    wassert(actual(clow.top_1x1).statsEqual(0, 45451, 2, 12.33, 42));
    wassert(actual(clow.neve_1x1).statsEqual(0, 13904, 1, 1, 1));
    wassert(actual(clow.corr_1x1).statsEqual(0, 65536, 0, 0, 0));
    wassert(actual(clow.conv_1x1).statsEqual(0, 65536, 0, 0, 0));

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
    //print_stats("cb", *cb, cerr);
    wassert(actual(*cb->qual).statsEqual(1, 58.81, 99));
    wassert(actual(cb->top).statsEqual(0, 97281, 1, 11.99, 76));
    wassert(actual(cb->first_level).statsEqual(0, 0.54, 3));
    wassert(actual(cb->first_level_static).statsEqual(0, 0.54, 3));
    wassert(actual(cb->bb_first_level).statsEqual(0, 0.48, 2));
    wassert(actual(cb->beam_blocking).statsEqual(0, 10.36, 51));
    wassert(actual(cb->dem).statsEqual(-9999, 47.59, 4070.26));
    wassert(actual(cb->quota).statsEqual(0, 0, 0));
    wassert(actual(cb->dato_corrotto).statsEqual(0, 0, 0));
    wassert(actual(*cb->calcolo_vpr->flag_vpr).statsEqual(0, 0.78, 1));
    wassert(actual(cb->calcolo_vpr->corr_polar).statsEqual(0, 0, 0));
    wassert(actual(cb->calcolo_vpr->neve).statsEqual(0, 0, 0));

    Cart cart(cb->volume.max_beam_size());
    wassert(actual(cart.max_bin) == 512);
    cart.creo_cart(*cb);
//    print_stats("cart", cart, cerr);
    wassert(actual(cart.cart).statsEqual(0, 228715, 1, 70.46, 255));
    wassert(actual(cart.cartm).statsEqual(0, 1048576, 0, 0, 0));
    wassert(actual(cart.topxy).statsEqual(0, 770367, 2, 12.18, 76));
    wassert(actual(cart.qual_Z_cart).statsEqual(0, 228715, 1, 32.79, 98));
    wassert(actual(cart.quota_cart).statsEqual(0, 1048576, 0, 0, 0));
    wassert(actual(cart.dato_corr_xy).statsEqual(0, 1048576, 0, 0, 0));
    wassert(actual(cart.beam_blocking_xy).statsEqual(0, 729618, 1, 46.64, 51));
    wassert(actual(cart.elev_fin_xy).statsEqual(0, 847848, 1, 1.15, 3));
    wassert(actual(cart.neve_cart).statsEqual(0, 228715, 1, 1, 1));
    wassert(actual(cart.corr_cart).statsEqual(0, 1048576, 0, 0, 0));
    wassert(actual(cart.conv_cart).statsEqual(0, 1048565, 100, 100, 100));

    CartLowris clow(cb->do_medium ? 512: 256, *cb, cart);
    wassert(actual(clow.ZLR_N_ELEMENTARY_PIXEL) == 4);
    wassert(actual(clow.ZLR_OFFSET) == 0);
    clow.creo_cart_z_lowris();
//     print_stats("clow", clow, cerr);
    wassert(actual(clow.z_out).statsEqual(0, 13904, 1, 84.01, 255));
    wassert(actual(clow.qual_Z_1x1).statsEqual(0, 13904, 1, 32.05, 97));
    wassert(actual(clow.quota_1x1).statsEqual(0, 0, 128, 128, 128));
    wassert(actual(clow.dato_corr_1x1).statsEqual(0, 65536, 0, 0, 0));
    wassert(actual(clow.elev_fin_1x1).statsEqual(0, 55671, 1, 1.16, 3));
    wassert(actual(clow.beam_blocking_1x1).statsEqual(0, 44590, 1, 46.32, 51));
    wassert(actual(clow.top_1x1).statsEqual(0, 45451, 2, 12.33, 42));
    wassert(actual(clow.neve_1x1).statsEqual(0, 13904, 1, 1, 1));
    wassert(actual(clow.corr_1x1).statsEqual(0, 65536, 0, 0, 0));
    wassert(actual(clow.conv_1x1).statsEqual(0, 65533, 100, 100, 100));

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
    VolumeStats stats;
    cb->volume.compute_stats(stats);
//    stats.print(stdout);
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
    wassert(actual(stats.sum_others[0]) == 19979741);
    wassert(actual(stats.sum_others[1]) == 15866483);
    wassert(actual(stats.sum_others[2]) == 12519289);
    wassert(actual(stats.sum_others[3]) == 10597453);
    wassert(actual(stats.sum_others[4]) ==  9178563);
   // print_stats("cb", *cb, cerr);

    Cart cart(cb->volume.max_beam_size());
    wassert(actual(cart.max_bin) == 849);
    cart.creo_cart(*cb);
//    print_stats("cart", cart, cerr);
    wassert(actual(cart.cart).statsEqual(0, 624947, 1, 36.03, 255));
    wassert(actual(cart.cartm).statsEqual(0, 2883204, 0, 0, 0));
    wassert(actual(cart.topxy).statsEqual(0, 2551711, 2, 15.55, 76));
    wassert(actual(cart.qual_Z_cart).statsEqual(0, 2883204, 0, 0, 0));
    wassert(actual(cart.quota_cart).statsEqual(0, 2883204, 0, 0, 0));
    wassert(actual(cart.dato_corr_xy).statsEqual(0, 2883204, 0, 0, 0));
    wassert(actual(cart.beam_blocking_xy).statsEqual(0, 2883204, 0, 0, 0));
    wassert(actual(cart.elev_fin_xy).statsEqual(0, 2883204, 0, 0, 0));
    wassert(actual(cart.neve_cart).statsEqual(0, 2883204, 0, 0, 0));
    wassert(actual(cart.corr_cart).statsEqual(0, 2883204, 0, 0, 0));
    wassert(actual(cart.conv_cart).statsEqual(0, 2883204, 0, 0, 0));

    CartLowris clow(cb->do_medium ? 512: 256, *cb, cart);
    wassert(actual(clow.ZLR_N_ELEMENTARY_PIXEL) == 4);
    wassert(actual(clow.ZLR_OFFSET) == -175);
    clow.creo_cart_z_lowris();
//    print_stats("clow", clow, cerr);
    wassert(actual(clow.z_out).statsEqual(0, 120384, 1, 45.09, 255));
    wassert(actual(clow.qual_Z_1x1).statsEqual(0, 262144, 0, 0, 0));
    wassert(actual(clow.quota_1x1).statsEqual(0, 262144, 0, 0, 0));
    wassert(actual(clow.dato_corr_1x1).statsEqual(0, 262144, 0, 0, 0));
    wassert(actual(clow.elev_fin_1x1).statsEqual(0, 262144, 0, 0, 0));
    wassert(actual(clow.beam_blocking_1x1).statsEqual(0, 262144, 0, 0, 0));
    wassert(actual(clow.top_1x1).statsEqual(0, 236585, 2, 16.96, 76));
    wassert(actual(clow.neve_1x1).statsEqual(0, 262144, 0, 0, 0));
    wassert(actual(clow.corr_1x1).statsEqual(0, 262144, 0, 0, 0));
    wassert(actual(clow.conv_1x1).statsEqual(0, 262144, 0, 0, 0));
    // TODO: scrivo_out_file_bin

    delete cb;
}

template<> template<>
void to::test<10>()
{
	LOG_CATEGORY("Test");
	LOG_INFO ("Start test 6");

    // versione BB che corrisponde al parametro algo_corto
    static const char* fname = "testdata/2014-05-09-12-40-00.itgat.PVOL.0.h5";
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
    cb->read_odim_volume(fname, 0);
    cb->setup_elaborazione(fname);
    wassert(actual(cb->calcolo_vpr) != (void*)0);

    cb->elabora_dato();

    VolumeStats stats;
    cb->volume.compute_stats(stats);
    wassert(actual(cb->volume).statsEqual( MISSING_DB, 0, -19.69 ,-16.63 , 67.37));

    wassert(actual(cb->beam_blocking).statsEqual(0, 15.83, 51));

    unsigned stats_size = cb->grid_stats.size_az * cb->grid_stats.size_beam;

    ArrayStats<unsigned> stat_anap_stats;
    stat_anap_stats.fill(cb->grid_stats.stat_anap, stats_size);
    //stat_anap_stats.print();
    wassert(actual((unsigned)stat_anap_stats.all_missing).isfalse());
    wassert(actual((unsigned)stat_anap_stats.min) == 0);
    wassert(actual((unsigned)stat_anap_stats.max) == 9);

    ArrayStats<unsigned> stat_anap_tot_stats;
    stat_anap_tot_stats.fill(cb->grid_stats.stat_tot, stats_size);
//    stat_anap_tot_stats.print();
    wassert(actual((unsigned)stat_anap_tot_stats.all_missing).isfalse());
    wassert(actual((unsigned)stat_anap_tot_stats.min) ==  0);
    wassert(actual((unsigned)stat_anap_tot_stats.max) == 1000);

    ArrayStats<unsigned> stat_bloc_stats;
    stat_bloc_stats.fill(cb->grid_stats.stat_bloc, stats_size);
//    stat_bloc_stats.print();
    wassert(actual((unsigned)stat_bloc_stats.all_missing).isfalse());
    wassert(actual((unsigned)stat_bloc_stats.min) == 0);
    wassert(actual((unsigned)stat_bloc_stats.max) == 0);

    ArrayStats<unsigned> stat_elev_stats;
    stat_elev_stats.fill(cb->grid_stats.stat_elev, stats_size);
//    stat_elev_stats.print();
    wassert(actual((unsigned)stat_elev_stats.all_missing).isfalse());
    wassert(actual((unsigned)stat_elev_stats.min) == 0);
    wassert(actual((unsigned)stat_elev_stats.max) == 9);

    wassert(actual(cb->dato_corrotto).statsEqual(0, 0.0, 1));

    cb->caratterizzo_volume();
    wassert(actual(cb->calcolo_vpr) != (void*)0);

//     print_stats("cb", *cb, cerr);
    wassert(actual(*cb->qual).statsEqual(1, 55.1, 99));
    wassert(actual(cb->top).statsEqual(0, 319186, 1, 9.01, 108));
    wassert(actual(cb->first_level).statsEqual(0, 0.33, 3));
    wassert(actual(cb->first_level_static).statsEqual(0, 0.33, 3));
    wassert(actual(cb->bb_first_level).statsEqual(0, 0.41, 2));
    wassert(actual(cb->beam_blocking).statsEqual(0, 15.83, 51));
    wassert(actual(cb->dem).statsEqual(-9999, 47.59, 4070.26));
    wassert(actual(cb->quota).statsEqual(1, 1953.17, 5813));
    wassert(actual(cb->dato_corrotto).statsEqual(0, 0, 1));
    wassert(actual(*cb->calcolo_vpr->flag_vpr).statsEqual(0, 0.89, 1));
    wassert(actual(cb->calcolo_vpr->corr_polar).statsEqual(0, 0, 0));
    wassert(actual(cb->calcolo_vpr->neve).statsEqual(0, 0, 0));

    Cart cart(cb->volume.max_beam_size());
    wassert(actual(cart.max_bin) == 849);
    cart.creo_cart(*cb);
//     print_stats("cart", cart, cerr);
    wassert(actual(cart.cart).statsEqual(0, 624947, 1, 7.23, 255));
    wassert(actual(cart.cartm).statsEqual(0, 2883204, 0, 0, 0));
    wassert(actual(cart.topxy).statsEqual(0, 2800578, 1, 15.8, 108));
    wassert(actual(cart.qual_Z_cart).statsEqual(0, 624947, 1, 19.59, 98));
    wassert(actual(cart.quota_cart).statsEqual(0, 624947, 3, 2618.84, 5813));
    wassert(actual(cart.dato_corr_xy).statsEqual(0, 2883203, 1, 1, 1));
    wassert(actual(cart.beam_blocking_xy).statsEqual(0, 1768227, 1, 49.92, 51));
    wassert(actual(cart.elev_fin_xy).statsEqual(0, 2671545, 1, 1.14, 3));
    wassert(actual(cart.neve_cart).statsEqual(0, 624947, 1, 1, 1));
    wassert(actual(cart.corr_cart).statsEqual(0, 2883204, 0, 0, 0));
    wassert(actual(cart.conv_cart).statsEqual(0, 2883204, 0, 0, 0));

    CartLowris clow(cb->do_medium ? 512: 256, *cb, cart);
    wassert(actual(clow.ZLR_N_ELEMENTARY_PIXEL) == 4);
    wassert(actual(clow.ZLR_OFFSET) == 337);
    clow.creo_cart_z_lowris();
//     print_stats("clow", clow, cerr);
    wassert(actual(clow.z_out).statsEqual(0, 0, 1, 18.9, 255));
    wassert(actual(clow.qual_Z_1x1).statsEqual(0, 0, 1, 32.46, 97));
    wassert(actual(clow.quota_1x1).statsEqual(0, 0, 128, 144.3, 185));
    wassert(actual(clow.dato_corr_1x1).statsEqual(0, 65536, 0, 0, 0));
    wassert(actual(clow.elev_fin_1x1).statsEqual(0, 53780, 1, 1.1, 3));
    wassert(actual(clow.beam_blocking_1x1).statsEqual(0, 42842, 1, 50.14, 51));
    wassert(actual(clow.top_1x1).statsEqual(0, 60002, 1, 13.21, 108));
    wassert(actual(clow.neve_1x1).statsEqual(0, 0, 1, 1, 1));
    wassert(actual(clow.corr_1x1).statsEqual(0, 65536, 0, 0, 0));
    wassert(actual(clow.conv_1x1).statsEqual(0, 65536, 0, 0, 0));

    // TODO: scrivo_out_file_bin

    delete cb;
}

template<> template<>
void to::test<11>()
{
	LOG_CATEGORY("Test");
	LOG_INFO ("Start test 6");

    // versione BB che corrisponde al parametro algo_corto
    static const char* fname = "esplosione/2014-03-01-09-15-00.itgat.PVOL.0.h5";
    unsetwork();
    setenv("DIR_OUT_PP_BLOC", "esplosione", 1);
    setenv("VPR0_FILE"		, "esplosione/vpr_GAT", 1);
    setenv("LAST_VPR"    	, "esplosione/last_vpr_GAT", 1);
    setenv("VPR_HMAX"    	, "esplosione/vpr_hmax_GAT",1); 
    setenv("VPR_HEATING" 	, "esplosione/vpr_heat_GAT",1);
    setenv("LOG_VPR"		, "esplosione/log_VPR_${SITO}",1);
    setenv("TEST_VPR"		, "esplosione/test_vpr",1);
    setenv("FILE_T"		, "esplosione/temperature.txt",1);
    setenv("FIRST_LEVEL_FILE", "../dati/FIRST_LEVEL_corto_GAT_PRI-EST_2011", 1);
    setenv("FILE_ZERO_TERMICO"	, "esplosione/0termico.prev", 1);
	printwork();

    CUM_BAC* cb = new CUM_BAC("GAT",false, 512);
    cb->do_readStaticMap=true;
    cb->do_clean= true;
    cb->do_devel=true;
    cb->do_beamblocking = true;
    cb->do_quality = true;
    cb->do_bloccorr = true;
    cb->do_class = true;
    cb->do_vpr = true;
    cb->read_odim_volume(fname, 0);

    cb->setup_elaborazione(fname);
    cb->elabora_dato();
    cb->caratterizzo_volume();
    wassert(actual(cb->calcolo_vpr) != (void*)0);
    cb->calcolo_vpr->classifica_rain();
    cb->calcolo_vpr->esegui_tutto();
    VolumeStats stats;
    cb->volume.compute_stats(stats);

    delete cb;
}

}
