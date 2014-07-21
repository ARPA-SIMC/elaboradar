#include <wibble/tests.h>
#include "cum_bac.h"
#include "config.h"
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

    Config cfg;
    CUM_BAC* cb = new CUM_BAC(cfg, "GAT");
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
//    stats.print(stdout);
    wassert(actual(stats.count_zeros[0]) == 197201);
    wassert(actual(stats.count_zeros[1]) == 197337);
    wassert(actual(stats.count_zeros[2]) == 197314);
    wassert(actual(stats.count_zeros[3]) == 197459);
    wassert(actual(stats.count_zeros[4]) == 197428);
    wassert(actual(stats.count_zeros[5]) == 196701);
    wassert(actual(stats.count_ones[0]) == 2);
    wassert(actual(stats.count_ones[1]) == 2);
    wassert(actual(stats.count_ones[2]) == 2);
    wassert(actual(stats.count_ones[3]) == 5);
    wassert(actual(stats.count_ones[4]) == 2);
    wassert(actual(stats.count_ones[5]) == 6);
    wassert(actual(stats.count_others[0]) == 397);
    wassert(actual(stats.count_others[1]) == 261);
    wassert(actual(stats.count_others[2]) == 284);
    wassert(actual(stats.count_others[3]) == 136);
    wassert(actual(stats.count_others[4]) == 170);
    wassert(actual(stats.count_others[5]) == 893);
    wassert(actual(stats.sum_others[0]) == 25439);
    wassert(actual(stats.sum_others[1]) == 15892);
    wassert(actual(stats.sum_others[2]) == 22283);
    wassert(actual(stats.sum_others[3]) ==  6522);
    wassert(actual(stats.sum_others[4]) ==  9684);
    wassert(actual(stats.sum_others[5]) == 68041);
    cb->caratterizzo_volume();

//    print_stats("cb", *cb, cerr);
    wassert(actual(*cb->qual).statsEqual(1, 42.34, 99));
    wassert(actual(cb->top).statsEqual(0, 196731, 1, 3.38, 15));
    wassert(actual(cb->first_level).statsEqual(0, 0.86, 3));
    wassert(actual(cb->first_level_static).statsEqual(0, 0.54, 3));
    wassert(actual(cb->bb_first_level).statsEqual(0, 0.49, 2));
    wassert(actual(cb->beam_blocking).statsEqual(0, 3.67, 49));
    wassert(actual(cb->dem).statsEqual(-9999, 47.59, 4070.26));
    wassert(actual(cb->quota).statsEqual(1, 1543.48, 5813));
    wassert(actual(cb->dato_corrotto).statsEqual(0, 2, 3));
    wassert(actual(*cb->calcolo_vpr->flag_vpr).statsEqual(0, 0.69, 1));
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

    Config cfg;
    CUM_BAC* cb = new CUM_BAC(cfg, "GAT");
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
//    stats.print(stdout);
    wassert(actual(stats.count_zeros[0]) == 197201);
    wassert(actual(stats.count_zeros[1]) == 197337);
    wassert(actual(stats.count_zeros[2]) == 197314);
    wassert(actual(stats.count_zeros[3]) == 197459);
    wassert(actual(stats.count_zeros[4]) == 197428);
    wassert(actual(stats.count_zeros[5]) == 196701);
    wassert(actual(stats.count_ones[0]) == 2);
    wassert(actual(stats.count_ones[1]) == 2);
    wassert(actual(stats.count_ones[2]) == 2);
    wassert(actual(stats.count_ones[3]) == 5);
    wassert(actual(stats.count_ones[4]) == 2);
    wassert(actual(stats.count_ones[5]) == 6);
    wassert(actual(stats.count_others[0]) == 397);
    wassert(actual(stats.count_others[1]) == 261);
    wassert(actual(stats.count_others[2]) == 284);
    wassert(actual(stats.count_others[3]) == 136);
    wassert(actual(stats.count_others[4]) == 170);
    wassert(actual(stats.count_others[5]) == 893);
    wassert(actual(stats.sum_others[0]) == 25439);
    wassert(actual(stats.sum_others[1]) == 15892);
    wassert(actual(stats.sum_others[2]) == 22283);
    wassert(actual(stats.sum_others[3]) ==  6522);
    wassert(actual(stats.sum_others[4]) ==  9684);
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

    Config cfg;
    CUM_BAC* cb = new CUM_BAC(cfg, "GAT");
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
    // print_stats("cb", *cb, cerr);
    wassert(actual(*cb->qual).statsEqual(1, 45.07, 99));
    wassert(actual(cb->top).statsEqual(0, 196731, 1, 3.38, 15));
    wassert(actual(cb->first_level).statsEqual(0, 0.54, 3));
    wassert(actual(cb->first_level_static).statsEqual(0, 0.54, 3));
    wassert(actual(cb->bb_first_level).statsEqual(0, 0.49, 2));
    wassert(actual(cb->beam_blocking).statsEqual(0, 4.1, 51));
    wassert(actual(cb->dem).statsEqual(-9999, 47.59, 4070.26));
    wassert(actual(cb->quota).statsEqual(1, 1164.04, 5732));
    wassert(actual(cb->dato_corrotto).statsEqual(0, 1.99, 3));
    wassert(actual(*cb->calcolo_vpr->flag_vpr).statsEqual(0, 0.74, 1));
    wassert(actual(cb->calcolo_vpr->corr_polar).statsEqual(0, 0, 0));
    wassert(actual(cb->calcolo_vpr->neve).statsEqual(0, 0, 0));

    Cart cart(cb->volume.max_beam_size());
    wassert(actual(cart.max_bin) == 494);
    cart.creo_cart(*cb);
    //setenv("DIR_DEBUG", "testdata/", 1);
    //cart.write_out(*cb, cb->assets);
    // print_stats("cart", cart, cout);
    wassert(actual(cart.cart).statsEqual(0, 213127, 1, 1.41, 240));
    wassert(actual(cart.cartm).statsEqual(0, 976144, 0, 0, 0));
    wassert(actual(cart.topxy).statsEqual(0, 974254, 1, 4.79, 15));
    wassert(actual(cart.qual_Z_cart).statsEqual(0, 213127, 1, 1.06, 93));
    wassert(actual(cart.quota_cart).statsEqual(0, 213127, 3, 1519.53, 5732));
    wassert(actual(cart.dato_corr_xy).statsEqual(0, 215809, 2, 2, 3));
    wassert(actual(cart.beam_blocking_xy).statsEqual(0, 973693, 24, 50.71, 51));
    wassert(actual(cart.elev_fin_xy).statsEqual(0, 761197, 1, 1.15, 3));
    wassert(actual(cart.neve_cart).statsEqual(0, 213127, 1, 1, 1));
    wassert(actual(cart.corr_cart).statsEqual(0, 976144, 0, 0, 0));
    wassert(actual(cart.conv_cart).statsEqual(0, 976144, 0, 0, 0));

    CartLowris clow(cb->do_medium ? 512: 256, *cb, cart);
    wassert(actual(clow.CART_DIM_ZLR) == 256);
    wassert(actual(clow.ZLR_N_ELEMENTARY_PIXEL) == 4);
    wassert(actual(clow.ZLR_OFFSET) == -18);
    clow.creo_cart_z_lowris();
    // print_stats("clow", clow, cerr);
    wassert(actual(clow.z_out).statsEqual(0, 17489, 1, 2.23, 240));
    wassert(actual(clow.qual_Z_1x1).statsEqual(0, 17489, 1, 1.22, 92));
    wassert(actual(clow.quota_1x1).statsEqual(0, 0, 128, 138.78, 185));
    wassert(actual(clow.dato_corr_1x1).statsEqual(0, 17972, 2, 2, 3));
    wassert(actual(clow.elev_fin_1x1).statsEqual(0, 52304, 1, 1.14, 3));
    wassert(actual(clow.beam_blocking_1x1).statsEqual(0, 65121, 24, 50.67, 51));
    wassert(actual(clow.top_1x1).statsEqual(0, 65377, 1, 5.2, 15));
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
    Config cfg;
    CUM_BAC* cb = new CUM_BAC(cfg, "GAT");
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
    wassert(actual(cb->top).statsEqual(0, 187460, 1, 7.34, 36));
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
    // print_stats("cart", cart, cerr);
    wassert(actual(cart.cart).statsEqual(0, 213127, 1, 16.54, 242));
    wassert(actual(cart.cartm).statsEqual(0, 976144, 0, 0, 0));
    wassert(actual(cart.topxy).statsEqual(0, 938838, 1, 11.14, 36));
    wassert(actual(cart.qual_Z_cart).statsEqual(0, 213127, 1, 36.51, 98));
    wassert(actual(cart.quota_cart).statsEqual(0, 976144, 0, 0, 0));
    wassert(actual(cart.dato_corr_xy).statsEqual(0, 976144, 0, 0, 0));
    wassert(actual(cart.beam_blocking_xy).statsEqual(0, 689694, 1, 47.03, 51));
    wassert(actual(cart.elev_fin_xy).statsEqual(0, 765448, 1, 1.15, 3));
    wassert(actual(cart.neve_cart).statsEqual(0, 213127, 1, 1, 1));
    wassert(actual(cart.corr_cart).statsEqual(0, 976144, 0, 0, 0));
    wassert(actual(cart.conv_cart).statsEqual(0, 976144, 0, 0, 0));

    CartLowris clow(cb->do_medium ? 512: 256, *cb, cart);
    wassert(actual(clow.ZLR_N_ELEMENTARY_PIXEL) == 4);
    wassert(actual(clow.ZLR_OFFSET) == -18);
    clow.creo_cart_z_lowris();
    // print_stats("clow", clow, cerr);
    wassert(actual(clow.z_out).statsEqual(0, 17489, 1, 25.89, 242));
    wassert(actual(clow.qual_Z_1x1).statsEqual(0, 17489, 1, 35.54, 98));
    wassert(actual(clow.quota_1x1).statsEqual(0, 0, 128, 128, 128));
    wassert(actual(clow.dato_corr_1x1).statsEqual(0, 65536, 0, 0, 0));
    wassert(actual(clow.elev_fin_1x1).statsEqual(0, 54065, 1, 1.13, 3));
    wassert(actual(clow.beam_blocking_1x1).statsEqual(0, 46700, 1, 47.04, 51));
    wassert(actual(clow.top_1x1).statsEqual(0, 62214, 1, 11.55, 36));
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
    Config cfg;
    CUM_BAC* cb = new CUM_BAC(cfg, "GAT");
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

//         print_stats("cb", *cb, cerr);
    wassert(actual(*cb->qual).statsEqual(1, 53.07, 99));
    wassert(actual(cb->top).statsEqual(0, 196731, 1, 3.38, 15));
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
    // print_stats("cart", cart, cerr);
    wassert(actual(cart.cart).statsEqual(0, 213127, 1, 1.31, 229));
    wassert(actual(cart.cartm).statsEqual(0, 976144, 0, 0, 0));
    wassert(actual(cart.topxy).statsEqual(0, 974260, 1, 4.77, 15));
    wassert(actual(cart.qual_Z_cart).statsEqual(0, 213127, 1, 36.6, 98));
    wassert(actual(cart.quota_cart).statsEqual(0, 976144, 0, 0, 0));
    wassert(actual(cart.dato_corr_xy).statsEqual(0, 976144, 0, 0, 0));
    wassert(actual(cart.beam_blocking_xy).statsEqual(0, 691160, 1, 47.08, 51));
    wassert(actual(cart.elev_fin_xy).statsEqual(0, 761192, 1, 1.15, 3));
    wassert(actual(cart.neve_cart).statsEqual(0, 213127, 1, 1, 1));
    wassert(actual(cart.corr_cart).statsEqual(0, 976144, 0, 0, 0));
    wassert(actual(cart.conv_cart).statsEqual(0, 976144, 0, 0, 0));

    CartLowris clow(cb->do_medium ? 512: 256, *cb, cart);
    wassert(actual(clow.ZLR_N_ELEMENTARY_PIXEL) == 4);
    wassert(actual(clow.ZLR_OFFSET) == -18);
    clow.creo_cart_z_lowris();
    // print_stats("clow", clow, cerr);
    wassert(actual(clow.z_out).statsEqual(0, 17489, 1, 1.99, 229));
    wassert(actual(clow.qual_Z_1x1).statsEqual(0, 17489, 1, 36.2, 98));
    wassert(actual(clow.quota_1x1).statsEqual(0, 0, 128, 128, 128));
    wassert(actual(clow.dato_corr_1x1).statsEqual(0, 65536, 0, 0, 0));
    wassert(actual(clow.elev_fin_1x1).statsEqual(0, 52302, 1, 1.15, 3));
    wassert(actual(clow.beam_blocking_1x1).statsEqual(0, 47319, 1, 47.16, 51));
    wassert(actual(clow.top_1x1).statsEqual(0, 65380, 1, 5.09, 15));
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

    Config cfg;
    CUM_BAC* cb = new CUM_BAC(cfg, "GAT");
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
    wassert(actual(stats.count_zeros[0]) == 170006);
    wassert(actual(stats.count_zeros[1]) == 209229);
    wassert(actual(stats.count_zeros[2]) == 161479);
    wassert(actual(stats.count_zeros[3]) == 108897);
    wassert(actual(stats.count_zeros[4]) == 118263);
    wassert(actual(stats.count_ones[0]) == 0);
    wassert(actual(stats.count_ones[1]) == 0);
    wassert(actual(stats.count_ones[2]) == 0);
    wassert(actual(stats.count_ones[3]) == 0);
    wassert(actual(stats.count_ones[4]) == 0);
    wassert(actual(stats.count_others[0]) == 169594);
    wassert(actual(stats.count_others[1]) == 130371);
    wassert(actual(stats.count_others[2]) == 102521);
    wassert(actual(stats.count_others[3]) ==  88703);
    wassert(actual(stats.count_others[4]) ==  79337);
    wassert(actual(stats.sum_others[0]) == 22253126);
    wassert(actual(stats.sum_others[1]) == 16789233);
    wassert(actual(stats.sum_others[2]) == 12888472);
    wassert(actual(stats.sum_others[3]) == 10805989);
    wassert(actual(stats.sum_others[4]) ==  9340826);

    wassert(actual(cb->beam_blocking).statsEqual(0, 5.09, 51));

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
    wassert(actual((unsigned)dato_corrotto_stats.max) == 3);


    cb->caratterizzo_volume();
    wassert(actual(cb->calcolo_vpr) != (void*)0);
//    print_stats("cb", *cb, cerr);
    wassert(actual(*cb->qual).statsEqual(1, 51.9, 99));
    wassert(actual(cb->top).statsEqual(0, 222889, 1, 13.23, 76));
    wassert(actual(cb->first_level).statsEqual(0, 0.33, 3));
    wassert(actual(cb->first_level_static).statsEqual(0, 0.33, 3));
    wassert(actual(cb->bb_first_level).statsEqual(0, 0.48, 2));
    wassert(actual(cb->beam_blocking).statsEqual(0, 5.09, 51));
    wassert(actual(cb->dem).statsEqual(-9999, 47.59, 4070.26));
    wassert(actual(cb->quota).statsEqual(1, 2027.25, 5825));
    wassert(actual(cb->dato_corrotto).statsEqual(0, 1.03, 3));
    wassert(actual(*cb->calcolo_vpr->flag_vpr).statsEqual(0, 0.85, 1));
    wassert(actual(cb->calcolo_vpr->corr_polar).statsEqual(0, 0, 0));
    wassert(actual(cb->calcolo_vpr->neve).statsEqual(0, 0, 0));

    Cart cart(cb->volume.max_beam_size());
    wassert(actual(cart.max_bin) == 849);
    cart.creo_cart(*cb);
    // print_stats("cart", cart, cerr);
    wassert(actual(cart.cart).statsEqual(0, 624947, 1, 43.39, 255));
    wassert(actual(cart.cartm).statsEqual(0, 2883204, 0, 0, 0));
    wassert(actual(cart.topxy).statsEqual(0, 2516298, 2, 16.45, 76));
    wassert(actual(cart.qual_Z_cart).statsEqual(0, 624947, 1, 11.64, 98));
    wassert(actual(cart.quota_cart).statsEqual(0, 624947, 3, 2716.01, 5825));
    wassert(actual(cart.dato_corr_xy).statsEqual(0, 1369711, 2, 2.02, 3));
    wassert(actual(cart.beam_blocking_xy).statsEqual(0, 2647084, 1, 46.19, 51));
    wassert(actual(cart.elev_fin_xy).statsEqual(0, 2683727, 1, 1.15, 3));
    wassert(actual(cart.neve_cart).statsEqual(0, 624947, 1, 1, 1));
    wassert(actual(cart.corr_cart).statsEqual(0, 2883204, 0, 0, 0));
    wassert(actual(cart.conv_cart).statsEqual(0, 2883204, 0, 0, 0));

    CartLowris clow(cb->do_medium ? 512: 256, *cb, cart);
    wassert(actual(clow.ZLR_N_ELEMENTARY_PIXEL) == 4);
    wassert(actual(clow.ZLR_OFFSET) == 337);
    clow.creo_cart_z_lowris();
    // print_stats("clow", clow, cerr);
    wassert(actual(clow.z_out).statsEqual(0, 0, 1, 82.8, 255));
    wassert(actual(clow.qual_Z_1x1).statsEqual(0, 0, 1, 22.61, 97));
    wassert(actual(clow.quota_1x1).statsEqual(0, 0, 128, 144.45, 186));
    wassert(actual(clow.dato_corr_1x1).statsEqual(0, 39084, 2, 2.06, 3));
    wassert(actual(clow.elev_fin_1x1).statsEqual(0, 55923, 1, 1.16, 3));
    wassert(actual(clow.beam_blocking_1x1).statsEqual(0, 53233, 1, 44.89, 51));
    wassert(actual(clow.top_1x1).statsEqual(0, 42124, 2, 14.01, 60));
    wassert(actual(clow.neve_1x1).statsEqual(0, 0, 1, 1, 1));
    wassert(actual(clow.corr_1x1).statsEqual(0, 65536, 0, 0, 0));
    wassert(actual(clow.conv_1x1).statsEqual(0, 65536, 0, 0, 0));

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

    Config cfg;
    CUM_BAC* cb = new CUM_BAC(cfg, "GAT");
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

//    print_stats("cb", *cb, cerr);
    wassert(actual(*cb->qual).statsEqual(1, 52.92, 99));
    wassert(actual(cb->top).statsEqual(0, 222889, 1, 13.23, 76));
    wassert(actual(cb->first_level).statsEqual(0, 0.33, 3));
    wassert(actual(cb->first_level_static).statsEqual(0, 0.33, 3));
    wassert(actual(cb->bb_first_level).statsEqual(0, 0.48, 2));
    wassert(actual(cb->beam_blocking).statsEqual(0, 17.1, 51));
    wassert(actual(cb->dem).statsEqual(-9999, 47.59, 4070.26));
    wassert(actual(cb->quota).statsEqual(0, 0, 0));
    wassert(actual(cb->dato_corrotto).statsEqual(0, 0, 0));
    wassert(actual(*cb->calcolo_vpr->flag_vpr).statsEqual(0, 0.73, 1));
    wassert(actual(cb->calcolo_vpr->corr_polar).statsEqual(0, 0, 0));
    wassert(actual(cb->calcolo_vpr->neve).statsEqual(0, 0, 0));
    // TODO: cb->stampa_vpr()

    Cart cart(cb->volume.max_beam_size());
    wassert(actual(cart.max_bin) == 849);
    cart.creo_cart(*cb);
    // print_stats("cart", cart, cerr);
    wassert(actual(cart.cart).statsEqual(0, 624947, 1, 42.31, 255));
    wassert(actual(cart.cartm).statsEqual(0, 2883204, 0, 0, 0));
    wassert(actual(cart.topxy).statsEqual(0, 2516370, 2, 16.42, 76));
    wassert(actual(cart.qual_Z_cart).statsEqual(0, 624947, 1, 19.13, 98));
    wassert(actual(cart.quota_cart).statsEqual(0, 2883204, 0, 0, 0));
    wassert(actual(cart.dato_corr_xy).statsEqual(0, 2883204, 0, 0, 0));
    wassert(actual(cart.beam_blocking_xy).statsEqual(0, 1687528, 1, 45.08, 51));
    wassert(actual(cart.elev_fin_xy).statsEqual(0, 2683498, 1, 1.15, 3));
    wassert(actual(cart.neve_cart).statsEqual(0, 624947, 1, 1, 1));
    wassert(actual(cart.corr_cart).statsEqual(0, 2883204, 0, 0, 0));
    wassert(actual(cart.conv_cart).statsEqual(0, 2883204, 0, 0, 0));

    CartLowris clow(cb->do_medium ? 512: 256, *cb, cart);
    wassert(actual(clow.ZLR_N_ELEMENTARY_PIXEL) == 4);
    wassert(actual(clow.ZLR_OFFSET) == 337);
    clow.creo_cart_z_lowris();
    // print_stats("clow", clow, cerr);
    wassert(actual(clow.z_out).statsEqual(0, 0, 1, 81.14, 255));
    wassert(actual(clow.qual_Z_1x1).statsEqual(0, 0, 1, 28.47, 97));
    wassert(actual(clow.quota_1x1).statsEqual(0, 0, 128, 128, 128));
    wassert(actual(clow.dato_corr_1x1).statsEqual(0, 65536, 0, 0, 0));
    wassert(actual(clow.elev_fin_1x1).statsEqual(0, 55892, 1, 1.16, 3));
    wassert(actual(clow.beam_blocking_1x1).statsEqual(0, 37285, 1, 45.33, 51));
    wassert(actual(clow.top_1x1).statsEqual(0, 42190, 2, 13.88, 59));
    wassert(actual(clow.neve_1x1).statsEqual(0, 0, 1, 1, 1));
    wassert(actual(clow.corr_1x1).statsEqual(0, 65536, 0, 0, 0));
    wassert(actual(clow.conv_1x1).statsEqual(0, 65536, 0, 0, 0));

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

    Config cfg;
    CUM_BAC* cb = new CUM_BAC(cfg, "GAT");
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
//    print_stats("cb", *cb, cerr);
    wassert(actual(*cb->qual).statsEqual(1, 52.92, 99));
    wassert(actual(cb->top).statsEqual(0, 222889, 1, 13.23, 76));
    wassert(actual(cb->first_level).statsEqual(0, 0.33, 3));
    wassert(actual(cb->first_level_static).statsEqual(0, 0.33, 3));
    wassert(actual(cb->bb_first_level).statsEqual(0, 0.48, 2));
    wassert(actual(cb->beam_blocking).statsEqual(0, 17.1, 51));
    wassert(actual(cb->dem).statsEqual(-9999, 47.59, 4070.26));
    wassert(actual(cb->quota).statsEqual(0, 0, 0));
    wassert(actual(cb->dato_corrotto).statsEqual(0, 0, 0));
    wassert(actual(*cb->calcolo_vpr->flag_vpr).statsEqual(0, 0.73, 1));
    wassert(actual(cb->calcolo_vpr->corr_polar).statsEqual(0, 0, 0));
    wassert(actual(cb->calcolo_vpr->neve).statsEqual(0, 0, 0));

    Cart cart(cb->volume.max_beam_size());
    wassert(actual(cart.max_bin) == 849);
    cart.creo_cart(*cb);
    // print_stats("cart", cart, cerr);
    wassert(actual(cart.cart).statsEqual(0, 624947, 1, 42.3, 255));
    wassert(actual(cart.cartm).statsEqual(0, 2883204, 0, 0, 0));
    wassert(actual(cart.topxy).statsEqual(0, 2516370, 2, 16.42, 76));
    wassert(actual(cart.qual_Z_cart).statsEqual(0, 624947, 1, 19.13, 98));
    wassert(actual(cart.quota_cart).statsEqual(0, 2883204, 0, 0, 0));
    wassert(actual(cart.dato_corr_xy).statsEqual(0, 2883204, 0, 0, 0));
    wassert(actual(cart.beam_blocking_xy).statsEqual(0, 1687528, 1, 45.08, 51));
    wassert(actual(cart.elev_fin_xy).statsEqual(0, 2683498, 1, 1.15, 3));
    wassert(actual(cart.neve_cart).statsEqual(0, 624947, 1, 1, 1));
    wassert(actual(cart.corr_cart).statsEqual(0, 2883204, 0, 0, 0));
    wassert(actual(cart.conv_cart).statsEqual(0, 2882639, 100, 100, 100));

    CartLowris clow(cb->do_medium ? 512: 256, *cb, cart);
    wassert(actual(clow.ZLR_N_ELEMENTARY_PIXEL) == 4);
    wassert(actual(clow.ZLR_OFFSET) == 337);
    clow.creo_cart_z_lowris();
    // print_stats("clow", clow, cerr);
    wassert(actual(clow.z_out).statsEqual(0, 0, 1, 81.14, 255));
    wassert(actual(clow.qual_Z_1x1).statsEqual(0, 0, 1, 28.47, 97));
    wassert(actual(clow.quota_1x1).statsEqual(0, 0, 128, 128, 128));
    wassert(actual(clow.dato_corr_1x1).statsEqual(0, 65536, 0, 0, 0));
    wassert(actual(clow.elev_fin_1x1).statsEqual(0, 55892, 1, 1.16, 3));
    wassert(actual(clow.beam_blocking_1x1).statsEqual(0, 37285, 1, 45.33, 51));
    wassert(actual(clow.top_1x1).statsEqual(0, 42190, 2, 13.88, 59));
    wassert(actual(clow.neve_1x1).statsEqual(0, 0, 1, 1, 1));
    wassert(actual(clow.corr_1x1).statsEqual(0, 65536, 0, 0, 0));
    wassert(actual(clow.conv_1x1).statsEqual(0, 65529, 100, 100, 100));

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

    Config cfg;
    CUM_BAC* cb = new CUM_BAC(cfg, "GAT",true,1024);
    cb->do_medium= true;
    cb->do_readStaticMap=true;
//    cb->do_zlr_media=true; 
    cb->read_sp20_volume(fname, 0);
    cb->setup_elaborazione(fname);
    cb->assets.write_vpr_heating(0);

    cb->elabora_dato();
    VolumeStats stats;
    cb->volume.compute_stats(stats);
    //stats.print(stdout);
    wassert(actual(stats.count_zeros[0]) == 172263);
    wassert(actual(stats.count_zeros[1]) == 207286);
    wassert(actual(stats.count_zeros[2]) == 158971);
    wassert(actual(stats.count_zeros[3]) == 108187);
    wassert(actual(stats.count_zeros[4]) == 117678);
    wassert(actual(stats.count_ones[0]) == 4928);
    wassert(actual(stats.count_ones[1]) == 4928);
    wassert(actual(stats.count_ones[2]) == 3012);
    wassert(actual(stats.count_ones[3]) ==    0);
    wassert(actual(stats.count_ones[4]) ==    0);
    wassert(actual(stats.count_others[0]) == 162409);
    wassert(actual(stats.count_others[1]) == 127386);
    wassert(actual(stats.count_others[2]) == 102017);
    wassert(actual(stats.count_others[3]) ==  89413);
    wassert(actual(stats.count_others[4]) ==  79922);
    wassert(actual(stats.sum_others[0]) == 21242356);
    wassert(actual(stats.sum_others[1]) == 16379035);
    wassert(actual(stats.sum_others[2]) == 12831241);
    wassert(actual(stats.sum_others[3]) == 10859368);
    wassert(actual(stats.sum_others[4]) ==  9385431);
   // print_stats("cb", *cb, cerr);

    Cart cart(cb->volume.max_beam_size());
    wassert(actual(cart.max_bin) == 849);
    cart.creo_cart(*cb);
    // print_stats("cart", cart, cerr);
    wassert(actual(cart.cart).statsEqual(0, 624947, 1, 39.72, 255));
    wassert(actual(cart.cartm).statsEqual(0, 2883204, 0, 0, 0));
    wassert(actual(cart.topxy).statsEqual(0, 2522699, 2, 16.34, 76));
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
    // print_stats("clow", clow, cerr);
    wassert(actual(clow.z_out).statsEqual(0, 120384, 1, 49.81, 255));
    wassert(actual(clow.qual_Z_1x1).statsEqual(0, 262144, 0, 0, 0));
    wassert(actual(clow.quota_1x1).statsEqual(0, 262144, 0, 0, 0));
    wassert(actual(clow.dato_corr_1x1).statsEqual(0, 262144, 0, 0, 0));
    wassert(actual(clow.elev_fin_1x1).statsEqual(0, 262144, 0, 0, 0));
    wassert(actual(clow.beam_blocking_1x1).statsEqual(0, 262144, 0, 0, 0));
    wassert(actual(clow.top_1x1).statsEqual(0, 234005, 2, 18.09, 76));
    wassert(actual(clow.neve_1x1).statsEqual(0, 262144, 0, 0, 0));
    wassert(actual(clow.corr_1x1).statsEqual(0, 262144, 0, 0, 0));
    wassert(actual(clow.conv_1x1).statsEqual(0, 262144, 0, 0, 0));

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

    Config cfg;
    CUM_BAC* cb = new CUM_BAC(cfg, "GAT");
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
    wassert(actual(cb->volume).statsEqual( MISSING_DB, 0, -31.5,-25.86 , 71.25));

    wassert(actual(cb->beam_blocking).statsEqual(0, 1.25, 51));

    unsigned stats_size = cb->grid_stats.size_az * cb->grid_stats.size_beam;

    ArrayStats<unsigned> stat_anap_stats;
    stat_anap_stats.fill(cb->grid_stats.stat_anap, stats_size);
    //stat_anap_stats.print();
    wassert(actual((unsigned)stat_anap_stats.all_missing).isfalse());
    wassert(actual((unsigned)stat_anap_stats.min) == 0);
    wassert(actual((unsigned)stat_anap_stats.max) == 1);

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
    wassert(actual((unsigned)stat_elev_stats.max) == 1);

    wassert(actual(cb->dato_corrotto).statsEqual(0, 1.6, 3));

    cb->caratterizzo_volume();
    wassert(actual(cb->calcolo_vpr) != (void*)0);

//     print_stats("cb", *cb, cerr);
    wassert(actual(*cb->qual).statsEqual(1, 52.29, 99));
    wassert(actual(cb->top).statsEqual(0, 315507, 1, 9.38, 108));
    wassert(actual(cb->first_level).statsEqual(0, 0.33, 3));
    wassert(actual(cb->first_level_static).statsEqual(0, 0.33, 3));
    wassert(actual(cb->bb_first_level).statsEqual(0, 0.41, 2));
    wassert(actual(cb->beam_blocking).statsEqual(0, 1.25, 51));
    wassert(actual(cb->dem).statsEqual(-9999, 47.59, 4070.26));
    wassert(actual(cb->quota).statsEqual(1, 1953.16, 5813));
    wassert(actual(cb->dato_corrotto).statsEqual(0, 1.6, 3));
    wassert(actual(*cb->calcolo_vpr->flag_vpr).statsEqual(0, 0.82, 1));
    wassert(actual(cb->calcolo_vpr->corr_polar).statsEqual(0, 0, 0));
    wassert(actual(cb->calcolo_vpr->neve).statsEqual(0, 0, 0));

    Cart cart(cb->volume.max_beam_size());
    wassert(actual(cart.max_bin) == 849);
    cart.creo_cart(*cb);
    // print_stats("cart", cart, cerr);
    wassert(actual(cart.cart).statsEqual(0, 624947, 1, 9.39, 255));
    wassert(actual(cart.cartm).statsEqual(0, 2883204, 0, 0, 0));
    wassert(actual(cart.topxy).statsEqual(0, 2783442, 1, 16.05, 108));
    wassert(actual(cart.qual_Z_cart).statsEqual(0, 624947, 1, 4.63, 98));
    wassert(actual(cart.quota_cart).statsEqual(0, 624947, 3, 2618.04, 5813));
    wassert(actual(cart.dato_corr_xy).statsEqual(0, 786700, 2, 2.01, 3));
    wassert(actual(cart.beam_blocking_xy).statsEqual(0, 2834226, 1, 48.6, 51));
    wassert(actual(cart.elev_fin_xy).statsEqual(0, 2672926, 1, 1.14, 3));
    wassert(actual(cart.neve_cart).statsEqual(0, 624947, 1, 1, 1));
    wassert(actual(cart.corr_cart).statsEqual(0, 2883204, 0, 0, 0));
    wassert(actual(cart.conv_cart).statsEqual(0, 2883204, 0, 0, 0));

    CartLowris clow(cb->do_medium ? 512: 256, *cb, cart);
    wassert(actual(clow.ZLR_N_ELEMENTARY_PIXEL) == 4);
    wassert(actual(clow.ZLR_OFFSET) == 337);
    clow.creo_cart_z_lowris();
    // print_stats("clow", clow, cerr);
    wassert(actual(clow.z_out).statsEqual(0, 0, 1, 24.26, 255));
    wassert(actual(clow.qual_Z_1x1).statsEqual(0, 0, 1, 11.26, 97));
    wassert(actual(clow.quota_1x1).statsEqual(0, 0, 128, 144.25, 185));
    wassert(actual(clow.dato_corr_1x1).statsEqual(0, 11862, 2, 2.03, 3));
    wassert(actual(clow.elev_fin_1x1).statsEqual(0, 54038, 1, 1.09, 3));
    wassert(actual(clow.beam_blocking_1x1).statsEqual(0, 63117, 1, 46.77, 51));
    wassert(actual(clow.top_1x1).statsEqual(0, 58889, 1, 13.52, 108));
    wassert(actual(clow.neve_1x1).statsEqual(0, 0, 1, 1, 1));
    wassert(actual(clow.corr_1x1).statsEqual(0, 65536, 0, 0, 0));
    wassert(actual(clow.conv_1x1).statsEqual(0, 65536, 0, 0, 0));

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
    setenv("VPR_ARCH"           , "esplosione/201403010915_vpr_GAT",1);
	printwork();

    Config cfg;
    CUM_BAC* cb = new CUM_BAC(cfg, "GAT",false, 512);
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
