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

    wassert(actual(*cb->qual).statsEqual(0, 55.06, 99));
    wassert(actual(*cb->calcolo_vpr->flag_vpr).statsEqual(0, 0.88, 1));
    wassert(actual(cb->top).statsEqual(0, 0.01, 15));

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
    wassert(actual(*cb->qual).statsEqual(0, 55.07, 99));
    wassert(actual(*cb->calcolo_vpr->flag_vpr).statsEqual(0, 0.88, 1));
    //wassert(actual(*cb->calcolo_vpr->flag_vpr).statsEqual(0, 143635, 0, 5.27, 10));
    wassert(actual(cb->top).statsEqual(0, 0.01, 15));
    //wassert(actual(cb->top).statsEqual(0, 204765, 1, 0, 15));

    Cart cart(cb->volume.max_beam_size());
    cart.creo_cart(*cb);
    //print_stats("cart", cart, cerr);
    wassert(actual(cart.cart).statsEqual(0, 0.94, 227));
    wassert(actual(cart.cartm).statsEqual(0, 0, 0));
    wassert(actual(cart.topxy).statsEqual(0, 0.01, 15));
    wassert(actual(cart.qual_Z_cart).statsEqual(0, 31.78, 98));
    wassert(actual(cart.quota_cart).statsEqual(0, 1187.93, 5732));
    wassert(actual(cart.dato_corr_xy).statsEqual(0, 0, 1));
    wassert(actual(cart.beam_blocking_xy).statsEqual(0, 8.9, 51));
    wassert(actual(cart.elev_fin_xy).statsEqual(0, 0.25, 4));
    wassert(actual(cart.neve_cart).statsEqual(0, 0.78, 1));
    wassert(actual(cart.corr_cart).statsEqual(0, 0, 0));
    wassert(actual(cart.conv_cart).statsEqual(0, 0, 0));

    CartLowris clow(cb->do_medium ? 512: 256);
    clow.creo_cart_z_lowris(*cb, cart);
    //print_stats("clow", clow, cerr);
    wassert(actual(clow.z_out).statsEqual(0, 0.89, 227));
    wassert(actual(clow.qual_Z_1x1).statsEqual(0, 29.59, 97));
    wassert(actual(clow.quota_1x1).statsEqual(128, 138.7, 185));
    wassert(actual(clow.dato_corr_1x1).statsEqual(0, 0, 1));
    wassert(actual(clow.elev_fin_1x1).statsEqual(0, 0.24, 3));
    wassert(actual(clow.beam_blocking_1x1).statsEqual(0, 8.28, 51));
    wassert(actual(clow.top_1x1).statsEqual(0, 0.01, 15));
    wassert(actual(clow.neve_1x1).statsEqual(0, 0, 1));
    wassert(actual(clow.corr_1x1).statsEqual(0, 0, 0));
    wassert(actual(clow.conv_1x1).statsEqual(0, 0, 0));

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
    wassert(actual(*cb->qual).statsEqual(0, 55.09, 99));
    wassert(actual(*cb->calcolo_vpr->flag_vpr).statsEqual(/*0, 143320,*/ 0, 0.88, 1));
    wassert(actual(cb->top).statsEqual(/*203211, 4,*/ 0, 0.3, 36));

    // la combina_profili restituisce 1 se non riesce a costruire un profilo
    // perchè non piove o piove poco
    int ier = cb->calcolo_vpr->combina_profili();
    wassert(actual(ier) == 1);

    cb->calcolo_vpr->heating = cb->calcolo_vpr->profile_heating();
    wassert(actual(cb->calcolo_vpr->heating) == 0);

    //ier = cb->calcolo_vpr->corr_vpr();
    //wassert(actual(ier) == 1);

    // TODO: cb->stampa_vpr()

    Cart cart(cb->volume.max_beam_size());
    cart.creo_cart(*cb);
    //print_stats("cart", cart, cerr);
    wassert(actual(cart.cart).statsEqual(0, 10.68, 227));
    wassert(actual(cart.cartm).statsEqual(0, 0.00, 0));
    wassert(actual(cart.topxy).statsEqual(0, 0.34, 36));
    wassert(actual(cart.qual_Z_cart).statsEqual(0, 31.74, 98));
    wassert(actual(cart.quota_cart).statsEqual(0, 0.00, 0));
    wassert(actual(cart.dato_corr_xy).statsEqual(0, 0.00, 0));
    wassert(actual(cart.beam_blocking_xy).statsEqual(0, 8.92, 51));
    wassert(actual(cart.elev_fin_xy).statsEqual(0, 0.25, 3));
    wassert(actual(cart.neve_cart).statsEqual(0, 0.78, 1));
    wassert(actual(cart.corr_cart).statsEqual(0, 0.00, 0));
    wassert(actual(cart.conv_cart).statsEqual(0, 0.00, 0));

    CartLowris clow(cb->do_medium ? 512: 256);
    clow.creo_cart_z_lowris(*cb, cart);
    //print_stats("clow", clow, cerr);
    wassert(actual(clow.z_out).statsEqual(0, 16, 227));
    wassert(actual(clow.qual_Z_1x1).statsEqual(0, 29.55, 98));
    wassert(actual(clow.quota_1x1).statsEqual(128, 128.00, 128));
    wassert(actual(clow.dato_corr_1x1).statsEqual(0, 0.00, 0));
    wassert(actual(clow.elev_fin_1x1).statsEqual(0, 0.23, 3));
    wassert(actual(clow.beam_blocking_1x1).statsEqual(0, 8.30, 51));
    wassert(actual(clow.top_1x1).statsEqual(0, 0.32, 36));
    wassert(actual(clow.neve_1x1).statsEqual(0, 0.00, 1));
    wassert(actual(clow.corr_1x1).statsEqual(0, 0.00, 0));
    wassert(actual(clow.conv_1x1).statsEqual(0, 0.00, 0));

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
    wassert(actual(*cb->qual).statsEqual(0, 55.09, 99));
    wassert(actual(*cb->calcolo_vpr->flag_vpr).statsEqual(/*0, 143320,*/ 0, 0.88, 1));
    wassert(actual(cb->top).statsEqual(/*0, 204764,*/ 0, 0.01, 15));

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

    Cart cart(cb->volume.max_beam_size());
    cart.creo_cart(*cb);
    //print_stats("cart", cart, cerr);
    wassert(actual(cart.cart).statsEqual(0, 0.96, 227));
    wassert(actual(cart.cartm).statsEqual(0, 0.00, 0));
    wassert(actual(cart.topxy).statsEqual(0, 0.01, 15));
    wassert(actual(cart.qual_Z_cart).statsEqual(0, 31.77, 98));
    wassert(actual(cart.quota_cart).statsEqual(0, 0.00, 0));
    wassert(actual(cart.dato_corr_xy).statsEqual(0, 0.00, 0));
    wassert(actual(cart.beam_blocking_xy).statsEqual(0, 8.91, 51));
    wassert(actual(cart.elev_fin_xy).statsEqual(0, 0.25, 3));
    wassert(actual(cart.neve_cart).statsEqual(0, 0.78, 1));
    wassert(actual(cart.corr_cart).statsEqual(0, 0.00, 0));
    wassert(actual(cart.conv_cart).statsEqual(0, 0.00, 0));

    CartLowris clow(cb->do_medium ? 512: 256);
    clow.creo_cart_z_lowris(*cb, cart);
    //print_stats("clow", clow, cerr);
    wassert(actual(clow.z_out).statsEqual(0, 0.91, 227));
    wassert(actual(clow.qual_Z_1x1).statsEqual(0, 29.59, 97));
    wassert(actual(clow.quota_1x1).statsEqual(128, 128.00, 128));
    wassert(actual(clow.dato_corr_1x1).statsEqual(0, 0.00, 0));
    wassert(actual(clow.elev_fin_1x1).statsEqual(0, 0.24, 3));
    wassert(actual(clow.beam_blocking_1x1).statsEqual(0, 8.29, 51));
    wassert(actual(clow.top_1x1).statsEqual(0, 0.01, 15));
    wassert(actual(clow.neve_1x1).statsEqual(0, 0.00, 1));
    wassert(actual(clow.corr_1x1).statsEqual(0, 0.00, 0));
    wassert(actual(clow.conv_1x1).statsEqual(0, 0.00, 0));

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

    wassert(actual(cb->beam_blocking).statsEqual(0, 13.49, 51));

    unsigned stats_size = cb->grid_stats.size_az * cb->grid_stats.size_beam;

    ArrayStats<unsigned> stat_anap_stats;
    stat_anap_stats.fill(cb->grid_stats.stat_anap, stats_size);
//    stat_anap_stats.print();
    wassert(actual((unsigned)stat_anap_stats.all_missing).isfalse());
    wassert(actual((unsigned)stat_anap_stats.min) == 0);
    wassert(actual((unsigned)stat_anap_stats.max) == 75);

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
    wassert(actual((unsigned)stat_elev_stats.max) == 75);

    ArrayStats<unsigned char> dato_corrotto_stats;
    dato_corrotto_stats.fill(cb->dato_corrotto);
//    dato_corrotto_stats.print();
    wassert(actual((unsigned)dato_corrotto_stats.all_missing).isfalse());
    wassert(actual((unsigned)dato_corrotto_stats.min) == 0);
    wassert(actual((unsigned)dato_corrotto_stats.max) == 1);


    cb->caratterizzo_volume();
    wassert(actual(cb->calcolo_vpr) != (void*)0);
    //print_stats("cb", *cb, cerr);

    wassert(actual(*cb->qual).statsEqual(0, 58.80, 99));
    //wassert(actual((unsigned)stats_qual.count_zeros) == 0);
    //wassert(actual((unsigned)stats_qual.count_ones) == 162599);

    wassert(actual(*cb->calcolo_vpr->flag_vpr).statsEqual(0, 0.91, 1));
    //wassert(actual((unsigned)stats_flag_vpr.count_zeros) == 167079);
    //wassert(actual((unsigned)stats_flag_vpr.count_ones) == 1650521);

    wassert(actual(cb->top).statsEqual(0, 6.29, 76));
    //wassert(actual((unsigned)stats_top.count_zeros) == 112778);
    //wassert(actual((unsigned)stats_top.count_ones) == 159);

    Cart cart(cb->volume.max_beam_size());
    cart.creo_cart(*cb);
    //print_stats("cart", cart, cerr);
    wassert(actual(cart.cart).statsEqual(0, 54.82, 255));
    wassert(actual(cart.cartm).statsEqual(0, 0.00, 0));
    wassert(actual(cart.topxy).statsEqual(0, 3.23, 76));
    wassert(actual(cart.qual_Z_cart).statsEqual(0, 25.66, 98));
    wassert(actual(cart.quota_cart).statsEqual(0, 1201.10, 5825));
    wassert(actual(cart.dato_corr_xy).statsEqual(0, 0.00, 1));
    wassert(actual(cart.beam_blocking_xy).statsEqual(0, 14.12, 51));
    wassert(actual(cart.elev_fin_xy).statsEqual(0, 0.22, 3));
    wassert(actual(cart.neve_cart).statsEqual(0, 0.78, 1));
    wassert(actual(cart.corr_cart).statsEqual(0, 0.00, 0));
    wassert(actual(cart.conv_cart).statsEqual(0, 0.00, 0));

    CartLowris clow(cb->do_medium ? 512: 256);
    clow.creo_cart_z_lowris(*cb, cart);
    //print_stats("clow", clow, cerr);
    wassert(actual(clow.z_out).statsEqual(0, 65.91, 255));
    wassert(actual(clow.qual_Z_1x1).statsEqual(0, 25.65, 98));
    wassert(actual(clow.quota_1x1).statsEqual(128, 139.64, 186));
    wassert(actual(clow.dato_corr_1x1).statsEqual(0, 0.00, 1));
    wassert(actual(clow.elev_fin_1x1).statsEqual(0, 0.22, 3));
    wassert(actual(clow.beam_blocking_1x1).statsEqual(0, 14.12, 51));
    wassert(actual(clow.top_1x1).statsEqual(0, 3.23, 42));
    wassert(actual(clow.neve_1x1).statsEqual(0, 0.00, 1));
    wassert(actual(clow.corr_1x1).statsEqual(0, 0.00, 0));
    wassert(actual(clow.conv_1x1).statsEqual(0, 0.00, 0));

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

    Cart cart(cb->volume.max_beam_size());
    cart.creo_cart(*cb);
    //print_stats("cart", cart, cerr);
    CartLowris clow(cb->do_medium ? 512: 256);
    clow.creo_cart_z_lowris(*cb, cart);
    //print_stats("clow", clow, cerr);

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

    Cart cart(cb->volume.max_beam_size());
    cart.creo_cart(*cb);
    //print_stats("cart", cart, cerr);
    CartLowris clow(cb->do_medium ? 512: 256);
    clow.creo_cart_z_lowris(*cb, cart);
    //print_stats("clow", clow, cerr);

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

    Cart cart(cb->volume.max_beam_size());
    cart.creo_cart(*cb);
    //print_stats("cart", cart, cerr);
    wassert(actual(cart.cart).statsEqual(0, 28.22, 255));
    wassert(actual(cart.cartm).statsEqual(0, 0.00, 0));
    wassert(actual(cart.topxy).statsEqual(0, 1.75, 76));
    wassert(actual(cart.qual_Z_cart).statsEqual(0, 0.00, 0));
    wassert(actual(cart.quota_cart).statsEqual(0, 0.00, 0));
    wassert(actual(cart.dato_corr_xy).statsEqual(0, 0.00, 0));
    wassert(actual(cart.beam_blocking_xy).statsEqual(0, 0.00, 0));
    wassert(actual(cart.elev_fin_xy).statsEqual(0, 0.00, 0));
    wassert(actual(cart.neve_cart).statsEqual(0, 0.00, 0));
    wassert(actual(cart.corr_cart).statsEqual(0, 0.00, 0));
    wassert(actual(cart.conv_cart).statsEqual(0, 0.00, 0));

    CartLowris clow(cb->do_medium ? 512: 256);
    clow.creo_cart_z_lowris(*cb, cart);
    //print_stats("clow", clow, cerr);
    wassert(actual(clow.z_out).statsEqual(0, 24.45, 255));
    wassert(actual(clow.qual_Z_1x1).statsEqual(0, 0.00, 0));
    wassert(actual(clow.quota_1x1).statsEqual(0, 0.00, 0));
    wassert(actual(clow.dato_corr_1x1).statsEqual(0, 0.00, 0));
    wassert(actual(clow.elev_fin_1x1).statsEqual(0, 0.00, 0));
    wassert(actual(clow.beam_blocking_1x1).statsEqual(0, 0.00, 0));
    wassert(actual(clow.top_1x1).statsEqual(0, 0.00, 0));
    wassert(actual(clow.neve_1x1).statsEqual(0, 0.00, 0));
    wassert(actual(clow.corr_1x1).statsEqual(0, 0.00, 0));
    wassert(actual(clow.conv_1x1).statsEqual(0, 0.00, 0));

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
    wassert(actual(cb->volume).statsEqual(MISSING_DB, 0, -19.6863, -16.88, 65.01));

    wassert(actual(cb->beam_blocking).statsEqual(0, 21.36, 51));

    unsigned stats_size = cb->grid_stats.size_az * cb->grid_stats.size_beam;

    ArrayStats<unsigned> stat_anap_stats;
    stat_anap_stats.fill(cb->grid_stats.stat_anap, stats_size);
//    stat_anap_stats.print();
    wassert(actual((unsigned)stat_anap_stats.all_missing).isfalse());
    wassert(actual((unsigned)stat_anap_stats.min) == 0);
    wassert(actual((unsigned)stat_anap_stats.max) == 476);

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
    wassert(actual((unsigned)stat_elev_stats.max) == 476);

    wassert(actual(cb->dato_corrotto).statsEqual(0, 0.04, 1));

    cb->caratterizzo_volume();
    wassert(actual(cb->calcolo_vpr) != (void*)0);

    
    //print_stats("cb", *cb, cerr);
    wassert(actual(*cb->qual).statsEqual(0, 50.16, 99));
    wassert(actual(*cb->calcolo_vpr->flag_vpr).statsEqual(0, 0.87, 1));
    wassert(actual(cb->top).statsEqual(0, 0.54, 108));

    Cart cart(cb->volume.max_beam_size());
    cart.creo_cart(*cb);
    //print_stats("cart", cart, cerr);
    wassert(actual(cart.cart).statsEqual(0, 3.14, 255));
    wassert(actual(cart.cartm).statsEqual(0, 0.00, 0));
    wassert(actual(cart.topxy).statsEqual(0, 0.40, 108));
    wassert(actual(cart.qual_Z_cart).statsEqual(0, 9.35, 98));
    wassert(actual(cart.quota_cart).statsEqual(0, 2076.55, 7749));
    wassert(actual(cart.dato_corr_xy).statsEqual(0, 0.02, 1));
    wassert(actual(cart.beam_blocking_xy).statsEqual(0, 21.42, 51));
    wassert(actual(cart.elev_fin_xy).statsEqual(0, 0.10, 4));
    wassert(actual(cart.neve_cart).statsEqual(0, 0.00, 0));
    wassert(actual(cart.corr_cart).statsEqual(0, 185.63, 237));
    wassert(actual(cart.conv_cart).statsEqual(0, 0.00, 0));

    CartLowris clow(cb->do_medium ? 512: 256);
    clow.creo_cart_z_lowris(*cb, cart);
    //print_stats("clow", clow, cerr);
    wassert(actual(clow.z_out).statsEqual(0, 6.00, 243));
    wassert(actual(clow.qual_Z_1x1).statsEqual(0, 10.29, 97));
    wassert(actual(clow.quota_1x1).statsEqual(0, 147.17, 185));
    wassert(actual(clow.dato_corr_1x1).statsEqual(0, 0.01, 1));
    wassert(actual(clow.elev_fin_1x1).statsEqual(0, 0.13, 4));
    wassert(actual(clow.beam_blocking_1x1).statsEqual(0, 33.53, 51));
    wassert(actual(clow.top_1x1).statsEqual(0, 0.23, 49));
    wassert(actual(clow.neve_1x1).statsEqual(0, 0.00, 0));
    wassert(actual(clow.corr_1x1).statsEqual(0, 201.60, 237));
    wassert(actual(clow.conv_1x1).statsEqual(0, 0.00, 0));

    // TODO: scrivo_out_file_bin

    delete cb;
}


}
