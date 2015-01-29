#include <wibble/tests.h>
#include "cum_bac.h"
#include "config.h"
#include <elaboradar/logging.h>
#include <elaboradar/image.h>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include "site.h"
#include "setwork.h"
#include "cartproducts.h"
#include "test-utils.h"
#include <unistd.h>

using namespace wibble::tests;
using namespace elaboradar;
using namespace testradar;
using namespace std;

namespace tut {

struct process_shar {
    Logging logging;
};
TESTGRP(process);

namespace {

struct CBTest
{
    Config cfg;
    const Site& site;
    Volume<double> volume;
    bool do_medium;
    unsigned max_bin;

    CBTest(const char* site_name, bool do_medium, unsigned max_bin=512):
        site(Site::get(site_name)), do_medium(do_medium), max_bin(max_bin)
    {
    }

    void read_sp20(const char* fname, bool do_clean)
    {
        CUM_BAC::read_sp20_volume(volume, site, fname, do_clean, do_medium);
    }

    void read_odim(const char* fname, bool do_clean)
    {
        CUM_BAC::read_odim_volume(volume, site, fname, do_clean, do_medium);
    }

    unique_ptr<CUM_BAC> make_cumbac()
    {
        return unique_ptr<CUM_BAC>(new CUM_BAC(volume, cfg, site, do_medium, max_bin));
    }
};

}

template<> template<>
void to::test<1>()
{
    // Test elabora_dato, con tutti i do_* a true
    static const char* fname = "../testdata/DBP2_070120141530_GATTATICO";

    unsetwork();
    setenv("FIRST_LEVEL_FILE", "../dati/FIRST_LEVEL_corto_GAT_2006_INV", 1);
    setenv("DIR_OUT_PP_BLOC", "../testdata", 1);
    setenv("FILE_T", "../testdata/temperature.txt", 1);
    setenv("LAST_VPR","../testdata/last_vpr",1);
    setenv("FILE_ZERO_TERMICO","../testdata/zero_termico.txt",1);

    CBTest test("GAT", false);
    test.read_sp20(fname, true);
    auto cb = test.make_cumbac();
    cb->want_vpr();
    cb->do_quality = true;
    cb->do_beamblocking = true;
    cb->do_declutter = false;
    cb->do_bloccorr = true;
    cb->do_class = true;
    cb->do_readStaticMap=true;

    wassert(actual((int)(cb->calcolo_vpr->t_ground * 100)) == 1010);

    cb->declutter_anaprop();

    // Check results
    VolumeStats stats;
    cb->volume.compute_stats(stats);
    // stats.print(stdout);
    wassert(actual(stats.count_zeros[0]) == 139396);
    wassert(actual(stats.count_zeros[1]) == 181135);
    wassert(actual(stats.count_zeros[2]) == 192056);
    wassert(actual(stats.count_zeros[3]) == 196000);
    wassert(actual(stats.count_zeros[4]) == 195715);
    wassert(actual(stats.count_zeros[5]) == 195897);
    wassert(actual(stats.count_ones[0]) ==  97);
    wassert(actual(stats.count_ones[1]) == 217);
    wassert(actual(stats.count_ones[2]) == 119);
    wassert(actual(stats.count_ones[3]) ==  70);
    wassert(actual(stats.count_ones[4]) ==  66);
    wassert(actual(stats.count_ones[5]) ==  43);
    wassert(actual(stats.count_others[0]) == 58107);
    wassert(actual(stats.count_others[1]) == 16248);
    wassert(actual(stats.count_others[2]) ==  5425);
    wassert(actual(stats.count_others[3]) ==  1530);
    wassert(actual(stats.count_others[4]) ==  1819);
    wassert(actual(stats.count_others[5]) ==  1660);
    wassert(actual(stats.sum_others[0]) == 5475665);
    wassert(actual(stats.sum_others[1]) == 1130268);
    wassert(actual(stats.sum_others[2]) ==  332252);
    wassert(actual(stats.sum_others[3]) ==   52204);
    wassert(actual(stats.sum_others[4]) ==  100286);
    wassert(actual(stats.sum_others[5]) ==  101241);
    cb->caratterizzo_volume();

//    print_stats("cb", *cb, cerr);
    wassert(actual(*cb->qual).statsEqual(1, 61.77, 99));
    wassert(actual(cb->top).statsEqual(0, 187460, 1, 7.34, 36));
    wassert(actual(cb->first_level).statsEqual(0, 0.86, 3));
    wassert(actual(cb->first_level_static).statsEqual(0, 0.54, 3));
    wassert(actual(cb->bb_first_level).statsEqual(0, 0.49, 2));
    wassert(actual(cb->beam_blocking).statsEqual(0, 4.72, 49));
    wassert(actual(cb->dem).statsEqual(-9999, 47.59, 4070.26));
    wassert(actual(cb->anaprop.quota).statsEqual(0, 0, 0));
    wassert(actual(cb->anaprop.dato_corrotto).statsEqual(0, 0, 0));
    wassert(actual(*cb->calcolo_vpr->flag_vpr).statsEqual(0, 0.99, 1));
    wassert(actual(cb->calcolo_vpr->corr_polar).statsEqual(0, 0, 0));
    wassert(actual(cb->calcolo_vpr->neve).statsEqual(0, 0, 0));
}

template<> template<>
void to::test<2>()
{
    // Test elabora_dato, con tutti i do_* a true
    static const char* fname = "../testdata/DBP2_070120141530_GATTATICO";

    unsetwork();
    setenv("FIRST_LEVEL_FILE", "../dati/FIRST_LEVEL_corto_GAT_2006_INV", 1);
    setenv("DIR_OUT_PP_BLOC", "../testdata", 1);

    CBTest test("GAT", false);
    test.read_sp20(fname, true);
    auto cb = test.make_cumbac();
    cb->want_vpr();
    cb->do_quality = true;
    cb->do_beamblocking = true;
    cb->do_declutter = false;
    cb->do_bloccorr = true;
    cb->do_readStaticMap=true;

    cb->declutter_anaprop();

    // Check results
    VolumeStats stats;
    cb->volume.compute_stats(stats);
    //stats.print(stdout);
    wassert(actual(stats.count_zeros[0]) == 139396);
    wassert(actual(stats.count_zeros[1]) == 181135);
    wassert(actual(stats.count_zeros[2]) == 192056);
    wassert(actual(stats.count_zeros[3]) == 196000);
    wassert(actual(stats.count_zeros[4]) == 195715);
    wassert(actual(stats.count_zeros[5]) == 195897);
    wassert(actual(stats.count_ones[0]) ==  97);
    wassert(actual(stats.count_ones[1]) == 217);
    wassert(actual(stats.count_ones[2]) == 119);
    wassert(actual(stats.count_ones[3]) ==  70);
    wassert(actual(stats.count_ones[4]) ==  66);
    wassert(actual(stats.count_ones[5]) ==  43);
    wassert(actual(stats.count_others[0]) == 58107);
    wassert(actual(stats.count_others[1]) == 16248);
    wassert(actual(stats.count_others[2]) ==  5425);
    wassert(actual(stats.count_others[3]) ==  1530);
    wassert(actual(stats.count_others[4]) ==  1819);
    wassert(actual(stats.count_others[5]) ==  1660);
    wassert(actual(stats.sum_others[0]) == 5475665);
    wassert(actual(stats.sum_others[1]) == 1130268);
    wassert(actual(stats.sum_others[2]) ==  332252);
    wassert(actual(stats.sum_others[3]) ==   52204);
    wassert(actual(stats.sum_others[4]) ==  100286);
    wassert(actual(stats.sum_others[5]) ==  101241);
}

template<> template<>
void to::test<3>()
{
    // versione BB che corrisponde al parametro algo_corto
    static const char* fname = "../testdata/DBP2_070120141530_GATTATICO";
    unsetwork();
    setenv("FIRST_LEVEL_FILE", "../dati/FIRST_LEVEL_corto_GAT_2006_INV", 1);
    setenv("DIR_OUT_PP_BLOC", "../testdata", 1);
    setenv("VPR_HEATING", "../testdata/vpr_heat_GAT", 1);
    unlink("../testdata/vpr_heat_GAT");
    setenv("FILE_T", "../testdata/temperature.txt", 1);

    CBTest test("GAT", false);
    test.read_sp20(fname, true);
    auto cb = test.make_cumbac();
    cb->want_vpr();
    cb->do_quality = true;
    cb->do_beamblocking = true;
    cb->do_declutter =false ;
    cb->do_class = false;
    cb->do_bloccorr = false;
    cb->do_readStaticMap=true;

    cb->declutter_anaprop();

    cb->caratterizzo_volume();
//    print_stats("cb", *cb, cerr);
    wassert(actual(*cb->qual).statsEqual(1, 60.01, 99));
    wassert(actual(cb->top).statsEqual(0, 187460, 1, 7.34, 36));
    wassert(actual(cb->first_level).statsEqual(0, 0.54, 3));
    wassert(actual(cb->first_level_static).statsEqual(0, 0.54, 3));
    wassert(actual(cb->bb_first_level).statsEqual(0, 0.49, 2));
    wassert(actual(cb->beam_blocking).statsEqual(0, 10.44, 51));
    wassert(actual(cb->dem).statsEqual(-9999, 47.59, 4070.26));
    wassert(actual(cb->anaprop.quota).statsEqual(0, 0, 0));
    wassert(actual(cb->anaprop.dato_corrotto).statsEqual(0, 0, 0));
    wassert(actual(*cb->calcolo_vpr->flag_vpr).statsEqual(0, 0.96, 1));
    wassert(actual(cb->calcolo_vpr->corr_polar).statsEqual(0, 0, 0));
    wassert(actual(cb->calcolo_vpr->neve).statsEqual(0, 0, 0));

    Cart cart(cb->volume.max_beam_size());
    wassert(actual(cart.max_bin) == 494);
    cart.creo_cart(*cb);
    //setenv("DIR_DEBUG", "../testdata/", 1);
    //cart.write_out(*cb, cb->assets);
//     print_stats("cart", cart, cout);
    wassert(actual(cart.cart).statsEqual(0, 213127, 1, 27.13, 251));
    wassert(actual(cart.cartm).statsEqual(0, 976144, 0, 0, 0));
    wassert(actual(cart.topxy).statsEqual(0, 935454, 1, 11, 36));
    wassert(actual(cart.qual_Z_cart).statsEqual(0, 213127, 1, 36.89, 98));
    wassert(actual(cart.quota_cart).statsEqual(0, 976144, 0, 0, 0));
    wassert(actual(cart.dato_corr_xy).statsEqual(0, 976144, 0, 0, 0));
    wassert(actual(cart.beam_blocking_xy).statsEqual(0, 693564, 1, 47.08, 51));
    wassert(actual(cart.elev_fin_xy).statsEqual(0, 976144, 0, 0, 0));
    wassert(actual(cart.neve_cart).statsEqual(0, 213127, 1, 1, 1));
    wassert(actual(cart.corr_cart).statsEqual(0, 976144, 0, 0, 0));
    wassert(actual(cart.conv_cart).statsEqual(0, 976144, 0, 0, 0));

    CartLowris clow(cb->do_medium ? 512: 256, *cb, cart);
    wassert(actual(clow.CART_DIM_ZLR) == 256);
    wassert(actual(clow.ZLR_N_ELEMENTARY_PIXEL) == 4);
    wassert(actual(clow.ZLR_OFFSET) == -18);
    clow.creo_cart_z_lowris();
//     print_stats("clow", clow, cerr);
    wassert(actual(clow.z_out).statsEqual(0, 17489, 1, 40.13, 251));
    wassert(actual(clow.qual_Z_1x1).statsEqual(0, 17489, 1, 37.21, 97));
    wassert(actual(clow.quota_1x1).statsEqual(0, 0, 128, 128, 128));
    wassert(actual(clow.dato_corr_1x1).statsEqual(0, 65536, 0, 0, 0));
    wassert(actual(clow.elev_fin_1x1).statsEqual(0, 65536, 0, 0, 0));
    wassert(actual(clow.beam_blocking_1x1).statsEqual(0, 48079, 1, 47.09, 51));
    wassert(actual(clow.top_1x1).statsEqual(0, 60383, 1, 10.71, 36));
    wassert(actual(clow.neve_1x1).statsEqual(0, 17489, 1, 1, 1));
    wassert(actual(clow.corr_1x1).statsEqual(0, 65536, 0, 0, 0));
    wassert(actual(clow.conv_1x1).statsEqual(0, 65536, 0, 0, 0));

    CartProducts products(cb->volume, 256, 4);
    cb->generate_maps(products, true);
    wassert(actual(products.scaled.image_offset) == -18);
    // write_image(cb->volume[0], "/tmp/zout-raw.png", "png");
    // write_image(clow.z_out, "/tmp/zout-old.png", "png");
    // write_image(products.z_out, "/tmp/zout-new.png", "png");
    // print_stats("products", products, cerr);
    wassert(actual(products.z_out).statsEqual(0, 17240, 1, 40.23, 251));
    wassert(actual(products.qual_Z_1x1).statsEqual(0, 17240, 1, 37.14, 98));
    wassert(actual(products.quota_1x1).statsEqual(0, 0, 128, 128, 128));
    wassert(actual(products.dato_corr_1x1).statsEqual(0, 65536, 0, 0, 0));
    wassert(actual(products.elev_fin_1x1).statsEqual(0, 65536, 0, 0, 0));
    wassert(actual(products.beam_blocking_1x1).statsEqual(0, 47936, 1, 47.24, 51));
    wassert(actual(products.top_1x1).statsEqual(0, 60304, 1, 10.73, 36));
    wassert(actual(products.neve_1x1).statsEqual(0, 17240, 1, 1, 1));
    wassert(actual(products.corr_1x1).statsEqual(0, 65536, 0, 0, 0));
    wassert(actual(products.conv_1x1).statsEqual(0, 65536, 0, 0, 0));


    // Try again, but with averages
    cb->do_zlr_media = true;

    Cart cart1(cb->volume.max_beam_size());
    cart1.creo_cart(*cb);
    // print_stats("cart1", cart1, cerr);
    wassert(actual(cart1.cart).statsEqual(0, 213127, 1, 27.13, 251));
    wassert(actual(cart1.cartm).statsEqual(0, 213127, 0.01, 327.25, 749048.25));
    wassert(actual(cart1.topxy).statsEqual(0, 935454, 1, 11, 36));
    wassert(actual(cart1.qual_Z_cart).statsEqual(0, 213127, 1, 36.89, 98));
    wassert(actual(cart1.quota_cart).statsEqual(0, 976144, 0, 0, 0));
    wassert(actual(cart1.dato_corr_xy).statsEqual(0, 976144, 0, 0, 0));
    wassert(actual(cart1.beam_blocking_xy).statsEqual(0, 693564, 1, 47.08, 51));
    wassert(actual(cart1.elev_fin_xy).statsEqual(0, 976144, 0, 0, 0));
    wassert(actual(cart1.neve_cart).statsEqual(0, 213127, 1, 1, 1));
    wassert(actual(cart1.corr_cart).statsEqual(0, 976144, 0, 0, 0));
    wassert(actual(cart1.conv_cart).statsEqual(0, 976144, 0, 0, 0));

    CartLowris clow1(cb->do_medium ? 512: 256, *cb, cart1);
    wassert(actual(clow1.CART_DIM_ZLR) == 256);
    wassert(actual(clow1.ZLR_N_ELEMENTARY_PIXEL) == 4);
    wassert(actual(clow1.ZLR_OFFSET) == -18);
    clow1.creo_cart_z_lowris();
    //print_stats("clow1", clow1, cerr);
    wassert(actual(clow1.z_out).statsEqual(0, 17489, 1, 33.31, 235));
    wassert(actual(clow1.qual_Z_1x1).statsEqual(0, 17489, 1, 37.21, 97));
    wassert(actual(clow1.quota_1x1).statsEqual(0, 0, 128, 128, 128));
    wassert(actual(clow1.dato_corr_1x1).statsEqual(0, 65536, 0, 0, 0));
    wassert(actual(clow1.elev_fin_1x1).statsEqual(0, 65536, 0, 0, 0));
    wassert(actual(clow1.beam_blocking_1x1).statsEqual(0, 48079, 1, 47.09, 51));
    wassert(actual(clow1.top_1x1).statsEqual(0, 60383, 1, 10.71, 36));
    wassert(actual(clow1.neve_1x1).statsEqual(0, 17489, 1, 1, 1));
    wassert(actual(clow1.corr_1x1).statsEqual(0, 65536, 0, 0, 0));
    wassert(actual(clow1.conv_1x1).statsEqual(0, 65536, 0, 0, 0));

    CartProducts products1(cb->volume, 256, 4);
    cb->generate_maps(products1, true);
    wassert(actual(products1.scaled.image_offset) == -18);
    // write_image(cb->volume[0], "/tmp/zout-raw.png", "png");
    // write_image(clow.z_out, "/tmp/zout-old.png", "png");
    // write_image(products.z_out, "/tmp/zout-new.png", "png");
    // print_stats("products", products, cerr);
    wassert(actual(products1.z_out).statsEqual(0, 17489, 1, 40.13, 251));
    wassert(actual(products1.qual_Z_1x1).statsEqual(0, 17489, 1, 37.21, 97));
    wassert(actual(products1.quota_1x1).statsEqual(0, 0, 128, 128, 128));
    wassert(actual(products1.dato_corr_1x1).statsEqual(0, 65536, 0, 0, 0));
    wassert(actual(products1.elev_fin_1x1).statsEqual(0, 65536, 0, 0, 0));
    wassert(actual(products1.beam_blocking_1x1).statsEqual(0, 48079, 1, 47.09, 51));
    wassert(actual(products1.top_1x1).statsEqual(0, 60383, 1, 10.71, 36));
    wassert(actual(products1.neve_1x1).statsEqual(0, 17489, 1, 1, 1));
    wassert(actual(products1.corr_1x1).statsEqual(0, 65536, 0, 0, 0));
    wassert(actual(products1.conv_1x1).statsEqual(0, 65536, 0, 0, 0));
}

template<> template<>
void to::test<4>()
{
    // versione BB_VPR che corrisponde al parametro algo_corto_dev
	LOG_CATEGORY("Test");

    static const char* fname = "../testdata/DBP2_070120141530_GATTATICO";
    unsetwork();
    setenv("FIRST_LEVEL_FILE", "../dati/FIRST_LEVEL_corto_GAT_2006_INV", 1);
    setenv("DIR_OUT_PP_BLOC", "../testdata", 1);
    setenv("VPR_HEATING", "../testdata/vpr_heat_GAT", 1);
    setenv("VPR0_FILE", "../testdata/ultimo_vpr", 1);
    unlink("../testdata/vpr_heat_GAT");
    setenv("FILE_T", "../testdata/temperature.txt", 1);
    setenv("LAST_VPR","../testdata/last_vpr",1);
    printwork();

    CBTest test("GAT", false);
    test.read_sp20(fname, true);
    auto cb = test.make_cumbac();
    cb->want_vpr();
    cb->do_quality = true;
    cb->do_beamblocking = true;
    cb->do_declutter = true;
    cb->do_class = false;
    cb->do_bloccorr = false;
    cb->do_readStaticMap=true;
LOG_INFO ("Chiamo elabora_dato");
    cb->declutter_anaprop();
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
    wassert(actual(cb->anaprop.quota).statsEqual(0, 0, 0));
    wassert(actual(cb->anaprop.dato_corrotto).statsEqual(0, 0, 0));
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
//    print_stats("clow", clow, cerr);
    wassert(actual(clow.z_out).statsEqual(0, 17489, 1, 25.84, 242));
    wassert(actual(clow.qual_Z_1x1).statsEqual(0, 17489, 1, 35.66, 98));
    wassert(actual(clow.quota_1x1).statsEqual(0, 0, 128, 128, 128));
    wassert(actual(clow.dato_corr_1x1).statsEqual(0, 65536, 0, 0, 0));
    wassert(actual(clow.elev_fin_1x1).statsEqual(0, 53987, 1, 1.13, 3));
    wassert(actual(clow.beam_blocking_1x1).statsEqual(0, 46817, 1, 47.04, 51));
    wassert(actual(clow.top_1x1).statsEqual(0, 62100, 1, 11.39, 36));
    wassert(actual(clow.neve_1x1).statsEqual(0, 17489, 1, 1, 1));
    wassert(actual(clow.corr_1x1).statsEqual(0, 65536, 0, 0, 0));
    wassert(actual(clow.conv_1x1).statsEqual(0, 65536, 0, 0, 0));
}

template<> template<>
void to::test<5>()
{
    // versione BB_VPR_CLASS che corrisponde al parametro algo_corto_dev
	LOG_CATEGORY("Test");
	LOG_INFO ("Start test 5");

    static const char* fname = "../testdata/DBP2_070120141530_GATTATICO";
    unsetwork();
    setenv("FIRST_LEVEL_FILE", "../dati/FIRST_LEVEL_corto_GAT_2006_INV", 1);
    setenv("DIR_OUT_PP_BLOC", "../testdata", 1);
    setenv("VPR_HEATING", "../testdata/vpr_heat_GAT", 1);
    unlink("../testdata/vpr_heat_GAT");
    setenv("FILE_T", "../testdata/temperature.txt", 1);
    setenv("LAST_VPR","../testdata/last_vpr",1);
unlink("LAST_VPR");
    setenv("FILE_ZERO_TERMICO","../testdata/zero_termico.txt",1);
    printwork();

    CBTest test("GAT", false);
    test.read_sp20(fname, true);
    auto cb = test.make_cumbac();
    cb->want_vpr();
    cb->do_quality = true;
    cb->do_beamblocking = true;
    cb->do_declutter = true;
    cb->do_class = true;
    cb->do_bloccorr = false;
    cb->do_readStaticMap=true;

    cb->declutter_anaprop();

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
    wassert(actual(cb->top).statsEqual(0, 187460, 1, 7.34, 36));
    wassert(actual(cb->first_level).statsEqual(0, 0.54, 3));
    wassert(actual(cb->first_level_static).statsEqual(0, 0.54, 3));
    wassert(actual(cb->bb_first_level).statsEqual(0, 0.49, 2));
    wassert(actual(cb->beam_blocking).statsEqual(0, 10.44, 51));
    wassert(actual(cb->dem).statsEqual(-9999, 47.59, 4070.26));
    wassert(actual(cb->anaprop.quota).statsEqual(0, 0, 0));
    wassert(actual(cb->anaprop.dato_corrotto).statsEqual(0, 0, 0));
    wassert(actual(*cb->calcolo_vpr->flag_vpr).statsEqual(0, 0.7, 1));
    wassert(actual(cb->calcolo_vpr->corr_polar).statsEqual(0, 0, 0));
    wassert(actual(cb->calcolo_vpr->neve).statsEqual(0, 0, 0));

    // TODO: cb->stampa_vpr()

    cb->conversione_convettiva();

    Cart cart(cb->volume.max_beam_size());
    wassert(actual(cart.max_bin) == 494);
    cart.creo_cart(*cb);
//    print_stats("cart", cart, cerr);
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
//     print_stats("clow", clow, cerr);
    wassert(actual(clow.z_out).statsEqual(0, 17489, 1, 25.84, 242));
    wassert(actual(clow.qual_Z_1x1).statsEqual(0, 17489, 1, 35.66, 98));
    wassert(actual(clow.quota_1x1).statsEqual(0, 0, 128, 128, 128));
    wassert(actual(clow.dato_corr_1x1).statsEqual(0, 65536, 0, 0, 0));
    wassert(actual(clow.elev_fin_1x1).statsEqual(0, 53987, 1, 1.13, 3));
    wassert(actual(clow.beam_blocking_1x1).statsEqual(0, 46817, 1, 47.04, 51));
    wassert(actual(clow.top_1x1).statsEqual(0, 62100, 1, 11.39, 36));
    wassert(actual(clow.neve_1x1).statsEqual(0, 17489, 1, 1, 1));
    wassert(actual(clow.corr_1x1).statsEqual(0, 65536, 0, 0, 0));
    wassert(actual(clow.conv_1x1).statsEqual(0, 65536, 0, 0, 0));
    // TODO: scrivo_out_file_bin
}

template<> template<>
void to::test<6>()
{
	LOG_CATEGORY("Test");
	LOG_INFO ("Start test 6");

    // versione BB che corrisponde al parametro algo_corto
    static const char* fname = "../testdata/DBP2_060220140140_GATTATICO";
    unsetwork();
    setenv("FIRST_LEVEL_FILE", "../dati/FIRST_LEVEL_corto_GAT_2006_INV", 1);
    setenv("DIR_OUT_PP_BLOC", "../testdata", 1);
    setenv("VPR_HEATING", "../testdata/vpr_heat_GAT", 1);
    setenv("FILE_ZERO_TERMICO", "../testdata/20140206/0termico.prev", 1);
    unlink("../testdata/vpr_heat_GAT");
    setenv("VPR0_FILE", "../testdata/ultimo_vpr", 1);
    unlink("../testdata/ultimo_vpr");
    setenv("FILE_T", "../testdata/temperature.txt", 1);
	printwork();

    CBTest test("GAT", false);
    test.read_sp20(fname, true);
    auto cb = test.make_cumbac();
    cb->want_vpr();
    cb->do_quality = true;
    cb->do_beamblocking = true;
    cb->do_declutter = false;
    cb->do_class = false;
    cb->do_bloccorr = false;
    cb->do_readStaticMap=true;
    wassert(actual(cb->calcolo_vpr) != (void*)0);

    cb->declutter_anaprop();

    VolumeStats stats;
    cb->volume.compute_stats(stats);
    //stats.print(stdout);
    wassert(actual(stats.count_zeros[0]) == 165040);
    wassert(actual(stats.count_zeros[1]) == 207601);
    wassert(actual(stats.count_zeros[2]) == 161283);
    wassert(actual(stats.count_zeros[3]) == 108238);
    wassert(actual(stats.count_zeros[4]) == 117678);
    wassert(actual(stats.count_ones[0]) == 0);
    wassert(actual(stats.count_ones[1]) == 0);
    wassert(actual(stats.count_ones[2]) == 0);
    wassert(actual(stats.count_ones[3]) == 0);
    wassert(actual(stats.count_ones[4]) == 0);
    wassert(actual(stats.count_others[0]) == 174560);
    wassert(actual(stats.count_others[1]) == 131999);
    wassert(actual(stats.count_others[2]) == 102717);
    wassert(actual(stats.count_others[3]) ==  89362);
    wassert(actual(stats.count_others[4]) ==  79922);
    wassert(actual(stats.sum_others[0]) == 24480412);
    wassert(actual(stats.sum_others[1]) == 17301958);
    wassert(actual(stats.sum_others[2]) == 12929034);
    wassert(actual(stats.sum_others[3]) == 10855770);
    wassert(actual(stats.sum_others[4]) ==  9385431);

    wassert(actual(cb->beam_blocking).statsEqual(0, 17.1, 51));

    unsigned stats_size = cb->anaprop.grid_stats.size_az * cb->anaprop.grid_stats.size_beam;

    ArrayStats<unsigned> stat_anap_stats;
    stat_anap_stats.fill(cb->anaprop.grid_stats.stat_anap, stats_size);
    //stat_anap_stats.print();
    wassert(actual((unsigned)stat_anap_stats.all_missing).isfalse());
    wassert(actual((unsigned)stat_anap_stats.min) == 0);
    wassert(actual((unsigned)stat_anap_stats.max) == 0);

    ArrayStats<unsigned> stat_anap_tot_stats;
    stat_anap_tot_stats.fill(cb->anaprop.grid_stats.stat_tot, stats_size);
    //stat_anap_tot_stats.print();
    wassert(actual((unsigned)stat_anap_tot_stats.all_missing).isfalse());
    wassert(actual((unsigned)stat_anap_tot_stats.min) ==  0);
    wassert(actual((unsigned)stat_anap_tot_stats.max) == 0);

    ArrayStats<unsigned> stat_bloc_stats;
    stat_bloc_stats.fill(cb->anaprop.grid_stats.stat_bloc, stats_size);
    //stat_bloc_stats.print();
    wassert(actual((unsigned)stat_bloc_stats.all_missing).isfalse());
    wassert(actual((unsigned)stat_bloc_stats.min) == 0);
    wassert(actual((unsigned)stat_bloc_stats.max) == 0);

    ArrayStats<unsigned> stat_elev_stats;
    stat_elev_stats.fill(cb->anaprop.grid_stats.stat_elev, stats_size);
    //stat_elev_stats.print();
    wassert(actual((unsigned)stat_elev_stats.all_missing).isfalse());
    wassert(actual((unsigned)stat_elev_stats.min) == 0);
    wassert(actual((unsigned)stat_elev_stats.max) == 0);

    ArrayStats<unsigned char> dato_corrotto_stats;
    dato_corrotto_stats.fill(cb->anaprop.dato_corrotto);
    //dato_corrotto_stats.print();
    wassert(actual((unsigned)dato_corrotto_stats.all_missing).isfalse());
    wassert(actual((unsigned)dato_corrotto_stats.min) == 0);
    wassert(actual((unsigned)dato_corrotto_stats.max) == 0);


    cb->caratterizzo_volume();
    wassert(actual(cb->calcolo_vpr) != (void*)0);
//    print_stats("cb", *cb, cerr);
    wassert(actual(*cb->qual).statsEqual(1, 55.99, 99));
    wassert(actual(cb->top).statsEqual(0, 222889, 1, 13.23, 76));
    wassert(actual(cb->first_level).statsEqual(0, 0.33, 3));
    wassert(actual(cb->first_level_static).statsEqual(0, 0.33, 3));
    wassert(actual(cb->bb_first_level).statsEqual(0, 0.48, 2));
    wassert(actual(cb->beam_blocking).statsEqual(0, 17.1, 51));
    wassert(actual(cb->dem).statsEqual(-9999, 47.59, 4070.26));
    wassert(actual(cb->anaprop.quota).statsEqual(0, 0, 0));
    wassert(actual(cb->anaprop.dato_corrotto).statsEqual(0, 0, 0));
    wassert(actual(*cb->calcolo_vpr->flag_vpr).statsEqual(0, 0.94, 1));
    wassert(actual(cb->calcolo_vpr->corr_polar).statsEqual(0, 0, 0));
    wassert(actual(cb->calcolo_vpr->neve).statsEqual(0, 0, 0));

    Cart cart(cb->volume.max_beam_size());
    wassert(actual(cart.max_bin) == 849);
    cart.creo_cart(*cb);
//     print_stats("cart", cart, cerr);
    wassert(actual(cart.cart).statsEqual(0, 624947, 1, 46.62, 255));
    wassert(actual(cart.cartm).statsEqual(0, 2883204, 0, 0, 0));
    wassert(actual(cart.topxy).statsEqual(0, 2511651, 1, 16.37, 76));
    wassert(actual(cart.qual_Z_cart).statsEqual(0, 624947, 1, 18.42, 98));
    wassert(actual(cart.quota_cart).statsEqual(0, 2883204, 0, 0, 0));
    wassert(actual(cart.dato_corr_xy).statsEqual(0, 2883204, 0, 0, 0));
    wassert(actual(cart.beam_blocking_xy).statsEqual(0, 1692726, 1, 45.11, 51));
    wassert(actual(cart.elev_fin_xy).statsEqual(0, 2883204, 0, 0, 0));
    wassert(actual(cart.neve_cart).statsEqual(0, 624947, 1, 1, 1));
    wassert(actual(cart.corr_cart).statsEqual(0, 2883204, 0, 0, 0));
    wassert(actual(cart.conv_cart).statsEqual(0, 2883204, 0, 0, 0));

    CartLowris clow(cb->do_medium ? 512: 256, *cb, cart);
    wassert(actual(clow.ZLR_N_ELEMENTARY_PIXEL) == 4);
    wassert(actual(clow.ZLR_OFFSET) == 337);
    clow.creo_cart_z_lowris();
//     print_stats("clow", clow, cerr);
    wassert(actual(clow.z_out).statsEqual(0, 0, 1, 89.28, 255));
    wassert(actual(clow.qual_Z_1x1).statsEqual(0, 0, 1, 27.69, 97));
    wassert(actual(clow.quota_1x1).statsEqual(0, 0, 128, 128, 128));
    wassert(actual(clow.dato_corr_1x1).statsEqual(0, 65536, 0, 0, 0));
    wassert(actual(clow.elev_fin_1x1).statsEqual(0, 65536, 0, 0, 0));
    wassert(actual(clow.beam_blocking_1x1).statsEqual(0, 39117, 1, 45.59, 51));
    wassert(actual(clow.top_1x1).statsEqual(0, 40993, 3, 13.89, 57));
    wassert(actual(clow.neve_1x1).statsEqual(0, 0, 1, 1, 1));
    wassert(actual(clow.corr_1x1).statsEqual(0, 65536, 0, 0, 0));
    wassert(actual(clow.conv_1x1).statsEqual(0, 65536, 0, 0, 0));
}

template<> template<>
void to::test<7>()
{
	LOG_CATEGORY("Test");
	LOG_INFO ("Start test 7");

    // versione BB_VPR che corrisponde al parametro algo_corto_dev
    static const char* fname = "../testdata/DBP2_060220140140_GATTATICO";
    unsetwork();
    setenv("FIRST_LEVEL_FILE", "../dati/FIRST_LEVEL_corto_GAT_2006_INV", 1);
    setenv("DIR_OUT_PP_BLOC", "../testdata", 1);
    setenv("FILE_T", "../testdata/temperature.txt", 1);
    setenv("FILE_ZERO_TERMICO", "../testdata/20140206/0termico.prev", 1);
    setenv("VPR_HEATING", "../testdata/vpr_heat_GAT", 1);
    unlink("../testdata/vpr_heat_GAT");
    setenv("VPR0_FILE", "../testdata/ultimo_vpr", 1);
    unlink("../testdata/ultimo_vpr");
    setenv("LAST_VPR", "../testdata/last_vpr", 1);
    unlink("../testdata/last_vpr");
    setenv("VPR_HMAX", "../testdata/vpr_hmax", 1);
    setenv("TEST_VPR", "../testdata/test_vpr", 1);
    setenv("VPR_ARCH", "../testdata/vpr_arch",1);
	printwork();

    CBTest test("GAT", false);
    test.read_sp20(fname, true);
    auto cb = test.make_cumbac();
    cb->want_vpr();
    cb->do_quality = true;
    cb->do_beamblocking = true;
    cb->do_declutter = true;
    cb->do_class = false;
    cb->do_bloccorr = false;
    cb->do_readStaticMap=true;

    cb->declutter_anaprop();

    cb->caratterizzo_volume();

    // la combina_profili restituisce 1 se non riesce a costruire un profilo
    // perchè non piove o piove poco
   cb->calcolo_vpr->test_vpr=fopen("../testdata/test_vpr","a+");

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
    wassert(actual(cb->anaprop.quota).statsEqual(0, 0, 0));
    wassert(actual(cb->anaprop.dato_corrotto).statsEqual(0, 0, 0));
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
//     print_stats("clow", clow, cerr);
    wassert(actual(clow.z_out).statsEqual(0, 0, 1, 81.19, 255));
    wassert(actual(clow.qual_Z_1x1).statsEqual(0, 0, 1, 28.5, 97));
    wassert(actual(clow.quota_1x1).statsEqual(0, 0, 128, 128, 128));
    wassert(actual(clow.dato_corr_1x1).statsEqual(0, 65536, 0, 0, 0));
    wassert(actual(clow.elev_fin_1x1).statsEqual(0, 55911, 1, 1.16, 3));
    wassert(actual(clow.beam_blocking_1x1).statsEqual(0, 37368, 1, 45.35, 51));
    wassert(actual(clow.top_1x1).statsEqual(0, 42175, 2, 13.88, 57));
    wassert(actual(clow.neve_1x1).statsEqual(0, 0, 1, 1, 1));
    wassert(actual(clow.corr_1x1).statsEqual(0, 65536, 0, 0, 0));
    wassert(actual(clow.conv_1x1).statsEqual(0, 65536, 0, 0, 0));
}

template<> template<>
void to::test<8>()
{
    // versione BB_VPR_CLASS che corrisponde al parametro algo_corto_dev
    static const char* fname = "../testdata/DBP2_060220140140_GATTATICO";
    unsetwork();
    setenv("FIRST_LEVEL_FILE", "../dati/FIRST_LEVEL_corto_GAT_2006_INV", 1);
    setenv("DIR_OUT_PP_BLOC", "../testdata", 1);
    setenv("FILE_T", "../testdata/temperature.txt", 1);
    setenv("FILE_ZERO_TERMICO", "../testdata/20140206/0termico.prev", 1);
    setenv("VPR_HEATING", "../testdata/vpr_heat_GAT", 1);
    unlink("../testdata/vpr_heat_GAT");
    setenv("VPR0_FILE", "../testdata/ultimo_vpr", 1);
    unlink("../testdata/ultimo_vpr");
    setenv("LAST_VPR", "../testdata/last_vpr", 1);
    unlink("../testdata/last_vpr");
    setenv("VPR_HMAX", "../testdata/vpr_hmax", 1);
    setenv("VPR_ARCH", "../testdata/vpr_arch",1);

    CBTest test("GAT", false);
    test.read_sp20(fname, true);
    auto cb = test.make_cumbac();
    cb->want_vpr();
    cb->do_quality = true;
    cb->do_beamblocking = true;
    cb->do_declutter = true;
    cb->do_class = true;
    cb->do_bloccorr = false;
    cb->do_readStaticMap=true;

    cb->declutter_anaprop();
    cb->caratterizzo_volume();
    cb->calcolo_vpr->classifica_rain();

    // la combina_profili restituisce 1 se non riesce a costruire un profilo
    // perchè non piove o piove poco
    cb->calcolo_vpr->test_vpr=fopen("../testdata/test_vpr","a+");
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
    wassert(actual(cb->anaprop.quota).statsEqual(0, 0, 0));
    wassert(actual(cb->anaprop.dato_corrotto).statsEqual(0, 0, 0));
    wassert(actual(*cb->calcolo_vpr->flag_vpr).statsEqual(0, 0.66, 1));
    wassert(actual(cb->calcolo_vpr->corr_polar).statsEqual(0, 0, 0));
    wassert(actual(cb->calcolo_vpr->neve).statsEqual(0, 0, 0));

    Cart cart(cb->volume.max_beam_size());
    wassert(actual(cart.max_bin) == 849);
    cart.creo_cart(*cb);
//     print_stats("cart", cart, cerr);
    wassert(actual(cart.cart).statsEqual(0, 624947, 1, 41.64, 249));
    wassert(actual(cart.cartm).statsEqual(0, 2883204, 0, 0, 0));
    wassert(actual(cart.topxy).statsEqual(0, 2516548, 2, 16.42, 76));
    wassert(actual(cart.qual_Z_cart).statsEqual(0, 624947, 1, 19.13, 98));
    wassert(actual(cart.quota_cart).statsEqual(0, 2883204, 0, 0, 0));
    wassert(actual(cart.dato_corr_xy).statsEqual(0, 2883204, 0, 0, 0));
    wassert(actual(cart.beam_blocking_xy).statsEqual(0, 1687643, 1, 45.08, 51));
    wassert(actual(cart.elev_fin_xy).statsEqual(0, 2683482, 1, 1.15, 3));
    wassert(actual(cart.neve_cart).statsEqual(0, 624947, 1, 1, 1));
    wassert(actual(cart.corr_cart).statsEqual(0, 2883204, 0, 0, 0));
    wassert(actual(cart.conv_cart).statsEqual(0, 2717594, 100, 100, 100));

    CartLowris clow(cb->do_medium ? 512: 256, *cb, cart);
    wassert(actual(clow.ZLR_N_ELEMENTARY_PIXEL) == 4);
    wassert(actual(clow.ZLR_OFFSET) == 337);
    clow.creo_cart_z_lowris();
//     print_stats("clow", clow, cerr);
    wassert(actual(clow.z_out).statsEqual(0, 0, 1, 80.12, 249));
    wassert(actual(clow.qual_Z_1x1).statsEqual(0, 0, 1, 28.51, 97));
    wassert(actual(clow.quota_1x1).statsEqual(0, 0, 128, 128, 128));
    wassert(actual(clow.dato_corr_1x1).statsEqual(0, 65536, 0, 0, 0));
    wassert(actual(clow.elev_fin_1x1).statsEqual(0, 55917, 1, 1.16, 3));
    wassert(actual(clow.beam_blocking_1x1).statsEqual(0, 37379, 1, 45.34, 51));
    wassert(actual(clow.top_1x1).statsEqual(0, 42201, 2, 13.87, 57));
    wassert(actual(clow.neve_1x1).statsEqual(0, 0, 1, 1, 1));
    wassert(actual(clow.corr_1x1).statsEqual(0, 65536, 0, 0, 0));
    wassert(actual(clow.conv_1x1).statsEqual(0, 59430, 100, 100, 100));

    // TODO: scrivo_out_file_bin
}

template<> template<>
void to::test<9>()
{
	LOG_CATEGORY("Test");
	LOG_INFO ("Start test 9");

    // versione BB che corrisponde al parametro algo_corto
    static const char* fname = "../testdata/DBP2_060220140140_GATTATICO";
    unsetwork();
    setenv("FIRST_LEVEL_FILE", "../dati/FIRST_LEVEL_corto+medio_GAT_INV_2011", 1);
    printwork();

    CBTest test("GAT", true, 1024);
    test.read_sp20(fname, false);
    auto cb = test.make_cumbac();
    cb->do_medium= true;
    cb->do_readStaticMap=true;
//    cb->do_zlr_media=true; 
    // FIXME: we don't compute VPR here, does it make sense to write vpr heating?
    cb->assets.write_vpr_heating(0);

    cb->declutter_anaprop();
    VolumeStats stats;
    cb->volume.compute_stats(stats);
//    stats.print(stdout);
    wassert(actual(stats.count_zeros[0]) == 165040);
    wassert(actual(stats.count_zeros[1]) == 207601);
    wassert(actual(stats.count_zeros[2]) == 161283);
    wassert(actual(stats.count_zeros[3]) == 108238);
    wassert(actual(stats.count_zeros[4]) == 117678);
    wassert(actual(stats.count_ones[0]) ==    0);
    wassert(actual(stats.count_ones[1]) ==    0);
    wassert(actual(stats.count_ones[2]) ==    0);
    wassert(actual(stats.count_ones[3]) ==    0);
    wassert(actual(stats.count_ones[4]) ==    0);
    wassert(actual(stats.count_others[0]) == 174560);
    wassert(actual(stats.count_others[1]) == 131999);
    wassert(actual(stats.count_others[2]) == 102717);
    wassert(actual(stats.count_others[3]) ==  89362);
    wassert(actual(stats.count_others[4]) ==  79922);
    wassert(actual(stats.sum_others[0]) == 24480412);
    wassert(actual(stats.sum_others[1]) == 17301958);
    wassert(actual(stats.sum_others[2]) == 12929034);
    wassert(actual(stats.sum_others[3]) == 10855770);
    wassert(actual(stats.sum_others[4]) ==  9385431);
    //print_stats("cb", *cb, cerr);

    Cart cart(cb->volume.max_beam_size());
    wassert(actual(cart.max_bin) == 849);
    cart.creo_cart(*cb);
//    print_stats("cart", cart, cerr);
    wassert(actual(cart.cart).statsEqual(0, 624947, 1, 46.62, 255));
    wassert(actual(cart.cartm).statsEqual(0, 2883204, 0, 0, 0));
    wassert(actual(cart.topxy).statsEqual(0, 2511651, 1, 16.37, 76));
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
//     print_stats("clow", clow, cerr);
    wassert(actual(clow.z_out).statsEqual(0, 120384, 1, 57.36, 255));
    wassert(actual(clow.qual_Z_1x1).statsEqual(0, 262144, 0, 0, 0));
    wassert(actual(clow.quota_1x1).statsEqual(0, 262144, 0, 0, 0));
    wassert(actual(clow.dato_corr_1x1).statsEqual(0, 262144, 0, 0, 0));
    wassert(actual(clow.elev_fin_1x1).statsEqual(0, 262144, 0, 0, 0));
    wassert(actual(clow.beam_blocking_1x1).statsEqual(0, 262144, 0, 0, 0));
    wassert(actual(clow.top_1x1).statsEqual(0, 231493, 3, 17.92, 76));
    wassert(actual(clow.neve_1x1).statsEqual(0, 262144, 0, 0, 0));
    wassert(actual(clow.corr_1x1).statsEqual(0, 262144, 0, 0, 0));
    wassert(actual(clow.conv_1x1).statsEqual(0, 262144, 0, 0, 0));
}

template<> template<>
void to::test<10>()
{
	LOG_CATEGORY("Test");
	LOG_INFO ("Start test 6");

    // versione BB che corrisponde al parametro algo_corto
    static const char* fname = "../testdata/2014-05-09-12-40-00.itgat.PVOL.0.h5";
    unsetwork();
    setenv("FIRST_LEVEL_FILE", "../dati/FIRST_LEVEL_corto_GAT_2006_INV", 1);
    setenv("DIR_OUT_PP_BLOC", "../testdata", 1);
    setenv("VPR_HEATING", "../testdata/vpr_heat_GAT", 1);
    setenv("FILE_ZERO_TERMICO", "../testdata/20140206/0termico.prev", 1);
    unlink("../testdata/vpr_heat_GAT");
    setenv("VPR0_FILE", "../testdata/ultimo_vpr", 1);
    unlink("../testdata/ultimo_vpr");
    setenv("FILE_T", "../testdata/temperature.txt", 1);
	printwork();

    CBTest test("GAT", false);
    test.read_odim(fname, true);
    auto cb = test.make_cumbac();
    cb->want_vpr();
    cb->do_quality = true;
    cb->do_beamblocking = true;
    cb->do_declutter = false;
    cb->do_class = false;
    cb->do_bloccorr = false;
    cb->do_readStaticMap=true;
    wassert(actual(cb->calcolo_vpr) != (void*)0);

    cb->declutter_anaprop();

    VolumeStats stats;
    cb->volume.compute_stats(stats);
    wassert(actual(cb->volume).statsEqual(MISSING_DB, 0, -31.5, -25.15, 78.51));

    wassert(actual(cb->beam_blocking).statsEqual(0, 15.83, 51));

    unsigned stats_size = cb->anaprop.grid_stats.size_az * cb->anaprop.grid_stats.size_beam;

    ArrayStats<unsigned> stat_anap_stats;
    stat_anap_stats.fill(cb->anaprop.grid_stats.stat_anap, stats_size);
    //stat_anap_stats.print();
    wassert(actual((unsigned)stat_anap_stats.all_missing).isfalse());
    wassert(actual((unsigned)stat_anap_stats.min) == 0);
    wassert(actual((unsigned)stat_anap_stats.max) == 0);

    ArrayStats<unsigned> stat_anap_tot_stats;
    stat_anap_tot_stats.fill(cb->anaprop.grid_stats.stat_tot, stats_size);
    //stat_anap_tot_stats.print();
    wassert(actual((unsigned)stat_anap_tot_stats.all_missing).isfalse());
    wassert(actual((unsigned)stat_anap_tot_stats.min) ==  0);
    wassert(actual((unsigned)stat_anap_tot_stats.max) == 0);

    ArrayStats<unsigned> stat_bloc_stats;
    stat_bloc_stats.fill(cb->anaprop.grid_stats.stat_bloc, stats_size);
//    stat_bloc_stats.print();
    wassert(actual((unsigned)stat_bloc_stats.all_missing).isfalse());
    wassert(actual((unsigned)stat_bloc_stats.min) == 0);
    wassert(actual((unsigned)stat_bloc_stats.max) == 0);

    ArrayStats<unsigned> stat_elev_stats;
    stat_elev_stats.fill(cb->anaprop.grid_stats.stat_elev, stats_size);
//    stat_elev_stats.print();
    wassert(actual((unsigned)stat_elev_stats.all_missing).isfalse());
    wassert(actual((unsigned)stat_elev_stats.min) == 0);
    wassert(actual((unsigned)stat_elev_stats.max) == 0);

    wassert(actual(cb->anaprop.dato_corrotto).statsEqual(0, 0., 0));

    cb->caratterizzo_volume();
    wassert(actual(cb->calcolo_vpr) != (void*)0);

//     print_stats("cb", *cb, cerr);
    wassert(actual(*cb->qual).statsEqual(1, 58.2, 99));
    wassert(actual(cb->top).statsEqual(0, 315507, 1, 9.38, 108));
    wassert(actual(cb->first_level).statsEqual(0, 0.33, 3));
    wassert(actual(cb->first_level_static).statsEqual(0, 0.33, 3));
    wassert(actual(cb->bb_first_level).statsEqual(0, 0.41, 2));
    wassert(actual(cb->beam_blocking).statsEqual(0, 15.83, 51));
    wassert(actual(cb->dem).statsEqual(-9999, 47.59, 4070.26));
    wassert(actual(cb->anaprop.quota).statsEqual(0, 0, 0));
    wassert(actual(cb->anaprop.dato_corrotto).statsEqual(0, 0, 0));
    wassert(actual(*cb->calcolo_vpr->flag_vpr).statsEqual(0, 0.94, 1));
    wassert(actual(cb->calcolo_vpr->corr_polar).statsEqual(0, 0, 0));
    wassert(actual(cb->calcolo_vpr->neve).statsEqual(0, 0, 0));

    Cart cart(cb->volume.max_beam_size());
    wassert(actual(cart.max_bin) == 849);
    cart.creo_cart(*cb);
//     print_stats("cart", cart, cerr);
    wassert(actual(cart.cart).statsEqual(0, 624947, 1, 13.3, 255));
    wassert(actual(cart.cartm).statsEqual(0, 2883204, 0, 0, 0));
    wassert(actual(cart.topxy).statsEqual(0, 2780214, 1, 15.87, 108));
    wassert(actual(cart.qual_Z_cart).statsEqual(0, 624947, 1, 18.35, 98));
    wassert(actual(cart.quota_cart).statsEqual(0, 2883204, 0, 0, 0));
    wassert(actual(cart.dato_corr_xy).statsEqual(0, 2883204, 0, 0, 0));
    wassert(actual(cart.beam_blocking_xy).statsEqual(0, 1768756, 1, 49.92, 51));
    wassert(actual(cart.elev_fin_xy).statsEqual(0, 2883204, 0, 0, 0));
    wassert(actual(cart.neve_cart).statsEqual(0, 624947, 1, 1, 1));
    wassert(actual(cart.corr_cart).statsEqual(0, 2883204, 0, 0, 0));
    wassert(actual(cart.conv_cart).statsEqual(0, 2883204, 0, 0, 0));

    CartLowris clow(cb->do_medium ? 512: 256, *cb, cart);
    wassert(actual(clow.ZLR_N_ELEMENTARY_PIXEL) == 4);
    wassert(actual(clow.ZLR_OFFSET) == 337);
    clow.creo_cart_z_lowris();
//     print_stats("clow", clow, cerr);
    wassert(actual(clow.z_out).statsEqual(0, 0, 1, 32.63, 255));
    wassert(actual(clow.qual_Z_1x1).statsEqual(0, 0, 1, 30.37, 92));
    wassert(actual(clow.quota_1x1).statsEqual(0, 0, 128, 128, 128));
    wassert(actual(clow.dato_corr_1x1).statsEqual(0, 65536, 0, 0, 0));
    wassert(actual(clow.elev_fin_1x1).statsEqual(0, 65536, 0, 0, 0));
    wassert(actual(clow.beam_blocking_1x1).statsEqual(0, 43347, 1, 50.13, 51));
    wassert(actual(clow.top_1x1).statsEqual(0, 57731, 1, 12.99, 108));
    wassert(actual(clow.neve_1x1).statsEqual(0, 0, 1, 1, 1));
    wassert(actual(clow.corr_1x1).statsEqual(0, 65536, 0, 0, 0));
    wassert(actual(clow.conv_1x1).statsEqual(0, 65536, 0, 0, 0));
}

template<> template<>
void to::test<11>()
{
	LOG_CATEGORY("Test");
	LOG_INFO ("Start test 6");

    // versione BB che corrisponde al parametro algo_corto
    static const char* fname = "../esplosione/2014-03-01-09-15-00.itgat.PVOL.0.h5";
    unsetwork();
    setenv("DIR_OUT_PP_BLOC", "../esplosione", 1);
    setenv("VPR0_FILE"		, "../esplosione/vpr_GAT", 1);
    setenv("LAST_VPR"    	, "../esplosione/last_vpr_GAT", 1);
    setenv("VPR_HMAX"    	, "../esplosione/vpr_hmax_GAT",1); 
    setenv("VPR_HEATING" 	, "../esplosione/vpr_heat_GAT",1);
    setenv("LOG_VPR"		, "../esplosione/log_VPR_${SITO}",1);
    setenv("TEST_VPR"		, "../esplosione/test_vpr",1);
    setenv("FILE_T"		, "../esplosione/temperature.txt",1);
    setenv("FIRST_LEVEL_FILE", "../dati/FIRST_LEVEL_corto_GAT_PRI-EST_2011", 1);
    setenv("FILE_ZERO_TERMICO"	, "../esplosione/0termico.prev", 1);
    setenv("VPR_ARCH"           , "../esplosione/201403010915_vpr_GAT",1);
	printwork();

    CBTest test("GAT", false);
    test.read_odim(fname, true);
    auto cb = test.make_cumbac();
    cb->want_vpr();
    cb->do_readStaticMap=true;
    cb->do_devel=true;
    cb->do_beamblocking = true;
    cb->do_quality = true;
    cb->do_bloccorr = true;
    cb->do_class = true;

    cb->declutter_anaprop();
    cb->caratterizzo_volume();
    wassert(actual(cb->calcolo_vpr) != (void*)0);
    cb->calcolo_vpr->classifica_rain();
    cb->calcolo_vpr->esegui_tutto();
    VolumeStats stats;
    cb->volume.compute_stats(stats);
}


template<> template<>
void to::test<12>()
{
	LOG_CATEGORY("Test");
	LOG_INFO ("Start test 12");

    // versione BB che corrisponde al parametro algo_corto
    static const char* fname = "vpr/2014-03-01-01-35-00.itgat.PVOL.0.h5";
    unsetwork();
    setenv("DIR_OUT_PP_BLOC", "vpr", 1);
    setenv("VPR0_FILE"		, "vpr/vpr_GAT", 1);
    setenv("LAST_VPR"    	, "vpr/last_vpr_GAT", 1);
    setenv("VPR_HMAX"    	, "vpr/vpr_hmax_GAT",1); 
    setenv("VPR_HEATING" 	, "vpr/vpr_heat_GAT",1);
    setenv("LOG_VPR"		, "vpr/log_VPR_${SITO}",1);
    setenv("TEST_VPR"		, "vpr/test_vpr",1);
    setenv("FILE_T"		, "vpr/temperature.txt",1);
//setenv("FIRST_LEVEL_FILE", "../dati/FIRST_LEVEL_corto_GAT_INV_2011_15el", 1);
    setenv("FIRST_LEVEL_FILE", "../dati/FIRST_LEVEL_corto+medio_GAT_INV_2011", 1);
    setenv("FILE_ZERO_TERMICO"	, "vpr/0termico.prev", 1);
    setenv("OUTPUT_Z_LOWRIS_DIR","vpr",1);
    setenv("DIR_QUALITY","vpr",1);
    setenv("DIR_DEBUG","vpr",1);

	printwork();

    CBTest test("GAT", false);
    test.read_odim(fname, true);
    auto cb = test.make_cumbac();
    cb->want_vpr();
    cb->do_quality = true;
    cb->do_beamblocking = true;
    cb->do_declutter = false;
    cb->do_class = true;
    cb->do_bloccorr = true;
    cb->do_devel= true;
    cb->do_readStaticMap=true;
    
    cb->declutter_anaprop();
    if (cb->do_quality){
       cb->caratterizzo_volume();
       wassert(actual(cb->calcolo_vpr) != (void*)0);
       if (cb->do_class) cb->calcolo_vpr->classifica_rain();
       cb->calcolo_vpr->esegui_tutto();
    }
 //   if (cb->do_class)
 //       cb->conversione_convettiva();

    /*--------------------------------------------------
      | conversione di coordinate da polare a cartesiana |
      --------------------------------------------------*/
    LOG_INFO("Creazione Matrice Cartesiana");
    Cart cart_maker(cb->volume.max_beam_size());
    cart_maker.creo_cart(*cb);
    cb->assets.write_gdal_image(cart_maker.cart,"DIR_DEBUG","cart","PNG" );

    //-------------------Se definita Z_LOWRIS creo matrice 1X1  ZLR  stampo e stampo coeff MP (serve?)------------------

    LOG_INFO("Estrazione Precipitazione 1X1");
    CartLowris cart_low(cb->do_medium ? 512: 256, *cb, cart_maker);
    cart_low.creo_cart_z_lowris();

    LOG_INFO("Scrittura File Precipitazione 1X1\n");
    cart_low.write_out(*cb,cb->assets);
    cb->assets.write_gdal_image(cart_low.z_out,"DIR_DEBUG","ZLR","PNG" );
}

}
