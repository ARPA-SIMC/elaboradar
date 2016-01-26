#include "radarelab/utils/tests.h"
#include <radarelab/logging.h>
#include <radarelab/image.h>
#include "cum_bac.h"
#include "config.h"
#include <cstdio>
#include <cstdlib>
#include <vector>
#include "site.h"
#include "setwork.h"
#include "cartproducts.h"
#include "test-utils.h"
#include <unistd.h>

using namespace radarelab::utils::tests;
using namespace radarelab;
using namespace elaboradar;
using namespace testradar;
using namespace std;

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

class Tests : public TestCase
{
    using TestCase::TestCase;

    void register_tests() override;
} test("process");

void Tests::register_tests() {

add_method("elabora_all_true_caratterizzo", []() {
    // Test elabora_dato, con tutti i do_* a true
    static const char* fname = "../testdata/DBP2_070120141530_GATTATICO";

    unsetwork();
    setenv("FIRST_LEVEL_FILE", "../dati/FIRST_LEVEL_corto_GAT_2006_INV", 1);
    setenv("DIR_OUT_PP_BLOC", "../testdata", 1);
    setenv("FILE_T", "../testdata/temperature.txt", 1);
    setenv("LAST_VPR","../testdata/last_vpr",1);
    setenv("FILE_ZERO_TERMICO","../testdata/zero_termico.txt",1);
    setenv("FILE_DEM_GAT", "../testdata/dem_Gatta.txt", 1);
    setenv("FILE_DEM_SPC", "../testdata/dem_SanPi.txt", 1);

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
    wassert(actual(stats.count_zeros[0]) == 145426u);
    wassert(actual(stats.count_zeros[1]) == 184085u);
    wassert(actual(stats.count_zeros[2]) == 193145u);
    wassert(actual(stats.count_zeros[3]) == 196213u);
    wassert(actual(stats.count_zeros[4]) == 196085u);
    wassert(actual(stats.count_zeros[5]) == 196080u);
    wassert(actual(stats.count_ones[0]) == 104u);
    wassert(actual(stats.count_ones[1]) == 204u);
    wassert(actual(stats.count_ones[2]) == 104u);
    wassert(actual(stats.count_ones[3]) ==  58u);
    wassert(actual(stats.count_ones[4]) ==  61u);
    wassert(actual(stats.count_ones[5]) ==  38u);
    wassert(actual(stats.count_others[0]) == 52070u);
    wassert(actual(stats.count_others[1]) == 13311u);
    wassert(actual(stats.count_others[2]) ==  4351u);
    wassert(actual(stats.count_others[3]) ==  1329u);
    wassert(actual(stats.count_others[4]) ==  1454u);
    wassert(actual(stats.count_others[5]) ==  1482u);
    wassert(actual(stats.sum_others[0]) == 4759195u);
    wassert(actual(stats.sum_others[1]) ==  914997u);
    wassert(actual(stats.sum_others[2]) ==  257459u);
    wassert(actual(stats.sum_others[3]) ==   46349u);
    wassert(actual(stats.sum_others[4]) ==   78749u);
    wassert(actual(stats.sum_others[5]) ==   90563u);
    cb->caratterizzo_volume();

//    print_stats("cb", *cb, cerr);
    wassert(actual(*cb->qual).statsEqual(1, 61.81, 99));
    wassert(actual(cb->top).statsEqual(0, 189090, 1, 7.10, 36));
    wassert(actual(cb->first_level).statsEqual(0, 0.86, 3));
    wassert(actual(cb->first_level_static).statsEqual(0, 0.54, 3));
    wassert(actual(cb->bb_first_level).statsEqual(0, 0.49, 2));
    wassert(actual(cb->beam_blocking).statsEqual(0, 4.72, 49));
    wassert(actual(cb->dem).statsEqual(-9999, 47.59, 4070.26));
    wassert(actual(cb->anaprop.quota).statsEqual(0, 0.00, 0));
    wassert(actual(cb->anaprop.dato_corrotto).statsEqual(0, 0.00, 0));
    wassert(actual(*cb->calcolo_vpr->flag_vpr).statsEqual(0, 0.99, 1));
    wassert(actual(cb->calcolo_vpr->corr_polar).statsEqual(0, 0.00, 0));
    wassert(actual(cb->calcolo_vpr->neve).statsEqual(0, 0.00, 0));
});

add_method("elabora_all_true", []() {
    // Test elabora_dato, con tutti i do_* a true
    static const char* fname = "../testdata/DBP2_070120141530_GATTATICO";

    unsetwork();
    setenv("FIRST_LEVEL_FILE", "../dati/FIRST_LEVEL_corto_GAT_2006_INV", 1);
    setenv("DIR_OUT_PP_BLOC", "../testdata", 1);
    setenv("FILE_DEM_GAT", "../testdata/dem_Gatta.txt", 1);
    setenv("FILE_DEM_SPC", "../testdata/dem_SanPi.txt", 1);

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
    wassert(actual(stats.count_zeros[0]) == 145426u);
    wassert(actual(stats.count_zeros[1]) == 184085u);
    wassert(actual(stats.count_zeros[2]) == 193145u);
    wassert(actual(stats.count_zeros[3]) == 196213u);
    wassert(actual(stats.count_zeros[4]) == 196085u);
    wassert(actual(stats.count_zeros[5]) == 196080u);
    wassert(actual(stats.count_ones[0]) == 104u);
    wassert(actual(stats.count_ones[1]) == 204u);
    wassert(actual(stats.count_ones[2]) == 104u);
    wassert(actual(stats.count_ones[3]) ==  58u);
    wassert(actual(stats.count_ones[4]) ==  61u);
    wassert(actual(stats.count_ones[5]) ==  38u);
    wassert(actual(stats.count_others[0]) == 52070u);
    wassert(actual(stats.count_others[1]) == 13311u);
    wassert(actual(stats.count_others[2]) ==  4351u);
    wassert(actual(stats.count_others[3]) ==  1329u);
    wassert(actual(stats.count_others[4]) ==  1454u);
    wassert(actual(stats.count_others[5]) ==  1482u);
    wassert(actual(stats.sum_others[0]) == 4759195u);
    wassert(actual(stats.sum_others[1]) ==  914997u);
    wassert(actual(stats.sum_others[2]) ==  257459u);
    wassert(actual(stats.sum_others[3]) ==   46349u);
    wassert(actual(stats.sum_others[4]) ==   78749u);
    wassert(actual(stats.sum_others[5]) ==   90563u);
});

add_method("bb_algo_corto", []() {
    // versione BB che corrisponde al parametro algo_corto
    static const char* fname = "../testdata/DBP2_070120141530_GATTATICO";
    unsetwork();
    setenv("FIRST_LEVEL_FILE", "../dati/FIRST_LEVEL_corto_GAT_2006_INV", 1);
    setenv("DIR_OUT_PP_BLOC", "../testdata", 1);
    setenv("VPR_HEATING", "../testdata/vpr_heat_GAT", 1);
    unlink("../testdata/vpr_heat_GAT");
    setenv("FILE_T", "../testdata/temperature.txt", 1);
    setenv("FILE_DEM_GAT", "../testdata/dem_Gatta.txt", 1);
    setenv("FILE_DEM_SPC", "../testdata/dem_SanPi.txt", 1);

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
    wassert(actual(*cb->qual).statsEqual(1, 60.03, 99));
    wassert(actual(cb->top).statsEqual(0, 189090, 1, 7.10, 36));
    wassert(actual(cb->first_level).statsEqual(0, 0.54, 3));
    wassert(actual(cb->first_level_static).statsEqual(0, 0.54, 3));
    wassert(actual(cb->bb_first_level).statsEqual(0, 0.49, 2));
    wassert(actual(cb->beam_blocking).statsEqual(0, 10.44, 51));
    wassert(actual(cb->dem).statsEqual(-9999, 47.59, 4070.26));
    wassert(actual(cb->anaprop.quota).statsEqual(0, 0.00, 0));
    wassert(actual(cb->anaprop.dato_corrotto).statsEqual(0, 0.00, 0));
    wassert(actual(*cb->calcolo_vpr->flag_vpr).statsEqual(0, 0.96, 1));
    wassert(actual(cb->calcolo_vpr->corr_polar).statsEqual(0, 0.00, 0));
    wassert(actual(cb->calcolo_vpr->neve).statsEqual(0, 0.00, 0));

    Cart cart(cb->volume.max_beam_size());
    wassert(actual(cart.max_bin) == 494);
    cart.creo_cart(*cb);
    //setenv("DIR_DEBUG", "../testdata/", 1);
    //cart.write_out(*cb, cb->assets);
//     print_stats("cart", cart, cout);
    wassert(actual(cart.cart).statsEqual(0, 213127, 1, 23.95, 251));
    wassert(actual(cart.cartm).statsEqual(0.00, 976144, 0.00, 0.00, 0.00));
    wassert(actual(cart.topxy).statsEqual(0, 941551, 1, 10.84, 36));
    wassert(actual(cart.qual_Z_cart).statsEqual(0, 213127, 1, 36.94, 98));
    wassert(actual(cart.quota_cart).statsEqual(0, 976144, 0, 0.00, 0));
    wassert(actual(cart.dato_corr_xy).statsEqual(0, 976144, 0, 0.00, 0));
    wassert(actual(cart.beam_blocking_xy).statsEqual(0, 693556, 1, 47.06, 51));
    wassert(actual(cart.elev_fin_xy).statsEqual(0, 976144, 0, 0.00, 0));
    wassert(actual(cart.neve_cart).statsEqual(0, 213127, 1, 1.00, 1));
    wassert(actual(cart.corr_cart).statsEqual(0, 976144, 0, 0.00, 0));
    wassert(actual(cart.conv_cart).statsEqual(0, 976144, 0, 0.00, 0));

    CartLowris clow(cb->do_medium ? 512: 256, *cb, cart);
    wassert(actual(clow.CART_DIM_ZLR) == 256);
    wassert(actual(clow.ZLR_N_ELEMENTARY_PIXEL) == 4);
    wassert(actual(clow.ZLR_OFFSET) == -18);
    clow.creo_cart_z_lowris();
//     print_stats("clow", clow, cerr);
    wassert(actual(clow.z_out).statsEqual(0, 17489, 1, 36.57, 251));
    wassert(actual(clow.qual_Z_1x1).statsEqual(0, 17489, 1, 37.29, 97));
    wassert(actual(clow.quota_1x1).statsEqual(0, 0, 128, 128.00, 128));
    wassert(actual(clow.dato_corr_1x1).statsEqual(0, 65536, 0, 0.00, 0));
    wassert(actual(clow.elev_fin_1x1).statsEqual(0, 65536, 0, 0.00, 0));
    wassert(actual(clow.beam_blocking_1x1).statsEqual(0, 48089, 1, 47.04, 51));
    wassert(actual(clow.top_1x1).statsEqual(0, 61004, 1, 10.53, 36));
    wassert(actual(clow.neve_1x1).statsEqual(0, 17489, 1, 1.00, 1));
    wassert(actual(clow.corr_1x1).statsEqual(0, 65536, 0, 0.00, 0));
    wassert(actual(clow.conv_1x1).statsEqual(0, 65536, 0, 0.00, 0));

    CartProducts products(cb->volume, 256, 4);
    cb->generate_maps(products, true);
    wassert(actual(products.scaled.image_offset) == -18);
    // write_image(cb->volume[0], "/tmp/zout-raw.png", "png");
    // write_image(clow.z_out, "/tmp/zout-old.png", "png");
    // write_image(products.z_out, "/tmp/zout-new.png", "png");
//     print_stats("products", products, cerr);
    wassert(actual(products.z_out).statsEqual(0, 17240, 1, 36.72, 251));
    wassert(actual(products.qual_Z_1x1).statsEqual(0, 17240, 1, 37.21, 98));
    wassert(actual(products.quota_1x1).statsEqual(0, 0, 128, 128.00, 128));
    wassert(actual(products.dato_corr_1x1).statsEqual(0, 65536, 0, 0.00, 0));
    wassert(actual(products.elev_fin_1x1).statsEqual(0, 65536, 0, 0.00, 0));
    wassert(actual(products.beam_blocking_1x1).statsEqual(0, 47922, 1, 47.13, 51));
    wassert(actual(products.top_1x1).statsEqual(0, 60943, 1, 10.56, 36));
    wassert(actual(products.neve_1x1).statsEqual(0, 17240, 1, 1.00, 1));
    wassert(actual(products.corr_1x1).statsEqual(0, 65536, 0, 0.00, 0));
    wassert(actual(products.conv_1x1).statsEqual(0, 65536, 0, 0.00, 0));


    // Try again, but with averages
    cb->do_zlr_media = true;

    Cart cart1(cb->volume.max_beam_size());
    cart1.creo_cart(*cb);
//     print_stats("cart1", cart1, cerr);
    wassert(actual(cart1.cart).statsEqual(0, 213127, 1, 23.95, 251));
    wassert(actual(cart1.cartm).statsEqual(0.00, 213127, 0.01, 264.48, 749048.27));
    wassert(actual(cart1.topxy).statsEqual(0, 941551, 1, 10.84, 36));
    wassert(actual(cart1.qual_Z_cart).statsEqual(0, 213127, 1, 36.94, 98));
    wassert(actual(cart1.quota_cart).statsEqual(0, 976144, 0, 0.00, 0));
    wassert(actual(cart1.dato_corr_xy).statsEqual(0, 976144, 0, 0.00, 0));
    wassert(actual(cart1.beam_blocking_xy).statsEqual(0, 693556, 1, 47.06, 51));
    wassert(actual(cart1.elev_fin_xy).statsEqual(0, 976144, 0, 0.00, 0));
    wassert(actual(cart1.neve_cart).statsEqual(0, 213127, 1, 1.00, 1));
    wassert(actual(cart1.corr_cart).statsEqual(0, 976144, 0, 0.00, 0));
    wassert(actual(cart1.conv_cart).statsEqual(0, 976144, 0, 0.00, 0));

    CartLowris clow1(cb->do_medium ? 512: 256, *cb, cart1);
    wassert(actual(clow1.CART_DIM_ZLR) == 256);
    wassert(actual(clow1.ZLR_N_ELEMENTARY_PIXEL) == 4);
    wassert(actual(clow1.ZLR_OFFSET) == -18);
    clow1.creo_cart_z_lowris();
//    print_stats("clow1", clow1, cerr);
    wassert(actual(clow1.z_out).statsEqual(0, 17489, 1, 29.85, 235));
    wassert(actual(clow1.qual_Z_1x1).statsEqual(0, 17489, 1, 37.29, 97));
    wassert(actual(clow1.quota_1x1).statsEqual(0, 0, 128, 128.00, 128));
    wassert(actual(clow1.dato_corr_1x1).statsEqual(0, 65536, 0, 0.00, 0));
    wassert(actual(clow1.elev_fin_1x1).statsEqual(0, 65536, 0, 0.00, 0));
    wassert(actual(clow1.beam_blocking_1x1).statsEqual(0, 48089, 1, 47.04, 51));
    wassert(actual(clow1.top_1x1).statsEqual(0, 61004, 1, 10.53, 36));
    wassert(actual(clow1.neve_1x1).statsEqual(0, 17489, 1, 1.00, 1));
    wassert(actual(clow1.corr_1x1).statsEqual(0, 65536, 0, 0.00, 0));
    wassert(actual(clow1.conv_1x1).statsEqual(0, 65536, 0, 0.00, 0));

    CartProducts products1(cb->volume, 256, 4);
    cb->generate_maps(products1, true);
    wassert(actual(products1.scaled.image_offset) == -18);
    // write_image(cb->volume[0], "/tmp/zout-raw.png", "png");
    // write_image(clow1.z_out, "/tmp/zout-old.png", "png");
    // write_image(products1.z_out, "/tmp/zout-new.png", "png");
//     print_stats("products", products, cerr);
    wassert(actual(products.z_out).statsEqual(0, 17240, 1, 36.72, 251));
    wassert(actual(products.qual_Z_1x1).statsEqual(0, 17240, 1, 37.21, 98));
    wassert(actual(products.quota_1x1).statsEqual(0, 0, 128, 128.00, 128));
    wassert(actual(products.dato_corr_1x1).statsEqual(0, 65536, 0, 0.00, 0));
    wassert(actual(products.elev_fin_1x1).statsEqual(0, 65536, 0, 0.00, 0));
    wassert(actual(products.beam_blocking_1x1).statsEqual(0, 47922, 1, 47.13, 51));
    wassert(actual(products.top_1x1).statsEqual(0, 60943, 1, 10.56, 36));
    wassert(actual(products.neve_1x1).statsEqual(0, 17240, 1, 1.00, 1));
    wassert(actual(products.corr_1x1).statsEqual(0, 65536, 0, 0.00, 0));
    wassert(actual(products.conv_1x1).statsEqual(0, 65536, 0, 0.00, 0));
});

add_method("bb_algo_corto_dev", []() {
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
    setenv("FILE_DEM_GAT", "../testdata/dem_Gatta.txt", 1);
    setenv("FILE_DEM_SPC", "../testdata/dem_SanPi.txt", 1);
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
    wassert(actual(cb->top).statsEqual(0, 189090, 1, 7.10, 36));
    wassert(actual(cb->first_level).statsEqual(0, 0.54, 3));
    wassert(actual(cb->first_level_static).statsEqual(0, 0.54, 3));
    wassert(actual(cb->bb_first_level).statsEqual(0, 0.49, 2));
    wassert(actual(cb->beam_blocking).statsEqual(0, 10.44, 51));
    wassert(actual(cb->dem).statsEqual(-9999, 47.59, 4070.26));
    wassert(actual(cb->anaprop.quota).statsEqual(0, 0.00, 0));
    wassert(actual(cb->anaprop.dato_corrotto).statsEqual(0, 0.00, 0));
    wassert(actual(*cb->calcolo_vpr->flag_vpr).statsEqual(0, 0.70, 1));
    wassert(actual(cb->calcolo_vpr->corr_polar).statsEqual(0, 0.00, 0));
    wassert(actual(cb->calcolo_vpr->neve).statsEqual(0, 0.00, 0));

    //ier = cb->calcolo_vpr->corr_vpr();
    //wassert(actual(ier) == 1);

    // TODO: cb->stampa_vpr()

    Cart cart(cb->volume.max_beam_size());
    wassert(actual(cart.max_bin) == 494);
    cart.creo_cart(*cb);
//     print_stats("cart", cart, cerr);
    wassert(actual(cart.cart).statsEqual(0, 213127, 1, 14.12, 227));
    wassert(actual(cart.cartm).statsEqual(0.00, 976144, 0.00, 0.00, 0.00));
    wassert(actual(cart.topxy).statsEqual(0, 945040, 1, 11.02, 36));
    wassert(actual(cart.qual_Z_cart).statsEqual(0, 213127, 1, 36.52, 98));
    wassert(actual(cart.quota_cart).statsEqual(0, 976144, 0, 0.00, 0));
    wassert(actual(cart.dato_corr_xy).statsEqual(0, 976144, 0, 0.00, 0));
    wassert(actual(cart.beam_blocking_xy).statsEqual(0, 689773, 1, 47.01, 51));
    wassert(actual(cart.elev_fin_xy).statsEqual(0, 764742, 1, 1.15, 3));
    wassert(actual(cart.neve_cart).statsEqual(0, 213127, 1, 1.00, 1));
    wassert(actual(cart.corr_cart).statsEqual(0, 976144, 0, 0.00, 0));
    wassert(actual(cart.conv_cart).statsEqual(0, 976144, 0, 0.00, 0));

    CartLowris clow(cb->do_medium ? 512: 256, *cb, cart);
    wassert(actual(clow.ZLR_N_ELEMENTARY_PIXEL) == 4);
    wassert(actual(clow.ZLR_OFFSET) == -18);
    clow.creo_cart_z_lowris();
//    print_stats("clow", clow, cerr);
    wassert(actual(clow.z_out).statsEqual(0, 17489, 1, 22.69, 227));
    wassert(actual(clow.qual_Z_1x1).statsEqual(0, 17489, 1, 35.78, 98));
    wassert(actual(clow.quota_1x1).statsEqual(0, 0, 128, 128.00, 128));
    wassert(actual(clow.dato_corr_1x1).statsEqual(0, 65536, 0, 0.00, 0));
    wassert(actual(clow.elev_fin_1x1).statsEqual(0, 53710, 1, 1.14, 3));
    wassert(actual(clow.beam_blocking_1x1).statsEqual(0, 46906, 1, 46.97, 51));
    wassert(actual(clow.top_1x1).statsEqual(0, 62679, 1, 11.32, 36));
    wassert(actual(clow.neve_1x1).statsEqual(0, 17489, 1, 1.00, 1));
    wassert(actual(clow.corr_1x1).statsEqual(0, 65536, 0, 0.00, 0));
    wassert(actual(clow.conv_1x1).statsEqual(0, 65536, 0, 0.00, 0));
});

add_method("bb_vpr_class_algo_corto_dev", []() {
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
    setenv("FILE_DEM_GAT", "../testdata/dem_Gatta.txt", 1);
    setenv("FILE_DEM_SPC", "../testdata/dem_SanPi.txt", 1);
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
    wassert(actual(cb->top).statsEqual(0, 189090, 1, 7.10, 36));
    wassert(actual(cb->first_level).statsEqual(0, 0.54, 3));
    wassert(actual(cb->first_level_static).statsEqual(0, 0.54, 3));
    wassert(actual(cb->bb_first_level).statsEqual(0, 0.49, 2));
    wassert(actual(cb->beam_blocking).statsEqual(0, 10.44, 51));
    wassert(actual(cb->dem).statsEqual(-9999, 47.59, 4070.26));
    wassert(actual(cb->anaprop.quota).statsEqual(0, 0.00, 0));
    wassert(actual(cb->anaprop.dato_corrotto).statsEqual(0, 0.00, 0));
    wassert(actual(*cb->calcolo_vpr->flag_vpr).statsEqual(0, 0.70, 1));
    wassert(actual(cb->calcolo_vpr->corr_polar).statsEqual(0, 0.00, 0));
    wassert(actual(cb->calcolo_vpr->neve).statsEqual(0, 0.00, 0));

    // TODO: cb->stampa_vpr()

    cb->conversione_convettiva();

    Cart cart(cb->volume.max_beam_size());
    wassert(actual(cart.max_bin) == 494);
    cart.creo_cart(*cb);
//    print_stats("cart", cart, cerr);
    wassert(actual(cart.cart).statsEqual(0, 213127, 1, 14.12, 227));
    wassert(actual(cart.cartm).statsEqual(0.00, 976144, 0.00, 0.00, 0.00));
    wassert(actual(cart.topxy).statsEqual(0, 945040, 1, 11.02, 36));
    wassert(actual(cart.qual_Z_cart).statsEqual(0, 213127, 1, 36.52, 98));
    wassert(actual(cart.quota_cart).statsEqual(0, 976144, 0, 0.00, 0));
    wassert(actual(cart.dato_corr_xy).statsEqual(0, 976144, 0, 0.00, 0));
    wassert(actual(cart.beam_blocking_xy).statsEqual(0, 689773, 1, 47.01, 51));
    wassert(actual(cart.elev_fin_xy).statsEqual(0, 764742, 1, 1.15, 3));
    wassert(actual(cart.neve_cart).statsEqual(0, 213127, 1, 1.00, 1));
    wassert(actual(cart.corr_cart).statsEqual(0, 976144, 0, 0.00, 0));
    wassert(actual(cart.conv_cart).statsEqual(0, 976144, 0, 0.00, 0));

    CartLowris clow(cb->do_medium ? 512: 256, *cb, cart);
    wassert(actual(clow.ZLR_N_ELEMENTARY_PIXEL) == 4);
    wassert(actual(clow.ZLR_OFFSET) == -18);
    clow.creo_cart_z_lowris();
//     print_stats("clow", clow, cerr);
    wassert(actual(clow.z_out).statsEqual(0, 17489, 1, 22.69, 227));
    wassert(actual(clow.qual_Z_1x1).statsEqual(0, 17489, 1, 35.78, 98));
    wassert(actual(clow.quota_1x1).statsEqual(0, 0, 128, 128.00, 128));
    wassert(actual(clow.dato_corr_1x1).statsEqual(0, 65536, 0, 0.00, 0));
    wassert(actual(clow.elev_fin_1x1).statsEqual(0, 53710, 1, 1.14, 3));
    wassert(actual(clow.beam_blocking_1x1).statsEqual(0, 46906, 1, 46.97, 51));
    wassert(actual(clow.top_1x1).statsEqual(0, 62679, 1, 11.32, 36));
    wassert(actual(clow.neve_1x1).statsEqual(0, 17489, 1, 1.00, 1));
    wassert(actual(clow.corr_1x1).statsEqual(0, 65536, 0, 0.00, 0));
    wassert(actual(clow.conv_1x1).statsEqual(0, 65536, 0, 0.00, 0));
    // TODO: scrivo_out_file_bin
});

add_method("bb_algo_corto", []() {
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
    setenv("FILE_DEM_GAT", "../testdata/dem_Gatta.txt", 1);
    setenv("FILE_DEM_SPC", "../testdata/dem_SanPi.txt", 1);
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
//    stats.print(stdout);
    wassert(actual(stats.count_zeros[0]) == 172245);
    wassert(actual(stats.count_zeros[1]) == 210954);
    wassert(actual(stats.count_zeros[2]) == 163036);
    wassert(actual(stats.count_zeros[3]) == 109587);
    wassert(actual(stats.count_zeros[4]) == 118624);
    wassert(actual(stats.count_ones[0]) == 0);
    wassert(actual(stats.count_ones[1]) == 0);
    wassert(actual(stats.count_ones[2]) == 0);
    wassert(actual(stats.count_ones[3]) == 0);
    wassert(actual(stats.count_ones[4]) == 0);
    wassert(actual(stats.count_others[0]) == 167355);
    wassert(actual(stats.count_others[1]) == 128646);
    wassert(actual(stats.count_others[2]) == 100964);
    wassert(actual(stats.count_others[3]) ==  88013);
    wassert(actual(stats.count_others[4]) ==  78976);
    wassert(actual(stats.sum_others[0]) == 23175278);
    wassert(actual(stats.sum_others[1]) == 16731796);
    wassert(actual(stats.sum_others[2]) == 12611046);
    wassert(actual(stats.sum_others[3]) == 10596468);
    wassert(actual(stats.sum_others[4]) ==  9176701);

    wassert(actual(cb->beam_blocking).statsEqual(0, 17.1, 51));

    unsigned stats_size = cb->anaprop.grid_stats.size_az * cb->anaprop.grid_stats.size_beam;

    ArrayStats<unsigned> stat_anap_stats;
    stat_anap_stats.fill(cb->anaprop.grid_stats.stat_anap, stats_size);
//    stat_anap_stats.print();
    wassert(actual((unsigned)stat_anap_stats.all_missing).isfalse());
    wassert(actual((unsigned)stat_anap_stats.min) == 0);
    wassert(actual((unsigned)stat_anap_stats.max) == 0);

    ArrayStats<unsigned> stat_anap_tot_stats;
    stat_anap_tot_stats.fill(cb->anaprop.grid_stats.stat_tot, stats_size);
//    stat_anap_tot_stats.print();
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
    wassert(actual(*cb->qual).statsEqual(1, 56.26, 99));
    wassert(actual(cb->top).statsEqual(0, 227845, 1, 12.72, 76));
    wassert(actual(cb->first_level).statsEqual(0, 0.33, 3));
    wassert(actual(cb->first_level_static).statsEqual(0, 0.33, 3));
    wassert(actual(cb->bb_first_level).statsEqual(0, 0.48, 2));
    wassert(actual(cb->beam_blocking).statsEqual(0, 17.10, 51));
    wassert(actual(cb->dem).statsEqual(-9999, 47.59, 4070.26));
    wassert(actual(cb->anaprop.quota).statsEqual(0, 0.00, 0));
    wassert(actual(cb->anaprop.dato_corrotto).statsEqual(0, 0.00, 0));
    wassert(actual(*cb->calcolo_vpr->flag_vpr).statsEqual(0, 0.94, 1));
    wassert(actual(cb->calcolo_vpr->corr_polar).statsEqual(0, 0.00, 0));
    wassert(actual(cb->calcolo_vpr->neve).statsEqual(0, 0.00, 0));

    Cart cart(cb->volume.max_beam_size());
    wassert(actual(cart.max_bin) == 849);
    cart.creo_cart(*cb);
//     print_stats("cart", cart, cerr);
    wassert(actual(cart.cart).statsEqual(0, 624947, 1, 43.33, 255));
    wassert(actual(cart.cartm).statsEqual(0.00, 2883204, 0.00, 0.00, 0.00));
    wassert(actual(cart.topxy).statsEqual(0, 2539486, 1, 15.63, 76));
    wassert(actual(cart.qual_Z_cart).statsEqual(0, 624947, 1, 18.77, 98));
    wassert(actual(cart.quota_cart).statsEqual(0, 2883204, 0, 0.00, 0));
    wassert(actual(cart.dato_corr_xy).statsEqual(0, 2883204, 0, 0.00, 0));
    wassert(actual(cart.beam_blocking_xy).statsEqual(0, 1693879, 1, 45.13, 51));
    wassert(actual(cart.elev_fin_xy).statsEqual(0, 2883204, 0, 0.00, 0));
    wassert(actual(cart.neve_cart).statsEqual(0, 624947, 1, 1.00, 1));
    wassert(actual(cart.corr_cart).statsEqual(0, 2883204, 0, 0.00, 0));
    wassert(actual(cart.conv_cart).statsEqual(0, 2883204, 0, 0.00, 0));

    CartLowris clow(cb->do_medium ? 512: 256, *cb, cart);
    wassert(actual(clow.ZLR_N_ELEMENTARY_PIXEL) == 4);
    wassert(actual(clow.ZLR_OFFSET) == 337);
    clow.creo_cart_z_lowris();
//     print_stats("clow", clow, cerr);
    wassert(actual(clow.z_out).statsEqual(0, 0, 1, 85.06, 255));
    wassert(actual(clow.qual_Z_1x1).statsEqual(0, 0, 1, 28.35, 97));
    wassert(actual(clow.quota_1x1).statsEqual(0, 0, 128, 128.00, 128));
    wassert(actual(clow.dato_corr_1x1).statsEqual(0, 65536, 0, 0.00, 0));
    wassert(actual(clow.elev_fin_1x1).statsEqual(0, 65536, 0, 0.00, 0));
    wassert(actual(clow.beam_blocking_1x1).statsEqual(0, 39215, 1, 45.68, 51));
    wassert(actual(clow.top_1x1).statsEqual(0, 42404, 2, 13.40, 57));
    wassert(actual(clow.neve_1x1).statsEqual(0, 0, 1, 1.00, 1));
    wassert(actual(clow.corr_1x1).statsEqual(0, 65536, 0, 0.00, 0));
    wassert(actual(clow.conv_1x1).statsEqual(0, 65536, 0, 0.00, 0));
});

add_method("combina_profili", []() {
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
    setenv("FILE_DEM_GAT", "../testdata/dem_Gatta.txt", 1);
    setenv("FILE_DEM_SPC", "../testdata/dem_SanPi.txt", 1);
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
    wassert(actual(*cb->qual).statsEqual(1, 53.05, 99));
    wassert(actual(cb->top).statsEqual(0, 227845, 1, 12.72, 76));
    wassert(actual(cb->first_level).statsEqual(0, 0.33, 3));
    wassert(actual(cb->first_level_static).statsEqual(0, 0.33, 3));
    wassert(actual(cb->bb_first_level).statsEqual(0, 0.48, 2));
    wassert(actual(cb->beam_blocking).statsEqual(0, 17.10, 51));
    wassert(actual(cb->dem).statsEqual(-9999, 47.59, 4070.26));
    wassert(actual(cb->anaprop.quota).statsEqual(0, 0.00, 0));
    wassert(actual(cb->anaprop.dato_corrotto).statsEqual(0, 0.00, 0));
    wassert(actual(*cb->calcolo_vpr->flag_vpr).statsEqual(0, 0.73, 1));
    wassert(actual(cb->calcolo_vpr->corr_polar).statsEqual(0, 0.00, 0));
    wassert(actual(cb->calcolo_vpr->neve).statsEqual(0, 0.00, 0));
    // TODO: cb->stampa_vpr()

    Cart cart(cb->volume.max_beam_size());
    wassert(actual(cart.max_bin) == 849);
    cart.creo_cart(*cb);
//    print_stats("cart", cart, cerr);
    wassert(actual(cart.cart).statsEqual(0, 624947, 1, 39.18, 255));
    wassert(actual(cart.cartm).statsEqual(0.00, 2883204, 0.00, 0.00, 0.00));
    wassert(actual(cart.topxy).statsEqual(0, 2545764, 2, 15.70, 76));
    wassert(actual(cart.qual_Z_cart).statsEqual(0, 624947, 1, 19.30, 98));
    wassert(actual(cart.quota_cart).statsEqual(0, 2883204, 0, 0.00, 0));
    wassert(actual(cart.dato_corr_xy).statsEqual(0, 2883204, 0, 0.00, 0));
    wassert(actual(cart.beam_blocking_xy).statsEqual(0, 1688372, 1, 45.11, 51));
    wassert(actual(cart.elev_fin_xy).statsEqual(0, 2682479, 1, 1.15, 3));
    wassert(actual(cart.neve_cart).statsEqual(0, 624947, 1, 1.00, 1));
    wassert(actual(cart.corr_cart).statsEqual(0, 2883204, 0, 0.00, 0));
    wassert(actual(cart.conv_cart).statsEqual(0, 2883204, 0, 0.00, 0));

    CartLowris clow(cb->do_medium ? 512: 256, *cb, cart);
    wassert(actual(clow.ZLR_N_ELEMENTARY_PIXEL) == 4);
    wassert(actual(clow.ZLR_OFFSET) == 337);
    clow.creo_cart_z_lowris();
//     print_stats("clow", clow, cerr);
    wassert(actual(clow.z_out).statsEqual(0, 0, 1, 76.95, 255));
    wassert(actual(clow.qual_Z_1x1).statsEqual(0, 0, 1, 28.82, 97));
    wassert(actual(clow.quota_1x1).statsEqual(0, 0, 128, 128.00, 128));
    wassert(actual(clow.dato_corr_1x1).statsEqual(0, 65536, 0, 0.00, 0));
    wassert(actual(clow.elev_fin_1x1).statsEqual(0, 55707, 1, 1.16, 3));
    wassert(actual(clow.beam_blocking_1x1).statsEqual(0, 37495, 1, 45.41, 51));
    wassert(actual(clow.top_1x1).statsEqual(0, 43679, 2, 13.37, 57));
    wassert(actual(clow.neve_1x1).statsEqual(0, 0, 1, 1.00, 1));
    wassert(actual(clow.corr_1x1).statsEqual(0, 65536, 0, 0.00, 0));
    wassert(actual(clow.conv_1x1).statsEqual(0, 65536, 0, 0.00, 0));
});

add_method("combina_profili1", []() {
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
    setenv("FILE_DEM_GAT", "../testdata/dem_Gatta.txt", 1);
    setenv("FILE_DEM_SPC", "../testdata/dem_SanPi.txt", 1);

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
    wassert(actual(*cb->qual).statsEqual(1, 53.05, 99));
    wassert(actual(cb->top).statsEqual(0, 227845, 1, 12.72, 76));
    wassert(actual(cb->first_level).statsEqual(0, 0.33, 3));
    wassert(actual(cb->first_level_static).statsEqual(0, 0.33, 3));
    wassert(actual(cb->bb_first_level).statsEqual(0, 0.48, 2));
    wassert(actual(cb->beam_blocking).statsEqual(0, 17.10, 51));
    wassert(actual(cb->dem).statsEqual(-9999, 47.59, 4070.26));
    wassert(actual(cb->anaprop.quota).statsEqual(0, 0.00, 0));
    wassert(actual(cb->anaprop.dato_corrotto).statsEqual(0, 0.00, 0));
    wassert(actual(*cb->calcolo_vpr->flag_vpr).statsEqual(0, 0.73, 1));
    wassert(actual(cb->calcolo_vpr->corr_polar).statsEqual(0, 0.00, 0));
    wassert(actual(cb->calcolo_vpr->neve).statsEqual(0, 0.00, 0));

    Cart cart(cb->volume.max_beam_size());
    wassert(actual(cart.max_bin) == 849);
    cart.creo_cart(*cb);
//     print_stats("cart", cart, cerr);
    wassert(actual(cart.cart).statsEqual(0, 624947, 1, 39.18, 255));
    wassert(actual(cart.cartm).statsEqual(0.00, 2883204, 0.00, 0.00, 0.00));
    wassert(actual(cart.topxy).statsEqual(0, 2545764, 2, 15.70, 76));
    wassert(actual(cart.qual_Z_cart).statsEqual(0, 624947, 1, 19.30, 98));
    wassert(actual(cart.quota_cart).statsEqual(0, 2883204, 0, 0.00, 0));
    wassert(actual(cart.dato_corr_xy).statsEqual(0, 2883204, 0, 0.00, 0));
    wassert(actual(cart.beam_blocking_xy).statsEqual(0, 1688372, 1, 45.11, 51));
    wassert(actual(cart.elev_fin_xy).statsEqual(0, 2682479, 1, 1.15, 3));
    wassert(actual(cart.neve_cart).statsEqual(0, 624947, 1, 1.00, 1));
    wassert(actual(cart.corr_cart).statsEqual(0, 2883204, 0, 0.00, 0));
    wassert(actual(cart.conv_cart).statsEqual(0, 2881984, 100, 100.00, 100));

    CartLowris clow(cb->do_medium ? 512: 256, *cb, cart);
    wassert(actual(clow.ZLR_N_ELEMENTARY_PIXEL) == 4);
    wassert(actual(clow.ZLR_OFFSET) == 337);
    clow.creo_cart_z_lowris();
//     print_stats("clow", clow, cerr);
    wassert(actual(clow.z_out).statsEqual(0, 0, 1, 76.94, 255));
    wassert(actual(clow.qual_Z_1x1).statsEqual(0, 0, 1, 28.82, 97));
    wassert(actual(clow.quota_1x1).statsEqual(0, 0, 128, 128.00, 128));
    wassert(actual(clow.dato_corr_1x1).statsEqual(0, 65536, 0, 0.00, 0));
    wassert(actual(clow.elev_fin_1x1).statsEqual(0, 55707, 1, 1.16, 3));
    wassert(actual(clow.beam_blocking_1x1).statsEqual(0, 37495, 1, 45.41, 51));
    wassert(actual(clow.top_1x1).statsEqual(0, 43679, 2, 13.37, 57));
    wassert(actual(clow.neve_1x1).statsEqual(0, 0, 1, 1.00, 1));
    wassert(actual(clow.corr_1x1).statsEqual(0, 65536, 0, 0.00, 0));
    wassert(actual(clow.conv_1x1).statsEqual(0, 65477, 100, 100.00, 100));

    // TODO: scrivo_out_file_bin
});

add_method("test_9", []() {
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
    wassert(actual(stats.count_zeros[0]) == 172245);
    wassert(actual(stats.count_zeros[1]) == 210954);
    wassert(actual(stats.count_zeros[2]) == 163036);
    wassert(actual(stats.count_zeros[3]) == 109587);
    wassert(actual(stats.count_zeros[4]) == 118624);
    wassert(actual(stats.count_ones[0]) ==    0);
    wassert(actual(stats.count_ones[1]) ==    0);
    wassert(actual(stats.count_ones[2]) ==    0);
    wassert(actual(stats.count_ones[3]) ==    0);
    wassert(actual(stats.count_ones[4]) ==    0);
    wassert(actual(stats.count_others[0]) == 167355);
    wassert(actual(stats.count_others[1]) == 128646);
    wassert(actual(stats.count_others[2]) == 100964);
    wassert(actual(stats.count_others[3]) ==  88013);
    wassert(actual(stats.count_others[4]) ==  78976);
    wassert(actual(stats.sum_others[0]) == 23175278);
    wassert(actual(stats.sum_others[1]) == 16731796);
    wassert(actual(stats.sum_others[2]) == 12611046);
    wassert(actual(stats.sum_others[3]) == 10596468);
    wassert(actual(stats.sum_others[4]) ==  9176701);
    //print_stats("cb", *cb, cerr);

    Cart cart(cb->volume.max_beam_size());
    wassert(actual(cart.max_bin) == 849);
    cart.creo_cart(*cb);
//    print_stats("cart", cart, cerr);
    wassert(actual(cart.cart).statsEqual(0, 624947, 1, 43.33, 255));
    wassert(actual(cart.cartm).statsEqual(0.00, 2883204, 0.00, 0.00, 0.00));
    wassert(actual(cart.topxy).statsEqual(0, 2539486, 1, 15.63, 76));
    wassert(actual(cart.qual_Z_cart).statsEqual(0, 2883204, 0, 0.00, 0));
    wassert(actual(cart.quota_cart).statsEqual(0, 2883204, 0, 0.00, 0));
    wassert(actual(cart.dato_corr_xy).statsEqual(0, 2883204, 0, 0.00, 0));
    wassert(actual(cart.beam_blocking_xy).statsEqual(0, 2883204, 0, 0.00, 0));
    wassert(actual(cart.elev_fin_xy).statsEqual(0, 2883204, 0, 0.00, 0));
    wassert(actual(cart.neve_cart).statsEqual(0, 2883204, 0, 0.00, 0));
    wassert(actual(cart.corr_cart).statsEqual(0, 2883204, 0, 0.00, 0));
    wassert(actual(cart.conv_cart).statsEqual(0, 2883204, 0, 0.00, 0));

    CartLowris clow(cb->do_medium ? 512: 256, *cb, cart);
    wassert(actual(clow.ZLR_N_ELEMENTARY_PIXEL) == 4);
    wassert(actual(clow.ZLR_OFFSET) == -175);
    clow.creo_cart_z_lowris();
//     print_stats("clow", clow, cerr);
    wassert(actual(clow.z_out).statsEqual(0, 120384, 1, 53.52, 255));
    wassert(actual(clow.qual_Z_1x1).statsEqual(0, 262144, 0, 0.00, 0));
    wassert(actual(clow.quota_1x1).statsEqual(0, 262144, 0, 0.00, 0));
    wassert(actual(clow.dato_corr_1x1).statsEqual(0, 262144, 0, 0.00, 0));
    wassert(actual(clow.elev_fin_1x1).statsEqual(0, 262144, 0, 0.00, 0));
    wassert(actual(clow.beam_blocking_1x1).statsEqual(0, 262144, 0, 0.00, 0));
    wassert(actual(clow.top_1x1).statsEqual(0, 233808, 2, 17.16, 76));
    wassert(actual(clow.neve_1x1).statsEqual(0, 262144, 0, 0.00, 0));
    wassert(actual(clow.corr_1x1).statsEqual(0, 262144, 0, 0.00, 0));
    wassert(actual(clow.conv_1x1).statsEqual(0, 262144, 0, 0.00, 0));
});

add_method("test_10", []() {
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
    setenv("FILE_DEM_GAT", "../testdata/dem_Gatta.txt", 1);
    setenv("FILE_DEM_SPC", "../testdata/dem_SanPi.txt", 1);
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
    wassert(actual(cb->volume).statsEqual(MISSING_DB, 0, -31.5, -26.43, 77.85));

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
    wassert(actual(*cb->qual).statsEqual(1, 58.52, 99));
    wassert(actual(cb->top).statsEqual(0, 323585, 1, 8.03, 108));
    wassert(actual(cb->first_level).statsEqual(0, 0.33, 3));
    wassert(actual(cb->first_level_static).statsEqual(0, 0.33, 3));
    wassert(actual(cb->bb_first_level).statsEqual(0, 0.41, 2));
    wassert(actual(cb->beam_blocking).statsEqual(0, 15.83, 51));
    wassert(actual(cb->dem).statsEqual(-9999, 47.59, 4070.26));
    wassert(actual(cb->anaprop.quota).statsEqual(0, 0.00, 0));
    wassert(actual(cb->anaprop.dato_corrotto).statsEqual(0, 0.00, 0));
    wassert(actual(*cb->calcolo_vpr->flag_vpr).statsEqual(0, 0.94, 1));
    wassert(actual(cb->calcolo_vpr->corr_polar).statsEqual(0, 0.00, 0));
    wassert(actual(cb->calcolo_vpr->neve).statsEqual(0, 0.00, 0));

    Cart cart(cb->volume.max_beam_size());
    wassert(actual(cart.max_bin) == 849);
    cart.creo_cart(*cb);
//     print_stats("cart", cart, cerr);
    wassert(actual(cart.cart).statsEqual(0, 624947, 1, 8.56, 255));
    wassert(actual(cart.cartm).statsEqual(0.00, 2883204, 0.00, 0.00, 0.00));
    wassert(actual(cart.topxy).statsEqual(0, 2817169, 1, 14.25, 108));
    wassert(actual(cart.qual_Z_cart).statsEqual(0, 624947, 1, 18.97, 98));
    wassert(actual(cart.quota_cart).statsEqual(0, 2883204, 0, 0.00, 0));
    wassert(actual(cart.dato_corr_xy).statsEqual(0, 2883204, 0, 0.00, 0));
    wassert(actual(cart.beam_blocking_xy).statsEqual(0, 1769167, 1, 49.92, 51));
    wassert(actual(cart.elev_fin_xy).statsEqual(0, 2883204, 0, 0.00, 0));
    wassert(actual(cart.neve_cart).statsEqual(0, 624947, 1, 1.00, 1));
    wassert(actual(cart.corr_cart).statsEqual(0, 2883204, 0, 0.00, 0));
    wassert(actual(cart.conv_cart).statsEqual(0, 2883204, 0, 0.00, 0));

    CartLowris clow(cb->do_medium ? 512: 256, *cb, cart);
    wassert(actual(clow.ZLR_N_ELEMENTARY_PIXEL) == 4);
    wassert(actual(clow.ZLR_OFFSET) == 337);
    clow.creo_cart_z_lowris();
//     print_stats("clow", clow, cerr);
    wassert(actual(clow.z_out).statsEqual(0, 0, 1, 24.54, 255));
    wassert(actual(clow.qual_Z_1x1).statsEqual(0, 0, 1, 31.62, 92));
    wassert(actual(clow.quota_1x1).statsEqual(0, 0, 128, 128.00, 128));
    wassert(actual(clow.dato_corr_1x1).statsEqual(0, 65536, 0, 0.00, 0));
    wassert(actual(clow.elev_fin_1x1).statsEqual(0, 65536, 0, 0.00, 0));
    wassert(actual(clow.beam_blocking_1x1).statsEqual(0, 43435, 1, 50.17, 51));
    wassert(actual(clow.top_1x1).statsEqual(0, 59537, 1, 11.87, 108));
    wassert(actual(clow.neve_1x1).statsEqual(0, 0, 1, 1.00, 1));
    wassert(actual(clow.corr_1x1).statsEqual(0, 65536, 0, 0.00, 0));
    wassert(actual(clow.conv_1x1).statsEqual(0, 65536, 0, 0.00, 0));
});

add_method("test_11", []() {
    LOG_CATEGORY("Test");
    LOG_INFO ("Start test 6");

    // versione BB che corrisponde al parametro algo_corto
    static const char* fname = "../esplosione/2014-03-01-09-15-00.itgat.PVOL.0.h5";
    unsetwork();
    setenv("DIR_OUT_PP_BLOC", "../esplosione", 1);
    setenv("VPR0_FILE"      , "../esplosione/vpr_GAT", 1);
    setenv("LAST_VPR"       , "../esplosione/last_vpr_GAT", 1);
    setenv("VPR_HMAX"       , "../esplosione/vpr_hmax_GAT",1); 
    setenv("VPR_HEATING"    , "../esplosione/vpr_heat_GAT",1);
    setenv("LOG_VPR"        , "../esplosione/log_VPR_${SITO}",1);
    setenv("TEST_VPR"       , "../esplosione/test_vpr",1);
    setenv("FILE_T"     , "../esplosione/temperature.txt",1);
    setenv("FIRST_LEVEL_FILE", "../dati/FIRST_LEVEL_corto_GAT_PRI-EST_2011", 1);
    setenv("FILE_ZERO_TERMICO"  , "../esplosione/0termico.prev", 1);
    setenv("VPR_ARCH"           , "../esplosione/201403010915_vpr_GAT",1);
    setenv("FILE_DEM_GAT", "../testdata/dem_Gatta.txt", 1);
    setenv("FILE_DEM_SPC", "../testdata/dem_SanPi.txt", 1);
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
});

add_method("test_12", []() {
    LOG_CATEGORY("Test");
    LOG_INFO ("Start test 12");

    // versione BB che corrisponde al parametro algo_corto
    static const char* fname = "../testdata/vpr/2014-03-01-01-35-00.itgat.PVOL.0.h5";
    unsetwork();
    setenv("DIR_OUT_PP_BLOC", "../testdata/vpr", 1);
    setenv("VPR0_FILE"      , "../testdata/vpr/vpr_GAT", 1);
    setenv("LAST_VPR"       , "../testdata/vpr/last_vpr_GAT", 1);
    setenv("VPR_HMAX"       , "../testdata/vpr/vpr_hmax_GAT",1); 
    setenv("VPR_HEATING"    , "../testdata/vpr/vpr_heat_GAT",1);
    setenv("LOG_VPR"        , "../testdata/vpr/log_VPR_${SITO}",1);
    setenv("TEST_VPR"       , "../testdata/vpr/test_vpr",1);
    setenv("FILE_T"     , "../testdata/vpr/temperature.txt",1);
//setenv("FIRST_LEVEL_FILE", "../dati/FIRST_LEVEL_corto_GAT_INV_2011_15el", 1);
    setenv("FIRST_LEVEL_FILE", "../dati/FIRST_LEVEL_corto+medio_GAT_INV_2011", 1);
    setenv("FILE_ZERO_TERMICO"  , "../testdata/vpr/0termico.prev", 1);
    setenv("OUTPUT_Z_LOWRIS_DIR","../testdata/vpr",1);
    setenv("DIR_QUALITY","../testdata/vpr",1);
    setenv("DIR_DEBUG","../testdata/vpr",1);
    setenv("FILE_DEM_GAT", "../testdata/dem_Gatta.txt", 1);
    setenv("FILE_DEM_SPC", "../testdata/dem_SanPi.txt", 1);

    printwork();

    CBTest test("GAT", false, 1024);
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
});

}

}
