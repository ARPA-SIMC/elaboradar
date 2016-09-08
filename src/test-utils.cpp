#include <test-utils.h>
#include <cum_bac.h>
#include <cartproducts.h>

using namespace elaboradar;

namespace testradar {

void print_stats(const std::string& name, const elaboradar::CUM_BAC& cb, std::ostream& out)
{
    print_stats("*" + name + "->qual", *cb.qual, out);
    print_stats(name + "->top", cb.top, 0, out);
    print_stats(name + "->first_level", cb.first_level, out);
    print_stats(name + "->first_level_static", cb.first_level_static, out);
    print_stats(name + "->bb_first_level", cb.bb_first_level, out);
    print_stats(name + "->beam_blocking", cb.beam_blocking, out);
    print_stats(name + "->dem", cb.dem, out);
    print_stats(name + "->anaprop.quota", cb.anaprop.quota, out);
    print_stats(name + "->anaprop.dato_corrotto", cb.anaprop.dato_corrotto, out);
    print_stats(name + "->flag_vpr", cb.flag_vpr, out);
    print_stats(name + "->calcolo_vpr->corr_polar", cb.calcolo_vpr->corr_polar, out);
    print_stats(name + "->calcolo_vpr->neve", cb.calcolo_vpr->neve, out);
}

#if 0
void print_stats(const std::string& name, const elaboradar::Cart& cart, std::ostream& out)
{
    print_stats(name + ".cart", cart.cart, 0, out);
    print_stats(name + ".cartm", cart.cartm, 0.0, out);
    print_stats(name + ".topxy", cart.topxy, 0, out);
    print_stats(name + ".qual_Z_cart", cart.qual_Z_cart, 0, out);
    print_stats(name + ".quota_cart", cart.quota_cart, 0, out);
    print_stats(name + ".dato_corr_xy", cart.dato_corr_xy, 0, out);
    print_stats(name + ".beam_blocking_xy", cart.beam_blocking_xy, 0, out);
    print_stats(name + ".elev_fin_xy", cart.elev_fin_xy, 0, out);
    print_stats(name + ".neve_cart", cart.neve_cart, 0, out);
    print_stats(name + ".corr_cart", cart.corr_cart, 0, out);
    print_stats(name + ".conv_cart", cart.conv_cart, 0, out);
}

void print_stats(const std::string& name, const elaboradar::CartLowris& cart, std::ostream& out)
{
    print_stats(name + ".z_out", cart.z_out, 0, out);
    print_stats(name + ".qual_Z_1x1", cart.qual_Z_1x1, 0, out);
    print_stats(name + ".quota_1x1", cart.quota_1x1, 0, out);
    print_stats(name + ".dato_corr_1x1", cart.dato_corr_1x1, 0, out);
    print_stats(name + ".elev_fin_1x1", cart.elev_fin_1x1, 0, out);
    print_stats(name + ".beam_blocking_1x1", cart.beam_blocking_1x1, 0, out);
    print_stats(name + ".top_1x1", cart.top_1x1, 0, out);
    print_stats(name + ".neve_1x1", cart.neve_1x1, 0, out);
    print_stats(name + ".corr_1x1", cart.corr_1x1, 0, out);
    print_stats(name + ".conv_1x1", cart.conv_1x1, 0, out);
}

void print_stats(const std::string& name, const elaboradar::CartProducts& cart, std::ostream& out)
{
    print_stats(name + ".z_out", cart.z_out, 0, out);
    print_stats(name + ".qual_Z_1x1", cart.qual_Z_1x1, 0, out);
    print_stats(name + ".quota_1x1", cart.quota_1x1, 0, out);
    print_stats(name + ".dato_corr_1x1", cart.dato_corr_1x1, 0, out);
    print_stats(name + ".elev_fin_1x1", cart.elev_fin_1x1, 0, out);
    print_stats(name + ".beam_blocking_1x1", cart.beam_blocking_1x1, 0, out);
    print_stats(name + ".top_1x1", cart.top_1x1, 0, out);
    print_stats(name + ".neve_1x1", cart.neve_1x1, 0, out);
    print_stats(name + ".corr_1x1", cart.corr_1x1, 0, out);
    print_stats(name + ".conv_1x1", cart.conv_1x1, 0, out);
}
#endif

}
