#include <test-utils.h>
#include <cum_bac.h>

namespace testradar {

void print_stats(const std::string& name, const cumbac::CUM_BAC& cb, std::ostream& out)
{
    print_stats("*" + name + "->qual", *cb.qual, out);
    print_stats("*" + name + "->calcolo_vpr->flag_vpr", *cb.calcolo_vpr->flag_vpr, out);
    print_stats(name + "->top", cb.top, out);
}

void print_stats(const std::string& name, const cumbac::Cart& cart, std::ostream& out)
{
    print_stats(name + ".cart", cart.cart, out);
    print_stats(name + ".cartm", cart.cartm, out);
    print_stats(name + ".topxy", cart.topxy, out);
    print_stats(name + ".qual_Z_cart", cart.qual_Z_cart, out);
    print_stats(name + ".quota_cart", cart.quota_cart, out);
    print_stats(name + ".dato_corr_xy", cart.dato_corr_xy, out);
    print_stats(name + ".beam_blocking_xy", cart.beam_blocking_xy, out);
    print_stats(name + ".elev_fin_xy", cart.elev_fin_xy, out);
    print_stats(name + ".neve_cart", cart.neve_cart, out);
    print_stats(name + ".corr_cart", cart.corr_cart, out);
    print_stats(name + ".conv_cart", cart.conv_cart, out);
}

void print_stats(const std::string& name, const cumbac::CartLowris& cart, std::ostream& out)
{
    print_stats(name + ".z_out", cart.z_out, out);
    print_stats(name + ".qual_Z_1x1", cart.qual_Z_1x1, out);
    print_stats(name + ".quota_1x1", cart.quota_1x1, out);
    print_stats(name + ".dato_corr_1x1", cart.dato_corr_1x1, out);
    print_stats(name + ".elev_fin_1x1", cart.elev_fin_1x1, out);
    print_stats(name + ".beam_blocking_1x1", cart.beam_blocking_1x1, out);
    print_stats(name + ".top_1x1", cart.top_1x1, out);
    print_stats(name + ".neve_1x1", cart.neve_1x1, out);
    print_stats(name + ".corr_1x1", cart.corr_1x1, out);
    print_stats(name + ".conv_1x1", cart.conv_1x1, out);
}

}
