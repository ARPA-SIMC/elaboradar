#include "radarelab/utils/tests.h"
#include <radarelab/logging.h>
#include "cum_bac.h"
#include "config.h"
#include <cstdio>
#include <cstdlib>
#include <vector>
#include "site.h"
#include <radarelab/cart.h>
#include "test-utils.h"
#include <unistd.h>

using namespace radarelab::utils::tests;
using namespace radarelab;
using namespace elaboradar;
using namespace testradar;

using namespace std;

namespace {

static const auto mis = IndexMapping::missing;

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
};

class Tests : public TestCase
{
    using TestCase::TestCase;

    void register_tests() override;
} test("cart");

void Tests::register_tests() {

/// Check the basic mapping calculations
add_method("basic", []() {
    CoordinateMapping mapping(494);

    // Write out images in /tmp
    //setenv("DIR_DEBUG", "/tmp", 1);
    //Config cfg;
    //Assets assets(cfg);
    //assets.configure(Site::get("GAT"), time(NULL));
    //assets.write_gdal_image(mapping.map_azimuth, "DIR_DEBUG", "map_azimuth", "png");
    //assets.write_gdal_image(mapping.map_range, "DIR_DEBUG", "map_range", "png");

    wassert(actual(mapping.beam_size) == 494u);

    // Check range

    // Range in the middle (on top of the radar) should be 0
    wassert(actual(floor(mapping.map_range(493, 493))) == 0);
    wassert(actual(floor(mapping.map_range(493, 494))) == 0);
    wassert(actual(floor(mapping.map_range(494, 493))) == 0);
    wassert(actual(floor(mapping.map_range(494, 494))) == 0);

    // Range at map edges at 90° angles should be beam_size
    wassert(actual(floor(mapping.map_range(  0, 493))) == 493); // Middle left
    wassert(actual(floor(mapping.map_range(  0, 494))) == 493); // Middle left
    wassert(actual(floor(mapping.map_range(987, 493))) == 493); // Middle right
    wassert(actual(floor(mapping.map_range(987, 494))) == 493); // Middle right
    wassert(actual(floor(mapping.map_range(493,   0))) == 493); // Top middle
    wassert(actual(floor(mapping.map_range(494,   0))) == 493); // Top middle
    wassert(actual(floor(mapping.map_range(493, 987))) == 493); // Bottom middle
    wassert(actual(floor(mapping.map_range(494, 987))) == 493); // Bottom middle

    // Range at map corners should be missing (out of range)
    wassert(actual(floor(mapping.map_range(  0,   0))) == 697);
    wassert(actual(floor(mapping.map_range(  0, 987))) == 697);
    wassert(actual(floor(mapping.map_range(987, 987))) == 697);
    wassert(actual(floor(mapping.map_range(987,   0))) == 697);


    // Check azimuth

    // Azimuth at the map corners and edge middle points
    wassert(actual(floor(mapping.map_azimuth(  0,   0) * 400/360)) == 350); // Top left
    wassert(actual(floor(mapping.map_azimuth(  0, 493) * 400/360)) == 399); // Top middle
    wassert(actual(floor(mapping.map_azimuth(  0, 494) * 400/360)) ==   0); // Top middle
    wassert(actual(floor(mapping.map_azimuth(  0, 987) * 400/360)) ==  50); // Top right
    wassert(actual(floor(mapping.map_azimuth(493, 987) * 400/360)) ==  99); // Middle right
    wassert(actual(floor(mapping.map_azimuth(494, 987) * 400/360)) == 100); // Middle right
    wassert(actual(floor(mapping.map_azimuth(987, 987) * 400/360)) == 150); // Bottom right
    wassert(actual(floor(mapping.map_azimuth(987, 494) * 400/360)) == 199); // Bottom middle
    wassert(actual(floor(mapping.map_azimuth(987, 493) * 400/360)) == 200); // Bottom middle
    wassert(actual(floor(mapping.map_azimuth(987,   0) * 400/360)) == 250); // Bottom left
    wassert(actual(floor(mapping.map_azimuth(494,   0) * 400/360)) == 299); // Middle left
    wassert(actual(floor(mapping.map_azimuth(493,   0) * 400/360)) == 300); // Middle left

    // Azimuth at the 4 corners around the centre should be at 45° angles
    wassert(actual(floor(mapping.map_azimuth(493, 493) * 400/360)) == 350);
    wassert(actual(floor(mapping.map_azimuth(493, 494) * 400/360)) ==  50);
    wassert(actual(floor(mapping.map_azimuth(494, 494) * 400/360)) == 150);
    wassert(actual(floor(mapping.map_azimuth(494, 493) * 400/360)) == 250);
});

/// Check the actual mapping calculations on an image
add_method("image", []() {
    static const char* fname = "../testdata/DBP2_070120141530_GATTATICO";

    CBTest test("GAT", false);
    test.read_sp20(fname, true);

    // Create the mapping
    FullsizeIndexMapping mapping(test.volume[0].beam_size);
    mapping.map_max_sample(test.volume[0]);


    // Write out images in /tmp
    //setenv("DIR_DEBUG", "/tmp", 1);
    //Assets assets(test.cfg);
    //assets.configure(test.site, test.volume.load_info->acq_date);
    //assets.write_gdal_image(mapping.map_azimuth, "DIR_DEBUG", "map_azimuth", "png");
    //assets.write_gdal_image(mapping.map_range, "DIR_DEBUG", "map_range", "png");

    wassert(actual(test.volume.beam_count) == 400u);
    wassert(actual(mapping.width) == 988);
    wassert(actual(mapping.height) == 988);

    const auto mis = IndexMapping::missing;

    // Check range

//for (unsigned y=0; y< test.volume[0].beam_size*2; y++) printf (" %d - %d %d \n",y,  mapping.map_range(y, 493), mapping.map_range(y, 494));
    // Range in the middle (on top of the radar) should be 0
    wassert(actual(mapping.map_range(493, 493)) == 0);
    wassert(actual(mapping.map_range(493, 494)) == 0);
    wassert(actual(mapping.map_range(494, 493)) == 0);
    wassert(actual(mapping.map_range(494, 494)) == 0);

    // Range at map edges at 90° angles should be beam_size
    wassert(actual(mapping.map_range(  0, 493)) == 493); // Middle left
    wassert(actual(mapping.map_range(  0, 494)) == 493); // Middle left
    wassert(actual(mapping.map_range(987, 493)) == 493); // Middle right
    wassert(actual(mapping.map_range(987, 494)) == 493); // Middle right
    wassert(actual(mapping.map_range(493,   0)) == 493); // Top middle
    wassert(actual(mapping.map_range(494,   0)) == 493); // Top middle
    wassert(actual(mapping.map_range(493, 987)) == 493); // Bottom middle
    wassert(actual(mapping.map_range(494, 987)) == 493); // Bottom middle

    // Range at map corners should be missing (out of range)
    wassert(actual(mapping.map_range(  0,   0)) == mis);
    wassert(actual(mapping.map_range(  0, 987)) == mis);
    wassert(actual(mapping.map_range(987, 987)) == mis);
    wassert(actual(mapping.map_range(987,   0)) == mis);

    // Check azimuth

    // Azimuth at the map corners and edge middle points
//for (unsigned x=0; x< test.volume[0].beam_size*2; x++) printf (" map_az %d - %d \n",x,  mapping.map_azimuth(493,x));
    wassert(actual(mapping.map_azimuth(  0,   0)) == mis); // Top left
    wassert(actual(mapping.map_azimuth(  0, 493)) ==   0); // Top middle
    wassert(actual(mapping.map_azimuth(  0, 494)) ==   0); // Top middle
    wassert(actual(mapping.map_azimuth(  0, 987)) == mis); // Top right
    wassert(actual(mapping.map_azimuth(493, 987)) == 100); // Middle right
    wassert(actual(mapping.map_azimuth(494, 987)) == 100); // Middle right
    wassert(actual(mapping.map_azimuth(987, 987)) == mis); // Bottom right
    wassert(actual(mapping.map_azimuth(987, 494)) == 200); // Bottom middle
    wassert(actual(mapping.map_azimuth(987, 493)) == 200); // Bottom middle
    wassert(actual(mapping.map_azimuth(987,   0)) == mis); // Bottom left
    wassert(actual(mapping.map_azimuth(494,   0)) == 300); // Middle left
    wassert(actual(mapping.map_azimuth(493,   0)) == 300); // Middle left

    // Azimuth at the 4 corners around the centre should be at 45° angles
 //   wassert(actual(mapping.map_azimuth(493, 493)) == 349);
 //   wassert(actual(mapping.map_azimuth(493, 494)) ==  49);
 //   wassert(actual(mapping.map_azimuth(494, 494)) == 149);
 //   wassert(actual(mapping.map_azimuth(494, 493)) == 249);
    // Azimuth at the 4 corners around the centre should be at 45° angles, but it takes the az for the max value
    wassert(actual(mapping.map_azimuth(493, 493)) == 286);
    wassert(actual(mapping.map_azimuth(493, 494)) == 386);
    wassert(actual(mapping.map_azimuth(494, 494)) ==  86);
    wassert(actual(mapping.map_azimuth(494, 493)) == 186);


    // 1:4 scaled mapping, perfectly fitting
    ScaledIndexMapping smapping(test.volume[0].beam_size, 247, 4);
    smapping.map_max_sample(test.volume[0], mapping);
    wassert(actual(smapping.width) == 247);
    wassert(actual(smapping.height) == 247);
    wassert(actual(smapping.image_offset) == 0);

    // Write out images in /tmp
    //setenv("DIR_DEBUG", "/tmp", 1);
    //Assets assets(test.cfg);
    //assets.configure(test.site, test.volume.load_info->acq_date);
    //assets.write_gdal_image(smapping.map_azimuth, "DIR_DEBUG", "map_azimuth", "png");
    //assets.write_gdal_image(smapping.map_range, "DIR_DEBUG", "map_range", "png");

    // Check range

    // Range in the middle (on top of the radar) should be 0
    wassert(actual(smapping.map_range(123, 123)) == 1);

    // Range at map edges at 90° angles should be beam_size
    wassert(actual(smapping.map_range(  0, 123)) == 492); // Middle left
    wassert(actual(smapping.map_range(246, 123)) == 490); // Middle right
    wassert(actual(smapping.map_range(123,   0)) == 493); // Top middle
    wassert(actual(smapping.map_range(123, 246)) == 490); // Bottom middle

    // Range at map corners should be missing (out of range)
    wassert(actual(smapping.map_range(  0,   0)) == mis);
    wassert(actual(smapping.map_range(  0, 246)) == mis);
    wassert(actual(smapping.map_range(246, 246)) == mis);
    wassert(actual(smapping.map_range(246,   0)) == mis);

    // Check azimuth

    // Azimuth at the map corners and edge middle points
    wassert(actual(smapping.map_azimuth(  0,   0)) == mis); // Top left
    wassert(actual(smapping.map_azimuth(  0, 123)) ==   0); // Top middle
    wassert(actual(smapping.map_azimuth(  0, 246)) == mis); // Top right
    wassert(actual(smapping.map_azimuth(123, 246)) == 100); // Middle right
    wassert(actual(smapping.map_azimuth(246, 246)) == mis); // Bottom right
    wassert(actual(smapping.map_azimuth(246, 123)) == 200); // Bottom middle
    wassert(actual(smapping.map_azimuth(246,   0)) == mis); // Bottom left
    wassert(actual(smapping.map_azimuth(123,   0)) == 300); // Middle left

    // Azimuth at the 4 corners around the centre should be at roughly 45°
    // angles
    //wassert(actual(smapping.map_azimuth(122, 122)) == 364);
    //wassert(actual(smapping.map_azimuth(122, 124)) ==  72);
    //wassert(actual(smapping.map_azimuth(124, 124)) == 149);
    //wassert(actual(smapping.map_azimuth(124, 122)) == 261);
    wassert(actual(smapping.map_azimuth(122, 122)) == 346);
    wassert(actual(smapping.map_azimuth(122, 124)) ==  72);
    wassert(actual(smapping.map_azimuth(124, 124)) == 174);
    wassert(actual(smapping.map_azimuth(124, 122)) == 223);
});

}

}
