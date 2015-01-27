#include "image.h"
#include <memory>
#include <gdal.h>
#include <gdal_priv.h>

using namespace std;

namespace elaboradar {

void gdal_init_once()
{
    static bool initialized = false;
    if (initialized) return;
    GDALAllRegister();
    initialized = true;
}

template<typename T> GDALDataType get_gdal_datatype() { throw std::runtime_error("get_gdal_datatype called for unsupported type"); }
template<> GDALDataType get_gdal_datatype<unsigned char>() { return GDT_Byte; }
template<> GDALDataType get_gdal_datatype<unsigned short>() { return GDT_UInt16; }
template<> GDALDataType get_gdal_datatype<short>() { return GDT_Int16; }
template<> GDALDataType get_gdal_datatype<double>() { return GDT_Float64; }
template<> GDALDataType get_gdal_datatype<int>() { return GDT_Int32; }


template<typename T>
class MatrixDataset : public GDALDataset
{
public:
    const Matrix2D<T>& image;

    MatrixDataset(const Matrix2D<T>& image);

    //int GetRasterXSize() override { return image.cols(); }
    //int GetRasterYSize() override { return image.rows(); }
    //int GetRasterCount() override { return 1; }

    //const char* GetProjectionRef() override { }

    // CPLErr GetGeoTransform (double *) override
};

template<typename T>
class MatrixRasterBand : public GDALRasterBand
{
public:
    const Matrix2D<T>& image;

    MatrixRasterBand(MatrixDataset<T>& ds)
        : image(ds.image)
    {
        poDS = &ds;
        nBand = 1;
        nBlockXSize = image.cols();
        nBlockYSize = image.rows();
        // SetDescription(name.c_str());

        eDataType = get_gdal_datatype<T>();
    }

    // const char* GetUnitType() override {}

    CPLErr IReadBlock(int xblock, int yblock, void *buf) override
    {
        if (xblock != 0 || yblock != 0)
        {
            CPLError(CE_Failure, CPLE_AppDefined, "Invalid block number");
            return CE_Failure;
        }

        memcpy(buf, image.data(), image.size() * sizeof(T));

        return CE_None;
    }

    // double GetOffset(int* pbSuccess=NULL) override {}
    // double GetScale(int* pbSuccess=NULL) override {}
    // double GetNoDataValue(int* pbSuccess=NULL) override {};
};

template<typename T>
MatrixDataset<T>::MatrixDataset(const Matrix2D<T>& image)
    : image(image)
{
    nRasterXSize = image.cols();
    nRasterYSize = image.rows();
    SetBand(1, new MatrixRasterBand<T>(*this));
}

template<typename T>
void write_image(const Matrix2D<T>& image, const std::string& fname, const std::string& format)
{
    unique_ptr<MatrixDataset<T>> src(new MatrixDataset<T>(image));
    GDALDriver *driver = GetGDALDriverManager()->GetDriverByName(format.c_str());
    if (driver == NULL)
        throw std::runtime_error("driver not found for " + format);

    GDALDataset* dst = driver->CreateCopy(fname.c_str(), src.get(), false, NULL, NULL, NULL);
    if (dst == NULL)
        throw std::runtime_error("cannot create " + fname);
    GDALClose(dst);
}


template void write_image(const Matrix2D<unsigned char>&, const std::string&, const std::string&);
template void write_image(const Matrix2D<unsigned short>&, const std::string&, const std::string&);
template void write_image(const Matrix2D<double>&, const std::string&, const std::string&);
template void write_image(const Matrix2D<int>&, const std::string&, const std::string&);
template void write_image(const Matrix2D<short>&, const std::string&, const std::string&);

std::string gdal_extension_for_format(const std::string& format)
{
    GDALDriver *driver = GetGDALDriverManager()->GetDriverByName(format.c_str());
    if (driver == NULL)
        throw std::runtime_error("driver not found for " + format);

    const char* ext = driver->GetMetadataItem(GDAL_DMD_EXTENSION, NULL);
    if (ext == NULL)
        throw std::runtime_error("extension not found for format " + format);
    return ext;
}

}
