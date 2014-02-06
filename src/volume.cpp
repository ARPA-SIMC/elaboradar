#include "volume.h"
#include <cstring>

/// This needs to be a global variable, as it is expected by libsp20
int elev_array[NEL];

namespace cumbac {

Volume::Volume()
    : acq_date(0), size_cell(0), declutter_rsp(false)
{
    memset(vol_pol, 0, sizeof(vol_pol));
    memset(nbeam_elev, 0, sizeof(nbeam_elev));
}

void Volume::fill_beam(double theta, double alpha, unsigned size, const unsigned char* data)
{
    int teta = theta / FATT_MOLT_EL;

    int el_num = elevation_index_MDB(teta);
    if (el_num > NEL) return;

    int alfa = alpha / FATT_MOLT_AZ;
    if (alfa >= 4096) return;

    int az_num = azimut_index_MDB(alfa);

    if (az_num == 0)
    {
        printf("fbeam ϑ%f→%d α%f→%d %u", theta, el_num, alpha, az_num, size);
        for (unsigned i = 0; i < 20; ++i)
            printf(" %d", (int)data[i]);
        printf("\n");
    }

    merge_beam(&vol_pol[el_num][az_num], theta, alpha, az_num, el_num, size, data);
    if(az_num*0.9 - alpha < 0.)
    {
        int new_az_num = (az_num + 1) % 400;
        merge_beam(&vol_pol[el_num][new_az_num], theta, alpha, new_az_num, el_num, size, data);
    }
    else if(az_num*0.9 - alpha > 0.)
    {
        int new_az_num = (az_num -1+400) %400;
        merge_beam(&vol_pol[el_num][new_az_num], theta, alpha, new_az_num, el_num, size, data);
    }
}

void Volume::merge_beam(VOL_POL* raggio, double theta, double alpha, int az_num, int el_num, unsigned size, const unsigned char* dati)
{
    if (raggio->flag == 0)
    {
        for (unsigned i = 0; i < size; i++)
        {
            if(dati[i])
                raggio->ray[i] = dati[i];
            else
                raggio->ray[i] = 1;
        }
        nbeam_elev[el_num]++;
    }
    else
        for (unsigned i = 0; i < size; i++)
            if(raggio->ray[i]<dati[i])
                raggio->ray[i]=dati[i];

    raggio->flag=1;
    raggio->b_header.alfa =(short)(az_num*.9/FATT_MOLT_AZ);
    raggio->b_header.teta = elev_array[el_num];
    raggio->alfa_true = alpha / FATT_MOLT_AZ;
    raggio->teta_true = theta / FATT_MOLT_EL;
    raggio->b_header.tipo_gran = INDEX_Z;  // FIXME: to be changed when we load different quantities
    raggio->b_header.max_bin = size;
}



}
