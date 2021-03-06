#include <radarelab/cylindrical.h>
#include <cmath>

#define DTOR  M_PI/180. /* esternalizzo?*/ //fattore conversione gradi-radianti

using namespace std;

namespace radarelab {

CylindricalVolume::CylindricalVolume(const Volume<double>& volume, double missing_value, double x_res, double z_res)
{
    const double size_cell = volume[0].cell_size;
    double range_min=0.5 * size_cell/1000.;
    double range_maxLowestRay = (volume[0].beam_size - 0.5) * size_cell / 1000.;

    double xmin = floor(range_min*cos(volume.elevation_max()*DTOR)); // distanza orizzontale minima dal radar
    double zmin = volume[0].sample_height(0) / 1000. + volume.radarSite.getTotalHeight(); // quota  minima in prop standard
    double xmax = floor(range_maxLowestRay*cos(volume.elevation_min()*DTOR)); // distanza orizzontale massima dal radar
    double zmax = volume.back().sample_height(volume.back().beam_size - 1) / 1000. + volume.radarSite.getTotalHeight();//quota massima
    //LOG_DEBUG(" Range min maxL maxU  %7.3f %7.3f %7.3f  --  xmin %7.3f xmax %7.3f zmin %7.3f zmax %7.3f", range_min, range_maxLowestRay, range_maxUpperRay, xmin,xmax,zmin,zmax);

    x_size = (xmax - xmin) / x_res; //dimensione orizzontale
    // FIXME: usiamo volume.max_beam_size invece di MyMAX_BIN?
    if (x_size > volume.max_beam_size()) x_size=volume.max_beam_size();
    z_size = (zmax - zmin) / z_res; //dimensione verticale

    resol[0] = x_res;
    resol[1] = z_res;
    slices.reserve(volume.beam_count);
    for (unsigned i = 0; i < volume.beam_count; ++i)
        slices.push_back(new Matrix2D<double>(Matrix2D<double>::Constant(x_size, z_size, missing_value)));

    resample(volume);
}

void CylindricalVolume::resample(const Volume<double>& volume)
{
    unsigned max_bin = x_size;
    /* ---------------------------------- */
    /*           FASE 1 */
    /* ---------------------------------- */
    /*    Costruzione matrice generale per indicizzare il caricamento dei dati */
    /*    da coordinate radar a coordinate (X,Z) */
    /*  Calcolo distanza per ogni punto sul raggio */
    /*  xx contiene la distanza sulla superfice terrestre in funzione del range radar e dell'elevazione */
    /*  Metodo 1 - -Calcolo le coordinate di ogni punto del RHI mediante utilizzo di cicli */

    // estremi x e z (si procede per rhi)
#warning to compute scan by scan?
    double cell_size = volume.scan(0).cell_size;
    double range_min=0.5 * cell_size / 1000.;
    //double range_max=(max_bin-0.5) * size_cell/1000.;

    double xmin=floor(range_min*cos(volume.elevation_max()*DTOR)); // distanza orizzontale minima dal radar
    double zmin=volume[0].sample_height(0) / 1000. + volume.radarSite.getTotalHeight(); // quota  minima in prop standard


    //float w_size[2]={3.,1.5}; //dimensione della matrice pesi
    const double w_size[2]={3.,0.3}; //dimensione della matrice pesi

    //LOG_INFO("calcolati range_min e range_max , dimensione orizzontale e dimensione verticale range_min=%f  range_max=%f x_size=%d z_size=%d",range_min,range_max,x_size,z_size);

    int w_x_size=ceil((w_size[0]/resol[0])/2)*2+1; //dimensione x matrice pesi
    int w_z_size=ceil((w_size[1]/resol[1])/2)*2+1; //dimensione z matrice pesi

    if (w_x_size < 3) w_x_size=3;
    if (w_z_size < 3) w_z_size=3;

    int w_x_size_2=w_x_size/2;
    int w_z_size_2=w_z_size/2;

    Matrix2D<int> i_xx_min(Matrix2D<int>::Zero(max_bin, volume.size()));
    Matrix2D<int> i_zz_min(Matrix2D<int>::Zero(max_bin, volume.size()));
    Matrix2D<int> im(Matrix2D<int>::Zero(max_bin, volume.size()));
    Matrix2D<int> ix(Matrix2D<int>::Zero(max_bin, volume.size()));
    Matrix2D<int> jm(Matrix2D<int>::Zero(max_bin, volume.size()));
    Matrix2D<int> jx(Matrix2D<int>::Zero(max_bin, volume.size()));
    for (unsigned i = 0; i < max_bin; i++){
        double range = (i + 0.5) * cell_size/1000.;

          for (unsigned k=0; k < volume.size(); k++){
            double elev_rad = volume.scan(k).elevation * DTOR;
            double zz = volume.scan(k).sample_height(i) / 1000. + volume.radarSite.getTotalHeight();// quota
            double xx = range*cos(elev_rad); // distanza
            int i_zz=floor((zz - zmin)/resol[1]);// indice in z, nella proiezione cilindrica, del punto i,k  PPA
            int i_xx=floor((xx - xmin)/resol[0]);// indice in x, nella proiezione cilindrica, del punto i,k  PPA
            // Enrico RHI_ind[k][i]=i_xx+i_zz*x_size;
            //shift orizzontale negativo del punto di indice i_xx per costruire la finestra in x
            // se l'estremo minimo in x della finestra è negativo assegno come shift il massimo possibile e cioè la distanza del punto dall'origine
            i_xx_min(i, k)=i_xx;
            if (i_xx -  w_x_size_2 >= 0)
                i_xx_min(i, k)= w_x_size_2;

            //shift orizzontale positivo attorno al punto di indice i_xx per costruire la finestra in x
            int i_xx_max = x_size-i_xx-1;
            if (i_xx+w_x_size_2 < (int) x_size)
                i_xx_max = w_x_size_2;

            //shift verticale negativo attorno al punto di indice i_zz per costruire la finestra in z
            i_zz_min(i, k)=i_zz;
            if (i_zz_min(i, k) -  w_z_size_2 > 0)
                i_zz_min(i, k) = w_z_size_2;

            //shift verticale positivo attorno al punto di indice i_zz per costruire la finestra in z
            int i_zz_max = z_size-i_zz-1;
            if (i_zz+w_z_size_2 < z_size)
                i_zz_max = w_z_size_2;

            //indici minimo e massimo in x e z per definire la finestra sul punto
            im(i, k)=i_xx-i_xx_min(i, k);
            ix(i, k)=i_xx+i_xx_max;
            jm(i, k)=i_zz-i_zz_min(i, k);
            jx(i, k)=i_zz+i_zz_max;

        }
    }
    /*
       ;------------------------------------------------------------------------------
       ;          FASE 2
       ;------------------------------------------------------------------------------
       ;   Costruzione matrice pesi
       ;   Questa matrice contiene i pesi (in funzione della distanza) per ogni punto.
       ;-----------------------------------------------------------------------------*/

    vector<double> w_x(w_x_size);
    for (unsigned k=0;k<(unsigned) w_x_size;k++)
        w_x[k]=exp(-pow(k-w_x_size_2,2.)/pow(w_x_size_2/2.,2.));

    vector<double> w_z(w_z_size);
    for (unsigned k=0;k<(unsigned)w_z_size;k++)
        w_z[k]=exp(-pow(k-w_z_size_2,2.)/pow(w_z_size_2/2.,2.));

    Matrix2D<double> w_tot(w_x_size, w_z_size);
    for (unsigned i=0;i<(unsigned)w_x_size;i++){
        for (unsigned j=0;j<(unsigned)w_z_size;j++){
            w_tot(i, j)=w_x[i]*w_z[j];
        }
    }

    /* ;----------------------------------------------------------- */
    /* ;   Matrici per puntare sul piano cartesiano velocemente */
    /* ;---------------------------------- */
    /* ;          FASE 3 */
    /* ;---------------------------------- */
    /* ; Selezione dati per formare RHI */
    /* ;---------------------------------- */

/*     for(k=0;k<MAX_BIN;k++){
        beamXweight[k]=(float **) malloc(w_x_size*sizeof(float *));
        for(i=0;i<w_x_size;i++){
            beamXweight[k][i]=(float *) malloc(w_z_size*sizeof(float));
        }
    }
*/
    Matrix2D<double> RHI_beam(volume.size(), max_bin);
    for (unsigned iaz=0; iaz<volume.beam_count; iaz++)
    {
        Matrix2D<double>& rhi_cart = *slices[iaz];
        Matrix2D<double> rhi_weight(Matrix2D<double>::Zero(x_size, z_size));

        volume.read_vertical_slice(iaz, RHI_beam, MISSING_DB);
        /* ;---------------------------------- */
        /* ;          FASE 4 */
        /* ;---------------------------------- */
        /* ;   Costruzione RHI */
        /* ;---------------------------------- */

        // Enrico: non sforare se il raggio è piú lungo di MAX_BIN
        unsigned ray_size = volume.scan(0).beam_size;
        if (ray_size > max_bin)
            ray_size = max_bin;

        for (unsigned iel=0;iel<volume.size();iel++){
            for (unsigned ibin=0;ibin<ray_size;ibin++) {
                double beamXweight[w_x_size][w_z_size];

                for(unsigned kx=0;kx<(unsigned)w_x_size;kx++){
                    for(unsigned kz=0;kz<(unsigned)w_z_size;kz++){
//std::cout<<"ibin , kx, kz "<<ibin<<" "<<kx<<" "<<kz<<" "<<w_x_size<< " "<<w_z_size<<" "<<MAX_BIN<<std::endl;
//std::cout<<"beam "<<  beamXweight[ibin][kx][kz]<<std::endl;
//std::cout<<"RHI "<<  RHI_beam[iel][ibin]<<std::endl;
//std::cout<<"w_tot "<<  w_tot[kx][kz]<<std::endl;
                        beamXweight[kx][kz] = RHI_beam(iel, ibin) * w_tot(kx, kz);
                    }
                }
                unsigned imin = im(ibin, iel);
                unsigned imax = ix(ibin, iel);
                unsigned jmin = jm(ibin, iel);
                unsigned jmax = jx(ibin, iel);

                unsigned wimin=w_x_size_2-i_xx_min(ibin, iel);
                //wimax=w_x_size_2+i_xx_max[ibin][iel];
                unsigned wjmin=w_z_size_2-i_zz_min(ibin, iel);
                //wjmax=w_z_size_2+i_zz_max[ibin][iel];
                for (unsigned i=imin;i<=imax;i++) {
                    for (unsigned j=jmin;j<=jmax;j++) {
                        rhi_cart(i, j) = rhi_cart(i, j) + beamXweight[wimin+(i-imin)][wjmin+(j-jmin)];
                        rhi_weight(i, j) = rhi_weight(i, j)+w_tot(wimin+(i-imin), wjmin+(j-jmin));
                    }
                }
            }
        }
        for (unsigned i=0;i<x_size;i++) {
            for (unsigned j=0;j<z_size;j++) {
                if (rhi_weight(i, j) > 0.0)
                    rhi_cart(i, j)=rhi_cart(i, j)/rhi_weight(i, j);
                else {
                    rhi_cart(i, j)=MISSING_DB;

                }
            }
        }
    }
}

}
