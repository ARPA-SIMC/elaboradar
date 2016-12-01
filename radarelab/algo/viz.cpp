#include <radarelab/algo/viz.h>
#include <radarelab/par_class.h>
#include <radarelab/vpr_par.h>
#include <radarelab/cylindrical.h>

#define MISSING 0 /*valore mancante*/

using namespace std;

namespace radarelab {
namespace algo {

CalcoloVIZ::CalcoloVIZ(const CylindricalVolume& cil, double htbb, double hbbb, double t_ground)
    : cil(cil), x_size(cil.x_size), z_size(cil.z_size), htbb(htbb), hbbb(hbbb), t_ground(t_ground), res_vert_cil (cil.resol[1]),
      conv_VIZ(Matrix2D<unsigned char>::Constant(cil.slices.size(), x_size, MISSING)),
      stratiform(Matrix2D<unsigned char>::Constant(cil.slices.size(), x_size, MISSING))
{
    logging_category = log4c_category_get("radar.vpr");
}

void CalcoloVIZ::classifico_VIZ()
{
    float cil_Z,base;
    Matrix2D<double> Zabb(Matrix2D<double>::Zero(cil.slices.size(), x_size));
    Matrix2D<double> Zbbb(Matrix2D<double>::Zero(cil.slices.size(), x_size));
    double ext_abb,ext_bbb;
    float LIM_VERT= 8.;//questo l'ho messo io
    long int ncv = 0;

    unsigned kbbb=floor(hbbb/res_vert_cil);   //08/01/2013...MODIFICA, inserito questo dato
    unsigned ktbb=ceil(htbb/res_vert_cil);
    unsigned kmax=ceil(LIM_VERT/res_vert_cil);
    // kmax=ceil(z_size/res_vert_cil);

    if (t_ground < T_MAX_ML) kmax=0;/////se t suolo dentro t melting layer pongo kmax=00 e in tal modo non classifico
    if (ktbb>z_size) ktbb=z_size;
    LOG_DEBUG("kmax= %i \n kbbb= %i \n ktbb= %i \n  z_size= %i",kmax,kbbb,ktbb,z_size);

    //inizio l'integrazione
    for(unsigned i=0; i<cil.slices.size(); i++){
        for(unsigned j=0; j<x_size; j++)
        {
            ext_abb=0.;
            ext_bbb=0.;

            //modifica 08/01/2013 .. afggiungo questo ..per fare l'integrazione anche con i dati sotto la bright band
            for(unsigned k=0; k<kbbb; k++)
            {
                if (cil(i, j, k) > -19.){   // 08/01/2013..modifica, prendo fin dove ho un segnale
                    base=(cil(i, j, k))/10.;
                    cil_Z=pow(10.,base);
                    Zbbb(i, j) = Zbbb(i, j) + res_vert_cil*cil_Z;
                    ext_bbb=res_vert_cil+ext_bbb;
                }
            }
//std::cout<<"Z_size :"<<z_size<<" kbbb :"<<kbbb<<" ktbb "<<ktbb<<std::endl;
            for(unsigned k=kbbb; k<ktbb; k++)
            {
                if (k < 4 ){
                    if (cil(i, j, k)>10. &&  cil(i, j, k+4)> 5.){
                        if (cil(i, j, k) - cil(i, j, k+4) > 5.)
                            stratiform(i, j)=1;
                    }
                }
                else if (cil(i, j, k)>10. &&  cil(i, j, k+4)> 5. &&  cil(i, j, k-4) > 5.){
                    if (cil(i, j, k) - cil(i, j, k+4) > 5.&&   cil(i, j, k)- cil(i, j, k-4) > 5. )
                        stratiform(i, j)=1;

                }

                if (cil(i, j, k) - cil(i, j, k+4) > 5.)
                    stratiform(i, j)=1;

                for(k=ktbb; k<z_size; k++)
                {
                    if (cil(i, j, k) > -19.){    // 08/01/2013..modifica, prendo fin dove ho un segnale
                        base=(cil(i, j, k))/10.;
                        cil_Z=pow(10.,base);
                        Zabb(i, j) = Zabb(i, j) + res_vert_cil*cil_Z;
                        ext_abb=res_vert_cil+ext_abb;
                    }
                }


                //solo se l'estensione verticale del segnale sopra il top della bright band è maggiore di 0.8 Km classifico

                if (ext_bbb +  ext_abb>0.8) {
                    //if ( ext_abb>0.8) {
                    if ((Zabb(i, j) +Zbbb(i, j))/(ext_bbb+ext_abb) > THR_VIZ){
                        //if ((Zabb(i, j) /ext_abb) > THR_VIZ){
                        conv_VIZ(i, j)=CONV_VAL;
                        //lista_conv[ncv][0]= i;
                        //lista_conv[ncv][1]= j;
                        ncv=ncv+1;
                    }
                }

            }
        }
    }
    LOG_DEBUG("numero nuclei VIZ = %li",ncv);

    return;
}

}
}
