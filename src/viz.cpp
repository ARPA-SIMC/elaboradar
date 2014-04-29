#include "viz.h"
#include "par_class.h"
#include "vpr_par.h"
#include "cylindrical.h"

#define MISSING 0 /*valore mancante*/

using namespace std;

namespace cumbac {

CalcoloVIZ::CalcoloVIZ(const CylindricalVolume& cil, double htbb, double hbbb, double t_ground)
    : cil(cil), x_size(cil.x_size), z_size(cil.z_size), htbb(htbb), hbbb(hbbb), t_ground(t_ground),
      conv_VIZ(x_size, NUM_AZ_X_PPI, MISSING), stratiform(x_size, NUM_AZ_X_PPI, MISSING)
{
    logging_category = log4c_category_get("radar.vpr");
}

void CalcoloVIZ::classifico_VIZ()
{
    int i,j,k,kbbb=0,ktbb=0,kmax=0;
    float cil_Z,base;
    Matrix2D<double> Zabb(x_size, NUM_AZ_X_PPI, 0.);
    Matrix2D<double> Zbbb(x_size, NUM_AZ_X_PPI, 0.);
    double ext_abb,ext_bbb;
    float LIM_VERT= 8.;//questo l'ho messo io
    long int ncv = 0;

    kbbb=floor(hbbb/RES_VERT_CIL);   //08/01/2013...MODIFICA, inserito questo dato
    ktbb=ceil(htbb/RES_VERT_CIL);

    kmax=ceil(LIM_VERT/RES_VERT_CIL);
    // kmax=ceil(z_size/RES_VERT_CIL);
    if (t_ground < T_MAX_ML) kmax=0;/////se t suolo dentro t melting layer pongo kmax=00 e in tal modo non classifico
    if (ktbb>z_size) ktbb=z_size;
    LOG_DEBUG("kmax= %i \n kbbb= %i \n ktbb= %i \n  z_size= %i",kmax,kbbb,ktbb,z_size);

    //inizio l'integrazione
    for(i=0; i<NUM_AZ_X_PPI; i++){
        for(j=0; j<x_size; j++)
        {
            ext_abb=0.;
            ext_bbb=0.;

            //modifica 08/01/2013 .. afggiungo questo ..per fare l'integrazione anche con i dati sotto la bright band
            for(k=0; k<kbbb; k++)
            {
                if (cil[i][j][k] > -19.){   // 08/01/2013..modifica, prendo fin dove ho un segnale
                    base=(cil[i][j][k])/10.;
                    cil_Z=pow(10.,base);
                    Zbbb[i][j] = Zbbb[i][j] + RES_VERT_CIL*cil_Z;
                    ext_bbb=RES_VERT_CIL+ext_bbb;
                }
            }
//std::cout<<"Z_size :"<<z_size<<" kbbb :"<<kbbb<<" ktbb "<<ktbb<<std::endl;
            for(k=kbbb; k<ktbb; k++)
            {
                if (k < 4 ){
                    if (cil[i][j][k]>10. &&  cil[i][j][k+4]> 5.){
                        if (cil[i][j][k] - cil[i][j][k+4] > 5.)
                            stratiform[i][j]=1;
                    }
                }
                else if (cil[i][j][k]>10. &&  cil[i][j][k+4]> 5. &&  cil[i][j][k-4] > 5.){
                    if (cil[i][j][k] - cil[i][j][k+4] > 5.&&   cil[i][j][k]- cil[i][j][k-4] > 5. )
                        stratiform[i][j]=1;

                }


                if (cil[i][j][k] - cil[i][j][k+4] > 5.)
                    stratiform[i][j]=1;

                for(k=ktbb; k<kmax; k++)
                {
                    if (cil[i][j][k] > -19.){    // 08/01/2013..modifica, prendo fin dove ho un segnale
                        base=(cil[i][j][k])/10.;
                        cil_Z=pow(10.,base);
                        Zabb[i][j] = Zabb[i][j] + RES_VERT_CIL*cil_Z;
                        ext_abb=RES_VERT_CIL+ext_abb;
                    }

                }


                //solo se l'estensione verticale del segnale sopra il top della bright band è maggiore di 0.8 Km classifico

                if (ext_bbb +  ext_abb>0.8) {
                    //if ( ext_abb>0.8) {
                    if ((Zabb[i][j] +Zbbb[i][j])/(ext_bbb+ext_abb) > THR_VIZ){
                        //if ((Zabb[i][j] /ext_abb) > THR_VIZ){
                        conv_VIZ[i][j]=CONV_VAL;
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
