#ifndef ARCHIVIATORE_STEINER_CLASS_H
#define ARCHIVIATORE_STEINER_CLASS_H

#include "logging.h"
#include "volume.h"
#include "volume/elev_fin.h"

namespace elaboradar {

namespace steiner {

struct Point
{
    int azimut;
    int range;
    unsigned npoints;
    // Valore della Z di background in dB
    double bckgr;
    // Valore della Z di background in mm^6/m^3
    double Z_bckgr;
    double convective_radius;

    Point() : azimut(-999), range(-999), npoints(0), bckgr(0), Z_bckgr(0), convective_radius(0) {}
    Point(int azimut, int range) : azimut(azimut), range(range), npoints(0), bckgr(0), Z_bckgr(0), convective_radius(0) {}

    void add_sample(double sample);
    void finalize();
};

}

struct CalcoloSteiner
{
    log4c_category_t* logging_category;

    const Volume<double>& volume;
    const volume::ElevFin<double>& elev_fin;
    const unsigned max_bin;
    const double size_cell;

    Matrix2D<unsigned char> conv_STEINER;

    // Pixel precipitanti
    std::vector<steiner::Point> lista_bckg;

    CalcoloSteiner(const Volume<double>& volume, const volume::ElevFin<double>& elev_fin, unsigned max_bin);

    /**
     *  calcola valore di background per individuare pixel convettivo
     *
     *  @brief funzione  che calcola il background 
     *  @details la classificazione di Steiner non ha bisogno di ricampionamento cilindrco perciò  uso direttamente la matrice polare
     */
    void calcolo_background();

    /**
     *  @brief funzione  che classifica secondo STEINER
     *  @details segna come convettivi i punti che hanno valore superiore a 40 dBZ e differenza col background elevata, quindi ingrandisce i nuclei di un raggio variabile
     */
    void classifico_STEINER();

    /**
     *  @brief funzione  che ingrandisce i nuclei di Steiner
     *  @details ingrandisce i nuclei di Steiner di un valore pari al raggio convettivo
     *  @param[in] cr raggio convettivo
     *  @param[in] ja indice di azimut
     *  @param[in] kr indice di range
     */ 
    void ingrasso_nuclei(float cr,int ja,int kr);

    void add_sample(unsigned pos, unsigned azimut, unsigned range);
};

}

#endif
