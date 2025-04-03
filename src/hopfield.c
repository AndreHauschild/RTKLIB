/*------------------------------------------------------------------------------
* hopfield.c : Hopfield tropospheric model
*
* Copyright (C) 2025 by A.HAUSCHILD, All rights reserved.
*
* references:
*
* version : $Revision:$ $Date:$
* history : 2025/04/03 1.0 new
*-----------------------------------------------------------------------------*/
#include "rtklib.h"

/* troposphere model (Hopfield model) -----------------------------------------
* compute tropospheric delay by simplified Hopfield model
* args   : gtime_t time     I   time
*          double *pos      I   receiver position {lat,lon,h} (rad,m)
*          double *azel     I   azimuth/elevation angle {az,el} (rad)
*          double humi      I   relative humidity
* return : tropospheric delay (m)
*-----------------------------------------------------------------------------*/
extern double tropmodelHpf(gtime_t time, const double *pos, const double *azel,
                           double humi, double *zwd)
{
    double hgt,pres,temp,e,z,trph,trpw;

    if (pos[2]<-100.0||1E4<pos[2]||azel[1]<=0) return 0.0;

    /* standard atmosphere */
    hgt=pos[2]<0.0?0.0:pos[2];

    pres = 1010.0;
    temp = 291.1;
    e = 10.4;

    /* Hopfield model */
    trph = (77.6e-6 * (-613.3768/temp+148.98)*pres)/5.0;
    trpw = (77.6e-6 * 11000.0*4810.0*e/pow(temp,2))/5.0;

    if (zwd) *zwd=trpw;

    return trph;
}

/* troposphere mapping function for Hopfield model -----------------------------
* compute tropospheric mapping function by NMF
* args   : gtime_t t        I   time
*          double *pos      I   receiver position {lat,lon,h} (rad,m)
*          double *azel     I   azimuth/elevation angle {az,el} (rad)
*          double *mapfw    IO  wet mapping function (NULL: not output)
* return : dry mapping function
*-----------------------------------------------------------------------------*/
extern double tropmapfHpf(gtime_t time, const double pos[], const double azel[],
                          double *mapfw)
{
    trace(4,"tropmapfhpf: pos=%10.6f %11.6f %6.1f azel=%5.1f %4.1f\n",
          pos[0]*R2D,pos[1]*R2D,pos[2],azel[0]*R2D,azel[1]*R2D);

    double mapfh;

    if (pos[2]<-1000.0||pos[2]>20000.0) {
        if (mapfw) *mapfw=0.0;
        return 0.0;
    }

     mapfh = 1.0/(sin(sqrt(pow(azel[1],2)+pow(PI/ 72.0,2))));
    *mapfw = 1.0/(sin(sqrt(pow(azel[1],2)+pow(PI/120.0,2))));

    return mapfh;
}
