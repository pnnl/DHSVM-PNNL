/*
 * SUMMARY:      SeparateRadiation.c - Separate radiation
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-96
 * DESCRIPTION:  Separate observed downward solar radiation into diffuse and 
 *               direct beam radiation based on the clearness index (kt)
 * DESCRIP-END.
 * FUNCTIONS:    SeparateRadiation()
 * COMMENTS:
 * $Id: SeparateRadiation.c,v 1.4 2003/07/01 21:26:24 olivier Exp $     
 */

#include <math.h>
#include "settings.h"
#include "rad.h"

/*****************************************************************************
  Function name: SeparateRadiation()

  Purpose      : Separate observed downward solar radiation into diffuse and 
                 direct beam radiation based on the clearness index (kt)

  Required     :
    float TotalSolar - Total amount of solar (shortwave) radiation (W/m2)
    float ClearIndex - Clearness index
    float *Beam      - Direct beam radiation (W/m2)
    float *Diffuse   - Diffuse radiation (W/m2)

  Returns      : void

  Modifies     :
    float *Beam      - Direct beam radiation (W/m2)
    float *Diffuse   - Diffuse radiation (W/m2)

  Comments     : Based on: Erbs, D.G., S.A. Klein, and J.A. Duffie, 
                 "Estimation of the diffuse fraction for hourly, daily and 
                 monthly-average global radiation", Solar Energy, V.28, n.4, 
                 pp.293-302, 1982
*****************************************************************************/
void SeparateRadiation(float TotalSolar, float ClearIndex,
		       float *Beam, float *Diffuse)
{
  /* The following relationships were developed for hourly data, and are here
     used varying timesteps (Time.Dt).  It would probably be better to
     develop explicit relationships for different timesteps */
  /* if (ClearIndex <= 0.22) {
     *Diffuse = (1 - 0.09 * ClearIndex);
     } 
     else if (ClearIndex <= 0.80) {
     *Diffuse = 0.9511 - 0.1604 * ClearIndex + 
     4.388 * ClearIndex*ClearIndex -
     16.638 * ClearIndex*ClearIndex*ClearIndex + 
     12.366 * ClearIndex*ClearIndex*ClearIndex*ClearIndex;
     } else
     *Diffuse = 0.165;

     *Diffuse *= TotalSolar;
     *Beam     = TotalSolar - *Diffuse; */

  /* The following relationships were taken from Chen and Black or one
     of the papers that Chen and Black reference, for application in the PNW */
  /* The clear index is with respect to top of atmosphere radiation */
  if (ClearIndex > 0.8) {
    *Diffuse = TotalSolar * 0.13;
  }
  else {
    *Diffuse = TotalSolar * (0.943 + 0.734 * ClearIndex -
			     4.9 * ClearIndex * ClearIndex +
			     1.796 * ClearIndex * ClearIndex * ClearIndex +
			     2.058 * ClearIndex * ClearIndex * ClearIndex *
			     ClearIndex);
  }

  *Beam = TotalSolar - *Diffuse;

}
