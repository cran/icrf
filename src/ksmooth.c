/*  MODIFIED ksmooth!!!
 *
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1998-2016	The R Foundation
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 */

#include <math.h>
#include <R.h>			/* for NA_REAL, includes math.h */
#include <Rinternals.h>
//#define RF_DEBUG

double normCDF(double x) {
  return erf(x * 0.7071068) / 2.0 + 0.5;
}

// x: original time, y: original density, n: length(x),
// xp: grid timepoints, y: a density vector to be returned, np: length(xp)
// bw: bandwidth
void ksmooth(double *x, double *y, int n,
             double *xp, double *yp, int np,
             double bw)
{
  int imin = 0;
  double cutoff = 0.0, num, den, x0, w, lb, ub, lbInt;

  /* bandwidth is in units of half inter-quartile range. */
  //bw *= 0.3706506;
  cutoff = 4*bw;
  int smooth = bw > 0.0001;

  while(x[imin] < xp[0] - cutoff && imin < n) imin++;

  yp[0] = 0.0;  // f(0) = 0.0 by default.
  for(int j = 1; j < np; j++) {
    num = den = 0.0;
    x0 = xp[j];
    lb = smooth? x0 - cutoff: xp[j-1];
    ub = smooth? x0 + cutoff: x0;

#ifdef RF_DEBUG
    int m = 10;
    //if (j==0) Rprintf("j=%d, x0 = %2.2f, y[j] = %2.3f, (w, y[i], num, den) = \n", j, x0, y[j]);
    if (j>m & j<m+15) Rprintf("j = %d, imin = %d, lb = %2.3f, x0 = %2.3f, ub = %2.3f\n", j, imin, lb, x0, ub);
#endif

    for(int i = imin; i < n; i++) {

      if(x[i] < lb) imin = i;
      else {

        if(x[i] >= ub) {
          lbInt = (x[i-1] >= lb) ? x[i-1] : lb;
          w = smooth? normCDF(cutoff/bw) - normCDF((lbInt - x0)/bw) : ub - lbInt;
          num += y[i] * w;
          den += w;
#ifdef RF_DEBUG
          if (j>m & j<m+15)   Rprintf("  [60] x[%d] = %2.3f, y[i] = %2.3f w = %2.4f num = %2.4f\n", i, x[i], y[i], w, num);
#endif
          break;
        } else {
          if (x[i-1] < lb) {
            w = smooth? normCDF((x[i] - x0)/bw) - normCDF(-cutoff/bw): (x[i] - lb);
#ifdef RF_DEBUG
            if (j>m & j<m+15)  Rprintf("     [65] (x[%d] - x0)/bw = %2.3f, w = %2.3f, -cutoff/bw = %2.3f, normCDF((x[i] - x0)/bw) = %2.3f, normCDF(-cutoff/bw) = %2.3f\n",
                i, (x[i] - x0)/bw, w, -cutoff/bw, normCDF((x[i] - x0)/bw), normCDF(-cutoff/bw));
#endif
          } else {
            w = smooth? normCDF((x[i] - x0)/bw) - normCDF((x[i - 1] - x0)/bw): (x[i] - x[i - 1]);
          }
          num += y[i] * w;
#ifdef RF_DEBUG
          if (j>m & j<m+15)     Rprintf("  [70] x[%d] = %2.3f, y[i] = %2.3f w = %2.4f num = %2.4f\n", i, x[i], y[i], w, num);
#endif
        }
        //if (j==0) Rprintf("j=%d, imin = %d, (%2.3f,  %2.3f, %2.2f, %2.2f) = \n", j, imin, w, y[i], num, den);
        //den = ub - lb;
        den += w;
      }
    }
    if (den > 0) {
      yp[j] = num/den;
    } else if (num == 0){  // This condition has been added! (tail of surv curve is just zero.)
      yp[j] = 0.0;
    } else {
      yp[j] = 0.0/0.0;
    }
    //if (j<10) Rprintf("   (%2.3f / %2.3f = %2.3f) ", num, den, yp[j]);
#ifdef RF_DEBUG
    if (j>m & j<m+15)  Rprintf("   x0 = %2.2f, y[j]=%2.2f, num = %2.2f, den = %2.2f, yp[j] = %2.2f\n", j, x0, y[j], num, den, yp[j]);
#endif
  }
}

// For use in R.
// normalized ksmooth. prob(y) converted into density, then smoothed density converted into smoothed density (yp)
void ksmoothProb (double *x, double *y, int *n, double *xp, double *yp, int *np, double *bw,
                  int *xinf, int *xpinf, int *sumToOne) {
  int i;
  for (i = 0; i < *n; i++)
    if (xinf[i]) x[i] = 1.0 / 0.0; //convert prob to density.

  for (i = 0; i < *np; i++)
    if (xpinf[i]) xp[i] = 1.0 / 0.0; //convert prob to density.

  for (i = 1; i < *n; i++) {
    y[i] /= (x[i] - x[i-1]);
  }

  ksmooth (x, y, *n, xp, yp, *np, *bw);
  for (i = 1; i < *np; i++) {
    yp[i] *= (xp[i] - xp[i-1]); //convert density to prob.
  }

  if (*sumToOne) {
    // Getting a reference value to make sure density sums to one.
    double *yp2, ypCum = 0.0, yp2Cum = 0.0, ratio = 1.0; //non-smoothed vector with the same gridline.
    yp2   = (double *) S_alloc(*np, sizeof(double));
    ksmooth (x, y, *n, xp, yp2, *np, 0.0);
    for (i = 0; i < *np - 1; i++) { //only up to (*np - 1)th element.
      yp2Cum += yp2[i] *(xp[i] - xp[i-1]); //convert density to prob and cumsum.
      ypCum += yp[i]; // cumsum (already converted to density above.)
    }
    ratio = yp2Cum / ypCum;

    for (i = 0; i < *np - 1; i++) {
      yp[i] *= ratio;
    }
    yp[*np - 1] = 1.0 - yp2Cum;  //yp[last] is the residual.
  }
}
