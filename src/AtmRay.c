#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <R.h>
#include <Rmath.h>

void proplin(double *p, double *az, const double *zs, const double *zr, const double *c0, const double *gc, const double *wx0, const double *gwx, const double *wy0, const double *gwy, const double *rho0, const double *grho, double *x, double *y, double *t, double *A)
{
  //  Rprintf("p %f, az %f, zs %f, zr %f, c0 %f, gc %f, wx0 %f, gwx %f, wy0 %f, gwy %f, rho0 %f, grho %f, x %f, y %f, t %f, A %f", *p, *az, *zs, *zr, *c0, *gc, *wx0, *gwx, *wy0, *gwy, *rho0, *grho, *x, *y, *t, *A);

  double u0, v0, c0_eff, gu, gv, gc_eff, cs, cr, epss, epsr, rad, T, dpdx, E, lambdar, tan, us, ur, dthdx;
  double pi = M_PI;
  double p_eps = pow(10, -8);
  double gc_eps = pow(10, -5);

  // this fails for p == 0, so set p as 10^-8
  if(fabs(*p) < p_eps){
    *p = p_eps;
  }
  //Rprintf("1\n");
  // Define effective sound speed slope and intercept
  u0 = cos(*az * pi/180) * *wy0 + sin(*az * pi/180) * *wx0;
  v0 = sin(*az * pi/180) * *wy0 - cos(*az * pi/180) * *wx0;
  c0_eff = *c0 + u0;
  gu = cos(*az * pi/180) * *gwy + sin(*az * pi/180) * *gwx;
  gv = sin(*az * pi/180) * *gwy - cos(*az * pi/180) * *gwx;
  gc_eff = *gc + gu;

  // function crashes when gc = 0
  if(fabs(gc_eff) < gc_eps){
    gc_eff = -gc_eps;
  }
  //Rprintf("2\n");

  // Define source/receiver effective sound speed
  //  cs = c0_eff + *zs * gc_eff;
  //  cr = c0_eff + *zr * gc_eff;
  cs = *c0 + *zs * *gc;
  cr = *c0 + *zr * *gc;
  us = u0 + *zs * gu;
  ur = u0 + *zr * gu;
  //Rprintf("3\n");
  
  // These are useful variables to define separately
  //  epss = sqrt(1 - *p * *p * cs * cs);
  //  epsr = sqrt(1 - *p * *p * cr * cr);
  // Calculate displacement, time, amplitude
  // Equations from Jeff's notes for Villarrica paper from Slotnick
  // JFA 9/13/2012: YOU CANNOT USE EFFECTIVE SOUND SPEED IN 3D!!! Rewriting using equations from Garces and numeric integration (farther down).
  //  rad = (*zs - *zr) * (epss - epsr)/(cr - cs)/ *p; // radial distance traveled//BAD EQUATION!
  
  //Rprintf("4\n");
  //Rprintf("%f\n", (*zs - *zr)/(cs - cr));
  //Rprintf("%f\n", (cs * (1 + epsr))/(cr * (1 + epss)));
  //Rprintf("%f\n", log((cs * (1 + epsr))/(cr * (1 + epss))));
  //Rprintf("%f\n", (*zs - *zr)/(cs - cr) * log((cs * (1 + epsr))/(cr * (1 + epss))));
  //*t = fabs((*zs - *zr)/(cs - cr) * log((cs * (1 + epsr))/(cr * (1 + epss)))); // traveltime

  // numerical integration (trapezoidal) of equations from Garces et al.
  *t = 0;
  rad = 0;
  tan = 0;
  int nsteps = 100;
  double dz = (*zs - *zr)/nsteps;
  double zi1, zi2, s1, s2, u1, u2, sqrt1, sqrt2;
  double p2 = *p + pow(10, -9); // for use in computing rad2, which is used to find dp/dx
  double rad2 = 0;
  double wind_offset = 0;
  double wind_offset2 = 0;
  for(int i = 0; i < nsteps; i++){
    zi1 = *zr + i * dz;
    zi2 = *zr + (i+1)*dz;
    s1 = 1/(*c0 + *gc * zi1);
    s2 = 1/(*c0 + *gc * zi2);
    u1 = u0 + gu * zi1;
    u2 = u0 + gu * zi2;
    sqrt1 = sqrt(pow(s1, 2) - *p * *p/pow(1 - *p * u1, 2));
    sqrt2 = sqrt(pow(s2, 2) - *p * *p/pow(1 - *p * u2, 2));
    *t += 0.5 * (pow(s1, 2) / sqrt1 + pow(s2, 2) / sqrt2 ) * fabs(dz);
    rad += 0.5 * ( (*p/(1 - *p * u1) + pow(s1,2) * u1)/sqrt1 + (*p/(1 - *p * u2) + pow(s2,2) * u2)/sqrt2 ) * fabs(dz);
    tan += 0.5 * ( pow(s1, 2) * (v0 + gv * zi1)/sqrt1 + pow(s2, 2) * (v0 + gv * zi2)/sqrt2) * fabs(dz);
    wind_offset += 0.5 * (u1 * pow(s1, 2) / sqrt1 + u2 * pow(s2, 2) / sqrt2 ) * fabs(dz); // "drift" from wind, used for energy equation.
    sqrt1 = sqrt(pow(s1, 2) - p2 * p2/pow(1 - p2 * u1, 2));
    sqrt2 = sqrt(pow(s2, 2) - p2 * p2/pow(1 - p2 * u2, 2));
    rad2 += 0.5 * ( (p2/(1 - p2 * u1) + pow(s1,2) * u1)/sqrt1 + (p2/(1 - p2 * u2) + pow(s2,2) * u2)/sqrt2 ) * fabs(dz);
    wind_offset2 += 0.5 * (u1 * pow(s1, 2) / sqrt1 + u2 * pow(s2, 2) / sqrt2 ) * fabs(dz); // "drift" from wind, used for energy equation.
  }
    //  Rprintf("%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f\n", rad, tan, *t, zi1, zi2, s1, s2, u1, u2, sqrt1, sqrt2);
		    
  // Calculate amplitude (several steps)
  //Rprintf("%f, %f, %f, %f, %f, %f, %f\n", cs, cr, *p, epsr, epss, *zs, *zr);
  //dpdx = (cs - cr) * *p * *p * epsr * epss/ ((*zs - *zr) * (epss * epss * epsr - epsr * epsr * epss + epsr * *p * *p * cs * cs - epss * *p * *p * cr * cr)); //d(ray parameter)/d(radial distance): used to calculate energy
  //Rprintf("5\n");
    
  //  dpdx = (ps1 - ps2)/(rad - rad2);
  //E = *p * cs * cs/(rad * (1 - *p * *p * cs * cs)) * dpdx; // arrival energy densiy
  //  E = ps1 * cs * cs/(rad * (1 - ps1 * ps1 * cs * cs)) * dpdx; // arrival energy density

  //  double ths = atan2(cs * cs/(1/ *p - us) + us, cs * sqrt(1-pow(cs/(1/ *p - us), 2)));
  //  double thr = atan2(cr * cr/(1/ *p - ur) + ur, cr * sqrt(1-pow(cr/(1/ *p - ur), 2)));
  //  double ths2 = atan2(cs * cs/(1/ p2 - us) + us, cs * sqrt(1-pow(cs/(1/p2 - us), 2)));

  double ths = asin(cs/(1/ *p - us));
  double ths2 = asin(cs/(1/p2 - us));
  double thr = asin(cr/(1/ *p - ur));

  dthdx = (ths - ths2)/(rad - rad2 + wind_offset2 - wind_offset);
  E = sin(ths) * dthdx/(4 * pi * (rad - wind_offset) * cos(thr));
  
  //  Rprintf("%f, %f, %f, %f, %f, %f, %f, %f\n", ths, thr, ths2, dthdx, wind_offset, rad, rad2, E);

  //Rprintf("6\n");
  lambdar = pow(*c0 + *gc * *zr, 2) * (*rho0 + *grho * *zr); // the lame' parameter at the receiver
  *A = sqrt(lambdar * E); // arrival amplitude

  // calculate transverse motion (ray can travel perpendicular to azimuth with a transverse wind)
	  //  tan = gv/(gc_eff*gc_eff) * asin(cs * *p)/ *p + (c0_eff * gv - v0 * gc_eff)/(gc_eff * gc_eff) * log((1 / *p + sqrt(1/(*p * *p) - cs*cs))/cs) - (gv/(gc_eff * gc_eff) * asin(cr * *p)/ *p + (c0_eff*gv - v0*gc_eff)/(gc_eff * gc_eff) * log((1/ *p + sqrt(1/(*p * *p) - cr*cr))/cr));  // tangential motion: did this by hand over multiple days--int v dT = int v dT/dz dz
  //Rprintf("7\n");
 
  // Finally: rotate rad, tan to find outputs x and y
  *x = sin(*az * pi/180) * rad - cos(*az * pi/180) * tan;
  *y = sin(*az * pi/180) * tan + cos(*az * pi/180) * rad;
  //Rprintf("8\n");
  //  Rprintf("%f %f %f %f %f\n", *x, *y, *t, *A, *p);
  if(isnan(*t)){
    *x = 1.0/0.0;
    *y = 1.0/0.0;
    *t = 1.0/0.0;
    *A = 0.0/0.0;
  }
    //Rprintf("%f %f %f %f %f\n", *x, *y, *t, *A, *p);
      
}







/////////////////////////////////////////////
//////////////////////////////////////////////
///////////////////////////////////////////////
//////////////////////////////////////////////



void p4xlin(double *x, double *y, const double *zs, const double *zr, const double *c0, const double *gc, const double *wx0, const double *gwx, const double *wy0, const double *gwy, const double *rho0, const double *grho, double *maxerror, double *p, double *az, double *error){
 
  double cr, cs, A, B, C, D, Af, Bf, Cf, Df, Amp, t, mx, my, R0, R1, R0f, R1f, Mr, minmin;
  double pi = M_PI;
  // initial values for search variables
  double ptol = 0.00001;
  double aztol = 0.1;
  double phi = (sqrt(5.0)-1.0)/2.0;
  int n = 0;
  *error = 1.0/0.0;

  *az = atan2(*x, *y) * 180/pi;

  // we can use a shortcut if wind and dc/dz are zero: just propagate in a straight line
  double eps = 0.000001;
  if(fabs(*gc) < eps && fabs(*wx0) < eps && fabs(*wy0) < eps && fabs(*gwx) < eps && fabs(*gwy) < eps){
    *p = sqrt(*x * *x + *y * *y)/sqrt(*x * *x + *y * *y + pow(*zs - *zr, 2))/ *c0;
    proplin(p, az, zs, zr, c0, gc, wx0, gwx, wy0, gwy, rho0, grho, &mx, &my, &t, &Amp); // this call changes mx and my

    *error = sqrt(pow(mx - *x, 2) + pow(my - *y, 2));
    return;
  }

  // loop until error is acceptable
  while(*error > *maxerror){
    //Rprintf("%f, %f, %f, %f, %f, \n", *p, *az, *error, mx, my);
    /* Could be that the target is in a shadow zone, so there is no
       solution.  Break after 10 iterations to avoid infinite loops
       (no calculation should ever take that many).
    */
    if(n > 10){ // return NaNs
      *p = 0.0/0.0;
      *az = 0.0/0.0;
      *error = 1.0/0.0;
      return;
    }
    n++;
    // define search range for p
    cr = *c0 + *wx0 * sin(*az * pi/180) + *wy0 * cos(*az * pi/180) + *zr * (*gc + *gwx * sin(*az * pi/180) + *gwy * cos(*az * pi/180));
    cs = *c0 + *wx0 * sin(*az * pi/180) + *wy0 * cos(*az * pi/180) + *zs * (*gc + *gwx * sin(*az * pi/180) + *gwy * cos(*az * pi/180));

    // use a golden section search for each variable.
    // usually, the golden search method is described in terms of
    // 2 exterior points and 1 interior points.  I'm finding the
    // bookkeeping easier to keep track of 2 exterior and 2 interior
    // points, without loss of efficiency.  A is far left, then B, C,
    // and D.  If A is 0 and D is 1, B is .38 (1-phi) and C is 0.62 (phi).

    // Take the interior points with the higher functional value and make it
    // an exterior point.  Relabel the others accordingly, and create a new
    // interior point in the proper place.
    
    A = 0; // minimum p to consider
    // find the maximum p--1/cmax (anything greater can't propagate at least somewhere in the zs--zr range.
    if(cs > cr){
      D = 1/(cs + 1e-12);
    }else{
      D = 1/(cr + 1e-12);
    }
    //    Rprintf("cs %f, cr %f, D %f\n", cs, cr, D);
    // define the interior points
    B = A + (D-A) * (1-phi);
    C = A + (D-A) * phi;

    // find the functional values of interior points
    //proplin(&A, az, zs, zr, c0, gc, wx0, gwx, wy0, gwy, rho0, grho, &mx, &my, &t, &Amp); // this call changes *mx and *my
    //    Af = sqrt(pow(mx - *x, 2) + pow(my - *y, 2));  // don't need Af for the p search
    proplin(&B, az, zs, zr, c0, gc, wx0, gwx, wy0, gwy, rho0, grho, &mx, &my, &t, &Amp); // this call changes *mx and *my
    Bf = sqrt(pow(mx - *x, 2) + pow(my - *y, 2));
    proplin(&C, az, zs, zr, c0, gc, wx0, gwx, wy0, gwy, rho0, grho, &mx, &my, &t, &Amp); 
    Cf = sqrt(pow(mx - *x, 2) + pow(my - *y, 2));
    proplin(&D, az, zs, zr, c0, gc, wx0, gwx, wy0, gwy, rho0, grho, &mx, &my, &t, &Amp); // this call changes *mx and *my
    Df = sqrt(pow(mx - *x, 2) + pow(my - *y, 2));
    //Rprintf("%f %f %f %f\n", A, B, C, D);
    // R0, R1 are used for error constraint (see below)
    R0 = D;
    R0f = Df;
    R1 = 1/(1/D + 1e-11); // very slightly less than D
    proplin(&R1, az, zs, zr, c0, gc, wx0, gwx, wy0, gwy, rho0, grho, &mx, &my, &t, &Amp); // this call changes *mx and *my
    R1f = sqrt(pow(mx - *x, 2) + pow(my - *y, 2));

  // loop through p until diff p is acceptable
    *error = 1.0/0.0;
    while(*error > *maxerror){
      if(Cf < Bf){
        A = B;
        B = C;
        Bf = Cf;
        C = A + D - B;
	proplin(&C, az, zs, zr, c0, gc, wx0, gwx, wy0, gwy, rho0, grho, &mx, &my, &t, &Amp); 
	Cf = sqrt(pow(mx - *x, 2) + pow(my - *y, 2));

      }else{ // C is to the right of the minimum, so we can update R0, R1
        D = C;
	Df = Cf;
        C = B;
        Cf = Bf;
        B = A + D - C;
	proplin(&B, az, zs, zr, c0, gc, wx0, gwx, wy0, gwy, rho0, grho, &mx, &my, &t, &Amp); 
	Bf = sqrt(pow(mx - *x, 2) + pow(my - *y, 2));
        R0 = R1;
        R0f = R1f;
        R1 = D;
        R1f = Df;

      }
      // calculate maximum error.  This assumes that derror/dp increases monotonely on the right side (left side slope becomes steeper as you approach the correct p!), so that the true min cannot be below the right tangent line for x = B.
      Mr = (R0f - R1f)/(R0-R1);
      minmin = Mr * (A - R1) + R1f;
      *error = fmax2(Bf, Cf) - minmin;
    }
    // at this point, we've narrowed down on a p value
    *p = (B + C)/2;
    
    // for azimuth, we center our search interval at *az (what it would be without wind).  error should be small close to *az and almost maximized at *az +/- 180.  So, our search interval is (*az - 180, *az + 180).
    A = *az - 180;
    D = *az + 180;
    B = A + (D-A) * (1-phi);
    C = A + (D-A) * phi;
    // find the functional values of these points
    proplin(p, &B, zs, zr, c0, gc, wx0, gwx, wy0, gwy, rho0, grho, &mx, &my, &t, &Amp); // this call changes *mx and *my
    Bf = sqrt(pow(mx - *x, 2) + pow(my - *y, 2));
    proplin(p, &C, zs, zr, c0, gc, wx0, gwx, wy0, gwy, rho0, grho, &mx, &my, &t, &Amp); 
    Cf = sqrt(pow(mx - *x, 2) + pow(my - *y, 2));

    // loop through az until diff az is acceptable
    while(D-A > aztol){
      if(Cf < Bf){
        A = B;
        B = C;
        Bf = Cf;
        C = A + D - B;
	proplin(p, &C, zs, zr, c0, gc, wx0, gwx, wy0, gwy, rho0, grho, &mx, &my, &t, &Amp); 
	Cf = sqrt(pow(mx - *x, 2) + pow(my - *y, 2));

      }else{
        D = C;
        C = B;
        Cf = Bf;
        B = A + D - C;
	proplin(p, &B, zs, zr, c0, gc, wx0, gwx, wy0, gwy, rho0, grho, &mx, &my, &t, &Amp); // this call changes *mx and *my
	Bf = sqrt(pow(mx - *x, 2) + pow(my - *y, 2));
      }
    }
    
    *az = (C+B)/2;
    proplin(p, az, zs, zr, c0, gc, wx0, gwx, wy0, gwy, rho0, grho, &mx, &my, &t, &Amp); // this call changes mx and my
    //    Rprintf("%f, %f\n", mx, my);
    *error = sqrt(pow(mx - *x, 2) + pow(my - *y, 2));
  }
  //  Rprintf("%i\n", n);
}      













////////////////////////////////////////////////////    
////////////////////////////////////////////////////    
////////////////////////////////////////////////////    
////////////////////////////////////////////////////    
////////////////////////////////////////////////////    
////////////////////////////////////////////////////    




















void makearrivals(const double *xs, const double *ys, const double *zs, const double *xr, const double *yr, const double *zr, const double *dt, const int *nt, const double *timing, const double *c0, const double *gc, const double *wx0, const double *gwx, const double *wy0, const double *gwy, const double *rho0, const double *grho, int *nsrc, int *nrec, double *P){
  // P must be defined from calling function as rep(0, nt*nrec)
  //Rprintf("nrec %i, nsrc %i\n", *nrec, *nsrc);
  // initialize variables
  //  t = 1:nt * dt
  int src=0;
  int rec=0;
  double dx, dy, p, az, error, mx, my, t, Amp;
  double maxerror = 3;
  int ai;

  //maxt = (nt+1) * dt
  // loop through source indices
    // loop through receiver indices
  for(rec = 0; rec < *nrec; rec++){
    for(src = 0; src < *nsrc; src++){

      //Rprintf("src %i, rec %i, nrec %i, nsrc %i, rec<nrec %i\n", src, rec, *nrec, *nsrc, rec < *nrec);
  // calculate arrival times and amplitudes
      dx = xr[rec] - xs[src];
      dy = yr[rec] - ys[src];
      //Rprintf("before p4xlin: dx %f, dy %f, zs %f, zr %f, c0 %f, gc %f, wx0 %f, gwx %f, wy0 %f, gwy %f, rho0 %f, grho %f, maxerror %f, p %f, az %f, error %f\n", dx, dy, zs[src], zr[rec], *c0, *gc, *wx0, *gwx, *wy0, *gwy, *rho0, *grho, maxerror, p, az, error);
      p4xlin(&dx, &dy, zs +src, zr + rec, c0, gc, wx0, gwx, wy0, gwy, rho0, grho, &maxerror, &p, &az, &error);
      //Rprintf("before proplin: dx %f, dy %f, zs %f, zr %f, c0 %f, gc %f, wx0 %f, gwx %f, wy0 %f, gwy %f, rho0 %f, grho %f, maxerror %f, p %f, az %f, error %f\n", dx, dy, *zs, *zr, *c0, *gc, *wx0, *gwx, *wy0, *gwy, *rho0, *grho, maxerror, p, az, error);
       proplin(&p, &az, zs+src, zr+rec, c0, gc, wx0, gwx, wy0, gwy, rho0, grho, &mx, &my, &t, &Amp);
       //      Rprintf("t %f\n", t);
      
      // update P with the new arrival
      if(t < (*nt-1) * *dt){
	ai = 0.5 + (t+timing[src])/ *dt; // arrival index--this coerces it to an integer by truncation; added 0.5 because trunc(x+0.5) = round(x)
	P[ai + rec * *nt] += Amp;
	//Rprintf("t %f, ai %i\n", t, ai);
      }
    }
  }
}
