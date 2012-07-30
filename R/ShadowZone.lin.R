ShadowZone.lin = function(zr, zs, az = 0, ATM){

  c0eff = ATM$c0 + ATM$wx0 * sin(az * pi/180) + ATM$wy0 * cos(az * pi/180)
  gceff = ATM$gc + ATM$gwx * sin(az * pi/180) + ATM$gwy * cos(az * pi/180)
  cr = ATM$c0 + ATM$wx0 * sin(az * pi/180) + ATM$wy0 * cos(az * pi/180) + zr * (ATM$gc + ATM$gwx * sin(az * pi/180) + ATM$gwy * cos(az * pi/180))
  cs = ATM$c0 + ATM$wx0 * sin(az * pi/180) + ATM$wy0 * cos(az * pi/180) + zs * (ATM$gc + ATM$gwx * sin(az * pi/180) + ATM$gwy * cos(az * pi/180))

  p = 1/max(cr, cs)

  X = sqrt(abs(2*(zr-zs)/gceff/p) - (zr-zs)^2) # no energy will arrive beyond this distance (without turning)
  
  return(list(X=X, p = p))
}
