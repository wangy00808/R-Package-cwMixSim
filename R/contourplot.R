


contour2d = function(opt, method, xlim =NULL, ylim = NULL, mar = c(3.5,3.5,3.5,3.5), nlevels= 100, lwd = 0.2,  col = "red", yaxt = NULL, xaxt = NULL){
  
  if(method =="cwm") contour2d.cwm(opt, xlim =xlim, ylim = ylim, mar = mar, nlevels= nlevels, lwd = lwd,  col = col, yaxt = yaxt, xaxt = xaxt)
  if(method == "skew") contour2d.skew(opt, xlim =xlim, ylim = ylim, mar = mar, nlevels= nlevels, lwd = lwd, col = col, yaxt = yaxt, xaxt =xaxt)
  if(method == "skewcwm") contour2d.skewcwm(opt, xlim =xlim, ylim = ylim, mar = mar, nlevels= nlevels, lwd = lwd,  col = col, yaxt = yaxt, xaxt = xaxt)
    
}
  