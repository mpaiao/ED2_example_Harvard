#----- Close all devices and delete all variables. ----------------------------------------#
graphics.off()
rm(list=ls())
#------------------------------------------------------------------------------------------#




#------------------------------------------------------------------------------------------#
#       Settings for this run.                                                             #
#------------------------------------------------------------------------------------------#
#------ Prefix of the daily means. --------------------------------------------------------#
prefix.analysis = "/home/pecan/ed2ws.harvard/analy/harvard"
suffix.analysis = "g01.h5"
#------ First and last day to include in the plots (MM/DD/YYYY). --------------------------#
whena           = "06/10/2000"  # First day
whenz           = "06/20/2000"  # Last day
dwhen           = 1             # Time between outputs, in hours
#------ PFTs that we used. ----------------------------------------------------------------#
include.pft     = c(5,6,8,9,10,11)
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#       Some plot parameters.                                                              #
#------------------------------------------------------------------------------------------#
plot.path   = "/home/pecan/ed2ws.harvard/hourly_plots"
plot.prefix = "harvard"
plot.place  = "Harvard Forest"
plot.ptsz   = 16     # Default font size
plot.width  = 9.7    # Window width
plot.height = 6.0    # Window.height
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#      Handy constants.                                                                    #
#------------------------------------------------------------------------------------------#
Watts.2.MEin = 4.6     #  Watts/m2 -> umol/m2/s (MEin)
t3ple        = 273.16  #  Triple point of water         [K     ]
alvl3        = 2.5e6   #  Latent heat of vap @ t3ple    [J/kg  ]
cph2o        = 1859    #  Specific heat of water vapor  [J/kg/K]
cliq         = 4186    #  Specific heat of liquid water [J/kg/K]
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#       Load packages.                                                                     #
#------------------------------------------------------------------------------------------#
ok = require(chron); if (! ok) stop("Package chron is not available...")
ok = require(hdf5) ; if (! ok) stop("Package hdf5 is not available...")
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#       Source some useful functions.                                                      #
#------------------------------------------------------------------------------------------#
source("timeutils.r")
#------------------------------------------------------------------------------------------#


#----- Avoid beeps. -----------------------------------------------------------------------#
options(locatorBell=FALSE)
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#       Create a vector with all days to be included.                                      #
#------------------------------------------------------------------------------------------#
when   = seq(from=as.numeric(chron(whena)),to=as.numeric(chron(whenz)),by=dwhen/24)
when   = chron(when)
n.when = length(when) 
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#      Create the output path.                                                             #
#------------------------------------------------------------------------------------------#
if (! file.exists(plot.path)) dir.create(plot.path)
#------------------------------------------------------------------------------------------#





#------------------------------------------------------------------------------------------#
#       Set up the arrays we want to plot.                                                 #
#------------------------------------------------------------------------------------------#
#------ 1. Energy budget. -----------------------------------------------------------------#
rshort.net  = rep(NA,times=n.when)
rlong.down  = rep(NA,times=n.when)
rlong.up    = rep(NA,times=n.when)
sensible    = rep(NA,times=n.when)
latent      = rep(NA,times=n.when)
#------ 2. Water fluxes. ------------------------------------------------------------------#
vapor.ca    = rep(0,times=n.when)
vapor.lc    = rep(0,times=n.when)
vapor.wc    = rep(0,times=n.when)
vapor.gc    = rep(0,times=n.when)
transp      = rep(0,times=n.when)
vapor.st    = rep(0,times=n.when)
#------ 3. Light response curve. ----------------------------------------------------------#
gpp         = rep(0,times=n.when)
par.down    = rep(0,times=n.when)
daytime     = rep(FALSE,times=n.when)
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#       Loop over time.                                                                    #
#------------------------------------------------------------------------------------------#
for (w in 1:n.when){

   #------ Make the file name. ------------------------------------------------------------#
   year.now  = sprintf("%4.4i",numyears (when[w]))
   month.now = sprintf("%2.2i",nummonths(when[w]))
   day.now   = sprintf("%2.2i",numdays  (when[w]))
   hour.now  = sprintf("%2.2i",hours    (when[w]))
   minu.now  = sprintf("%2.2i",minutes  (when[w]))
   seco.now  = sprintf("%2.2i",seconds  (when[w]))
   file.now  = paste(prefix.analysis,"-I-",year.now,"-",month.now,"-",day.now,"-"
                    ,hour.now,minu.now,seco.now,"-",suffix.analysis,sep="")
   cat(" - Reading file :",file.now,"...","\n")
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #       Load the data.                                                                  #
   #---------------------------------------------------------------------------------------#
   now = hdf5load(file=file.now,load=FALSE,verbosity=0,tidy=TRUE)
   #---------------------------------------------------------------------------------------#


   #---------------------------------------------------------------------------------------#
   #     Grab the variables.                                                               #
   #---------------------------------------------------------------------------------------#
   #----- Energy. -------------------------------------------------------------------------#
   rshort.net [w] = now$AVG.RSHORT * (1. - now$AVG.ALBEDO)
   rlong.down [w] = now$AVG.RLONG
   rlong.up   [w] = now$AVG.RLONGUP
   sensible   [w] = now$AVG.SENSIBLE.GC + now$AVG.SENSIBLE.LC + now$AVG.SENSIBLE.WC
   latent     [w] = ( ( alvl3 + (cph2o-cliq) * (now$AVG.CAN.TEMP - t3ple) )
                    * ( now$AVG.VAPOR.GC + now$AVG.VAPOR.LC
                      + now$AVG.VAPOR.WC + now$AVG.TRANSP   ) )
   #----- Water. --------------------------------------------------------------------------#
   vapor.ca   [w] = - now$AVG.VAPOR.AC * 86400
   vapor.lc   [w] =   now$AVG.VAPOR.LC * 86400
   vapor.wc   [w] =   now$AVG.VAPOR.WC * 86400
   vapor.gc   [w] =   now$AVG.VAPOR.GC * 86400
   transp     [w] =   now$AVG.TRANSP   * 86400
   vapor.st   [w] = ( now$AVG.VAPOR.AC + now$AVG.VAPOR.LC + now$AVG.VAPOR.WC
                    + now$AVG.VAPOR.GC + now$AVG.TRANSP   ) * 86400
   #----- Light response curve. -----------------------------------------------------------#
   par.down   [w] = ( now$AVG.PAR.BEAM + now$AVG.PAR.DIFFUSE ) * Watts.2.MEin
   gpp        [w] =   now$AVG.GPP
   daytime    [w] =   now$AVG.RSHORT.DIFFUSE > 1.0
   #---------------------------------------------------------------------------------------#
}#end for
#------------------------------------------------------------------------------------------#





#------------------------------------------------------------------------------------------#
#     Plot the energy cycle.                                                               #
#------------------------------------------------------------------------------------------#
pdf(file=file.path(plot.path,"energy.pdf"),onefile=FALSE,pointsize=plot.ptsz
   ,width=plot.width,height=plot.height)


   #----- Plot annotation. ----------------------------------------------------------------#
   my.title      = paste("Energy fluxes at ",plot.place,sep="")
   my.xlab       = "Time [UTC]"
   my.ylab       = "Fluxes [W/m2]"
   #---------------------------------------------------------------------------------------#



   #----- Split the wind so the legend stays outside the plot. ----------------------------#
   layout(mat=rbind(1,2),heights=c(5,1))
   #---------------------------------------------------------------------------------------#



   #----- Find convenient scales for x and y, and nice time stamps. -----------------------#
   xlimit   = range(when)
   ylimit   = range(c(rshort.net,rlong.down,rlong.up,sensible,latent))
   whenplot = pretty.time(when,n=10)
   #---------------------------------------------------------------------------------------#



   #----- Main plot. ----------------------------------------------------------------------#
   par(mar=c(4.1,4.1,4.1,2.1),bg="white")
   plot.new()
   plot.window(xlim=xlimit,ylim=ylimit)
   title(main=my.title,xlab=my.xlab,ylab=my.ylab)
   box()
   axis(side=1,at=whenplot$levels,labels=whenplot$labels,padj=whenplot$padj)
   axis(side=2)
   abline(v=whenplot$levels,h=axTicks(side=2),col="grey50",lty="dotted")
   #----- Add components. -----------------------------------------------------------------#
   lines(x=when,y=rshort.net,lwd=2.5,col="darkgoldenrod")
   lines(x=when,y=rlong.down,lwd=2.5,col="slateblue3"   )
   lines(x=when,y=rlong.up  ,lwd=2.5,col="darkorange"   )
   lines(x=when,y=sensible  ,lwd=2.5,col="firebrick"    )
   lines(x=when,y=latent    ,lwd=2.5,col="royalblue4"   )
   #---------------------------------------------------------------------------------------#



   #----- Legend. -------------------------------------------------------------------------#
   par(mar=c(0.1,4.1,0.1,2.1))
   plot.new()
   plot.window(xlim=c(0,1),ylim=c(0,1))
   legend(x="bottom",inset=0,lwd=2.5,ncol=3
         ,legend=c("Net SW","Down LW","Up LW","Sensible","Latent")
         ,col=c("darkgoldenrod","slateblue3","darkorange","firebrick","royalblue4"))
   #---------------------------------------------------------------------------------------#
dev.off()
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#





#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#     Plot the water cycle.                                                                #
#------------------------------------------------------------------------------------------#
pdf(file=file.path(plot.path,"water.pdf"),onefile=FALSE,pointsize=plot.ptsz
   ,width=plot.width,height=plot.height)


   #----- Plot annotation. ----------------------------------------------------------------#
   my.title      = paste("Water fluxes at ",plot.place,sep="")
   my.xlab       = "Time [UTC]"
   my.ylab       = "Fluxes [kg/m2/day]"
   #---------------------------------------------------------------------------------------#



   #----- Split the wind so the legend stays outside the plot. ----------------------------#
   layout(mat=rbind(1,2),heights=c(5,1))
   #---------------------------------------------------------------------------------------#



   #----- Find convenient scales for x and y, and nice time stamps. -----------------------#
   xlimit   = range(when)
   ylimit   = range(c(vapor.ca,vapor.lc,vapor.wc,vapor.gc,transp,vapor.st))
   whenplot = pretty.time(when,n=10)
   #---------------------------------------------------------------------------------------#



   #----- Main plot. ----------------------------------------------------------------------#
   par(mar=c(4.1,4.1,4.1,2.1),bg="white")
   plot.new()
   plot.window(xlim=xlimit,ylim=ylimit)
   title(main=my.title,xlab=my.xlab,ylab=my.ylab)
   box()
   axis(side=1,at=whenplot$levels,labels=whenplot$labels,padj=whenplot$padj)
   axis(side=2)
   abline(v=whenplot$levels,h=axTicks(side=2),col="grey50",lty="dotted")
   #----- Add components. -----------------------------------------------------------------#
   lines(x=when,y=vapor.ca  ,lwd=2.5,col="royalblue3"   )
   lines(x=when,y=vapor.lc  ,lwd=2.5,col="green2"       )
   lines(x=when,y=vapor.wc  ,lwd=2.5,col="darkgoldenrod")
   lines(x=when,y=vapor.gc  ,lwd=2.5,col="sienna4"      )
   lines(x=when,y=transp    ,lwd=2.5,col="chartreuse4"  )
   lines(x=when,y=vapor.st  ,lwd=2.5,col="purple4"      )
   #---------------------------------------------------------------------------------------#



   #----- Legend. -------------------------------------------------------------------------#
   par(mar=c(0.1,4.1,0.1,2.1))
   plot.new()
   plot.window(xlim=c(0,1),ylim=c(0,1))
   legend(x="bottom",inset=0,lwd=2.5,ncol=3
         ,legend=c("CAS->Atm","Leaf->CAS","Wood->CAS","Ground->CAS"
                  ,"Transpiration","Storage")
         ,col=c("royalblue3","green2","darkgoldenrod","sienna4","chartreuse4","purple4"))
   #---------------------------------------------------------------------------------------#
dev.off()
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#





#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#     Plot the light response curve.                                                       #
#------------------------------------------------------------------------------------------#
pdf(file=file.path(plot.path,"light.pdf"),onefile=FALSE,pointsize=plot.ptsz
   ,width=plot.width,height=plot.height)


   #----- Plot annotation. ----------------------------------------------------------------#
   my.title      = paste("Stand-level light response curve at ",plot.place,sep="")
   my.xlab       = "PAR    [umol/m2/s]"
   my.ylab       = "GPP    [umol/m2/s]"
   #---------------------------------------------------------------------------------------#


   #----- Grab daytime data, and order them by PAR. ---------------------------------------#
   o           = order(par.down)
   dt.daytime  = daytime [o]
   dt.gpp      = gpp     [o]
   dt.par.down = par.down[o]
   dt.gpp      = dt.gpp     [dt.daytime]
   dt.par.down = dt.par.down[dt.daytime]
   #---------------------------------------------------------------------------------------#



   #----- Find convenient scales for x and y, and nice time stamps. -----------------------#
   xlimit   = c(0.,max(dt.par.down))
   ylimit   = c(0.,max(dt.gpp     ))
   whenplot = pretty.time(when,n=10)
   #---------------------------------------------------------------------------------------#



   #----- Main plot. ----------------------------------------------------------------------#
   par(bg="white")
   plot.new()
   plot.window(xlim=xlimit,ylim=ylimit)
   title(main=my.title,xlab=my.xlab,ylab=my.ylab)
   box()
   axis(side=1)
   axis(side=2)
   grid(col="grey50",lty="dotted")
   #----- Add components. -----------------------------------------------------------------#
   points(x=dt.par.down,y=dt.gpp,col="chartreuse",pch=16,cex=0.8)
   #---------------------------------------------------------------------------------------#


   #----- Fit a light response curve. -----------------------------------------------------#
   daytime.data    = data.frame( gpp = dt.gpp, par.down = dt.par.down)
   fit.light  = nls( formula = gpp ~ a1 + a2 * par.down / (a3 + par.down)
                   , data    = daytime.data
                   , start   = list(a1 = 1., a2 = -40., a3 = 500.)
                   )#end nls
   pred.light = predict(fit.light,data=daytime.data)
   lines(x=daytime.data$par.down,y=pred.light,lwd=3.0,col="darkgreen")
   #---------------------------------------------------------------------------------------#
dev.off()
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
