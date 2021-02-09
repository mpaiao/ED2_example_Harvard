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
#------ First and last month/year to include in the plots (MM/DD/YYYY). -------------------#
montha          = 06    # First month
yeara           = 2000  # First year
monthz          = 12    # Last month
yearz           = 2000  # Last year 
#------------------------------------------------------------------------------------------#


#------------------------------------------------------------------------------------------#
#       Some plot parameters.                                                              #
#------------------------------------------------------------------------------------------#
plot.path   = "/home/pecan/ed2ws.harvard/monthly_plots"
plot.place  = "Harvard Forest"
plot.ptsz   = 14     # Default font size
plot.width  = 9.7    # Window width
plot.height = 6.0    # Window.height
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#      Handy constants.                                                                    #
#------------------------------------------------------------------------------------------#
Watts.2.MEin      = 4.6     #  Watts/m2 -> umol/m2/s (MEin)
t3ple             = 273.16  #  Triple point of water             [     K]
alvl3             = 2.5e6   #  Latent heat of vap @ t3ple        [  J/kg]
cph2o             = 1859.   #  Specific heat of water vapor      [J/kg/K]
cliq              = 4186.   #  Specific heat of liquid water     [J/kg/K]
min.dbh.recruit   = 5.      #  Minimum DBH for demographic rates [    cm]
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
#       Find the number of months.                                                         #
#------------------------------------------------------------------------------------------#
n.when     = (yearz-yeara-1)*12+monthz+(12-montha+1)
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#      Create the output path.                                                             #
#------------------------------------------------------------------------------------------#
if (! file.exists(plot.path)) dir.create(plot.path)
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#      Names and colors for all PFTs, so it works for tropical and temperate.              #
#------------------------------------------------------------------------------------------#
pft.names = c("C4 grass"          ,"Early tropical"    ,"Mid tropical"      
             ,"Late tropical"     ,"Temperate C3 Grass","North Pine"        
             ,"South Pine"        ,"Late conifer"      ,"Early hardwood"    
             ,"Mid hardwood"      ,"Late hardwood"     ,"C3 crop"           
             ,"C3 pasture"        ,"C4 crop"           ,"C4 pasture"        
             ,"C3 grass"          ,"Araucaria"         ,"Total"             )
pft.cols  = c("gold"              ,"chartreuse"        ,"chartreuse4"       
             ,"#004E00"           ,"mediumpurple1"     ,"deepskyblue"       
             ,"mediumturquoise"   ,"royalblue4"        , "darkorange"       
             ,"orangered"         ,"firebrick4"         , "purple4"          
             ,"darkorchid1"       ,"darkgoldenrod"     ,   "khaki"          
             ,"lightgoldenrod3"   ,"steelblue3"        ,   "grey22"         )
n.pft     = length(pft.names) - 1
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#      DBH classes.                                                                        #
#------------------------------------------------------------------------------------------#
dbh.brks   <<- c(-Inf,10,25,40,Inf)
dbh.names  <<- c("< 10","10-25","25-40","> 40","Total")
dbh.cols   <<- c("royalblue3","chartreuse3","yellow3","orangered","grey22")
n.dbh      <<- length(dbh.names) - 1
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#       Set up the arrays we want to plot.                                                 #
#------------------------------------------------------------------------------------------#
when         = rep(NA,times=n.when)
#----- 1. Variables that we are going to plot by PFT. -------------------------------------#
lai.pft      = matrix(0,nrow=n.when,ncol=n.pft+1)
agb.pft      = matrix(0,nrow=n.when,ncol=n.pft+1)
bsa.pft      = matrix(0,nrow=n.when,ncol=n.pft+1)
#----- 2. Variables that we are going to plot by DBH. -------------------------------------#
gpp.dbh      = matrix(0,nrow=n.when,ncol=n.pft+1)
i.gpp.dbh    = matrix(0,nrow=n.when,ncol=n.pft+1)
#----- 3. Individual-based demographic rates. ---------------------------------------------#
growth.dbh   = rep(NA,times=n.when)
mortality    = rep(NA,times=n.when)
recruitment  = rep(NA,times=n.when)
#----- 4. Variables that we are going to use to plot in 3 dimensions. ---------------------#
coh.nplant   = list()
coh.age      = list()
coh.dbh      = list()
coh.pft      = list()
coh.height   = list()
coh.cbal     = list()
#------------------------------------------------------------------------------------------#



#------------------------------------------------------------------------------------------#
#       Loop over time.                                                                    #
#------------------------------------------------------------------------------------------#
w    = 0
for (yy in yeara:yearz){

   if (yy == yeara){
      month.begin = montha
   }else{
      month.begin = 1
   }#end if

   if (yy == yearz){
      month.end = monthz
   }else{
      month.end = 12
   }#end if

   for (mm in month.begin:month.end){
      w = w + 1

      #------ Make the file name. ---------------------------------------------------------#
      year.now  = sprintf("%4.4i",yy )
      month.now = sprintf("%2.2i",mm )
      day.now   = sprintf("%2.2i",0)
      hour.now  = sprintf("%2.2i",0)
      minu.now  = sprintf("%2.2i",0)
      seco.now  = sprintf("%2.2i",0)

      when [w]  = chron(dates=paste(mm,1,yy,sep="/"),times=paste(0,0,0,sep=":"))

      file.now  = paste(prefix.analysis,"-Q-",year.now,"-",month.now,"-",day.now,"-"
                       ,hour.now,minu.now,seco.now,"-",suffix.analysis,sep="")
      cat(" - Reading file :",file.now,"...","\n")
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #       Load the data.                                                               #
      #------------------------------------------------------------------------------------#
      now = hdf5load(file=file.now,load=FALSE,verbosity=0,tidy=TRUE)
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #        Define the global number of patches and cohorts.                            #
      #------------------------------------------------------------------------------------#
      npatches        = now$SIPA.N
      ncohorts        = now$PACO.N
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #      In case there is no cohort (i.e. desert), skip the data loading.              #
      #------------------------------------------------------------------------------------#
      if (sum(ncohorts) != 0){

         #---------------------------------------------------------------------------------#
         #     Extend the area and age of each patch so it has the same length as the      #
         # cohorts.                                                                        #
         #---------------------------------------------------------------------------------#
         area       = now$AREA * rep(now$AREA.SI,times=npatches)
         area       = rep(area    ,times=ncohorts)
         age        = rep(now$AGE ,times=ncohorts)
         #---------------------------------------------------------------------------------#



         #---------------------------------------------------------------------------------#
         #---------------------------------------------------------------------------------#
         #      Grab other cohort level variables.                                         #
         #---------------------------------------------------------------------------------#
         ipft       = now$PFT
         dbh        = now$DBH
         idbh       = as.integer(cut(x=now$DBH,dbh.brks))
         nplant     = now$NPLANT       * area
         lai        = now$LAI.CO       * area
         height     = now$HITE
         gpp        = now$MMEAN.GPP.CO
         agb        = now$AGB.CO
         bsa        = now$BA.CO
         cbal       = ( ( now$MMEAN.GPP.CO
                        - now$MMEAN.LEAF.RESP.CO     - now$MMEAN.ROOT.RESP.CO
                        - now$MMEAN.GROWTH.RESP.CO   - now$MMEAN.STORAGE.RESP.CO 
                        - now$MMEAN.VLEAF.RESP.CO    - now$MMEAN.LEAF.MAINTENANCE
                        - now$MMEAN.ROOT.MAINTENANCE - now$MMEAN.LEAF.DROP.CO     ) 
                      / now$AGB.CO )
         #---------------------------------------------------------------------------------#
         #---------------------------------------------------------------------------------#





         #---------------------------------------------------------------------------------#
         #---------------------------------------------------------------------------------#
         #     Aggregate the data by PFT.                                                  #
         #---------------------------------------------------------------------------------#
         #----- 1. LAI, sum is fine. ------------------------------------------------------#
         lai.pft.now        = tapply(X=lai,INDEX=ipft,FUN=sum)
         idx                = as.numeric(names(lai.pft.now))
         lai.pft[w,idx   ]  = lai.pft.now
         lai.pft[w,n.pft+1] = sum(lai)
         #----- 2. AGB, must scale by demographic density. --------------------------------#
         agb.pft.now        = tapply(X=agb*nplant,INDEX=ipft,FUN=sum)
         idx                = as.numeric(names(agb.pft.now))
         agb.pft[w,idx   ]  = agb.pft.now
         agb.pft[w,n.pft+1] = sum(agb*nplant)
         #----- 3. Basal area, must scale by demographic density. -------------------------#
         bsa.pft.now        = tapply(X=bsa*nplant,INDEX=ipft,FUN=sum)
         idx                = as.numeric(names(bsa.pft.now))
         bsa.pft[w,idx   ]  = bsa.pft.now
         bsa.pft[w,n.pft+1] = sum(bsa*nplant)
         #---------------------------------------------------------------------------------#
         #---------------------------------------------------------------------------------#





         #---------------------------------------------------------------------------------#
         #---------------------------------------------------------------------------------#
         #     Aggregate the data by DBH.                                                  #
         #---------------------------------------------------------------------------------#
         #----- 1. GPP in kgC/m2/yr, we must correct by demographic density. --------------#
         gpp.dbh.now          = tapply(X=gpp*nplant,INDEX=idbh,FUN=sum)
         idx                  = as.numeric(names(gpp.dbh.now))
         gpp.dbh[w,idx   ]    = gpp.dbh.now
         gpp.dbh[w,n.dbh+1]   = sum(gpp*nplant)
         #----- 2. GPP in kgC/plant/yr, we must do a weighted average. --------------------#
         tmp                  = data.frame(x=gpp,w=nplant)
         tmp                  = split(tmp,f=idbh)
         i.gpp.dbh.now        = sapply(X=tmp,FUN=weighted.mean)
         idx                  = as.numeric(names(i.gpp.dbh.now))
         i.gpp.dbh[w,idx]     = i.gpp.dbh.now
         i.gpp.dbh[w,n.dbh+1] = weighted.mean(x=gpp,w=nplant)
         #---------------------------------------------------------------------------------#
         #---------------------------------------------------------------------------------#





         #---------------------------------------------------------------------------------#
         #---------------------------------------------------------------------------------#
         #     Find the demographic rates.                                                 #
         #---------------------------------------------------------------------------------#
         #----- 0. Select trees that are large enough to be surveyed. ---------------------#
         is.rec = dbh >= min.dbh.recruit
         #---------------------------------------------------------------------------------#
         #  1.  Growth rates (%DBH/month).  We use the exponential form, then find the     #
         #      weighted mean for all cohorts:                                             #
         #                                                                                 #
         #             ln(DBH_now_n) - ln(DBH_previous_n)                                  #
         #      g_n = ------------------------------------                                 #
         #                          Delta t                                                #
         #                                                                                 #
         #---------------------------------------------------------------------------------#
         #----- 1. Growth rates, save results in %dbh/month. ------------------------------#
         dbh.before     = pmin(dbh,dbh - now$DDBH.DT / 12)
         growth.now     = log(dbh/dbh.before)
         growth.dbh [w] = 100. * weighted.mean(x=growth.now[is.rec],w=nplant[is.rec]) / 12.
         #---------------------------------------------------------------------------------#
         #  2.  Mortality rates (%population/month).  We use the exponential form:         #
         #                                                                                 #
         #           ln(N) - ln(S)       N -- individuals alive in the previous time       #
         #      m = ---------------      S -- individuals alive this time                  #
         #              Delta t                                                            #
         #                                                                                 #
         #---------------------------------------------------------------------------------#
         if (sum(ncohorts) == 0){
            mort.coh   = NA
         }else if (sum(ncohorts) == 1){
            mort.coh   = sum(now$MMEAN.MORT.RATE.CO)
         }else{
            mort.coh   = rowSums(now$MMEAN.MORT.RATE.CO)
         }#end if (sum(nohorts) == 0)
         S             = sum(nplant[is.rec])
         N             = sum(nplant[is.rec] * exp(mort.coh[is.rec]/12))
         mortality [w] = 100. * log(N/S)    # mortality in %population/month
         #---------------------------------------------------------------------------------#
         #  3.  Recruitment rates (%population/month).  We use the exponential form:       #
         #                                                                                 #
         #           ln(N) - ln(E)       E -- individuals that were recruits last time     #
         #      r = ---------------      N -- individuals that are recruits this time      #
         #              Delta t                                                            #
         #                                                                                 #
         #---------------------------------------------------------------------------------#
         last.dbh       = dbh * exp(-growth.now / 12)
         is.est         = is.rec & ( last.dbh >= min.dbh.recruit )
         E              = sum(nplant[is.est])
         N              = sum(nplant[is.rec])
         recruitment[w] = 100. * log(N/E)
         #---------------------------------------------------------------------------------#


         #---------------------------------------------------------------------------------#
         #      Copy cohort variables to the lists.                                        #
         #---------------------------------------------------------------------------------#
         coh.nplant [[w]] = nplant
         coh.age    [[w]] = age
         coh.dbh    [[w]] = dbh
         coh.pft    [[w]] = ipft
         coh.height [[w]] = height
         coh.cbal   [[w]] = cbal
         #---------------------------------------------------------------------------------#
         #---------------------------------------------------------------------------------#
      }else{
         #---------------------------------------------------------------------------------#
         #---------------------------------------------------------------------------------#
         #     Assign a single NA to cohort variables so we know it became a desert.       #
         #---------------------------------------------------------------------------------#
         coh.nplant [[w]] = NA
         coh.age    [[w]] = NA
         coh.dbh    [[w]] = NA
         coh.pft    [[w]] = NA
         coh.height [[w]] = NA
         coh.cbal   [[w]] = NA
         #---------------------------------------------------------------------------------#
      }#end if
      #------------------------------------------------------------------------------------#
   }#end for
}#end for
#------------------------------------------------------------------------------------------#



#----- Make sure that when is still a chron object. ---------------------------------------#
when = chron(when)
#------------------------------------------------------------------------------------------#



#----- Keep only the PFTs of interest. ----------------------------------------------------#
pft.use = which(colMeans(agb.pft) > 0.)
#------------------------------------------------------------------------------------------#


#----- Find the time-related annotation. --------------------------------------------------#
xlimit   = range(when)
whenplot = pretty.time(when,n=10)
#------------------------------------------------------------------------------------------#



#==========================================================================================#
#==========================================================================================#
#     Plot the time series by PFT.                                                         #
#------------------------------------------------------------------------------------------#
pdf(file=file.path(plot.path,"pft-plot.pdf"),onefile=FALSE,pointsize=plot.ptsz
   ,width=plot.width,height=plot.height)

   #----- Plot annotation. ----------------------------------------------------------------#
   my.title = plot.place
   #---------------------------------------------------------------------------------------#



   #----- Split the window into 4. --------------------------------------------------------#
   par   (oma=c(0,0,2,0),bg="white")
   layout(mat=rbind(c(1,2),c(3,4)))
   #---------------------------------------------------------------------------------------#



   #----- First plot: Above-ground biomass by PFT. ----------------------------------------#
   par(mar=c(3.1,4.1,1.1,1.1))
   plot.new()
   plot.window(xlim=xlimit,ylim=c(0,max(agb.pft[,pft.use])))
   title(ylab="Above-ground biomass [kgC/m2]")
   axis(side=1,at=whenplot$levels,labels=whenplot$labels,padj=whenplot$padj)
   axis(side=2)
   box()
   abline(v=whenplot$levels,h=axTicks(side=2),col="grey50",lty="dotted")
   for (p in pft.use) lines(x=when,y=agb.pft[,p],lwd=2.5,col=pft.cols[p])
   #---------------------------------------------------------------------------------------#



   #----- Second plot: Basal area by PFT. -------------------------------------------------#
   par(mar=c(3.1,4.1,1.1,1.1))
   plot.new()
   plot.window(xlim=xlimit,ylim=c(0,max(bsa.pft[,pft.use])))
   title(ylab="Basal area [cm2/m2]")
   axis(side=1,at=whenplot$levels,labels=whenplot$labels,padj=whenplot$padj)
   axis(side=2)
   box()
   abline(v=whenplot$levels,h=axTicks(side=2),col="grey50",lty="dotted")
   for (p in pft.use) lines(x=when,y=bsa.pft[,p],lwd=2.5,col=pft.cols[p])
   #---------------------------------------------------------------------------------------#



   #----- Third plot: Leaf area index by PFT. ---------------------------------------------#
   par(mar=c(3.1,4.1,1.1,1.1))
   plot.new()
   plot.window(xlim=xlimit,ylim=c(0,max(lai.pft[,pft.use])))
   title(ylab="Leaf area index [m2/m2]")
   axis(side=1,at=whenplot$levels,labels=whenplot$labels,padj=whenplot$padj)
   axis(side=2)
   box()
   abline(v=whenplot$levels,h=axTicks(side=2),col="grey50",lty="dotted")
   for (p in pft.use) lines(x=when,y=lai.pft[,p],lwd=2.5,col=pft.cols[p])
   #---------------------------------------------------------------------------------------#



   #----- Fourth plot: the legend. --------------------------------------------------------#
   par(mar=c(0.1,0.1,0.1,0.1))
   plot.new()
   plot.window(xlim=c(0,1),ylim=c(0,1))
   legend(x="center",inset=0.0,legend=pft.names[pft.use],col=pft.cols[pft.use],lwd=2.5
         ,ncol=2,title=expression(bold("Plant functional types")))
   #---------------------------------------------------------------------------------------#


   #----- Title. --------------------------------------------------------------------------#
   mtext(side=3,text=paste("Time series by PFT - ",plot.place,sep=""),outer=TRUE
        ,font=2,cex=1.2)
   #---------------------------------------------------------------------------------------#

dev.off()
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#     Plot GPP by DBH class.                                                               #
#------------------------------------------------------------------------------------------#
pdf(file=file.path(plot.path,"gpp-by-dbh.pdf"),onefile=FALSE,pointsize=plot.ptsz
   ,width=plot.width,height=plot.height)


   #----- Plot annotation. ----------------------------------------------------------------#
   my.title = plot.place
   #---------------------------------------------------------------------------------------#



   #----- Split the window into 4. --------------------------------------------------------#
   par   (oma=c(0,0,2,0),bg="white")
   layout(mat=rbind(c(1,2),c(3,3)),heights=c(9,2))
   #---------------------------------------------------------------------------------------#



   #----- First plot: GPP [kgC/m2/yr]. ----------------------------------------------------#
   par(mar=c(3.1,4.1,1.1,1.1))
   plot.new()
   plot.window(xlim=xlimit,ylim=c(0,max(gpp.dbh)))
   title(ylab="Gross Primary Productivity [kgC/m2/yr]")
   axis(side=1,at=whenplot$levels,labels=whenplot$labels,padj=whenplot$padj)
   axis(side=2)
   box()
   abline(v=whenplot$levels,h=axTicks(side=2),col="grey50",lty="dotted")
   for (d in sequence(n.dbh+1)) lines(x=when,y=gpp.dbh[,d],lwd=2.5,col=dbh.cols[d])
   #---------------------------------------------------------------------------------------#



   #----- First plot: GPP [kgC/plant/yr]. -------------------------------------------------#
   par(mar=c(3.1,4.1,1.1,1.1))
   plot.new()
   plot.window(xlim=xlimit,ylim=c(0,max(i.gpp.dbh)))
   title(ylab="Gross Primary Productivity [kgC/plant/yr]")
   axis(side=1,at=whenplot$levels,labels=whenplot$labels,padj=whenplot$padj)
   axis(side=2)
   box()
   abline(v=whenplot$levels,h=axTicks(side=2),col="grey50",lty="dotted")
   for (d in sequence(n.dbh+1)) lines(x=when,y=i.gpp.dbh[,d],lwd=2.5,col=dbh.cols[d])
   #---------------------------------------------------------------------------------------#



   #----- Third plot: the legend. ---------------------------------------------------------#
   par(mar=c(0.1,0.1,0.1,0.1))
   plot.new()
   plot.window(xlim=c(0,1),ylim=c(0,1))
   legend(x="bottom",inset=0.0,legend=dbh.names,col=dbh.cols,lwd=2.5,ncol=3
         ,title=expression(bold("DBH classes [cm]")))
   #---------------------------------------------------------------------------------------#


   #----- Title. --------------------------------------------------------------------------#
   mtext(side=3,text=paste("Time series by DBH - ",plot.place,sep=""),outer=TRUE
        ,font=2,cex=1.1)
   #---------------------------------------------------------------------------------------#
dev.off()
#==========================================================================================#
#==========================================================================================#





#==========================================================================================#
#==========================================================================================#
#     Plot demographic rates.                                                              #
#------------------------------------------------------------------------------------------#
pdf(file=file.path(plot.path,"demography.pdf"),onefile=FALSE,pointsize=plot.ptsz
   ,width=plot.width,height=plot.height)


   #----- Split the window into 4. --------------------------------------------------------#
   par   (bg="white")
   layout(mat=rbind(1,2),heights=c(5,1))
   #---------------------------------------------------------------------------------------#



   #----- Plot the three demographic rates together. --------------------------------------#
   par(mar=c(3.1,4.1,4.1,2.1))
   plot.new()
   plot.window(xlim=xlimit,ylim=c(0,max(c(mortality,recruitment,growth.dbh))))
   title(main=paste("Demographic rates - ",plot.place),ylab="Rates [%/month]")
   axis(side=1,at=whenplot$levels,labels=whenplot$labels,padj=whenplot$padj)
   axis(side=2)
   box()
   abline(v=whenplot$levels,h=axTicks(side=2),col="grey50",lty="dotted")
   lines(x=when,y=recruitment,lwd=2.5,col="chartreuse3" )
   lines(x=when,y=mortality  ,lwd=2.5,col="purple3"     )
   lines(x=when,y=growth.dbh ,lwd=2.5,col="deepskyblue3")
   #---------------------------------------------------------------------------------------#


   #----- Plot the legend. ----------------------------------------------------------------#
   par(mar=c(1.1,0.1,0.1,0.1))
   plot.new()
   plot.window(xlim=c(0,1),ylim=c(0,1))
   legend(x="bottom",inset=0.0,legend=c("Recruitment","Mortality","DBH growth")
         ,col=c("chartreuse3","purple3","deepskyblue3"),lwd=2.5,ncol=3)
   #---------------------------------------------------------------------------------------#
dev.off()
#==========================================================================================#
#==========================================================================================#




#==========================================================================================#
#==========================================================================================#
#     Create the size-and-age structure of NPP.                                            #
#------------------------------------------------------------------------------------------#
#----- Define fixed ranges for the grid. --------------------------------------------------#
coh.var  = lapply(X=coh.cbal,FUN="*",100)
var.pref = "cbal"
var.desc = "Relative carbon balance"
var.unit = "%AGB/yr"
popmin   = min(unlist(coh.nplant),na.rm=TRUE)
popmax   = max(unlist(coh.nplant),na.rm=TRUE)
#------------------------------------------------------------------------------------------#


#----- Define some properties for perspective plot. ---------------------------------------#
theta           = 315.                    # Azimuth for perspective projection
phi             = 30.                     # Vertical angle for perspective projection
ltheta          = -210.                   # Azimuth angle for light
shade           = 0.125                   # Shade intensity
expz            = 0.5                     # Expansion factor for Z axis
gcol            = c("lightblue","white")  # Colours for the 50's style floor
cexmin          = 0.5                     # Minimum "head" size of the lollipop
cexmax          = 2.0                     # Maximum "head" size of the lollipop
#------------------------------------------------------------------------------------------#


#----- Define the grid information for the 3-D plot. --------------------------------------#
ageaxis   = pretty(unlist(coh.age),n=20)
dbhaxis   = pretty(unlist(coh.dbh),n=20)
xlimit    = range(ageaxis)
ylimit    = range(dbhaxis)
hlimit    = range(unlist(coh.height), na.rm=TRUE)
vlimit    = range(unlist(coh.var)   , na.rm=TRUE)
#------------------------------------------------------------------------------------------#


#----- Define the floor location. ---------------------------------------------------------#
if ((vlimit[1] > 0) != (vlimit[2] > 0)){
   floor3d = 0.
}else if (vlimit[1] > 0){
   floor3d = vlimit[1]
}else{
   floor3d = vlimit[2]
}#end if
flooraxis  = matrix(floor3d,nrow=length(ageaxis),ncol=length(dbhaxis))
zeroaxis   = matrix(0.,nrow=length(ageaxis),ncol=length(dbhaxis))
#------------------------------------------------------------------------------------------#

for (w in 1:n.when){

   #----- Find which year we are plotting. ------------------------------------------------#
   cyear    = sprintf("%4.4i",numyears (when[w]))
   cmonth   = sprintf("%2.2i",nummonths(when[w]))
   mon.nice = month.name[nummonths(when[w])]
   cat(" - Plot SAS view of ",var.desc," (",mon.nice,", ",cyear,")","\n")
   #---------------------------------------------------------------------------------------#



   #---------------------------------------------------------------------------------------#
   #     We only plot the SAS figures when the polygon is not an absolute desert.          #
   #---------------------------------------------------------------------------------------#
   if (any (! is.na(coh.var[[w]]))){

      #----- Expand the lines to make the "trees". ----------------------------------------#
      n.coh.now  = length(coh.var[[w]])
      ageww      = rep(coh.age[[w]],each=3)
      dbhww      = rep(coh.dbh[[w]],each=3)
      pftww      = rep(coh.pft[[w]],each=3)
      hgtww      = as.vector(rbind(rep(0.,times=n.coh.now)
                                  ,coh.height[[w]]
                                  ,rep(NA,times=n.coh.now)))
      varww      = as.vector(rbind(rep(floor3d,times=n.coh.now)
                                  ,coh.var[[w]]
                                  ,rep(NA,times=n.coh.now)))
      pchww      = rep(c(NA,16,NA),times=n.coh.now)
      #------------------------------------------------------------------------------------#


      #----- Find the scale for the heads (bigger for more abundant cohorts). -------------#
      cexww    = ( cexmin + (cexmax - cexmin) 
                          * log(coh.nplant[[w]]/popmin) / log(popmax/popmin))
      cexww    = rep(cexww,each=3)
      #------------------------------------------------------------------------------------#


      #----- Find the colours of the PFTs. ------------------------------------------------#
      colww    = pft.cols[pftww]
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     Find the PFTs that are still available at this time.                           #
      #------------------------------------------------------------------------------------#
      pftin   = sort(unique(coh.pft[[w]]))
      colleg  = pft.cols  [pftin]
      pftleg  = pft.names [pftin]
      #------------------------------------------------------------------------------------#




      #----- Open the graphics. -----------------------------------------------------------#
      out.file = file.path(plot.path
                          ,paste("sas_",var.pref,"-",cyear,"-",cmonth,".pdf",sep=""))
      pdf(file=out.file,onefile=FALSE,pointsize=plot.ptsz
         ,width=plot.width,height=plot.height)
      #------------------------------------------------------------------------------------#




      #------------------------------------------------------------------------------------#
      #       Shared annotation.                                                           #
      #------------------------------------------------------------------------------------#
      mytitle = paste(plot.place," - ",mon.nice,"/",cyear,sep="")
      myvlab  = paste(var.desc," [",var.unit,"]",sep="")
      myhlab  = paste("Height [m]")
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #      Split the window into 3.                                                      #
      #------------------------------------------------------------------------------------#
      par(oma=c(0,0,2,0),bg="white")
      layout(mat=rbind(c(1,2),c(3,3)),heights=c(9,2))
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #      First plot: Height SAS.                                                       #
      #------------------------------------------------------------------------------------#
      par(mar=c(3.1,4.1,1.1,1.1))
      #----- The box. ---------------------------------------------------------------------#
      pout = persp(x=ageaxis,y=dbhaxis,z=zeroaxis,xlim=xlimit,ylim=ylimit
                  ,zlim=hlimit,theta=theta,phi=phi,col=gcol,expand=expz
                  ,ticktype="detailed",border=NA,xlab="Gap age [yr]"
                  ,ylab="DBH [cm]",zlab=myhlab,shade=shade,ltheta=ltheta)
      #----- Second plot, the actual data. ------------------------------------------------#
      lines (trans3d(x=ageww,y=dbhww,z=hgtww,pmat=pout),type="l",col="grey29",lwd=2)
      points(trans3d(x=ageww,y=dbhww,z=hgtww,pmat=pout),type="p",pch=pchww
                    ,col=colww,cex=cexww)
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #      Second plot: The variable.                                                    #
      #------------------------------------------------------------------------------------#
      par(mar=c(3.1,4.1,1.1,1.1))
      #----- The box. ---------------------------------------------------------------------#
      pout = persp(x=ageaxis,y=dbhaxis,z=flooraxis,xlim=xlimit,ylim=ylimit
                  ,zlim=vlimit,theta=theta,phi=phi,col=gcol,expand=expz
                  ,ticktype="detailed",border=NA,xlab="Gap age [yr]"
                  ,ylab="DBH [cm]",zlab=myvlab,shade=shade,ltheta=ltheta)
      #----- Second plot, the actual data. ------------------------------------------------#
      lines (trans3d(x=ageww,y=dbhww,z=varww,pmat=pout),type="l",col="grey29",lwd=2)
      points(trans3d(x=ageww,y=dbhww,z=varww,pmat=pout),type="p",pch=pchww
                    ,col=colww,cex=cexww)
      #------------------------------------------------------------------------------------#



      #------------------------------------------------------------------------------------#
      #     The legend.                                                                    #
      #------------------------------------------------------------------------------------#
      par(mar=c(0.1,0.1,0.1,0.1))
      plot.new()
      plot.window(xlim=c(0,1),ylim=c(0,1))
      legend(x="center",inset=0.0,legend=pftleg,fill=colleg,ncol=3
            ,title=expression(bold("Plant functional types")))
      #------------------------------------------------------------------------------------#


      #----- Title. -----------------------------------------------------------------------#
      mtext(side=3,text=mytitle,outer=TRUE,font=2,cex=1.1)
      #------------------------------------------------------------------------------------#


      dev.off()
   }#end if is.na(varww)
}#end for nameco
#==========================================================================================#
#==========================================================================================#


