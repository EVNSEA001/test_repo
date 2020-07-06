########################################################
#####                                              #####
#####             Hunting time method              #####
#####  Optimised automated Broken stick algorithm  #####
#####                                              #####
########################################################

## Karine Heerah - February 2014
## A new method to quantify within dive foraging behaviour in marine predators - K. Heerah, M. Hindell, C. Guinet, J-B, charrassin
## karine.heerah@hotmail.fr
## LOCEAN - UMR 7159, CNRS/UPMC/IRD/MNHN, 4 place Jussieu 75252 Paris Cedex 05, France.
 
rm(list=ls())
 
## To setwd() press cntr + Shift + H

#data_path = "C:/Users/Sean/Documents/sean/MS/ Marion project/Masters project/Computing/Dive data for Sean Evans" # path where your TDR files are stored (after ZOC)
#fig_path="C:/Users/Sean/Documents/sean/MSc Marion project/Masters project/Computing/Dive data for Sean Evans/Figs/4_VDB_2009W_0890439_FB515" # path for figures to be stored
getwd()
setwd(data_path)






################ Fur seals #####################
dt <- read.table("4_VDB_2009W_0890439_FB515.tab",header=TRUE) ## load TDR files sep = '\t',header = F, skip = 2,quote='', comment=''
#head(dt)

############### Fur seals ##################

## Renaming columns and changing time format
names(dt)[names(dt)=="Depth"] <- "depth" #rename the depth column
seal <- "5" #insert seal unique id number here for when you merge files
dt$seal <- seal #add that as a column to your file
dt$daytime <- ISOdatetime(1970,1,1,0,0,0, tz="GMT")+dt$Time #turn that ISOdate number into something you can understand
#range(dt$daytime)

#class(dt$daytime)

#dplot <- dt[,c(2,7)]      
#head(dplot)
# #Plot the dive depths
#par(mfrow=c(2,1))
#theYlim <- c(max(dplot$depth,na.rm=TRUE),-10) ##set the y axis scale
#with(dplot, plot(depth, ylim=theYlim,pch=19, type="h"))# cex=.5, pch=19, type="h"))
#with(dplot, hist(depth[dplot$depth>2]))


#write.table(fn,"27_VDB_2013S_10L0056_FB535.tab", sep = "\t")



####################################################################
##########Part 2: Analysing fur seal dives.##############

install.packages("maps")
install.packages("mapdata")
install.packages("fields")
install.packages("lattice")
install.packages("effects")
install.packages("trip")
install.packages("Matrix")

require(maps)
require(mapdata)
require(fields)
require(lattice)
require(effects)
require(trip)
require(Matrix)


dives <- dt[,c("daytime", "depth", "seal")]
#head(dives)
str(dives)

samp <-1 ##sets the sampling interval 1 second here

##indentify the first and last records
(sta <- dives$daytime[1])
(en  <- dives$daytime[nrow(dives)])

# #Plot the dive depths
# par(mfrow=c(2,1))
# theYlim <- c(max(dives$depth,na.rm=TRUE),-10) ##set the y axis scale
# with(dives, plot(gmt,depth, ylim=theYlim, cex=.5, pch=19, type="h"))
# with(dives, hist(depth[dives$depth>2]))


##Step 1
##The first step is to correct for any drift in the depth sensor over time - called the 
## zero-off set correction (ZOC)
## This is function to do that - you need to run it first.... 
zoc <- function(x, zval = 10, blim = -5){
  h <- hist(x[x < zval & x > blim], breaks = seq(blim, zval, 1))
  # print(h$breaks[which.max(h$counts)])
  
  # print(length(x - h$breaks[which.max(h$counts)]))
  #scan("", 1)
  x - h$breaks[which.max(h$counts)]
}

## ..and then use it to Zoc your data file

ct <- cut(dives$daytime, "2 hours") ## breaks the data in to 2 hour segments (This is not at 1 sec intervals?????????) ???

zd <- unlist(tapply(dives$depth, ct, zoc)) ## runs the zoc function on each 2 h segment

#sets surface intervals within the "surface noise" value to 0
sn <- 1
names(zd) <- NULL #Names of ZD is the timestamp...
nzd <- zd # So this is just then the zoc'ed depth
nzd[zd < sn] <- 0

##......and adds these data to the dataframe
(l <- length(nzd)-length(dives$depth))
if (l==0) dives$cor.depth <- nzd
if (l!=0) dives$cor.depth <- c(nzd,0)

#head(dives)
##assigns numbers to each dive based on a user defined "surface noise" (sn) value 
##Gives each surface or dive section a consecutive number

###############################################
#diff(nzd == 5)

num <- ceiling(cumsum(abs(c(1, diff(nzd == 6))))/2) #ceiling method is used to return the smallest integer value which is greater than or equal to an individual number
dives$num <- num

##re-write dives to include only the corrected depth data
dives <- dives[,c("daytime", "cor.depth", "num", "seal")]
colnames(dives) <- c("daytime", "depth", "num", "seal") #rename the variables
#head(dives)
## calculate some basic statistics: 
#(i) what is the maximum depth reached by your seal? 
#(ii) how many dives are there?
#(iii) what is the maximum depth for each dive

#max(dives$depth)
#max(dives$num)
max.d <- aggregate(depth~num, data=dives, max)
#max.d[1:50,] ##look at the first 50 dives
##add to dives
dives <- merge(dives, max.d, by.x="num", by.y="num", all.x=T)
colnames(dives) <- c("num", "daytime", "depth", "seal", "max.dep") #rename the variables
#head(dives)

##extract dives less than 5m maximum depth
# identify dives with max.dep > 5m 
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Should I only use dives >2m or perhaps >6m...
# Want to compare with binned data so sticking with 5m is best I think.
d.list <- which(dives$max.dep > 10)
# #retain only these records
dives_10m <- dives[d.list,]
#dives_5m[1:50,] ##look at the first 50 dives

#How many dives are we keeping from our analysis?
length(unique(dives_10m$num))
# length(unique(dives2$num))



#####################################################################
################# Heerah Broken stick model #########################

                                         ## POSIXct format
## column names: date-hour = "daytime", id = "seal", depth = "depth", dive number = "num", temperature = "temp", 
## ambient light = "light" etc. 
## The only variables needed for the algorithm are: daytime, depth and dive number. All the others are optionnal

#dt <- dt[,c(1,2,5)] ## We only keep "daytime", "depth" and "num" columns.

#-----------------------------------------------------------------------------------------------------------------------------------------

##### Optimised Broken Stick Algorithm

# 1. Creation of output dataframes

d_env <- dives_10m ## put your dataframe in the d_env variable which is used troughout the script

num = unique(d_env$num) ## number of dives you have
num.list <- num

dbs <- data.frame("num"=rep(0,1), "all.dur"=0,"start"=0,"end"=0,"depth_start"=0,"depth_end"=0,"seg"=0,"npoints"=0,
      "dur"=0,"dur.per"=0,"coef"=0, "mean_depth"=0,"max.depth"=0,"sinuosity"=0,"mean_err"=0,"foraging"=0) ## Broken stick dataframe

## num = dive number, all.dur = total dive duration, start = date of segment start, end = date of segment start, 
## depth_start = depth of segment start, seg = broken stick segment number, npoints = number of points summarising the dive
## dur = duration of each segment, dur.per = % of total duration, coef = slope coefficient of the segment, 
## mean_depth = mean depth of segment, max.depth = maximum dive depth, sinuosity = vertical sinuosity associated to this part of the dive,
## mean_err = mean distance between original dive profile and the reconstructed one for the optimal number of broken stick points,
## foraging = behaviour according to vertical sinuosity threshold

ncdv = data.frame("daytime"=rep(0,1) ,"depth"=0,"num"=0) ## dataframe in which the dives for which the fit doesn't work will be stored.
                                                         ## Needs to have a same column names exactly
d_env <- d_env[,c("daytime","depth", "num")]
#-----------------------------------------------------------------------------------------------------------------------------------------
# loop for each dive

for(d in 1:length(num.list)){
  print(d)
  dt <- d_env[d_env$num==num.list[d],]
  if(nrow(dt) > 60 ) {  ##consider dives of more than 60s as the resolution of the dataset is 1s
    ndive=num.list[d]
    
    #plot(as.numeric(dt$daytime), dt$depth, ylim=c(max(dt$depth),0), t="l", ylab="depth (m)", xlab="",xaxt="n") 
    #Use only if you want to check your dive
    #axis.POSIXct(1,x=dt$daytime, format="%H:%M:%S", labels = TRUE,cex.lab=0.5)
    #idem
    
    np <- c(3:30)  ## number of broken stick iterations to see which optimal number of points summarise your dive
    npe=rep(NA,28) ## vector where the average distance between original and reconstructed dive profile is stored
    npo=rep(NA,28) ## vector where the number of broken stick points describing the dive profile is stored
    

### Finding the optimal number of Broken points for each dive

# 1.Loop to define the mean distance depending on the number of broken stick points 
    
  for (k in 1:length(np)){
      npp = np[k] # selection of the number of iteration: from 3 to 30
      # 2 lines below: selection of the depth and time for the 2 surface points and the maximum depth point
      ref <- c(dt$depth[1],max(dt$depth),dt$depth[nrow(dt)])
      tim <-c(as.numeric(dt$daytime[1]),as.numeric(dt$daytime[dt$depth==max(dt$depth)][1]),as.numeric(dt$daytime[nrow(dt)]))
      
      for (i in 1:npp){
        #plot(as.numeric(dt$daytime), dt$depth, ylim=c(max(dt$depth),0), t="l", ylab="depth (m)", xlab="",xaxt="n")
        # plot only if you want to see how the broken stick algorithm is working
        #points(tim,ref, pch=19, cex=1, col="red")
        #idem
        interp <- approx(tim,ref,xout=dt$daytime,method="linear") #linear interpolation between broken stick points at TDR time interval
        #lines(interp,col="red")
        #idem
        dif_x <- as.numeric(interp$x - dt$daytime) # time differences between original and reconstructed profiles
        dif_y <- interp$y - dt$depth # depth differences between original and reconstructed profiles
        dst <- sqrt(dif_x^2 + dif_y^2) # calculate distances between original and reconstructed profiles
        
        ii <- which(dst==max(dst))[1] # index of the data point of maximum difference between original and reconstructed profiles
        #points(dt$daytime[ii],dt$depth[ii],col="blue",pch=19,cex=1)
        #idem
        tim <- c(as.numeric(tim),as.numeric(dt$daytime[ii])) # add new broken stick point time
        tim <- ISOdatetime(1970,1,1, 0,0,0, tz="gmt") + tim 
        ref <- c(ref,dt$depth[ii]) # add new broken stick point depth
      }
      npe[k] = mean(dst) # average distance between original and reconstructed dive profiles
      npo[k] = length(tim) # number of broken stick points describing the dive profile
    }
    

# 2. Defining the optimal number of broken stick points
    
    f <- data.frame(npe=npe, npo=npo)
    #plot(f$npo, f$npe,xlab="nb of points", ylab="mean error") #plot of mean distance between original and reconstructed dive profiles
                                                               # according to the number of broken stick points describing the dive
                                                               #activate only if you want to check
    
    # Use of a gompertz model to find the curve which best fit our data
    Asym <- 0; b2 <- -5; b3 <- 0.9
    fm1 <- -999
    try(fm1 <- nls(npe ~ SSgompertz(npo, Asym, b2, b3), data=f, control=nls.control(maxiter=500)),TRUE) #gompertz model to fit an asymptote
    #curve to the mean distance between original and reconstructed dive profiles plot
    
    if (class(fm1) == "nls"){ # if the model converged, we can go to the next steps
      #summary(fm1)
      tt <-predict(fm1, f$npe)
      
      
      # plot of the mean distance 
      png(paste(fig_path,"WED_BS_",ndive,"_",substr(dt$daytime[1],1,10),".png", sep=""),pointsize=12*1.5,height=480*1.5,width=480*1.5)
      
      par(mfrow=c(2,1),mar=c(4,4,2,2))
      tit=paste("BS_WED08_",ndive,"_",substr(dt$daytime[1],1,10))
      plot(f$npo, f$npe,xlab="nb of points", ylab="mean error",main=tit)
      lines(na.omit(f$npo),tt[1:28],col="red")
      
      # Plot the linear approximation between the first and last point of the fitted curve
      t <- data.frame(npe=c(f$npe[1], f$npe[28]), npo=c(f$npo[1], f$npo[28]))
      interp <- approx(c(f$npo[1],f$npo[28]),c(tt[1],tt[28]), xout=f$npo,method="linear")
      interp$x <- interp$x[!is.na(interp$x)]
      interp$y <- interp$y[!is.na(interp$y)]
      lines(interp$x, interp$y,col="blue")
      
      # Looking for the inflexion point which is the furthest point between the fitted curve and the approximation 
      dif_x <- interp$x - na.omit(f$npo)
      dif_y <- interp$y - tt[1:28]
      dst <- sqrt(dif_x^2 + dif_y^2)
      dm <- f$npo[which(dst==max(dst))]
      
      points(f$npo[which(dst==max(dst))], f$npe[which(dst==max(dst))], pch=19, col="red") ## inflexion point
      
      
#3. optimal broken stick method for each dive 
      
      # The two lines below select the optimal number of broken stick points (in their order of appearance in the BS iteration)
      # example: surface start point, max. depth point, surface end point + x other points
      tim= tim[1:dm] 
      ref=ref[1:dm]
      
      tim2 <- sort(tim)
      dep_tim <- as.data.frame(cbind(ref,tim))
      dep_tim <- dep_tim[order(tim),]
      
      dbs2 <- data.frame("num"=rep(0,(nrow(dep_tim)-1)) , "all.dur"=0,"start"=0,"end"=0,"depth_start"=0,"depth_end"=0,"seg"=0,"npoints"=0, "dur"=0,"dur.per"=0,"coef"=0, "mean_depth"=0, "max.depth"=0,"sinuosity"=0,"mean_err"=0)
      
      # Loop to calculate the different metrics for each broken stick segments
      for (n in 1:(nrow(dep_tim)-1)){
        x1= dep_tim$tim[n] # start of BS segment
        x2= dep_tim$tim[n+1] #end of BS segment
        dbs2$num[n]=ndive
        dbs2$all.dur[n]=difftime(dt$daytime[nrow(dt)], dt$daytime[1], tz,units = c("secs")) #dive duration
        dbs2$start[n]=x1 
        dbs2$end[n]=x2
        dbs2$depth_start[n]= dep_tim$ref[n] # depth of start of BS segment
        dbs2$depth_end[n]= dep_tim$ref[n+1] # depth of end of BS segment
        dbs2$seg[n]=n #segment number
        dbs2$npoints[n]=nrow(dep_tim) # optimal BS points summarising the original dive profile
        dbs2$dur[n]= difftime(tim2[n+1], tim2[n], tz,units = c("secs")) #duration of the segment in sec.
        dbs2$dur.per[n]=(dbs2$dur[n]/dbs2$all.dur[n])*100 #% of segment duration according to total dive duration
        dbs2$coef[n]=(dep_tim$ref[n+1] - dep_tim$ref[n])/(x2 - x1) # slope coefficient of the segment
        dbs2$mean_depth[n]=mean(dt$depth[which(as.numeric(dt$daytime)==x1):which(as.numeric(dt$daytime)==x2)]) #mean depth of the segment 
        # calculated from original profile depths
        dbs2$max.depth[n]= max(dt$depth) # dive max. depth
        
        #Calculation of vertical sinuosity
        deuc= abs(dep_tim$ref[n+1] - dep_tim$ref[n]) # Vertical distance swum between 2 BS points
        dobs=sum(abs(diff(dt$depth[which(dt$daytime==x1):which(dt$daytime==x2)]))) # sum of all the vertical distances from the original
        #profile between the two corresponding BS depth points
        
        dbs2$sinuosity[n]=deuc/dobs # vertical sinuosity index
        dbs2$mean_err[n]=f$npe[which(dst==max(dst))] # mean distance between original and reconstructed dive profiles for the optimal
        #number of BS points summarising the dive.
      }
      
      #-----------------------------------------------------------------------------------------------------------------------------------
      # IMPORTANT:
      #-----------
      # Attribution of behaviour according to vertical sinuosity -- Remind that the sinuosity threshold used here was determined according
      # to the histogram/density plot of vertical sinuosity for every BS segments of every dive
      # so, before setting your threshold at 0.9, check if it suits your dataset (i.e after running the BS on all your dives)
      
      dbs2$foraging <- 2 ## 2 stands for "hunting" mode
      dbs2$foraging[dbs2$sinuosity >=0.9 & dbs2$sinuosity <=1] <-1 ## 1 stands for "transit" mode
      #-----------------------------------------------------------------------------------------------------------------------------------      

      # Dive plot: original dive profile and Broken stick reconstructed profile
      sg <- unique(dbs2$seg)
      cl <- c("blue","red")
      dbs2$code[dbs2$sinuosity]

      plot(as.numeric(dt$daytime), dt$depth, ylim=c(max(dt$depth),0),t="l",ylab="depth (m)", xlab="",xaxt="n")
      points(tim,ref, pch=19, cex=1, col="black")
      lines(approx(tim,ref,xout=dt$daytime,method="linear"),col="black")
      for(i in 1:length(sg)){
        lines(c(dbs2$start[dbs2$seg==sg[i]],dbs2$end[dbs2$seg==sg[i]]),c(dbs2$depth_start[dbs2$seg==sg[i]],
        dbs2$depth_end[dbs2$seg==sg[i]]),col=cl[dbs2$foraging][dbs2$seg==sg[i]],lwd=2.5)
        }
      axis.POSIXct(1,x=dt$daytime, format="%H:%M:%S",labels = TRUE,cex.lab=0.5)
      dev.off()
      dbs <-rbind(dbs,dbs2)
    } else {ncdv<-rbind(ncdv,dt)} # allows to keep somewhere the data for which the fit of the Gompertz model didn't work ? 
  } #end of if loop for 60 s
  #save(dbs, file="BS_fitmet_WED_08_samp.RData")
} #end of for loop for dive number

dbs<- dbs[-1,]

hist(dbs$sinuosity,xlab="sinuosity",breaks=seq(0,1,0.1),main="") ## See line 191
abline(v=0.9,col="red",lwd=2)

ncdv <- ncdv[-1,]
#save(ncdv, file="BS_err_WED_08_samp.RData")
#save(dbs, file="BS_fitmet_WED_08_samp.RData")
mean_for_depth <- dbs[,c("dur","mean_depth","mean_err","coef","foraging")]
mds.list <- which(mean_for_depth$foraging==2 & mean_for_depth$coef<0.2 
                  & mean_for_depth$dur>20 & mean_for_depth$mean_depth>6)
mds.list <- which(dbs$start >"15:00:00" & dbs$start<"00:00:00")
#mfd.list <- mean_for_depth[,"foraging"==2]
mean_for_depth <- mean_for_depth[mds.list,]
#names(dt)[names(dt)=="Depth"] <- "depth"
hist(mean_for_depth$mean_depth,xlab="mean_depth",breaks=seq(6,80,1),main="") ## See line 191
#abline(v=0.9,col="red",lwd=2)
length(mean_for_depth$mean_depth)
dbs$start <- ISOdatetime(1970,1,1,0,0,0, tz="GMT")+dbs$start
dbs$end <- ISOdatetime(1970,1,1,0,0,0, tz="GMT")+dbs$end



























############### Ellies ##################

rm(list = ls())
fig_path="C:/Users/Sean/Documents/sean/MSc Marion project/Masters project/Computing/MEOP data/Figs/ct99-LB037-13_ODV" # path for figures to be stored

## Ellies
dt <- read.table("ct99-LB037-13_ODV.txt",header=TRUE, sep = "\t", dec = ".")
#head(dt)

names(dt)[names(dt)=="Depth..m."] <- "depth" #rename the depth column
names(dt)[names(dt)=="Temperature...C."] <- "Temperature" #rename temp
#seal <- "13" #insert seal unique id number here for when you merge files
#dt$seal <- seal #add that as a column to your file
#dt$daytime <- ISOdatetime(1970,1,1,0,0,0, tz="GMT")+dt$Time
#class(dt$mon.day.yr)
#as.Date(dt$mon.day.yr)

library("lubridate")
require("lubridate")
parse_date_time(dt$mon.day.yr, orders="%m/%d/%y", tz = "GMT") -> dt$mon.day.yr
#class(dt$hh.mm)

#parse_date_time(dt$hh.mm, orders="%H:%M", tz = "GMT") -> dt$hh.mm

#format(dt$hh.mm, "%H:%M") 
#dt$mon.day.yr <- format(dt$mon.day.yr, "%m/%d/%y") + format(dt$hh.mm, "%H:%M") 

#as.POSIXct(paste(dt$mon.day.yr, dt$hh.mm), format="%Y-%m-%d %H:%M")
dt$daytime <- with(dt, ymd(mon.day.yr) + hm(hh.mm))

dives <- dt[,c(8,14)]



dives$daytime <- as_datetime(dives$daytime, tz = "GMT")

####################################################################
##########Part 2: Analysing Ellies dives.##############


install.packages("maps")
install.packages("mapdata")
install.packages("fields")
install.packages("lattice")
install.packages("effects")
install.packages("trip")
install.packages("Matrix")

require(maps)
require(mapdata)
require(fields)
require(lattice)
require(effects)
require(trip)
require(Matrix)

###############################################


num <- ceiling(cumsum(abs(c(0, diff(dives$depth == 10))))/2) #ceiling method is used to return the smallest integer value which is greater than or equal to an individual number
dives$num <- num

## calculate some basic statistics: 
#(i) what is the maximum depth reached by your seal? 
#(ii) how many dives are there?
#(iii) what is the maximum depth for each dive

max(dives$depth)
max(dives$num)
max.d <- aggregate(depth~num, data=dives, max)
#max.d[1:50,] ##look at the first 50 dives
##add to dives
dives <- merge(dives, max.d, by.x="num", by.y="num", all.x=T)
colnames(dives) <- c("num", "depth","daytime", "max.dep") #rename the variables
#head(dives)

##extract dives less than 5m maximum depth
# identify dives with max.dep > 5m 
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Should I only use dives >2m or perhaps >6m...
# Want to compare with binned data so sticking with 5m is best I think.
d.list <- which(dives$max.dep > 10)
# #retain only these records
dives_10m <- dives[d.list,]
#dives_5m[1:50,] ##look at the first 50 dives

#How many dives are we keeping from our analysis?
length(unique(dives_10m$num))
# length(unique(dives2$num))



#####################################################################
################# Heerah Broken stick model #########################

## POSIXct format
## column names: date-hour = "daytime", id = "seal", depth = "depth", dive number = "num", temperature = "temp", 
## ambient light = "light" etc. 
## The only variables needed for the algorithm are: daytime, depth and dive number. All the others are optionnal

#dt <- dt[,c(1,2,5)] ## We only keep "daytime", "depth" and "num" columns.

#-----------------------------------------------------------------------------------------------------------------------------------------

##### Optimised Broken Stick Algorithm

# 1. Creation of output dataframes

d_env <- dives_10m ## put your dataframe in the d_env variable which is used troughout the script

num = unique(d_env$num) ## number of dives you have
num.list <- num

dbs <- data.frame("num"=rep(0,1), "all.dur"=0,"start"=0,"end"=0,"depth_start"=0,"depth_end"=0,"seg"=0,"npoints"=0,
                  "dur"=0,"dur.per"=0,"coef"=0, "mean_depth"=0,"max.depth"=0,"sinuosity"=0,"mean_err"=0,"foraging"=0) ## Broken stick dataframe

## num = dive number, all.dur = total dive duration, start = date of segment start, end = date of segment start, 
## depth_start = depth of segment start, seg = broken stick segment number, npoints = number of points summarising the dive
## dur = duration of each segment, dur.per = % of total duration, coef = slope coefficient of the segment, 
## mean_depth = mean depth of segment, max.depth = maximum dive depth, sinuosity = vertical sinuosity associated to this part of the dive,
## mean_err = mean distance between original dive profile and the reconstructed one for the optimal number of broken stick points,
## foraging = behaviour according to vertical sinuosity threshold

ncdv = data.frame("daytime"=rep(0,1) ,"depth"=0,"num"=0) ## dataframe in which the dives for which the fit doesn't work will be stored.
## Needs to have a same column names exactly
d_env <- d_env[,c("daytime","depth", "num")]
#-----------------------------------------------------------------------------------------------------------------------------------------
# loop for each dive

for(d in 1:length(num.list)){
  print(d)
  dt <- d_env[d_env$num==num.list[d],]
  if(nrow(dt) > 60 ) {  ##consider dives of more than 60s as the resolution of the dataset is 1s
    ndive=num.list[d]
    
    #plot(as.numeric(dt$daytime), dt$depth, ylim=c(max(dt$depth),0), t="l", ylab="depth (m)", xlab="",xaxt="n") 
    #Use only if you want to check your dive
    #axis.POSIXct(1,x=dt$daytime, format="%H:%M:%S", labels = TRUE,cex.lab=0.5)
    #idem
    
    np <- c(3:30)  ## number of broken stick iterations to see which optimal number of points summarise your dive
    npe=rep(NA,28) ## vector where the average distance between original and reconstructed dive profile is stored
    npo=rep(NA,28) ## vector where the number of broken stick points describing the dive profile is stored
    
    
    ### Finding the optimal number of Broken points for each dive
    
    # 1.Loop to define the mean distance depending on the number of broken stick points 
    
    for (k in 1:length(np)){
      npp = np[k] # selection of the number of iteration: from 3 to 30
      # 2 lines below: selection of the depth and time for the 2 surface points and the maximum depth point
      ref <- c(dt$depth[1],max(dt$depth),dt$depth[nrow(dt)])
      tim <-c(as.numeric(dt$daytime[1]),as.numeric(dt$daytime[dt$depth==max(dt$depth)][1]),as.numeric(dt$daytime[nrow(dt)]))
      
      for (i in 1:npp){
        #plot(as.numeric(dt$daytime), dt$depth, ylim=c(max(dt$depth),0), t="l", ylab="depth (m)", xlab="",xaxt="n")
        # plot only if you want to see how the broken stick algorithm is working
        #points(tim,ref, pch=19, cex=1, col="red")
        #idem
        interp <- approx(tim,ref,xout=dt$daytime,method="linear") #linear interpolation between broken stick points at TDR time interval
        #lines(interp,col="red")
        #idem
        dif_x <- as.numeric(interp$x - dt$daytime) # time differences between original and reconstructed profiles
        dif_y <- interp$y - dt$depth # depth differences between original and reconstructed profiles
        dst <- sqrt(dif_x^2 + dif_y^2) # calculate distances between original and reconstructed profiles
        
        ii <- which(dst==max(dst))[1] # index of the data point of maximum difference between original and reconstructed profiles
        #points(dt$daytime[ii],dt$depth[ii],col="blue",pch=19,cex=1)
        #idem
        tim <- c(as.numeric(tim),as.numeric(dt$daytime[ii])) # add new broken stick point time
        tim <- ISOdatetime(1970,1,1, 0,0,0, tz="gmt") + tim 
        ref <- c(ref,dt$depth[ii]) # add new broken stick point depth
      }
      npe[k] = mean(dst) # average distance between original and reconstructed dive profiles
      npo[k] = length(tim) # number of broken stick points describing the dive profile
    }
    
    
    # 2. Defining the optimal number of broken stick points
    
    f <- data.frame(npe=npe, npo=npo)
    #plot(f$npo, f$npe,xlab="nb of points", ylab="mean error") #plot of mean distance between original and reconstructed dive profiles
    # according to the number of broken stick points describing the dive
    #activate only if you want to check
    
    # Use of a gompertz model to find the curve which best fit our data
    Asym <- 0; b2 <- -5; b3 <- 0.9
    fm1 <- -999
    try(fm1 <- nls(npe ~ SSgompertz(npo, Asym, b2, b3), data=f, control=nls.control(maxiter=500)),TRUE) #gompertz model to fit an asymptote
    #curve to the mean distance between original and reconstructed dive profiles plot
    
    if (class(fm1) == "nls"){ # if the model converged, we can go to the next steps
      #summary(fm1)
      tt <-predict(fm1, f$npe)
      
      
      # plot of the mean distance 
      png(paste(fig_path,"WED_BS_",ndive,"_",substr(dt$daytime[1],1,10),".png", sep=""),pointsize=12*1.5,height=480*1.5,width=480*1.5)
      
      par(mfrow=c(2,1),mar=c(4,4,2,2))
      tit=paste("BS_WED08_",ndive,"_",substr(dt$daytime[1],1,10))
      plot(f$npo, f$npe,xlab="nb of points", ylab="mean error",main=tit)
      lines(na.omit(f$npo),tt[1:28],col="red")
      
      # Plot the linear approximation between the first and last point of the fitted curve
      t <- data.frame(npe=c(f$npe[1], f$npe[28]), npo=c(f$npo[1], f$npo[28]))
      interp <- approx(c(f$npo[1],f$npo[28]),c(tt[1],tt[28]), xout=f$npo,method="linear")
      interp$x <- interp$x[!is.na(interp$x)]
      interp$y <- interp$y[!is.na(interp$y)]
      lines(interp$x, interp$y,col="blue")
      
      # Looking for the inflexion point which is the furthest point between the fitted curve and the approximation 
      dif_x <- interp$x - na.omit(f$npo)
      dif_y <- interp$y - tt[1:28]
      dst <- sqrt(dif_x^2 + dif_y^2)
      dm <- f$npo[which(dst==max(dst))]
      
      points(f$npo[which(dst==max(dst))], f$npe[which(dst==max(dst))], pch=19, col="red") ## inflexion point
      
      
      #3. optimal broken stick method for each dive 
      
      # The two lines below select the optimal number of broken stick points (in their order of appearance in the BS iteration)
      # example: surface start point, max. depth point, surface end point + x other points
      tim= tim[1:dm] 
      ref=ref[1:dm]
      
      tim2 <- sort(tim)
      dep_tim <- as.data.frame(cbind(ref,tim))
      dep_tim <- dep_tim[order(tim),]
      
      dbs2 <- data.frame("num"=rep(0,(nrow(dep_tim)-1)) , "all.dur"=0,"start"=0,"end"=0,"depth_start"=0,"depth_end"=0,"seg"=0,"npoints"=0, "dur"=0,"dur.per"=0,"coef"=0, "mean_depth"=0, "max.depth"=0,"sinuosity"=0,"mean_err"=0)
      
      # Loop to calculate the different metrics for each broken stick segments
      for (n in 1:(nrow(dep_tim)-1)){
        x1= dep_tim$tim[n] # start of BS segment
        x2= dep_tim$tim[n+1] #end of BS segment
        dbs2$num[n]=ndive
        dbs2$all.dur[n]=difftime(dt$daytime[nrow(dt)], dt$daytime[1], tz,units = c("secs")) #dive duration
        dbs2$start[n]=x1 
        dbs2$end[n]=x2
        dbs2$depth_start[n]= dep_tim$ref[n] # depth of start of BS segment
        dbs2$depth_end[n]= dep_tim$ref[n+1] # depth of end of BS segment
        dbs2$seg[n]=n #segment number
        dbs2$npoints[n]=nrow(dep_tim) # optimal BS points summarising the original dive profile
        dbs2$dur[n]= difftime(tim2[n+1], tim2[n], tz,units = c("secs")) #duration of the segment in sec.
        dbs2$dur.per[n]=(dbs2$dur[n]/dbs2$all.dur[n])*100 #% of segment duration according to total dive duration
        dbs2$coef[n]=(dep_tim$ref[n+1] - dep_tim$ref[n])/(x2 - x1) # slope coefficient of the segment
        dbs2$mean_depth[n]=mean(dt$depth[which(as.numeric(dt$daytime)==x1):which(as.numeric(dt$daytime)==x2)]) #mean depth of the segment 
        # calculated from original profile depths
        dbs2$max.depth[n]= max(dt$depth) # dive max. depth
        
        #Calculation of vertical sinuosity
        deuc= abs(dep_tim$ref[n+1] - dep_tim$ref[n]) # Vertical distance swum between 2 BS points
        dobs=sum(abs(diff(dt$depth[which(dt$daytime==x1):which(dt$daytime==x2)]))) # sum of all the vertical distances from the original
        #profile between the two corresponding BS depth points
        
        dbs2$sinuosity[n]=deuc/dobs # vertical sinuosity index
        dbs2$mean_err[n]=f$npe[which(dst==max(dst))] # mean distance between original and reconstructed dive profiles for the optimal
        #number of BS points summarising the dive.
      }
      
      #-----------------------------------------------------------------------------------------------------------------------------------
      # IMPORTANT:
      #-----------
      # Attribution of behaviour according to vertical sinuosity -- Remind that the sinuosity threshold used here was determined according
      # to the histogram/density plot of vertical sinuosity for every BS segments of every dive
      # so, before setting your threshold at 0.9, check if it suits your dataset (i.e after running the BS on all your dives)
      
      dbs2$foraging <- 2 ## 2 stands for "hunting" mode
      dbs2$foraging[dbs2$sinuosity >=0.9 & dbs2$sinuosity <=1] <-1 ## 1 stands for "transit" mode
      #-----------------------------------------------------------------------------------------------------------------------------------      
      
      # Dive plot: original dive profile and Broken stick reconstructed profile
      sg <- unique(dbs2$seg)
      cl <- c("blue","red")
      dbs2$code[dbs2$sinuosity]
      
      plot(as.numeric(dt$daytime), dt$depth, ylim=c(max(dt$depth),0),t="l",ylab="depth (m)", xlab="",xaxt="n")
      points(tim,ref, pch=19, cex=1, col="black")
      lines(approx(tim,ref,xout=dt$daytime,method="linear"),col="black")
      for(i in 1:length(sg)){
        lines(c(dbs2$start[dbs2$seg==sg[i]],dbs2$end[dbs2$seg==sg[i]]),c(dbs2$depth_start[dbs2$seg==sg[i]],
                                                                         dbs2$depth_end[dbs2$seg==sg[i]]),col=cl[dbs2$foraging][dbs2$seg==sg[i]],lwd=2.5)
      }
      axis.POSIXct(1,x=dt$daytime, format="%H:%M:%S",labels = TRUE,cex.lab=0.5)
      dev.off()
      dbs <-rbind(dbs,dbs2)
    } else {ncdv<-rbind(ncdv,dt)} # allows to keep somewhere the data for which the fit of the Gompertz model didn't work ? 
  } #end of if loop for 60 s
  #save(dbs, file="BS_fitmet_WED_08_samp.RData")
} #end of for loop for dive number

dbs<- dbs[-1,]

hist(dbs$sinuosity,xlab="sinuosity",breaks=seq(0,1,0.1),main="") ## See line 191
abline(v=0.9,col="red",lwd=2)

ncdv <- ncdv[-1,]
#save(ncdv, file="BS_err_WED_08_samp.RData")
#save(dbs, file="BS_fitmet_WED_08_samp.RData")
mean_for_depth <- dbs[,c("dur","mean_depth","mean_err","coef","foraging")]
mds.list <- which(mean_for_depth$foraging==2 & mean_for_depth$coef<0.2 
                  & mean_for_depth$dur>20 & mean_for_depth$mean_depth>6)
mds.list <- which(dbs$start >"15:00:00" & dbs$start<"00:00:00")
#mfd.list <- mean_for_depth[,"foraging"==2]
mean_for_depth <- mean_for_depth[mds.list,]
#names(dt)[names(dt)=="Depth"] <- "depth"
hist(mean_for_depth$mean_depth,xlab="mean_depth",breaks=seq(6,80,1),main="") ## See line 191
#abline(v=0.9,col="red",lwd=2)
length(mean_for_depth$mean_depth)
dbs$start <- ISOdatetime(1970,1,1,0,0,0, tz="GMT")+dbs$start
dbs$end <- ISOdatetime(1970,1,1,0,0,0, tz="GMT")+dbs$end