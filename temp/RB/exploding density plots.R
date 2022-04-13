
##------ LOAD LIBRARIES

library(rgeos) #--in addition, need some of the libraries from  OutputCompilaton.R

##------ DEFINE THE FUNCTIONS

ExplodeRasterPolygon<-function(r,poly,poly.name.field,center.xy,fac,type="internal"){
  poly.field.list<-sort(as.character(unique(poly@data[,poly.name.field])))
  r.list<-poly.list<-shifted.r.list<-shifted.poly.list<-list()
  r<-crop(r,poly)
  r<-mask(r,poly)
  # plot(r)
  # plot(poly,add=TRUE)
  # i<-poly.field.list[1]
  i<-10
  tick<-0
  for(i in 1:length(poly)){
    try(
      {
        tick<-tick+1
        poly.list[[tick]]<-poly[tick,]#poly@data[,poly.name.field]%in%i,]
        r.temp<-crop(r,poly.list[[tick]])
        r.list[[tick]]<-mask(r.temp,poly.list[[tick]])
        
      },
      silent=TRUE)
    
  }
  
  to.keep<-which(unlist(lapply(r.list,function(x)!is.null(x))))
  r.list<-r.list[to.keep]
  poly.list<-poly.list[to.keep]
  
  # #center.xy<-apply(coordinates(r)[!is.na(r[]),],2,median)
  # 
  # #center.xy<-c(mean(extent(r)[c(1,2)]), mean(extent(r)[c(3,4)]))
  # abline(v=mean(extent(r)[c(1,2)]))
  # abline(h=mean(extent(r)[c(3,4)]))
  # points(center.xy,pch=19,cex=2)
  # x<-r.list[[1]]
  
  
  center.poly.xy<-do.call(rbind,lapply(poly.list,function(x){
    coordinates(gCentroid(x))
  }))
  
  # poly.field.list<-which(unlist(lapply(r.list,function(x)!is.null(x))))
  # plot(poly)
  # plot(extent(r),add=TRUE,col="red")
  #i<- poly.field.list[1]
  i<-1
  for(i in 1:length(poly.list)){
    
    temp<-r.list[[i]]
    proj4string( temp)<-proj4string(poly)
    temp.xy<-apply(coordinates(temp),2,mean)
    
    if(type=="internal"){
      a<-temp.xy[1]-center.xy[1]
      b<-temp.xy[2]-center.xy[2]
      sign.a<-ifelse(a>=0,1,-1)
      sign.b<-ifelse(b>=0,1,-1)
      c<-sqrt(a^2+b^2)
      C<-c+c#*fac
      
      a.shift<-  sign.a*abs((C-c)*a)/abs(a)*fac[1]#sign.a*abs((abs(a) - abs(a*C/sqrt(a^2+b^2))))
      b.shift<- sign.b*abs((C-c)*b)/abs(b)*fac[2]#sign.b*abs(a.shift)/abs(a)
      print(i)
      print(a.shift)
      print(b.shift)
    }
    
    
    if(type=="uniform"){
      
      
      
      
      
      x.bounds<-   c(center.xy[,1] - (center.xy[,1]-range(center.poly.xy[,1])[1] ) *fac[1],
                     center.xy[,1] + (range(center.poly.xy[,1])[2]- center.xy[,1]) *fac[1])
      
      
      y.bounds<-   c(center.xy[,2] - (center.xy[,2]-range(center.poly.xy[,2])[1] ) *fac[2],
                     center.xy[,2] + (range(center.poly.xy[,2])[2]- center.xy[,2]) *fac[2])
      
      #y.bounds<-  center.xy[,2] + (diff(range(center.poly.xy[,2]))*fac[2]*c(-0.5,0.5))
      
      a<-DoScale(center.poly.xy[,1],x.bounds[1],x.bounds[2])
      b<-DoScale(center.poly.xy[,2],y.bounds[1],y.bounds[2])
      
      #a<-seq(x.bounds[1],x.bounds[2],length.out=length(poly.field.list) )   
      #b<-seq(y.bounds[1],y.bounds[2],length.out=length(poly.field.list) )   
      
      a.shift<-  a[i]-center.poly.xy[i,1]# sign.a*abs((C-c)*a)/abs(a)#sign.a*abs((abs(a) - abs(a*C/sqrt(a^2+b^2))))
      b.shift<-  b[i]-center.poly.xy[i,2]#sign.b*abs((C-c)*b)/abs(b)#sign.b*abs(a.shift)/abs(a)
      
      # plot(poly.sp)
      # points(center.poly.xy,pch=19,col="pink")
      # points(b~a,pch=19,col="blue")
      # text(a,b,1:length(poly))
      # 
      # plot(b~a,pch=19,col="blue")
      # points(center.poly.xy,pch=19,col="red")
      # points(center.xy,pch=19,col="red",cex=2)
      # plot(poly.sp,add=TRUE)
    }
    
    
    #seq(-10,10,length.out=2)
    if(type=="h"){
      new.loc<-seq(center.xy[1]-fac,center.xy[1]+fac,length.out=length(poly.list))[i]
      a.shift<-new.loc-center.xy[1]
      b.shift<-0
    }
    
    #--SHIFT RASTER
    temp<-shift(temp,dx=a.shift,dy=b.shift)
    shifted.r.list[[i]]<-temp
    # plot(poly.sp)
    # plot(temp,add=TRUE,legend=FALSE)
    # points(center.poly.xy[i,2]~center.poly.xy[i,1],pch=19,col="red")
    # points(b[i]~a[i],pch=19,col="blue")
    #print(i)
    #locator(1)
    
    #--SHIFT POLYGON
    temp<-shift(poly.list[[i]],dx=a.shift,dy=b.shift)
    shifted.poly.list[[i]]<-temp
    #plot(temp,add=TRUE,legend=FALSE)
    
  }
  #plot(poly,add=TRUE)
  #points(center.xy[2]~center.xy[1],cex=4,pch=19)  
  out<-list(shifted.r.list=shifted.r.list,shifted.poly.list=shifted.poly.list)
  return(out)
  
}


DoScale<-function(data,l=0,u=1){
  x<-data+min(data,na.rm=TRUE)
  
  if(all(is.na(x)))
  {
    scaled.x<-rep(NA,length(data))
  }
  else
  {
    min.x<-min(x,na.rm=TRUE)
    max.x<-max(x,na.rm=TRUE)
    scaled.x<-(u-l)*(x-min.x)/(max.x-min.x)+l
  }
  
  return(scaled.x)
  
}


##------  INPUT FORMATTING

regions.sp<-as(regions,"Spatial")
regions.center.xy<-gCentroid(regions.sp,byid =TRUE)


alps.sp<-as(alps,"Spatial")
regions.sp$dummy<-"dummy"
studyarea.sp<- aggregate(regions.sp,by="dummy")
notalps.sp<-gDifference(studyarea.sp,alps.sp)
notalps.sf <- st_cast(st_as_sf(notalps.sp), "POLYGON")
notalps.sp<-as(notalps.sf,"Spatial")
notalps.sp@data$id<-1:length(notalps.sp)
notalps.sp<-notalps.sp[raster::area(notalps.sp)>1e+6,]
alps.sp$id<-"alps"
allalps.sp<-raster::bind(list(notalps.sp,alps.sp))
allalps.sp$id[allalps.sp$id!="alps"]<-"not.alps"
allalps.sp<-raster::aggregate(allalps.sp,by="id")
allalps.center.xy<-gCentroid(allalps.sp,byid =TRUE)


studyarea.center.xy<-gCentroid(studyarea.sp,byid =TRUE)


##------ COLOR DEFINITIONS 
  
    max <- max(meanDensity.R[],na.rm=TRUE)
    cuts <- seq(0,max,length.out = 100)   #set breaks
    colfunc<- colorRampPalette(c("black","orange","yellow","white"))
    col <- colfunc(100)
    
    
   
##---- REGIONAL EXPLOSION   
    
    graphics.off()
    
    
    
    
    path<-file.path(analysisDir, modelName,"inkscape/Fig_region_exploded.png")
    
     png(file = path,
         width = 15, height = 10, units="in",  pointsize = 12,res=300)

    
    
    head(poly.sp)
    
    out<- ExplodeRasterPolygon(r=meanDensity.R
                               ,
                               poly=regions.sp
                               ,
                               poly.name.field="DEN_UTS"
                               ,
                               center.xy=coordinates(regions.center.xy)
                               ,
                                 fac=0.95* c(35,8)
                               
                               ,
                               type="uniform"#
    )
    
    
    #---to fix some weird offset
    temp<-gCentroid(do.call(raster::bind,out$shifted.poly.list),byid =FALSE)
    dx<- -(coordinates(temp)[,1] - coordinates(studyarea.center.xy)[,1])
    dy<--(coordinates(temp)[,2] - coordinates(studyarea.center.xy)[,2])
    shifted.poly.list<-lapply(out$shifted.poly.list,function(x)shift(x,dx=dx,dy=dy))
    shifted.raster.list<-lapply(out$shifted.r.list,function(x)shift(x,dx=dx,dy=dy))
    
    exploded.poly.sp <- do.call(raster::bind,shifted.poly.list)
 
    #plot(regions.sp,border="orange")
    #plot(exploded.poly.sp,add=T,col="red")
    
    #  SETUP A REFERENCE POLYGON FOR ALL PLOTS
    
    e <- extent(exploded.poly.sp)
    e.poly <- as(e, 'SpatialPolygons')  
    proj4string(e.poly)<-proj4string(exploded.poly.sp)
    
    
    par(bg="black")
    plot(e.poly)
    
  lapply(shifted.raster.list,image,add=TRUE,legend=FALSE,breaks=c(cuts,max(cuts)+1000), col = col,legend=FALSE)
    
  plot(exploded.poly.sp,add=TRUE,border=grey(0.65),lwd=0.75)

graphics.off()
  



##----ZONAL EXPLOSION   
  

path<-file.path(analysisDir, modelName,"inkscape/Fig_zone_exploded.png")

png(file = path,
    width = 15, height = 10, units="in",  pointsize = 12,res=300)

  head(allalps.sp)
  
  out<- ExplodeRasterPolygon(r=meanDensity.R
                             ,
                             poly=allalps.sp
                             ,
                             poly.name.field="id"
                             ,
                             center.xy=coordinates(allalps.sp)
                             ,
                             fac=c(-0.3,2.5)
                             
                             ,
                             type="uniform"#"internal"#
  )
  
  
  #---to fix some weird offset
  temp<-gCentroid(do.call(raster::bind,out$shifted.poly.list),byid =FALSE)
  dx<- -(coordinates(temp)[,1] - coordinates(studyarea.center.xy)[,1])
  dy<--(coordinates(temp)[,2] - coordinates(studyarea.center.xy)[,2])
  shifted.poly.list<-lapply(out$shifted.poly.list,function(x)shift(x,dx=dx,dy=dy))
  shifted.raster.list<-lapply(out$shifted.r.list,function(x)shift(x,dx=dx,dy=dy))
  
  exploded.poly.sp <- do.call(raster::bind,shifted.poly.list)
  

  par(bg="black")
  plot(e.poly)
  
  lapply(shifted.raster.list,image,add=TRUE,legend=FALSE,breaks=c(cuts,max(cuts)+1000), col = col,legend=FALSE)
  
  plot(exploded.poly.sp,add=TRUE,border=grey(0.65),lwd=0.75)
  
  graphics.off()
  
  
  
  
##---- TOTAL DENSITY MAP
  

  path<-file.path(analysisDir, modelName,"inkscape/Fig_unexploded.png")
  
  png(file = path,
      width = 15, height = 10, units="in",  pointsize = 12,res=300)
  
  par(bg="black")
  
  plot(e.poly)
  
  image(meanDensity.R,add=TRUE,legend=FALSE,breaks=c(cuts,max(cuts)+1000), col = col,legend=FALSE)
  
  #--because there are "holes" in the study area (geometry issues)
  
  temp <- RemoveHolesSp(studyarea.sp)
  plot(temp,add=TRUE,border=grey(0.65),lwd=0.75)
  
  graphics.off()
  
  
##---- LEGEND
  
  
  path<-file.path(analysisDir, modelName,"inkscape/Fig_legend.png")
  
  png(file = path,
      width = 1.75, height = 5, units="in",  pointsize = 12,res=300)
  
  par(mar=c(1,1,1,4),las=1,cex.axis=1.6,xaxs="i",yaxs="i")#,bg="black")
  plot(1,type="n",axes=FALSE,ylim=c(1,100),xlim=c(0,1),xlab="",ylab="")
  resolution<-prod(res(meanDensity.R)/1000)
  num<-seq(min(cuts)/resolution,max(cuts)/resolution,0.05)
  lab<-DoScale(num,1,100)
  segments(0,1:length(col),1,1:length(col),col=col,lwd=5,lend=1)
  axis(4,lab,num*100)
  box()
  graphics.off()
  
  
  
  