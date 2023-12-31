#generates matrices with different resolutions for vegetation maps, from 1/6lat*1/4lon to coarser resolution
#2/21/04 chg
#
#  $Id: gen_veg_mat.r,v 1.2 2007-06-27 11:54:03 skoerner Exp $
#---------------------------------------------------------------------------------------------------

library(foreign) #need to use data.restore...

if(exists("vegpath")){
  data.restore(paste(vegpath,"vegflux.dmp",sep=""),print=T)
  } else {stop("gen_veg_mat: please define 'vegpath' for reading/writing flux and vegetation grids")}

#vegflux.dmp contains veg.1, veg.2, ..., veg.31 (INPUT FIELDS FOR THIS SCRIPT)
#all at high resolution (1/4 lon, 1/6 lat)
#veg.1, ... , veg.17 are IGBP vegetation classes
#[1]  Evergreen needleleaf forest        
#[2]  Evergreen broadleaf forest        
#[3]  Deciduous needleleaf forest        
#[4]  Deciduous broadleaf forest        
#[5]  Mixed forest                       
#[6]  Closed shrublands                  
#[7]  Open shrublands                    
#[8]  Woody savannas
#[9]  Savannas                           
#[10] Grasslands                        
#[11] Permanent wetlands                 
#[12] Croplands                         
#[13] Urban and built-up                 
#[14] Cropland/Natural vegetation mosaic
#[15] Snow and ice                       
#[16] Barren or sparsely vegetated      
#[17] Water bodies                      
#[18] CO2 fossil fuel emission  (GEIA 1 deg)
#[19] CO fossil fuel emission (NAPAP 20 km res., GEIA 1 deg)
#[20-31] CH4 emission for 12 months (JAN-DEC) (EDGAR/GISS, provided by John Miller, NOAA)

#OUTPUT FIELDS generated by this script:
#veg.vv.xyrr
#Names for different vegetation types and different flux fields:
#vv= 1-31 for veg. classes and flux fields (see above under INPUT)
#xy= 1-3, 1-3 (for grids which are divided up): xpart can be west, middle, east; ypart can be south, middle, north
#rr: resolution 1-16 according to table below
#rr:	  (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)
#coarsex_c(1,1,2,2,2,4,4,4,8, 8, 8,16,16,16,32,32) #factors by which grids have been made coarser
#coarsey_c(1,2,1,2,4,2,4,8,4, 8,16, 8,16,32,16,32)#e.g., '4' means grid is 4 times coarser
#    e.g. veg.19.312 is CO emission inventory, xpart=3 (east), ypart=1 (south), resolution rr=2 (i.e. coarser by factor 2 in y, that is 2*1/6 lat=1/3 lat))


#create 16 element vectors for x and y reduction
coarsex<-c(1,1,2,2,2,4,4,4,8,8,8,16,16,16,32,32)#factors by which grids have been made coarser
coarsey<-c(1,2,1,2,4,2,4,8,4,8,16,8,16,32,16,32)#e.g., '4' means grid is 4 times coarser
#run all vegetations
#for(veg in 1:1){
for(veg in 1:31){
	veg.raw<-get(paste("veg.",as.character(veg),sep=""))
#run all different resolutions
	for (n in 16:1){
	lat.old<-as.numeric(dimnames(veg.raw)[[1]])
	lon.old<-as.numeric(dimnames(veg.raw)[[2]])
	gty.ngm<-1:324			#get NGM area only
	lat.ngm<-(gty.ngm-1)/6+11	#south-west corner represents gridcell
	gtx.ngm<-1:376			#get NGM area only
	lon.ngm<-(gtx.ngm-1)/4-145	#south-west corner represents gridcell
	sel.lat<-(lat.old>=min(lat.ngm))&(lat.old<=max(lat.ngm))
	sel.lon<-(lon.old>=min(lon.ngm))&(lon.old<=max(lon.ngm))
	veg.raw<-veg.raw[sel.lat,sel.lon]	#get NGM area only

	#create matrix for results; assumes south-west corner to represent grid
	gty.new<-floor((gty.ngm+coarsey[n]-1)/coarsey[n])*coarsey[n]-coarsey[n]+1
	gtx.new<-floor((gtx.ngm+coarsex[n]-1)/coarsex[n])*coarsex[n]-coarsex[n]+1
	
	result<-matrix(0,nrow=length(unique(gty.new)),ncol=length(unique(gtx.new)))
	minx<-unique(gtx.new);miny<-unique(gty.new)
#add shifted matrices together and take average
#shift matrices
for (dy in 1:coarsey[n]){
	for (dx in 1:coarsex[n]){
#first make sure not to run over the grid limit
gtx.shift<-minx+dx-1; if(gtx.shift[length(gtx.shift)]>376)gtx.shift[length(gtx.shift)]<-376
gty.shift<-miny+dy-1; if(gty.shift[length(gty.shift)]>324)gty.shift[length(gty.shift)]<-324
result<-result+veg.raw[gty.shift,gtx.shift]
}#for gtx
}#for gty
result<-result/coarsex[n]/coarsey[n] #turn sum into average

#check whether result is identical in orientation and pattern
#motif()
#image(t(result))
#motif()
#image(t(veg.raw))
dimnames(result)<-list(miny,minx)#south-west corner is reference for grid cell
#first generate complete fields, indicated by "00"
assignr(paste("veg.",as.character(veg),".00",as.character(n),sep=""),result,path=vegpath,printTF=T)

#generate divided fields, indicated by "xpart ypart" in name
if(n==1){xparts<-3;yparts<-3}
if(n==2){xparts<-3;yparts<-1}
if(n==3){xparts<-1;yparts<-3}
if(n<=3){
for (xpart in 1:xparts){
for (ypart in 1:yparts){
selx<-minx>(xpart-1)*376/(xparts+1)&minx<=(xpart+1)*376/(xparts+1)
#selx<-minx>=(xpart-1)*376/(xparts+1)&minx<=(xpart+1)*376/(xparts+1)
sely<-miny>(ypart-1)*324/(yparts+1)&miny<=(ypart+1)*324/(yparts+1)
#sely<-miny>=(ypart-1)*324/(yparts+1)&miny<=(ypart+1)*324/(yparts+1)
assignr(paste("veg.",as.character(veg),".",as.character(xpart),as.character(ypart),as.character(n),sep=""),result[sely,selx],path=vegpath,printTF=T)
}#for ypart
}#for xpart
}#for if(n<=3)
print(memory.size())
	}#of run all different resolutions
}#of run all vegetations
