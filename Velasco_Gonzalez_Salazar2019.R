# Genarate a virtual species
library(raster)
library(rgdal)
library(ecospat)
library(dismo)
library(ENMeval)
library(rJava)
source("pkgs.R") # functions to generate ellipsoids and rasters
source("fun_tutorial.R") # functions to generate ellipsoids and rasters

# load a shapefile
shape <- readOGR(dsn = ".", layer = "aplodontia_rufa_P")
apl_ruf_p<-shape@data
apl_ruf_p.c<-subset(apl_ruf_p, bio_10 > 0)
pts<-apl_ruf_p.c # points
pts2<-pts[,c(2:3)]

# load climatic layers
varclim<-list.files(path="~/Julian_Velasco/InPrep/ModelComplex/Capas&Species/capas",  pattern="asc", full.names=T)
varclim<-stack(varclim)

# to convert layers to points l
bios_points <- rasterToPoints(varclim[[1]], fun=NULL, spatial=TRUE)
bios_all <- extract(varclim, bios_points)
bios_all <- data.frame(coordinates(bios_points),bios_all)

# select only 3 climatic variables
treeVars <- varclim[[c("bio_1","bio_6","bio_12")]]

#to generate a fundamental realized ecological niche

# to calculate covariance matrix from 3 variables
cent_cov_sp <- cov_cent(sp_clim,mve=T,level=0.99,nvars=3, interactive=T)

# fit an ellipsoid model
suits_sp <-ellipsoidfit(data=treeVars,level=0.999,threshold=0.01, plot=T,size=2.2, 
                            centroid=cent_cov_sp$centroid,covar=cent_cov_sp$covarice)
# plot in a map                  
plot(suits_sp$suitRaster)
points(sp_clim[,c("lon", "lat")], pch=19, cex=0.1, col="black")

# to project the true niche in the geography
suits_sp_points <- suits_sp$suits
suits_sp_raster <- suits_sp$suitRaster

# to reclassify the true niche in a binary map (presence-absence)
bin_sp_raster<-reclassify(suits_sp_raster, c(-Inf,0,0,0,Inf,1))

# to select presences (suitability > 0)
sample_sp<-subset(suits_sp_points, suitability > 0)
#sample.size <- floor(0.1 * nrow(sample_sp)) # to take a sample of 10%
sample.size <- 200 # sample of 200 presence points
set.seed(123)
train <- sample(seq_len(nrow(sample_sp)), size = sample.size)
train_sp1 <- sample_sp[train, ] # training 
test_sp1 <- sample_sp[-train, ] # testing
test_sp1<-sample_sp[sample(nrow(test_sp1), 200),] # to take 200 presences to validation

# calculate prevalence to estimate the number of absences
virtual.pres<-freq(bin_sp_raster, value=1)
virtual.abs<-freq(bin_sp_raster, value=0)
preval<-virtual.pres/(virtual.pres+virtual.abs)
abs.numb<-((1-preval)*10)*200 # absences number

# to select absences from the true distribution (binary map)
points.bin<-rasterToPoints(bin_sp_raster, fun=NULL, spatial=TRUE)
points.presabs<-extract(bin_sp_raster,points.bin)
points.presabs<-data.frame(coordinates(points.bin), points.presabs)
ausencias_sp<-subset(points.presabs, points.presabs==0)
abs_sp1<-ausencias_sp[sample(nrow(ausencias_sp), abs.numb),] # extract absences according to prevalence

test_sp1["PresAbs"] <- 1
abs_sp1["PresAbs"] <- 0
test_sp1 <- test_sp1[,c(1:2,8)]
abs_sp1 <- abs_sp1[,c(1:2,4)]
testing.data<-rbind(test_sp1, abs_sp1) # dataset to validate estimated models


# to estimated models with different features and regularization values
pts<-train_sp1[,1:2]
coordinates(pts)<-~x+y
make.args=c(extrapolate=FALSE) # to set Maxent to avoid extrapolation

# to generate models and calculate AICc values
sp1enm.rkf <-ENMevaluate(pts, varclim, bg.coords=NULL, occ.grp=NULL, 
                     bg.grp=NULL, RMvalues=seq(0.5, 4, 0.5), 
                     fc=c("L","LQ","H","LQH","LQHP","LQHPT"),
                     categoricals=NULL, n.bg=10000, method="randomkfold", 
                     overlap=TRUE, aggregation.factor=c(1,1), 
                     kfolds=2, bin.output=TRUE, clamp=TRUE, rasterPreds=TRUE, make.args)

pred.sp1<-sp1enm.rkf@predictions
results.sp1<-sp1enm.rkf@results

# function to standardized models (0 to 1)
stand.biomod<-function(raster)
{
  ras.max=cellStats(raster,"max")
  ras.std=raster/ras.max
  return(ras.std)
}

list_r<-c("L_0_5","LQ_0_5","H_0_5","LQH_0_5","LQHP_0_5","LQHPT_0_5","L_1","LQ_1",
          "H_1","LQH_1","LQHP_1","LQHPT_1","L_1_5","LQ_1_5","H_1_5","LQH_1_5",
          "LQHP_1_5","LQHPT_1_5","L_2","LQ_2","H_2","LQH_2","LQHP_2","LQHPT_2",
          "L_2_5","LQ_2_5","H_2_5","LQH_2_5","LQHP_2_5","LQHPT_2_5","L_3",
          "LQ_3", "H_3","LQH_3","LQHP_3","LQHPT_3","L_3_5","LQ_3_5","H_3_5",
          "LQH_3_5","LQHP_3_5","LQHPT_3_5","L_4","LQ_4","H_4","LQH_4","LQHP_4","LQHPT_4")


# standardized all models (0 to 1)
pred.sp1.s<-stack(pred.s<-sapply(1:48,function(x) stand.biomod(pred.sp1[[x]]))) 

# to check all models have the same size and resolution that the true map
pred.sp1.s.crop<-crop(pred.sp1.s, suits_sp_raster)

# reclassify estimated models using the same threshold that real species
mtp<-0.0100006146948381 # Aplodontia rufa
mtp<-0.01009563 # Arborimus albipes

pred.bin.same<-stack(pred.bin<-sapply(1:48,function(x) 
      reclassify(pred.sp1.s.crop[[x]], c(-Inf,mtp,0, mtp,Inf,1)))) 
names(pred.bin.same)<-list_r

# save all rasters
writeRaster(stack(pred.bin.same), names(pred.bin.same), bylayer=TRUE, format="ascii", overwrite=TRUE)	

# to calculate values for each threshold criteria
test.data<-testing.data
coordinates(test.data)<-~x+y
test.data.suits<- as.data.frame(sapply(1:48, function(x) extract(pred.sp1.s.crop[[x]], test.data)))
test.suits<-cbind(testing.data, test.data.suits)
abs.pred<-subset(test.suits, PresAbs==0)
pres.pred<-subset(test.suits, PresAbs==1)
abs.pred<-abs.pred[,c(-1:-3)]
pres.pred<-pres.pred[,c(-1:-3)]
valid<-sapply(1:48, function(x) evaluate(p=pres.pred[,x], a=abs.pred[,x]))
thresholds<-sapply(1:48, function(x) threshold(valid[[x]]))

# Minimum training presence for all
mtp.all<-as.data.frame(sapply(1:48,function(x) min(extract(pred.sp1.s.crop[[x]], occs))))

# Maximum training sensitivity plus specificity
spec_sens<-as.data.frame(as.numeric(sapply(1:48, function(x) thresholds[2,x])))

# Equal sensitivity and specificity
equal_sens_spec<-as.data.frame(as.numeric(sapply(1:48, function(x) thresholds[5,x])))

# Kappa is highest
kappa_high<-as.data.frame(as.numeric(sapply(1:48, function(x) thresholds[1,x])))

# to generate binary predictions using each threshold value

# minimum training presence
pred.bin.mtp<-stack(pred.bin<-sapply(1:48,function(x) 
  reclassify(pred.sp1.s.crop[[x]], c(-Inf,mtp.all[x,],0, mtp.all[x,],Inf,1)))) 
names(pred.bin.mtp)<-list_r
writeRaster(stack(pred.bin.mtp), names(pred.bin.mtp), bylayer=TRUE, format="ascii", overwrite=TRUE)

# maximize especificidad y sensibilidad
pred.bin.spec_sens<-stack(pred.bin<-sapply(1:48,function(x) 
  reclassify(pred.sp1.s.crop[[x]], c(-Inf,spec_sens[x,],0, spec_sens[x,],Inf,1)))) 
names(pred.bin.spec_sens)<-list_r
writeRaster(stack(pred.bin.spec_sens), names(pred.bin.spec_sens), bylayer=TRUE, format="ascii", overwrite=TRUE)

# igual sensibilidad y especificidad
pred.bin.equal_sens_spec<-stack(pred.bin<-sapply(1:48,function(x) 
  reclassify(pred.sp1.s.crop[[x]], c(-Inf,equal_sens_spec[x,],0, equal_sens_spec[x,],Inf,1)))) 
names(pred.bin.equal_sens_spec)<-list_r
writeRaster(stack(pred.bin.equal_sens_spec), names(pred.bin.equal_sens_spec), bylayer=TRUE, format="ascii", overwrite=TRUE)

# kappa maximo
pred.bin.kappa<-stack(pred.bin<-sapply(1:48,function(x) 
  reclassify(pred.sp1.s.crop[[x]], c(-Inf,kappa_high[x,],0, kappa_high[x,],Inf,1)))) 
names(pred.bin.kappa)<-list_r
writeRaster(stack(pred.bin.kappa), names(pred.bin.kappa), bylayer=TRUE, format="ascii", overwrite=TRUE)

# VALIDATION METRICS

# Validation metrics - Threshold same virtual
validation<-testing.data
# function to calculate all performance metrics
allMetrics <-  lapply(1:48, function(x){
  
  return(confu_mat_con(rasterBinary = pred.bin.same[[x]],valData = validation ,x = "x",y="y"))
  
})

allMetricsDF <- do.call("rbind.data.frame",allMetrics)
results_all_same_virtual<-cbind(results.sp1, allMetricsDF)
write.csv(results_all, file="results_sp2_t_same_virtual.csv")

# Validation metrics - MTP
validation<-testing.data
# function to calculate all performance metrics
allMetrics <-  lapply(1:48, function(x){
  
  return(confu_mat_con(rasterBinary = pred.bin.sp1.crop.mtp[[x]],valData = validation ,x = "x",y="y"))
  
})

allMetricsDF <- do.call("rbind.data.frame",allMetrics)
results_all_same_mtp<-cbind(results.sp1, allMetricsDF)
write.csv(results_all, file="results_sp2_t_mtp.csv")

# Validation metrics - same Espec + Sens
validation<-testing.data
# function to calculate all performance metrics
allMetrics <-  lapply(1:48, function(x){
  
  return(confu_mat_con(rasterBinary = pred.bin.sp1.crop.mtp[[x]],valData = validation ,x = "x",y="y"))
  
})

allMetricsDF <- do.call("rbind.data.frame",allMetrics)
results_all_same_Esp_Sens<-cbind(results.sp1, allMetricsDF)
write.csv(results_all, file="results_sp2_t_equal_spec_sens.csv")

# Validation metrics - max Espec + Sens
validation<-testing.data
# function to calculate all performance metrics
allMetrics <-  lapply(1:48, function(x){
  
  return(confu_mat_con(rasterBinary = pred.bin.sp1.crop.mtp[[x]],valData = validation ,x = "x",y="y"))
  
})

allMetricsDF <- do.call("rbind.data.frame",allMetrics)
results_all_max_Esp_Sens<-cbind(results.sp1, allMetricsDF)
write.csv(results_all, file="results_sp2_t_max_sens_spec.csv")

# Validation metrics - max Kappa
validation<-testing.data
# function to calculate all performance metrics
allMetrics <-  lapply(1:48, function(x){
  
  return(confu_mat_con(rasterBinary = pred.bin.sp1.crop.mtp[[x]],valData = validation ,x = "x",y="y"))
  
})

allMetricsDF <- do.call("rbind.data.frame",allMetrics)
results_all_max_kappa<-cbind(results.sp1, allMetricsDF)
write.csv(results_all, file="results_sp2_t_max_kappa.csv")






# inter pixel x pixel
real<-real.bin
bin.est<-list.files("C:\Users\fcupul\Documents\Julian_Velasco\InPrep\ModelComplex\Umbrales",full.names = T)

# Ejecutar la funcion para contar interseccion pixel a pixel
results <- lapply(bin.est, function(x){
  modelo <- raster(x)
  return(compare_maps(mapReal = real,mapMod = modelo))
  
})





# Comparaciones nichos (espacio E) de real vs. estimada

real.p <- rasterToPoints(real.bin, fun=NULL, spatial=TRUE)
real.pts <- extract(real.bin, real.p)
real.pts <- data.frame(coordinates(real.p),real.pts)

pred.bin.p<-rasterToPoints(pred.bin.sp1.crop[[1]], fun=NULL, spatial=TRUE)
pred.pts <- extract(pred.bin.sp1.crop, pred.bin.p)
pred.pts <- data.frame(coordinates(pred.bin.p),pred.pts)

bad.pred.pts<-pred.pts[,c(1,2,8)] # LQHPT 0.5 highest AICc
good.pred.pts<-pred.pts[,c(1,2,36)] # LQH 3 lowest AICc

poor.clim<-extract(varclim, bad.pred.pts[,-3])
poor.pts.clim<-cbind(bad.pred.pts, poor.clim)
poor.pts.clim.pres<-subset(poor.pts.clim, LQHPT_0_5==1)

best.clim<-extract(varclim, good.pred.pts[,-3])
best.pts.clim<-cbind(good.pred.pts, best.clim)
best.pts.clim.pres<-subset(best.pts.clim, LQH_3==1)

real.clim<-extract(varclim, real.pts[,-3])
real.pts.clim<-cbind(real.pts, real.clim)
real.pts.clim.pres<-subset(real.pts.clim, real.pts==1)

pdf("Aplondotia rufa P.pdf")
par(mfrow=c(2,2))
plot(best.pts.clim$bio_1, best.pts.clim$bio_12, col="gray", main="Lowest AICc", pch=15, xlab="bio1", ylab="bio12")
points(best.pts.clim.pres$bio_1, best.pts.clim.pres$bio_12, col="red", pch=15)
plot(poor.pts.clim$bio_1, poor.pts.clim$bio_12, col="gray", main="Highest AICc", pch=15, xlab="bio1", ylab="bio12")
points(poor.pts.clim.pres$bio_1, poor.pts.clim.pres$bio_12, col="red", pch=15)
plot(real.pts.clim$bio_1, real.pts.clim$bio_12, col="gray", main="Real", pch=15, xlab="bio1", ylab="bio12")
points(real.pts.clim.pres$bio_1, real.pts.clim.pres$bio_12, col="red", pch=15)
dev.off()

# Suitability comparisons (Broennimann)

bios_points <- rasterToPoints(varclim[[1]], fun=NULL, spatial=TRUE)
bios_all<-extract(varclim,bios_points)
bios_all<-data.frame(coordinates(bios_points), bios_all)

# densidad de ocurrencias
best.pres<-best.pts.clim.pres[,1:2]
poor.pres<-poor.pts.clim.pres[,1:2]
real.pres<-real.pts.clim.pres[,1:2]

best.pres["species"]<-"best_model"
poor.pres["species"]<-"worst_model"
real.pres["species"]<-"real_model"

best.pres.s<-best.pres[sample(nrow(best.pres), 1000),] # sacar 500 pts
real.pres.s<-real.pres[sample(nrow(real.pres), 1000),] # sacar 500 pts
poor.pres.s<-poor.pres

pair1<- rbind(best.pres.s, real.pres.s)
pair2<- rbind(poor.pres, real.pres.s)

write.csv(pair1, file="pair1.csv")
pair1<-read.csv("pair1.csv")
head(pair1)

write.csv(pair2, file="pair2.csv")
pair2<-read.csv("pair2.csv")
head(pair2)

sp1.N<-subset(pair1, species=="best_model")
sp2.N<-subset(pair1, species=="real_model")
sp3.N<-subset(pair2, species=="worst_model")

#eliminar puntos a menos de 10 km de distancia
occ.sp1.N <- ecospat.occ.desaggregation(dfvar=sp1.N,colxy=2:3,colvar=NULL, min.dist=0.05,plot=FALSE)
occ.sp2.N <-ecospat.occ.desaggregation(dfvar=sp2.N,colxy=2:3,colvar=NULL, min.dist=0.05,plot=FALSE)
occ.sp3.N <-ecospat.occ.desaggregation(dfvar=sp3.N,colxy=2:3,colvar=NULL, min.dist=0.05,plot=FALSE)


occ.sp_test <- na.exclude(ecospat.sample.envar(dfsp=pair2,colspxy=2:3,colspkept=1:4,dfvar=bios_all,
                                               colvarxy=1:2,colvar="all",resolution=0.1))

occ.sp_test<-occ.sp_test[,-1]
# lista de especies
sp.list<-levels(occ.sp_test[,3])

sp.nbocc<-c()
for (i in 1:length(sp.list)){sp.nbocc<-c(sp.nbocc,length(which(occ.sp_test[,3] == sp.list[i])))} 

nb.sp<-length(sp.list) #nb of species

# Que variables incluir en los analisis?
Xvar<-c(3:21)
nvar<-length(Xvar)
#numero de iteraciones
iterations<-100
#resolucion de espacio ambiental
R=100

### PCA-ENVIRONMENT ###
data<-rbind(occ.sp_test[,Xvar+1],bios_all[,Xvar])
w1<-c(rep(0,nrow(occ.sp_test)),rep(1,nrow(bios_all)))

#vector of weight, 0 for the occurences, 1 for the sites of the study area
pca.cal <-dudi.pca(data, row.w = w1, center = T, scale = T, scannf = F, nf = 2)

# the pca is calibrated on all the sites of the study area
# occurences are not used for the calibration, but their scores are calculated

sp.combn1<-combn(1:2,2)
for(i in 1:ncol(sp.combn1)) {
  row.sp1.p1<-which(occ.sp_test[,3] == sp.list[sp.combn1[1,i]]) # rows in data corresponding to sp1
  row.sp2.p1<-which(occ.sp_test[,3] == sp.list[sp.combn1[2,i]]) # rows in data corresponding to sp2
  name.sp1.p1<-sp.list[sp.combn1[1,i]]
  name.sp2.p1<-sp.list[sp.combn1[2,i]]
  # predict the scores on the axes
  scores.clim.p1<- pca.cal1$li[(nrow(occ.sp_test)+1):nrow(data),] #scores for global climate
  scores.sp1.p1<- pca.cal1$li[row.sp1.p1,] #scores for sp1
  scores.sp2.p1<- pca.cal1$li[row.sp2.p1,] #scores for sp2
}


#Dynamic Occurrence Densities Grid ecospat.grid.clim.dyn
z1.dyn.p1<-ecospat.grid.clim.dyn (scores.clim.p1, scores.clim.p1, scores.sp1.p1, R=100)
z2.dyn.p2<-ecospat.grid.clim.dyn (scores.clim.p1, scores.clim.p1, scores.sp2.p1, R=100)

# test of niche equivalency and similarity according to Warren et al. 2008
a.dyn.p1 <- ecospat.niche.equivalency.test(z1.dyn.p1, z2.dyn.p2, rep = 100) #niche equivalency
b.dyn.p1 <- ecospat.niche.similarity.test(z1.dyn.p1, z2.dyn.p2, rep = 100) #niche similarity 1 to 2
c.dyn.p1 <- ecospat.niche.similarity.test(z2.dyn.p2, z1.dyn.p1, rep = 100) #niche similarity 2 to 1

#Plot Species Density ecospat.plot.niche for one species
z1<- ecospat.grid.clim.dyn(scores.clim.p1,scores.clim.p1,scores.sp1.p1,R)
z2<- ecospat.grid.clim.dyn(scores.clim.p1,scores.clim.p1,scores.sp2.p1,R)


pdf("Aplodontia_rufa_worst_vs_real.pdf", width=12, height=7)
par(mfrow=c(3,3))
ecospat.plot.niche (z1, title="worst", name.axis1="PC1", name.axis2="PC2", cor=F)
ecospat.plot.niche (z2, title="real", name.axis1="PC1", name.axis2="PC2", cor=F)
ecospat.plot.contrib(pca.cal$co,pca.cal$eig)
ecospat.plot.overlap.test(a.dyn.p1, "I", "Equivalency")
ecospat.plot.overlap.test(b.dyn.p1, "I", "Similarity 1 -> 2")
ecospat.plot.overlap.test(c.dyn.p1, "I", "Similarity 2 -> 1")
ecospat.plot.overlap.test(a.dyn.p1, "D", "Equivalency")
ecospat.plot.overlap.test(b.dyn.p1, "D", "Similarity 1 -> 2")
ecospat.plot.overlap.test(c.dyn.p1, "D", "Similarity 2 -> 1")
dev.off()




# plots validation
quartz()
par(mfrow=c(2,2))
plot(results_all$tss, results_all$AICc, xlab="TSS", ylab="AICc")
plot(results_all$tss, results_all$delta.AICc, xlab="TSS", ylab="Delta AICc")
plot(results_all$tss, results_all$w.AIC, xlab="TSS", ylab="weights")
plot(results_all$tss, results_all$nparam, xlab="TSS", ylab="nparam")

quartz()
par(mfrow=c(2,2))
plot(results_all$Mean.AUC, results_all$AICc, xlab="mAUC", ylab="AICc")
plot(results_all$Mean.AUC, results_all$delta.AICc, xlab="mAUC", ylab="Delta AICc")
plot(results_all$Mean.AUC, results_all$w.AIC, xlab="mAUC", ylab="weights")
plot(results_all$Mean.AUC, results_all$nparam, xlab="mAUC", ylab="nparam")

quartz()
par(mfrow=c(2,2))
plot(results_all$kappa, results_all$AICc, xlab="Kappa", ylab="AICc")
plot(results_all$kappa, results_all$delta.AICc, xlab="Kappa", ylab="Delta AICc")
plot(results_all$kappa, results_all$w.AIC, xlab="Kappa", ylab="weights")
plot(results_all$kappa, results_all$nparam, xlab="Kappa", ylab="nparam")


save.image("Aplondotia_rufa_P.RData")



