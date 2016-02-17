
# Load libraries
library(TMB)
library(INLA)
library(SpatialDeltaGLMM)
library(ThorsonUtilities)

file.sources = list.files(path="/Users/jiecao/Desktop/Google Drive/work/Postdoc research/GOM_shrimp_spatial/GOM_Shrimp_Geostatistical_Delta_GLMM/R",pattern="*.R",full.names=TRUE)
for (f in file.sources) {
  source(f)
  print(f)}

###############
# Settings
###############
covariate=1

Data_Set = c("GOM_shrimp")
Sim_Settings = list("Species_Set"=1:100, "Nyears"=10, "Nsamp_per_year"=600, "Depth_km"=-1, "Depth_km2"=-1, "Dist_sqrtkm"=0, "SigmaO1"=0.5, "SigmaO2"=0.5, "SigmaE1"=0.5, "SigmaE2"=0.5, "SigmaVY1"=0.05, "Sigma_VY2"=0.05, "Range1"=1000, "Range2"=500, "SigmaM"=1)
Version = "geo_index_v3l"
n_x = c(100, 250, 500, 1000, 2000)[2] # Number of stations
FieldConfig = c("Omega1"=0, "Epsilon1"=0, "Omega2"=1, "Epsilon2"=1) # 1=Presence-absence; 2=Density given presence
CovConfig = c("SST"=0, "RandomNoise"=0) # DON'T USE DURING REAL-WORLD DATA FOR ALL SPECIES (IT IS UNSTABLE FOR SOME)
Q_Config = c("Pass"=0)
VesselConfig = c("Vessel"=0, "VesselYear"=0)
ObsModel = 5  # 0=normal (log-link); 1=lognormal; 2=gamma; 4=ZANB; 5=ZINB; 11=lognormal-mixture; 12=gamma-mixture
Kmeans_Config = list( "randomseed"=1, "nstart"=100, "iter.max"=1e3 )     # Samples: Do K-means on trawl locs; Domain: Do K-means on extrapolation grid

# Determine region
Region = switch( Data_Set,"GOM_shrimp"="West_GOM", "SAWC_jacopever"="South_Africa", "Sim"="California_current")

strata.limits = list('All_areas'=c(4100, 4080, 4060, 4030, 4010, 4050, 4040, 4070, 4020, 4090, 4110, 4120),
                     'Six_strata'=c(4100, 4080, 4060, 4030, 4050, 4070))

# This is where all runs will be located
if(covariate==1){
  DateFile = paste('/Users/jiecao/Desktop/Google Drive/work/Postdoc research/GOM_shrimp_spatial/Geostatistical_model_results','/',Sys.Date(),'-','covariates','/',sep='')
}else{
  DateFile = paste('/Users/jiecao/Desktop/Google Drive/work/Postdoc research/GOM_shrimp_spatial/Geostatistical_model_results','/',Sys.Date(),'-','no-covariates','/',sep='')
}
dir.create(DateFile)


# Compile TMB software
TmbDir = system.file("executables", package="SpatialDeltaGLMM")
TmbDir = "/Users/jiecao/Desktop/Google Drive/work/Postdoc research/GOM_shrimp_spatial/Geostatistical_model_results"
setwd( TmbDir )
compile( paste(Version,".cpp",sep="") )

################
# Prepare data
# (THIS WILL VARY FOR DIFFERENT DATA SETS) 
################
setwd('/Users/jiecao/Desktop/Google Drive/work/Postdoc research/GOM_shrimp_spatial/GOM_Shrimp_Geostatistical_Delta_GLMM')

# Read in shrimp summer survey data
load("West_GOM_Shrimp_survey.rda")
#shrimp_catch <- shrimp_data_new$B0 + shrimp_data_new$B1 + shrimp_data_new$B2 + shrimp_data_new$B3 + shrimp_data_new$B4 + shrimp_data_new$B5 + shrimp_data_new$B6 +shrimp_data_new$B7
shrimp_catch <- shrimp_data_new$N0 + shrimp_data_new$N1 + shrimp_data_new$N2 + shrimp_data_new$N3 + shrimp_data_new$N4 + shrimp_data_new$N5 + shrimp_data_new$N6 +shrimp_data_new$N7

Data_Geostat = data.frame( "Catch_KG"=shrimp_catch, "Year"=shrimp_data_new[,'Year'], 
                           "Vessel"="missing", "AreaSwept_km2"=0.93*0.03, 
                           "Lat"=shrimp_data_new[,'Lat'], 
                           "Lon"=shrimp_data_new[,'Lon'],"Depth"=shrimp_data_new[,'Depth'],
                           "Bottom_temp"=shrimp_data_new[,'Temp'])
Year_Set = sort(unique(Data_Geostat[,'Year']))

# Get extrapolation data

if( Region == "West_GOM" ){
  Extrapolation_List = Prepare_WGOM_Extrapolation_Data_Fn( strata.limits=strata.limits )
}

#########################################################################################
#########################################################################################

# Convert to an Eastings-Northings coordinate system
tmpUTM = Convert_LL_to_UTM_Fn( Lon=Data_Geostat[,'Lon'], Lat=Data_Geostat[,'Lat'], zone=Extrapolation_List[["zone"]] )                                                         #$
Data_Geostat = cbind( Data_Geostat, 'E_km'=tmpUTM[,'X'], 'N_km'=tmpUTM[,'Y'])

# Calculate k-means centroids (but only once for all species)
Kmeans = Calc_Kmeans(n_x=n_x, loc_orig=Data_Geostat[,c("E_km", "N_km")], randomseed=Kmeans_Config[["randomseed"]], nstart=Kmeans_Config[["nstart"]], iter.max=Kmeans_Config[["iter.max"]], DirPath=DateFile)
loc_x = Kmeans$centers
Data_Geostat = cbind( Data_Geostat, "knot_i"=Kmeans$cluster )

# Calc design matrix and areas
PolygonList = Calc_Polygon_Areas_and_Polygons_Fn(loc_x=loc_x, Data_Extrap=Extrapolation_List[["Data_Extrap"]], Covariates=c("none"), a_el=Extrapolation_List[["a_el"]])
a_xl = PolygonList[["a_xl"]]
NN_Extrap = PolygonList[["NN_Extrap"]]

# Calc average covariates
load("static_covariates.rda")
load("dynamic_covariates.rda")
if (covariate==1){
  CovariateList = Calc_Covariates_Knot_Fn(loc_x=loc_x,static_covariates=static_covariates, dynamic_covariates=dynamic_covariates,Data_Extrap=Extrapolation_List[["Data_Extrap"]])
  X_xj = as.matrix(CovariateList[["X_xj"]])
  #X_xj = NULL
  #X_xtp = NULL
  X_xtp = array(CovariateList[["X_xtp"]][,,2], c(nrow(loc_x),dim(dynamic_covariates)[2],1))
  #X_xtp = CovariateList[["X_xtp"]]
}else{
  X_xj = NULL
  X_xtp = NULL
}

# Make mesh and info for anisotropy
MeshList = Calc_Anisotropic_Mesh(loc_x=loc_x)

################
# Make and Run TMB model
# (THIS WILL BE SIMILAR FOR EVERY DATA SET) 
################

# Save options for future records
Record = bundlelist( c("Data_Set","Sim_Settings","Version","n_x","FieldConfig","CovConfig","Q_Config","VesselConfig","ObsModel","Kmeans_Config") )
capture.output( Record, file=paste0(DateFile,"Record.txt"))

# Make TMB data list
TmbData = Data_Fn("Version"=Version, "FieldConfig"=FieldConfig, "ObsModel"=ObsModel, "b_i"=Data_Geostat[,'Catch_KG'], "a_i"=Data_Geostat[,'AreaSwept_km2'], "v_i"=as.numeric(Data_Geostat[,'Vessel'])-1, "s_i"=Data_Geostat[,'knot_i']-1, "t_i"=Data_Geostat[,'Year']-min(Data_Geostat[,'Year']), "a_xl"=a_xl, "X_xj"=X_xj,"X_xtp"=X_xtp, "MeshList"=MeshList)

# Make TMB object
TmbList = Build_TMB_Fn(TmbData, TmbDir=TmbDir, Version=Version, VesselConfig=VesselConfig, loc_x=loc_x, DiagnosticDir=DateFile)
Obj = TmbList[["Obj"]]

# Run first time -- marginal likelihood
Start_time = Sys.time()
Obj$fn(Obj$par)
# Run first time -- gradient with respect to fixed effects
Obj$gr(Obj$par)

#Obj$he(Obj$par)

# Run model
Opt = nlminb(start=Obj$par, objective=Obj$fn, gradient=Obj$gr, lower=TmbList[["Lower"]], upper=TmbList[["Upper"]], control=list(eval.max=1e4, iter.max=1e4, trace=1))  # , rel.tol=1e-20
Opt[["final_diagnostics"]] = data.frame( "Name"=names(Opt$par), "Lwr"=TmbList[["Lower"]], "Est"=Opt$par, "Upr"=TmbList[["Upper"]], "Gradient"=Obj$gr(Opt$par) )
Opt[["total_time_to_run"]] = Sys.time() - Start_time
Opt[["number_of_coefficients"]] = c("Total"=length(unlist(Obj$env$parameters)), "Fixed"=length(Obj$par), "Random"=length(unlist(Obj$env$parameters))-length(Obj$par) )
capture.output( Opt, file=paste0(DateFile,"Opt.txt"))

# Reports
Report = Obj$report()                                      
Sdreport = sdreport(Obj, bias.correct=TRUE)

# Save stuff
Save = list("Opt"=Opt, "Report"=Report, "Sdreport"=Sdreport, "ParHat"=Obj$env$parList(Opt$par), "TmbData"=TmbData)
save(Save, file=paste0(DateFile,"Save.RData"))
capture.output( Opt, file=paste0(DateFile,"Opt.txt"))
capture.output( Sdreport, file=paste0(DateFile,"Sdreport.txt"))
file.copy( from=paste0(system.file("executables", package="SpatialDeltaGLMM"),"/",Version,".cpp"), to=paste0(DateFile,Version,".cpp"), overwrite=TRUE)

################
# Make diagnostic plots
################

Plot_data_and_knots (Data_Extrap=Extrapolation_List[["Data_Extrap"]], Extrap_Area_km2=Extrapolation_List[["Area_km2_x"]], 
                     loc_x=loc_x, Data_Geostat=Data_Geostat, Plot_name=paste0(DateFile,"/Data_and_knots.png"))

# Plot Anisotropy  
if( TmbData$Options_vec['Aniso']==1 ){
  PlotAniso_Fn( FileName=paste0(DateFile,"Aniso.png"), Report=Report )
}

# Plot surface
Dim = c("Nrow"=ceiling(sqrt(length(Year_Set)))); Dim = c(Dim,"Ncol"=ceiling(length(Year_Set)/Dim['Nrow']))
par( mfrow=Dim )
MapDetails_List = MapDetails_Fn( "Region"=Region, "NN_Extrap"=NN_Extrap, "Extrapolation_List"=Extrapolation_List )
PlotResultsOnMap_Fn(plot_set=1:3, MappingDetails=MapDetails_List[["MappingDetails"]], Report=Report, PlotDF=MapDetails_List[["PlotDF"]], MapSizeRatio=MapDetails_List[["MapSizeRatio"]], Xlim=MapDetails_List[["Xlim"]], Ylim=MapDetails_List[["Ylim"]], FileName=paste0(DateFile,"Field_"), Year_Set=Year_Set, Rotate=MapDetails_List[["Rotate"]], mfrow=Dim, mar=c(0,0,2,0), oma=c(3.5,3.5,0,0), Cex=MapDetails_List[["Cex"]])

# Plot index
PlotIndex_Fn( DirName=DateFile, TmbData=TmbData, Sdreport=Sdreport, Year_Set=Year_Set, strata_names=ls(strata.limits) )

# Positive catch rate Q-Q plot
Q = QQ_Fn( TmbData=TmbData, Report=Report, FileName_PP=paste0(DateFile,"Posterior_Predictive.jpg"), FileName_Phist=paste0(DateFile,"Posterior_Predictive-Histogram.jpg"), FileName_QQ=paste0(DateFile,"Q-Q_plot.jpg"), FileName_Qhist=paste0(DateFile,"Q-Q_hist.jpg"))

# Covariate effect
#PlotCov_Fn(Report=Report, NN_Extrap=NN_Extrap, X_xj=X_xj, FileName=paste0(DateFile,"Cov_"))


Plot_range_quantiles( Data_Extrap=Extrapolation_List[["Data_Extrap"]], Report=Report, a_xl=a_xl, NN_Extrap=NN_Extrap, 
                      Year_Set=NULL, Prob_vec=c(0.10,0.5,0.90), FileName_Quantiles=paste0(DateFile,"/Range_boundary.png") )
  
Plot_range_shifts(Sdreport, Report, TmbData, Year_Set=NULL, FileName_COG=paste0(DateFile,"/center_of_gravity.png"), FileName_Area=paste0(DateFile,"/Area.png"), Znames=rep("",ncol(Report$mean_Z_tm)),use_biascorr=FALSE)









