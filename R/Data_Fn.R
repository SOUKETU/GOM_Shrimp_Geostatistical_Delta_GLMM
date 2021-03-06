Data_Fn <-
function( Version, FieldConfig, ObsModel, b_i, a_i, v_i, s_i, t_i, a_xl, X_xj=NULL, X_xtp=NULL, Q_ik=NULL, MeshList, Aniso=1, R2_interpretation=0, RhoConfig=c("Beta1"=0,"Beta2"=0,"Epsilon1"=0,"Epsilon2"=0), Options=c('SD_site_density'=0,'SD_site_logdensity'=0,'Calculate_Range'=0,'Calculate_evenness'=0,'Calculate_effective_area'=0), CheckForErrors=TRUE, Alpha=2 ){

  # Determine dimensions
  n_t = max(t_i) + 1
  n_v = max(v_i) + 1
  n_i = length(b_i)
  n_x = nrow(a_xl)
  n_l = ncol(a_xl)

  # Covariates and defaults
  if( is.null(X_xj) ) X_xj = matrix(0, nrow=n_x, ncol=1)
  if( is.null(X_xtp) ) X_xtp = array(0, dim=c(n_x,n_t,1))
  if( is.null(Q_ik) ) Q_ik = matrix(0, nrow=n_i, ncol=1)
  n_j = ncol(X_xj)
  n_p = dim(X_xtp)[3]
  n_k = ncol(Q_ik)

  # by default, add nothing as Z_xl
  if( Options['Calculate_Range']==FALSE ){
    Z_xm = matrix(0, nrow=nrow(a_xl), ncol=ncol(a_xl) ) # Size so that it works for Version 3g-3j
  }
  if(Options['Calculate_Range']==TRUE ){
    Z_xm = MeshList$loc_x
    message( "Calculating range shift for stratum #1:",colnames(a_xl[1]))
  }

  # Check for bad data entry
  if( CheckForErrors==TRUE ){
    if( !is.matrix(a_xl) | !is.matrix(X_xj) | !is.matrix(Q_ik) ) stop("a_xl, X_xj, and Q_ik should be matrices")
    if( max(s_i)-1 > MeshList$mesh$n | min(s_i)<0 ) stop("s_i exceeds bounds in MeshList")
  }

  # Check for bad data entry
  if( CheckForErrors==TRUE ){
    if( length(b_i)!=n_i | length(a_i)!=n_i | length(v_i)!=n_i | length(s_i)!=n_i | length(t_i)!=n_i ) stop("b_i, a_i, v_i, s_i, or t_i doesn't have length n_i")
    if( nrow(a_xl)!=n_x | ncol(a_xl)!=n_l ) stop("a_xl has wrong dimensions")
    if( nrow(X_xj)!=n_x | ncol(X_xj)!=n_j ) stop("X_xj has wrong dimensions")
    if( nrow(Q_ik)!=n_i | ncol(Q_ik)!=n_k ) stop("Q_ik has wrong dimensions")
    if( length(unique(t_i))!=n_t ) stop("Please renumber t_i consecutively from 0 to final year minus 1")
    if( length(unique(v_i))!=n_v ) stop("Please renumber v_i consecutively from 0 to final vessel minus 1")
  }

  # Output tagged list
  if(Version=="geo_index_v3a"){
    Return = list( "n_i"=n_i, "n_s"=MeshList$spde$n.spde, "n_x"=n_x, "n_t"=n_t, "n_v"=n_v, "n_j"=n_j, "n_k"=n_k, "Aniso"=Aniso, "FieldConfig"=FieldConfig, "ObsModel"=ObsModel, "Options"=Options, "b_i"=b_i, "a_i"=a_i, "v_i"=v_i, "s_i"=s_i, "t_i"=t_i, "a_x"=a_xl[,1], "X_xj"=X_xj, "Q_ik"=Q_ik, "n_tri"=nrow(MeshList$mesh$graph$tv), "Tri_Area"=MeshList$Tri_Area, "E0"=MeshList$E0, "E1"=MeshList$E1, "E2"=MeshList$E2, "TV"=MeshList$TV-1, "G0_inv"=inla.as.dgTMatrix(solve(MeshList$spde$param.inla$M0)), "G0"=MeshList$spde$param.inla$M0, "G1"=MeshList$spde$param.inla$M1, "G2"=MeshList$spde$param.inla$M2 )
  }
  if(Version=="geo_index_v3b"){
    Return = list( "n_i"=n_i, "n_s"=MeshList$spde$n.spde, "n_x"=n_x, "n_t"=n_t, "n_v"=n_v, "n_j"=n_j, "n_k"=n_k, "n_l"=n_l, "Aniso"=Aniso, "FieldConfig"=FieldConfig, "ObsModel"=ObsModel, "Options"=Options, "b_i"=b_i, "a_i"=a_i, "v_i"=v_i, "s_i"=s_i, "t_i"=t_i, "a_xl"=a_xl, "X_xj"=X_xj, "Q_ik"=Q_ik, "n_tri"=nrow(MeshList$mesh$graph$tv), "Tri_Area"=MeshList$Tri_Area, "E0"=MeshList$E0, "E1"=MeshList$E1, "E2"=MeshList$E2, "TV"=MeshList$TV-1, "G0_inv"=inla.as.dgTMatrix(solve(MeshList$spde$param.inla$M0)), "G0"=MeshList$spde$param.inla$M0, "G1"=MeshList$spde$param.inla$M1, "G2"=MeshList$spde$param.inla$M2 )
  }
  if(Version%in%c("geo_index_v3d","geo_index_v3c")){
    Return = list( "n_i"=n_i, "n_s"=MeshList$spde$n.spde, "n_x"=n_x, "n_t"=n_t, "n_v"=n_v, "n_j"=n_j, "n_k"=n_k, "n_l"=n_l, "Aniso"=Aniso, "FieldConfig"=FieldConfig, "ObsModel"=ObsModel, "Options"=Options, "b_i"=b_i, "a_i"=a_i, "v_i"=v_i, "s_i"=s_i, "t_i"=t_i, "a_xl"=a_xl, "X_xj"=X_xj, "Q_ik"=Q_ik, "spde"=list(), "G0"=MeshList$spde$param.inla$M0, "G1"=MeshList$spde$param.inla$M1, "G2"=MeshList$spde$param.inla$M2 )
    Return[['spde']] = list("n_s"=MeshList$spde$n.spde, "n_tri"=nrow(MeshList$mesh$graph$tv), "Tri_Area"=MeshList$Tri_Area, "E0"=MeshList$E0, "E1"=MeshList$E1, "E2"=MeshList$E2, "TV"=MeshList$TV-1, "G0"=MeshList$spde$param.inla$M0, "G0_inv"=inla.as.dgTMatrix(solve(MeshList$spde$param.inla$M0)) )
  }
  if(Version%in%c("geo_index_v3e")){
    Return = list( "n_i"=n_i, "n_s"=MeshList$spde$n.spde, "n_x"=n_x, "n_t"=n_t, "n_v"=n_v, "n_j"=n_j, "n_k"=n_k, "n_l"=n_l, "Options_vec"=c(Aniso,R2_interpretation), "FieldConfig"=FieldConfig, "ObsModel"=ObsModel, "Options"=Options, "b_i"=b_i, "a_i"=a_i, "v_i"=v_i, "s_i"=s_i, "t_i"=t_i, "a_xl"=a_xl, "X_xj"=X_xj, "Q_ik"=Q_ik, "spde"=list(), "G0"=MeshList$spde$param.inla$M0, "G1"=MeshList$spde$param.inla$M1, "G2"=MeshList$spde$param.inla$M2 )
    Return[['spde']] = list("n_s"=MeshList$spde$n.spde, "n_tri"=nrow(MeshList$mesh$graph$tv), "Tri_Area"=MeshList$Tri_Area, "E0"=MeshList$E0, "E1"=MeshList$E1, "E2"=MeshList$E2, "TV"=MeshList$TV-1, "G0"=MeshList$spde$param.inla$M0, "G0_inv"=inla.as.dgTMatrix(solve(MeshList$spde$param.inla$M0)) )
  }
  if(Version%in%c("geo_index_v3g","geo_index_v3f")){
    Options_vec = c("Aniso"=Aniso, "R2_interpretation"=R2_interpretation, "Rho_betaTF"=ifelse(RhoConfig[["Beta1"]]|RhoConfig[["Beta2"]],1,0) )
    Return = list( "n_i"=n_i, "n_s"=MeshList$spde$n.spde, "n_x"=n_x, "n_t"=n_t, "n_v"=n_v, "n_j"=n_j, "n_k"=n_k, "n_l"=n_l, "Options_vec"=Options_vec, "FieldConfig"=FieldConfig, "ObsModel"=ObsModel, "Options"=Options, "b_i"=b_i, "a_i"=a_i, "v_i"=v_i, "s_i"=s_i, "t_i"=t_i, "a_xl"=a_xl, "X_xj"=X_xj, "Q_ik"=Q_ik, "Z_xl"=Z_xm, "spde"=list(), "G0"=MeshList$spde$param.inla$M0, "G1"=MeshList$spde$param.inla$M1, "G2"=MeshList$spde$param.inla$M2 )
    Return[['spde']] = list("n_s"=MeshList$spde$n.spde, "n_tri"=nrow(MeshList$mesh$graph$tv), "Tri_Area"=MeshList$Tri_Area, "E0"=MeshList$E0, "E1"=MeshList$E1, "E2"=MeshList$E2, "TV"=MeshList$TV-1, "G0"=MeshList$spde$param.inla$M0, "G0_inv"=inla.as.dgTMatrix(solve(MeshList$spde$param.inla$M0)) )
  }
  if(Version%in%c("geo_index_v3i","geo_index_v3h")){
    Options_vec = c("Aniso"=Aniso, "R2_interpretation"=R2_interpretation, "Rho_betaTF"=ifelse(RhoConfig[["Beta1"]]|RhoConfig[["Beta2"]],1,0), "Alpha"=Alpha )
    Return = list( "n_i"=n_i, "n_s"=MeshList$spde$n.spde, "n_x"=n_x, "n_t"=n_t, "n_v"=n_v, "n_j"=n_j, "n_k"=n_k, "n_l"=n_l, "Options_vec"=Options_vec, "FieldConfig"=FieldConfig, "ObsModel"=ObsModel, "Options"=Options, "b_i"=b_i, "a_i"=a_i, "v_i"=v_i, "s_i"=s_i, "t_i"=t_i, "a_xl"=a_xl, "X_xj"=X_xj, "Q_ik"=Q_ik, "Z_xl"=Z_xm, "spde"=list(), "spde_aniso"=list(), "G0"=MeshList$spde$param.inla$M0, "G1"=MeshList$spde$param.inla$M1, "G2"=MeshList$spde$param.inla$M2 )
    Return[['spde']] = inla.spde2.matern(MeshList$mesh)$param.inla[c("M0","M1","M2")]
    Return[['spde_aniso']] = list("n_s"=MeshList$spde$n.spde, "n_tri"=nrow(MeshList$mesh$graph$tv), "Tri_Area"=MeshList$Tri_Area, "E0"=MeshList$E0, "E1"=MeshList$E1, "E2"=MeshList$E2, "TV"=MeshList$TV-1, "G0"=MeshList$spde$param.inla$M0, "G0_inv"=inla.as.dgTMatrix(solve(MeshList$spde$param.inla$M0)) )
  }
  if(Version%in%c("geo_index_v3j")){
    Options_vec = c("Aniso"=Aniso, "R2_interpretation"=R2_interpretation, "Rho_betaTF"=ifelse(RhoConfig[["Beta1"]]|RhoConfig[["Beta2"]],1,0), "Alpha"=Alpha, "AreaAbundanceCurveTF"=0 )
    Return = list( "n_i"=n_i, "n_s"=MeshList$spde$n.spde, "n_x"=n_x, "n_t"=n_t, "n_v"=n_v, "n_j"=n_j, "n_k"=n_k, "n_l"=n_l, "Options_vec"=Options_vec, "FieldConfig"=FieldConfig, "ObsModel"=ObsModel, "Options"=Options, "b_i"=b_i, "a_i"=a_i, "v_i"=v_i, "s_i"=s_i, "t_i"=t_i, "a_xl"=a_xl, "X_xj"=X_xj, "Q_ik"=Q_ik, "Z_xl"=Z_xm, "spde"=list(), "spde_aniso"=list() )
    Return[['spde']] = inla.spde2.matern(MeshList$mesh)$param.inla[c("M0","M1","M2")]
    Return[['spde_aniso']] = list("n_s"=MeshList$spde$n.spde, "n_tri"=nrow(MeshList$mesh$graph$tv), "Tri_Area"=MeshList$Tri_Area, "E0"=MeshList$E0, "E1"=MeshList$E1, "E2"=MeshList$E2, "TV"=MeshList$TV-1, "G0"=MeshList$spde$param.inla$M0, "G0_inv"=inla.as.dgTMatrix(solve(MeshList$spde$param.inla$M0)) )
  }
  if(Version%in%c("geo_index_v3k")){
    Options_vec = c("Aniso"=Aniso, "R2_interpretation"=R2_interpretation, "Rho_betaTF"=ifelse(RhoConfig[["Beta1"]]|RhoConfig[["Beta2"]],1,0), "Alpha"=Alpha, "AreaAbundanceCurveTF"=0 )
    Return = list( "n_i"=n_i, "n_s"=MeshList$spde$n.spde, "n_x"=n_x, "n_t"=n_t, "n_v"=n_v, "n_j"=n_j, "n_k"=n_k, "n_l"=n_l, "n_m"=ncol(Z_xm), "Options_vec"=Options_vec, "FieldConfig"=FieldConfig, "ObsModel"=ObsModel, "Options"=Options, "b_i"=b_i, "a_i"=a_i, "v_i"=v_i, "s_i"=s_i, "t_i"=t_i, "a_xl"=a_xl, "X_xj"=X_xj, "Q_ik"=Q_ik, "Z_xm"=Z_xm, "spde"=list(), "spde_aniso"=list() )
    Return[['spde']] = inla.spde2.matern(MeshList$mesh)$param.inla[c("M0","M1","M2")]
    Return[['spde_aniso']] = list("n_s"=MeshList$spde$n.spde, "n_tri"=nrow(MeshList$mesh$graph$tv), "Tri_Area"=MeshList$Tri_Area, "E0"=MeshList$E0, "E1"=MeshList$E1, "E2"=MeshList$E2, "TV"=MeshList$TV-1, "G0"=MeshList$spde$param.inla$M0, "G0_inv"=inla.as.dgTMatrix(solve(MeshList$spde$param.inla$M0)) )
  }
  if(Version%in%c("geo_index_v3l")){
    Options_vec = c("Aniso"=Aniso, "R2_interpretation"=R2_interpretation, "Rho_betaTF"=ifelse(RhoConfig[["Beta1"]]|RhoConfig[["Beta2"]],1,0), "Alpha"=Alpha, "AreaAbundanceCurveTF"=0 )
    Return = list( "n_i"=n_i, "n_s"=MeshList$spde$n.spde, "n_x"=n_x, "n_t"=n_t, "n_v"=n_v, "n_j"=n_j, "n_p"=n_p, "n_k"=n_k, "n_l"=n_l, "n_m"=ncol(Z_xm), "Options_vec"=Options_vec, "FieldConfig"=FieldConfig, "ObsModel"=ObsModel, "Options"=Options, "b_i"=b_i, "a_i"=a_i, "v_i"=v_i, "s_i"=s_i, "t_i"=t_i, "a_xl"=a_xl, "X_xj"=X_xj, "X_xtp"=X_xtp, "Q_ik"=Q_ik, "Z_xm"=Z_xm, "spde"=list(), "spde_aniso"=list() )
    Return[['spde']] = inla.spde2.matern(MeshList$mesh)$param.inla[c("M0","M1","M2")]
    Return[['spde_aniso']] = list("n_s"=MeshList$spde$n.spde, "n_tri"=nrow(MeshList$mesh$graph$tv), "Tri_Area"=MeshList$Tri_Area, "E0"=MeshList$E0, "E1"=MeshList$E1, "E2"=MeshList$E2, "TV"=MeshList$TV-1, "G0"=MeshList$spde$param.inla$M0, "G0_inv"=inla.as.dgTMatrix(solve(MeshList$spde$param.inla$M0)) )
  }

  # Check for NAs
  if( CheckForErrors==TRUE ){
    if( any(sapply(Return, FUN=function(Array){any(is.na(Array))})==TRUE) ) stop("Please find and eliminate the NA from your inputs") 
  }

  # Return
  return( Return )
}
