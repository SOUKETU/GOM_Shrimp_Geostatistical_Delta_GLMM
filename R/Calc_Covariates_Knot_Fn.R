
Calc_Covariates_Knot_Fn <-
  function( loc_x, static_covariates, dynamic_covariates, Data_Extrap ){

    # Nearest extrapolation grid for each knot
    NN_Extrap = nn2( data=loc_x[,c('E_km','N_km')], query=Data_Extrap[,c('E_km','N_km')], k=1 )
    
    # Calculate covariate for each knot
    # Covariate j at location x is the average value in the extrapolation_grid for all grid cells s that are nearest to location x
    X_xj = matrix(NA, ncol=ncol(static_covariates)-2, nrow=nrow(loc_x))
      for(j in 1:ncol(X_xj)){
        X_xj[,j] = tapply(static_covariates[,2+j], INDEX=factor(NN_Extrap$nn.idx,levels=1:nrow(loc_x)), FUN=mean, na.rm=TRUE)
      }
    
    X_xtp = array(NA, c(nrow(loc_x),dim(dynamic_covariates)[2],dim(dynamic_covariates)[3]))
    for (z in 1:dim(dynamic_covariates)[3]){
      for (i in 1:dim(dynamic_covariates)[2]){
        X_xtp[,i,z] = tapply(dynamic_covariates[,i,z], INDEX=factor(NN_Extrap$nn.idx,levels=1:nrow(loc_x)), FUN=mean, na.rm=TRUE)
      }
    }
    
    # Return stuff
    Return = list( "X_xj"=X_xj, "X_xtp"=X_xtp )
  }
