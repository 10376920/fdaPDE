# TEST
# Invoking smooth.FEM.basis with covariates != NULL in some cases causes
# errors and/or warnings.

library(fdaPDE)

order = 1
mesh = create.MESH.2D(nodes=rbind(c(0, 0), c(0, 1), c(0.5, 0.5), c(1, 1), c(1, 0)),
                      segments=rbind(c(1, 2), c(2, 3), c(3, 4), c(4, 5), c(5,1)),
                      order = order)
FEMbasis = create.FEM.basis(mesh)

locations = rbind(c(0, 0), c(0, 1), c(0.5, 0.5), c(1, 1), c(1, 0))
observations = c(1,2,1,2,1)
covariates1 = cbind(c(1,2,3,4,5))
covariates2 = cbind(c(1,2,3,4,5), c(5,4,3,2,1))

lambda1 = c(1)
lambda2 = c(1,2)

# ABBREVIATIONS OF ERRORS AND WARNINGS
# (F) fnhat:
#       "Error in fnhat[, i] : incorrect number of dimensions"
# (B) betahat:
#       "In betahat[i] = as.vector(lm.fit(covariates, as.vector(observations -  :
#           number of items to replace is not a multiple of replacement length"
# (N) NaNs:
#       "In sqrt(stderr2) : NaNs produced"

# SUMMARY
#
# - 1 covariate:
#               |   loc not     |   loc         |
#               |   on nodes    |   on nodes    |
#   ---------------------------------------------
#   1 lambda    |       -       |       F       |
#   ---------------------------------------------
#   2 lambdas   |       -       |       -       |
#   ---------------------------------------------
#
# - 2 covariates:
#               |   loc not     |   loc         |
#               |   on nodes    |   on nodes    |
#   ---------------------------------------------
#   1 lambda    |       BN      |      FBN      |
#   ---------------------------------------------
#   2 lambdas   |       BN      |      BN       |
#   ---------------------------------------------


# errors: none
# warnings: none
print("locations not on nodes, 1 lambda, 1 covariate")
output1 =
    smooth.FEM.basis(observations = observations,
                     locations=locations,
                     FEMbasis = FEMbasis,
                     lambda = lambda1,
                     covariates=covariates1,
                     GCV = TRUE,
                     CPP_CODE = TRUE)

# errors: none
# warnings: betahat, NaNs
print("locations not on nodes, 1 lambda, 2 covariates")
    smooth.FEM.basis(observations = observations,
                     locations = locations,
                     FEMbasis = FEMbasis,
                     lambda = lambda1,
                     covariates=covariates2,
                     GCV = TRUE,
                     CPP_CODE = TRUE)

# errors: none
# warnings: none
print("locations not on nodes, 2 lambdas, 1 covariate")
output3 =
    smooth.FEM.basis(observations = observations,
                     locations = locations,
                     FEMbasis = FEMbasis,
                     lambda = lambda2,
                     covariates=covariates1,
                     GCV = TRUE,
                     CPP_CODE = TRUE)

# errors: none
# warnings: betahat, Nans
print("locations not on nodes, 2 lambdas, 2 covariate")
output4 =
    smooth.FEM.basis(observations = observations,
                     locations = locations,
                     FEMbasis = FEMbasis,
                     lambda = lambda2,
                     covariates=covariates2,
                     GCV = TRUE,
                     CPP_CODE = TRUE)

# errors: fnhat
# warnings: none
print("locations on nodes, 1 lambda, 1 covariate")
output5 =
    smooth.FEM.basis(observations = observations,
                     FEMbasis = FEMbasis,
                     lambda = lambda1,
                     covariates=covariates1,
                     GCV = TRUE,
                     CPP_CODE = TRUE)

# errors: fnhat
# warnings: betahat, NaNs
print("locations on nodes, 1lambda, 2 covariates")
output6 =
    smooth.FEM.basis(observations = observations,
                     FEMbasis = FEMbasis,
                     lambda = lambda1,
                     covariates=covariates2,
                     GCV = TRUE,
                     CPP_CODE = TRUE)

# errors: none
# warnings: none
print("locations on nodes, 2 lambdas, 1 covariate")
output7 =
    smooth.FEM.basis(observations = observations,
                     FEMbasis = FEMbasis,
                     lambda = lambda2,
                     covariates=covariates1,
                     GCV = TRUE,
                     CPP_CODE = TRUE)

# errors: none
# warnings: betahat, NaNs
print("locations on nodes, 2 lambdas, 2 covariates")
output8 =
    smooth.FEM.basis(observations = observations,
                     FEMbasis = FEMbasis,
                     lambda = lambda2,
                     covariates=covariates2,
                     GCV = TRUE,
                     CPP_CODE = TRUE)
