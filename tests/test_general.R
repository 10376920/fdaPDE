pathnames = list.files(pattern="[.]R$", path="./R/", full.names=TRUE);
sapply(pathnames, FUN=source);
source("tests/mesh.R")

libraries =
    c("fdaPDE_original",
      "fdaPDE_stochastic",
      "fdaPDE_temp",
      "fdaPDE_woodbury_whole",
      "fdaPDE_MUMPS_whole")

################################################################################
# Edit here

# Index (according to the order of the vector libraries) of the versions to be
# tested
idx_libs_to_test = c(4)
# A vector containing the Ns of the grids to be used
N = c(40)
# The number of observations to be generated (same length as N)
n_observations = c(40)
# The "true" coefficients of the covariates
beta = rbind(4.5, -2.0, 3.2, 5.9)
# The lambda to be used
lambda = c(10)
# The order of FEM
order = 1
################################################################################

n_meshes = length(N)
n_covariates = length(beta)
n_libs_to_test = length(idx_libs_to_test)

dyn.load(paste("fdaPDE_original", .Platform$dynlib.ext, sep = ""))
FEMbasis = vector("list", n_meshes)
locations = vector("list", n_meshes)
covariates = vector("list", n_meshes)
covariates_on_nodes = vector("list", n_meshes)
observations = vector("list", n_meshes)
observations_on_nodes = vector("list", n_meshes)
output1 = vector("list", n_libs_to_test*n_meshes)
output2 = vector("list", n_libs_to_test*n_meshes)
output3 = vector("list", n_libs_to_test*n_meshes)
output4 = vector("list", n_libs_to_test*n_meshes)
set.seed(7)
for (i in 1:n_meshes) {
    grid = mesh_quadratounitario(N[i])
    mesh = create.MESH.2D(nodes=grid$nodes, order = order)
    FEMbasis[[i]] = create.FEM.basis(mesh)
    locations[[i]] = cbind(cbind(runif(n_observations[i],0,1)),
                                cbind(runif(n_observations[i],0,1)))
    covariates[[i]] = matrix(runif(n_observations*n_covariates,-100,100),
                             nrow = n_observations,
                             ncol = n_covariates)
    covariates_on_nodes[[i]] =
        matrix(runif(nrow(mesh$nodes)*n_covariates,-100,100),
               nrow = nrow(mesh$nodes),
               ncol = n_covariates)
    observations[[i]] = sin(2*pi*locations[[i]][,1])
                        + covariates[[i]] %*% beta
                        + rnorm(n = nrow(locations[[i]]), sd = 0.1)
    observations_on_nodes[[i]] = sin(2*pi*mesh$nodes[,1])
                                 + covariates_on_nodes[[i]] %*% beta
                                 + rnorm(n = nrow(mesh$nodes), sd = 0.1)
}

# COVARIATES, LOC NOT ON NODES
if (1) {
    cat("\nCOVARIATES, LOC NOT ON NODES\n\n")
    for (k in 1:n_libs_to_test) {
        dyn.load(paste(libraries[idx_libs_to_test[k]],
                       .Platform$dynlib.ext, sep = ""))
        cat("LIBRARY: ", libraries[idx_libs_to_test[k]], "\n")
        for (i in 1:n_meshes) {
            cat("grid: ", N, "x", N, "nodes\n")
            output_CPP =
            smooth.FEM.basis(observations = observations[[i]],
                             locations=locations[[i]],
                             FEMbasis = FEMbasis[[i]],
                             lambda = lambda,
                             covariates=covariates[[i]],
                             GCV = TRUE,
                             CPP_CODE = TRUE,
                             nrealizations = 1000)
            cat("edf = ", output_CPP$edf, "\n")
            output1[[(k-1)*n_meshes + i]] = output_CPP
        }
    }
}

# COVARIATES, LOC ON NODES
if (0) {
    cat("\nCOVARIATES, LOC ON NODES\n\n")
    for (k in 1:n_libs_to_test) {
        dyn.load(paste(libraries[idx_libs_to_test[k]],
                       .Platform$dynlib.ext, sep = ""))
        cat("LIBRARY: ", libraries[idx_libs_to_test[k]], "\n")
        for (i in 1:n_meshes) {
            cat("grid: ", N, "x", N, "nodes\n")
            output_CPP =
            smooth.FEM.basis(observations = observations_on_nodes[[i]],
                             FEMbasis = FEMbasis[[i]],
                             lambda = lambda,
                             covariates=covariates_on_nodes[[i]],
                             GCV = TRUE,
                             CPP_CODE = TRUE)
            cat("edf = ", output_CPP$edf, "\n")
            output2[[(k-1)*n_meshes + i]] = output_CPP
        }
    }
}

# NO COVARIATES, LOC NOT ON NODES
if (0) {
    cat("\nNO COVARIATES, LOC NOT ON NODES\n\n")
    for (k in 1:n_libs_to_test) {
        dyn.load(paste(libraries[idx_libs_to_test[k]],
                       .Platform$dynlib.ext, sep = ""))
        cat("LIBRARY: ", libraries[idx_libs_to_test[k]], "\n")
        for (i in 1:n_meshes) {
            cat("grid: ", N, "x", N, "nodes\n")
            output_CPP =
            smooth.FEM.basis(observations = observations[[i]],
                             locations=locations[[i]],
                             FEMbasis = FEMbasis[[i]],
                             lambda = lambda,
                             GCV = TRUE,
                             CPP_CODE = TRUE)
            cat("edf = ", output_CPP$edf, "\n")
            output3[[(k-1)*n_meshes + i]] = output_CPP
        }
    }
}

# NO COVARIATES, LOC ON NODES
if (0) {
    cat("\nNO COVARIATES, LOC ON NODES\n\n")
    for (k in 1:n_libs_to_test) {
        dyn.load(paste(libraries[idx_libs_to_test[k]],
                       .Platform$dynlib.ext, sep = ""))
        cat("LIBRARY: ", libraries[idx_libs_to_test[k]], "\n")
        for (i in 1:n_meshes) {
        cat("grid: ", N, "x", N, "nodes\n")
            output_CPP =
            smooth.FEM.basis(observations = observations[[i]],
                             FEMbasis = FEMbasis[[i]],
                             lambda = lambda,
                             GCV = TRUE,
                             CPP_CODE = TRUE)
            cat("edf = ", output_CPP$edf, "\n")
            output4[[(k-1)*n_meshes + i]] = output_CPP
        }
    }
}
