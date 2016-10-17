library(Rmpi)

pathnames = list.files(pattern="[.]R$", path="./R/", full.names=TRUE);
sapply(pathnames, FUN=source);
source("tests/mesh.R")

libraries =
    c("fdaPDE_original",
      "fdaPDE_stochastic",
      "fdaPDE_temp",
      "fdaPDE_woodbury_whole",
      "fdaPDE_MUMPS_whole",
      "fdaPDE_woodbury_decomposeQ")

################################################################################
# Edit here

# Index (according to the order of the vector libraries) of the versions to be
# tested

idx_libs_to_test = c(4)

# A vector containing the Ns of the grids to be used
N = c(30)
# The number of observations to be generated (same length as N)
n_observations = c(30)
# The "true" coefficients of the covariates
beta = rbind(0.2, -0.4, 0.7, -0.05)
# Functions to be used to generate the covariates
f = vector("list", length(beta-1))
# Specify a function for each covariate -1, which is random
#f[[1]] <- function (x,y){x}
#f[[2]] <- function (x,y){y}
#f[[3]] <- function (x,y){x*y}
f[[1]] <- function (x,y){sin(2*pi*x*y)}
f[[2]] <- function (x,y){sin(2*pi*x)*sin(2*pi*y)}
f[[3]] <- function (x,y){sin(3*pi*x)*cos(4*pi*y)}
# The lambda to be used
lambda = c(1) #seq(1,20,1)
# The order of FEM
order = 1
################################################################################

n_meshes = length(N)
n_covariates = length(beta)
n_libs_to_test = length(idx_libs_to_test)

dyn.load(paste("fdaPDE_woodbury_whole", .Platform$dynlib.ext, sep = ""))
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
    fun = runif(n_observations[i],-1,1)
    for (j in 1:(length(beta)-1)){
        fun = cbind(fun,f[[j]](locations[[i]][,1],locations[[i]][,2]))
    }
    covariates[[i]] = matrix(fun,
                             nrow = n_observations[i],
                             ncol = n_covariates)
    fun_on_nodes = runif(nrow(mesh$nodes),-1,1)
    for (j in 1:(length(beta)-1)){
        fun_on_nodes = cbind(fun_on_nodes,f[[j]](mesh$nodes[,1],mesh$nodes[,2]))
    }
    covariates_on_nodes[[i]] =
        matrix(fun_on_nodes,
               nrow = nrow(mesh$nodes),
               ncol = n_covariates)
    observations[[i]] = sin(2*pi*locations[[i]][,1])
                        + covariates[[i]] %*% beta
                        + rnorm(n = nrow(locations[[i]]), sd = 0.1)
    observations_on_nodes[[i]] = sin(2*pi*mesh$nodes[,1])
                                 + covariates_on_nodes[[i]] %*% beta
                                 + rnorm(n = nrow(mesh$nodes), sd = 0.1)
    indeces_to_cut= sample(1:length(observations_on_nodes[[i]]), N[i] - n_observations[i], replace=F)
    observations_on_nodes[[i]][indeces_to_cut[]]=NaN
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

            if (k == 1){
                  edf_Eardi1=output_CPP$edf
                  gcv_Eardi1=output_CPP$GCV
            }
            if (k == 2){
                edf_Woodbury1=output_CPP$edf
                gcv_Woodbury1=output_CPP$GCV
            }
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

            if (k == 1){
                  edf_Eardi2=output_CPP$edf
                  gcv_Eardi2=output_CPP$GCV
            }
            if (k == 2){
                edf_Woodbury2=output_CPP$edf
                gcv_Woodbury2=output_CPP$GCV
            }
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

            if (k == 1){
                  edf_Eardi3=output_CPP$edf
                  gcv_Eardi3=output_CPP$GCV
            }
            if (k == 2){
                edf_Woodbury3=output_CPP$edf
                gcv_Woodbury3=output_CPP$GCV
            }
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

            if (k == 1){
                  edf_Eardi4=output_CPP$edf
                  gcv_Eardi4=output_CPP$GCV
            }
            if (k == 2){
                edf_Woodbury4=output_CPP$edf
                gcv_Woodbury4=output_CPP$GCV
            }
        }
    }
}

mpi.quit()
#################################################################################################
#                            ###PLOT###
##TEST1
## to save
##jpeg('test1_edf')
#plot(lambda, edf_Eardi1,ylim=c(4.5,5.5),type="s",col="black")
#points(lambda, edf_Woodbury1,type="s",col="red")
##dev.off()

##TEST2
#plot(lambda, edf_Eardi2,type="s",col="black")
#points(lambda, edf_Woodbury2,type="s",col="red")

##TEST3
#plot(lambda, edf_Eardi3,ylim=c(0.7,1.6),type="s",col="black")
#points(lambda, edf_Woodbury3,type="s",col="red")

#plot(lambda, gcv_Eardi3,type="s",col="black")
#points(lambda, gcv_Woodbury3,type="s",col="red")

##TEST4
#plot(lambda, edf_Eardi4,ylim= c(0.7,1.5),type="s",col="black")
#points(lambda, edf_Woodbury4,type="s",col="red")

#plot(lambda, gcv_Eardi4,ylim=c(0.37,0.385),type="s",col="black")
#points(lambda, gcv_Woodbury4,type="s",col="red")
