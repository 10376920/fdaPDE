library = "../fdaPDE_schur"

dyn.load(paste(library, .Platform$dynlib.ext, sep = ""))
pathnames <- list.files(pattern="[.]R$", path="../R/", full.names=TRUE);
sapply(pathnames, FUN=source);

source("mesh.R")

order = 1
output <- mesh_quadratounitario(60)
mesh<-create.MESH.2D(nodes=output$nodes,
                     segments=output$segments, order = order)

FEMbasis = create.FEM.basis(mesh)

observations = sin(2*pi*mesh$nodes[,1]) + rnorm(n = nrow(mesh$nodes), sd = 0.1)

lambda = 7

print("Test: no covariates, locations on nodes")
#output_R_00 =
#    smooth.FEM.basis(observations = observations, 
#                     FEMbasis = FEMbasis,
#                     lambda = lambda, 
#                     GCV = TRUE,
#                     CPP_CODE = FALSE)


output_CPP_00 =
    smooth.FEM.basis(observations = observations, 
                     FEMbasis = FEMbasis,
                     lambda = lambda, 
                     GCV = TRUE,
                     CPP_CODE = TRUE)
