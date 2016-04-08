library(fdaPDE)

order = 1
for (n in 1:N) {
    output <- mesh_quadratounitario(10*n)
    mesh<-create.MESH.2D(nodes=output$nodes, order = order)
    FEMbasis = create.FEM.basis(mesh)
    observations = sin(2*pi*mesh$nodes[,1]) + rnorm(n = nrow(mesh$nodes), sd = 0.1)
    lambda = 7
}

mesh<-create.MESH.2D(nodes=rbind(c(0, 0), c(0, 1), c(0.5, 0.5), c(1, 1), c(1, 0)),
                     segments=rbind(c(1, 2), c(2, 3), c(3, 4), c(4, 5), c(5, 1)), order = order)

FEMbasis = create.FEM.basis(mesh)

lambda = c(1,2,3)

locations = rbind(c(0, 0), c(0, 1), c(0.5, 0.5), c(1, 1), c(1, 0))
observations = c(1,2,1,2,1)
data = c(1,2,1,2,1)
covariates = cbind(c(1, 2, 3, 4, 5))
BC = NULL

output_CPP = smooth.FEM.basis(locations  = as.matrix(locations), 
                              observations = data, 
                              FEMbasis = FEMbasis, lambda = lambda, 
                              covariates = covariates, 
                              GCV = TRUE,
                              CPP_CODE = TRUE)
