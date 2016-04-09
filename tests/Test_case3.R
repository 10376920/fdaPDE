#TEST caso 3: CON COVARIATE LOCATION ON NODES

librerie=c("fdaPDE_original","fdaPDE_schur")

for( k in 1:2){
	print("Test: covariates, locations on nodes")
	cat("Library = ", librerie[k])
	cat("\n")
	library = librerie[k]
	dyn.load(paste(library, .Platform$dynlib.ext, sep = ""))
	pathnames <- list.files(pattern="[.]R$", path="./R/", full.names=TRUE);
	sapply(pathnames, FUN=source);
	source("tests/mesh.R")

	order = 1
	lambda = 1

	for( n in 1:7){
		output <- mesh_quadratounitario(n*10)
		mesh<-create.MESH.2D(nodes=output$nodes, order = order)
		observations_on_nodes= sin(2*pi*mesh$nodes[,1]) + rnorm(n = nrow(mesh$nodes), sd = 0.1)
		#locations <- cbind(cbind(runif(n*100*n,0,1)),cbind(runif(n*100*n,0,1)))
		#observations_not_on_nodes = sin(2*pi*locations[,1]) + rnorm(n = nrow(locations), sd = 0.1)
		covariates = cbind(cbind(runif(n*100*n,-100,100)),cbind(runif(n*100*n,-100,100)))

		FEMbasis = create.FEM.basis(mesh)

		cat("n = ", n*10)

		output_CPP_00 =
		    smooth.FEM.basis(observations = observations_on_nodes, 
		                     FEMbasis = FEMbasis,
		                     lambda = lambda, 
		                  	 covariates = covariates,
		                     GCV = TRUE,
		                     CPP_CODE = TRUE)
	}
}