#' generate_boston_map
#'
#' @param plot A boolean value indicating whether or not to plot the resulting `sf` object of Boston.
#'
#' @return An `sf` object containing geometries for 5 counties in or near Boston.
#' @export
#'
#' @examples
#'
#' generate_boston_map(plot = TRUE)
#'

generate_boston_map <- function(plot = FALSE){

  boston_area <- c("essex","middlesex","norfolk","plymouth","suffolk")
  boston_area_map <- maps::map(database = "county",
                               regions = paste("massachusetts",boston_area,sep = ","),
                               fill = TRUE,
                               plot = FALSE)
  boston_area_sf <- sf::st_as_sf(boston_area_map)

  if(plot){
    p <- ggplot2::ggplot(boston_area_sf)+ggplot2::geom_sf()
    print(p)
  }

  return(boston_area_sf)

}

#' generate_grid
#'
#' Given an sf object, generate an underlying grid sf object. Options exist to rotate the grid.
#'
#' @param sf_object An sf object over which to make the grid
#' @param plot A logical value representing whether or not to plot the result
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#'
#' boston_map <- generate_boston_map()
#' generate_grid(boston_map,plot = TRUE,n=c(10,10),theta = -pi/4)
#'
generate_grid <- function(sf_object,theta = 0,plot = FALSE,...){

  #cellsize = base_cellsize*scaling
  grid <- sf::st_make_grid(sf_object,...)
  #grid <- st_as_sf(grid)

  center <- sf::st_centroid(sf::st_union(grid))
  m = matrix(c(cos(theta),sin(theta),-sin(theta),cos(theta)), 2)
  grid_temp <- sf::st_geometry(grid-center)*m+center
  grid <- sf::st_set_crs(grid_temp,sf::st_crs(grid))

  grid <- sf::st_as_sf(grid)

  if(plot){
   p <-  ggplot2::ggplot(sf_object)+
     ggplot2::geom_sf(lwd = 1)+
     ggplot2::geom_sf(data = grid, alpha = 0)
   print(p)
  }

  return(grid)

}


#' generate_gaussian_random_field
#'
#' @param l Larger l means slower decay in covarinace (higher covariance between nearby points)
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#'
#' boston_map <- generate_boston_map()
#' gaussian_random_field <- generate_gaussian_random_field(boston_map,l = 1,n = c(30,30), theta = -pi/4, plot = TRUE)
#'
generate_gaussian_random_field <- function(sf_object,l = 1,plot = FALSE,...){

  fine_grid <- generate_grid(sf_object,...)
  centroids <- sf::st_centroid(fine_grid)

  n_grid <- dim(fine_grid)[1]

  distances <- sf::st_distance(centroids,centroids)
  standardized_distances <- distances/sd(distances)

  covariances <- exp(-matrix(standardized_distances,nrow = n_grid,ncol = n_grid)/l)
  value <- crossprod(chol(covariances),rnorm(n_grid,0,1))
  fine_grid$value <- value

  if(plot){
    p <-   ggplot2::ggplot()+
      ggplot2::geom_sf(data = sf_object)+
      ggplot2::geom_sf(data = fine_grid, alpha = 0.9,ggplot2::aes(fill = value))
    print(p)
  }


  return(fine_grid)

}

#' compute_intersection_area_percentages
#'
#' Intersect two sf geometries and return the fraction of area of the intersection out of its parent in `sf_object_small`.
#'
#' @param sf_object_small
#' @param sf_object_big
#'
#' @return
#' @export
#'
#' @examples
#'
#' boston_map <- generate_boston_map()
#' fine_grid <- generate_grid(boston_map,n = c(4,4))
#' coarse_grid <- generate_grid(boston_map,n = c(2,2),theta = -pi/4)
#' compute_intersection_area_fractions(fine_grid,coarse_grid,plot = TRUE)
#'
#'
compute_intersection_area_fractions <- function(sf_object_small,sf_object_big,plot = FALSE){

  small_areas <- sf::st_area(sf_object_small)

  intersection <- sf::st_intersection(sf_object_small,sf_object_big)
  intersection_centroids <- sf::st_centroid(intersection)
  intersection_area <- sf::st_area(intersection)
  intersection_area_outer <- crossprod(sf::st_contains(sf_object_small,intersection_centroids,sparse = FALSE),small_areas)

  fraction_of_area <- as.vector(intersection_area/intersection_area_outer) #as.vector removes units
  intersection$fraction_of_area <- fraction_of_area

  if(plot){
    p <- ggplot2::ggplot()+
      ggplot2::geom_sf(data = intersection)+
      ggplot2::geom_sf(data = intersection_centroids)+
      ggplot2::geom_sf_label(data = intersection, ggplot2::aes(label = round(fraction_of_area,2)))
    print(p)
  }

  return(intersection)


}

#' compute_intersection_areas
#'
#' Intersect two sf geometries and return the area of the intersection
#' @param sf_object_small
#' @param sf_object_big
#'
#' @return
#' @export
#'
#' @examples
#'
#' boston_map <- generate_boston_map()
#' fine_grid <- generate_grid(boston_map,n = c(4,4))
#' coarse_grid <- generate_grid(boston_map,n = c(2,2),theta = -pi/4)
#' compute_intersection_areas(fine_grid,coarse_grid,plot = TRUE)
#'
#'
compute_intersection_areas <- function(sf_object_small,sf_object_big,plot = FALSE){

  small_areas <- sf::st_area(sf_object_small)

  intersection <- sf::st_intersection(sf_object_small,sf_object_big)
  intersection_centroids <- sf::st_centroid(intersection)
  intersection_area <- sf::st_area(intersection)

  intersection$area <- as.vector(intersection_area) #as.vector to remove the units

  if(plot){
    p <- ggplot2::ggplot()+
      ggplot2::geom_sf(data = intersection)+
      ggplot2::geom_sf(data = intersection_centroids)+
      ggplot2::geom_sf_label(data = intersection, ggplot2::aes(label = round(area,2)))
    print(p)
  }

  return(intersection)


}

#' compute_distance_covariance_matrix
#'
#' Given an sf object, compute a covariance matrix based on distances between
#' centroids using an exponential kernel.
#'
#' @param sf_object Object of class `sf`
#' @param l Smoothness parameter: bigger `l` represents less correlation between
#' neighboring points.
#'
#' @return a matrix representing the covariance between points at different distances.
#' @export
#'
#' @examples
compute_distance_covariance_matrix <- function(sf_object,l){

  n <- dim(sf_object)[1]

  centroids <- st_centroid(sf_object)
  distances <- sf::st_distance(centroids,centroids)
  standardized_distances <- distances/sd(distances)

  covariances <- exp(-matrix(standardized_distances,nrow = n,ncol = n)/l)

  return(covariances)

}

#' generate_coarse_grid_value
#'
#' Given an sf object associated with a `value` and a second sf object, compute a value
#' for the second sf object based on the first objects' values over intersections.
#'
#' @param fine_grid A sf object thought to have smaller geometries; requires a `value` field
#' @param coarse_grid A sf object thought to have larger geometries
#'
#' @return An sf object equal to `coarse_grid` with a `value` field attached.
#' @export
#'
#' @examples
generate_coarse_grid_value <- function(fine_grid,coarse_grid){

  intersection <- st_intersection(coarse_grid,fine_grid)

  intersection_centroids <- sf::st_centroid(intersection)
  fine_grid_containment_matrix <- sf::st_contains(fine_grid,intersection_centroids,sparse = FALSE)
  coarse_grid_containment_matrix <- sf::st_contains(coarse_grid,intersection_centroids,sparse = FALSE)

  coarse_grid_areas <- as.vector(st_area(coarse_grid))

  intersection_areas <- compute_intersection_areas(fine_grid,coarse_grid,plot = FALSE)
  fine_area_mat <- apply(fine_grid_containment_matrix,1,function(x){intersection_areas$area*x})
  coarse_area_weight_mat <- apply(coarse_grid_containment_matrix,2,function(x) 1/coarse_grid_areas*x)

  A <- coarse_area_weight_mat %*% fine_area_mat

  coarse_grid_value <- A %*% fine_grid$value

  coarse_grid$value <- coarse_grid_value

  return(coarse_grid)

}

#' estimate_fine_grid
#'
#' Uses the values associated with geometries in a coarse grid to estimate values
#' for geometries in a finer grid. This is done using regularization that assumes
#' a gaussian prior with a given covariance structure.#'
#'
#' @param fine_grid A finer grid which we would like to estimate values for
#' @param coarse_grid A coarser grid whose values are known and in a `value` field
#' @param covariance_matrix A covariance matrix with one row/col per element of fine_grid
#' @param lambda A scaling parameter controlling the regularization
#'
#' @return
#' @export
#'
#' @examples
estimate_fine_grid <- function(fine_grid,coarse_grid,covariance_matrix,lambda){

  intersection <- st_intersection(coarse_grid,fine_grid)

  intersection_centroids <- sf::st_centroid(intersection)
  fine_grid_containment_matrix <- sf::st_contains(fine_grid,intersection_centroids,sparse = FALSE)
  coarse_grid_containment_matrix <- sf::st_contains(coarse_grid,intersection_centroids,sparse = FALSE)

  coarse_grid_areas <- as.vector(st_area(coarse_grid))

  intersection_areas <- compute_intersection_areas(fine_grid,coarse_grid,plot = FALSE)
  fine_area_mat <- apply(fine_grid_containment_matrix,1,function(x){intersection_areas$area*x})
  coarse_area_weight_mat <- apply(coarse_grid_containment_matrix,2,function(x) 1/coarse_grid_areas*x)

  A <- coarse_area_weight_mat %*% fine_area_mat

  estimated_fine_values <- solve(crossprod(A)+lambda^2*crossprod(solve(covariance_matrix)))%*%crossprod(A,coarse_grid$value)

  fine_grid$estimated_value <- estimated_fine_values

  return(fine_grid)


}

