% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/grid_simulation.R
\name{compute_intersection_areas}
\alias{compute_intersection_areas}
\title{compute_intersection_areas}
\usage{
compute_intersection_areas(sf_object_small, sf_object_big, plot = FALSE)
}
\arguments{
\item{sf_object_big}{}
}
\description{
Intersect two sf geometries and return the area of the intersection
}
\examples{

boston_map <- generate_boston_map()
fine_grid <- generate_grid(boston_map,n = c(4,4))
coarse_grid <- generate_grid(boston_map,n = c(2,2),theta = -pi/4)
compute_intersection_areas(fine_grid,coarse_grid,plot = TRUE)


}
