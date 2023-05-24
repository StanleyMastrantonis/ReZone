#' Calculate pairwise distances between two matrices
#'
#' This function calculates the pairwise Euclidean distances between two numeric matrices.
#' It compares each row in the first matrix (\code{m1}) with each row in the second matrix (\code{m2}),
#' and returns a matrix of distances.
#'
#' @param m1 A numeric matrix.
#' @param m2 A numeric matrix.
#'
#' @return A numeric matrix containing the pairwise distances between rows of \code{m1} and \code{m2}.
#' The resulting matrix has dimensions \code{nrow(m1)} x \code{nrow(m2)}.
#'
#' @details The number of columns in \code{m1} must be equal to the number of columns in \code{m2}.
#' If the dimensions are not compatible, an error is thrown.
#'
#' @examples
#' m1 <- matrix(c(1, 2, 3, 4), nrow = 2)
#' m2 <- matrix(c(5, 6, 7, 8), nrow = 2)
#' crossdist(m1, m2)
#'
#' @export
cppFunction('NumericMatrix crossdist(NumericMatrix m1, NumericMatrix m2) {
  int nrow1 = m1.nrow();
  int nrow2 = m2.nrow();
  int ncol = m1.ncol();
  if (ncol != m2.ncol()) {
    throw std::runtime_error("Incompatible number of dimensions");
  }
  NumericMatrix out(nrow1, nrow2);
  for (int r1 = 0; r1 < nrow1; r1++) {
    for (int r2 = 0; r2 < nrow2; r2++) {
      double total = 0;
      for (int c12 = 0; c12 < ncol; c12++) {
        total += pow(m1(r1, c12) - m2(r2, c12), 2);
      }
      out(r1,r2) = sqrt(total);
    }
  }
  return out;
}')

#' Calculate pairwise distances within a matrix
#'
#' This function calculates the pairwise Euclidean distances within a numeric matrix.
#' It compares each row in the matrix (\code{m1}) with every other row,
#' and returns a vector of distances.
#'
#' @param m1 A numeric matrix.
#'
#' @return A numeric vector containing the pairwise distances between rows of \code{m1}.
#' The length of the resulting vector is \code{nrow(m1)}.
#'
#' @examples
#' m1 <- matrix(c(1, 2, 3, 4), nrow = 2)
#' pairdist(m1)
#'
#' @export
cppFunction('NumericVector pairdist(NumericMatrix m1, NumericMatrix m2) {
  int nrow1 = m1.nrow();
  int ncol = m1.ncol();
  NumericVector out(nrow1);
  for (int r1 =0 ; r1 < nrow1; r1++) {
      double total = 0;
      for (int c12 = 0; c12 < ncol; c12++) {
        total += pow(m1(r1, c12) - m2(r1, c12), 2);
      }
      out(r1) = sqrt(total);
    
  }
  return out;
}')


#' Find the minimum distance between two matrices
#'
#' This function calculates the minimum pairwise distance between rows of two numeric matrices.
#' It uses the \code{crossdist} function to compute the pairwise distances, and returns a vector
#' containing the minimum distances for each column.
#'
#' @param x A numeric matrix.
#' @param y A numeric matrix.
#'
#' @return A numeric vector containing the minimum distances between rows of \code{x} and \code{y},
#' for each column.
#'
#' @examples
#' x <- matrix(c(1, 2, 3, 4), nrow = 2)
#' y <- matrix(c(5, 6, 7, 8), nrow = 2)
#' mindist_func(x, y)
#'
#' @export
mindist_func <- function(x, y) {
  d_list <- list()
  cd <- crossdist(as.matrix(x), as.matrix(y))
  for (i in 1:ncol(cd)) {
    d_list[[i]] <- cd[, i]
  }
  mind <- unlist(lapply(d_list, min))
  mind
}

#' Find the maximum distance between two matrices
#'
#' This function calculates the maximum pairwise distance between rows of two numeric matrices.
#' It uses the \code{crossdist} function to compute the pairwise distances, and returns a vector
#' containing the maximum distances for each column.
#'
#' @param x A numeric matrix.
#' @param y A numeric matrix.
#'
#' @return A numeric vector containing the maximum distances between rows of \code{x} and \code{y},
#' for each column.
#'
#' @examples
#' x <- matrix(c(1, 2, 3, 4), nrow = 2)
#' y <- matrix(c(5, 6, 7, 8), nrow = 2)
#' maxdist_func(x, y)
#'
#' @export
maxdist_func <- function(x, y) {
  d_list <- list()
  cd <- crossdist(as.matrix(x), as.matrix(y))
  for (i in 1:ncol(cd)) {
    d_list[[i]] <- cd[, i]
  }
  mind <- unlist(lapply(d_list, max))
  mind
}

#' Find the index of the minimum distance between two matrices
#'
#' This function calculates the index of the row in the second matrix (\code{y})
#' that has the minimum distance to each row in the first matrix (\code{x}).
#' It uses the \code{crossdist} function to compute the pairwise distances,
#' and returns a vector containing the indices of the minimum distances for each column.
#'
#' @param x A numeric matrix.
#' @param y A numeric matrix.
#'
#' @return A numeric vector containing the indices of the rows in \code{y} with the minimum distances
#' to each row in \code{x}, for each column.
#'
#' @examples
#' x <- matrix(c(1, 2, 3, 4), nrow = 2)
#' y <- matrix(c(5, 6, 7, 8), nrow = 2)
#' whichmindist_func(x, y)
#'
#' @export
whichmindist_func <- function(x, y) {
  d_list <- list()
  cd <- crossdist(as.matrix(x), as.matrix(y))
  for (i in 1:ncol(cd)) {
    d_list[[i]] <- cd[, i]
  }
  which_min <- unlist(lapply(d_list, which.min))
  which_min
}

#' Find the indices of distances within a threshold
#'
#' This function calculates the indices of the rows in the second matrix (\code{y})
#' that have distances within a specified threshold to each row in the first matrix (\code{x}).
#' It uses the \code{crossdist} function to compute the pairwise distances,
#' and returns a vector containing the indices of the rows in \code{y} that satisfy the distance criterion
#' for each column.
#'
#' @param x A numeric matrix.
#' @param y A numeric matrix.
#' @param dist The distance threshold.
#'
#' @return A numeric vector containing the indices of the rows in \code{y} that have distances
#' within the threshold (\code{dist}) to each row in \code{x}, for each column.
#'
#' @examples
#' x <- matrix(c(1, 2, 3, 4), nrow = 2)
#' y <- matrix(c(5, 6, 7, 8), nrow = 2)
#' whichdist_func(x, y, dist = 5)
#'
#' @export
whichdist_func <- function(x, y, dist) {
  d_list <- list()
  cd <- crossdist(as.matrix(x), as.matrix(y))
  for (i in 1:ncol(cd)) {
    d_list[[i]] <- cd[, i]
  }
  which <- unlist(lapply(d_list, function(x) { which(x <= dist) }))
  which
}

#' Find the indices of distances within a threshold (as a list)
#'
#' This function calculates the indices of the rows in the second matrix (\code{y})
#' that have distances within a specified threshold to each row in the first matrix (\code{x}).
#' It uses the \code{crossdist} function to compute the pairwise distances,
#' and returns a list containing the indices of the rows in \code{y} that satisfy the distance criterion
#' for each column.
#'
#' @param x A numeric matrix.
#' @param y A numeric matrix.
#' @param dist The distance threshold.
#'
#' @return A list containing the indices of the rows in \code{y} that have distances
#' within the threshold (\code{dist}) to each row in \code{x}, for each column.
#'
#' @examples
#' x <- matrix(c(1, 2, 3, 4), nrow = 2)
#' y <- matrix(c(5, 6, 7, 8), nrow = 2)
#' whichdist_func_list(x, y, dist = 5)
#'
#' @export
whichdist_func_list <- function(x, y, dist) {
  d_list <- list()
  cd <- crossdist(as.matrix(x), as.matrix(y))
  for (i in 1:ncol(cd)) {
    d_list[[i]] <- cd[, i]
  }
  which <- lapply(d_list, function(x) { which(x <= dist) })
  which
}

#' Scale values between 0 and 1
#'
#' This function scales the values of a numeric vector between 0 and 1 using min-max normalization.
#' The scaled values are calculated as \code{(x - min(x)) / (max(x) - min(x))}.
#'
#' @param x A numeric vector.
#'
#' @return A numeric vector with values scaled between 0 and 1.
#'
#' @examples
#' x <- c(1, 2, 3, 4)
#' scale_func(x)
#'
#' @export
scale_func <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

#' Scale values within a specified range
#'
#' This function scales the values of a numeric vector to a specified range \code{(a, b)}.
#' The scaled values are calculated as \code{(b - a) * ((x - min(x)) / (max(x) - min(x))) + a}.
#'
#' @param x A numeric vector.
#' @param a The lower bound of the desired range.
#' @param b The upper bound of the desired range.
#'
#' @return A numeric vector with values scaled to the specified range.
#'
#' @examples
#' x <- c(1, 2, 3, 4)
#' scale_func_range(x, 0, 10)
#'
#' @export
scale_func_range <- function(x, a, b) {
  (b - a) * ((x - min(x)) / (max(x) - min(x))) + a
}

#' Scale values within a raster to the range [0, 1]
#'
#' This function scales the values of a raster to the range [0, 1] using min-max normalization.
#' The minimum and maximum values of the raster are obtained using the \code{cellStats} function.
#' The scaled values are calculated as \code{(x - min(x)) / (max(x) - min(x))}.
#'
#' @param x A raster object.
#'
#' @return A raster object with values scaled to the range [0, 1].
#'
#' @examples
#' r <- raster(matrix(c(1, 2, 3, 4), nrow = 2))
#' scale_func_ras(r)
#'
#' @export
scale_func_ras <- function(x) {
  (x - cellStats(x, stat = 'min')) / (cellStats(x, stat = 'max') - cellStats(x, stat = 'min'))
}

#' Scale values within a raster to a specified range
#'
#' This function scales the values of a raster to a specified range \code{(a, b)}.
#' The minimum and maximum values of the raster are obtained using the \code{cellStats} function.
#' The scaled values are calculated as \code{(b - a) * ((x - min(x)) / (max(x) - min(x))) + a}.
#'
#' @param x A raster object.
#' @param a The lower bound of the desired range.
#' @param b The upper bound of the desired range.
#'
#' @return A raster object with values scaled to the specified range.
#'
#' @examples
#' r <- raster(matrix(c(1, 2, 3, 4), nrow = 2))
#' scale_func_range_ras(r, 0, 10)
#'
#' @export
scale_func_range_ras <- function(x, a, b) {
  (b - a) * (x - cellStats(x, stat = 'min')) / (cellStats(x, stat = 'max') - cellStats(x, stat = 'min')) + a
}

#' Sigmoid function
#'
#' This function calculates the sigmoid function for a given input vector \code{x}.
#' The sigmoid function is defined as \code{1 / (1 + exp((x - d) / coef))}.
#'
#' @param x A numeric vector.
#' @param coef Coefficient of the sigmoid function.
#' @param d Value at which the sigmoid function changes rapidly.
#'
#' @return A numeric vector with sigmoid function values.
#'
#' @examples
#' x <- c(1, 2, 3, 4)
#' sigmoid_func(x, 1, 2)
#'
#' @export
sigmoid_func <- function(x, coef, d) {
  1 / (1 + exp((x - d) / coef))
}

#' Sigmoid function with maximum value adjustment
#'
#' This function calculates the sigmoid function for a given input vector \code{x},
#' but with an additional maximum value adjustment \code{max.add}.
#' The sigmoid function is defined as \code{max.add / (1 + exp((x - d) / coef))}.
#'
#' @param x A numeric vector.
#' @param coef Coefficient of the sigmoid function.
#' @param d Value at which the sigmoid function changes rapidly.
#' @param max.add Maximum value adjustment.
#'
#' @return A numeric vector with sigmoid function values.
#'
#' @examples
#' x <- c(1, 2, 3, 4)
#' sigmoid_func_2(x, 1, 2, 10)
#'
#' @export
sigmoid_func_2 <- function(x, coef, d, max.add) {
  max.add / (1 + exp((x - d) / coef))
}

#' Calculate the geometric mean of positive values
#'
#' This function calculates the geometric mean of a numeric vector, considering only the positive values.
#' The function ignores zero and negative values while computing the geometric mean.
#' The geometric mean is calculated as \code{prod(x[x > 0])^(1 / length(x[x > 0]))}.
#'
#' @param x A numeric vector.
#'
#' @return The geometric mean of the positive values in \code{x}.
#'
#' @examples
#' x <- c(1, 2, 3, 4, -5)
#' geo_mean(x)
#'
#' @export
geo_mean <- function(x) {
  prod(x[x > 0])^(1 / length(x[x > 0]))
}

#' Calculate the sum of positive values
#'
#' This function calculates the sum of a numeric vector, considering only the positive values.
#' If the sum is zero or \code{NA}, it returns 0. Otherwise, it returns the sum.
#'
#' @param x A numeric vector.
#'
#' @return The sum of the positive values in \code{x}, or 0 if the sum is zero or \code{NA}.
#'
#' @examples
#' x <- c(1, 2, -3, 4)
#' sum_func(x)
#'
#' @export
sum_func <- function(x) {
  if (is.na(sum(x)) || sum(x) <= 0) {
    return(0)
  } else {
    x_sub <- x[x > 0]
    y <- sum(x_sub)
    sum <- ifelse(y > 1, 1, y)
    return(sum)
  }
}

#' Calculate the mean of positive values
#'
#' This function calculates the mean of a numeric vector, considering only the positive values.
#' If the mean is zero or \code{NA}, it returns 0. Otherwise, it returns the mean.
#'
#' @param x A numeric vector.
#'
#' @return The mean of the positive values in \code{x}, or 0 if the mean is zero or \code{NA}.
#'
#' @examples
#' x <- c(1, 2, -3, 4)
#' mean_func(x)
#'
#' @export
mean_func <- function(x) {
  if (is.na(sum(x)) || sum(x) <= 0) {
    return(0)
  } else {
    x_sub <- x[x > 0]
    y <- mean(x_sub)
    return(y)
  }
}

#' Find the index of the value closest to a given sequence
#'
#' This function finds the index of the value in a vector \code{vector} that is closest to a given sequence \code{seq}.
#' If multiple values are equally close, it returns the index of the first occurrence.
#'
#' @param seq A numeric sequence.
#' @param vector A numeric vector.
#'
#' @return The index of the value in \code{vector} closest to \code{seq}.
#'
#' @examples
#' seq <- c(0.2, 0.4, 0.6, 0.8)
#' vector <- c(0.15, 0.35, 0.65, 0.85)
#' find_func(seq, vector)
#'
#' @export
find_func <- function(seq, vector) {
  y <- which(abs(seq - vector) == min(abs(seq - vector)))
  y[1]
}

#' Find all indices of values closest to a given sequence
#'
#' This function finds all the indices of the values in a vector \code{vector} that are closest to a given sequence \code{seq}.
#' If multiple values are equally close, it returns the indices of all occurrences.
#'
#' @param seq A numeric sequence.
#' @param vector A numeric vector.
#'
#' @return A numeric vector with the indices of the values in \code{vector} closest to \code{seq}.
#'
#' @examples
#' seq <- c(0.2, 0.4, 0.6, 0.8)
#' vector <- c(0.15, 0.35, 0.65, 0.85)
#' find_func_all(seq, vector)
#'
#' @export
find_func_all <- function(seq, vector) {
  which(abs(seq - vector) == min(abs(seq - vector)))
}




#' Calculate inverse distance weighting (IDW) with kernel adjustment
#'
#' This function calculates the inverse distance weighting (IDW) with a kernel adjustment. It takes a numeric vector of weights and a numeric vector of distances as arguments and returns a numeric value representing the IDW calculation.
#' 
#' The IDW calculation is adjusted based on the length of the weight vector. If the length exceeds a specified maximum kernel size, it is truncated to the maximum kernel size. The IDW value is then multiplied by the ratio of the adjusted length to the maximum kernel size.
#'
#' @param weight A numeric vector of weights.
#' @param distance A numeric vector of distances.
#'
#' @return A numeric value representing the IDW calculation with kernel adjustment.
#'
#' @examples
#' weight <- c(1, 2, 3, 4)
#' distance <- c(0.5, 0.7, 0.9, 1.2)
#' idw_func_k(weight, distance)
#'
#' @export

idw_func_k <- function(weight, distance){
  if(length(weight) < 1){
    return(0)
  } else {
    max_kern <- 100  # Maximum kernel size
    x_len <- ifelse(length(weight) > max_kern, max_kern, length(weight))  # Set the length of the kernel
    xw <- weight / distance  # Calculate weighted values
    xw_x <- sum(xw)  # Sum of the weighted values
    pd <- sum(1 / distance)  # Sum of the inverse distances
    dw <- (xw_x / pd) * (x_len / max_kern)  # Calculate the IDW value with kernel adjustment
    return(dw)
  }
}

#' Calculate the sum of distances
#'
#' This function takes a numeric vector of weights and a numeric vector of distances as arguments and returns a numeric value representing the sum of distances.
#' 
#' If the length of the weight vector is less than 1, it returns 0. Otherwise, it calculates the weighted values using the formula: 1 - (weight / (distance + weight)). It then calculates the product of the weighted values using the `prod()` function and subtracts it from 1 to obtain the sum of distances.
#'
#' @param weight A numeric vector of weights.
#' @param distance A numeric vector of distances.
#'
#' @return A numeric value representing the sum of distances.
#'
#' @examples
#' weight <- c(1, 2, 3, 4)
#' distance <- c(0.5, 0.7, 0.9, 1.2)
#' d_sum_func(weight, distance)
#'
#' @export

d_sum_func <- function(weight, distance){
  if(length(weight) < 1){
    return(0)
  } else {
    xw <- 1 - (weight / (distance + weight))  # Calculate weighted values
    xw_x <- 1 - prod(xw)  # Product of the weighted values
    return(xw_x)
  }
}

#' Calculate the sum of distances (alternative implementation)
#'
#' This function takes a numeric vector of weights and a numeric vector of distances as arguments and returns a numeric value representing the sum of distances.
#' 
#' If the length of the weight vector is less than 1, it returns 0. Otherwise, it calculates the weighted values using the formula: 1 - (1 - weight). It then calculates the product of the weighted values using the `prod()` function and subtracts it from 1 to obtain the sum of distances.
#'
#' @param weight A numeric vector of weights.
#' @param distance A numeric vector of distances.
#'
#' @return A numeric value representing the sum of distances.
#'
#' @examples
#' weight <- c(1, 2, 3, 4)
#' distance <- c(0.5, 0.7, 0.9, 1.2)
#' d_sum_func_2(weight, distance)
#'
#' @export

d_sum_func_2 <- function(weight, distance){
  if(length(weight) < 1){
    return(0)
  } else {
    xw <- 1 - (1 - weight)  # Calculate weighted values
    xw_x <- 1 - prod(xw)  # Product of the weighted values
    return(xw_x)
  }
}




#' Calculate inverse distance weighting (IDW)
#'
#' This function takes a numeric vector of weights and a numeric vector of distances as arguments and returns a numeric value representing the IDW calculation.
#' 
#' If the length of the weight vector is less than 1, it returns 0. Otherwise, it calculates the weighted values using the formula: weight / distance. It then calculates the sum of the weighted values and the sum of the inverse distances using the `sum()` function. Finally, it divides the sum of the weighted values by the sum of the inverse distances to obtain the IDW value.
#'
#' @param weight A numeric vector of weights.
#' @param distance A numeric vector of distances.
#'
#' @return A numeric value representing the IDW calculation.
#'
#' @examples
#' weight <- c(1, 2, 3, 4)
#' distance <- c(0.5, 0.7, 0.9, 1.2)
#' idw_func(weight, distance)
#'
#' @export

idw_func <- function(weight, distance){
  if(length(weight) < 1){
    return(0)
  } else {
    xw <- weight / distance  # Calculate weighted values
    xw_x <- sum(xw)  # Sum of the weighted values
    pd <- sum(1 / distance)  # Sum of the inverse distances
    dw <- xw_x / pd  # Calculate the IDW value
    return(dw)
  }
}


#' Calculate the area of a patch
#'
#' This function takes a numeric vector representing a patch as an argument and returns a numeric value representing the area of the patch.
#' 
#' It calculates the sum of positive values in the patch using the `sum()` function. It also calculates the length of the patch and divides the sum of positive values by the length to obtain the proportion of positive values.
#'
#' @param patch A numeric vector representing a patch.
#'
#' @return A numeric value representing the area of the patch.
#'
#' @examples
#' patch <- c(0.1, 0.2, 0, 0.4, 0.3)
#' area_func(patch)
#'
#' @export

area_func <- function(patch){
  a_p <- sum(patch[patch > 0])  # Sum of positive values in the patch
  area <- length(patch)  # Length of the patch
  prop <- a_p / area  # Calculate the proportion of positive values
  return(prop)
}



#' Get topologies based on area and resolution
#'
#' This function takes an area and a resolution as arguments and returns a grid of topologies based on the specified area and resolution.
#' 
#' It calculates the bounding box using the `bbox()` function and rounds it to the nearest multiple of the resolution. It then creates a `GridTopology` object using the rounded bounding box, resolution, and dimensions of the cells. The dimensions of the cells are calculated based on the difference between the coordinates of the bounding box and the resolution. 
#'
#' Next, it creates a `SpatialPoints` object using the `GridTopology` object. It uses the `over()` function to extract the topology information from the `study_area` object for each point in the `SpatialPoints` object. It then creates a data frame called `topo_coords_study` by combining the coordinates of the points and the topology information. Any rows with missing topology information are removed from the data frame.
#'
#' Finally, it converts the `topo_coords_study` data frame into a matrix called `grid` and returns it.
#'
#' @param area A numeric value representing the area.
#' @param res A numeric value representing the resolution.
#'
#' @return A matrix representing the grid of topologies based on the specified area and resolution.
#'
#' @examples
#' area <- 100
#' res <- 0.1
#' get_topologies(area, res)
#'
#' @export

get_topologies <- function(area, res){
    bb = res*round(bbox(study_area)/res)
    gt = GridTopology(cellcentre.offset = bb[,1], 
                      cellsize = c(res, res),
                      cells.dim = (c(diff(bb[1,]), diff(bb[2,]))/res) + 1)
    
    sp = SpatialPoints(gt)
    sp_topo_study = over(sp, study_area)
    topo_coords_study = data.frame(cbind(coordinates(sp), sp_topo_study)[!is.na(sp_topo_study),1:2])
    
    grid = as.matrix(topo_coords_study, rownames.force = FALSE)
    return(grid)
}


#' Create a data frame with sigmoidal values
#'
#' This function takes a minimum value and a maximum value as arguments and creates a data frame with sigmoidal values calculated using the `sigmoid_func_2` function.
#' 
#' It assigns the maximum value to `maxdist_base` and the minimum value to `mindist_base`. It then generates a sequence of values from `mindist_base` to `maxdist_base` with an increment of 0.1 using the `seq()` function. The `scale_func` function is applied to the sequence to scale it. The scaled sequence is then passed to the `sigmoid_func_2` function along with parameters `0.08`, `0.5`, and `1` to calculate the sigmoidal values.
#'
#' Next, it creates a data frame called `value_df` by combining the sequence and the sigmoidal values. The column names of the data frame are set to 'value_seq' and 'value'.
#'
#' Finally, it returns the `value_df` data frame.
#'
#' @param min A numeric value representing the minimum value.
#' @param max A numeric value representing the maximum value.
#'
#' @return A data frame with sigmoidal values.
#'
#' @examples
#' min_val <- 0
#' max_val <- 10
#' sdf(min_val, max_val)
#'
#' @export

sdf <- function(min, max){
  maxdist_base = max
  mindist_base = min
  value = sigmoid_func_2(scale_func(seq(mindist_base,maxdist_base,0.1)), 0.08,0.5,1)
  value_df = data.frame(seq(mindist_base,maxdist_base,0.1),value)
  colnames(value_df) = c('value_seq','value')
  return(value_df)
}




#' Calculate the maximum kernel size
#'
#' This function takes two resources as input, `resource1` and `resource2`, and calculates the maximum kernel size based on the distance between them.
#'
#' It creates a data frame called `kern_d_list` by computing the pairwise distances between `resource1` and `resource2` using the `crossdist()` function. The distances are stored in the 'distance' column of the data frame.
#' The 'index' column is created by repeating a sequence from 1 to the number of rows in `resource1` for the number of rows in `resource2`.
#'
#' Next, the 'inside' column is created by checking if the distances in the 'distance' column are less than or equal to a variable named `forage`.
#' The 'inside' column is used to filter the rows in `kern_d_list` to include only the distances that are within the range of `forage`.
#'
#' The variable `max_kern` is computed by taking the ceiling of the median of the table of frequencies for the 'index' column in `kern_list`.
#'
#' Finally, the function returns the value of `max_kern`.
#'
#' @param resource1 A data frame or matrix representing the first resource.
#' @param resource2 A data frame or matrix representing the second resource.
#'
#' @return A numeric value representing the maximum kernel size.
#'
#' @examples
#' res1 <- data.frame(x = c(1, 2, 3), y = c(4, 5, 6))
#' res2 <- data.frame(x = c(7, 8, 9), y = c(10, 11, 12))
#' max_kern_func(res1, res2)
#'
#' @export

max_kern_func <- function(resource1, resource2){
  kern_d_list = data.frame(matrix(crossdist(resource1, resource2), ncol =1))
  kern_d_list$index = rep(seq(1,nrow(resource1),1), times = nrow(resource2))
  colnames(kern_d_list) =  c('distance','index' )
  kern_d_list$inside  = kern_d_list$distance <= forage
  kern_list = kern_d_list[kern_d_list$inside == TRUE,]
  max_kern = ceiling(median(table(kern_list$index)))
  return(max_kern)
}




#' Calculate resource values
#'
#' This function takes three inputs: `resource1`, `grid`, and `ac`. It calculates the resource values based on the given inputs.
#'
#' It starts by converting the `ac` raster object to a matrix using the `rasterToPoints()` function and assigns column names to the resulting matrix.
#' The `ac_coords` matrix is created by extracting the x and y coordinates from `ac_mat`.
#'
#' The `dist_list` is computed by calling the `whichdist_func_list()` function with `ac_coords`, `resource1`, and `forage` as arguments.
#' The `match_list` data frame is created by iterating over the length of `dist_list` and combining the relevant values from `resource1`, `ac_mat`, and `dist_list`.
#' Column names are assigned to `match_list` to indicate the corresponding variables.
#' The `ac_x` and `ac_y` columns are populated with the x and y coordinates from `ac_coords` based on the `ac_ind` values in `match_list`.
#' The `dist` column is computed by calling the `pairdist()` function on selected columns of `match_list`.
#' The unique values of `index` in `match_list` are stored in `ind_loop`.
#'
#' The `idw_val`, `x`, and `y` vectors are initialized as empty.
#' A loop is executed over the indices in `ind_loop`, and for each index, a subset of `match_list` is created based on the index value and non-zero distances.
#' The `idw_val` is computed using the `idw_func()` function on the `ac` and `dist` values in the subset, and the first values of `x` and `y` are stored.
#' The resulting `x`, `y`, and `idw_val` vectors are combined into the `df_agg` data frame, with appropriate column names.
#'
#' Finally, the function returns the `ac` values from `df_agg`.
#'
#' @param resource1 A data frame or matrix representing the resource.
#' @param grid A grid object representing the spatial grid.
#' @param ac A raster object representing the variable 'ac'.
#'
#' @return A numeric vector representing the resource values.
#'
#' @examples
#' res1 <- data.frame(x = c(1, 2, 3), y = c(4, 5, 6))
#' grid <- GridTopology(cellcentre.offset = c(0, 0), cellsize = c(1, 1), cells.dim = c(3, 3))
#' ac <- raster(matrix(1:9, nrow = 3))
#' res_value_function(res1, grid, ac)
#'
#' @export

res_value_function <- function(resource1, grid, ac){
  
  ac_mat = rasterToPoints(ac)
  colnames(ac_mat) = c('x','y','val')
  ac_coords = as.matrix(ac_mat[,1:2])
  dist_list = whichdist_func_list(ac_coords, resource1, forage)
  match_list = data.frame(do.call(rbind,lapply(1:length(dist_list), function(xi){
    cbind(resource1[xi,1],resource1[xi,2],xi,ac_mat[dist_list[[xi]],3])
  }
  )))
  
  match_list$ac_ind = unlist(dist_list)
  colnames(match_list) = c('x','y','index','ac','ac_ind')
  match_list$ac_x = ac_coords[match_list$ac_ind,1]
  match_list$ac_y = ac_coords[match_list$ac_ind,2]
  match_list$dist = pairdist(as.matrix(match_list[,1:2]), as.matrix(match_list[,6:7]))
  ind_loop = unique(match_list$index)
  idw_val = NULL
  x = NULL
  y = NULL
  
  for (i in 1:length(ind_loop)){
    sub = match_list[match_list$index == ind_loop[i] & match_list$dist > 0,]
    idw_val[i] = ifelse(is.null(dim(sub)),0,idw_func(sub$ac,sub$dist))
    x[i] = sub[1,1]
    y[i] = sub[1,2]
  }
  
  df_agg = data.frame(x,y,idw_val)
  colnames(df_agg) =  c('x','y','ac')
  return(df_agg$ac)
  
}


#' Calculate resource values based on area
#'
#' This function takes three inputs: `resource1`, `grid`, and `ac`. It calculates the resource values based on the given inputs and the area of the patch.
#'
#' It starts by converting the `ac` raster object to a matrix using the `rasterToPoints()` function and assigns column names to the resulting matrix.
#' The `ac_coords` matrix is created by extracting the x and y coordinates from `ac_mat`.
#'
#' The `dist_list` is computed by calling the `whichdist_func_list()` function with `ac_coords`, `resource1`, and `forage` as arguments.
#' The `match_list` data frame is created by iterating over the length of `dist_list` and combining the relevant values from `resource1`, `ac_mat`, and `dist_list`.
#' Column names are assigned to `match_list` to indicate the corresponding variables.
#' The `ac_x` and `ac_y` columns are populated with the x and y coordinates from `ac_coords` based on the `ac_ind` values in `match_list`.
#' The `dist` column is computed by calling the `pairdist()` function on selected columns of `match_list`.
#' The unique values of `index` in `match_list` are stored in `ind_loop`.
#'
#' The `area_val`, `x`, and `y` vectors are initialized as empty.
#' A loop is executed over the indices in `ind_loop`, and for each index, a subset of `match_list` is created based on the index value and non-zero distances.
#' The `area_val` is computed using the `area_func()` function on the `ac` values in the subset, and the first values of `x` and `y` are stored.
#' The resulting `x`, `y`, and `area_val` vectors are combined into the `df_agg` data frame, with appropriate column names.
#'
#' Finally, the function returns the `ac` values from `df_agg`.
#'
#' @param resource1 A data frame or matrix representing the resource.
#' @param grid A grid object representing the spatial grid.
#' @param ac A raster object representing the variable 'ac'.
#'
#' @return A numeric vector representing the resource values based on area.
#'
#' @examples
#' res1 <- data.frame(x = c(1, 2, 3), y = c(4, 5, 6))
#' grid <- GridTopology(cellcentre.offset = c(0, 0), cellsize = c(1, 1), cells.dim = c(3, 3))
#' ac <- raster(matrix(1:9, nrow = 3))
#' res_value_function_area(res1, grid, ac)
#'
#' @export

res_value_function_area <- function(resource1, grid, ac){
  
  ac_mat = rasterToPoints(ac)
  colnames(ac_mat) = c('x','y','val')
  ac_coords = as.matrix(ac_mat[,1:2])
  dist_list = whichdist_func_list(ac_coords, resource1, forage)
  match_list = data.frame(do.call(rbind,lapply(1:length(dist_list), function(xi){
    cbind(resource1[xi,1],resource1[xi,2],xi,ac_mat[dist_list[[xi]],3])
  }
  )))
  
  match_list$ac_ind = unlist(dist_list)
  colnames(match_list) = c('x','y','index','ac','ac_ind')
  match_list$ac_x = ac_coords[match_list$ac_ind,1]
  match_list$ac_y = ac_coords[match_list$ac_ind,2]
  match_list$dist = pairdist(as.matrix(match_list[,1:2]), as.matrix(match_list[,6:7]))
  ind_loop = unique(match_list$index)
  area_val = NULL
  x = NULL
  y = NULL
  
  for (i in 1:length(ind_loop)){
    sub = match_list[match_list$index == ind_loop[i] & match_list$dist > 0,]
    area_val[i] = ifelse(is.null(dim(sub)),0,area_func(sub$ac))
    x[i] = sub[1,1]
    y[i] = sub[1,2]
  }
  
  df_agg = data.frame(x,y,area_val)
  colnames(df_agg) =  c('x','y','ac')
  return(df_agg$ac)
  
}


#' Calculate landscape value
#'
#' This function takes three inputs: `resource1`, `resource2`, and `ac_val`. It calculates the landscape value based on the given inputs.
#'
#' It starts by creating the `value_list` data frame, which contains the pairwise distances between `resource1` and `resource2`.
#' The `index` column is populated with values representing the indices of `resource1`.
#' Column names are assigned to `value_list` to indicate the corresponding variables.
#' The `inside` column is computed based on whether the distance values are less than or equal to `forage`.
#' The `value` column is computed by calling the `find_func` function on the `value_df$value_seq` values in `value_list$distance`.
#' The unique values of `index` in `value_list` are stored in `ind_loop`.
#' The `idw_val` vector is initialized as empty.
#'
#' A loop is executed over the indices in `ind_loop`, and for each index, a subset of `value_list` is created based on the index value.
#' The `idw_val` is computed using the `idw_func` function on the relevant subset of `value_list` where `inside` is TRUE.
#' The resulting `x`, `y`, and `idw_val` vectors are combined into the `agg_value` data frame, with appropriate column names.
#' The `landscape_total_value` is computed as the sum of the element-wise multiplication of `agg_value$value` and `ac_val`.
#'
#' Finally, the function returns the `landscape_total_value`.
#'
#' @param resource1 A data frame or matrix representing the first resource.
#' @param resource2 A data frame or matrix representing the second resource.
#' @param ac_val A numeric vector representing the variable 'ac'.
#'
#' @return A numeric value representing the landscape value.
#'
#' @examples
#' res1 <- data.frame(x = c(1, 2, 3), y = c(4, 5, 6))
#' res2 <- data.frame(x = c(2, 3, 4), y = c(5, 6, 7))
#' ac_val <- c(0.1, 0.2, 0.3)
#' landscape_value(res1, res2, ac_val)
#'
#' @export

landscape_value <- function(resource1, resource2, ac_val){
  
  value_list = data.frame(matrix(crossdist(resource1, resource2), ncol = 1))
  value_list$index = rep(seq(1, nrow(resource1), 1), times = nrow(resource2))
  colnames(value_list) = c('distance', 'index')
  value_list$inside = value_list$distance <= forage
  value_list$value = value_df$value[sapply(value_list$distance, find_func, value_df$value_seq)]
  ind_loop = unique(value_list$index)
  idw_val = NULL
  
  for (i in 1:length(ind_loop)){
    sub = value_list[value_list$index == ind_loop[i], ]
    idw_val[i] = ifelse(all(sub$inside == FALSE), 0, idw_func(sub[sub$inside == TRUE, 4], sub[sub$inside == TRUE, 1]))
  }
  
  agg_value = data.frame(resource1[, 1], resource1[, 2], idw_val)
  colnames(agg_value) = c('x', 'y', 'value')
  landscape_total_value = sum(agg_value$value * ac_val)
  return(landscape_total_value)
  
}


#' Calculate the sum of landscape values
#'
#' This function takes three inputs: `resource1`, `resource2`, and `ac_val`. It calculates the sum of landscape values based on the given inputs.
#'
#' It starts by creating the `value_list` data frame, which contains the pairwise distances between `resource1` and `resource2`.
#' The `index` column is populated with values representing the indices of `resource1`.
#' Column names are assigned to `value_list` to indicate the corresponding variables.
#' The `inside` column is computed based on whether the distance values are less than or equal to `forage`.
#' The `value` column is computed by calling the `find_func` function on the `value_df$value_seq` values in `value_list$distance`.
#' The unique values of `index` in `value_list` are stored in `ind_loop`.
#' The `sum_val` vector is initialized as empty.
#'
#' A loop is executed over the indices in `ind_loop`, and for each index, a subset of `value_list` is created based on the index value.
#' The `sum_val` is computed using the `d_sum_func_2` function on the relevant subset of `value_list` where `inside` is TRUE.
#' The resulting `x`, `y`, and `sum_val` vectors are combined into the `agg_value` data frame, with appropriate column names.
#' The `landscape_total_value` is computed as the sum of the element-wise multiplication of `agg_value$value` and `ac_val`.
#'
#' Finally, the function returns the `landscape_total_value`.
#'
#' @param resource1 A data frame or matrix representing the first resource.
#' @param resource2 A data frame or matrix representing the second resource.
#' @param ac_val A numeric vector representing the variable 'ac'.
#'
#' @return A numeric value representing the sum of landscape values.
#'
#' @examples
#' res1 <- data.frame(x = c(1, 2, 3), y = c(4, 5, 6))
#' res2 <- data.frame(x = c(2, 3, 4), y = c(5, 6, 7))
#' ac_val <- c(0.1, 0.2, 0.3)
#' landscape_value_sum(res1, res2, ac_val)
#'
#' @export

landscape_value_sum <- function(resource1, resource2, ac_val){
  
  value_list = data.frame(matrix(crossdist(resource1, resource2), ncol = 1))
  value_list$index = rep(seq(1, nrow(resource1), 1), times = nrow(resource2))
  colnames(value_list) = c('distance', 'index')
  value_list$inside = value_list$distance <= forage
  value_list$value = value_df$value[sapply(value_list$distance, find_func, value_df$value_seq)]
  ind_loop = unique(value_list$index)
  sum_val = NULL
  
  for (i in 1:length(ind_loop)){
    sub = value_list[value_list$index == ind_loop[i], ]
    sum_val[i] = ifelse(all(sub$inside == FALSE), 0, d_sum_func_2(sub[sub$inside == TRUE, 4], sub[sub$inside == TRUE, 1]))
  }
  
  agg_value = data.frame(resource1[, 1], resource1[, 2], sum_val)
  colnames(agg_value) = c('x', 'y', 'value')
  landscape_total_value = sum(agg_value$value * ac_val)
  return(landscape_total_value)
  
}



#' Calculate the cell values
#'
#' This function takes three inputs: `resource1`, `resource2`, and `ac_val`. It calculates the cell values based on the given inputs.
#'
#' It starts by computing the `maxvalue` using the `landscape_value_sum` function with `resource1`, `resource2`, and `ac_val` as arguments.
#' The `cell_value`, `maxval_new`, `xxi`, and `yyi` vectors are initialized as empty.
#'
#' A loop is executed over the rows of `resource1`, and for each row, a subset of `resource1` and `ac_val` is created by removing the current row.
#' The `maxvalue_new` is computed using the `landscape_value_sum` function on the updated `res1_new` and `resource2` with corresponding `val_ind`.
#' The difference between `maxvalue` and `maxvalue_new` is stored in `cell_value` for the current row, and the corresponding x and y coordinates are stored in `xxi` and `yyi`.
#'
#' The resulting `xxi`, `yyi`, and `cell_value` vectors are combined into the `res1_value` data frame, with appropriate column names.
#' The `cell_value`, `xxi`, and `yyi` vectors are reset as empty.
#'
#' Another loop is executed over the rows of `resource2`, and for each row, a subset of `resource2` is created by removing the current row.
#' The `maxvalue_new` is computed using the `landscape_value_sum` function on `resource1` and the updated `res2_new` with `ac_val`.
#' The difference between `maxvalue` and `maxvalue_new` is stored in `cell_value` for the current row, and the corresponding x and y coordinates are stored in `xxi` and `yyi`.
#'
#' The resulting `xxi`, `yyi`, and `cell_value` vectors are combined into the `res2_value` data frame, with appropriate column names.
#' The `res1_value` and `res2_value` data frames are combined using `rbind` to create the `study_value` data frame, which contains the cell values for both resources.
#'
#' Finally, the function returns the `study_value`.
#'
#' @param resource1 A data frame or matrix representing the first resource.
#' @param resource2 A data frame or matrix representing the second resource.
#' @param ac_val A numeric vector representing the variable 'ac'.
#'
#' @return A data frame with three columns: 'x', 'y', and 'val', representing the x and y coordinates and the corresponding cell values.
#'
#' @examples
#' res1 <- data.frame(x = c(1, 2, 3), y = c(4, 5, 6))
#' res2 <- data.frame(x = c(2, 3, 4), y = c(5, 6, 7))
#' ac_val <- c(0.1, 0.2, 0.3)
#' cellvalue(res1, res2, ac_val)
#'
#' @export

cellvalue <- function(resource1, resource2, ac_val){
  maxvalue = landscape_value_sum(resource1, resource2, ac_val)
  cell_value = maxval_new = xxi = yyi = NULL
  for (n in 1:nrow(resource1)) {
    res1_new = resource1[-n, ]
    val_ind = ac_val[-n]
    maxvalue_new = landscape_value_sum(res1_new, resource2, val_ind)
    cell_value[n] = maxvalue - maxvalue_new
    xxi[n] = resource1[n, 1]
    yyi[n] = resource1[n, 2]
  }
  res1_value = data.frame(xxi, yyi, cell_value)
  colnames(res1_value) = c('x', 'y', 'val')
  cell_value = xxi = yyi = NULL
  for (m in 1:nrow(resource2)){
    res2_new = resource2[-m, ]
    maxvalue_new = landscape_value_sum(resource1, res2_new, ac_val)
    cell_value[m] = ifelse((maxvalue - maxvalue_new) > 1, 1, maxvalue - maxvalue_new)
    xxi[m] = resource2[m, 1]
    yyi[m] = resource2[m, 2]
  }
  res2_value = data.frame(xxi, yyi, cell_value)
  colnames(res2_value) = c('x', 'y', 'val')
  study_value = rbind(res1_value, res2_value)
  return(study_value)
}



#' Calculate raster values
#'
#' This function takes two inputs: `cellvalue` and `grid`. It calculates raster values based on the given inputs.
#'
#' It starts by creating the `buff_seq` sequence from 1 to `buffer` with a step of 1 and calculates the corresponding `y_range` using a power transformation.
#' The `y_scale` is computed using the `scale_func` function on `y_range`.
#' The `buffer_df` data frame is created with `buff_seq` and `y_scale` as columns.
#' The `cellvalue_mat` is converted to a matrix, containing the x and y coordinates of the cell values.
#'
#' The `distbuffer` data frame is initialized with the pairwise distances between `grid` and `cellvalue_mat` using the `crossdist` function.
#' It includes columns for the index of the corresponding cell value, the value itself, the index of the grid cell, and the x and y coordinates of the grid cell.
#' Column names are assigned accordingly.
#' The `inside` column is created to determine if the distance is within the buffer range.
#' The `val_in` column is set to 0 for distances outside the buffer and to the corresponding cell value for distances inside the buffer.
#' The `distbuffer` is subsetted to include only rows where `inside` is TRUE.
#' The `weight` column is computed by finding the corresponding weights from `buffer_df$y_scale` based on the distances using the `find_func_all` function.
#' The `w_val` column is the element-wise product of `val_in` and `weight`.
#' The `agg_buff` data frame is created by aggregating the weighted values by the x and y coordinates of the grid cells using the `sum_func`.
#'
#' The `ras` raster is created using the `rasterFromXYZ` function with the aggregated `agg_buff` data frame.
#' Any NA values in `ras` are replaced with 0.
#' The resulting raster is returned.
#'
#' @param cellvalue A data frame with three columns: 'x', 'y', and 'val', representing the x and y coordinates and the corresponding cell values.
#' @param grid A data frame or matrix representing the grid cells.
#'
#' @return A raster object representing the calculated raster values.
#'
#' @examples
#' cellval <- data.frame(x = c(1, 2, 3), y = c(4, 5, 6), val = c(0.1, 0.2, 0.3))
#' grid <- expand.grid(x = 1:5, y = 1:5)
#' raster_value_func(cellval, grid)
#'
#' @export

raster_value_func <- function(cellvalue, grid){
  buff_seq = seq(1, buffer, 1)
  y_range = -1 * (buff_seq^0.2) + buffer
  y_scale = scale_func(y_range)
  buffer_df = data.frame(buff_seq, y_scale)
  cellvalue_mat = as.matrix(cellvalue[, 1:2])

  distbuffer = data.frame(matrix(crossdist(grid, cellvalue_mat), ncol = 1))
  distbuffer$index = rep(1:nrow(cellvalue_mat), each = nrow(grid))
  distbuffer$value = cellvalue[distbuffer$index, 3]
  distbuffer$grid_index = rep(1:nrow(grid), times = nrow(cellvalue))
  distbuffer$x = grid[distbuffer$grid_index, 1]
  distbuffer$y = grid[distbuffer$grid_index, 2]
  colnames(distbuffer) = c('dist', 'index', 'value', 'grid_index', 'x', 'y')
  distbuffer$inside = distbuffer$dist <= buffer
  distbuffer$val_in = ifelse(distbuffer$inside == FALSE, 0, distbuffer$value)
  distbuffer = distbuffer[distbuffer$inside == TRUE, ]
  distbuffer$weight = buffer_df$y_scale[sapply(distbuffer$dist, find_func_all, buffer_df$buff_seq)]
  distbuffer$w_val = distbuffer$val_in * distbuffer$weight
  agg_buff = aggregate(distbuffer$w_val, list(distbuffer$x, distbuffer$y), sum_func)

  ras = rasterFromXYZ(agg_buff)
  ras[is.na(ras)] = 0
  return(ras)
}



#' Calculate passage values
#'
#' This function takes five inputs: `res1`, `res2`, `value`, `rasval`, and `perm`. It calculates passage values based on the given inputs.
#'
#' It starts by splitting the `value` data frame into two parts: `p_1` and `p_2` based on the number of rows in `res1`.
#' Random sampling is performed on `p_1` to obtain a subset `p_1_samp` of size 20, with sampling probabilities determined by the third column of `p_1`.
#' `p_2_samp` is set to `p_2` for now.
#' The transition matrix `tran` is created using the `transition` function, where the transition probabilities are set to the inverse of the mean value.
#' The `tran_g` is the geometric correction of `tran`.
#' An empty vector `lraster` is created to store the passage rasters.
#' A loop is initiated for each row in `p_1`, where an origin point `opoint` is randomly sampled from `p_1_samp` and a goal point `gpoint` is randomly sampled from `p_2_samp`.
#' The `passage` function is then used to calculate the passage raster based on `tran_g`, `opoint`, `gpoint`, and a threshold of 0.5.
#' The resulting passage raster is appended to `lraster`.
#' The loop progress is printed for each iteration.
#'
#' The sum of all passage rasters in `lraster` is computed using the `Reduce` function with addition as the binary operator.
#' The `p_v` is set to `value`.
#' A `SpatialPointsDataFrame` `p` is created using the coordinates from `p_v` and the corresponding data.
#' The cells in `sum_rast` corresponding to the points in `p` are set to 0.
#' The scaling function `scale_func_ras` is applied to `sum_rast` to scale the values.
#' The `pass_ras` is created as a stack of `pass_ras` and `rasval` rasters.
#' The sum of all layers in `ras_list` is calculated using the `calc` function with the `sum` aggregation method.
#' Any values greater than 1 in `ras_sum` are set to 1.
#' The cells in `ras_sum` corresponding to the points in `p` are set to the corresponding values in `rasval`.
#' The `fix_points` matrix is created by combining the x and y coordinates from `value` with the values from `ras_sum` at the corresponding cell locations.
#' Column names are assigned accordingly.
#' The `ras_points` matrix is created by converting `ras_sum` to points using the `rasterToPoints` function.
#' The `fix_list` data frame is created with the pairwise distances between `ras_points` and `fix_points`, including columns for the index, value, and inside (indicating distances within 2 units).
#' The subset `fix_list_sub` includes rows where `inside` is TRUE and distance is greater than 0.
#' The `fix_agg` data frame is created by aggregating the values from `fix_list_sub` by the index, taking the mean value.
#' The cells in `ras_sum` corresponding to the points in `value` are updated with the values from `fix_agg`.
#' The resulting `ras_sum` raster is returned.
#'
#' @param res1 A data frame representing resource 1.
#' @param res2 A data frame representing resource 2.
#' @param value A data frame containing values associated with the resources.
#' @param rasval A raster representing the initial values of the landscape.
#' @param perm A permutation representing the direction of movement.
#' @return A raster representing the passage values.
#'
#' @examples
#' # Example usage of passage_func
#' passage_func(res1, res2, value, rasval, perm)
#'
#' @export
passage_func = function(res1, res2, value, rasval, perm) {
  p_1 = value[1:nrow(res1), ]
  p_2 = value[nrow(res1):nrow(value), ]
  
  p_1_samp = p_1[sample(1:nrow(p_1), size = 20, prob = p_1[, 3]), ]
  p_2_samp = p_2
  tran = transition(perm, function(x) 1 / mean(x), direction = 16, symm = FALSE)
  tran_g = geoCorrection(tran)
  
  lraster = c()
  
  for (i in 1:nrow(p_1)) {
    opoint = as.matrix(p_1[sample(1:nrow(p_1), size = 1, prob = p_1[, 3]), ][, 1:2], ncol = 2)
    gpoint = as.matrix(p_2[sample(1:nrow(p_2), size = 1, prob = p_2[, 3]), ][, 1:2], ncol = 2)
    passraster = passage(tran_g, opoint, gpoint, 0.5)
    lraster = append(lraster, passraster)
    print(i)
  }
  
  sum_rast = Reduce("+", lraster)
  p_v = value
  sum_rast[cellFromXY(sum_rast, p_v@coords)] = 0
  
  pass_ras = scale_func_ras(sum_rast)
  
  ras_list = stack(pass_ras, rasval)
  ras_sum = calc(ras_list, sum)
  ras_sum[ras_sum > 1] = 1
  
  fix_points = as.matrix(cbind(value[, 1:2], ras_sum[cellFromXY(ras_sum, value[, 1:2])]))
  colnames(fix_points) = c('x', 'y', 'value')
  ras_points = as.matrix(rasterToPoints(ras_sum))
  
  fix_list = data.frame(matrix(crossdist(ras_points[, 1:2], fix_points[, 1:2]), ncol = 1))
  fix_list$index = rep(seq(1, nrow(fix_points), 1), each = nrow(ras_points))
  fix_list$val = rep(ras_points[, 3], times = nrow(fix_points))
  colnames(fix_list) = c('distance', 'index', 'val')
  fix_list$inside = fix_list$distance <= 2
  fix_list_sub = fix_list[fix_list$inside == TRUE & fix_list$distance > 0, ]
  fix_agg = aggregate(fix_list_sub$val, list(fix_list_sub$index), mean)
  ras_sum[cellFromXY(ras_sum, value[, 1:2])] = fix_agg$x
  
  return(ras_sum)
}


#' Calculate Pareto Efficiency Frontier
#'
#' This function takes two inputs: `ras` and `conf`. It calculates the Pareto Efficiency Frontier based on the given inputs.
#'
#' The `ras` raster and `conf` raster are updated by replacing NA values with 0.
#' The `conf` and `ras` rasters are stacked together.
#' The `conf_grid` is created by converting the `conf` raster to points using the `rasterToPoints` function.
#' The `val_grid` is created by converting the `ras` raster to points using the `rasterToPoints` function.
#' Column names of `val_grid` are updated to 'x', 'y', and 'pass'.
#' The `merge_land` is created by merging `val_grid` and `conf_grid` based on the 'x' and 'y' columns, with all values from `val_grid` and matching values from `conf_grid`.
#' The `seq_loop` is created as a sequence from the minimum value of 'pass' in `merge_land` to the maximum value of 'pass' in `merge_land`, with a step of 0.01.
#' 
#' Various variables are initialized for further calculations: `noconf_baseline`, `conf_baseline`, `conf_extent`, `conf_amount`, `value_conf`, `landstate_df`, `subset_conf_0`, `conf_0_sum`, `value_calc`, `range_new`, `value_range`, `range_count`, `counter_1`, `counter_2`, `sub_choice`, and `conf_amount_ls`.
#' 
#' A loop is initiated for each value in `seq_loop`. Within the loop, the following steps are performed:
#'   - If the current value in `range_count` is 0, a subset of `merge_land` is created where 'pass' is 0 and 'conflict' is 1.
#'   - If the current value in `range_count` is not 0, a subset of `merge_land` is created where 'pass' is less than or equal to the current `range_count`, greater than `range_new`, and 'conflict' is 1.
#'   - The subset of `merge_land` is appended to `landstate_df`.
#'   - The sum of 'conflict' in the subset is calculated and stored in `conf_sum`.
#'   - The sum of 'pass' in the subset is calculated and stored in `value_sum`.
#'   - The current `range_count` becomes the new `range_new`.
#'   - The current `value_sum` becomes the new `value_range`.
#'   - The new `value_calc` is calculated by subtracting `value_range` from the previous `value_calc`.
#'   - The new `value_calc` is stored in `value_conf` at the corresponding index `counter_1`.
#'   - The `conf_extent` is stored in `conf_amount` at the corresponding index `counter_2`.
#'   - The `conf_sum` and `conf_extent` are added together and stored in `conf_amount` at the corresponding index `counter_1`.
#'   - The `conf_amount` is stored in `frontier_value_ind`.
#'   - The `landstate_df` and `conf_amount` are stored in `sub_choice` and `conf_amount_ls`, respectively.
#'   - The `counter_1` and `counter_2` are incremented by 1.
#'
#' The `df_frontier` data frame is created with columns `value_conf`, `conf_amount`, `percentage`, and `area`. The `percentage` column represents the percentage of `value_conf` relative to `noconf_baseline`, and the `area` column represents the percentage of `conf_amount` relative to `conf_baseline`.
#'
#' The Pareto Efficiency Frontier plot is generated using the `plot` function, displaying the `area` values on the x-axis and `percentage` values on the y-axis. The plot includes a red line connecting the points, labeled as "Pareto Frontier".
#'
#' The user is prompted to click on the plot to choose a land use allocation, which is stored in the `choice` variable.
#'
#' The `choice_ls` is extracted from `sub_choice` based on the chosen `choice` value, and the `choice_ras` raster is created using the `rasterFromXYZ` function.
#'
#' The `choice_ras` raster is returned as the result of the function.
#'
#' @param ras A raster representing the landscape.
#' @param conf A raster representing the conflict values.
#' @return A raster representing the Pareto Efficiency Frontier.
#'
#' @examples
#' # Example usage of pareto_func
#' pareto_func(ras, conf)
#'
#' @export
pareto_func = function(ras,conf){
    ras = ras
    conf[is.na(conf)] = 0
    ras[is.na(ras)] = 0

    stack(conf,ras)

    conf_grid = rasterToPoints(conf)
    val_grid = rasterToPoints(ras)
    colnames(val_grid) = c('x','y','pass')

    merge_land = merge(val_grid,conf_grid, by = c('x','y'), all.x = TRUE)
    seq_loop = seq(min(merge_land$pass),max(merge_land$pass), 0.01)


    noconf_baseline = sum(merge_land$pass)
    conf_baseline = sum(merge_land$conflict)
    conf_extent = conf_amount = value_conf = NULL
    landstate_df = setNames(data.frame(matrix(ncol = 4, nrow = 0)), names(merge_land))
    subset_conf_0 = subset(merge_land, merge_land$conflict == 0)
    conf_0_sum = sum(subset_conf_0$conflict)
    value_conf[1] = conf_extent[1] = value_calc =  noconf_baseline
    conf_amount[1] = conf_0_sum
    range_new = value_range = range_count = 0  
    counter_1 = 2
    counter_2 = 1

    sub_choice = conf_amount_ls = list()

    for (range_count in seq_loop){
      value_calc_new = value_calc

      if (range_count == 0){
        subset_landscape = subset(merge_land, merge_land$pass == 0 
                                  & merge_land$conflict == 1 )
      } else {
        subset_landscape = subset(merge_land, merge_land$pass <= range_count &
                                    merge_land$pass > range_new &
                                    merge_land$conflict == 1 )
      }

      landstate_df = rbind(landstate_df,subset_landscape)
      conf_sum = sum(subset_landscape$conflict)
      value_sum = sum(subset_landscape$pass)
      range_new = range_count
      value_range = value_sum
      value_calc = value_calc_new - value_range
      value_conf[counter_1] = value_calc
      conf_extent = conf_amount[counter_2]
      conf_amount[counter_1] = conf_sum+conf_extent
      frontier_value_ind = conf_amount
      sub_choice[[counter_2]] = landstate_df
      conf_amount_ls[[counter_2]] = conf_amount
      counter_1= counter_1+1
      counter_2= counter_2+1

    }


    df_frontier = data.frame(value_conf,conf_amount)
    df_frontier$percentage = df_frontier$value_conf/noconf_baseline*100
    df_frontier$area = df_frontier$conf_amount/conf_baseline*100


    options(scipen = 999)
    # par(mfrow = c(4,2), oma = c(1, 1, 1, 1), mar = c(4, 4, 1, 3) + 0.1)
    plot(df_frontier$area,df_frontier$percentage, ylab ='Total Habitat Value (%)',
         type = 'o',xlab = "", col = 'red', lwd = 4, cex = 0)
    title(xlab = "Area Cleared (%)")
    title('Pareto Efficiency Frontier', line = 0.2)

    #######IMPORTANT#######################################################
    #######click on the fontier plot to choose a land use allocation 'finish' in top right of plot
    choice = identify(df_frontier$area,
                      df_frontier$percentage, pos = FALSE, tolerance = 5,
                      atpen = FALSE, labels = '')


    ##########################################################################
    v = df_frontier[choice,4]
    abline(v =  v, lty = 5, lwd = 2)
    legend("topright", legend=c("Pareto Frontier", "Decision"),
           col=c("red", "black"), lty=c(1,5),lwd = 2, cex=1)

    choice_ls = sub_choice[[choice-1]][,c(1,2,4)]
    choice_ras = rasterFromXYZ(choice_ls)


    # up_ras = ras
    # # plot(up_ras)
    # 
    # up_ras[cellFromXY(up_ras, choice_ls[,1:2])] = 0
    # # plot(up_ras)
    # 
    # plot(up_ras, col = cols_out, legend = TRUE,
    #      asp =0, axes = T,xlab = 'x' , interpolate = TRUE,
    #      legend.width=2, legend.shrink= 1, legend.args = list(text = '', side = 4, 
    #                                                           font = 2, line = 2.5, cex = 1))
    # #Habitat Value
    # # title('Optimised Conservation Outcome', line = 0.2)
    # title(ylab = 'y')
    # 
    # plot(choice_ras, legend = FALSE, col= 'red', add = TRUE)
    # legend(1,300, legend = 'Clearing Extent', col = 'red', pch = 15, cex = 1.3)

    return(choice_ras)
}

