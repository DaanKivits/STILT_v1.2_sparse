create_sparse_matrix <- function(array){
## Purpose: Create a sparse matrix from a 2D array in R
## Input: 2D array 
## Output: Sparse matrix
## NOTE: This function is used to create sparse footprints for STILT specifically,
## therefore, index_x is the column index and index_y is the row index, since this is
## how STILT stores the footprints. 

indices <- which(array != 0, arr.ind = T)
values <- array[indices]
index_x = indices[, 'col']
index_y = indices[, 'row']    

return(list("values" = values, "index_x" = index_x, "index_y" = index_y))
}