#' @title Create Boolean Matrix
#' @export
# more options: https://stackoverflow.com/questions/17431524/create-a-binary-indicator-matrix-boolean-matrix-in-r
create_boolean_matrix <- function(kegg_merged) {
	kegg_merged = as.matrix(kegg_merged)
	rowVals <- sort(unique(unlist(kegg_merged[,1])))
	colVals <- sort(unique(unlist(kegg_merged[,2])))
	boolean_matrix <- matrix(0, nrow = length(rowVals), ncol = length(colVals), dimnames = list(rowVals, colVals))
	boolean_matrix[kegg_merged[,c(1:2)]] <- table(kegg_merged[,1], kegg_merged[,2])[kegg_merged[,c(1:2)]]

	boolean_matrix[boolean_matrix > 1] = 1
	return(boolean_matrix)
}

#' @title Create Square Matrix
#' @export
create_square_matrix <- function(data_table) {
	# https://stackoverflow.com/questions/46562173/creating-a-square-matrix-from-a-data-frame
	# get names for row and columns
	nameVals <- sort(unique(unlist(data_table[,c(1:2)])))
	# construct 0 matrix of correct dimensions with row and column names
	square_matrix <- matrix(0, length(nameVals), length(nameVals), dimnames = list(nameVals, nameVals))
	# fill in the matrix with matrix indexing on row and column names
	square_matrix[as.matrix(data_table[,c(1:2)])] <- data_table[[3]] #data_table[["combined_score"]]
	diag(square_matrix) = 1000
	return(square_matrix)
}

#' @title  Convert Matrix to Data.table
#' @export
matrix_to_collapsed_dt <- function(boolean_matrix) {
	collapsed_dt = data.table(matrix(data = rownames(boolean_matrix), nrow = dim(boolean_matrix)[1], ncol = 2))
	names(collapsed_dt) = c("Keys", "Values")
	setkey(collapsed_dt, Keys)

	for (row_index in 1:dim(boolean_matrix)[1]) {
		m = dimnames(boolean_matrix[row_index, boolean_matrix[row_index,] == 1, drop = F])
		#m_key = m[[1]]
		m_val = m[[2]] #as.list(m[[2]])

		collapsed_dt[row_index, 2] = paste(unlist(m_val), collapse = ", ")
	}
	return(collapsed_dt)
}
