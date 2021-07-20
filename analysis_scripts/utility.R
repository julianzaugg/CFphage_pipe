
# Function to convert a dataframe to a matrix. First column becomes row names
df2m <- function(mydataframe){
  mymatrix <- mydataframe
  rownames(mymatrix) <- mydataframe[,1]
  mymatrix[,1] <- NULL
  mymatrix <- as.matrix(mymatrix)
  mymatrix
}

# Function to convert a matrix to a dataframe. Rows become first column.
# Specify 'column_name' to set the name of the new column (defaults to "Row_variable")
m2df <- function(mymatrix, column_name = "Row_variable"){
  mydf <- as.data.frame(mymatrix)
  cur_names <- names(mydf)
  mydf[, column_name] <- rownames(mydf)
  rownames(mydf) <- NULL
  mydf <- mydf[,c(column_name,cur_names)]
  return(mydf)
}

# Make row names of a dataframe the value of a defined column (default 1st) and then remove the column
df_column_to_rownames <- function(mydf, rowname_col = 1){
  my_clean.df <- mydf
  rownames(my_clean.df) <- my_clean.df[,rowname_col]
  my_clean.df[,1] <- NULL
  return(my_clean.df)
}

# Assign colours to a dataframe for specified discrete columns.
assign_colours_to_df <- function(dataframe.df,
                                 columns,
                                 my_palette = NULL,
                                 auto_assign =T,
                                 my_default_palette=NULL){
  # TODO - make option for default continuous palette
  # Palettes
  colour_palette_soft_8 <- c("#8b90bc","#76cc5d","#9558b7","#d2c351","#cd5f88","#89cab7","#d06842","#858658")
  colour_palette_10 <- c("#8eec45","#0265e8","#f6a800","#bf6549","#486900","#c655a0","#00d1b6","#ff4431","#aeb85c","#7e7fc8")
  colour_palette_15 <- c("#77b642","#7166d9","#cfa240","#b351bb","#4fac7f","#d44891","#79843a","#c68ad4","#d15a2c","#5ba7d9","#ce4355","#6570ba",
                         "#b67249","#9b4a6f","#df8398")
  colour_palette_20 <- c("#0057b4","#7fff56","#d600bc","#d8d500","#e76eff","#019932","#9f8fff","#ffc730","#007fac","#a20019","#06fefd","#ff6782",
                         "#00774c","#e0c8ff","#717a00","#4b2952","#e2ed7d","#46321e","#ffbd76","#ffb4c6")
  colour_palette_30 <- c("#009348","#f579fe","#4fe16e","#b40085","#4d7e00","#4742b4","#f0c031","#016dd9","#d45200","#7499ff","#ef4d2d","#01c9c8",
                         "#f8394b","#88d7a6","#d20063","#c8cc5d","#882986","#fdb95d","#404f8f","#917300","#f3aefc","#5c5800","#ff75c3","#00674a",
                         "#ba001c","#979760","#8b354c","#ff875f","#943105","#cf9478")
  colour_palette_206 <- c("#cfefb4","#7d8b00","#a70079","#552155","#632900","#ffb173","#fbdcf2","#015a6a","#43fdf7","#ff443a","#008186","#3b8aff",
                          "#8b5fff","#ff9777","#4200a9","#85f6fd","#c96000","#36218a","#d28900","#0137d7","#30325b","#ff836b","#008b4f","#21ff9d",
                          "#00794d","#870052","#e9ec4b","#ce006b","#6e0044","#8a6500","#006971","#432e4b","#ca8dff","#f20059","#44ffe2","#00be5c",
                          "#a0d2ff","#1914ab","#4d284e","#59d7ff","#ab9aff","#0151d9","#1de740","#e24500","#9fc400","#610769","#0a4600","#1e365b",
                          "#018f3f","#b15fff","#009c5e","#005290","#506100","#f49aff","#0187c1","#ffb5f4","#daf100","#70081d","#ff9890","#c1baff",
                          "#ffbe5a","#1b3466","#ff2a7f","#ff5d3c","#e47800","#ac6bff","#1f6000","#006627","#4f4000","#dcd6ff","#ffd7c1","#ed2de4",
                          "#a50038","#a5a8ff","#0f2f7f","#b11700","#00e06b","#ffabb8","#015780","#82eaff","#1b2a88","#6f1600","#d3ef9c","#746e00",
                          "#01d851","#625300","#01d799","#96fd6c","#ff5ca1","#7b0017","#004c2b","#baf678","#f8aaff","#007c1b","#01a88a","#a71ed8",
                          "#fb8cff","#840079","#276d00","#556655","#02b0de","#c0efd7","#63193e","#8e9984","#017ac9","#ff925f","#ff63d7","#294100",
                          "#28baff","#5b2523","#35ab00","#69132e","#8a3b00","#a67700","#7fff6a","#002f96","#681a0b","#4d3003","#ff7de6","#0190d8",
                          "#a69700","#ff6282","#d3f266","#ffc4cf","#ffac3c","#d064ff","#d07aff","#c3005d","#9d0067","#0167c1","#8cfe82","#ffd68f",
                          "#8cfcaf","#f50096","#00c2a2","#aa5e00","#02c16d","#4e4bf6","#ffd962","#004793","#93d800","#462a58","#323a03","#4f9eff",
                          "#2b3a25","#2defff","#02edd6","#864e00","#ffc59f","#e7e9ab","#014cc4","#437bff","#00afba","#ff7d82","#8a1ed4","#ff48b3",
                          "#acf7ab","#005550","#7600a6","#bc0028","#00adab","#02dfbf","#ba004c","#004760","#ebc5ff","#0162d7","#9b3900","#5869ff",
                          "#ff6160","#87b6ff","#ff6796","#ff8422","#ff8440","#b500a8","#937fff","#0132bd","#f48e00","#1e8800","#462370","#3e3614",
                          "#9ca800","#efe5bf","#aeb6a0","#d9aaff","#d8ef89","#cec800","#ffb8b3","#4a2c42","#01715b","#b8ebff","#ff9ec0","#ff93ec",
                          "#ffe0aa","#65b300","#6a8b00","#f6e77c","#ff85c0","#5de522","#a5f6ca","#c70077","#5a4149","#a3b700","#ff63c4","#63fecd",
                          "#93f6e7","#01b4a4")
  # Create palette list for each column
  palette_list <- list()
  for (col in columns){
    column_values <- as.character(unique(dataframe.df[,col]))
    column_values <- column_values[!is.na(column_values)]
    n_values <- length(column_values)
    if (is.null(my_default_palette)){
      if (n_values <= 10){
        # default_colour_palette <- setNames(colorRampPalette(colour_palette_10)(n_values), column_values)
        default_colour_palette <- colour_palette_10
      } else if (n_values > 10 & n_values <= 15){
        # default_colour_palette <- setNames(colorRampPalette(colour_palette_15)(n_values), column_values)
        default_colour_palette <- colour_palette_15
      } else if (n_values > 15 & n_values <= 20){
        # default_colour_palette <- setNames(colorRampPalette(colour_palette_20)(n_values), column_values)
        default_colour_palette <- colour_palette_20
      } else if (n_values > 20 & n_values <= 30){
        # default_colour_palette <- setNames(colorRampPalette(colour_palette_30)(n_values), column_values)
        default_colour_palette <- colour_palette_30
      } else{
        # default_colour_palette <- setNames(colorRampPalette(colour_palette_206)(n_values), column_values)
        default_colour_palette <- colour_palette_206
      }
    } else{
      default_colour_palette <- my_default_palette
    }

    if (typeof(my_palette) == "character"){
      if (length(my_palette) >= n_values){
        palette_list[[col]] <- my_palette
      } else{
        stop(paste0("Provided palette has ", length(my_palette), " colours, but column ", col, " has ", n_values, ". Provide a longer palette."))
      }
    }
    else if (is.null(my_palette) & auto_assign == T){
      palette_list[[col]] <- default_colour_palette
    }
    else if (typeof(my_palette) == "list"){
      if (! col %in% names(my_palette)){
        if (auto_assign == T){
          warning(paste0("Column ", col, " does not have a specified palette. Assigning a default palette"))
          palette_list[[col]] <- default_colour_palette
        } else{
          warning(paste0("Column ", col, " does not have a specified palette. Skipping."))
        }
      } else{
        palette_list[[col]] <- my_palette[[col]]
      }
    }
    else {
      stop("Provided palette should be a character vector of colours or a list of separate character vectors, or set auto_assign=T")
    }
    # print(col)
    # print(palette_list)
    # print(palette_list[[col]])
    # print(typeof(palette_list[[col]]))
    # print(palette_list[col])
    # print(typeof(palette_list[col]))
    # print("*******")
    if (!is.null(names(palette_list[[col]]))){
      # print(col)
      col_colours <- palette_list[[col]]
    } else{
      # print(col)
      col_colours <- setNames(palette_list[[col]][1:n_values], column_values)
    }

    dataframe.df[,paste0(col, "_colour")] <- as.character(lapply(as.character(dataframe.df[,col]), function(x) as.character(col_colours[x])))
  }
  dataframe.df
}

darken <- function(color, factor=1.4){
  col <- col2rgb(color)
  col <- col/factor
  col <- rgb(t(col), maxColorValue=255)
  col
}

lighten <- function(color, factor=1.4){
  col <- col2rgb(color)
  col <- col*factor
  col <- rgb(t(col), maxColorValue=255)
  col
}



first_resolved_taxonomy <- function(x) {
  # Get the lowest taxonomy level (or 'first') that has been resolved from a taxonomy string
  ranks.v <- unlist(strsplit(x, split = ';'))
  # print(x)
  # print(ranks.v)
  for (i in length(ranks.v):1) {
    split.v <- unlist(strsplit(ranks.v[i], split = '__'))
    if (!is.na(split.v[2]) & split.v[2] != "") {
      if (!grepl(split.v[2], pattern = "Unassigned|uncultured")) {
        return(ranks.v[i])
      }
    }
  }
  return(ranks.v[1])
}
