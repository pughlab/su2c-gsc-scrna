library(Matrix)
library(Seurat)
library(harmony)
library(scran)
library (AUCell)
library(ggplot2)
library(SummarizedExperiment)
library(R.devices)
library(reticulate)
library(GSA)
library(reshape2)
library(scales)
library(cowplot)
library(mclust)

LM <- reticulate::import('sklearn.linear_model')
SKLEARN_NEIGHBORS <- reticulate::import('sklearn.neighbors')
MODEL_SELECTION <- reticulate::import('sklearn.model_selection')
METRICS <- reticulate::import('sklearn.metrics')

subset_seurat <- function(inp_seurat_obj, subset_size = 5000) {
  set.seed(12345)
  cell_names <- rownames(inp_seurat_obj@meta.data)
  cells_subs <- sample(cell_names, size = subset_size, replace = FALSE)
  inp_seurat_obj <- SubsetData(inp_seurat_obj, cells = cells_subs, do.clean = TRUE, subset.raw = TRUE)
}

AUCell_batch <- function(inp_data, genesets, num_batches = 100) {
  ## Scores a data matrix with AUCell in batches. Idea is to limit memory consumption when
  ## scoring with AUCell
  ## INPUTS:
  ##    inp_data = input data, either a dxn matrix of d features, n samples or a Seurat object
  ##                containing such a matrix
  ##    genesets = named list of character vectors, each consisting of a set of gene symbols
  ##    num_batches = number of batches to run AUCell for. More batches = fewer cells (observations)
  ##                  for each batch used for scoring
  ##    slot = slot to use if using a Seurat object
  ##    assay = assay to use if using a Seurat object
  ## RETURNS:
  ##  either an nxp matrix (samples x scores)
  if (is.matrix(inp_data) || is(inp_data, 'dgCMatrix')) {
    num_cells <- ncol(inp_data)
    batch_size <- ceiling(num_cells/num_batches)
    score_mat <- c()
    print('Running AUCell scoring')
    Sys.time()
    for (i in 1:num_batches) {
      print(paste('batch', i, Sys.time()))
      ind1 <- (i-1)*batch_size + 1
      ind2 <- i*batch_size
      if (ind2 > num_cells) {
        ind2 <- num_cells
      }
      gene_rankings <- AUCell::AUCell_buildRankings(inp_data[,ind1:ind2], plotStats = FALSE)
      score_mat_i <- AUCell::AUCell_calcAUC(geneSets = genesets, rankings = gene_rankings)
      score_mat_i <- t(assay(score_mat_i, 'AUC'))
      score_mat <- rbind(score_mat, score_mat_i)
      gc(full = TRUE, verbose = TRUE)
    }
    print('Finished Scoring')
    print(Sys.time())
    return(score_mat)
  } else if (class(inp_data) == 'seurat') {
    score_mat <- AUCell_batch(inp_data = GetAssayData(inp_data, slot = 'raw.data', assay = 'RNA'), 
                              genesets = genesets, 
                              num_batches = num_batches)
    colnames(score_mat) <- paste0(colnames(score_mat), '_AUC')
    inp_data <- AddMetaData(inp_data, as.data.frame(score_mat))
  }
}

calc_feat_diff <- function(seurat_obj, var1, var2, new_feat_name, scale = FALSE) {
  ## Calculate the difference between 2 numeric features
  ## INPUTS:
  ##  seurat_obj = seurat object
  ##  var1 = first feature, a character vector of length 1. can be in metadata
  ##          or in data slot for seurat object
  ##  var2 = second feature, a character vector of length 1
  ##  new_feat_name = character vector of length 1, denotes new name of difference
  ##                  vector between var1 and var2, to be added to metadata of seurat_obj
  ## RETURNS:
  ##  Seurat object with difference vector added to metadata
  
  scale_data <- function(x) {
    return((x - min(x))/(max(x) - min(x)))
  }
  
  fetched_data <- Seurat::FetchData(seurat_obj, vars.all = c(var1, var2))
  stopifnot(is.numeric(fetched_data[,var1]))
  stopifnot(is.numeric(fetched_data[,var2]))
  feat_vect1 <- fetched_data[,var1]
  feat_vect2 <- fetched_data[,var2]
  if (scale) {
    feat_vect1 <- scale_data(feat_vect1)
    feat_vect2 <- scale_data(feat_vect2)
  }
  diff_vect <-  feat_vect1 - feat_vect2
  names(diff_vect) <- rownames(fetched_data)
  seurat_obj <- Seurat::AddMetaData(seurat_obj, metadata = diff_vect, col.name = new_feat_name)
  return(seurat_obj)
}

regress_data <- function(seurat_obj, vars_regress, genes_use = NULL) {
  # routine for getting scaled data for Seurat, regressing out specified variables,
  # scaling specified features, or by default Variable Genes
  # get nUMI.
  nUMI <- Matrix::colSums(seurat_obj@raw.data)
  seurat_obj@meta.data$nUMI <- nUMI
  # calculate percent mito as sum of mitochondrial features / sum of all features
  mito_genes <- grep('^MT-', rownames(GetAssayData(seurat_obj, assay = 'RNA', slot = 'raw.data')))
  mito_counts <- Matrix::colSums(GetAssayData(seurat_obj, assay = 'RNA', slot = 'raw.data')[mito_genes,])
  seurat_obj@meta.data$percent.mito <- mito_counts/nUMI
  # calculate difference in expression of s-phase and G2M genes
  # read in a list of cell cycle markers, from Tirosh et al, 2015.
  # We regress out cell cycle difference as suggested in alternate workflow
  # outlined by the Satija Lab: https://satijalab.org/seurat/v2.4/cell_cycle_vignette.html
  # Change filepath for genelist if necessary
  cc_genes <- readLines(con = "~/projects/su2c_v2/data/raw/GeneSets/regev_lab_cell_cycle_tirosh_2015/regev_lab_cell_cycle_genes.txt")
  s_genes <- cc_genes[1:43]
  g2m_genes <- cc_genes[44:97]
  seurat_obj <- Seurat::CellCycleScoring(seurat_obj, 
                                         g2m.genes = g2m_genes, 
                                         s.genes = s_genes)
  seurat_obj <- calc_feat_diff(seurat_obj = seurat_obj, var1 = 'G2M.Score', var2 = 'S.Score', new_feat_name = 'CC.Difference')
  seurat_obj <- Seurat::ScaleData(seurat_obj, 
                                  vars.to.regress = vars_regress,
                                  do.scale = TRUE, 
                                  do.center = TRUE,
                                  genes.use = genes_use)
}

seurat_pipeline <- function(seurat_obj, vars_regress = c('percent.mito', 'nUMI', 'CC.Difference')) {
  ## INPUTS:
  ##  seurat_obj = seurat v2.3.4 object
  ##  vars_regress = variables whose effects are to be regressed from data via linear model
  ##                  note that cell cycle scores are recalculated as they are relative
  ##                  scores among the population of cells in input data
  ## RETURNS:
  ##  seurat object with or without batch correction performed
  # get log normalized data
  seurat_obj <- Seurat::NormalizeData(seurat_obj,
                                      assay.type = "RNA", 
                                      normalization.method = "LogNormalize",
                                      scale.factor = 10000,
                                      display.progress = TRUE)
  # find variable genes
  seurat_obj <- Seurat::FindVariableGenes(object = seurat_obj, 
                                          mean.function = ExpMean, 
                                          dispersion.function = LogVMR, 
                                          x.low.cutoff = 0.0125, 
                                          x.high.cutoff = 3, 
                                          y.cutoff = 0.5,
                                          do.plot = FALSE)
  # correct for covariates not of interest
  seurat_obj <- regress_data(seurat_obj = seurat_obj, vars_regress = vars_regress)
  # run pca without ribosomal genes
  ribo.genes <- c(rownames(seurat_obj@data)[grep("^RP[1:9]", rownames(seurat_obj@data))],
                  rownames(seurat_obj@data)[grep("^RP[L,S]", rownames(seurat_obj@data))]
  )
  seurat_obj <- Seurat::RunPCA(object = seurat_obj, 
                               pc.genes = rownames(seurat_obj@data)[!rownames(seurat_obj@data) %in% ribo.genes], #all genes except ribo
                               do.print = FALSE, 
                               pcs.compute = min(100, nrow(seurat_obj@meta.data) - 1)
                               )
  return(seurat_obj)
}

seurat_subroutine <- function(raw_data, meta_data, min_cells = 20, min_genes = 200) {
  print('Creating Seurat Object')
  seurat_obj <- Seurat::CreateSeuratObject(raw.data = raw_data, project = "SeuratProject", min.cells = min_cells,
                                           min.genes = min_genes, is.expr = 0, do.center = FALSE, meta.data = meta_data)
  print('initial dims after filtering for min 20 cells per gene, min 200 genes per cell')
  dim(seurat_obj@raw.data)
  print('running Seurat Pipeline')
  print(Sys.time())
  seurat_obj <- seurat_pipeline(seurat_obj)
  seurat_obj <- Seurat::RunTSNE(seurat_obj)
  # seurat_obj <- Seurat::RunDiffusion(seurat_obj)
  print(Sys.time())
  print('Clustering')
  seurat_obj <- Seurat::FindClusters(seurat_obj, print.output = FALSE)
  return(seurat_obj)
}

scoring_subroutine <- function(seurat_obj, genesets, preproc_data_dir, prefix, min_length = 3) {
  print('Filtering Genesets')
  orig_length <- unlist(lapply(genesets, length))
  orig_gensets_condensed <- unlist(lapply(genesets, FUN = function(x) {paste(x, collapse = ';')}))
  genesets_filt <- lapply(genesets, FUN = function(x) {return(x[x %in% rownames(seurat_obj@raw.data)])})
  filt_length <- unlist(lapply(genesets_filt, length))
  filt_gensets_condensed <- unlist(lapply(genesets_filt, FUN = function(x) {paste(x, collapse = ';')}))
  kept <- filt_length >= min_length
  genesets_summary <- data.frame(geneset_name = names(genesets), 
                                 original_length = orig_length, 
                                 original_genes = orig_gensets_condensed,
                                 filtered_length = filt_length, 
                                 filtered_genes = filt_gensets_condensed,
                                 kept = kept)
  genesets_summary_file <- paste0(prefix, '_genesets_summary.csv')
  write.csv(genesets_summary, file = file.path(preproc_data_dir, genesets_summary_file))
  if (any(kept)) {
    genesets_filt <- genesets_filt[which(kept)]
  } else {
    stop('0 genesets kept as none with >= 5 genes')
  }
  print('Running Scoring')
  print(Sys.time())
  # Run Scoring
  seurat_obj <- AUCell_batch(inp_data = seurat_obj, genesets = genesets_filt, num_batches = 1)
  print(Sys.time())
  return(seurat_obj)
}

filter_df <- function(df, filter_by = NULL, filter_class = NULL) {
  ## filter by value for a particular class
  ## INPUTS:
  ##    df = data.frame
  ##    filter_by = NULL or string denoting class by which data is to be subsetted
  ##    filter_class = set of valid classes in filter_by categorical variable to use in subset
  ##    Note that if filter_by is set to NULL, no filtering is done
  ##    If filter_class is left unspecified and filter_by is specified as a valid column name,
  ##    an error will be thrown
  ## RETURNS:
  ##  either the same data.frame or if valid non-NULL arguments provided to filter_by and
  ##  filter_class, a subset of that data consisting of all members of the valid classes
  ##  specified in filter_class
  if (is.character(filter_by)) {
    if (is.character(filter_class)) {
      filter_vect <- df[,filter_by]
      if (is.numeric(filter_vect)) {
        warning('specified filter feature is of numeric value, coercing to character')
      }
      filter_vect <- as.character(filter_vect)
      if (all(filter_class %in% filter_vect)) {
        valid_cell_inds <- filter_vect %in% filter_class
        df <- df[which(valid_cell_inds),]
      } else {
        stop('Not all specified classes are present in the specified filter_by column')
      }
    } else {
      warning('filter_class left unspecified, continuing to plot without filtering data')
    }
  } else if (!is.null(filter_by)) {
    warning('filter_by must be specified as string (character vector of length 1) if non-NULL\
    in order to filter. Ignoring argument')
  }
  return(df)
}

get_lims <- function(x) {
  ## Get axis limits for a set of points
  ## INPUTS:
  ##  x = numeric vector
  ## RETURNS:
  range_x <- max(x) - min(x)
  buff <- 0.05*range_x
  return(c(min(x) - buff, max(x) + buff))
}

get_density_grid <- function(df, var1, var2, eps = 3.0, grid_dist = 1.0, density_ceiling = 0.25) {
  ## INPUTS:
  ##  df = data.frame
  ##  var1 = character, denotes variable for x coordinates
  ##  var2 = character, denotes variable for y coordiantes
  ##  eps = radius within which to include 'neighbors' for calculating density.
  ##        density is calculated as number of neighbors in a given radius divided by area of circle
  ##        of radius 'eps'
  ##  grid_dist = distance to increment when making grid. lower values = more points = longer run time
  ##  density_ceiling = percentile density at which to apply ceiling
  ## RETURNS:
  ##  a data.frame with x and y coordinates in grid and their corresponding density values.
  ##  number of grid coordinates is ~ ((xmax-xmin)/grid_dist)*((ymax-ymin)/grid_dist)
  vars_get <- unique(c(var1, var2))
  df <- df[,vars_get]
  ## make a grid of points
  x_lims <- get_lims(df[,var1])
  y_lims <- get_lims(df[,var2])
  x_seq <- seq(x_lims[1], x_lims[2], grid_dist)
  y_seq <- seq(y_lims[1], y_lims[2], grid_dist)
  grid_coords <- c()
  for (i in 1:length(x_seq)) {
    grid_coords <- rbind(grid_coords, cbind(rep(x_seq[i], length(y_seq)), y_seq))
  }
  grid_density <- apply(grid_coords, MARGIN = 1, FUN = function(x) {
    # calculate pairwise distances for each point in grid (row)
    # to points in actual data. Take density as # of cells within distance eps/area
    x_rep <- matrix(rep(x, nrow(df)), ncol = 2, byrow = TRUE)
    diff_mat <- x_rep - as.matrix(df[,c(var1, var2)])
    distances <- sqrt(rowSums(diff_mat*diff_mat))
    n_points <- length(which(distances < eps))
    circle_area <- pi*(eps^2)
    return(n_points/circle_area)
  })
  grid_coords <- as.data.frame(grid_coords)
  colnames(grid_coords) <- c(var1, var2)
  grid_coords$density <- grid_density
  if (!is.null(density_ceiling)) {
    if (is.numeric(density_ceiling)) {
      if (0 < density_ceiling && density_ceiling < 1.0) {
        pct_ceiling <- quantile(grid_coords$density, probs = density_ceiling)
        grid_coords$density[grid_coords$density > pct_ceiling] <- pct_ceiling
      } else {
        stop('density_ceiling must indicate a value between 0 and 1.0, corresponding to a percentile at which to cut off density')
      }
    } else {
      stop('non numeric argument for density_ceiling')
    }
  }
  return(grid_coords)
}

discrete_color_mapping <- function(categorical_vals) {
  unique_vals <- unique(categorical_vals) 
  unique_vals <- unique_vals[order(unique_vals)]
  col_map <- hue_pal()(length(unique_vals))
  names(col_map) <- unique_vals
  return(col_map)
}

make_contour_plot <- function(inp_data, var1 = NULL, var2 = NULL, split_by = NULL, 
                              dim1 = NULL, dim2 = NULL, reduction_use = 'pca', eps = 3.0, 
                              grid_dist = 1.0, bins = 10, density_ceiling = 0.25, 
                              legend_pt_size = 5) {
  ## INPUTS:
  ##  inp_data = data.frame or Seurat object
  ##  split_by = categorical variable by which to split data. if not NULL,
  ##            will make 3 separate contour plots for each class in the
  ##            categorical variable
  ##  var1, var2, epse, grid_dist, density_ceiling: see get_density_grid
  ##  dim1, dim2 = dimensions to use from reduced dims if using Seurat input
  ##  
  ##  bins = number of bins for contour plot (passed to geom_contour)
  ##  legend_pt_size = size of points in legend
  ## RETURNS:
  ##  ggplot2 object with contours plotted for all datapoints.
  if (class(inp_data) == 'seurat') {
    if (!is.null(var1) || !is.null(var2)) {
      warning('ignoring var1 and var2 for seurat input')
    }
    if (!is.numeric(dim1) || !(is.numeric(dim2))) {
      stop('dim1 or dim2 must be specified to integer value')
    }
    embeddings <- Seurat::GetCellEmbeddings(inp_data, reduction.type = reduction_use)[,c(dim1, dim2)]
    fetched_data <- Seurat::FetchData(inp_data, vars.all = c(split_by))
    input_df <- cbind(as.data.frame(embeddings), fetched_data)
    p <- make_contour_plot(inp_data = input_df, var1 = colnames(embeddings)[1], var2 = colnames(embeddings)[2], split_by = split_by, 
                           eps = eps, grid_dist = grid_dist, bins = bins, density_ceiling = density_ceiling, 
                           legend_pt_size = legend_pt_size)
  } else if (is.data.frame(inp_data)) {
    if (!is.character(var1) || !is.character(var2)) {
      stop('for data.frame input, var1 and var2 must be specified')
    }
    if (is.null(split_by)) {
      # just do contour plot using density computed using all datapoints
      grid_coords <- get_density_grid(inp_data, var1, var2, eps = eps, 
                                      grid_dist = grid_dist, density_ceiling = density_ceiling)
      p <- ggplot(data = grid_coords, mapping = aes_string(x = var1, y = var2, z = 'density'))
      p <- p + geom_contour(bins = bins)
    } else {
      # separately make a contour plot for each class in categorical variable 'split_by'
      if (is.numeric(inp_data[,split_by])) {
        stop(paste('argument \'split_by\' should correspond to a categorical variable'))
      }
      discrete_col_mapping <- discrete_color_mapping(inp_data[,split_by])
      for (i in 1:length(discrete_col_mapping)) {
        class_i <- names(discrete_col_mapping)[i]
        df_subs <- filter_df(inp_data, filter_by = split_by, filter_class = class_i)
        grid_coords_subs <- get_density_grid(df_subs, var1, var2, eps = eps, 
                                             grid_dist = grid_dist, density_ceiling = density_ceiling)
        grid_coords_subs[,split_by] <- rep(class_i, nrow(grid_coords_subs))
        if (i == 1) {
          p <- ggplot(data = grid_coords_subs, mapping = aes_string(x = var1, y = var2, z = 'density', color = split_by))
          p <- p + geom_contour(bins = bins)
        } else {
          p <- p + geom_contour(data = grid_coords_subs, mapping = aes_string(x = var1, y = var2, z = 'density', color = split_by), 
                                bins = bins)
        }
      }
      p <- p + guides(color = guide_legend(override.aes = list(color = discrete_col_mapping,
                                                               size = legend_pt_size,
                                                               shape = 15,
                                                               labels = names(discrete_col_mapping))))
    }
  } else {
    stop(paste('unrecognized argument of class', class(inp_data), 'for inp_data'))
  }
  return(p)
}

scatter_generic <- function(df, var1, var2, color_by, 
                            filter_by = NULL, filter_class = NULL,
                            pt.size = 0.5, pt.shape =  16, alpha = 0.4, 
                            legend_pt_size = 20, legend_pt_shape = 15,
                            xlim = NULL, ylim = NULL, legend_title = FALSE) {
  ## generic scatterplot function for data frame
  ## INPUTS:
  ##  df = data.frame
  ##  var1 = x axis variable
  ##  var2 = y axis variable
  ##  color_by = variable to color points by
  ##    filter_by = subset data by a particular categorical variable
  ##    filter_class = set of allowed classes for subset
  ##    pt.size, alpha, legend_pt_size, legend_pt_shape, xlim, ylim: arguments
  ##    passed to ggplot2 plotting functions. Note that if xlim and ylim are set to null,
  ##    by default the limits are taken as (min(x) - 0.05*range(x), max(x) + 0.05*range(x))
  ##    for the x limits and the same for the y limits but replacing x with y.
  ##    the limits are calculated with the whole set of points prior to filtering
  ##    Note that argument alpha is ignored if we have a continuous scale AND the range
  ##    of values includes negative and positive values (i.e crosses 0)
  ## RETURNS:
  ##  ggplot2 object
  # set x and y limits if NULL
  
  if (is.null(xlim)) {
    xlim <- get_lims(df[,var1])
  }
  if (is.null(ylim)) {
    ylim <- get_lims(df[,var2])
  }
  
  tryCatch({df[,color_by]},
           error = function(e) {
             # since warnings not going to stderr file for Rmd, if feature is not found in
             # df we specify an error that says the name of the column not found
             stop(paste('could not find feature', color_by, 'in metadata'))
           })
  # decide if dealing with categorical data, set appropriate settings
  is_categorical <- !is.numeric(df[,color_by])
  
  # make plot
  if (is_categorical) {
    # categorical values
    discrete_col_mapping <- discrete_color_mapping(df[,color_by])
    # filter data as desired after getting appropriate color mapping
    df <- filter_df(df = df, filter_by = filter_by, filter_class = filter_class)
    p <- ggplot(data = df, mapping = aes_(x = as.name(var1), y = as.name(var2))) +
      geom_point(aes_(color = as.name(color_by)), size = pt.size, alpha = alpha, pch = pt.shape)
    p <- p + scale_color_manual(values = discrete_col_mapping)
    p <- p + guides(colour = guide_legend(override.aes = list(size=legend_pt_size,
                                                              shape=legend_pt_shape),
                                          title = ifelse(legend_title, color_by, '')))
  } else {
    # numeric values
    qtl_vals <- quantile(df[,color_by], seq(0,1, 0.05))
    # set values for color scale
    min_val <- qtl_vals['5%']
    max_val <- qtl_vals['95%']
    if (sign(min_val) == -sign(max_val)) {
      # continuous_colors <- c('blue', 'white', 'red')
      continuous_breaks <- c(qtl_vals['5%'], 0, qtl_vals['95%'])
      zero_center <- TRUE
    } else {
      continuous_colors <- c('red', 'yellow')
      continuous_breaks <- c(qtl_vals['5%'], qtl_vals['95%'])
      zero_center <- FALSE
    }
    # set color values
    df$color_vals <- df[,color_by]
    df[df[,color_by] >= max_val, 'color_vals'] <- max_val
    df[df[,color_by] <= min_val, 'color_vals'] <- min_val
    # if (zero_center) {
    #   # set alpha values to absolute value of feature of interest, using top 5%/bottom 25% as
    #   # max/min values for alpha, then scale to be between 0 and 1
    #   alpha_vals <- abs(df[,color_by])
    #   alpha_quantiles <- quantile(alpha_vals, seq(0,1, 0.05))
    #   pct5_alpha <- alpha_quantiles['5%']
    #   pct95_alpha <- alpha_quantiles['95%']
    #   alpha_vals[alpha_vals > pct95_alpha] <- pct95_alpha
    #   alpha_vals[alpha_vals < pct5_alpha] <- pct5_alpha
    #   alpha_vals <- min(alpha_vals/pct95_alpha, 0.5*pct95_alpha)
    #   df$alpha_vals <- alpha_vals
    # }
    # filter data as desired after deciding appropriate color mapping
    df <- filter_df(df = df, filter_by = filter_by, filter_class = filter_class)
    # do the plot
    p <- ggplot(data = df, mapping = aes_(x = as.name(var1), y = as.name(var2)))
    if (zero_center) {
      p <- p + geom_point(aes_(color = df$color_vals), size = pt.size, alpha = alpha, pch = pt.shape)
      # p <- p + geom_point(aes_(color = df$color_vals), size = pt.size, alpha = alpha, pch = pt.shape)
      # p <- p + scale_color_gradient2(midpoint = 0, limits = c(min_val, max_val), low = 'turquoise', mid = 'white', high = 'orange')
      p <- p + scale_color_gradient2(midpoint = 0, low = 'blue', mid = 'darkgrey', high = 'red')
      # p <- p + theme_dark()
      p <- p + theme_light()
    } else {
      p <- p + geom_point(aes_(color = df$color_vals), size = pt.size, alpha = alpha, pch = pt.shape)
      # p <- p + scale_color_gradient2(low = 'turquoise', mid = 'white', high = 'orange',
      #                                midpoint = mean(df$color_vals),
      #                                limits = c(min_val, max_val))
      p <- p + scale_color_gradient2(low = 'darkblue', mid = 'darkgrey', high = 'red', midpoint = mean(df$color_vals))
      # p <- p + theme_dark()
      p <- p + theme_light()
    }
    if (!legend_title) {
      p <- p + guides(color = guide_colorbar(title = NULL))
    } else {
      p <- p + guides(color = guide_colorbar(title = color_by))
    }
    # p <- p + scale_color_gradientn(colors = continuous_colors, breaks = continuous_breaks)
    
  }
  p <- p + coord_cartesian(xlim = xlim, ylim = ylim)
  
}

single_dr_plot <- function(obj, dim1, dim2, reduction.use, feat,
                           filter_by = NULL, filter_class = NULL,
                           pt.size = 0.5, pt.shape = 16, alpha = 0.4, 
                           legend_pt_size = 20, legend_pt_shape = 15,
                           xlim = NULL, ylim = NULL, legend_title = FALSE) {
  ## make dimensionality reduction plot for Seurat given a dim reduction to use,
  ##  numeric features, categorical features, and specified dimensions
  ## INPUTS: 
  ##    obj = seurat object
  ##    dim1, dim2 = dimensions from reduced dims to plot
  ##    reduction.use = reduced dims to use. assume they're already calculated
  ##    feat = feature to plot
  ##    ... other arguments passed to scatter_generic
  ##  RETURNS: 
  ##    ggplot2 object
  # get reduced dims
  dim_red_emb <- 	Seurat::GetCellEmbeddings(object = obj, reduction.type = reduction.use)
  dim_red_emb <- dim_red_emb[,c(dim1, dim2)]
  dim_red_names <- colnames(dim_red_emb)
  dim_red_emb <- as.data.frame(dim_red_emb)
  # get data frame of metadata
  feats_fetch <- feat
  if (any(feats_fetch %in% dim_red_names)) {
    # if reduced dims in numeric features to plot, remove
    # as these are already fetched and will be added
    # to plot_data
    feats_fetch <- setdiff(feats_fetch, dim_red_names)
  }
  if (is.character(filter_by)) {
    feats_fetch <- union(feats_fetch, filter_by)
  }
  if (length(feats_fetch) > 0) {
    not_found_warnings <- capture.output(plot_data <- Seurat::FetchData(obj, vars.all = feats_fetch), type = 'message')
    if (is.matrix(plot_data)) {
      plot_data <- as.data.frame(plot_data)
    } else if (!is.data.frame(plot_data)) {
      stop('Expect Seurat::FetchData to produced data.frame or matrix with row and column names')
    }
    stopifnot(all(rownames(dim_red_emb) == rownames(plot_data)))
    plot_data <- cbind(dim_red_emb, plot_data)
    if (length(not_found_warnings) >= 1) {
      for(msg in 1:length(not_found_warnings)) {
        warning(msg)
      }
    }
  } else {
    plot_data <- dim_red_emb
  }
  p <- scatter_generic(df = plot_data, 
                      var1 = names(dim_red_emb)[1], 
                      var2 = names(dim_red_emb)[2],
                      color_by = feat,
                      filter_by = filter_by, 
                      filter_class = filter_class,
                      pt.size = pt.size, 
                      alpha = alpha, 
                      legend_pt_size = legend_pt_size, 
                      legend_pt_shape = legend_pt_shape,
                      xlim = xlim, 
                      ylim = ylim, 
                      legend_title = legend_title) 
  p <- p + ggtitle(feat)
  return(p)
}

make_dr_plots <- function(obj, dim1, dim2, reduction.use, feats_plot, return_plotlist = FALSE,
                          filter_by = NULL, filter_class = NULL, pt.size = 0.5, pt.shape = 16, alpha = 0.4, 
                          legend_pt_size = 20, legend_pt_shape = 15,
                          xlim = NULL, ylim = NULL, legend_title = FALSE) {
  ## Make multiple dimensionality reduction plots
  ##  INPUTS:
  ##    obj = Seurat object
  ##    dim1, dim2 = dimensions from reduced dims to plot
  ##    reduction.use = type of reduced dims to use, assume they're already calculated
  ##    feats_plot = character vector of features to plot
  ##    return_plotlist = return list of ggplot2 objects instead of printing plots
  ##    ... = arguments passed to scatter_generic
  ##  RETURNS:
  ##    either prints all plots or returns them as a list of ggplot2 objects
  plotlist <- list()
  for (feat in feats_plot) {
    p <- single_dr_plot(obj = obj, dim1 = dim1, dim2 = dim2, 
                        reduction.use = reduction.use, 
                        feat = feat, 
                        filter_by = filter_by, 
                        filter_class = filter_class,
                        pt.size = pt.size, 
                        pt.shape = pt.shape,
                        alpha = alpha, 
                        legend_pt_size = legend_pt_size, 
                        legend_pt_shape = legend_pt_shape,
                        xlim = xlim, 
                        ylim = ylim, 
                        legend_title = legend_title)
    if (return_plotlist) {
      plotlist[[feat]] <- p
    } else {
      print(p)
    }
  }
  if (return_plotlist) {
    return(plotlist)
  }
}

corplot_2_var <- function(obj, var1, var2, color_by, slot = 'data', method = 'spearman', 
                          filter_by = NULL, filter_class = NULL, pt.size = 1, alpha = 0.4, 
                          legend_pt_size = 20, legend_pt_shape = 15,
                          xlim = NULL, ylim = NULL, legend_title = TRUE) {
  ## Show correlation of two numeric variables for either a dataframe or seurat object
  ## INPUTS:
  ##  obj = Seurat object or data.frame
  ##  var1, var2 = variables to plot on x and y axis, respectively
  ##  slot = slot to pull any gene expression data from if input is a Seurat object
  ##  pt_size = point size for plot
  if (is(obj, 'seurat')) {
    vars_fetch <- union(c(var1, var2), color_by)
    vars_fetch <- union(vars_fetch, filter_by)
    df <- Seurat::FetchData(obj, vars.all = vars_fetch)
    p <- corplot_2_var(obj = df, var1 = var1, var2 = var2, color_by = color_by, 
                       slot = slot, method = method, filter_by = filter_by, 
                       filter_class = filter_class, pt.size = pt.size, alpha = alpha, 
                       legend_pt_size = legend_pt_size, legend_pt_shape = legend_pt_shape, 
                       xlim = xlim, ylim = ylim,
                       legend_title = legend_title)
    return(p)
  } else if (is.data.frame(obj)) {
    df_subs <- filter_df(obj, filter_by = filter_by, filter_class = filter_class)
    # calculate correlation
    cor_val <- cor(df_subs[,var1], df_subs[,var2], method = method)
    cor_val <- round(cor_val, 2)
    cor_pval <- cor.test(df_subs[,var1], df_subs[,var2], method = method)$p.value
    if (cor_pval < .Machine$double.eps) {
      cor_pval <- paste0(' < ', signif(.Machine$double.eps, 2))
    } else {
      cor_pval <- paste0(' = ', signif(cor_pval, 2))
    }
    # make plot
    p <- scatter_generic(df = obj, 
                         var1 = var1, 
                         var2 = var2, 
                         color_by = color_by,
                         filter_by = filter_by, 
                         filter_class = filter_class,
                         pt.size = pt.size, 
                         alpha = alpha, 
                         legend_pt_size = legend_pt_size, 
                         legend_pt_shape = legend_pt_shape,
                         xlim = xlim, 
                         ylim = ylim, 
                         legend_title = legend_title)
    title_use <- paste(var1, 'vs', var2)
    subtitle_use <- paste(method, 'correlation:', cor_val, '; p', cor_pval)
    p <- p + labs(title = title_use, subtitle = subtitle_use)
    return(p)
  } else {
    stop(paste('unexpected argument obj of class', class(obj)))
  }
}

# make function for stacked bar charts
stacked_barcharts <- function(inp_data, split_by, color_by, title = '', show_proportion = FALSE, filter_by = NULL, filter_class = NULL) {
  inp_data <- filter_df(inp_data, filter_by = filter_by, filter_class = filter_class)
  if (show_proportion) {
    # show stacked bar chart with proportions per class
    p <- ggplot(data = inp_data, aes_(x = as.name(split_by)))
    p <- p + geom_bar(aes_(fill = as.name(color_by)), position = position_fill()) + theme(axis.text.x = element_text(size = 7)) +
      ylab('proportion') + ggtitle(title)
  } else {
    # show stacked bar chart with counts per class
    p <- ggplot(data = inp_data, aes_(x = as.name(split_by)))
    p <- p + geom_bar(aes_(fill = as.name(color_by))) + theme(axis.text.x = element_text(size = 7)) + 
      ggtitle(title)
  }
  return(p)
}
# show numeric features split by a group, and potentially further split into subgroups
sc_jitter_plot <- function(seurat_obj, y_axis, split_by = 'SampleType', color_by = 'SampleType', 
                          filter_by = NULL, filter_class = NULL, comparisons = NULL, mode = 'jitter', size = 2.0, 
                           dodge_width = 0.75, jitter_width = 0.2) {
  
  # seurat_obj = seurat object
  # y_axis = variable to plot on y axis
  # split_by = groups to split by
  # color_by = same as split_by if you want to compare split by groups, different if you want to do within group comparisons
  # filter_by = subset data by this particular variable
  # filter_class = class filtered for within filter_by variable
  # comparisons = 'all', or list of comparisons you'd like to perform. if color_by != split by, must set to all, and 
  #   color_by variavle must have exactly 2 classes. NULL if you don't want any comparisons. will add wilcoxon rank sum p-values to plot.
  #   note that this does not address the issue of multiple testing, so if you want to be certain about something, 
  #   run a identical wilcoxon ranksum test for all features of interest, and adjust p-values as desired.
  # mode = jitter for jitter plot, violin for violin plot, jitter-dodge for violin plot with jittered points,
  #   boxplot for boxplot
  # returns: ggplot2 object with desired plot.
  
  if (is.null(color_by)) {
    color_by <- split_by
  }
  if (!is.null(filter_by)) {
    input_df <- Seurat::FetchData(seurat_obj, vars.all = unique(c(split_by, color_by, y_axis, filter_by)))
    input_df <- filter_df(df = input_df, filter_by = filter_by, filter_class = filter_class)
  } else {
    input_df <- Seurat::FetchData(seurat_obj, vars.all = unique(c(split_by, color_by, y_axis)))
  }
  stopifnot(is.character(input_df[,split_by]) || is.factor(input_df[,split_by])
            || is.logical(input_df[,split_by]))
  if (is.logical(input_df[,split_by])) {
    input_df[,split_by] <- factor(input_df[,split_by])
  }
  if (!is.numeric(input_df[,y_axis])) {
    stop(paste('could not make plot as variable for y_axis', y_axis, 'is non-numeric'))
  }
  p <- ggplot(data = input_df, aes_(as.name(split_by), as.name(y_axis)))
  if (mode == 'jitter') {
    p <- p + geom_jitter(aes_(color = as.name(color_by)), size = size)
  } else if (mode == 'violin') {
    p <- p + geom_violin(aes_(fill = as.name(color_by)), size = size)
  } else if (mode == 'jitter-dodge') {
    p <- p + geom_violin(aes_(color = as.name(color_by)), position = position_dodge(width = dodge_width)) + 
      geom_point(aes_(color = as.name(color_by)), pch = 21, 
                 position = position_jitterdodge(jitter.width = jitter_width, dodge.width = dodge_width))
    
  } else if (mode == 'boxplot') {
    p <- p + geom_boxplot(aes_(fill = as.name(color_by)))
  } else {
    stop('unexpected mode argument')
  }
  p <- p + ylab(y_axis) + xlab(split_by) + guides(fill = guide_legend(title = color_by))
  
  if (!is.null(comparisons)) {
    tryCatch({
      if (is.list(comparisons)) {
        stopifnot(all(as.logical(lapply(comparisons, FUN = function(x) {is.character(x) && length(x) == 2}))))
      } else if (is.character('comparisons')) {
        stopifnot(length(comparisons) == 1 && comparisons == 'all')
      }
    }, error = function(e) {
      stop('comparisons must be NULL, \'all\', or a list of character vectors of length 2')
    })
    
    # internal function to map p-values to dots
    pval_2_dot <- function(x) {
      if (is.na(x)) {
        return('NA')
      }
      if (x < 0.001) {
        return('***')
      } else if (x < 0.01) {
        return('**')
      } else if (x < 0.05) {
        return('*')
      } else {
        return('NS')
      }
    }
    
    if (color_by == split_by) {
      # comparisons should compare groups among split_by variable
      levels_use <- unique(input_df[,split_by])
      if (length(levels_use) < 2) {
        stop('only one class detected in split_by, comparisons should be set to NULL')
      }
      if (!is.list(comparisons)) {
        comparisons <- list()
        for (i in 1:(length(levels_use) - 1)) {
          for (j in (i+1):length(levels_use)) {
            level_i <- as.character(levels_use[i])
            level_j <- as.character(levels_use[j])
            compar_name <- paste0(level_i, '_vs_', level_j)
            comparisons[[compar_name]] <- c(level_i, level_j)
          }
        }
      }
      stopifnot(is.list(comparisons))
      ext_plot_data <- data.frame(grp = character(0), height_bar = numeric(0), dots = character(0))
      for (compar in comparisons) {
        gr_a <- as.character(compar[1])
        gr_b <- as.character(compar[2])
        bad_levels <- setdiff(c(gr_a, gr_b), levels_use)
        if (length(bad_levels) > 0) {
          stop(paste('found', length(bad_levels), 'categories not in split_by:', paste(bad_levels, collapse = ';')))
        }
        gr_a_vals <- input_df[ input_df[,split_by] == gr_a , y_axis ]
        gr_b_vals <- input_df[ input_df[,split_by] == gr_b , y_axis ]
        p_val <- wilcox.test(gr_a_vals, gr_b_vals, alternative = 'two.sided')$p.value
        pct_95 <- quantile(c(gr_a_vals, gr_b_vals), probs = seq(0, 1, 0.05))['95%']
        height_bar <- 0.01 + 2*pct_95
        height_dots <- height_bar*1.10
        dot_char <- pval_2_dot(p_val)
        ext_plot_data <- rbind(ext_plot_data, 
                               data.frame(grp = c(gr_a, gr_b), 
                                          height_bar = rep(height_bar, 2),
                                          height_dots = rep(height_dots, 2), 
                                          dots = c(dot_char, ''),
                                          stringsAsFactors = FALSE)
                               )
      }
      ext_plot_data$grp <- factor(ext_plot_data$grp)
      p <- p + geom_path(data = ext_plot_data, mapping = aes(grp, y = height_bar, group = 1))
      p <- p + geom_text(data = ext_plot_data, mapping = aes(grp, y = height_dots, label = dots), nudge_x = 0.5, color = 'red')
    } else {
      # pairwise comparisons between color by classes within groups specified by split_by
      levels_use <- unique(input_df[,color_by])
      if (length(levels_use) > 2) {
        stop('only pairwise comparisons supported if color_by != split_by (i.e. for comparisons within split_by groups)')
      }
      if (comparisons != 'all') {
        stop('comparisons should be set to \'all\' if color_by != split_by')
      }
      split_grs <- unique(input_df[,split_by])
      compar_grs <- unique(input_df[,color_by])
      ext_plot_data <- data.frame(split_gr_id = character(0), dot_height = numeric(0), dot_char = character(0))
      for (split_gr_id in split_grs) {
        data_subs <- input_df[ input_df[,split_by] == split_gr_id , ]
        if (!all(compar_grs %in% data_subs[,color_by])) {
          dot_char <- 'NA'
          dot_height <-  0.01 + 2*quantile(data_subs[,y_axis], probs = seq(0, 1, 0.05))['95%']
        } else {
          gr_a <- levels_use[1]
          gr_b <- levels_use[2]
          gr_a_vals <- data_subs[ data_subs[,color_by] == gr_a , y_axis]
          gr_b_vals <- data_subs[ data_subs[,color_by] == gr_b , y_axis]
          p_val <- wilcox.test(gr_a_vals, gr_b_vals, alternative = 'two.sided')$p.value
          dot_char <- pval_2_dot(p_val)
          dot_height <- 0.01 + 2*quantile(data_subs[,y_axis], probs = seq(0, 1, 0.05))['95%']
        }
        ext_plot_data <- rbind(ext_plot_data, 
                               data.frame(split_gr_id = split_gr_id,
                                          dot_height = dot_height,
                                          dot_char = dot_char,
                                          stringsAsFactors = FALSE)
                               )
      }
      ext_plot_data$split_gr_id <- factor(ext_plot_data$split_gr_id)
      p <- p + geom_text(data = ext_plot_data, aes(split_gr_id, dot_height, label = dot_char), color = 'red') 
    }
    # legend for p-values
    p <- p + ggtitle('Wilcoxon RS Test', sub = '*** = p < 0.001\n** = p < 0.01\n* = p < 0.05\nNS = Not Significant\tNA = no p-val')
  }
  
  return(p)
}

run_fisher_pathways <- function(genesets, GSA_gmt, master.set, FDR = 0.1, filter.by.fdr = TRUE) {
  # INPUTS:
  #   genesets = named list of genesets
  #   GSA_gmt: output from GSA.read.gmt (contains list of pathway genes, list of pathway names)
  # RETURNS:
  #   data frame with fisher's exact test result with columns (geneset, pathway, FDR, passes_fdr)
  # run fisher's exact test on all genesets for all pathways, storing adjusted p-values
  # in a matrix. p-values adjsuted among all tests corresponding to same geneset
  result_mat <- c()
  FDR <- 0.10
  for (gs_name in names(genesets)) {
    gs_i <- genesets[[gs_name]]
    p_vals <- as.numeric(lapply(GSA_gmt$genesets, FUN = function(x) {
      test_res <- fisher.test.sets(set1 = gs_i, set2 = x, master.set = master.set)
      return(test_res$p.value)
    }))
    adj_pvals <- p.adjust(p_vals, method = 'BH')
    result_mat <- rbind(result_mat, adj_pvals)
  }
  rownames(result_mat) <- names(genesets)
  colnames(result_mat) <- GSA_gmt$geneset.names
  result_df <- melt(result_mat, value.name = 'FDR', varnames = c('geneset', 'pathway'))
  pathway_description <- GSA_gmt$geneset.descriptions
  names(pathway_description) <- GSA_gmt$geneset.names
  result_df$description <- pathway_description[result_df$pathway]
  result_df$passes_fdr <- as.logical(result_df$FDR <= FDR)
  if (filter.by.fdr) {
    result_df <- result_df[which(result_df$passes_fdr),]
  }
  return(result_df)
}

get_de_pathways <- function(obj, ident.1, ident.2, class, gene.fwer = 0.05, pathway.fdr = 0.1, GSA_gmt, ...) {
  ## given a seurat object, run differential expression between ident.1 and ident.2 for given class,
  ## and run pathway enrichment analysis on upreg and downreg genes
  meta_data <- obj@meta.data
  success <- FALSE
  try ({
    ident_vect <- meta_data[,class]
    success <- TRUE
  })
  if (!success) {
    stop(paste('could not find class', class, 'in metadata columns'))
  }
  obj <- Seurat::SetIdent(obj, ident.use = ident_vect)
  markers_df <- Seurat::FindMarkers(object = obj, ident.1 = ident.1, ident.2 = ident.2, ...)
  markers_df <- markers_df[markers_df$p_val_adj <= gene.fwer,]
  markers_df <- markers_df[order(markers_df$avg_logFC, decreasing = TRUE), ]
  gene_names <- rownames(markers_df)
  upreg_genes <- gene_names[markers_df$avg_logFC > 0]
  downreg_genes <- gene_names[markers_df$avg_logFC < 0]
  all.genes <- rownames(GetAssayData(obj, assay = 'RNA', slot = 'data'))
  upreg_pathways <- run_fisher_pathways(list(upreg_genes = upreg_genes), GSA_gmt, 
                                        master.set = all.genes, FDR = pathway.fdr,
                                        filter.by.fdr = TRUE)
  downreg_pathways <- run_fisher_pathways(list(downreg_genes = downreg_genes), GSA_gmt, 
                                          master.set = all.genes, FDR = pathway.fdr,
                                          filter.by.fdr = TRUE)
  return(list(diff_exp_results = markers_df, upreg_genes = upreg_genes, upreg_pathways = upreg_pathways,
              downreg_genes = downreg_genes, downreg_pathways = downreg_pathways))
}


hm_w_numbers <- function(num.mat, color.mat, color.value.name, order.samples, text_size = 2.0,
                         axis_label_size = 3.0) {
  num.mat <- num.mat[order.samples, order.samples]
  color.mat <- color.mat[order.samples, order.samples]
  mat.melt <- melt(num.mat, value.name = 'num.value')
  color.melt <- melt(color.mat, value.name = color.value.name)
  p <- ggplot(data = color.melt, aes_(x=as.name('Var1'), y=as.name('Var2'), fill=as.name(color.value.name)))
  p <- p +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = axis_label_size),
          axis.text.y = element_text(angle = 0, hjust = 1, size = axis_label_size)) +
    geom_tile(data = color.melt) + scale_fill_gradient(low = 'white', high = 'red') +
    geom_text(aes_(label = round(mat.melt[,'num.value'], 2)), size = text_size)
  return(p)
}
#
do_pickle <- function(output_dir, outfile, inp_data) {
  # convenience function for pickling python data
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  outfile_path <- file.path(output_dir, outfile)
  reticulate::py_save_object(inp_data, outfile_path, pickle = 'pickle')
  return(outfile_path)
}

fit_logistic <- function(X, y, random_state = 32L, penalty = 'l2') {
  ## fit logistic regression model with CV, outputting best model
  ## INPUTS:
  ##  X = input data
  ##  y = vector of outputs
  ##  random_state = random seed
  ##  penalty = penalty used when fitting model
  ## RETURNS:
  ##  list of class logistic_cv_result which contains
  ##  the model generated by LogisticRegressionCV, the average score (accuracy)
  ##  obtained for each lambda, and features with non-zero coefficients
  output.model <- LM$LogisticRegressionCV(cv = 5L,
                                          Cs = 50L,
                                          solver = switch(penalty,
                                                          l1 = 'liblinear', 
                                                          l2 = 'sag', 
                                                          elasticnet = 'saga'),
                                          penalty = penalty,
                                          random_state = random_state,
                                          l1_ratios = 1:100/100)

  output.model$fit(X = X, y = y)
  avg.score <- lapply(output.model$scores_, FUN = function(x) {
    apply(x, MARGIN = 2, mean)
  })
  n_class <- length(unique(y))
  coefs <- output.model$coef_
  if (n_class == 1) {
    coefs <- as.numeric(coefs)
    names(coefs) <- colnames(X)
    selected.features <- names(coefs)[which(abs(coefs) > 0)]
  } else {
   for (j in 1:nrow(coefs)) {
     coefs_j <- coefs[j,]
     names(coefs_j) <- colnames(X)
     nonzero_coefs_j <- names(coefs_j)[which(abs(coefs_j) > 0)]
     if (j == 1) {
       selected.features <- nonzero_coefs_j
     } else {
       selcted.features <- c(selected.features, nonzero_coefs_j)
     }
   } 
  }

  # check for non-zero coefs. these are selected features
  
  selected.features <- selected.features[order(selected.features)]
  
  # if (write_file) {
  #  outfile_path <- do_pickle(output_dir, outfile, output.model) 
  # } else {
  #   outfile_path <- NULL
  # }
  logistic_cv_result = list(model = output.model,
                            avg.score = avg.score,
                            selected.features = selected.features)
  class(logistic_cv_result) <- 'logistic_cv_result'
  
  return(logistic_cv_result)
}

# plot_cv_avg <- function(logistic_cv_result, class.use = 1L, xlab = 'log10(C)', ylab = 'score') {
#   # class.use = either integer for indexing avg.score list
#   # or a string indicating name of class.
#   c.vals <- as.numeric(logistic_cv_result$model$Cs_)
#   best.c <- as.numeric(logistic_cv_result$model$C_)
#   avg.score <- logistic_cv_result$avg.score[[class.use]]
#   df <- data.frame(log10c = log10(c.vals),
#                    avg.score = avg.score,
#                    best.model = c.vals == best.c)
#   cols <- c("FALSE" = "blue", "TRUE" = "red")
#   p <- ggplot(df, aes(log10c, avg.score, best.model)) +
#     geom_line(aes(x = log10c, y = avg.score), color = 'blue') +
#     geom_point(aes(x = log10c, y = avg.score, color = best.model)) +
#     scale_color_manual(values = cols) +
#     xlab(xlab) + ylab(ylab)
#   return(p)
# }

fit_logistic_train_test <- function(X, y, split_fraction = 0.25, random_state = 12345L, ...) {
  ## INPUTS:
  ##  X = input data matrix
  ##  y = numeric vector of results
  ##  split_fraction = fraction to hold out as test set
  ##  random_state = random seed
  ## RETURNS:
  ##  list of class logistic_train_test_result, which contains the final model fit on full data,
  ##  training accuracy for said model, testing accuracy for said model, and coefficients found on
  ##  the training data
  if(!is.character(colnames(X))) {
    colnames(X) <- 1:ncol(X)
  }
  res <- MODEL_SELECTION$train_test_split(X, y, test_size = split_fraction, stratify = y,
                                          random_state = random_state)
  X_train <- res[[1]]
  X_test <- res[[2]]
  y_train <- res[[3]]
  y_test <- res[[4]]
  rm(res)
  colnames(X_train) <- colnames(X_test) <- colnames(X)
  rm(X)
  gc(full = TRUE)
  print('Fitting Model on Training Data')
  print(Sys.time())
  logistic_cv_result <- fit_logistic(X_train, y_train, random_state = random_state, ...)
  print('Calculating Model Accuracy')
  print(Sys.time())
  y_pred_train <- logistic_cv_result$model$predict(X_train)
  y_pred_test <- logistic_cv_result$model$predict(X_test)
  train_accuracy <- METRICS$accuracy_score(y_true = y_train, y_pred = y_pred_train)
  test_accuracy <- METRICS$accuracy_score(y_true = y_test, y_pred = y_pred_test)
  logistic_train_test_result <- list(final_model = logistic_cv_result,
                                    train_accuracy = train_accuracy,
                                    test_accuracy = test_accuracy,
                                     train_coefs = logistic_cv_result$model$coef_)
  class(logistic_train_test_result) <- 'logistic_train_test_result'
  print('Finished')
  print(Sys.time())
  return(logistic_train_test_result)
}

fit_logistic_multi <- function(X, y, num_splits = 5, split_fraction = 0.25, base_seed = 12345L, ...) {
  ## INPUTS:
  ##  X, y, = see fit_logistic
  ##  num_splits = number of random splits to do
  ##  split_fraction = see fit_logistic_train_test
  ##  base_seed = smallest of a sequence of random seeds to pass to fit_logistic_train_test
  ## RETURNS:
  ##  list with train and test accuracy values for multiple splits, a matrix of model coefficients
  ##  obtained on training data (if a binary classification problem) or a list of such matrices,
  ##  with one for each class, as well as a logistic_cv_result 
  train_accuracy_vals <- c()
  test_accuracy_vals <- c()
  model_list <- list()
  
  n_class <- length(unique(y))
  if (n_class > 2) {
    multiclass <- TRUE
    mat_list <- list()
  } else {
    multiclass <- FALSE
    coefs_mat <- c()
  }
  for (i in 1:num_splits) {
    print(paste('Running for split', i, 'of', num_splits))
    print(Sys.time())
    logistic_train_test_result <- fit_logistic_train_test(X, y, split_fraction = split_fraction,
                                                          random_state = as.integer(base_seed*i), ...)
    train_accuracy_vals <- c(train_accuracy_vals, logistic_train_test_result$train_accuracy)
    test_accuracy_vals <- c(test_accuracy_vals, logistic_train_test_result$test_accuracy)
    if (multiclass) {
      for (j in 1:n_class) {
        if (i == 1) {
          mat_list[[j]] <- logistic_train_test_result$train_coefs[j,]
        } else {
          mat_list[[j]] <- rbind(mat_list[[j]], logistic_train_test_result$train_coefs[j,])
        }
      }
    } else {
      coefs_mat <- rbind(coefs_mat, logistic_train_test_result$train_coefs)
    }
    model_list[[i]] <- logistic_train_test_result$final_model
  }
  if (multiclass) {
    for (j in 1:n_class) {
      rownames(mat_list[[j]]) <- paste0('split_', 1:num_splits)
    }
    coefs_mat <- mat_list
    rm(mat_list)
  } else {
    rownames(coefs_mat) <- paste0('split_', 1:num_splits)
  }
  
  model_select_index <- which(as.numeric(test_accuracy_vals) == max(test_accuracy_vals))
  if (length(model_select_index) > 1) {
    warning('multiple models have test accuracy equal to maximum. selecting first one')
    model_select_index <- model_select_index[1]
  }
  final_model <- model_list[[model_select_index]]
  
  multi_fit_result <- list(train_accuracy_vals = train_accuracy_vals, 
                           test_accuracy_vals = test_accuracy_vals,
                           coefs_mat = coefs_mat,
                           final_model = final_model)
  class(multi_fit_result) <- 'multi_fit_result'
  print('Finished')
  print(Sys.time())
  return(multi_fit_result)
}

fit_knn_cv <- function(X, y, k_vals = 3:20) {
  ## Fit KNN classifier with 5 fold CV
  ## INPUTS:
  ##  X = input data
  ##  y = class labels
  ##  k_vals = range of k values to test
  ## OUTPUTS:
  ##  model fit with best performing value of k, as well as average 
  tryCatch({
    cv_avg <- c()
    cv_sdev <- c()
    for (k in k_vals) {
      knn_model <- SKLEARN_NEIGHBORS$KNeighborsClassifier(n_neighbors = k, metric = 'minkowski', p = 2L)
      cv_scores <- MODEL_SELECTION$cross_val_score(estimator = knn_model, X = X, y = y, cv = 5L, scoring = 'accuracy')
      cv_avg <- c(cv_avg, mean(cv_scores))
      cv_sdev <- c(cv_sdev, sd(cv_scores))
    }
    # take value of k with maximum CV score average. break ties by taking smaller value k
    k_select <- k_vals[which(cv_avg == max(cv_avg))]
    k_select <- k_select[order(k_select, decreasing = FALSE)][1]
    model_output<- SKLEARN_NEIGHBORS$KNeighborsClassifier(n_neighbors = k_select, metric = 'minkowski', p = 2)
    model_output$fit(X = X, y = y)
    names(cv_sdev) <- names(cv_avg) <- k_vals
    output <- list(model = model_output, cv_avg = cv_avg, cv_sdev = cv_sdev, k = k_select)
    class(output) <- 'knn_cv_result'
    return(output)
  }, error = function(e) {
    msg <- paste('failure of knn_cv:', e)
    class_X <- class(reticulate::r_to_py(X))
    class_y <- class(reticulate::r_to_py(y))
    class_k_vals <- class(reticulate::r_to_py(k_vals))
    msg <- paste(msg, '\n', 'Checking classes of X and y, and k_vals')
    msg <- paste(msg, '\n X:', class_X)
    msg <- paste(msg, '\n y:', class_y)
    msg <- paste(msg, '\n k_vals:', class_k_vals)
    stop(msg)
  })
  
}

fit_knn_cv_train_test <- function(X, y, split_fraction = 0.25, k_vals = 3:20, random_state = 42L) {
  if(!is.character(colnames(X))) {
    colnames(X) <- 1:ncol(X)
  }
  res <- MODEL_SELECTION$train_test_split(X, y, test_size = split_fraction, stratify = y,
                                          random_state = random_state)
  X_train <- res[[1]]
  X_test <- res[[2]]
  y_train <- res[[3]]
  y_test <- res[[4]]
  rm(res)
  colnames(X_train) <- colnames(X_test) <- colnames(X)
  rm(X)
  gc(full = TRUE)
  print('Fitting Model on Training Data')
  print(Sys.time())
  knn_cv_result <- fit_knn_cv(X_train, y_train, k_vals = k_vals)
  print('Calculating Model Accuracy')
  print(Sys.time())
  y_pred_train <- knn_cv_result$model$predict(X_train)
  y_pred_test <- knn_cv_result$model$predict(X_test)
  train_accuracy <- METRICS$accuracy_score(y_true = y_train, y_pred = y_pred_train)
  test_accuracy <- METRICS$accuracy_score(y_true = y_test, y_pred = y_pred_test)
  knn_cv_train_test_result <- list(final_model = knn_cv_result,
                                   train_accuracy = train_accuracy,
                                   test_accuracy = test_accuracy)
  class(knn_cv_train_test_result) <- 'knn_cv_train_test_result'
  print('Finished')
  print(Sys.time())
  return(knn_cv_train_test_result)
}

fit_knn_cv_multi <- function(X, y, num_splits = 5, split_fraction = 0.25, base_seed = 12345L, ...) {
  ## INPUTS:
  ##  X, y, = see fit_logistic
  ##  num_splits = number of random splits to do
  ##  split_fraction = see knn_cv_train_test
  ##  base_seed = smallest of a sequence of random seeds to pass to knn_cv_train_test
  ## RETURNS:
  ##  list with train and test accuracy values for multiple splits, a numeric vector of 
  ##  values k selected for each split, as well as a knn_cv_train_test_result for best performing
  ##  model
  train_accuracy_vals <- c()
  test_accuracy_vals <- c()
  model_list <- list()
  
  n_class <- length(unique(y))
  k_selected <- c()
  for (i in 1:num_splits) {
    print(paste('Running for split', i, 'of', num_splits))
    print(Sys.time())
    knn_cv_train_test_result <- fit_knn_cv_train_test(X, y, split_fraction = split_fraction,
                                                      random_state = as.integer(base_seed*i), ...)
    train_accuracy_vals <- c(train_accuracy_vals, knn_cv_train_test_result$train_accuracy)
    test_accuracy_vals <- c(test_accuracy_vals, knn_cv_train_test_result$test_accuracy)
    k_selected <- c(k_selected, knn_cv_train_test_result$final_model$k)
    model_list[[i]] <- knn_cv_train_test_result$final_model
  }
  # select model with k = median value
  model_select_index <- which(k_selected == floor(median(k_selected)))[1]
  final_model <- model_list[[model_select_index]]
  multi_fit_result <- list(train_accuracy_vals = train_accuracy_vals, 
                           test_accuracy_vals = test_accuracy_vals,
                           k_vals = k_selected,
                           final_model = final_model)
  class(multi_fit_result) <- 'multi_fit_knn'
  print('Finished')
  print(Sys.time())
  return(multi_fit_result)
}

# load_model <- function(inp_obj) {
#   # load previously calculated logistic regression model into
#   # object's model slot. meant to deal with fact that
#   # saving R objects with pointers to python objects
#   # via reticulate package will not result in python objects
#   # being saved with the R objects
#   tryCatch({
#     stopifnot(class(inp_obj) %in% c('multi_fit_result', 'logistic_cv_result', 'logistic_train_test_result'))
#   }, error = function(e) {
#     stop(paste('inappropriate object of class', class(inp_obj)))
#   })
#   
#   if (class(inp_obj) %in% c('multi_fit_result', 'logistic_train_test_result')) {
#     if (!class(inp_obj$final_model) == 'logistic_cv_result') {
#       stop('no final model stored')
#     }
#     inp_obj$final_model <- load_model(inp_obj$final_model)
#   } else if (class(inp_obj) == 'logistic_cv_result') {
#     if (!is.character(inp_obj$outfile_path) || !file.exists(inp_obj$outfile_path)) {
#       stop('could not load in model for object. either model was not saved in previous R session\
#     or object\'s filepath does not specify a valid relative path from current running script')
#     }
#     inp_obj$model <- reticulate::py_load_object(filename = inp_obj$outfile_path)
#   }
#   return(inp_obj)
# }
# result <- fit_logistic_multi(X, y)

make_table <- function(df, var1, var2, filter_by = NULL, filter_class = NULL) {
  ## INPUTS:
  ##  df = data.frame
  ##  var1 = categorical variable
  ##  var2 = 2 class categorical variable
  ##  filter_by = filter samples by a certain categorical variable
  ##  filter_class = class(es) acceptable for filter
  ## RETURNS:
  ##  data.frame representing contingency table of classes, 
  ##  + proportions of second class of var2 in each class denoted by var1
  df <- filter_df(df, filter_by, filter_class)
  if (!length(unique(df[,var2])) == 2) {
    stop('expect var2 to denote a categorical variable with 2 classes')
  }
  new_table <- table(df[,c(var1, var2)])
  new_df <- data.frame(new_table[,1], new_table[,2])
  colnames(new_df) <- colnames(new_table)
  proportion <- new_df[,2]/rowSums(new_df)
  new_df$proportion <- proportion
  return(new_df)
}

calcAvgExp <- function(object, genes.use, signature.name = 'signature', assay.use = 1) {
  ## object: Seurat object, or matrix. If one of the first 2, extracts data matrix
  ##          from said object.
  ## genes.use: genes to use
  ## signature.name: name for signature if adding avg expression to metadata. will be added to colnames of metadata
  ## assay.use: if object is SummarizedExperiment, which assay to index.
  
  ## return: either a data matrix or object of same class as input with avg expression for geneset added to sample metadata
  if (class(object) == 'seurat') {
    avg.signature <- calcAvgExp(as.matrix(object@scale.data), genes.use, signature.name)
    object <- Seurat::AddMetaData(object, metadata = as.data.frame(avg.signature), col.name = signature.name)
  } else if (is(object, 'SummarizedExperiment')) {
    avg.signature <- calcAvgExp(assay(object, assay.use), genes.use, signature.name)
    new.metadata <- DataFrame(cbind(as.data.frame(colData(object)), as.data.frame(avg.signature)))
    colnames(new.metadata)[ncol(new.metadata)] <- signature.name
    object@colData <- new.metadata
    
  } else if (class(object) == 'matrix') {
    stopifnot(is.character(genes.use))
    stopifnot(is.character(signature.name))
    genes.drop <- which(!genes.use %in% rownames(object))
    if (length(genes.drop) > 0) {
      warning(paste('using', length(genes.use) - length(genes.drop), 'of', length(genes.use), 'genes for signature', signature.name))
      genes.use <- genes.use[-genes.drop]
    }
    if (length(genes.use) < 1) {
      stop(paste('no genes for signature', signature))
    }
    avg.signature <- apply(matrix(object[genes.use,], nrow = length(genes.use)), MARGIN = 2, FUN = mean)
    return(avg.signature)
  } else {
    stop('unsupported input for object argument')
  }
}
calcAvgChrMat <- function(object, BM.mapping, chr.use = NULL, std.norm = T, min.genes = 50, assay.use = 'scale.data') {
  # calculate average expression of each chromosome.
  # standard normalize if so specified

  # object: either a matrix or a seurat object (from which scaled data matrix retrieved)
  # BM.mapping: data.frame with genes in 1st column chromosomes in second
  # std.norm: avg chromosome expression is standard normalized within sample if set to TRUE
  # min.genes: minimum genes per chromosome
  # chr.use: character vector denoting cromosomes to use
  # return: list containing output matrix (matrix), data frame containing final set of genes
  # used in calculating chromsome averages and the corresponding chromosomes (chr.mapping), and a
  # data frame showing chromsomes, number of genes on each chromosome,  and
  # number of genes used in calculating average (chr.summary)

  if (class(object) == 'seurat') {
    calcAvgChrMat(object = as.matrix(Seurat::GetAssayData(object, assay.type = 'RNA', slot = assay.use)),
                  BM.mapping = BM.mapping,
                  chr.use = chr.use,
                  std.norm = std.norm,
                  min.genes = min.genes)
  } else if (class(object) == 'matrix') {
    stopifnot(is.data.frame(BM.mapping))
    duplicate.genes <- duplicated(BM.mapping[,1])
    if (any(duplicate.genes)) {
      warning(paste('removing', sum(duplicate.genes), 'duplicate genes in provided input for BM.mapping'))
      BM.mapping <- BM.mapping[-which(duplicate.genes)]
    }
    if (is.null(chr.use)) {
      chr.use <- unique(BM.mapping[,2])
    }
    BM.mapping <- BM.mapping[BM.mapping[,2] %in% chr.use,]
    genes.use <- intersect(BM.mapping[,1], rownames(object))
    if (!all(chr.use %in% BM.mapping[,2])) {
      stop(paste('Not all chromosomes in chr.use represented in BM.mapping'))
    }
    object <- object[genes.use,]
    # number of genes used per chromsome will be updated as these genes are found
    chr.summary <- data.frame(chr = chr.use, n.genes = sapply(chr.use, FUN = function(x) {
      length(which(BM.mapping[,2] == x))
    }), n.genes.used = numeric(length(chr.use)))

    output.mat <- c()
    chr.mapping <- c()
    skipped.chr <- c()
    for (i in 1:length(chr.use)) {
      chr <- chr.use[i]
      chr.genes <- as.character(BM.mapping[which(BM.mapping[,2] == chr),1])
      chr.genes <- chr.genes[which(chr.genes %in% genes.use)]
      if (length(chr.genes) < min.genes) {
        warning(paste0(paste('Fewer genes than min.genes for region', chr, '. min.genes = ',
                             min.genes, '. genes on chromsome in data matrix:', length(chr.genes),'.\nskipping\n')))
        skipped.chr <- c(skipped.chr, i)
        next
      }
      chr.summary[which(chr.summary$chr == chr), 'n.genes.used'] <- length(chr.genes)
      chr.mapping <- rbind(chr.mapping, cbind(chr.genes, rep(chr, length(chr.genes))))
      output.mat <- rbind(output.mat,
                          calcAvgExp(object = object,
                                     genes.use = chr.genes,
                                     signature.name = as.character(chr))
      )
    }
    if (length(skipped.chr) > 0) {
      chr.use <- chr.use[-skipped.chr]
    }
    rownames(output.mat) <- chr.use
    colnames(output.mat) <- colnames(object)
    if (std.norm) {
      output.mat <- scale((output.mat), center = T, scale = T)
    }

    colnames(chr.mapping) <- c('chr', 'gene')
    chr.mapping <- as.data.frame(chr.mapping)

    output.list <- list(output.mat = output.mat, chr.mapping = chr.mapping, chr.summary = chr.summary)
    return(output.list)
  } else {
    stop()
  }
}

identify.glioma <- function(seurat_obj, clusters.use = 'RNA') {
  # seurat_obj = Seurat object, with z-scored average chromosome expression
  #           in colData. Expect columns with names chr.x, for x in 1:22. Expect a column
  #           'cluster' in colData denoting RNA based clusters
  # clusters.use = character. If 'RNA', classify individual clusters as glioma, other,
  #                 or unclassified based on distribution of probabilities for
  #                 glioma and non glioma clusters, then take set of all clusters classified
  #                 as glioma and label the cells as glioma. If 'chromosome', assign 
  #                 as glioma or other based on chromosome based cluster found with max
  #                 probability
  # details: Note that chromosome based clusters are found based on a Gaussian Mixture
  #           Model fit in space of chromosomes 7 and 10. The glioma cluster is determined
  #           using a crude classifier score with the chromosome 7 and 10 z scores z-normalized
  #           by sample (so that features are weighted equally) and calculated as chr.7 - chr.10.
  #           glioma.cluster is determined as the cluster with thte higher distribution of this
  #           score by means of a wilcoxon rank sum test
  # returns: Seurat object with probabilities of coming from glioma cluster or
  #           other cluster ('glioma.score', 'other.score'), call for glioma identify ('is.glioma'),
  #           classifier.score, all added to colData
  #
  
  # cluster based off of z-scored avg chromosome expression for chromosome 7, chromosome 10
  meta_data <- seurat_obj@meta.data
  chr_cols <- c('chr.7', 'chr.10')
  tryCatch({
    stopifnot(all(chr_cols %in% colnames(meta_data)))
  }, error = function(e) {
    stop(paste('must have columns', paste(chr_cols, collapse = ','), 'in data'))
  })
  
  chr.mat <- as.matrix(meta_data[,chr_cols])
  chr.GMM <- Mclust(data = chr.mat,
                    G = 2)
  # chr.GMM <- Mclust(data = chr.mat,
  #                   G = 2)
  chr.clust.labs <- apply(chr.GMM$z, MARGIN = 1, FUN = function(x) {which(x == max(x))})
  
  # chr.clust.labs <- cutree(hclust(dist(chr.mat)), k = 2)
  meta_data$chr.clust <- chr.clust.labs
  
  # make a crude 'classifier score' to score cells for 'glioma ness'
  # scale chr 7 score, and chr 10 score to unit variance
  # glioma classifier score is glioma stem cell score + chr 7 score - chr 10 score
  # recall that chr7 and chr10 score are within cell z-scores of average chromosome expression
  
  glioma.features <- chr.mat[,c('chr.7', 'chr.10')]
  glioma.features[,'chr.10'] <- -1*glioma.features[,'chr.10']
  glioma.features <- scale(glioma.features, center = T, scale = T)
  classifier.scores <- rowSums(glioma.features)
  meta_data$classifier.score <- classifier.scores
  p.vals <- numeric(2)
  p.val.thresh <- 0.05
  for (i in 1:2) {
    chr.clust.memb <- which(meta_data$chr.clust == i)
    wilcox.result <- wilcox.test(meta_data$classifier.score[chr.clust.memb],
                                 meta_data$classifier.score[-chr.clust.memb],
                                 alternative = 'greater')
    p.vals[i] <- wilcox.result$p.value
  }
  
  if (!any(p.vals < p.val.thresh)) {
    stop('could not confidently identify chromosome based cluster with lowest p-value')
  }
  cancer.clust <- which(p.vals == min(p.vals))
  
  meta_data$glioma.score <- chr.GMM$z[,cancer.clust]
  meta_data$other.score <- chr.GMM$z[,-cancer.clust]
  
  
  # divide RNA based clusters into 'glioma' clusters and 'non glioma' clusters by
  # glioma score
  meta_data$is.glioma <- character(nrow(meta_data))
  if (clusters.use == 'chromosome') {
    meta_data$is.glioma <- ifelse(meta_data$chr.clust  == cancer.clust, 'glioma', 'other')
  } else if (clusters.use == 'RNA') {
    for (i in unique(meta_data$cluster)) {
      cell.inds <- meta_data$cluster == i
      glioma.scores.i <- meta_data[cell.inds,'glioma.score']
      other.scores.i <- meta_data[cell.inds,'other.score']
      p.val.glioma <- wilcox.test(glioma.scores.i, other.scores.i, alternative = 'greater')$p.value
      p.val.other <- wilcox.test(glioma.scores.i, other.scores.i, alternative = 'less')$p.value
      if (p.val.glioma < p.val.thresh) {
        is.glioma.i <- 'glioma'
      } else if (p.val.other < p.val.thresh) {
        is.glioma.i <- 'other'
      } else {
        is.glioma.i <- 'unclassified'
      }
      meta_data[cell.inds,'is.glioma'] <- is.glioma.i
    } 
  } else {
    stop('clusters.use must be one of RNA, chromosome')
  }
  meta_data$is.glioma <- factor(meta_data$is.glioma)
  seurat_obj <- Seurat::AddMetaData(seurat_obj, metadata = meta_data)
  return(seurat_obj)
}