#' A Function That Computes Transcriptogram of a Gene.
#'
#' This function computes  RNA-Seq signal matrix of a gene from coverage signal.
#' @param exons_bed a bed-like dataframe of exons.
#' @param coverage_bwfiles a vector of bigwig files.
#' @param bin_size an integer specifying the bin to discretize RNA-Seq signal.
#' @param tmp_dir a string specifying temporary directory.
#' @export
#' @importFrom utils write.table
#' @importFrom utils read.table
get_transcriptogram = function(exons_bed, coverage_bwfiles, bin_size=5, tmp_dir="tmp") {
  bar = lapply(1:nrow(exons_bed), function(i) {
    before_region_start_length = 0
    after_region_start_length = exons_bed[i,3] - exons_bed[i,2]
    nb_bin = ceiling((before_region_start_length + after_region_start_length) / bin_size)
    after_region_start_length = nb_bin * bin_size
    dir.create(tmp_dir, recursive=TRUE, showWarnings=FALSE)
    matrix_out_filename = paste0(tmp_dir, "/matrix_", i, "_RNA.txt.gz")
    if (!file.exists(matrix_out_filename)) {
      regions_filename = paste0("tmp/exons_bed_", i, ".bed")
      write.table(exons_bed[i,1:6], file=regions_filename, sep="\t", quote=FALSE,row.names=FALSE, col.names=FALSE)      
      out_filename = matrix_out_filename
      matrix_out_filename = dt_compute_matrix(regions_filename=regions_filename, score_filename=coverage_bwfiles, out_filename=out_filename, bin_size=bin_size, before_region_start_length=before_region_start_length, after_region_start_length=after_region_start_length, FORCE_EXEC=TRUE)
      # matrix_file = matrix_out_filename
      # out_filename = "hm.png"
      # dt_plot_heatmap(matrix_file, out_filename)    
    }

    foo = read.table(gzfile(matrix_out_filename), skip=1, stringsAsFactors=FALSE)
    # head(foo)
    # print(dim(foo))
    data = matrix(as.numeric(unlist(t(foo[,-(1:6)]))), nb_bin*nrow(foo))
    rownames(data) = paste0("exon_", rep(i, nb_bin), "_", 1:nb_bin)
    colnames(data) = do.call(rbind, strsplit(do.call(rbind,lapply(strsplit(coverage_bwfiles, "/"), rev))[,1], "_"))[,1]
  
    return(data)
  })
  bar = do.call(rbind,bar)
  return(bar)  
}



#' Binding of deepTools `plotHeatmap` Function.
#'
#' This function binds deepTools `plotHeatmap` function.
#' @param matrix_file           deepTools `plotHeatmap` parameter.  
#' @param out_filename        deepTools `plotHeatmap` parameter.  
#' @param y_max               deepTools `plotHeatmap` parameter.  
#' @param y_min               deepTools `plotHeatmap` parameter.  
#' @param sort_regions        deepTools `sortRegions` parameter.  
#' @param FORCE_EXEC a boolean that force call to deepTools `computeMatrix` function.
#' @export
dt_plot_heatmap = function (matrix_file, out_filename, y_max, y_min, sort_regions="no", FORCE_EXEC=FALSE) {
  command = "plotHeatmap"  
  args = paste(
    c(
      "--matrixFile"                , matrix_file, 
      "--outFileName"             , out_filename,     
      "--colorMap"                , "bwr", 
      "--sortRegions"             , sort_regions
    ), 
  collapse = " ")
  if (!missing(y_max)) {
    args = paste(c(args, "--yMax", y_max), collapse=" ")
  }
  if (!missing(y_min)) {
    args = paste(c(args, "--yMin", y_min), collapse=" ")
  }
  print(paste(command, args))
  if (substr(Sys.info()["nodename"],1,4)=="luke" | FORCE_EXEC) {
    system2(command=command, args=args)
  } else {
    print(paste("Skiping deep tools calls on ", Sys.info()["nodename"]))
  }
  return(out_filename)
}

#' Binding of deepTools `computeMatrix` Function.
#'
#' This function binds deepTools `computeMatrix` function.
#' @param regions_filename           deepTools `computeMatrix` parameter.  
#' @param score_filename             deepTools `computeMatrix` parameter.  
#' @param out_filename               deepTools `computeMatrix` parameter.  
#' @param bin_size                   deepTools `computeMatrix` parameter.  
#' @param before_region_start_length deepTools `computeMatrix` parameter.  
#' @param after_region_start_length  deepTools `computeMatrix` parameter.  
#' @param number_of_processors       deepTools `computeMatrix` parameter.  
#' @param blacklist_filename         deepTools `computeMatrix` parameter.  
#' @param FORCE_EXEC a boolean that force call to deepTools `computeMatrix` function.
#' @export
dt_compute_matrix = function (regions_filename, score_filename, out_filename, bin_size=50, before_region_start_length=5000, after_region_start_length=5000, number_of_processors=12, blacklist_filename=NULL, FORCE_EXEC=FALSE) {
  command = "computeMatrix"  
  args = paste(
    c("reference-point", 
      "-R"                        , regions_filename, 
      "-S"                        , score_filename, 
      "--outFileName"             , out_filename,     
      "--referencePoint"          , "TSS", 
      "--binSize"                 , bin_size,  
      "--beforeRegionStartLength" , before_region_start_length,
      "--afterRegionStartLength"  , after_region_start_length,
      "--numberOfProcessors"      , number_of_processors,  
      "--sortRegions"             , "keep" 

    ), 
  collapse = " ")
  if (!is.null(blacklist_filename)) {
    args = paste(args, "--blackListFileName" , blacklist_filename) 
  }
  print(paste(command, args))
  if (substr(Sys.info()["nodename"],1,4)=="luke" | FORCE_EXEC) {
    system2(command=command, args=args)
  } else {
    print(paste("Skiping deep tools calls on ", Sys.info()["nodename"]))
  }
  return(out_filename)
}


#' Binding of deepTools `computeMatrix` Function.
#'
#' This function binds deepTools `computeMatrix` function.
#' @param regions_filename           deepTools `computeMatrix` parameter.  
#' @param score_filename             deepTools `computeMatrix` parameter.  
#' @param out_filename               deepTools `computeMatrix` parameter.  
#' @param bin_size                deepTools `computeMatrix` parameter.  
#' @param before_region_start_length deepTools `computeMatrix` parameter.  
#' @param after_region_start_length  deepTools `computeMatrix` parameter.  
#' @param region_body_length         deepTools `computeMatrix` parameter.  
#' @param number_of_processors       deepTools `computeMatrix` parameter.  
#' @param blacklist_filename         deepTools `computeMatrix` parameter.  
#' @param FORCE_EXEC a boolean that force call to deepTools `computeMatrix` function.
#' @export
dt_compute_matrix_scale = function (regions_filename, score_filename, out_filename, bin_size=50, before_region_start_length=0, after_region_start_length=0, region_body_length, number_of_processors=12, blacklist_filename=NULL, FORCE_EXEC=FALSE) {
  command = "computeMatrix"  
  args = paste(
    c("scale-regions", 
      "-R"                        , regions_filename, 
      "-S"                        , score_filename, 
      "--outFileName"             , out_filename,     
      "--startLabel"              , "win_start", 
      "--endLabel"                , "win_end", 
      "--regionBodyLength"        , region_body_length,  
      "--binSize"                 , bin_size,  
      "--beforeRegionStartLength" , before_region_start_length,
      "--afterRegionStartLength"  , after_region_start_length,
      "--numberOfProcessors"      , number_of_processors,  
      "--sortRegions"             , "keep" 

    ), 
  collapse = " ")
  if (!is.null(blacklist_filename)) {
    args = paste(args, "--blackListFileName" , blacklist_filename) 
  }
  print(paste(command, args))
  if (substr(Sys.info()["nodename"],1,4)=="luke" | FORCE_EXEC) {
    system2(command=command, args=args)
  } else {
    print(paste("Skiping deep tools calls on ", Sys.info()["nodename"]))
  }
  return(out_filename)
}


