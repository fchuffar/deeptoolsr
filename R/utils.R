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


