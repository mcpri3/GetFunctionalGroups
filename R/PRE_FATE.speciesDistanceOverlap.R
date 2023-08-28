### HEADER #####################################################################
##' @title Computation of niche overlap distances between species 
##' 
##' @name PRE_FATE.speciesDistanceOverlap
##'
##' @author Maya Guéguen, modified by Marie-Caroline Prima
##' 
##' @description This script is designed to create a distance matrix between 
##' species, based on co-occurrence of species. 
##'              
##' 
##' @param mat.overlap.object two options, depending on the value of 
##' \code{mat.overlap.option} :
##' \itemize{
##'   \item (\code{PCA} option) a \code{list} with 2 elements :
##'   \describe{
##'     \item{\code{tab.dom.PA}}{a \code{matrix} or \code{data.frame} with 
##'     sites in rows and species in columns, containing either \code{NA}, 
##'     \code{0} or \code{1} (see \code{\link{PRE_FATE.selectDominant}})}
##'     \item{\code{tab.env}}{a \code{matrix} or \code{data.frame} with 
##'     sites in rows and environmental variables in columns}
##'   }
##' }
##' 
##' @details 
##' 
##' This function allows to obtain a \strong{distance matrix between species} 
##' (1 - Schoeners D), based on niche overlap information :
##' 
##' \itemize{
##'   \item The degree of niche overlap will 
##'   be computed using the \code{\link[ecospat]{ecospat.niche.overlap}}. 
##' }
##' 
##' 
##' @return A \code{matrix} containing overlap distances between each pair 
##' of species, calculated as \code{1 - Schoeners D}.
##'  
##' @keywords niche overlap, Schoeners D
##' 
##' @seealso \code{\link[ecospat]{ecospat.niche.overlap}} 
##' 
##' @examples
##' 
##' @export
##' 
##' @importFrom stats as.dist
##' @importFrom methods as
##' @importFrom ade4 dudi.hillsmith suprow
##' @importFrom raster raster extension
##' @importFrom phyloclim niche.overlap
##' 
## END OF HEADER ###############################################################


PRE_FATE.speciesDistanceOverlap = function(mat.overlap.object, nme.pca
){
  library(utils)
  library(foreach)
  library(ecospat)
  #############################################################################
  
  cat("\n\n #------------------------------------------------------------#")
  cat("\n # PRE_FATE.speciesDistanceOverlap")
  cat("\n #------------------------------------------------------------# \n\n")
  
  
      mat.overlap = mat.overlap.object
      tab.dom.PA = mat.overlap[[1]]
      tab.env = mat.overlap[[2]]
      ind_rownames = sort(unique(intersect(rownames(tab.dom.PA), rownames(tab.env))))
      tab.dom.PA = tab.dom.PA[ind_rownames, ]
      tab.env = tab.env[ind_rownames, ]
        
        ## Calculate PCA for all environment
        pca.env = dudi.hillsmith(tab.env, scannf = F, nf = 2)
        saveRDS(pca.env, here::here(paste0('data/derived-data/ePCA/ePCA_', nme.pca)))
        
        scores.env = pca.env$li
 
        ## Calculate overlap matrix
        PROGRESS = txtProgressBar(min = 0, max = ncol(tab.dom.PA), style = 3)
        grid.list = foreach(ii = 1:ncol(tab.dom.PA)) %do%
          {
            setTxtProgressBar(pb = PROGRESS, value = ii)
            si.01 = rownames(tab.dom.PA)[which(!is.na(tab.dom.PA[, ii]))]
            si.1 = rownames(tab.dom.PA)[which(tab.dom.PA[, ii] > 0)]
            if (length(si.1) > 5)
            {
              ind.01 = which(rownames(tab.env) %in% si.01)
              ind.1 = which(rownames(tab.env) %in% si.1)
              scores.sp1.01 = suprow(pca.env, tab.env[ind.01, ])$li
              scores.sp1.1 = suprow(pca.env, tab.env[ind.1, ])$li

              openxlsx::write.xlsx(scores.sp1.1, here::here(paste0('data/derived-data/ePCA/SpeciesCoords_',colnames(tab.dom.PA)[ii], 
                                                                   '_ePCA_',nme.pca, '.xlsx')))
              grid.clim.sp1 = ecospat.grid.clim.dyn(glob = scores.env
                                                    , glob1 = scores.sp1.01
                                                    , sp = scores.sp1.1
                                                    , R = 100, th.sp = 0)
              return(grid.clim.sp1)
            } else { return(NULL) }
          }
        close(PROGRESS)
        
        n.sel = ncol(tab.dom.PA)
        mat.overlap = matrix(NA, nrow = n.sel, ncol = n.sel
                             , dimnames = list(colnames(tab.dom.PA), colnames(tab.dom.PA)))
        PROGRESS = txtProgressBar(min = 0, max = n.sel, style = 3)
        for (ii in 1:(n.sel-1))
        {
          setTxtProgressBar(pb = PROGRESS, value = ii)          
          if (!is.null(grid.list[[ii]]))
          {
            for(jj in (ii+1):n.sel)
            {
              if (!is.null(grid.list[[jj]]))
              {
                res = ecospat.niche.overlap(grid.list[[ii]], grid.list[[jj]], cor = TRUE)$D
                mat.overlap[ii, jj] = res
              }
            }
          }
        }
        close(PROGRESS)
        
        mat.overlap[lower.tri(mat.overlap, diag = FALSE)] = t(mat.overlap)[lower.tri(mat.overlap, diag = FALSE)]
        diag(mat.overlap) = 1
  
  
  ## Remove species with no overlap
  no_NA_values = apply(mat.overlap, 2, function(x) sum(is.na(x)))
  ind_NA_values = which(no_NA_values >= nrow(mat.overlap) - 1)
  if (length(ind_NA_values) > 0)
  {
    warning(paste0("Missing data!\n `mat.overlap` contains some species with no overlap values : "
                   , paste0(colnames(mat.overlap)[ind_NA_values], collapse = ", ")
                   , "\nThese species will not be taken into account ! \n\n"
    ))
    mat.overlap = mat.overlap[-ind_NA_values, -ind_NA_values]
    # names_species.overlap = sort(unique(colnames(mat.overlap)))
    # if (nrow(mat.overlap) <= 1)
    # {
    #   stop("Wrong dimension(s) of data!\n `mat.overlap` does not have the appropriate number of rows (>=2)")
    # }
  }
  
  ## Transform into dissimilarity distances (instead of similarity)
  mat.OVERLAP = (1 - mat.overlap)
  
  cat("\n> Done!\n")
  
  return(mat.OVERLAP)
  
}