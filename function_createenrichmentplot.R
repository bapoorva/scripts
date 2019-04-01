library(DOSE)
library(stats)
library(qvalue)
library(fgsea)
library(enrichplot)
creategseaobj = function(geneList, geneSets, exponent=1,
                         nPerm=1000,
                         minGSSize = 10,
                         maxGSSize = 500,
                         pvalueCutoff=0.05,
                         pAdjustMethod="BH"){
  
  tmp_res <- fgsea::fgsea(pathways=geneSets,
                   stats=geneList,
                   nperm=nPerm,
                   minSize=minGSSize,
                   maxSize=maxGSSize,
                   gseaParam=exponent,
                   nproc = 0)
  
  p.adj <- p.adjust(tmp_res$pval, method=pAdjustMethod)
  qvalues <- calculate_qvalue(tmp_res$pval)
  
  params <- list(pvalueCutoff = pvalueCutoff,
                 nPerm = nPerm,
                 pAdjustMethod = pAdjustMethod,
                 exponent = exponent,
                 minGSSize = minGSSize,
                 maxGSSize = maxGSSize
  )
  res <- res[!is.na(res$pvalue),]
  res <- res[ res$pvalue <= pvalueCutoff, ]
  res <- res[ res$p.adjust <= pvalueCutoff, ]
  idx <- order(res$pvalue, decreasing = FALSE)
  res <- res[idx, ]
  
  row.names(res) <- res$ID
  observed_info <- lapply(geneSets[res$ID], function(gs)
    gseaScores(geneSet=gs,
               geneList=geneList,
               exponent=exponent)
  )
  
  if (verbose)
    message("leading edge analysis...")
  
  ledge <- leading_edge(observed_info)
  
  res$rank <- ledge$rank
  res$leading_edge <- ledge$leading_edge
  res$core_enrichment <- sapply(ledge$core_enrichment, paste0, collapse='/')
  
  if (verbose)
    message("done...")
  
  new_res=new("gseaResult",
      result     = res,
      geneSets   = geneSets,
      geneList   = geneList,
      params     = params,
      readable   = FALSE
  )
  
  return(new_res)
}