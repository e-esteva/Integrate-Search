IntegrateSearch=function(ref,query,project_name=NULL,refcellID=NULL,out_dir=getwd(),plot.convergence=F,k=15,s=NULL){
  suppressPackageStartupMessages({
    library(harmony)
    library(Seurat)
    library(umap)
    library(ggplot2)
    library(dplyr)
    library(glue)
  })
  
  if(is.null(s)){
    s = k
  }
  source('/gpfs/data/reizislab/ee699/SingleCellAnalysisTools/IntegrateSearch/SeuratHarmony.R')  
  findNNinRef=function(knn,queryRange,cellID,search_depth){
    df.subset=knn[queryRange,seq(search_depth)]
    
    ref.hits=lapply(seq(dim(df.subset)[1]),function(x) df.subset[x,][df.subset[x,] >  max(queryRange)])
    
    nnREF=suppressWarnings(sapply(seq(dim(df.subset)[1]),function(x) ifelse(length(ref.hits[[x]]) > 0,ref.hits[[x]][1],NA)))
    
    
    celltype.assignments=sapply(seq(length(nnREF)),function(x) ifelse(!is.na(nnREF[x]),cellID[nnREF[x]],'noHit'))
    return(celltype.assignments)
  }
  
  if(class(ref) == 'Seurat' && class(query) =='Seurat'){
    
    features=SelectIntegrationFeatures(list(ref,query))
    
    if(is.null(refcellID)){
      cellID=Idents(ref)
    }else{
      cellID=ref[[refcellID]][[1]]
    }
    
    ref=ref[features]
    query=query[features]
    
    rm(features)
    
    integrated_ = SeuratHarmonyIntegration(object_list=list(query,ref),vars=c('query','ref'))
    
    rm(ref)    
    
    Idents(integrated_) = integrated_$dataset
    message('Running UMAP Seurat Native (UWOT) for visual inspection of batch correction')
    integrated_ = RunUMAP(integrated_,reduction='harmony',dims=1:20,return.model=T)
    p=DimPlot(integrated_)
    ggsave(glue('{out_dir}/batch-corrected-UMAP.pdf'),p,device='pdf')
    message('Running UMAP using umap-learn for KNN graph')
    integrated_knn <- integrated_@reductions$harmony@cell.embeddings %>% umap(n_neighbors=k,verbose=T)
    for(i in s){
     
     
    
     indx=integrated_knn$knn$indexes
     layout=data.frame(integrated_knn$layout)
     query_range=grep('query',integrated_$dataset)
    
     cellID.mod=c(rep('query',length(query_range)),cellID)
     celltypes_q = findNNinRef(knn=indx,queryRange=query_range,cellID=cellID.mod,search_depth=i)
     labels=c(paste(celltypes_q,'query',sep='-'),paste(cellID,'ref',sep='-'))
     layout$celltype=labels
    
     p=ggplot(layout,aes(X1,X2,col=celltype))+geom_point()
     out_dir=ifelse(is.null(out_dir),getwd(),out_dir)
     ggsave(glue('{out_dir}/integrated-UMAP-k={k}_s={i}.pdf'),p,device='pdf',width=20,height=20)

     colName=glue('predicted_cell_id_k={k}_s={i}')

     query[[colName]]=celltypes_q
    }
    project = ifelse(is.null(project_name),'integrate_search',project_name)
    integrated_$celltype=c(celltypes_q,cellID)    
    saveRDS(integrated_knn,glue('{out_dir}/{project}-umap.rds'))

    return(integrated_)    
    
  }
}
