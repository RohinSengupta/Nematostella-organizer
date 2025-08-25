# This is a function for exporting seurat object as h5ad (anndata format)
# to be analysed by SAMap. Use source('seurat_to_anndata.R') to load the
# function and then use the function:
# seurat_to_anndata(your_seurat_object,'exported file prefix')

seurat_to_anndata<-function(seurat_object,
                            h5ad_name = "exported_seurat"){
  
  require(Seurat)
  require(zellkonverter)
  require(Rhdf5lib)
  require(HDF5Array)
  require(scater)
  
  stopifnot(class(seurat_object)[1]=="Seurat")
  
  message("cell types found:")
  print(levels(Idents(seurat_object)))
  
  message("saving as cell_type variable")
  seurat_object$cell_type<-as.character(Idents(seurat_object))
  
  message("removing excess data from object")
  
  seurat_object <- DietSeurat(
    seurat_object,
    layers = 'counts',
    dimreducs = 'umap',
    graphs = NULL
  )
  
  sce_object<-as.SingleCellExperiment(seurat_object)
  
  assay(sce_object,'logcounts')<-NULL
  assay(sce_object,'scaledata')<-NULL
  
  message("converting to anndata format")
  
  rowData(sce_object)<-DataFrame("gene_id"=rownames(sce_object),'id_type'='gene_id')
  
  writeH5AD(
    sce_object,
    file = paste0(h5ad_name,'.h5ad'),
    colData = TRUE,
    metadata = FALSE,
    version = "0.7.6"
  )
  
  if(file.exists(paste0(h5ad_name,'.h5ad'))){
    message(paste0("object exported as: ",h5ad_name,'.h5ad'))
  }
}
