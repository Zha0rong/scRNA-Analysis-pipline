$HOSTNAME = ""
params.outdir = 'results'  


if (!params.H5FileInput){params.H5FileInput = ""} 

Channel.fromPath(params.H5FileInput, type: 'any').map{ file -> tuple(file.baseName, file) }.set{g_1_h5_file0_g_0}


process Load_Data_h5 {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /${name}.output.rds$/) "RDSFileOutPut/$filename"}
input:
 set val(name), file(Input) from g_1_h5_file0_g_0

output:
 set val(name),file("${name}.output.rds")  into g_0_rdsFile00_g_3


script:"""

#!/usr/bin/env Rscript
library('Seurat')


Data=Read10X_h5('${Input}')
Data=CreateSeuratObject(Data)
Data[['sample']]='${name}'
saveRDS(Data,'${name}.output.rds')


"""


}

if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 400
    $CPU  = 1
    $MEMORY = 100
    $QUEUE = "long"
}
process Filter_Seurat_Object_Markdown {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /${name}.filtered_seurat.rds$/) "RDSOutput/$filename"}
publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /.*.filtering_report.html$/) "HTML_Output/$filename"}
input:
 set val(name), file(seurat_obj) from g_0_rdsFile00_g_3

output:
 set val(name), file("${name}.filtered_seurat.rds")  into g_3_rdsFile00_g_6
 file "*.filtering_report.html"  into g_3_outputFileHTML11

label 'scrna_seurat'

shell:
minUMI = params.Filter_Seurat_Object_Markdown.minUMI
maxUMI = params.Filter_Seurat_Object_Markdown.maxUMI
nFeature_RNA = params.Filter_Seurat_Object_Markdown.nFeature_RNA
mitoRatio = params.Filter_Seurat_Object_Markdown.mitoRatio
'''
#!/usr/bin/env perl

my $script = <<'EOF';

---
title: "scRNAseq Filtering Reports"
output: 
  html_document:
    theme: united
    highlight: tango
    toc: true
    number_sections: true
    code_folding: hide
    toc_float:
      collapsed: false
      smooth_scroll: false
---

```{r setup, include=TRUE}
knitr::opts_chunk$set(echo = FALSE)
```

```{r error=FALSE, message=FALSE, warning=FALSE, cache=FALSE, results = FALSE}
suppressPackageStartupMessages({
library(dplyr)
library(Seurat)
library(ggplot2)
if(!require(remotes)) install.packages("remotes",repos = "http://cran.us.r-project.org")
remotes::install_github("chris-mcginnis-ucsf/DoubletFinder")
library(DoubletFinder)

})
```

```{r load, message=FALSE, warning=FALSE,include=TRUE}
# Load data from saved seurat objects
obj <- readRDS("!{seurat_obj}")
system(paste0("rm -rf ","!{seurat_obj}"))

if (any(grepl("^MT-",rownames(obj)))|any(grepl("^mt-",rownames(obj)))) {
	if (any(grepl("^MT-",rownames(obj)))) {
	obj[["percent.mt"]]=PercentageFeatureSet(obj,pattern="^MT-")
}
	else if (any(grepl("^mt-",rownames(obj)))) {
	obj[["percent.mt"]]=PercentageFeatureSet(obj,pattern="^mt-")
}} else {
	obj[["percent.mt"]]=0
}

```

## Filtering Cells
We removed cells with:

- Fewer genes than selected cutoff (nFeature_RNA < !{nFeature_RNA})
- Fewer than minimum number of unique transcript molecules (minTranscripts > !{minUMI}).
- More than maximum number of unique transcript molecules (maxTranscripts < !{maxUMI}).
- More than mitocondrial ratio (mitoRatio < %!{mitoRatio})

This will be highly sample dependent. If you find that you lose too many cells filtering with these thresholds, we recommend filtering with lower thresholding (50 genes/100 transcripts) and then removing low quality clusters.

```{r precellcount, warning=FALSE, message=FALSE,fig.width=10,fig.height=10,include=TRUE}

## determine how many cells are there before filtering
beforeFiltering <- ncol(obj)
VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)

plot1 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2 + NoLegend()

print("Summary Statistics for number of genes per cell")

print(summary(obj$nFeature_RNA))

print("Summary Statistics for number of UMIs per cell")

print(summary(obj$nCount_RNA))

print("Summary Statistics for mitochondrial contents per cell")

print(summary(obj$percent.mt))


```

We use DoubletFinder to remove potential multiplets (more than 1 cell in droplets).

```{r doublet removal, warning=FALSE, message=FALSE,results=FALSE,include=TRUE}

obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
obj <- ScaleData(obj)
obj <- RunPCA(obj)

annotations <- obj@meta.data$orig.ident

homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
nExp_poi <- round(0.05*nrow(obj@meta.data))  ## Assuming 5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))



obj <- doubletFinder_v3(obj, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
doublet.classification.name=colnames(obj@meta.data)[ncol(obj@meta.data)]
obj$Doublet.Classification=obj@meta.data[,doublet.classification.name]

VlnPlot(obj,"nCount_RNA",group.by = "Doublet.Classification")
VlnPlot(obj,"nFeature_RNA",group.by = "Doublet.Classification")


```



```{r filtering, warning=FALSE, message=FALSE,include=TRUE}
## filter features

obj <- subset(obj, subset= (nFeature_RNA > !{nFeature_RNA}) &
                          (nCount_RNA > !{minUMI}) & 
                          (nCount_RNA < !{maxUMI}) & 
                          (percent.mt < !{mitoRatio})&
                          (Doublet.Classification=="Singlet"))

```

```{r postfilteringcellcount, warning=FALSE, message=FALSE,include=TRUE}

## determine how many cells are recovered after filtering
afterFiltering <- ncol(obj)

data <- data.frame(matrix(c("BeforeFiltering", "AfterFiltering", beforeFiltering, afterFiltering), nrow=2, ncol=2, dimnames = list(c("Before", "After"), c("Filtering", "Value") )))
data$Filtering=factor(data$Filtering,levels = c("BeforeFiltering", "AfterFiltering"))
data$Value=as.numeric(data$Value)
ggplot(data=data, aes(x=Filtering, y=Value)) +
     geom_bar(stat="identity", fill="steelblue") +
     labs( x="", y = "The # of Cells")+
     geom_text(aes(label=Value), vjust=1.6, color="white", size=3.5)+
     theme_minimal()

# save Seurat object
saveRDS(obj, file="!{name}.filtered_seurat.rds")

```

EOF

open OUT, ">filtering_rmark.rmd";
print OUT $script;
close OUT;

runCommand("Rscript -e 'rmarkdown::render(\\"filtering_rmark.rmd\\",\\"html_document\\", output_file = \\"!{name}.filtering_report.html\\")'   ");

sub runCommand {
            my ($com) = @_;
            my $error = system($com);
            if   ($error) { die "Command failed: $error $com\\n"; }
            else          { print "Command successful: $com\\n"; }
          }

'''



}

//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 20
}
//* platform
//* platform
//* autofill

process Merge_Seurat_Objects {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /merged_filtered_seurat.rds$/) "Merged_Object/$filename"}
input:
 file seurat_obj from g_3_rdsFile00_g_6.collect()

output:
 file "merged_filtered_seurat.rds"  into g_6_rdsFile00_g_9

label 'scrna_seurat'

shell:

'''
#!/usr/bin/env Rscript

list_of_samples <- list.files(pattern = "*.filtered_seurat.rds")
if (length(list_of_samples)==1) {
	obj = readRDS(list_of_samples[1])
} else {
list_of_seurat <- list()
for(i in 1:length(list_of_samples)){
  # print name
  print(list_of_samples[i])
  list_of_seurat[[i]] <- readRDS(list_of_samples[i])
}

obj <- merge(list_of_seurat[[1]],list_of_seurat[-1])
}
saveRDS(obj, file="merged_filtered_seurat.rds")

'''


}


process Normalization_and_Feature_Selection {

publishDir params.outdir, mode: 'copy', saveAs: {filename -> if (filename =~ /Normalized.rds$/) "Normalized_Data/$filename"}
input:
 file seurat_object from g_6_rdsFile00_g_9

output:
 file "Normalized.rds"  into g_9_rdsFile00

label 'scrna_seurat'

shell:
varFeatures = params.Normalization_and_Feature_Selection.varFeatures
n_genes = params.Normalization_and_Feature_Selection.n_genes
normMethod = params.Normalization_and_Feature_Selection.normMethod
selmethod = params.Normalization_and_Feature_Selection.selmethod






shell:
'''
#!/usr/bin/env Rscript

# libraries
library(Seurat)
library(dplyr)

Data=readRDS("!{seurat_object}")

normMethod <- "!{normMethod}"
selmethod <- "!{selmethod}"
varFeatures <- "!{varFeatures}"
n_genes <- "!{n_genes}"

if (normMethod=="SCT") {
	if (length(unique(Data$sample))==1) {
		Data <- SCTransform(Data,variable.features.n=varFeatures,vars.to.regress=ifelse("percent.mt"%in%colnames(Data@meta.data),
		yes="percent.mt",
		no=NULL),verbose = FALSE)
		
	} else {
		Data <- SplitObject(Data, split.by = "sample")
		Data <- lapply(X = Data, FUN = function(x) {
    		x <- SCTransform(x, verbose = FALSE,variable.features.n=varFeatures,vars.to.regress=ifelse("percent.mt"%in%colnames(Data@meta.data),
		yes="percent.mt",
		no=NULL))
		
		})
		variable.features=SelectIntegrationFeatures(object.list = Data, nfeatures = 3000)
		Data <- merge(Data[[1]],Data[-1])
		VariableFeatures(Data) <- variable
		Data <- ScaleData(Data,vars.to.regress=ifelse("percent.mt"%in%colnames(Data@meta.data),
		yes="percent.mt",
		no=NULL))


	}
} else {
	Data <- NormalizeData(object = Data, normalization.method = normMethod, scale.factor = 10000)
	Data <- FindVariableFeatures(Data, selection.method = selmethod, nfeatures = varFeatures)
	Data <- ScaleData(Data,vars.to.regress=ifelse("percent.mt"%in%colnames(Data@meta.data),
		yes="percent.mt",
		no=NULL))
}

saveRDS(Data,"Normalized.rds")

'''


}


workflow.onComplete {
println "##Pipeline execution summary##"
println "---------------------------"
println "##Completed at: $workflow.complete"
println "##Duration: ${workflow.duration}"
println "##Success: ${workflow.success ? 'OK' : 'failed' }"
println "##Exit status: ${workflow.exitStatus}"
}
