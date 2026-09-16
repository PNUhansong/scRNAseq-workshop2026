# scRNA-seq analysis workshop with Seurat v5

 ## Preparation
### Download dataset
### (1) 검색 후 다운로드
[GEO GSE244515](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE244515)에서 아래 네 sample을 선택해 count matrix를 다운로드합니다. 연구 전체에는 당뇨군도 있지만 이번에는 **Healthy 2명과 Periodontitis 2명**의 PBMC를 사용합니다.

```
GEO database에 공개된 accession number GSE244515 데이터셋을 다운로드합니다.

[https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi)

위 링크 접속 후 GSE244515 검색
```
이번 실습에서는 실습시간을 고려하여 총 4명의 일반인과 치주염 환자로부터 획득한 PBMC 데이터를 사용합니다.
- 'healthy control 1', 'healthy control 2', 'PD1', 'PD2', for practice.
 'custom' 버튼을 이용해 필요한 샘플만 다운로드 받을 수 있습니다.

### (2) 압축 풀기
| 폴더 이름 | GEO sample | 연구 집단 |
|---|---|---|
| H1 | GSM7818495 | Healthy |
| H2 | GSM7818496 | Healthy |
| PD1 | GSM7818506 | Periodontitis |
| PD2 | GSM7818507 | Periodontitis |

`.tar` 파일은 압축을 풀고, **`.gz` 파일은 풀지 않습니다.** 

### (3) 형식에 맞추어 파일 생성
Sample마다 폴더를 만들고 파일 앞의 sample 이름을 지웁니다.
각 폴더에는 **barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz** 세 파일이 있어야 합니다.

```text
data/GSE244515_4sample/
├─ H1/
├─ H2/
├─ PD1/
└─ PD2/
```


## 실습 목표와 실행 방법:
**QC → SCTransform → integration → clustering → UMAP → annotation**을 진행합니다. 이어서 T cell, module score, cell composition을 살펴봅니다.

- 별도 R 파일 없이 이 문서의 코드 블록을 **하나씩 순서대로** 실행합니다.
- 기본 실습은 1–10절, 추가 실습은 11–13절입니다. 마지막에 14절에서 저장합니다.
- `선택`으로 표시한 접힌 부분은 필요할 때만 실행합니다.
- 오류가 나면 다음 블록으로 넘어가지 말고 강사와 함께 확인합니다.
- `<-`는 결과에 이름을 붙이는 기호입니다. 예를 들어 `rawdata`는 합친 데이터를 담는 이름입니다.
- 목표 환경은 **R 4.6.1 / Seurat 5.5.1**입니다. 


## 1. R / RStudio와 패키지 준비


처음 한 번만 설치합니다. 이미 설치했다면 건너뜁니다. 아래 명령은 실행 시점 CRAN 버전을 설치하므로 정확한 버전을 고정하지는 않습니다.

<details>
<summary>의존성 패키지 다운로드</summary>
 
### Rtools 먼저 설치 필요 - 시간 오래 소요
browseURL("https://cran.r-project.org/bin/windows/Rtools/")

### remotes가 설치되어 있지 않다면 먼저 설치 - 시간 오래 소요
```
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
install.packages(c(
  "cluster",
  "cowplot",
  "fastDummies",
  "fitdistrplus",
  "future",
  "future.apply",
  "generics",
  "ggplot2",
  "ggrepel",
  "ggridges",
  "httr",
  "ica",
  "igraph",
  "irlba",
  "jsonlite",
  "KernSmooth",
  "lifecycle",
  "lmtest",
  "MASS",
  "Matrix",
  "matrixStats",
  "miniUI",
  "patchwork",
  "pbapply",
  "plotly",
  "png",
  "progressr",
  "RANN",
  "RColorBrewer",
  "Rcpp",
  "RcppAnnoy",
  "RcppHNSW",
  "reticulate",
  "rlang",
  "ROCR",
  "RSpectra",
  "Rtsne",
  "scales",
  "scattermore",
  "sctransform",
  "shiny",
  "sp",
  "spam",
  "spatstat.explore",
  "spatstat.geom",
  "tibble",
  "uwot"
))

# 원하는 버전으로 설치
remotes::install_version("SeuratObject", version = "5.4.0")
remotes::install_version("Seurat",       version = "5.5.1")
remotes::install_version("sctransform",  version = "0.4.3")
remotes::install_version("patchwork",    version = "1.3.2")
remotes::install_version("scales",       version = "1.4.0")
remotes::install_version("Matrix",       version = "1.7-5")
remotes::install_version("future",       version = "1.70.0")
remotes::install_version("future.apply", version = "1.20.2")
remotes::install_version("uwot",         version = "0.2.4")
remotes::install_version("RcppAnnoy",    version = "0.0.23")
remotes::install_version("irlba",        version = "2.3.7")
remotes::install_version("igraph",       version = "2.3.0")
```
</details>

또는 환경에 따라

```r
install.packages(c("ggplot2", "dplyr", "patchwork"))
```

설치 후 **Session → Restart R**를 선택하고 버전을 확인합니다. Seurat 5.x가 아니거나 강사의 환경과 다르면 먼저 확인합니다.

```r
R.version.string
packageVersion("Seurat")
```


설치가 끝나면 **Session → Restart R**로 R session을 다시 시작합니다. 다음 블록부터 분석을 시작합니다. R을 다시 설치하거나 새 package library를 사용하는 경우에도 이 설치 블록을 사용합니다.


## 2. 분석 준비

패키지를 불러옵니다. R session을 다시 시작하면 이 블록도 다시 실행합니다.

```r
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
set.seed(12345)
```

난수를 사용하는 분석 함수에도 `12345`를 넣습니다. 함수가 자체 seed를 사용하는 경우가 있어 `set.seed()`만으로는 충분하지 않을 수 있습니다. 버전과 환경이 다르면 결과가 완전히 같지는 않을 수 있습니다.

데이터 경로와 저장 폴더를 정합니다.

역슬래쉬('\\')는 슬래쉬('/')로 바꾸어주어야 R에서 인식됩니다.

```r
data_dir <- "  ** 데이터가 저장된 경로 **  "   # 본인 설정에 맞추어 수정. 예시: data_dir <- "C:/Users/Desktop/GSE244515_RAW"
getwd()  # 현재 작업 디렉토리 출력
setwd("  ** 데이터 저장할 경로 ** ")    # 본인 설정에 맞추어 수정. 예시: setwd("C:/Users/")
dir.create("results", showWarnings = FALSE)
```

**확인:** 원하는 폴더에 `results/`가 보이나요?


## 3. 10x matrix에서 Seurat object 만들기

### 3-1. 네 sample 읽기

반복문 대신 sample마다 한 줄씩 읽습니다. 이번 파일은 gene expression matrix를 기준으로 합니다.

```r
H1_counts <- Read10X(file.path(data_dir, "H1"))
H2_counts <- Read10X(file.path(data_dir, "H2"))
PD1_counts <- Read10X(file.path(data_dir, "PD1"))
PD2_counts <- Read10X(file.path(data_dir, "PD2"))
```

행은 gene, 열은 cell barcode입니다.

```r
dim(H1_counts)
H1_counts[1:5, 1:5]

```

### 3-2. Sample별 object 만들기

`min.cells = 3`은 해당 sample에서 최소 3개 cell에 검출된 gene을 남깁니다.

```r
H1 <- CreateSeuratObject(H1_counts, project = "H1", min.cells = 3)
H2 <- CreateSeuratObject(H2_counts, project = "H2", min.cells = 3)
PD1 <- CreateSeuratObject(PD1_counts, project = "PD1", min.cells = 3)
PD2 <- CreateSeuratObject(PD2_counts, project = "PD2", min.cells = 3)
```

### 3-3. 네 sample 합치기

Cell 이름에 sample 접두어를 붙여 같은 barcode를 구분합니다. 아직 normalization 전이므로 count만 합칩니다.

```r
rawdata <- merge(
  H1, y = list(H2, PD1, PD2),
  add.cell.ids = c("H1", "H2", "PD1", "PD2"),
  merge.data = FALSE
)
table(rawdata$orig.ident)
Layers(rawdata[["RNA"]])
```

`orig.ident`는 출신 sample, layer는 sample별 측정값을 보관하는 칸입니다. 질환 정보도 추가합니다.

```r
condition_map <- c(H1 = "Healthy", H2 = "Healthy",
                   PD1 = "Periodontitis", PD2 = "Periodontitis")
rawdata$condition <- unname(condition_map[as.character(rawdata$orig.ident)])
table(rawdata$condition)
```

**확인:** H1·H2·PD1·PD2가 모두 있나요?


## 4. 세포별 QC와 필터링

| 지표 | 의미 |
|---|---|
| nFeature_RNA | 검출 gene 종류 수 |
| nCount_RNA | 전체 UMI 수 |
| percent.mt | mitochondrial gene UMI 비율 |

사람의 mitochondrial gene은 보통 `MT-`로 시작합니다. 아래 첫 결과가 0이면 gene 이름을 강사와 확인한 뒤 진행합니다.

```r
sum(grepl("^MT-", rownames(rawdata)))
rawdata[["percent.mt"]] <- PercentageFeatureSet(rawdata, pattern = "^MT-")
```

### 4-1. 필터링 전 분포

```r
qc_violin <- VlnPlot(
  rawdata, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  group.by = "orig.ident", layer = "counts", pt.size = 0, ncol = 3
)
qc_violin

VlnPlot(
  rawdata, features = c("nFeature_RNA"),
  group.by = "orig.ident", layer = "counts", pt.size = 0
) + geom_hline(yintercept = c(600,5000))

VlnPlot(
  rawdata, features = c("nCount_RNA"),
  group.by = "orig.ident", layer = "counts", pt.size = 0
) + geom_hline(yintercept = c(25000))

VlnPlot(
  rawdata, features = c("percent.mt"),
  group.by = "orig.ident", layer = "counts", pt.size = 0
) + geom_hline(yintercept = c(20))

```

### 4-2. 기준에 맞는 cell 남기기

아래 기준은 이 데이터의 **예시용 기준**입니다. 모든 조직에 그대로 적용하지 않습니다. UMI가 높다는 이유만으로 doublet을 확정할 수 없고, mitochondrial 비율도 cell type에 따라 달라집니다.

```r
filtdata <- subset(
  rawdata,
  subset = nFeature_RNA > 600 & nFeature_RNA < 5000 &
    nCount_RNA < 25000 & percent.mt < 20
)
```

Sample별 남은 cell 수를 비교합니다.

```r
qc_summary <- data.frame(
  before = table(factor(rawdata$orig.ident, levels = names(condition_map))),
  after = as.vector(table(factor(filtdata$orig.ident, levels = names(condition_map))))
)
colnames(qc_summary) <- c("sample", "before", "after")
qc_summary
```

**확인:** 특정 sample만 크게 줄었나요? Sample이 사라지거나 cell이 수십 개만 남았다면 다음 계산 전에 강사와 확인합니다. 이 필터만으로 doublet과 ambient RNA가 모두 제거되지는 않습니다.

## 5. SCTransform과 PCA - 약 5분 소요

(지금 Sample 별로 layer가 분리되어 있으므로 생략 가능) Sample별 layer를 준비합니다.

```r
filtdata[["RNA"]] <- JoinLayers(filtdata[["RNA"]])
filtdata[["RNA"]] <- split(filtdata[["RNA"]], f = filtdata$orig.ident)
Layers(filtdata[["RNA"]])
```

SCTransform으로 기술적 변이를 모델링하고 PCA로 주요 발현 차이를 요약합니다. Seurat v5 기본값인 SCT v2, variable genes 3,000개, PCA 50개는 코드에서 생략했습니다. `conserve.memory = TRUE`는 메모리 사용을 줄이기 위해 남겼습니다.

```r
filtdata <- SCTransform(filtdata, conserve.memory = TRUE, seed.use = 12345)
filtdata <- RunPCA(filtdata, seed.use = 12345)
ElbowPlot(filtdata, ndims = 50)
```

**확인:** ElbowPlot이 보이나요? 이번 실습은 이후 계산에 PC 1–30을 사용합니다.

<img width="400" height="400" alt="image" src="https://github.com/user-attachments/assets/19c502f4-9dd4-4adc-a922-638ee38dc476" />

### 데이터 저장: 
데이터의 크기가 크고 스크립트 실행시간이 길다면 실행된 데이터는 만약을 위해 **꼭 저장**하는 습관을 길러야합니다.
만약을 위해 github 페이지에 filtdata.RData 파일을 업로드해두었습니다.
(단, 세포 수를 조절해 용량을 줄인 데이터를 업로드 하였음.)
```
#save(filtdata, file = ' ** 원하는 경로 ** ')
#load(' ** 저장한 경로 ** ')
save(filtdata, file = 'C:/Users/results/filtdata.RData')
load('C:/Users/results/filtdata.RData')
```


## 6. Integration 전 UMAP - 생략

```r
filtdata <- RunUMAP(
  filtdata, dims = 1:30,
  reduction.name = "umap.unintegrated", seed.use = 12345
)
umap_before <- DimPlot(filtdata, reduction = "umap.unintegrated", group.by = "orig.ident")
umap_before
```

점 하나는 cell이고 색은 sample입니다. Sample별 분리는 기술적·생물학적 차이 모두에서 생길 수 있습니다. UMAP 축은 실제 조직 좌표가 아닙니다.

<img width="600" height="600" alt="image" src="https://github.com/user-attachments/assets/184a58ca-5009-41bc-a588-e8925dc221fb" />


## 7. RPCA integration

Sample 사이에 공유되는 구조를 정렬합니다. `SCT`를 사용했다는 설정과 결과 이름은 분석에 필요하므로 남깁니다.

```r
intdata <- IntegrateLayers(
  filtdata, method = RPCAIntegration,
  assay = "SCT", normalization.method = "SCT",
  new.reduction = "integrated.rpca", dims = 1:30
)
Reductions(intdata)
```

**확인:** `integrated.rpca`가 있나요? 원래 RNA count를 지우는 과정은 아닙니다.

## 8. Clustering과 integration 후 UMAP

### 8-1. 비슷한 이웃끼리 묶기

`resolution`은 cluster를 나누는 세밀함에 영향을 줍니다. 이번에는 `0.2`를 사용합니다.

```r
intdata <- FindNeighbors(intdata, reduction = "integrated.rpca", dims = 1:30)
intdata <- FindClusters(intdata, resolution = 0.2,
                        cluster.name = "rpca_clusters", random.seed = 12345)
```

### 8-2. 지도 그리기

```r
intdata <- RunUMAP(
  intdata, reduction = "integrated.rpca", dims = 1:30,
  reduction.name = "umap.rpca", seed.use = 12345
)
```

Sample 색과 cluster 색으로 나란히 봅니다.

```r
umap_sample <- DimPlot(intdata, reduction = "umap.rpca", group.by = "orig.ident")
umap_cluster <- DimPlot(intdata, reduction = "umap.rpca",
                        group.by = "rpca_clusters", label = TRUE)
umap_after <- umap_sample + umap_cluster
umap_after
```

**질문:** Sample이 섞였나요? 다음 절에서 marker도 유지되는지 확인합니다. Cluster 번호는 계산 결과의 이름표입니다.

## 9. Marker로 cell type 후보 찾기

### 9-1. RNA 발현값 준비

지도에는 integrated reduction을, marker 확인에는 RNA normalized data를 사용합니다.

```r
intdata[["RNA"]] <- JoinLayers(intdata[["RNA"]])
DefaultAssay(intdata) <- "RNA"
intdata <- NormalizeData(intdata)
Idents(intdata) <- "rpca_clusters"
```

### 9-2. Marker 목록

```r
marker_panels <- list(
  `T cell` = c("CD3D", "CD3E", "TRAC"),
  `CD4 T candidate` = c("IL7R", "LTB", "CCR7"),
  `CD8 T / NK` = c("CD8A", "NKG7", "GNLY", "PRF1"),
  `B cell` = c("MS4A1", "CD79A", "CD37"),
  Monocyte = c("LYZ", "S100A8", "FCGR3A", "LILRB1"),
  mDC = c("CD1C", "CLEC10A", "FCER1A"),
  pDC = c("GZMB", "IL3RA", "CLEC4C"),
  Platelet = c("PPBP", "PF4"),
  Erythroid = c("HBA1", "HBB"),
  Proliferating = c("MKI67", "TOP2A")
)
```

### 9-3. DotPlot 읽기

```r
marker_dotplot <- DotPlot(intdata, features = marker_panels) + RotatedAxis()
marker_dotplot
```

점 크기는 발현 cell 비율, 색은 gene별로 표준화한 cluster 평균입니다. 다른 gene끼리 색만 보고 절대 발현량을 비교하지 않습니다. Marker가 없다는 경고가 나오면 gene 이름을 확인합니다.

**활동:** 각 cluster의 후보 이름과 근거 marker 두 개를 적어봅니다. IL7R 하나로 CD4 T를, NKG7 하나로 NK를 확정하지 않습니다.

<details>
<summary>선택: cluster별 후보 marker 계산하기</summary>

시간이 남을 때만 실행합니다. 질환군 간 검정을 대신하는 분석은 아닙니다. `min.pct`와 `logfc.threshold`는 후보 선정 기준이므로 이전 실습값을 명시합니다.

```r
cluster_markers <- FindAllMarkers(intdata, only.pos = TRUE,
                                 min.pct = 0.25, logfc.threshold = 0.25,
                                 random.seed = 12345)
top_markers <- cluster_markers |>
  group_by(cluster) |>
  slice_max(avg_log2FC, n = 10, with_ties = FALSE)
top_markers
```

```r
write.csv(cluster_markers, "results/cluster_markers_all.csv", row.names = FALSE)
write.csv(top_markers, "results/cluster_markers_top10.csv", row.names = FALSE)
```

</details>

## 10. 직접 annotation 붙이기

먼저 모든 cluster를 `Unassigned`로 둡니다. **이 초기화 블록은 한 번만 실행합니다.**

```r
cluster_ids <- levels(Idents(intdata))
celltype_map <- setNames(rep("Unassigned", length(cluster_ids)), cluster_ids)
celltype_map
```

아래는 **작성 형식 예시이며 정답이 아닙니다.** Marker를 보고 번호와 이름을 수정한 뒤 앞의 `#`를 지워 실행합니다. 불확실한 cluster는 `Unassigned`로 남겨도 됩니다.

```r
# celltype_map["0"] <- "T cell"
# celltype_map["1"] <- "Monocyte"
```

이름표를 적용합니다. 이름을 수정했다면 아래 블록부터 다시 실행합니다.

```r
intdata$celltype <- unname(celltype_map[as.character(intdata$rpca_clusters)])
table(intdata$rpca_clusters, intdata$celltype)
annotated_umap <- DimPlot(intdata, reduction = "umap.rpca",
                          group.by = "celltype", label = TRUE)
annotated_umap
```

**확인:** 표에 이름이 올바르게 연결되었나요? 모두 `Unassigned`라면 marker를 다시 확인합니다. Annotation을 바꾼 뒤에는 추가 분석도 다시 계산합니다.

## 11. 추가 실습: T cell만 확대하기

### 11-1. T cell 선택

`T_CELL_LABELS`를 10절에서 사용한 실제 이름과 맞춥니다. 이름이 없으면 10절로 돌아갑니다. 이 절을 생략할 때는 12절도 건너뛰고 13절로 이동합니다.

```r
T_CELL_LABELS <- c("T cell", "CD4 T", "CD8 T", "Treg")
table(intdata$celltype)
tcell <- subset(intdata, subset = celltype %in% T_CELL_LABELS)
table(tcell$orig.ident)
```

Sample이 하나뿐이거나 cell이 매우 적으면 강사와 확인합니다. 전체 PBMC의 PCA를 확대하는 대신 T cell 내부의 차이를 다시 계산합니다.

```r
tcell[["RNA"]] <- JoinLayers(tcell[["RNA"]])
tcell[["RNA"]] <- split(tcell[["RNA"]], f = tcell$orig.ident)
tcell <- SCTransform(tcell, new.assay.name = "SCT_T",
                      conserve.memory = TRUE, seed.use = 12345)
tcell <- RunPCA(tcell, npcs = 30, reduction.name = "pca.tcell", seed.use = 12345)
```

### 11-2. 다시 통합하고 묶기

T cell에서는 PC 1–20, resolution 0.4를 사용합니다. 작은 subset을 고려해 `k.weight = 50`을 유지합니다.

```r
tcell <- IntegrateLayers(
  tcell, method = RPCAIntegration, assay = "SCT_T",
  normalization.method = "SCT", orig.reduction = "pca.tcell",
  new.reduction = "integrated.rpca.tcell", dims = 1:20, k.weight = 50
)
```

```r
tcell <- FindNeighbors(tcell, reduction = "integrated.rpca.tcell", dims = 1:20)
tcell <- FindClusters(tcell, resolution = 0.4,
                      cluster.name = "tcell_clusters", random.seed = 12345)
tcell <- RunUMAP(tcell, reduction = "integrated.rpca.tcell", dims = 1:20,
                 reduction.name = "umap.tcell", seed.use = 12345)
tcell_umap <- DimPlot(tcell, reduction = "umap.tcell",
                      group.by = "tcell_clusters", label = TRUE)
tcell_umap
```

### 11-3. Type과 state 구분하기

Marker를 볼 RNA 발현값을 준비합니다.

```r
tcell[["RNA"]] <- JoinLayers(tcell[["RNA"]])
DefaultAssay(tcell) <- "RNA"
tcell <- NormalizeData(tcell)
```

```r
tcell_markers <- list(
  `Naive / memory-like` = c("CCR7", "TCF7", "LEF1", "IL7R", "LTB"),
  Cytotoxicity = c("CD8A", "CCL5", "NKG7", "PRF1", "GZMB"),
  Treg = c("FOXP3", "IL2RA", "CTLA4"),
  Proliferating = c("MKI67", "TOP2A", "STMN1"),
  `IFN response` = c("ISG15", "IFIT1", "IFIT3", "MX1"),
  `Exhaustion-associated` = c("PDCD1", "LAG3", "HAVCR2", "TOX", "TIGIT")
)
tcell_dotplot <- DotPlot(tcell, features = tcell_markers,
                         group.by = "tcell_clusters") + RotatedAxis()
tcell_dotplot
```

**질문:** IFN response는 여러 subtype에 걸쳐 나타나나요? PDCD1·TIGIT 하나만으로 exhaustion을 확정할 수는 없습니다. TIGIT는 Treg에서도 나타납니다. 여러 marker와 질환 맥락을 함께 봅니다.

## 12. 추가 실습: Module score

11절에서 만든 `tcell`이 필요합니다. 여러 gene의 발현 경향을 상대 점수로 요약합니다.

```r
tcell_programs <- list(
  Cytotoxicity = c("NKG7", "CCL5", "PRF1", "GZMB"),
  IFN_response = c("ISG15", "IFIT1", "IFIT3", "MX1"),
  Exhaustion_associated = c("PDCD1", "LAG3", "HAVCR2", "TOX", "TIGIT")
)
```

데이터에 없는 gene을 제외하고 각 program에 남은 gene 수를 확인합니다. `lapply`는 목록 각각에 같은 처리를 적용합니다. **2개 미만인 항목이 있으면 계산 전에 강사와 확인합니다.**

```r
tcell_programs <- lapply(tcell_programs, intersect, y = rownames(tcell))
lengths(tcell_programs)
tcell <- AddModuleScore(tcell, features = tcell_programs,
                        name = "Program", seed = 12345)
```

목록 순서에 따라 `Program1`, `Program2`, `Program3`이 생깁니다. 반복문 대신 이름을 한 줄씩 붙입니다.

```r
tcell$Cytotoxicity_score <- tcell$Program1
tcell$IFN_response_score <- tcell$Program2
tcell$Exhaustion_associated_score <- tcell$Program3
```

```r
program_plot <- FeaturePlot(
  tcell, reduction = "umap.tcell",
  features = c("Cytotoxicity_score", "IFN_response_score", "Exhaustion_associated_score"),
  ncol = 3
)
program_plot
```

점수는 control gene과 비교한 상대값입니다. 그림마다 색 범위가 다를 수 있으며, 다른 program끼리 점수 크기를 직접 비교하지 않습니다. 각 program 내 패턴을 marker·sample 정보와 함께 해석합니다.

T cell 추가 실습을 마쳤다면 저장합니다.

```r
saveRDS(tcell, "results/GSE244515_Tcell_Seurat5.rds")
```

## 13. 추가 실습: Sample별 cell composition

10절까지의 `intdata`만 있으면 됩니다. **QC 후 남은 PBMC 전체**를 sample별 분모로 사용하며 `Unassigned`도 포함합니다.

```r
composition_table <- as.data.frame(table(
  sample = intdata$orig.ident, celltype = intdata$celltype
))
colnames(composition_table)[3] <- "cell_count"
composition_table <- composition_table |>
  group_by(sample) |>
  mutate(cell_proportion = cell_count / sum(cell_count)) |>
  ungroup()
composition_table$condition <- unname(condition_map[as.character(composition_table$sample)])
composition_table
```

막대 하나는 sample 하나입니다. 0.5는 50%를 뜻합니다.

```r
composition_plot <- ggplot(composition_table,
                           aes(sample, cell_proportion, fill = celltype)) +
  geom_col() +
  theme_classic() +
  labs(x = "Sample", y = "Cell proportion")
composition_plot
```

**질문:** PD1과 PD2에서 같은 방향의 차이가 보이나요? Healthy 2명과 Periodontitis 2명의 탐색적 결과입니다. Cell 수천 개가 독립 환자 수천 명을 뜻하지 않습니다. 채취·분리·QC에 따른 편향도 고려합니다.

## 14. 결과 저장과 재개

기본 실습 결과와 실행 환경을 저장합니다. 같은 파일 이름으로 저장하면 이전 결과를 갱신합니다.

```r
saveRDS(intdata, "results/GSE244515_4sample_Seurat5.rds")
writeLines(capture.output(sessionInfo()), "results/sessionInfo.txt")
```

그림은 RStudio의 **Plots → Export**로 저장할 수 있습니다. 아래 자동 저장은 선택입니다.

<details>
<summary>선택: 기본 실습 그림과 표 저장</summary>

1–10절을 마친 뒤 실행합니다. Width와 height는 inch 단위의 그림 크기입니다.

```r
ggsave("results/01_QC_violin_before_filtering.png", qc_violin, width = 14, height = 5)
ggsave("results/02_QC_scatter_before_filtering.png", qc_scatter, width = 10, height = 8)
ggsave("results/03_UMAP_before_integration.png", umap_before, width = 8, height = 6)
ggsave("results/04_UMAP_after_integration.png", umap_after, width = 14, height = 6)
ggsave("results/05_marker_DotPlot.png", marker_dotplot, width = 15, height = 8)
ggsave("results/06_UMAP_annotated.png", annotated_umap, width = 9, height = 7)
```

```r
write.csv(qc_summary, "results/QC_cell_numbers.csv", row.names = FALSE)
annotation_table <- data.frame(cluster = names(celltype_map), celltype = unname(celltype_map))
write.csv(annotation_table, "results/cluster_annotation.csv", row.names = FALSE)
```

</details>

<details>
<summary>선택: T cell 그림 저장 — 11–12절 완료 후</summary>

```r
ggsave("results/08_Tcell_UMAP.png", tcell_umap, width = 9, height = 7)
ggsave("results/09_Tcell_marker_DotPlot.png", tcell_dotplot, width = 15, height = 8)
ggsave("results/10_Tcell_program_scores.png", program_plot, width = 15, height = 5)
```

</details>

<details>
<summary>선택: Composition 저장 — 13절 완료 후</summary>

```r
ggsave("results/07_celltype_composition.png", composition_plot, width = 9, height = 6)
write.csv(composition_table, "results/celltype_composition_by_sample.csv", row.names = FALSE)
```

</details>

<details>
<summary>선택: 새 R session에서 저장한 object 불러오기</summary>

같은 Project를 열고 2절의 패키지를 불러온 뒤 실행합니다. 처음 분석할 때는 생략합니다.

```r
intdata <- readRDS("results/GSE244515_4sample_Seurat5.rds")
condition_map <- c(H1 = "Healthy", H2 = "Healthy",
                   PD1 = "Periodontitis", PD2 = "Periodontitis")
table(intdata$celltype)
```

11절 또는 13절로 이동합니다. 모두 `Unassigned`라면 9–10절에서 먼저 annotation을 붙입니다.

</details>

## 막혔을 때 확인하기

- **파일을 못 찾음:** Project 위치와 sample 폴더, 세 파일의 이름을 확인합니다. `.gz`는 풀지 않습니다.
- **object가 없음:** 앞 블록을 실행했는지 확인합니다. Session을 재시작하면 메모리의 object가 사라집니다.
- **MT- gene이 0개:** 사람 gene symbol 대신 다른 ID를 읽었는지 확인합니다.
- **Read10X 결과가 matrix가 아닌 list:** 이번 실습용 Gene Expression 파일을 받았는지 강사와 확인합니다.
- **메모리 부족:** 새 session에서 필요한 단계만 실행하거나 검증된 checkpoint를 사용합니다. `conserve.memory = TRUE`를 유지합니다.
- **data layers are not joined:** 9절의 RNA `JoinLayers()`를 실행합니다. SCT assay에는 그대로 적용하지 않습니다.
- **anchor / k.weight 오류:** Sample별 남은 cell 수와 공유 population을 먼저 확인합니다. 강사와 `k.weight`를 조정한 뒤 marker 구조도 다시 확인합니다.
- **T cell을 선택하지 못함:** 실제 annotation 이름과 `T_CELL_LABELS`를 맞춥니다. NK를 cytotoxic marker 하나로 T cell에 포함하지 않습니다.
- **module score gene이 부족함:** `lengths(tcell_programs)`와 gene 이름을 확인합니다.

## 참고 자료

- [Seurat v5 integration](https://satijalab.org/seurat/articles/seurat5_integration.html)
- [SCTransform 기본값](https://satijalab.org/seurat/reference/sctransform)
- [RunPCA 기본값과 seed](https://satijalab.org/seurat/reference/runpca)
- [RPCAIntegration](https://satijalab.org/seurat/reference/rpcaintegration)
- [AddModuleScore](https://satijalab.org/seurat/reference/addmodulescore)
- [GSE244515](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE244515)
- Lee H, Joo J, Song J, et al. *Immunological link between periodontitis and type 2 diabetes deciphered by single-cell RNA analysis*. Clinical and Translational Medicine. 2023. [doi:10.1002/ctm2.1503](https://doi.org/10.1002/ctm2.1503)
