# scRNA-seq analysis workshop with Seurat v5

사람 PBMC 4개 sample로 **QC → SCTransform v2 → RPCA integration → clustering → UMAP → marker annotation**을 실습합니다. 이어서 T cell을 확대하고, functional program과 sample별 cell composition을 살펴봅니다.

**이 README의 R 코드 블록을 위에서 아래로 실행합니다. 별도의 `.R` 파일이나 `source()` 명령은 필요하지 않습니다.** 데이터 파일은 별도로 준비합니다. 코드는 한 블록씩 RStudio Console에 붙여넣거나, RStudio에서 새 R Script를 열어 복사한 뒤 실행해도 됩니다.

- 워크숍 목표 환경: **R 4.6.1 / Seurat 5.5.1**. 아래 설치 코드는 Seurat 5.5.1 이상인 **5.x**를 확인하며 정확한 버전을 고정하는 설치는 아닙니다.
- `install.packages()`는 실행 시점의 CRAN 버전을 설치하므로, 수업 전에 강사와 학생의 설치 버전을 확인합니다.
- 기본 실습: 1–10절. 추가 실습: 11–13절. 마지막에 14절로 결과를 저장합니다.
- **10절 annotation에서는 잠시 멈추고 marker를 해석합니다.** 나머지 계산은 순서대로 실행하되, cluster 번호별 정답을 미리 복사하지 않습니다.
- 결과는 `results/`에 저장합니다. 같은 이름으로 다시 실행하면 해당 결과 파일을 갱신합니다.

각 단계의 완료 기준을 확인한 뒤 다음 단계로 이동합니다. 계산 시간은 cell 수와 컴퓨터 성능에 따라 달라집니다.


## 1. R / RStudio와 패키지 준비

RStudio에서 **File → New Project → Existing Directory**로 실습 폴더를 선택합니다. `getwd()`로 현재 폴더를 확인할 수 있습니다.

아래 설치 블록은 처음 한 번 실행합니다. 설치 중에는 분석 코드를 실행하지 않습니다.
```r
# scRNA-seq workshop: package installation
# Target environment: R 4.6.1, Seurat 5.5.1 or later in the Seurat 5 series

cran_repo <- "https://cloud.r-project.org"
required_packages <- c("Seurat", "ggplot2", "dplyr", "patchwork", "scales", "future")

if (getRversion() < package_version("4.6.1")) {
  stop(
    "이 실습은 R 4.6.1 이상을 기준으로 작성했습니다. 현재 R 버전: ",
    R.version.string
  )
}

installed <- installed.packages()
packages_to_install <- setdiff(required_packages, rownames(installed))

# Seurat가 설치되어 있어도 실습 기준보다 오래된 버전이면 업데이트합니다.
if ("Seurat" %in% rownames(installed) &&
    package_version(installed["Seurat", "Version"]) < package_version("5.5.1")) {
  packages_to_install <- c(packages_to_install, "Seurat")
}

if (length(packages_to_install) > 0) {
  install.packages(unique(packages_to_install), repos = cran_repo)
}

seurat_version <- packageVersion("Seurat")

if (seurat_version < package_version("5.5.1") ||
    seurat_version >= package_version("6.0.0")) {
  stop(
    "Seurat 5.5.1 이상인 Seurat 5 버전을 설치해 주세요. 현재 버전: ",
    as.character(seurat_version)
  )
}

cat("설치 확인 완료\n")
cat("R:", R.version.string, "\n")
cat("Seurat:", as.character(packageVersion("Seurat")), "\n")
cat("SeuratObject:", as.character(packageVersion("SeuratObject")), "\n")
```

설치가 끝나면 **Session → Restart R**로 R session을 다시 시작합니다. 다음 블록부터 분석을 시작합니다. R을 다시 설치하거나 새 package library를 사용하는 경우에도 이 설치 블록을 사용합니다.

## 2. 데이터와 분석 설정

[GSE244515](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE244515)는 치주염과 당뇨를 다룬 연구입니다. 이번 실습은 그중 **Healthy 2개와 Periodontitis 2개**만 사용합니다. 당뇨군 전체를 비교하는 실습은 아닙니다.

| 실습 이름 | GEO sample | 실습 condition |
|---|---|---|
| H1 | GSM7818495 | Healthy |
| H2 | GSM7818496 | Healthy |
| PD1 | GSM7818506 | Periodontitis |
| PD2 | GSM7818507 | Periodontitis |

실습 폴더 아래에 다음 구조로 파일을 준비합니다.

```text
data/
└─ GSE244515_4sample/
   ├─ H1/
   ├─ H2/
   ├─ PD1/
   └─ PD2/
```

각 sample 폴더에는 `barcodes.tsv.gz`, `features.tsv.gz`, `matrix.mtx.gz`가 있어야 합니다. GEO 파일 앞의 sample 접두어를 제거하여 이름을 맞추고 **`.gz`는 풀지 않습니다.** 예를 들어 `GSM7818495_H1_matrix.mtx.gz`는 `H1/matrix.mtx.gz`로 둡니다.

아래 코드는 패키지를 불러오고, 경로·QC 기준·분석 차원을 정의합니다. R session을 다시 시작했다면 이 블록부터 필요한 단계를 다시 실행합니다.
```r
# 실습 중 같은 결과를 얻을 수 있도록 난수 시작점을 고정합니다.
set.seed(20260909)

# Seurat v5 형식의 Assay를 만들도록 설정합니다.
options(Seurat.object.assay.version = "v5")

# 학생 PC에서는 순차 계산을 사용합니다. maxSize는 RAM을 늘리는 설정이 아닙니다.
future::plan(future::sequential)
options(future.globals.maxSize = 8 * 1024^3)

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

if (getRversion() < package_version("4.6.1")) {
  stop("R 4.6.1 이상이 필요합니다. 현재 버전: ", R.version.string)
}

seurat_version <- packageVersion("Seurat")

if (seurat_version < package_version("5.5.1") ||
    seurat_version >= package_version("6.0.0")) {
  stop(
    "Seurat 5.5.1 이상인 Seurat 5 버전이 필요합니다. 현재 버전: ",
    as.character(seurat_version)
  )
}

cat("R:", R.version.string, "\n")
cat("Seurat:", as.character(packageVersion("Seurat")), "\n")


# 0. 경로와 분석 설정 ---------------------------------------------------------

# RStudio Project의 최상위 폴더에서 이 스크립트를 실행한다고 가정합니다.
# 데이터는 data/GSE244515_4sample/H1 같은 구조로 둡니다.
data_dir <- file.path("data", "GSE244515_4sample")
result_dir <- "results"
dir.create(result_dir, showWarnings = FALSE, recursive = TRUE)

sample_dirs <- c(
  H1 = file.path(data_dir, "H1"),
  H2 = file.path(data_dir, "H2"),
  PD1 = file.path(data_dir, "PD1"),
  PD2 = file.path(data_dir, "PD2")
)

# 이 cutoff는 네 샘플을 이용한 교육용 출발점입니다.
# 다른 데이터에는 그대로 복사하지 말고 샘플별 분포와 세포 유형을 확인해야 합니다.
min_features <- 600
max_features <- 5000
max_counts <- 25000
max_percent_mt <- 20

# PCA와 통합 분석에서 사용할 차원입니다.
dims_use <- 1:30
```

완료 기준: R / Seurat 버전이 출력되고, 실습 폴더에 `results/`가 생깁니다.

## 3. 10x matrix에서 Seurat object 만들기

행은 gene, 열은 cell barcode입니다. 세 파일이 모두 있는지 확인한 뒤 sample별 object를 만들고 합칩니다. Cell 이름 앞에 sample 이름을 붙여 서로 같은 barcode를 구분합니다.

`min.cells = 3`은 sample 내 최소 3개 cell에서 검출된 gene을 남기는 설정입니다. 세포별 QC는 다음 절에서 수행합니다.
```r
# 1. 10x count matrix 읽기 ----------------------------------------------------

required_files <- c("barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz")

for (sample_id in names(sample_dirs)) {
  sample_path <- sample_dirs[[sample_id]]
  missing_files <- required_files[
    !file.exists(file.path(sample_path, required_files))
  ]

  if (length(missing_files) > 0) {
    stop(
      sample_id, " 폴더에서 다음 파일을 찾을 수 없습니다: ",
      paste(missing_files, collapse = ", "),
      "\n확인한 경로: ", normalizePath(sample_path, mustWork = FALSE)
    )
  }
}

read_gene_expression <- function(data_dir) {
  counts <- Read10X(data.dir = data_dir)

  # Feature Barcode 자료는 Gene Expression과 다른 feature type을 list로 돌려줄 수 있습니다.
  if (is.list(counts)) {
    if (!"Gene Expression" %in% names(counts)) {
      stop("Read10X 결과에 'Gene Expression' matrix가 없습니다: ", data_dir)
    }
    counts <- counts[["Gene Expression"]]
  }

  counts
}

counts_list <- lapply(sample_dirs, read_gene_expression)

# 각 count matrix에서 gene이 행, cell barcode가 열에 놓입니다.
cat("H1 count matrix 크기 (gene x cell):\n")
print(dim(counts_list$H1))
print(counts_list$H1[1:5, 1:5])

sample_objects <- Map(
  f = function(counts, sample_id) {
    CreateSeuratObject(
      counts = counts,
      project = sample_id,
      min.cells = 3,
      min.features = 0
    )
  },
  counts = counts_list,
  sample_id = names(counts_list)
)

rawdata <- merge(
  x = sample_objects[[1]],
  y = sample_objects[-1],
  add.cell.ids = names(sample_objects),
  project = "GSE244515_4sample",
  merge.data = FALSE
)

stopifnot(anyDuplicated(colnames(rawdata)) == 0)
condition_map <- c(H1 = "Healthy", H2 = "Healthy",
                   PD1 = "Periodontitis", PD2 = "Periodontitis")
rawdata$condition <- unname(condition_map[as.character(rawdata$orig.ident)])
stopifnot(!anyNA(rawdata$condition))

cat("합친 object 크기 (gene x cell):\n")
print(dim(rawdata))
print(table(rawdata$orig.ident))

# Seurat v5에서는 한 assay 안에 샘플별 count layer를 보관할 수 있습니다.
cat("RNA assay의 layer:\n")
print(Layers(rawdata[["RNA"]]))
```

완료 기준: H1 matrix 크기, 네 sample의 cell 수, RNA layer 이름이 출력됩니다. Layer는 같은 assay 안에서 sample별 측정값을 나누어 보관하는 칸입니다.

## 4. 세포별 QC와 필터링

| 지표 | 뜻 | 함께 생각할 점 |
|---|---|---|
| nFeature_RNA | 검출 gene 종류 수 | 너무 적으면 정보 부족, 유독 많으면 doublet 가능성 |
| nCount_RNA | 전체 UMI 수 | 높은 값만으로 doublet을 확정할 수 없음 |
| percent.mt | mitochondrial gene UMI 비율 | stress·RNA 손실의 단서, cell type별 차이도 고려 |

교육용 시작 기준은 **600 < nFeature < 5000, nCount < 25000, percent.mt < 20**입니다. 다른 조직이나 데이터의 정답으로 사용하지 않습니다. 아래 블록은 QC 그림을 그리고 `filtdata`를 실제로 만듭니다.
```r
# 2. Quality control ----------------------------------------------------------

# 사람의 미토콘드리아 유전자 이름은 보통 MT-로 시작합니다.
mt_genes <- grep("^MT-", rownames(rawdata), value = TRUE)
if (length(mt_genes) == 0) {
  stop("MT- gene이 없습니다. features.tsv.gz가 사람 gene symbol을 사용하는지 확인하세요.")
}
rawdata[["percent.mt"]] <- PercentageFeatureSet(rawdata, features = mt_genes)

qc_violin <- VlnPlot(
  rawdata,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  group.by = "orig.ident",
  layer = "counts",
  pt.size = 0,
  ncol = 3
)
print(qc_violin)
ggsave(
  filename = file.path(result_dir, "01_QC_violin_before_filtering.png"),
  plot = qc_violin,
  width = 14,
  height = 5,
  dpi = 200
)

qc_scatter <- rawdata[[]] |>
  ggplot(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt)) +
  geom_point(size = 0.35, alpha = 0.55) +
  scale_color_gradient(low = "grey80", high = "#0055FF") +
  facet_wrap(vars(orig.ident), ncol = 2) +
  geom_hline(yintercept = c(min_features, max_features), linetype = "dashed") +
  geom_vline(xintercept = max_counts, linetype = "dashed") +
  theme_classic() +
  labs(
    title = "세포별 QC 지표",
    subtitle = "점 하나는 cell barcode 하나입니다",
    color = "percent.mt"
  )
print(qc_scatter)
ggsave(
  filename = file.path(result_dir, "02_QC_scatter_before_filtering.png"),
  plot = qc_scatter,
  width = 10,
  height = 8,
  dpi = 200
)

qc_before <- rawdata[[]] |>
  count(orig.ident, name = "cells_before")

filtdata <- subset(
  rawdata,
  subset = nFeature_RNA > min_features &
    nFeature_RNA < max_features &
    nCount_RNA < max_counts &
    percent.mt < max_percent_mt
)

qc_after <- filtdata[[]] |>
  count(orig.ident, name = "cells_after")

qc_summary <- left_join(qc_before, qc_after, by = "orig.ident") |>
  mutate(
    cells_after = coalesce(cells_after, 0L),
    retained_percent = round(100 * cells_after / cells_before, 1)
  )

print(qc_summary)
write.csv(
  qc_summary,
  file = file.path(result_dir, "QC_cell_numbers.csv"),
  row.names = FALSE
)

# 주의: 이 기본 filtering만으로 doublet과 ambient RNA 문제가 모두 해결되지는 않습니다.

if (any(qc_summary$cells_after <= 51)) {
  stop("QC 후 너무 적은 cell이 남은 sample이 있습니다. QC 표와 데이터 입력을 먼저 확인하세요.")
}
```

완료 기준: QC 그림 두 개와 `QC_cell_numbers.csv`가 저장됩니다. Sample별로 몇 %가 남았는지 비교하세요. 이 필터링만으로 doublet과 ambient RNA 문제가 모두 해결되지는 않습니다.

## 5. SCTransform v2와 PCA

RNA layer를 sample별로 나눈 뒤 SCTransform으로 기술적 변이를 모델링합니다. Variable gene을 선택하고 PCA로 주요 발현 변이를 요약합니다. 모든 세포의 발현을 같게 만드는 과정은 아닙니다.
```r
# 3. Seurat v5 layer 준비와 SCTransform v2 ------------------------------------

DefaultAssay(filtdata) <- "RNA"

# merge 결과의 layer 이름이 환경에 따라 달라도, 한 번 합친 뒤 sample별로 다시 나누면
# 이후의 normalization과 integration에서 각 sample을 독립된 batch로 인식할 수 있습니다.
filtdata[["RNA"]] <- JoinLayers(filtdata[["RNA"]])
filtdata[["RNA"]] <- split(filtdata[["RNA"]], f = filtdata$orig.ident)

cat("sample별로 나눈 RNA layer:\n")
print(Layers(filtdata[["RNA"]]))

# Seurat v5에서는 SCT v2가 기본입니다.
# SCTransform은 normalization, variance stabilization, variable feature 선택을 수행합니다.
filtdata <- SCTransform(
  object = filtdata,
  assay = "RNA",
  new.assay.name = "SCT",
  vst.flavor = "v2",
  variable.features.n = 3000,
  conserve.memory = TRUE,
  verbose = FALSE
)

filtdata <- RunPCA(
  object = filtdata,
  assay = "SCT",
  npcs = 50,
  verbose = FALSE
)

print(ElbowPlot(filtdata, ndims = 50))
```

완료 기준: `SCT` assay와 `pca` reduction이 생성되고 ElbowPlot이 보입니다. 여기서는 30개 PC를 사용하지만 실제 분석에서는 데이터의 구조를 함께 점검합니다.

## 6. Integration 전 UMAP

통합 전 지도를 남겨 이후 결과와 비교합니다. 점 하나는 cell 하나이고 색은 sample입니다. Sample별 분리는 기술적 차이뿐 아니라 생물학적 차이에서도 생길 수 있습니다.
```r
# 4. Integration 전 결과 ------------------------------------------------------

# integration 전 UMAP을 남겨 두면 sample별 batch effect가 얼마나 보이는지 비교할 수 있습니다.
filtdata <- RunUMAP(
  object = filtdata,
  reduction = "pca",
  dims = dims_use,
  reduction.name = "umap.unintegrated",
  reduction.key = "UMAPunintegrated_",
  seed.use = 20260909,
  verbose = FALSE
)

umap_before <- DimPlot(
  filtdata,
  reduction = "umap.unintegrated",
  group.by = "orig.ident"
) +
  ggtitle("Integration 전")
print(umap_before)
ggsave(
  filename = file.path(result_dir, "03_UMAP_before_integration.png"),
  plot = umap_before,
  width = 8,
  height = 6,
  dpi = 200
)
```

완료 기준: `03_UMAP_before_integration.png`가 저장됩니다. UMAP 축은 gene 발현량이나 실제 조직 좌표가 아닙니다.

## 7. RPCA integration

Sample 사이에 공유되는 구조를 정렬합니다. 결과는 `integrated.rpca` reduction에 보관하고, 원래 RNA count는 유지합니다.
```r
# 5. Seurat v5 IntegrateLayers ------------------------------------------------

# RPCA integration은 sample 사이에서 공유되는 구조를 찾아 저차원 공간을 보정합니다.
# 원래 RNA count를 덮어쓰지 않고 integrated.rpca reduction을 새로 만듭니다.
intdata <- IntegrateLayers(
  object = filtdata,
  method = RPCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.rpca",
  assay = "SCT",
  normalization.method = "SCT",
  dims = dims_use,
  verbose = FALSE
)

cat("사용 가능한 dimensional reduction:\n")
print(Reductions(intdata))
```

완료 기준: `Reductions(intdata)`에 `integrated.rpca`가 나타납니다. 질병 차이를 무조건 제거하는 과정으로 해석하지 않습니다.

## 8. Clustering과 integration 후 UMAP

이웃 관계를 이용해 cluster를 만들고 UMAP에 표시합니다. Cluster 번호는 이름표이며 숫자의 크기에 생물학적 순서는 없습니다.
```r
# 6. Clustering과 UMAP --------------------------------------------------------

intdata <- FindNeighbors(
  object = intdata,
  reduction = "integrated.rpca",
  dims = dims_use,
  verbose = FALSE
)

intdata <- FindClusters(
  object = intdata,
  resolution = 0.2,
  cluster.name = "rpca_clusters",
  random.seed = 20260909,
  verbose = FALSE
)

intdata <- RunUMAP(
  object = intdata,
  reduction = "integrated.rpca",
  dims = dims_use,
  reduction.name = "umap.rpca",
  reduction.key = "UMAPRPCA_",
  seed.use = 20260909,
  verbose = FALSE
)

umap_sample <- DimPlot(
  intdata,
  reduction = "umap.rpca",
  group.by = "orig.ident"
) +
  ggtitle("RPCA integration 후: sample")

umap_cluster <- DimPlot(
  intdata,
  reduction = "umap.rpca",
  group.by = "rpca_clusters",
  label = TRUE,
  repel = TRUE
) +
  NoLegend() +
  ggtitle("RPCA integration 후: cluster")

print(umap_sample + umap_cluster)
ggsave(
  filename = file.path(result_dir, "04_UMAP_after_integration.png"),
  plot = umap_sample + umap_cluster,
  width = 14,
  height = 6,
  dpi = 200
)
```

완료 기준: sample 색과 cluster 색으로 그린 UMAP이 저장됩니다. 통합 전후를 비교하되, sample 혼합 정도와 다음 절의 marker 보존을 함께 확인합니다.

## 9. Marker로 cell type 후보 찾기

`integrated.rpca`는 이웃 탐색과 지도에, `RNA`의 normalized data는 marker 발현 확인에 사용합니다.

DotPlot의 점 크기는 발현 cell 비율입니다. 기본 설정의 색은 **gene별로 cluster 평균을 표준화한 값**이므로 서로 다른 gene의 절대 발현량을 비교하는 색이 아닙니다. `IL7R` 하나로 CD4 T를, `NKG7` 하나로 NK를 확정하지 말고 CD3 계열을 포함한 여러 marker를 함께 봅니다.

`RUN_FIND_ALL_MARKERS`는 기본 `FALSE`입니다. 후보 marker 계산도 해보고 싶으면 `TRUE`로 바꾸어 해당 블록을 실행합니다. 이 검정은 cluster 특징 탐색용이며 환자 간 질병 차이 검정을 대신하지 않습니다.
```r
# 7. Marker 확인과 cell type annotation ---------------------------------------

# RNA assay의 sample별 layer를 다시 합친 뒤 log-normalized data layer를 만듭니다.
# 이렇게 하면 marker 발현을 한 RNA layer에서 비교할 수 있습니다.
intdata[["RNA"]] <- JoinLayers(intdata[["RNA"]])
intdata <- NormalizeData(
  intdata,
  assay = "RNA",
  normalization.method = "LogNormalize",
  verbose = FALSE
)
DefaultAssay(intdata) <- "RNA"
Idents(intdata) <- "rpca_clusters"

marker_panels <- list(
  `T cell` = c("CD3D", "CD3E", "TRAC"),
  `CD4 T` = c("IL7R", "LTB", "CCR7"),
  `CD8 T / NK` = c("CD8A", "NKG7", "GNLY", "PRF1"),
  `B cell` = c("MS4A1", "CD79A", "CD37"),
  Monocyte = c("LYZ", "S100A8", "FCGR3A", "LILRB1"),
  mDC = c("CD1C", "CLEC10A", "FCER1A"),
  pDC = c("GZMB", "IL3RA", "CLEC4C"),
  Platelet = c("PPBP", "PF4"),
  Erythroid = c("HBA1", "HBB"),
  Proliferating = c("MKI67", "TOP2A")
)

# 데이터에 없는 gene은 자동으로 제외합니다.
marker_panels <- lapply(marker_panels, intersect, y = rownames(intdata))
marker_panels <- marker_panels[lengths(marker_panels) > 0]

marker_dotplot <- DotPlot(
  intdata,
  features = marker_panels,
  group.by = "rpca_clusters",
  assay = "RNA",
  dot.scale = 7
) +
  RotatedAxis() +
  labs(
    title = "Cluster annotation을 위한 marker 확인",
    subtitle = "점 크기: 발현 세포 비율, 색: gene별 표준화 평균 발현"
  )
print(marker_dotplot)
ggsave(
  filename = file.path(result_dir, "05_marker_DotPlot.png"),
  plot = marker_dotplot,
  width = 15,
  height = 8,
  dpi = 200
)

```

### 9-2. 선택: cluster별 후보 marker 계산하기

기본 DotPlot 확인만 진행하려면 아래 기본값 FALSE를 유지합니다.

```r
# 선택 실습: 각 cluster의 후보 marker를 계산합니다.
# 세포 수가 많으면 시간이 걸릴 수 있으므로 기본값은 FALSE입니다.
RUN_FIND_ALL_MARKERS <- FALSE

if (RUN_FIND_ALL_MARKERS) {
  cluster_markers <- FindAllMarkers(
    intdata,
    assay = "RNA",
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.25,
    verbose = FALSE
  )

  top_markers <- cluster_markers |>
    group_by(cluster) |>
    slice_max(order_by = avg_log2FC, n = 10, with_ties = FALSE) |>
    ungroup()

  write.csv(
    cluster_markers,
    file = file.path(result_dir, "cluster_markers_all.csv"),
    row.names = FALSE
  )
  write.csv(
    top_markers,
    file = file.path(result_dir, "cluster_markers_top10.csv"),
    row.names = FALSE
  )
}
```

완료 기준: `05_marker_DotPlot.png`가 저장됩니다. 각 cluster에 대해 후보 cell type과 근거 marker 두 개 이상을 적어보세요.

## 10. 직접 작성한 annotation 적용하기

**여기는 학생이 marker 근거로 이름표를 편집하는 단계입니다.** 모든 번호에 가짜 정답을 넣는 대신 `Unassigned`로 시작합니다. 불확실한 cluster는 그대로 남겨도 됩니다.

초기화 블록은 한 번만 실행합니다. 이름표를 수정한 후에는 `# 이름표를 모두 확인했으면` 아래부터 다시 실행하면 됩니다. 이후 annotation을 수정했다면 composition과 T cell 분석도 새 annotation으로 다시 계산해야 합니다.
```r
# 현재 cluster 번호를 확인합니다.
Idents(intdata) <- "rpca_clusters"
cluster_ids <- levels(Idents(intdata))
print(cluster_ids)

# 처음에는 모든 cluster를 미확정으로 둡니다.
celltype_map <- setNames(rep("Unassigned", length(cluster_ids)), cluster_ids)
print(celltype_map)

# 위 DotPlot을 보고 이 위치에 직접 작성합니다.
# 아래 줄은 형식 예시이며, 특정 cluster의 정답이 아닙니다.
# celltype_map["0"] <- "T cell"
# celltype_map["1"] <- "Monocyte"

```

초기화와 이름표 편집을 마쳤으면 다음 블록으로 annotation을 적용합니다.

```r
# 이름표를 모두 확인했으면 여기부터 실행합니다.
if (anyDuplicated(names(celltype_map)) ||
    !setequal(names(celltype_map), cluster_ids) ||
    anyNA(celltype_map) || any(!nzchar(trimws(celltype_map)))) {
  stop("celltype_map에 모든 cluster 번호와 유효한 이름이 있어야 합니다.")
}
intdata$celltype <- unname(celltype_map[as.character(intdata$rpca_clusters)])
if (all(intdata$celltype == "Unassigned")) {
  message("아직 모든 cell이 Unassigned입니다. marker를 확인하고 이름표를 수정하세요.")
}
print(table(intdata$rpca_clusters, intdata$celltype))

annotated_umap <- DimPlot(
  intdata, reduction = "umap.rpca", group.by = "celltype",
  label = TRUE, repel = TRUE
) + NoLegend() + ggtitle("Marker 근거로 작성한 cell type annotation")
print(annotated_umap)
ggsave(file.path(result_dir, "06_UMAP_annotated.png"),
       annotated_umap, width = 9, height = 7, dpi = 200)
write.csv(data.frame(cluster = names(celltype_map), celltype = unname(celltype_map)),
          file.path(result_dir, "cluster_annotation.csv"), row.names = FALSE)

# 기본 실습의 checkpoint: 심화 실습을 하지 않아도 저장합니다.
saveRDS(intdata, file.path(result_dir, "GSE244515_4sample_Seurat5.rds"),
        compress = FALSE)
writeLines(capture.output(sessionInfo()), file.path(result_dir, "sessionInfo.txt"))
```

완료 기준: cluster와 cell type의 대응표 및 annotated UMAP이 저장됩니다. T cell 심화는 T cell로 해석한 annotation이 하나 이상 있어야 진행할 수 있습니다.

## 11. T cell만 확대해서 다시 분석하기

전체 PBMC에서는 큰 계통 차이가 먼저 보입니다. T cell만 선택하여 variable gene, PCA, integration, clustering, UMAP을 다시 계산합니다. 앞 절의 annotation이 필요하며 이 절 전체는 추가 실습입니다. 생략할 경우 13절로 이동할 수 있습니다.

T cell subset에는 `k.weight = 50`을 사용합니다. Cell 수 조건을 통과해도 공유되는 population이 부족하면 anchor 오류가 날 수 있습니다. Sample별 cell 수와 marker 구조를 먼저 확인합니다.
```r
T_CELL_LABELS <- c("T cell", "CD4 T", "CD8 T", "Treg")
# 위 이름을 10절에서 실제로 사용한 T cell annotation과 맞춥니다.
available_tcell_labels <- intersect(T_CELL_LABELS, unique(intdata$celltype))
if (length(available_tcell_labels) == 0) {
  stop("T cell 이름표가 없습니다. 10절 annotation과 T_CELL_LABELS를 확인하세요.")
}
tcell_cells <- rownames(intdata[[]])[
  as.character(intdata$celltype) %in% available_tcell_labels
]

if (length(tcell_cells) < 100) {
  stop(
    "T cell로 선택된 cell이 100개 미만입니다. T_CELL_LABELS와 annotation을 확인하세요."
  )
}

tcell <- subset(intdata, cells = tcell_cells)
DefaultAssay(tcell) <- "RNA"

cat("T cell subset에 포함된 annotation:\n")
print(table(tcell$celltype))
cat("T cell subset의 sample별 cell 수:\n")
print(table(tcell$orig.ident))
if (length(unique(tcell$orig.ident)) < 2 || any(table(tcell$orig.ident) <= 31)) {
  stop("T cell integration에 사용할 sample 수와 sample별 cell 수를 확인하세요.")
}

# 전체 PBMC에서 사용한 PCA는 큰 cell type 차이를 잘 설명합니다.
# T cell 내부의 작은 차이를 보기 위해 T cell만으로 SCT, PCA와 clustering을 다시 수행합니다.
tcell[["RNA"]] <- JoinLayers(tcell[["RNA"]])
tcell[["RNA"]] <- split(tcell[["RNA"]], f = tcell$orig.ident)

```

### 11-2. T cell 내부의 발현 변이를 다시 계산하기

위에서 선택한 T cell로 SCT, PCA, integration과 지도를 다시 계산합니다.

```r
tcell <- SCTransform(
  object = tcell,
  assay = "RNA",
  new.assay.name = "SCT_T",
  vst.flavor = "v2",
  variable.features.n = 3000,
  conserve.memory = TRUE,
  verbose = FALSE
)

tcell <- RunPCA(
  object = tcell,
  assay = "SCT_T",
  npcs = 30,
  reduction.name = "pca.tcell",
  reduction.key = "PCATCELL_",
  verbose = FALSE
)

dims_tcell <- 1:20

tcell <- IntegrateLayers(
  object = tcell,
  method = RPCAIntegration,
  orig.reduction = "pca.tcell",
  new.reduction = "integrated.rpca.tcell",
  assay = "SCT_T",
  normalization.method = "SCT",
  dims = dims_tcell,
  k.weight = 50,
  verbose = FALSE
)

tcell <- FindNeighbors(
  object = tcell,
  reduction = "integrated.rpca.tcell",
  dims = dims_tcell,
  graph.name = c("tcell_nn", "tcell_snn"),
  verbose = FALSE
)

tcell <- FindClusters(
  object = tcell,
  graph.name = "tcell_snn",
  resolution = 0.4,
  cluster.name = "tcell_clusters",
  random.seed = 20260909,
  verbose = FALSE
)

tcell <- RunUMAP(
  object = tcell,
  reduction = "integrated.rpca.tcell",
  dims = dims_tcell,
  reduction.name = "umap.tcell",
  reduction.key = "UMAPTCELL_",
  seed.use = 20260909,
  verbose = FALSE
)

tcell_umap <- DimPlot(
  tcell,
  reduction = "umap.tcell",
  group.by = "tcell_clusters",
  label = TRUE,
  repel = TRUE
) +
  NoLegend() +
  ggtitle("T cell subset: 다시 계산한 cluster")
print(tcell_umap)

ggsave(
  filename = file.path(result_dir, "08_Tcell_UMAP.png"),
  plot = tcell_umap,
  width = 9,
  height = 7,
  dpi = 200
)

```

### 11-3. T cell의 type과 state marker 확인하기

새 cluster를 marker로 해석합니다. RNA normalized data에서 발현을 확인합니다.

```r
# Marker는 type과 state를 구분하여 해석합니다.
tcell_marker_panels <- list(
  `Naive / memory-like` = c("CCR7", "TCF7", "LEF1", "IL7R", "LTB"),
  Cytotoxicity = c("CD8A", "CCL5", "NKG7", "PRF1", "GZMB"),
  Treg = c("FOXP3", "IL2RA", "CTLA4"),
  Proliferating = c("MKI67", "TOP2A", "STMN1"),
  `IFN response` = c("ISG15", "IFIT1", "IFIT3", "MX1"),
  `Exhaustion-associated` = c("PDCD1", "LAG3", "HAVCR2", "TOX", "TIGIT")
)

tcell[["RNA"]] <- JoinLayers(tcell[["RNA"]])
tcell <- NormalizeData(tcell, assay = "RNA", verbose = FALSE)
DefaultAssay(tcell) <- "RNA"

tcell_marker_panels <- lapply(
  tcell_marker_panels,
  intersect,
  y = rownames(tcell)
)
tcell_marker_panels <- tcell_marker_panels[lengths(tcell_marker_panels) > 0]

tcell_marker_dotplot <- DotPlot(
  tcell,
  features = tcell_marker_panels,
  group.by = "tcell_clusters",
  assay = "RNA",
  dot.scale = 7
) +
  RotatedAxis() +
  labs(
    title = "T cell 안에서도 type과 state가 다릅니다",
    subtitle = "Exhaustion-associated marker 하나만으로 exhaustion을 확정하지 않습니다"
  )
print(tcell_marker_dotplot)

ggsave(
  filename = file.path(result_dir, "09_Tcell_marker_DotPlot.png"),
  plot = tcell_marker_dotplot,
  width = 15,
  height = 8,
  dpi = 200
)
```

완료 기준: `08_Tcell_UMAP.png`, `09_Tcell_marker_DotPlot.png`가 저장됩니다.

- Type/subtype 후보: Naive/memory-like, Treg 등.
- State/program 후보: Proliferation, IFN response, exhaustion-associated program.
- PDCD1 또는 TIGIT 하나만으로 exhaustion을 확정하지 않습니다. TIGIT는 Treg에서도 나타날 수 있습니다.

질문: Cytotoxicity가 높은 cluster는 CD8A도 높은가요? IFN response는 여러 subtype에 걸쳐 보이나요? Treg와 cytotoxic T cell의 기능은 암·자가면역·감염에서 어떻게 달라질까요?

## 12. Functional program을 module score로 비교하기

11절에서 만든 `tcell`이 필요합니다. 여러 gene의 발현을 control gene과 비교하여 상대 점수로 요약합니다. 데이터에 남은 gene 목록을 확인하고, 2개 미만인 program은 계산하지 않습니다.

서로 다른 program의 점수 크기를 직접 비교하여 어느 기능이 더 강하다고 결론내리지 않습니다. 각 program 내에서 cell 간 패턴을 보고 marker·sample 정보와 함께 해석합니다. 아래 그림도 program별 색 범위를 사용합니다.
```r
tcell_programs <- list(
  Cytotoxicity = c("NKG7", "CCL5", "PRF1", "GZMB"),
  IFN_response = c("ISG15", "IFIT1", "IFIT3", "MX1"),
  Exhaustion_associated = c("PDCD1", "LAG3", "HAVCR2", "TOX", "TIGIT")
)

tcell_programs <- lapply(tcell_programs, intersect, y = rownames(tcell))
tcell_programs <- tcell_programs[lengths(tcell_programs) >= 2]

print(tcell_programs)
if (length(tcell_programs) == 0) message("계산 가능한 program이 없습니다. gene 이름을 확인하세요.")
if (length(tcell_programs) > 0) {
  tcell <- AddModuleScore(
    object = tcell,
    features = tcell_programs,
    assay = "RNA",
    name = "Program",
    seed = 20260909
  )

  raw_score_names <- paste0("Program", seq_along(tcell_programs))
  score_names <- paste0(names(tcell_programs), "_score")

  for (i in seq_along(raw_score_names)) {
    tcell[[score_names[[i]]]] <- tcell[[raw_score_names[[i]]]][, 1]
  }

  program_score_plot <- FeaturePlot(
    tcell,
    features = score_names,
    reduction = "umap.tcell",
    ncol = length(score_names),
    keep.scale = "feature"
  )
  print(program_score_plot)

  ggsave(
    filename = file.path(result_dir, "10_Tcell_program_scores.png"),
    plot = program_score_plot,
    width = 5 * length(score_names),
    height = 5,
    dpi = 200
  )
}

saveRDS(tcell, file.path(result_dir, "GSE244515_Tcell_Seurat5.rds"), compress = FALSE)
```

완료 기준: 계산 가능한 program이 있으면 `10_Tcell_program_scores.png`가 저장됩니다. 점수는 절대 활성도나 임상 진단값이 아닙니다.

## 13. Cell composition은 sample별로 비교하기

T cell 심화 실습을 생략해도 10절의 `intdata`로 실행할 수 있습니다. 각 sample의 **QC 후 남은 PBMC 전체**를 분모로 cell type별 비율을 계산합니다. `Unassigned`도 표시하므로 미확정 cell을 숨겨 비율을 바꾸지 않습니다.
```r
# 10절 annotation이 필요합니다. 아직 미확정인 cell도 분모에 포함합니다.
# 존재하지 않는 sample-cell type 조합도 cell_count = 0으로 기록합니다.
composition_table <- as.data.frame(table(
  orig.ident = factor(intdata$orig.ident, levels = names(sample_dirs)),
  celltype = factor(intdata$celltype)
), responseName = "cell_count") |>
  mutate(condition = unname(condition_map[as.character(orig.ident)])) |>
  group_by(orig.ident) |>
  mutate(cell_proportion = cell_count / sum(cell_count)) |>
  ungroup()

print(composition_table)
write.csv(
  composition_table,
  file = file.path(result_dir, "celltype_composition_by_sample.csv"),
  row.names = FALSE
)

composition_plot <- ggplot(
  composition_table,
  aes(x = orig.ident, y = cell_proportion, fill = celltype)
) +
  geom_col(color = "white", linewidth = 0.15) +
  scale_y_continuous(labels = scales::percent) +
  theme_classic() +
  labs(
    title = "Cell composition은 sample별로 비교합니다",
    subtitle = "세포 수가 많아도 biological replicate는 sample입니다",
    x = "Sample",
    y = "Cell proportion",
    fill = "Cell type"
  )
print(composition_plot)

ggsave(
  filename = file.path(result_dir, "07_celltype_composition.png"),
  plot = composition_plot,
  width = 9,
  height = 6,
  dpi = 200
)
```

완료 기준: `07_celltype_composition.png`와 `celltype_composition_by_sample.csv`가 저장되고 sample별 비율 합은 1입니다.

PD1·PD2가 같은 방향을 보이는지, 한 sample이 결과를 주도하는지 비교하세요. Healthy 2명과 Periodontitis 2명의 탐색적 결과입니다. Cell이 수천 개라도 독립적인 환자가 수천 명인 것은 아닙니다. 채취·분리·QC에 따른 회수 편향도 composition에 영향을 줍니다.

## 14. 결과 저장과 재개

기본·추가 실습을 마친 뒤 실행합니다. T cell 분석을 생략해도 저장할 수 있습니다. RDS에는 계산한 object가 저장되고 `sessionInfo.txt`에는 실행 환경이 기록됩니다.
```r
# 11. 결과 저장 ---------------------------------------------------------------

saveRDS(
  intdata,
  file = file.path(result_dir, "GSE244515_4sample_Seurat5.rds"),
  compress = FALSE
)

if (exists("tcell") && inherits(tcell, "Seurat")) {
  saveRDS(
    tcell,
    file = file.path(result_dir, "GSE244515_Tcell_Seurat5.rds"),
    compress = FALSE
  )
}

writeLines(
  capture.output(sessionInfo()),
  con = file.path(result_dir, "sessionInfo.txt")
)
```

주요 결과:

| 단계 | 결과 파일 |
|---|---|
| QC | `01_QC_violin_before_filtering.png`, `02_QC_scatter_before_filtering.png`, `QC_cell_numbers.csv` |
| 통합 전후 | `03_UMAP_before_integration.png`, `04_UMAP_after_integration.png` |
| Marker / annotation | `05_marker_DotPlot.png`, `06_UMAP_annotated.png`, `cluster_annotation.csv` |
| Composition | `07_celltype_composition.png`, `celltype_composition_by_sample.csv` |
| T cell 추가 실습 | `08_Tcell_UMAP.png`, `09_Tcell_marker_DotPlot.png`, `10_Tcell_program_scores.png` |
| Object / 환경 | `GSE244515_4sample_Seurat5.rds`, `GSE244515_Tcell_Seurat5.rds`, `sessionInfo.txt` |

`FindAllMarkers`를 선택하면 `cluster_markers_all.csv`와 `cluster_markers_top10.csv`도 저장됩니다. 실행하지 않은 추가 단계의 결과는 생성되지 않습니다.

<details>
<summary>새 R session에서 저장한 PBMC object로 추가 실습 재개하기</summary>

2절 설정 블록을 먼저 실행한 뒤 다음 코드를 실행합니다. 이 블록은 처음 분석할 때는 생략합니다.

```r
intdata <- readRDS(file.path(result_dir, "GSE244515_4sample_Seurat5.rds"))
condition_map <- c(H1 = "Healthy", H2 = "Healthy",
                   PD1 = "Periodontitis", PD2 = "Periodontitis")
stopifnot("celltype" %in% colnames(intdata[[]]))
table(intdata$celltype)
```

11절 또는 13절로 이동합니다. 모두 `Unassigned`이면 먼저 9절 marker 확인과 10절 annotation을 진행합니다.
</details>

## 자주 생기는 문제

- **파일을 찾지 못함:** Project 폴더와 `data/GSE244515_4sample/H1` 등 네 경로, 세 파일 이름을 확인합니다. `.gz`는 풀지 않습니다.
- **object를 찾지 못함:** 해당 object를 만드는 앞 절을 실행했는지 확인합니다. R session을 재시작하면 메모리의 object는 사라집니다.
- **Seurat가 없음 / 버전 불일치:** 1절 설치 블록을 현재 R에서 실행하고 session을 재시작합니다. 목표 5.x가 아닌 버전이 설치되면 수업용 버전을 강사와 확인합니다.
- **data layers are not joined:** RNA marker 분석 전에 9절의 `JoinLayers()`와 `NormalizeData()`를 실행합니다. SCT assay에 그대로 적용하지 않습니다.
- **메모리 부족:** 새 session에서 필요 단계만 실행하고 `conserve.memory = TRUE`를 유지합니다. 계산을 나누거나 강사가 검증한 checkpoint를 사용합니다. `future.globals.maxSize`를 높이는 것만으로 RAM 부족이 해결되지는 않습니다.
- **Number of anchor cells is less than k.weight:** sample별 cell 수·공통 population을 확인합니다. 필요 시 해당 `IntegrateLayers()`의 `k.weight`를 확보된 anchor 수보다 작게 조정하되, 통합 후 marker 구조를 다시 점검합니다.
- **T cell 이름표가 없음:** 10절의 실제 cell type 이름과 11절 `T_CELL_LABELS`를 맞춥니다. NK를 cytotoxic marker만 보고 T cell에 포함하지 않습니다.
- **그림의 cluster 번호가 예시와 다름:** 버전·seed·설정에 따라 달라질 수 있습니다. 현재 결과의 marker로 해석합니다.

## 해석할 때 기억할 점

- Cluster는 계산 결과이고 cell type은 marker와 생물학을 이용한 해석입니다.
- UMAP의 거리·섬 크기만으로 기능 차이나 조직 내 위치를 판단하지 않습니다.
- Cell type과 cell state를 나누어 해석하고 한 marker로 기능을 확정하지 않습니다.
- 질환군 비교의 biological replicate는 환자/sample입니다.

## 참고 자료

- [Seurat v5 Essential Commands](https://satijalab.org/seurat/articles/seurat5_essential_commands.html)
- [Seurat v5 Integrative Analysis](https://satijalab.org/seurat/articles/seurat5_integration.html)
- [Using sctransform in Seurat](https://satijalab.org/seurat/articles/sctransform_vignette)
- [Seurat AddModuleScore](https://satijalab.org/seurat/reference/addmodulescore)
- [NCBI GEO GSE244515](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE244515)
- Lee H, Joo J, Song J, et al. *Immunological link between periodontitis and type 2 diabetes deciphered by single-cell RNA analysis*. Clinical and Translational Medicine. 2023. [doi:10.1002/ctm2.1503](https://doi.org/10.1002/ctm2.1503)


