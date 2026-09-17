# scRNA-seq analysis workshop with Seurat v5

## Preparation: Healthy와 severe COVID-19 PBMC 비교

**GSE149689의 Healthy 1개와 severe COVID-19 1개**, 총 두 사람의 PBMC를 사용합니다. 두 sample 모두 **10x Chromium Single Cell 3′ v3의 droplet 방식**으로 만들어졌습니다. 

**GSE149689 링크**: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149689

| 실습 이름 | GEO sample | 원래 sample 이름 | 조건 | 나이 / 성별 | 원본 barcode 접미사 |
|---|---|---|---|---|---|
| H1 | [GSM4509015](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4509015) | Sample 5 / Normal 1 | Healthy | 63세 / 여성 | `-5` |
| COVID1 | [GSM4509011](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4509011) | Sample 1 / nCoV 1 | Severe_COVID | 63세 / 남성 | `-1` |
<br>
<br>
### 실습용 파일 준비

아래 두 파일을 내려받아 폴더에 둡니다. GitHub 파일 화면에서는 **Download raw file**을 선택합니다. 학생은 원본 전체 matrix를 받을 필요가 없습니다. 

- [H1_counts.rds — Healthy, 약 10.5 MB](data/GSE149689_2sample/H1_counts.rds)
- [COVID1_counts.rds — Severe COVID-19, 약 9.3 MB](data/GSE149689_2sample/COVID1_counts.rds)

이 링크는 README와 함께 `data/GSE149689_2sample/` 폴더를 같은 저장소에 올렸을 때 작동합니다.

```text
실습폴더/
└─ data/
   └─ GSE149689_2sample/
      ├─ H1_counts.rds
      └─ COVID1_counts.rds
```

두 RDS는 **원본 count matrix에서 해당 sample의 열만 추출한 분석 시작 자료**입니다. 

## 실습 목표와 실행 방법:
**QC → SCTransform → integration → clustering → UMAP → annotation**을 진행합니다. 이어서 T cell, module score, cell composition을 살펴봅니다.

- 별도 R 파일 없이 이 문서의 코드 블록을 **하나씩 순서대로** 실행합니다.
- 기본 실습은 1–10절, 추가 실습은 11–13절입니다. 마지막에 14절에서 저장합니다.
- `선택`으로 표시한 접힌 부분은 필요할 때만 실행합니다.
- 오류가 나면 다음 블록으로 넘어가지 말고 강사와 함께 확인합니다.
- `<-`는 결과에 이름을 붙이는 기호입니다. 

<br>
<br>

## 1. R / RStudio와 패키지 준비

강사가 준비한 R / Seurat v5 환경을 우선 사용합니다. 설치가 안 된 컴퓨터에서만 아래 코드를 한 번 실행합니다. 필요한 의존성 패키지도 함께 설치됩니다.

```r
install.packages(c("Seurat", "ggplot2", "dplyr", "patchwork"))
```

설치 후 **Session → Restart R**를 선택하고 버전을 확인합니다.

```r
R.version.string
packageVersion("Seurat")
```

<br>
<br>

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

```r
getwd()
setwd("data 폴더를 다온로드 받은 디렉토리 설정")  
data_dir <- "data/GSE149689_2sample"
dir.create("results", showWarnings = FALSE)
list.files(data_dir)
```

**확인:** 위에서 준비한 `.rds` 파일 두 개가 보이나요? 다른 곳에 저장했다면 `data_dir`만 실제 경로로 바꿉니다. Windows 경로는 `C:/Users/...`처럼 `/`를 사용합니다.

<br>
<br>

## 3. 두 sample에서 Seurat object 만들기

### 3-1. Count matrix 읽기

`readRDS()`로 저장된 matrix를 읽습니다. 행은 gene, 열은 cell barcode, 숫자는 UMI count입니다.

```r
H1_counts <- readRDS(file.path(data_dir, "H1_counts.rds"))
COVID1_counts <- readRDS(file.path(data_dir, "COVID1_counts.rds"))
dim(H1_counts)
dim(COVID1_counts)
H1_counts[1:5, 1:5]
```

`[1:5, 1:5]`는 첫 다섯 gene과 첫 다섯 cell을 선택합니다.

Sparse matrix의 `.`은 0을 뜻합니다. 0은 이번 측정에서 검출되지 않았다는 뜻이며, 그 세포에 해당 RNA가 절대 없다는 뜻은 아닙니다.

### 3-2. Sample별 object 만들기

`min.cells = 3`은 해당 sample에서 최소 세 cell에 검출된 gene을 남깁니다.

```r
H1 <- CreateSeuratObject(H1_counts, project = "H1", min.cells = 3)
COVID1 <- CreateSeuratObject(COVID1_counts, project = "COVID1", min.cells = 3)
```

### 3-3. 합치고 조건 정보 붙이기

```r
rawdata <- merge(H1, y = COVID1, add.cell.ids = c("H1", "COVID1"))
head(rawdata@meta.data)

condition_map <- c(H1 = "Healthy", COVID1 = "Severe_COVID")
rawdata$condition <- unname(condition_map[as.character(rawdata$orig.ident)])
head(rawdata@meta.data)

table(rawdata$orig.ident, rawdata$condition)
Layers(rawdata[["RNA"]])
```

`orig.ident`는 sample 이름, `condition`은 비교 조건입니다. 

H1은 Healthy, COVID1은 Severe_COVID에만 속해야 합니다. 

RNA의 sample별 layer는 서로 다른 sample의 측정값을 보관하는 칸입니다.

<br>
<br>

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

아래는 **분포 확인 후 조정할 교육용 QC 시작 기준**입니다. GSE149689에 최적화되거나 원 논문을 재현하는 기준은 아닙니다. 강사는 두 sample의 분포와 남는 cell 수를 확인한 후 수업용 값을 정합니다. 모든 조직에 그대로 적용하지 않습니다. UMI가 높다는 이유만으로 doublet을 확정할 수 없고, mitochondrial 비율도 cell type에 따라 달라집니다.

필터링 후 Sample별 남은 cell 수를 비교합니다.

```r
filtdata <- subset(
  rawdata,
  subset = nFeature_RNA > 600 & nFeature_RNA < 5000 &
    nCount_RNA < 25000 & percent.mt < 20
)

qc_summary <- data.frame(
  before = table(factor(rawdata$orig.ident, levels = names(condition_map))),
  after = as.vector(table(factor(filtdata$orig.ident, levels = names(condition_map))))
)
colnames(qc_summary) <- c("sample", "before", "after")
qc_summary
```

<br>
<br>

## 5. SCTransform과 PCA

(지금 Sample 별로 layer가 분리되어 있으므로 생략 가능) Sample별 layer를 준비합니다.

```r
filtdata[["RNA"]] <- JoinLayers(filtdata[["RNA"]])
filtdata[["RNA"]] <- split(filtdata[["RNA"]], f = filtdata$orig.ident)
Layers(filtdata[["RNA"]])
```

SCTransform으로 기술적 변이를 모델링하고 PCA로 주요 발현 차이를 요약합니다. 

```r
BiocManager::install('glmGamPoi')   # SCTransform 계산 속도를 향상

filtdata <- SCTransform(filtdata, conserve.memory = TRUE, seed.use = 12345)
filtdata <- RunPCA(filtdata, seed.use = 12345)
ElbowPlot(filtdata, ndims = 50)
```

**확인:** ElbowPlot이 보이나요? 이번 실습은 이후 계산에 PC 1–30을 사용합니다.


### 중간 결과 저장

계산 시간은 컴퓨터와 cell 수에 따라 달라집니다. 아래는 **지금 분석한 COVID 데이터**를 저장하는 코드입니다. 이전 데이터의 checkpoint는 사용하지 않습니다.

```r
saveRDS(filtdata, "results/GSE149689_SCT_PCA.rds")
```

<details>
<summary>선택: 나중에 새 session에서 이 단계부터 재개하기</summary>

2절의 패키지를 먼저 불러온 뒤 실행합니다. 방금 저장한 경우 바로 실행할 필요는 없습니다.

```r
filtdata <- readRDS("results/GSE149689_SCT_PCA.rds")
condition_map <- c(H1 = "Healthy", COVID1 = "Severe_COVID")
```

이후 6절부터 진행합니다. `rawdata`와 QC 그림은 이 checkpoint에 들어 있지 않습니다.
</details>

## 6. Integration 전 UMAP

```r
filtdata <- RunUMAP(
  filtdata, dims = 1:30,
  reduction.name = "umap.unintegrated", seed.use = 12345
)
umap_before <- DimPlot(filtdata, reduction = "umap.unintegrated", group.by = "orig.ident")
umap_before
```

점 하나는 cell이고 색은 sample입니다. Sample별 분리는 기술적·생물학적 차이 모두에서 생길 수 있습니다. UMAP 축은 실제 조직 좌표가 아닙니다.

<br>
<br>

## 7. RPCA integration

Sample 사이에 공유되는 구조를 정렬합니다. 모든 sample이 완전히 섞이는 것이 목표는 아닙니다. 질환에 특이적인 cell state까지 지워지지 않았는지 integration 전후와 marker를 함께 봅니다. `SCT`를 사용했다는 설정과 결과 이름은 분석에 필요하므로 남깁니다.

```r
intdata <- IntegrateLayers(
  filtdata, method = RPCAIntegration,
  assay = "SCT", normalization.method = "SCT",
  new.reduction = "integrated.rpca", dims = 1:30
)
Reductions(intdata)
```

**확인:** `integrated.rpca`가 있나요? 원래 RNA count를 지우는 과정은 아닙니다.

<br>
<br>

## 8. Clustering과 integration 후 UMAP

### 8-1. 비슷한 이웃끼리 묶고 지도그리기

`resolution`은 cluster를 나누는 세밀함에 영향을 줍니다. 이번에는 `0.2`를 사용합니다.

```r
intdata <- FindNeighbors(intdata, reduction = "integrated.rpca", dims = 1:30)
intdata <- FindClusters(intdata, resolution = 0.2,
                        cluster.name = "rpca_clusters", random.seed = 12345)
intdata <- RunUMAP(
  intdata, reduction = "integrated.rpca", dims = 1:30,
  reduction.name = "umap.rpca", seed.use = 12345
)

umap_after <- DimPlot(intdata, reduction = "umap.rpca", group.by = "orig.ident") +
  DimPlot(intdata, reduction = "umap.rpca", group.by = "rpca_clusters", label = TRUE)
umap_after
DimPlot(intdata, reduction = "umap.rpca", group.by = "condition")
```

**질문:** Sample이 섞였나요? 다음 절에서 marker도 유지되는지 확인합니다. Cluster 번호는 계산 결과의 이름표입니다.

<br>
<br>

## 9. Marker로 cell type 후보 찾기

Marker는 원래 RNA의 normalized expression으로 확인합니다. Sample별 RNA layer를 합치고, cluster 번호를 현재 이름표로 지정합니다. RPCA 좌표는 cell을 묶는 데 사용하고 발현량 자체로 사용하지 않습니다.

```r
intdata[["RNA"]] <- JoinLayers(intdata[["RNA"]])
DefaultAssay(intdata) <- "RNA"
intdata <- NormalizeData(intdata)
Idents(intdata) <- "rpca_clusters"
```

```r
marker_panels <- list(
  `T cell` = c("CD3D", "CD3E", "TRAC", "TRDC"),
  `NK` = c("NKG7", "GNLY", "PRF1"),
  `B cell` = c("MS4A1", "CD79A", "CD37"),
  Monocyte = c("LYZ", "S100A8", "CD14", "FCGR3A"),
  Plasmablast = c("MZB1", "JCHAIN", "IGJ", "XBP1"),
  `Neutrophil-like` = c("LTF", "ELANE", "MPO"),
  mDC = c("CD1C", "CLEC10A", "FCER1A"),
  pDC = c("GZMB", "IL3RA", "CLEC4C"),
  Platelet = c("PPBP", "PF4"),
  Proliferating = c("MKI67", "TOP2A")
)
marker_panels <- lapply(marker_panels, intersect, y = rownames(intdata))
lengths(marker_panels)
marker_dotplot <- DotPlot(intdata, features = marker_panels) + RotatedAxis()
marker_dotplot
```

각 panel에 marker가 남아 있는지 먼저 확인합니다. `lapply`는 목록의 각 항목에 같은 작업을 적용하며, `intersect`는 실제 데이터에 있는 gene만 남깁니다. JCHAIN과 이전 symbol IGJ 중 실제 matrix에 있는 이름을 사용합니다. Neutrophil-like는 검토할 후보 이름이며 lineage 전환을 뜻하지 않습니다.

점 크기는 발현 cell 비율, 색은 gene별로 표준화한 cluster 평균입니다. 다른 gene끼리 색만 보고 절대 발현량을 비교하지 않습니다. Marker가 없다는 경고가 나오면 gene 이름을 확인합니다.

**활동:** 각 cluster의 후보 이름과 근거 marker 두 개를 적어봅니다. IL7R 하나로 CD4 T를, NKG7 하나로 NK를 확정하지 않습니다.

TRAC은 αβ T cell, TRDC는 γδ T cell을 구분하는 데 도움이 됩니다. Cytotoxic T cell과 NK는 NKG7·GNLY 등을 공유하므로 CD3D·CD3E와 함께 확인합니다.

<details>
<summary>선택: cluster별 후보 marker 계산하기</summary>

클러스터 별로 높게 발현되는 유전자를 찾습니다.

```r
cluster_markers <- FindAllMarkers(intdata, only.pos = TRUE,
                                 min.pct = 0.25, logfc.threshold = 0.25,
                                 random.seed = 12345)
top_markers <- cluster_markers |>
  group_by(cluster) |>
  slice_max(avg_log2FC, n = 10, with_ties = FALSE)
top_markers
```
</details>

<br>
<br>

## 10. 직접 annotation 붙이기

**이전 데이터의 cluster 번호를 재사용하지 않습니다.** 번호는 실행 환경과 분석 설정에 따라 달라질 수 있습니다. 9절 DotPlot에서 여러 marker를 확인한 뒤 이름을 붙입니다.

먼저 모든 cluster를 `Unassigned`로 둡니다. 초기화 블록은 한 번만 실행합니다.

```r
cluster_ids <- levels(Idents(intdata))
celltype_map <- setNames(rep("Unassigned", length(cluster_ids)), cluster_ids)
celltype_map
```

아래는 **입력 형식만 보여주는 예시**입니다. 관찰한 cluster 번호와 이름으로 수정하고, 실행할 줄의 맨 앞 `#`를 지웁니다. 여러 cluster가 같은 cell type이어도 됩니다.

```r
# celltype_map["0"] <- "T cell"
# celltype_map["1"] <- "Monocyte"
# celltype_map["2"] <- "B cell"
# celltype_map["3"] <- "Plasmablast"
```

가능한 이름: `T cell`, `NK`, `B cell`, `Monocyte`, `Plasmablast`, `Neutrophil-like`, `mDC`, `pDC`, `Platelet`. Proliferating은 별도의 lineage가 아니라 state일 수 있습니다. 근거가 부족하면 `Unassigned`로 남깁니다.

<details>
<summary>검증 실행의 참고 annotation — 직접 marker를 확인한 뒤 펼치세요</summary>

전체 QC 통과 cell 7,931개, resolution 0.2, 아래에 명시한 검증 환경에서 나온 **넓은 cell type 수준의 참고 분류**입니다. 다운샘플링했거나 cluster 수·marker가 다르면 이 번호를 그대로 적용하지 않습니다. 각 cluster 내부에 더 작은 subtype이나 혼합 population이 남아 있을 수 있습니다.

| Cluster | 참고 이름 | 확인한 근거 |
|---|---|---|
| 0 | Monocyte | LYZ, CD14, S100A8/S100A9 |
| 1 | T cell | CD3D/CD3E, TRAC, CD8A/CD8B |
| 2 | T cell | CD3D/CD3E, TRAC |
| 3 | B cell | MS4A1, CD79A |
| 4 | T cell | CD3D/CD3E, TRDC; γδ T-like |
| 5 | B cell | MS4A1, CD79A |
| 6 | NK | GNLY, NKG7, PRF1; 상대적으로 낮은 CD3D/TRAC |
| 7 | Monocyte | LYZ, FCGR3A; CD16 Monocyte-like |
| 8 | Platelet | PPBP, PF4 |
| 9 | mDC | CD1C, FCER1A, CLEC10A |

```r
celltype_map <- c(
  "0" = "Monocyte", "1" = "T cell", "2" = "T cell",
  "3" = "B cell", "4" = "T cell", "5" = "B cell",
  "6" = "NK", "7" = "Monocyte", "8" = "Platelet", "9" = "mDC"
)
```

Marker 목록에 있다고 모든 cell type이 별도 cluster로 검출되는 것은 아닙니다. 이번 설정에서는 Plasmablast나 pDC라는 이름을 억지로 붙이지 않았습니다.
</details>

작성한 이름을 cell별로 붙입니다.

```r
intdata$celltype <- unname(celltype_map[as.character(intdata$rpca_clusters)])
table(intdata$rpca_clusters, intdata$celltype)
annotated_umap <- DimPlot(intdata, reduction = "umap.rpca",
                          group.by = "celltype", label = TRUE)
annotated_umap
```

**완료 기준:** 주요 cluster에 근거 marker 두 개 이상을 설명할 수 있나요? 전부 `Unassigned`라면 annotation을 마친 후 다음으로 넘어갑니다. T cell을 하나 이상 확인해야 12–13절을 진행할 수 있습니다.

<br>
<br>

## 11. 추가 실습: Sample별 cell composition

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
  facet_grid(~condition, scales = "free_x", space = "free_x") +
  theme_classic() +
  labs(x = "Sample", y = "Cell proportion")
composition_plot
```

**질문:** H1과 COVID1에서 어떤 cell type의 비율이 다르게 보이나요? Healthy 한 명과 Severe_COVID 한 명의 탐색적 결과입니다. Cell 수천 개가 독립 환자 수천 명을 뜻하지 않습니다. 채취·분리·QC에 따른 편향도 고려합니다.

참고 annotation을 사용한 검증 결과에서 Monocyte 비율은 **H1 12.1%, COVID1 44.4%**, T cell 비율은 **H1 59.9%, COVID1 25.5%**였습니다. 이는 QC를 통과한 포획 cell 안에서의 비율이며 혈액 내 절대 세포 수를 뜻하지 않습니다.


---

### 11-1. 같은 Monocyte 안에서 cell state 비교하기

전체 PBMC 평균의 차이가 세포 구성 때문인지, 같은 종류의 세포 내부 차이인지 구분해 봅니다. 아래 코드는 10절에서 `Monocyte`라고 붙인 cell을 사용합니다. 두 sample 모두 Monocyte가 있는지 먼저 확인합니다.

```r
mono <- subset(intdata, subset = celltype == "Monocyte")
table(mono$orig.ident)
mono_genes <- c("CD14", "FCGR3A", "ISG15", "IFIT1", "IFIT3", "MX1",
                "S100A8", "S100A9", "IL1B", "NFKBIA")
mono_genes <- intersect(mono_genes, rownames(mono))
mono_dotplot <- DotPlot(mono, features = mono_genes,
                        group.by = "orig.ident", scale = FALSE) + RotatedAxis()
mono_dotplot
```

`scale = FALSE`로 두 sample 사이의 작은 차이를 gene별 표준화 색상으로 과장하지 않도록 합니다. 같은 gene의 두 sample을 비교하세요. Gene마다 원래 발현량이 달라 서로 다른 gene의 색을 직접 비교하지 않습니다.

**질문:** IFN-response gene과 염증 관련 gene은 두 sample에서 어떤 패턴인가요? 원 연구는 중증 COVID-19의 classical Monocyte에서 IFN response와 염증 반응을 보고했습니다. 아래 gene 목록은 이를 탐색하기 위한 짧은 교육용 예시이며, 해당 경로의 완전한 정의는 아닙니다.

검증 결과에서 COVID1에서 높게 관찰한 세 gene을 분포로 확인합니다.

```r
mono_inflammation_plot <- VlnPlot(
  mono, features = c("S100A8", "S100A9", "IL1B"),
  group.by = "orig.ident", pt.size = 0, ncol = 3
)
mono_inflammation_plot
```

이 세 gene의 발현만으로 모든 염증 경로나 실제 cytokine 분비량을 확정하지 않습니다.

여러 IFN-response gene을 함께 요약해 봅니다.

```r
mono_ifn <- list(c("ISG15", "IFIT1", "IFIT3", "MX1"))
mono_ifn <- lapply(mono_ifn, intersect, y = rownames(mono))
lengths(mono_ifn)
```

두 개 이상의 gene이 남았는지 확인한 뒤 실행합니다.

```r
mono <- AddModuleScore(mono, features = mono_ifn, name = "IFN", seed = 12345)
VlnPlot(mono, features = "IFN1", group.by = "orig.ident", pt.size = 0)
```

IFN1은 control gene과 비교한 상대 점수입니다. Monocyte 안에서도 CD14/CD16 subtype 구성에 따라 값이 달라질 수 있습니다. 두 사람만으로 질환 효과를 확정하지 않으며, 이번 실습에서는 질환군 간 p-value를 계산하지 않습니다. 관찰한 방향과 차이의 크기를 실제 결과로 설명합니다.

검증 실행에서 평균 IFN1은 **H1 약 0.100, COVID1 약 −0.047**이었습니다. 음수는 IFN RNA가 없다는 뜻이 아니라 control gene과 비교한 상대값입니다. 이번 비교에서는 염증 관련 gene과 IFN score가 같은 방향으로 변하지 않았습니다. 원 논문 전체의 경향을 이 한 쌍의 정답으로 강요하지 않습니다.

## 12. 추가 실습: T cell만 확대하기

### 12-1. T cell 선택

`T_CELL_LABELS`를 10절에서 사용한 실제 이름과 맞춥니다. 이름이 없으면 10절로 돌아갑니다. 이 절을 생략할 때는 13절도 건너뛰고 14절로 이동합니다.

```r
T_CELL_LABELS <- c("T cell", "CD4 T", "CD8 T", "Treg")
table(intdata$celltype)
tcell <- subset(intdata, subset = celltype %in% T_CELL_LABELS)
table(tcell$orig.ident)
```

Sample이 하나뿐이거나 cell이 매우 적으면 강사와 확인합니다. 전체 PBMC의 PCA를 확대하는 대신 T cell 내부의 차이를 다시 계산합니다.

아래의 `tcell[["SCT"]] <- NULL`은 **T cell object 안에 복사된 이전 SCT 결과만** 지웁니다. 원래 `intdata`는 유지됩니다. 이후 T cell의 RNA로 SCT를 새로 계산합니다.

```r
tcell[["RNA"]] <- JoinLayers(tcell[["RNA"]])
tcell[["RNA"]] <- split(tcell[["RNA"]], f = tcell$orig.ident)
DefaultAssay(tcell) <- "RNA"
tcell[["SCT"]] <- NULL
tcell <- SCTransform(tcell, conserve.memory = TRUE, seed.use = 12345)
tcell <- RunPCA(tcell, npcs = 30, reduction.name = "pca.tcell", seed.use = 12345)
```

### 12-2. 다시 통합하고 묶기

T cell에서는 PC 1–20, resolution 0.4를 사용합니다. 작은 subset을 고려해 `k.weight = 50`을 유지합니다.

```r
tcell <- IntegrateLayers(
  tcell, method = RPCAIntegration, assay = "SCT",
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

### 12-3. Type과 state 구분하기

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

## 13. 추가 실습: Module score

12절에서 만든 `tcell`이 필요합니다. 여러 gene의 발현 경향을 상대 점수로 요약합니다.

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

IFN response와 Exhaustion-associated score는 탐색용입니다. Monocyte에서 본 변화가 T cell에서도 똑같이 나타난다고 가정하지 않습니다. 높은 Exhaustion-associated score만으로 기능적 exhaustion을 확정할 수 없습니다.

점수는 control gene과 비교한 상대값입니다. 그림마다 색 범위가 다를 수 있으며, 다른 program끼리 점수 크기를 직접 비교하지 않습니다. 각 program 내 패턴을 marker·sample 정보와 함께 해석합니다.

같은 cell type 안에서도 **환자별**로 IFN response를 비교합니다.

```r
VlnPlot(tcell, features = "IFN_response_score", group.by = "orig.ident", pt.size = 0)
```

T cell subtype 비율이 다르면 전체 T cell 점수도 달라질 수 있습니다. 조건별 차이가 보이면 같은 subtype 안에서도 확인해야 합니다.

T cell 추가 실습을 마쳤다면 저장합니다.

```r
saveRDS(tcell, "results/GSE149689_Tcell_Seurat5.rds")
```


## 14. 결과 저장과 재개

기본 실습 결과와 실행 환경을 저장합니다. 같은 파일 이름으로 저장하면 이전 결과를 갱신합니다.

```r
saveRDS(intdata, "results/GSE149689_2sample_Seurat5.rds")
writeLines(capture.output(sessionInfo()), "results/sessionInfo.txt")
```

그림은 RStudio의 **Plots → Export**로 저장할 수 있습니다. 아래 자동 저장은 선택입니다.

<details>
<summary>선택: 기본 실습 그림과 표 저장</summary>

1–10절을 처음부터 마친 뒤 실행합니다. Checkpoint에서 시작했다면 만들지 않은 QC 그림과 `qc_summary`의 저장 줄은 건너뜁니다. Width와 height는 inch 단위의 그림 크기입니다.

```r
ggsave("results/01_QC_violin_before_filtering.png", qc_violin, width = 14, height = 5, bg = "white")
ggsave("results/03_UMAP_before_integration.png", umap_before, width = 8, height = 6, bg = "white")
ggsave("results/04_UMAP_after_integration.png", umap_after, width = 14, height = 6, bg = "white")
ggsave("results/05_marker_DotPlot.png", marker_dotplot, width = 15, height = 8, bg = "white")
ggsave("results/06_UMAP_annotated.png", annotated_umap, width = 9, height = 7, bg = "white")
```

```r
write.csv(qc_summary, "results/QC_cell_numbers.csv", row.names = FALSE)
annotation_table <- data.frame(cluster = names(celltype_map), celltype = unname(celltype_map))
write.csv(annotation_table, "results/cluster_annotation.csv", row.names = FALSE)
```

</details>

<details>
<summary>선택: T cell 그림 저장 — 12–13절 완료 후</summary>

```r
ggsave("results/08_Tcell_UMAP.png", tcell_umap, width = 9, height = 7, bg = "white")
ggsave("results/09_Tcell_marker_DotPlot.png", tcell_dotplot, width = 15, height = 8, bg = "white")
ggsave("results/10_Tcell_program_scores.png", program_plot, width = 15, height = 5, bg = "white")
```

</details>

<details>
<summary>선택: Composition 저장 — 11절 완료 후</summary>

```r
ggsave("results/07_celltype_composition.png", composition_plot, width = 9, height = 6, bg = "white")
write.csv(composition_table, "results/celltype_composition_by_sample.csv", row.names = FALSE)
```

</details>

<details>
<summary>선택: 새 R session에서 저장한 object 불러오기</summary>

같은 Project를 열고 2절의 패키지를 불러온 뒤 실행합니다. 처음 분석할 때는 생략합니다.

```r
intdata <- readRDS("results/GSE149689_2sample_Seurat5.rds")
condition_map <- c(H1 = "Healthy", COVID1 = "Severe_COVID")
table(intdata$celltype)
```

</details>

## 막혔을 때 확인하기

- **파일을 못 찾음:** Project 위치와 `data_dir`, 두 학생용 `.rds` 파일명을 확인합니다.
- **object가 없음:** 앞 블록을 실행했는지 확인합니다. Session을 재시작하면 메모리의 object가 사라집니다.
- **MT- gene이 0개:** 사람 gene symbol 대신 다른 ID를 읽었는지 확인합니다.
- **입력 파일 혼동:** 학생용 RDS는 count matrix입니다. Seq-Well의 exon/intron list나 이전 데이터의 분석 완료 object를 사용하지 않습니다.
- **원본 sample 선택 오류:** barcode 접미사 `-1$`과 `-5$`로 선택합니다. `$`가 빠지면 `-1`이 `-10` 등과 함께 선택될 수 있습니다.
- **메모리 부족:** 새 session에서 필요한 단계만 실행하거나 검증된 checkpoint를 사용합니다. `conserve.memory = TRUE`를 유지합니다.
- **data layers are not joined:** 9절의 RNA `JoinLayers()`를 실행합니다. SCT assay에는 그대로 적용하지 않습니다.
- **PCA에서 일부 feature가 scaled되지 않았다는 경고:** sample별 SCT 처리에서 일부 gene이 공통 PCA 입력에 남지 않을 수 있습니다. 검증 실행에서도 발생했으나 분석은 완료되었습니다. Marker는 별도로 RNA assay에서 확인합니다.
- **PC_ key가 이미 있다는 경고:** T cell PCA의 내부 이름이 자동으로 바뀌는 안내입니다. 코드에서 지정한 `pca.tcell` reduction 이름은 유지됩니다.
- **anchor / k.weight 오류:** Sample별 남은 cell 수와 공유 population을 먼저 확인합니다. 강사와 `k.weight`를 조정한 뒤 marker 구조도 다시 확인합니다.
- **T cell을 선택하지 못함:** 실제 annotation 이름과 `T_CELL_LABELS`를 맞춥니다. NK를 cytotoxic marker 하나로 T cell에 포함하지 않습니다.
- **module score gene이 부족함:** `lengths(tcell_programs)`와 gene 이름을 확인합니다.

## 참고 자료

- [Seurat v5 integration](https://satijalab.org/seurat/articles/seurat5_integration.html)
- [SCTransform 기본값](https://satijalab.org/seurat/reference/sctransform)
- [RunPCA 기본값과 seed](https://satijalab.org/seurat/reference/runpca)
- [RPCAIntegration](https://satijalab.org/seurat/reference/rpcaintegration)
- [AddModuleScore](https://satijalab.org/seurat/reference/addmodulescore)
- [GSE149689](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149689)
- 원 연구와 샘플 메타데이터: [GSE149689](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149689), [Healthy sample](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4509015), [Severe COVID-19 sample](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4509011).

