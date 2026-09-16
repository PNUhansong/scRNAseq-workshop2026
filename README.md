
# scRNA-seq analysis workshop with Seurat v5

## Preparation: GSE150728 다운로드

이번 실습은 **Healthy 2명과 severe COVID-19 2명**, 총 네 명의 PBMC를 비교합니다.
환자는 논문 Figure 1a에서 채혈 당시 **ARDS와 mechanical ventilation**이 확인되는 서로 다른 환자 C3·C4를 선택했습니다.

| 실습 sample 이름 | 논문 donor | GEO sample | GEO의 sample 이름 | 비교 조건 |
|---|---|---|---|---|
| H1 | H1 | [GSM4557334](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4557334) | HIP002 | Healthy |
| H2 | H2 | [GSM4557335](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4557335) | HIP015 | Healthy |
| C3 | C3 | [GSM4557330](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4557330) | covid_557 | Severe_COVID |
| C4 | C4 | [GSM4557331](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4557331) | covid_558 | Severe_COVID |

선정 근거: [GEO의 sample ID 대응표](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE150728), [원 논문 Figure 1a](https://pmc.ncbi.nlm.nih.gov/articles/PMC7382903/figure/F1/).
같은 환자의 두 시점을 서로 다른 환자로 세지 않습니다. 두 군은 연령·성별 등을 맞춘 임상 비교군이 아니므로 교육용 탐색으로 해석합니다.

### (1) 네 파일 다운로드

위 GEO sample 링크의 **Supplementary file → (http)**에서 각각 다운로드합니다.
전체 `GSE150728_RAW.tar`를 받았다면 tar를 풀고 아래 네 파일만 사용합니다.

### (2) `.gz`를 한 번 풀어 `.rds`로 준비

**이번 데이터는 이전 10x 데이터와 파일 구조가 다릅니다.** GEO의 `.rds.gz`는 압축된 RDS 파일을 다시 gzip으로 감싼 형태입니다. 7-Zip 등의 프로그램으로 바깥 `.gz`를 한 번 풀어 아래 이름의 `.rds` 파일을 만듭니다. `.rds` 내부 압축은 R이 처리하므로 더 풀지 않습니다.

RStudio에서 실습 폴더를 Project로 열고 아래처럼 배치합니다. Sample별 하위 폴더나 `matrix.mtx` 파일은 필요하지 않습니다.

```text
실습폴더/
├─ workshop.Rproj
└─ data/
   └─ GSE150728_4sample/
      ├─ GSM4557334_HIP002_cell.counts.matrices.rds
      ├─ GSM4557335_HIP015_cell.counts.matrices.rds
      ├─ GSM4557330_557_cell.counts.matrices.rds
      └─ GSM4557331_558_cell.counts.matrices.rds
```

Seq-Well은 작은 well에 세포와 barcode bead를 담는 방법입니다. 이 실습은 이미 만들어진 UMI count matrix에서 시작하며 FASTQ 처리나 alignment는 포함하지 않습니다.

## 실습 목표와 실행 방법:
**QC → SCTransform → integration → clustering → UMAP → annotation**을 진행합니다. 이어서 T cell, module score, cell composition을 살펴봅니다.

- 별도 R 파일 없이 이 문서의 코드 블록을 **하나씩 순서대로** 실행합니다.
- 기본 실습은 1–10절, 추가 실습은 11–13절입니다. 마지막에 14절에서 저장합니다.
- `선택`으로 표시한 접힌 부분은 필요할 때만 실행합니다.
- 오류가 나면 다음 블록으로 넘어가지 말고 강사와 함께 확인합니다.
- `<-`는 결과에 이름을 붙이는 기호입니다. 예를 들어 `rawdata`는 합친 데이터를 담는 이름입니다.
- 강의 목표 환경은 **R 4.6.1 / Seurat 5.5.1**입니다. 이 문서의 로컬 검증 환경은 R 4.5.2 / Seurat 5.3.1이므로, 수업 전 강사의 목표 환경에서 다시 실행해 확인합니다.

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

이 명령은 실행 시점 CRAN 버전을 설치하므로 특정 버전을 고정하지 않습니다. 강사와 다른 버전이면 수업 전에 맞춥니다. Windows에서 source package를 직접 컴파일할 때는 사용 중인 R에 맞는 Rtools가 필요할 수 있습니다.

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

데이터 경로와 저장 폴더를 정합니다. RStudio Project를 열었다면 상대경로를 그대로 사용할 수 있습니다.

```r
getwd()
data_dir <- "data/GSE150728_4sample"
dir.create("results", showWarnings = FALSE)
list.files(data_dir)
```

**확인:** 위에서 준비한 `.rds` 파일 네 개가 보이나요? 다른 곳에 저장했다면 `data_dir`만 실제 경로로 바꿉니다. Windows 경로는 `C:/Users/...`처럼 `/`를 사용합니다.

## 3. Seq-Well RDS에서 Seurat object 만들기

### 3-1. 네 sample 읽기

`readRDS()`는 저장된 R object를 읽습니다. 파일 하나에는 `exon`, `intron`, `spanning`이라는 세 count matrix가 있습니다. 이번 유전자 발현 실습은 **exon matrix**를 사용합니다. RNA velocity 분석은 하지 않으며 세 matrix를 임의로 더하지 않습니다.

```r
H1_data <- readRDS(file.path(data_dir, "GSM4557334_HIP002_cell.counts.matrices.rds"))
H2_data <- readRDS(file.path(data_dir, "GSM4557335_HIP015_cell.counts.matrices.rds"))
C3_data <- readRDS(file.path(data_dir, "GSM4557330_557_cell.counts.matrices.rds"))
C4_data <- readRDS(file.path(data_dir, "GSM4557331_558_cell.counts.matrices.rds"))
names(H1_data)
```

`$exon`은 목록에서 exon matrix를 꺼낸다는 뜻입니다. **행은 gene, 열은 cell barcode**이며 각 숫자는 해당 gene의 UMI 수입니다.

```r
H1_counts <- H1_data$exon
H2_counts <- H2_data$exon
C3_counts <- C3_data$exon
C4_counts <- C4_data$exon
dim(H1_counts)
H1_counts[1:5, 1:5]
```

### 3-2. Sample별 object 만들기

`min.cells = 3`은 해당 sample에서 최소 3개 cell에 검출된 gene을 남깁니다.

```r
H1 <- CreateSeuratObject(H1_counts, project = "H1", min.cells = 3)
H2 <- CreateSeuratObject(H2_counts, project = "H2", min.cells = 3)
C3 <- CreateSeuratObject(C3_counts, project = "C3", min.cells = 3)
C4 <- CreateSeuratObject(C4_counts, project = "C4", min.cells = 3)
```

### 3-3. 네 sample 합치기

Cell 이름에 sample 접두어를 붙여 같은 barcode를 구분합니다.

```r
rawdata <- merge(H1, y = list(H2, C3, C4),
                 add.cell.ids = c("H1", "H2", "C3", "C4"))
table(rawdata$orig.ident)
Layers(rawdata[["RNA"]])
```

`orig.ident`는 출신 sample, layer는 sample별 측정값을 보관하는 칸입니다. 이 실습에서는 sample 하나가 donor 한 명입니다.

```r
condition_map <- c(H1 = "Healthy", H2 = "Healthy",
                   C3 = "Severe_COVID", C4 = "Severe_COVID")
rawdata$condition <- unname(condition_map[as.character(rawdata$orig.ident)])
table(rawdata$orig.ident, rawdata$condition)
```

네 개의 sample이 두 조건에 올바르게 배치되었는지 확인합니다. 처음 읽은 중간 object를 비워 메모리를 확보합니다.

```r
rm(H1_data, H2_data, C3_data, C4_data,
   H1_counts, H2_counts, C3_counts, C4_counts, H1, H2, C3, C4)
gc()
```

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
) + geom_hline(yintercept = c(1000, 15000))

VlnPlot(
  rawdata, features = c("percent.mt"),
  group.by = "orig.ident", layer = "counts", pt.size = 0
) + geom_hline(yintercept = c(20))

```

### 4-2. 기준에 맞는 cell 남기기

아래는 **워크숍용 QC 시작 기준**입니다. 원 논문의 UMI 1,000–15,000 및 mitochondrial 20% 기준을 참고하고, 실습에서는 nFeature 600–5,000을 함께 사용합니다. 원 논문의 rRNA·complexity 필터와 후속 수동 제거까지 재현하는 분석은 아닙니다. 먼저 sample별 분포를 보고 기준을 점검합니다. 모든 조직에 그대로 적용하지 않습니다. UMI가 높다는 이유만으로 doublet을 확정할 수 없고, mitochondrial 비율도 cell type에 따라 달라집니다.

```r
filtdata <- subset(
  rawdata,
  subset = nFeature_RNA > 600 & nFeature_RNA < 5000 &
    nCount_RNA >= 1000 & nCount_RNA <= 15000 & percent.mt < 20
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

아래는 공개 파일의 exon matrix에 위 기준을 적용해 확인한 값입니다. QC 전 열에는 낮은 UMI를 가진 barcode도 포함되어 있어 모두 양질의 세포라는 뜻은 아닙니다.

| Sample | QC 전 barcode 수 | QC 후 cell 수 |
|---|---:|---:|
| H1 | 9,070 | 2,054 |
| H2 | 5,330 | 2,257 |
| C3 | 16,536 | 6,892 |
| C4 | 21,215 | 3,182 |
| 합계 | 52,151 | 14,385 |

<details>
<summary>선택: 노트북 실습을 위해 sample당 1,000개씩 사용하기</summary>

시간이나 메모리가 부족할 때 **5절 전에 한 번만** 실행합니다. 각 donor에서 무작위로 1,000개씩 골라 총 4,000개를 사용합니다. 기본 분석은 QC를 통과한 전체 cell을 사용합니다.

```r
Idents(filtdata) <- "orig.ident"
set.seed(12345)
filtdata <- subset(filtdata, downsample = 1000)
table(filtdata$orig.ident)
```

Downsampling을 하면 cluster 번호와 그림이 달라지고 드문 cell type이 줄어들 수 있습니다. 11절의 composition은 선택된 cell 안에서의 비율입니다. 강사와 학생은 전체 분석 또는 4,000개 분석 중 같은 방식을 사용합니다.
</details>


## 5. SCTransform과 PCA

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


### 중간 결과 저장

계산 시간은 컴퓨터와 cell 수에 따라 달라집니다. 아래는 **지금 분석한 COVID 데이터**를 저장하는 코드입니다. 이전 치주염 실습의 `filtdata.RData`는 사용하지 않습니다.

```r
saveRDS(filtdata, "results/GSE150728_SCT_PCA.rds")
```

<details>
<summary>선택: 나중에 새 session에서 이 단계부터 재개하기</summary>

2절의 패키지를 먼저 불러온 뒤 실행합니다. 방금 저장한 경우 바로 실행할 필요는 없습니다.

```r
filtdata <- readRDS("results/GSE150728_SCT_PCA.rds")
condition_map <- c(H1 = "Healthy", H2 = "Healthy",
                   C3 = "Severe_COVID", C4 = "Severe_COVID")
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
  `T cell` = c("CD3D", "CD3E", "TRAC"),
  `NK` = c("NKG7", "GNLY", "PRF1"),
  `B cell` = c("MS4A1", "CD79A", "CD37"),
  Monocyte = c("LYZ", "S100A8", "CD14", "FCGR3A"),
  Plasmablast = c("MZB1", "IGJ", "XBP1"),
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

각 panel에 marker가 남아 있는지 먼저 확인합니다. `lapply`는 목록의 각 항목에 같은 작업을 적용하며, `intersect`는 실제 데이터에 있는 gene만 남깁니다. 이 데이터에서는 JCHAIN의 이전 symbol인 `IGJ`를 사용합니다. Neutrophil-like는 검토할 후보 이름이며 lineage 전환을 뜻하지 않습니다.

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
</details>

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

작성한 이름을 cell별로 붙입니다.

```r
intdata$celltype <- unname(celltype_map[as.character(intdata$rpca_clusters)])
table(intdata$rpca_clusters, intdata$celltype)
annotated_umap <- DimPlot(intdata, reduction = "umap.rpca",
                          group.by = "celltype", label = TRUE)
annotated_umap
```

**완료 기준:** 주요 cluster에 근거 marker 두 개 이상을 설명할 수 있나요? 전부 `Unassigned`라면 annotation을 마친 후 다음으로 넘어갑니다. T cell을 하나 이상 확인해야 12–13절을 진행할 수 있습니다.

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

**질문:** C3과 C4에서 같은 방향의 차이가 보이나요? Healthy 2명과 Severe_COVID 2명의 탐색적 결과입니다. Cell 수천 개가 독립 환자 수천 명을 뜻하지 않습니다. 채취·분리·QC에 따른 편향도 고려합니다.


---

### 11-1. 같은 Monocyte 안에서 발현 비교하기

전체 PBMC 평균의 차이가 세포 구성 때문인지, 같은 종류의 세포 내부 변화인지 구분해 봅니다. 아래 코드는 10절에서 `Monocyte`라고 붙인 cell을 사용합니다.

```r
mono <- subset(intdata, subset = celltype == "Monocyte")
table(mono$orig.ident)
DotPlot(mono, features = c("CD14", "FCGR3A", "HLA-DRA", "HLA-DPB1",
                           "HLA-DMA", "ISG15", "IFIT3", "IFI27"),
        group.by = "orig.ident") + RotatedAxis()
```

**질문:** HLA class II gene과 IFN-response gene의 패턴이 C3·C4에서 같은가요? Monocyte 내부에서도 CD14/CD16 구성 차이가 영향을 줄 수 있습니다. 네 sample 모두 Monocyte가 있는지 확인하고, 없다면 annotation부터 점검합니다.

원 연구에서 관찰한 HLA-II 감소, 환자별 IFN 반응 차이를 탐색하는 활동입니다. 이번 네 donor의 재분석이 논문 전체 결과와 같을 필요는 없습니다. 그룹당 두 명이므로 환자별 패턴을 보여주며, 세포를 독립 환자로 취급하는 질환 DEG 검정은 하지 않습니다.

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

```r
tcell[["RNA"]] <- JoinLayers(tcell[["RNA"]])
tcell[["RNA"]] <- split(tcell[["RNA"]], f = tcell$orig.ident)
tcell <- SCTransform(tcell, new.assay.name = "SCT_T",
                      conserve.memory = TRUE, seed.use = 12345)
tcell <- RunPCA(tcell, npcs = 30, reduction.name = "pca.tcell", seed.use = 12345)
```

### 12-2. 다시 통합하고 묶기

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

이 연구에서 IFN response는 환자마다 달랐습니다. Severe_COVID 두 명 모두 같은 방향일 것이라고 가정하지 않습니다. Exhaustion-associated score는 탐색용이며 이 연구에서 CD8 T cell exhaustion 증가가 확정된 것은 아닙니다.

점수는 control gene과 비교한 상대값입니다. 그림마다 색 범위가 다를 수 있으며, 다른 program끼리 점수 크기를 직접 비교하지 않습니다. 각 program 내 패턴을 marker·sample 정보와 함께 해석합니다.

같은 cell type 안에서도 **환자별**로 IFN response를 비교합니다.

```r
VlnPlot(tcell, features = "IFN_response_score", group.by = "orig.ident", pt.size = 0)
```

T cell subtype 비율이 다르면 전체 T cell 점수도 달라질 수 있습니다. 조건별 차이가 보이면 같은 subtype 안에서도 확인해야 합니다.

T cell 추가 실습을 마쳤다면 저장합니다.

```r
saveRDS(tcell, "results/GSE150728_Tcell_Seurat5.rds")
```


## 14. 결과 저장과 재개

기본 실습 결과와 실행 환경을 저장합니다. 같은 파일 이름으로 저장하면 이전 결과를 갱신합니다.

```r
saveRDS(intdata, "results/GSE150728_4sample_Seurat5.rds")
writeLines(capture.output(sessionInfo()), "results/sessionInfo.txt")
```

그림은 RStudio의 **Plots → Export**로 저장할 수 있습니다. 아래 자동 저장은 선택입니다.

<details>
<summary>선택: 기본 실습 그림과 표 저장</summary>

1–10절을 처음부터 마친 뒤 실행합니다. Checkpoint에서 시작했다면 만들지 않은 QC 그림과 `qc_summary`의 저장 줄은 건너뜁니다. Width와 height는 inch 단위의 그림 크기입니다.

```r
ggsave("results/01_QC_violin_before_filtering.png", qc_violin, width = 14, height = 5)
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
<summary>선택: T cell 그림 저장 — 12–13절 완료 후</summary>

```r
ggsave("results/08_Tcell_UMAP.png", tcell_umap, width = 9, height = 7)
ggsave("results/09_Tcell_marker_DotPlot.png", tcell_dotplot, width = 15, height = 8)
ggsave("results/10_Tcell_program_scores.png", program_plot, width = 15, height = 5)
```

</details>

<details>
<summary>선택: Composition 저장 — 11절 완료 후</summary>

```r
ggsave("results/07_celltype_composition.png", composition_plot, width = 9, height = 6)
write.csv(composition_table, "results/celltype_composition_by_sample.csv", row.names = FALSE)
```

</details>

<details>
<summary>선택: 새 R session에서 저장한 object 불러오기</summary>

같은 Project를 열고 2절의 패키지를 불러온 뒤 실행합니다. 처음 분석할 때는 생략합니다.

```r
intdata <- readRDS("results/GSE150728_4sample_Seurat5.rds")
condition_map <- c(H1 = "Healthy", H2 = "Healthy",
                   C3 = "Severe_COVID", C4 = "Severe_COVID")
table(intdata$celltype)
```

</details>

## 막혔을 때 확인하기

- **파일을 못 찾음:** Project 위치와 `data_dir`, 네 `.rds` 파일명을 확인합니다.
- **object가 없음:** 앞 블록을 실행했는지 확인합니다. Session을 재시작하면 메모리의 object가 사라집니다.
- **MT- gene이 0개:** 사람 gene symbol 대신 다른 ID를 읽었는지 확인합니다.
- **unknown input format:** 원본 `.rds.gz`의 바깥 gzip을 한 번 풀고 `.rds`를 `readRDS()`로 읽습니다. 이 데이터에 `Read10X()`를 사용하지 않습니다.
- **RDS 결과가 list:** 정상입니다. `names(H1_data)`를 확인하고 `$exon`을 선택합니다.
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
- [GSE150728](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE150728)
- Wilk AJ, et al. *A single-cell atlas of the peripheral immune response in patients with severe COVID-19*. Nature Medicine. 2020. [doi:10.1038/s41591-020-0944-y](https://doi.org/10.1038/s41591-020-0944-y)
