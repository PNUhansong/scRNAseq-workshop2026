# scRNA-seq analysis workshop with Seurat v5

single-cell RNA-seq 실습 자료입니다. \
사람 PBMC 4개 샘플을 읽고, 세포별 QC부터 Seurat v5의 layer 기반 integration, UMAP, clustering, marker 기반 annotation을 진행해보겠습니다.

### GEO database (치주염-당뇨 환자의 PBMC scRNAseq)
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE244515

\
### 실습 환경
- 기준 환경: **R 4.6.1**
- 기준 패키지: **Seurat 5.5.1**
- 실습 데이터: **GSE244515**의 H1, H2, PD1, PD2
- 분석 흐름: `10x matrix → QC → SCTransform v2 → RPCA integration → clustering → UMAP → marker 확인`

> 이 자료는 2026년 9월의 CRAN 및 Seurat 공식 문서를 기준으로 작성했습니다. 

## 파일 구성

```text
├─ data/                         # 직접 준비, GitHub에는 올리지 않음
   └─ GSE244515_4sample/
      ├─ H1/
      ├─ H2/
      ├─ PD1/
      └─ PD2/
```

각 샘플 폴더에는 아래 파일 세 개가 있어야 합니다. `.gz` 파일은 압축을 풀지 않습니다.

```text
barcodes.tsv.gz
features.tsv.gz
matrix.mtx.gz
```

## 1. R과 RStudio 준비

R 4.6.1을 설치한 뒤 RStudio를 실행합니다. 이전 R 버전에서 설치한 package library를 그대로 복사하기보다 R 4.6.1에서 패키지를 다시 설치하는 편이 안전합니다.

RStudio에서 새 Project를 만들고 이 저장소의 최상위 폴더를 선택합니다. 이후 코드가 `data/`와 `results/`를 상대 경로로 찾기 때문에 `setwd("D:/...")`를 직접 수정할 필요가 없습니다.

처음 한 번만 다음 파일을 실행합니다. Seurat가 설치되어 있더라도 5.5.1보다 오래된 버전이면 이 스크립트가 업데이트합니다.

```r
source("00_install_packages.R")
```

설치가 끝난 뒤 R session을 다시 시작하고 버전을 확인합니다.

```r
R.version.string
packageVersion("Seurat")
packageVersion("SeuratObject")
```

2026년 9월 CRAN의 Seurat 최신 버전은 5.5.1이며 R 4.0.0 이상을 요구합니다. 이 실습은 R 4.6.1과 Seurat 5.5.1 조합을 목표로 하며, 설치 스크립트는 Seurat 6이 설치된 경우에도 멈춰 버전 차이를 확인하도록 합니다.

## 2. 데이터 준비

[NCBI GEO GSE244515](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE244515)에서 다음 네 샘플의 10x count matrix 파일을 준비합니다.

| 실습 이름 | GEO sample | 연구 집단 |
|---|---|---|
| H1 | GSM7818495 | healthy control 1 |
| H2 | GSM7818496 | healthy control 2 |
| PD1 | GSM7818506 | periodontitis 1 |
| PD2 | GSM7818507 | periodontitis 2 |

GEO에서 받은 파일은 샘플 이름이 앞에 붙어 있습니다. 각 샘플 폴더로 옮긴 뒤 `Read10X()`가 읽을 수 있도록 다음처럼 이름을 맞춥니다.

```text
GSM7818495_H1_barcodes.tsv.gz  → H1/barcodes.tsv.gz
GSM7818495_H1_features.tsv.gz  → H1/features.tsv.gz
GSM7818495_H1_matrix.mtx.gz    → H1/matrix.mtx.gz
```

H2, PD1, PD2도 같은 방법으로 정리합니다.

## 3. 전체 코드 실행

RStudio에서 `01_scRNAseq_workshop_Seurat5.R`을 엽니다. 처음에는 한 번에 `Source`하기보다 section별로 실행하고, 각 단계에서 object와 그림을 확인하는 것을 권장합니다.

전체 파일을 실행하려면 다음 명령을 사용합니다.

```r
source("01_scRNAseq_workshop_Seurat5.R", echo = TRUE)
```

## 4. Step 1: 10x matrix와 Seurat object

`Read10X()`가 읽은 sparse matrix에서 행은 gene, 열은 cell barcode입니다.

```r
counts_list <- lapply(sample_dirs, read_gene_expression)
dim(counts_list$H1)
counts_list$H1[1:5, 1:5]
```

샘플마다 object를 만든 뒤 `merge()`합니다. Seurat v5는 합쳐진 RNA assay 안에 샘플별 count를 여러 layer로 보관할 수 있습니다.

```r
Layers(rawdata[["RNA"]])
```

여기서 layer는 같은 종류의 측정값을 샘플별로 나누어 둔 칸이라고 생각하면 됩니다.

## 5. Step 2: 세포별 quality control

세 가지 기본 지표를 함께 봅니다.

| 지표 | 의미 | 주의해서 볼 경우 |
|---|---|---|
| `nFeature_RNA` | 한 cell barcode에서 검출된 gene 종류 수 | 너무 낮으면 정보가 부족할 수 있음 |
| `nCount_RNA` | 한 cell barcode의 전체 UMI 수 | 지나치게 높으면 doublet 가능성을 점검 |
| `percent.mt` | 전체 UMI 중 mitochondrial gene UMI 비율 | 높으면 stress 또는 cytoplasmic RNA 손실의 단서가 될 수 있음 |

실습 코드의 cutoff는 다음과 같습니다.

```r
nFeature_RNA > 600
nFeature_RNA < 5000
nCount_RNA < 25000
percent.mt < 20
```

이 값은 네 PBMC 샘플을 이용한 교육용 출발점입니다. `percent.mt = 20` 같은 숫자를 모든 조직의 정답처럼 사용하면 안 됩니다. 정상적인 범위는 species, tissue, cell type, 실험 방법에 따라 달라질 수 있습니다.

또한 이 filtering만으로 doublet이나 ambient RNA를 모두 제거할 수 없습니다. 이번 기초 실습에서는 세 가지 QC 지표의 의미에 집중하고, doublet detection과 ambient RNA correction은 후속 분석 주제로 남깁니다.

## 6. Step 3: SCTransform v2

먼저 RNA assay를 샘플별 layer로 나눕니다.

```r
filtdata[["RNA"]] <- JoinLayers(filtdata[["RNA"]])
filtdata[["RNA"]] <- split(filtdata[["RNA"]], f = filtdata$orig.ident)
Layers(filtdata[["RNA"]])
```

이후 `SCTransform()`을 실행합니다.

```r
filtdata <- SCTransform(
  filtdata,
  vst.flavor = "v2",
  variable.features.n = 3000,
  conserve.memory = TRUE
)
```

SCTransform은 세포마다 검출된 전체 분자 수 같은 기술적 차이를 모델링하면서, 생물학적 변이를 분석하기 좋은 형태로 바꿉니다. 모든 세포의 발현을 똑같게 만드는 과정은 아닙니다. Seurat v5에서는 SCT v2가 기본입니다.

## 7. Step 4: integration 전 지도

PCA 결과를 사용해 integration 전 UMAP을 만듭니다.

```r
filtdata <- RunPCA(filtdata, assay = "SCT", npcs = 50)
filtdata <- RunUMAP(
  filtdata,
  reduction = "pca",
  dims = 1:30,
  reduction.name = "umap.unintegrated"
)
```

UMAP에서 점 하나는 측정된 cell 하나입니다. 가까운 점은 발현 패턴이 비슷한 경향이 있지만, UMAP 1과 UMAP 2는 특정 gene의 발현량이나 실제 조직 좌표가 아닙니다.

## 8. Step 5: Seurat v5 integration

```r
intdata <- IntegrateLayers(
  object = filtdata,
  method = RPCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.rpca",
  assay = "SCT",
  normalization.method = "SCT",
  dims = 1:30
)
```

RPCA integration은 sample마다 공통으로 나타나는 cell state와 cell type 구조를 찾아 비교할 수 있도록 정렬합니다. 원래 RNA count를 지우는 과정은 아니며, 결과를 `integrated.rpca`라는 새로운 dimensional reduction에 저장합니다.

Integration은 질병에 따른 모든 차이를 없애는 자동 보정이 아닙니다. 결과가 sample별로 잘 섞였는지와 알려진 marker 패턴이 유지되는지를 함께 확인해야 합니다.

## 9. Step 6: clustering과 UMAP

```r
intdata <- FindNeighbors(
  intdata,
  reduction = "integrated.rpca",
  dims = 1:30
)
intdata <- FindClusters(
  intdata,
  resolution = 0.2,
  cluster.name = "rpca_clusters"
)
intdata <- RunUMAP(
  intdata,
  reduction = "integrated.rpca",
  dims = 1:30,
  reduction.name = "umap.rpca"
)
```

`FindClusters()`는 비슷한 이웃 관계를 가진 cell들을 묶습니다. cluster 번호는 이름표일 뿐이며 숫자의 크기에 생물학적 순서가 있는 것은 아닙니다.

## 10. Step 7: marker로 cell type 추론

Integration에 사용한 reduction과 marker 발현값의 역할을 구분해야 합니다.

- `integrated.rpca`: 이웃 탐색, clustering, UMAP에 사용
- `RNA`: 실제 marker 발현을 확인할 때 사용

RNA의 sample별 layer를 합치고 normalized data layer를 만듭니다.

```r
intdata[["RNA"]] <- JoinLayers(intdata[["RNA"]])
intdata <- NormalizeData(intdata, assay = "RNA")
DefaultAssay(intdata) <- "RNA"
```

`DotPlot()`에서는 점의 크기가 해당 cluster에서 gene이 검출된 cell 비율을, 색이 평균 발현량을 나타냅니다. marker 하나만으로 cell type을 결정하지 말고 여러 marker의 조합과 알려진 생물학을 함께 확인합니다.

스크립트의 `celltype_map`은 의도적으로 비어 있습니다.

```r
celltype_map <- character(0)
```

DotPlot을 확인한 뒤 현재 결과에 맞게 작성합니다.

```r
celltype_map <- c(
  "0" = "CD4 T",
  "1" = "Monocyte"
  # 나머지 cluster도 marker 근거를 확인한 뒤 추가합니다.
)
```


## 11. T cell만 확대해서 다시 분석하기

큰 PBMC 지도에서는 T cell과 B cell처럼 서로 다른 계통의 차이가 먼저 보입니다. T cell 내부의 작은 차이를 보기 위해서는 T cell만 분리한 뒤 SCT, PCA, integration, clustering과 UMAP을 다시 계산합니다.

먼저 `celltype_map`에서 실제로 사용한 이름에 맞게 다음 항목을 수정합니다.

```r
T_CELL_LABELS <- c("T cell", "CD4 T", "CD8 T", "Treg")
RUN_TCELL_DEEP_DIVE <- TRUE
```

스크립트는 해당 annotation의 cell만 선택하고 `tcell` object를 만듭니다. 전체 PBMC에서 계산한 PCA를 그대로 확대하는 것이 아니라 T cell 안에서 variable gene과 주요 변이를 다시 찾습니다.

T cell subset은 전체 PBMC보다 작기 때문에 RPCA integration에서 `k.weight = 50`을 사용합니다. 오류가 발생하면 sample별 T cell 수와 공유되는 population을 먼저 확인합니다.

T cell marker는 다음처럼 두 종류로 나누어 해석합니다.

| 구분 | 예시 | 대표 marker |
|---|---|---|
| 비교적 안정적인 type 또는 subtype | Naive/memory-like, cytotoxic T, Treg | `CCR7`, `TCF7`, `CCL5`, `FOXP3`, `IL2RA` |
| 변화할 수 있는 state 또는 program | Proliferation, IFN response, exhaustion-associated program | `MKI67`, `ISG15`, `PDCD1`, `TOX` |

`PDCD1`이나 `TIGIT` 하나가 검출되었다고 exhausted T cell로 확정하지 않습니다. activation과 exhaustion-associated state가 marker를 공유할 수 있고, scRNA-seq에는 dropout이 있기 때문입니다. 여러 marker의 조합과 질환 맥락을 함께 확인해야 합니다.

학생들에게 다음 질문을 제시할 수 있습니다.

1. Cytotoxicity가 높은 cluster는 `CD8A`도 함께 높은가?
2. IFN response는 하나의 T cell subtype에만 나타나는가, 여러 subtype에 걸쳐 나타나는가?
3. Treg의 면역 억제 기능은 자가면역과 암에서 각각 어떤 결과를 만들 수 있는가?
4. Cytotoxic T cell의 기능은 감염 제거와 조직 손상 중 어느 한쪽으로만 설명할 수 있는가?

## 12. Functional program을 module score로 비교하기

`AddModuleScore()`는 여러 gene의 발현 경향을 cell별 하나의 상대적인 점수로 요약합니다. 이 실습에서는 다음 세 program을 계산합니다.

```r
tcell_programs <- list(
  Cytotoxicity = c("NKG7", "CCL5", "PRF1", "GZMB"),
  IFN_response = c("ISG15", "IFIT1", "IFIT3", "MX1"),
  Exhaustion_associated = c("PDCD1", "LAG3", "HAVCR2", "TOX", "TIGIT")
)
```

Module score는 절대적인 활성도나 임상적 진단값이 아닙니다. 비슷한 평균 발현량을 가진 control gene과 비교한 상대 점수입니다. 따라서 점수 하나로 세포 상태를 확정하기보다 UMAP 위치, marker DotPlot과 sample 정보를 함께 봅니다.

## 13. Cell composition은 sample별로 비교하기

각 sample에서 cell type별 cell 수와 비율을 계산합니다.

```r
composition_table <- intdata[[]] |>
  count(orig.ident, condition, celltype, name = "cell_count") |>
  group_by(orig.ident) |>
  mutate(cell_proportion = cell_count / sum(cell_count))
```

막대 하나는 cell 하나가 아니라 **환자 또는 sample 하나**를 나타냅니다. PD1과 PD2에서 같은 방향의 변화가 나타나는지 확인하고, 한 sample이 결과를 주도하지 않는지 비교합니다.

이 데이터는 Healthy 2명과 Periodontitis 2명만 사용하므로 cell composition 그림은 탐색적 결과입니다. cell이 수천 개 있더라도 biological replicate가 수천 개가 되는 것은 아닙니다.

## 14. 결과 파일

코드를 실행하면 `results/`에 다음 파일이 생성됩니다.

```text
01_QC_violin_before_filtering.png
02_QC_scatter_before_filtering.png
QC_cell_numbers.csv
03_UMAP_before_integration.png
04_UMAP_after_integration.png
05_marker_DotPlot.png
GSE244515_4sample_Seurat5.rds
sessionInfo.txt
```

`celltype_map`을 작성하면 다음 결과도 추가됩니다.

```text
06_UMAP_annotated.png
07_celltype_composition.png
celltype_composition_by_sample.csv
08_Tcell_UMAP.png
09_Tcell_marker_DotPlot.png
10_Tcell_program_scores.png
GSE244515_Tcell_Seurat5.rds
```

T cell 관련 결과는 `T_CELL_LABELS`와 일치하는 annotation이 있을 때 생성됩니다.

## 자주 생기는 문제

### `Read10X()`가 파일을 찾지 못하는 경우

각 샘플 폴더 안에 `barcodes.tsv.gz`, `features.tsv.gz`, `matrix.mtx.gz`가 정확히 있는지 확인합니다. `.gz`는 풀지 않습니다.

```r
list.files("data/GSE244515_4sample/H1")
```

### R을 업데이트했는데 Seurat가 보이지 않는 경우

R 4.6.1은 이전 R과 다른 package library를 사용할 수 있습니다. 새 R session에서 `00_install_packages.R`을 다시 실행합니다.

### `data layers are not joined`와 비슷한 메시지가 나오는 경우

여러 RNA layer를 대상으로 marker 또는 differential expression을 실행하기 전에 다음 코드를 사용합니다.

```r
intdata[["RNA"]] <- JoinLayers(intdata[["RNA"]])
```

### 메모리가 부족한 경우

다른 큰 object를 닫고 R session을 다시 시작합니다. `SCTransform(..., conserve.memory = TRUE)`를 유지하고, 강사가 미리 만든 filtered object부터 시작하는 방법도 사용할 수 있습니다.

### `Number of anchor cells is less than k.weight` 오류가 나는 경우

QC 후 특정 sample에 남은 cell 수가 지나치게 적거나 sample 사이에 공유되는 cell population이 부족한지 먼저 확인합니다. 의도적으로 아주 작은 연습 데이터를 사용했다면 `IntegrateLayers()`에 `k.weight = 50`처럼 기본값보다 작은 값을 추가할 수 있습니다. 오류 메시지에 표시된 anchor 수보다 작게 정하되, 단순히 오류를 없애기 위한 숫자로 사용하지 말고 integration 결과에서 marker 구조가 보존되는지도 확인합니다.

## 해석할 때 기억할 점

- UMAP은 고차원 발현 패턴을 2차원으로 요약한 그림입니다.
- cluster는 계산 결과이고 cell type은 marker와 생물학적 지식을 이용한 해석입니다.
- `proliferating`, `interferon response`, `Exhaustion` 같은 표현은 cell type과 다른 **cell state**일 수 있습니다.
- 질환군 비교에서는 여러 cell을 각각 독립된 환자처럼 취급하면 안 됩니다. 환자 또는 sample 단위의 biological replicate를 고려해야 합니다.

## 참고 자료

- [Seurat v5 Essential Commands](https://satijalab.org/seurat/articles/seurat5_essential_commands.html)
- [Seurat v5 Integrative Analysis](https://satijalab.org/seurat/articles/seurat5_integration.html)
- [Using sctransform in Seurat](https://satijalab.org/seurat/articles/sctransform_vignette)
- [Seurat AddModuleScore](https://satijalab.org/seurat/reference/addmodulescore)
- [NCBI GEO GSE244515](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE244515)
- Lee H, Joo J, Song J, et al. *Immunological link between periodontitis and type 2 diabetes deciphered by single-cell RNA analysis*. Clinical and Translational Medicine. 2023. [doi:10.1002/ctm2.1503](https://doi.org/10.1002/ctm2.1503)



- [NCBI GEO GSE244515](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE244515)
- Lee H, Joo J, Song J, et al. *Immunological link between periodontitis and type 2 diabetes deciphered by single-cell RNA analysis*. Clinical and Translational Medicine. 2023. [doi:10.1002/ctm2.1503](https://doi.org/10.1002/ctm2.1503)
