# synthMEs - Synthetic ME Toolkit

A compact set of R functions to build **synthetic module eigengenes (MEs)** in a target cohort from a template network, **compare** eigengenes across cohorts with optional covariate adjustment, and **clean/repair** module color assignments using **kME-based** rules.

Functions included:

* `getSynthMEs()` — construct synthetic MEs in a target dataset using a template network’s modules and kME ranking
* `compareMEs()` — run one-way (or covariate-adjusted) ANOVA per ME in template and target, and create comparative boxplot + beeswarm PDFs
* `kMEcleanup()` — iteratively repair/lock module color assignments based on kME intramodule strength and max-minus-intra gaps (with optional hub preservation)

These algorithms have been published for use in inferring module eigengenes in -omic data of tissues or cells similar to other omics data which already has a coexpression network constructed, now made available as functions for easier reproducbility.

---

## Installation & Requirements

* **R >= 4.1** (recommended)
* Packages:

  * **WGCNA** (eigengenes, kME, color utilities)
  * **beeswarm** (only for `compareMEs()` plotting)
  * Base `stats`, `utils`, `grDevices`

Install dependencies:

```r
install.packages(c("WGCNA","beeswarm"))
```

---

## Data Conventions

* Expression matrices: **rows = features/assays/proteins, columns = samples**
* Colors/vectors are aligned to **rows of the expression matrix**
* MEs: **rows = samples, columns = modules**; columns may be named like `"MEturquoise"` or `"turquoise"` (functions handle either)

---

## 1) `getSynthMEs()`

### Purpose

Build **synthetic eigengenes** in a **target** cohort using a **template** cohort’s module structure. For each template module, pick top-kME members (with a rescue rule), intersect with features present in the target cohort, and compute eigengenes on the target.

### Signature

```r
getSynthMEs(
  minimumMEmembers = 4,
  topPercent = 0.20,
  minKmeThreshToRescue = 0.70,
  cleanDat.template,
  netColors.template,
  cleanDat.target,
  matchCSV = "Synthetic_eigengene_exact_match_members.csv"
)
```

### Arguments

* `minimumMEmembers` (int, default 4): Minimum markers per module to compute a synthetic ME in target
* `topPercent` (0–1, default 0.20): Start with top X% kME members from template module
* `minKmeThreshToRescue` (0–1, default 0.70): If topPercent yields too few markers, **rescue** by taking all members with kME >= this threshold
* `cleanDat.template` (matrix/df): Template expression (rows=features, cols=samples)
* `netColors.template` (vector): Template module colors; length = `nrow(cleanDat.template)`
* `cleanDat.target` (matrix/df): Target expression (rows=features, cols=samples)
* `matchCSV` (string): CSV filename to write **per-module feature members** actually used in the synthetic ME computation

### Returns

* A **data.frame** `MEsSynth` (rows = target samples, columns = modules) ordered to match the template’s non-grey modules.

### Side effects

* Prints dimension comparison of `MEsSynth` vs template MEs
* Prints whether **all** modules generated synthetic MEs (or how many are missing)
* **Writes** `matchCSV` listing marker IDs used per module (columns)

### Notes & Behavior

* Drops `grey` from the template’s ME list before matching/ordering
* Feature matching uses **rownames** of `cleanDat.target`
* If **no modules** meet criteria, returns an **empty** `MEsSynth`, prints messages, and writes an **empty** CSV

### Minimal Example

```r
MEsSynth <- getSynthMEs(
  minimumMEmembers = 4,
  topPercent = 0.20,
  minKmeThreshToRescue = 0.70,
  cleanDat.template = cleanDatMEGA,
  netColors.template = netMEGA$colors,
  cleanDat.target = cleanDat,
  matchCSV = "Synthetic_eigengene_exact_match_members.csv"
)
```

---

## 2) `compareMEs()`

### Purpose

Compare **template** vs **target** eigengenes by performing **one-way ANOVA** (or ANOVA adjusted for optional covariates `Sex`, `PMI`) **per ME**, then draw template/target **boxplots** with **beeswarm** overlays and highlight significant p-values.

### Signature

```r
compareMEs(
  MEs.template,
  MEs.target,
  template.name,
  target.name,
  Grouping.template,
  Grouping.target,
  sex.template = NULL,
  PMI.template = NULL,
  sex.target = NULL,
  PMI.target = NULL
)
```

### Arguments

* `MEs.template`, `MEs.target`: data.frame/matrix (rows = samples, cols = modules). Column names may or may not have `"ME"` prefix; `"MEgrey"` is removed automatically
* `template.name`, `target.name`: labels for plot titles and PDF filename
* `Grouping.template`, `Grouping.target`: group labels (factor/character) for each cohort’s samples (length must match rows of respective MEs)
* `sex.template`, `PMI.template`, `sex.target`, `PMI.target` (optional): covariates; if provided, model is `MEs ~ Group + Sex (+ PMI)`, else `MEs ~ Group`

### Returns

* `list(template = p_template, target = p_target)` — named numeric vectors of ANOVA p-values per ME for each cohort

### Side effects

* Creates a PDF named:
  `"{template.name}(template)-{target.name} Synthetic Eigenproteins.pdf"`
* For each ME:

  * **Left panel:** template cohort boxplot+beeswarm; title red if p < 0.05
  * **Right panel:** target cohort boxplot+beeswarm; title red if p < 0.05

### Notes & Behavior

* Uses an **MLM** (`lm(data.matrix(MEs) ~ ...)`) and extracts the **model F-test p-value** for each response (ME)
* `Sex` and `PMI` are **optional**; if omitted, they are **not** included in the model
* Grouping levels are ordered with simple **priority rules** to put likely controls first (e.g., “Control”, “WT”)

### Minimal Example

```r
pvals <- compareMEs(
  MEs.template = MEs.UPenn,
  MEs.target   = MEsSynth.ROSMAP,
  template.name = "UPenn.vasc",
  target.name   = "ROSMAP400.vasc",
  Grouping.template = numericMeta.UPenn$Group,
  Grouping.target   = numericMeta.ROSMAP$Dx,
  sex.template = numericMeta.UPenn$msex,
  PMI.template = numericMeta.UPenn$PMI
)
```

---

## 3) `kMEcleanup()`

### Purpose

Iteratively **repair** module color assignments (`NETcolors`) using **kME intramodule strength** and the **gap to the global kME max**. Designed to reconcile noisy or partial initial assignments (e.g., when mapping synthetic ME members) and converge to a stable network.

### Signature

```r
kMEcleanup(
  cleanDat,
  MEs,
  NETcolors = NULL,
  synthMEmembersCSV = NULL,
  kMEmaxDiff = 0.10,
  reassignIfGT = 0.30,
  greyIfLT = 0.30,
  keepColors = TRUE
)
```

### Arguments

* `cleanDat` (matrix/df): expression (rows = features, cols = samples)
* `MEs` (matrix/df): eigengenes (rows = samples, cols = modules). Column names may be `"ME<color>"` or `"<color>"`. A `grey` column (if present) is dropped internally for kME calculations
* `NETcolors` (named character vector, optional): starting color per **row** of `cleanDat`. If **missing**, it is **derived** from `synthMEmembersCSV`
* `synthMEmembersCSV` (string, optional): if `NETcolors` missing, CSV with **one column per color** and rows listing member assay IDs (rownames in `cleanDat`) used to build each synthetic ME
* `kMEmaxDiff` (default 0.10): if (kME.max - kME.intra) >= this value -> **reassign** to module of kME.max
* `reassignIfGT` (default 0.30): **grey** features with kME.max > threshold are reassigned to the kME.max module
* `greyIfLT` (default 0.30): features with **intra-module kME < threshold** are set to **grey**
* `keepColors` (logical, default **TRUE**): if TRUE, **preserve color identities** and protect **rolling hubs**

  * **Rolling hub protection:** at each iteration, the **top kME member** (hub) in each color is **pinned** to remain in that color for the next iteration; hubs are **excluded** from reassign/greying and **excluded** from convergence checks

### How it works (per iteration)

1. **Optionally enforce** last iteration’s hubs (when `keepColors=TRUE`)
2. Grey out low intra-module kME (< `greyIfLT`)
3. Reassign from grey to max-kME color if kME.max > `reassignIfGT`
4. Reassign non-grey rows to their max-kME color if (kME.max - kME.intra) >= `kMEmaxDiff`
5. **Recompute** eigengenes and kME on the updated colors (messages suppressed)
6. **Update** rolling hubs (if `keepColors=TRUE`)
7. **Converge** when:

   * All **non-grey**, **non-hub** rows have intra-module kME ≥ `greyIfLT` **and**
   * Their (kME.max - kME.intra) ≤ `kMEmaxDiff`
     (Hubs are intentionally excluded from the criteria when `keepColors=TRUE`.)

Hard stop after **30 iterations** with a warning.

### Returns

`list(kMEtable, NETcolors, MEs)` where:

* `kMEtable`: per-feature table with color, all kMEs, and a sortable key (module rank + intramodule kME or grey average)
* `NETcolors`: **final** color assignment per feature
* `MEs`: **final, unprefixed**, non-grey module eigengenes (rows = samples, cols = colors)

### Side effects

* Console messages: iteration status, **“clean in X iterations”** or **“max iterations reached”**, module sizes, grey counts (before & after)

### Minimal Example

```r
res <- kMEcleanup(
  cleanDat = cleanDat.ROSMAP,
  MEs = MEs.initial,
  # either provide NETcolors directly…
  NETcolors = initialColors,
  # …or build from a CSV of synthetic ME members:
  # synthMEmembersCSV = "Synthetic_eigengene_exact_match_members.csv",
  kMEmaxDiff = 0.10,
  reassignIfGT = 0.30,
  greyIfLT = 0.30,
  keepColors = TRUE
)

# Final outputs:
res$NETcolors   # cleaned colors
res$MEs         # final eigengenes (non-grey)
res$kMEtable    # detailed per-feature kME table (for export)
```

### Tips & Troubleshooting

* If convergence **stalls** with `keepColors=TRUE`, hubs are **intentionally pinned**. This is expected; hubs are excluded from the stop criteria, so the loop will still end once non-hub members stabilize.
* If you prefer classic WGCNA **relabeling by size** (e.g., turquoise, blue, brown ordering), set `keepColors=FALSE`.
* Ensure `rownames(cleanDat)` are **exactly** the assay IDs used in the synthetic member CSV (if you use it); any missing IDs are reported as warnings.
* For very small modules or sparse overlap, consider lowering `greyIfLT` or `kMEmaxDiff`, or raising `reassignIfGT`.

---

## Quick Start: End-to-End

```r
# 1) Build synthetic MEs in the target cohort
MEsSynth <- getSynthMEs(
  cleanDat.template = cleanDatMEGA,
  netColors.template = netMEGA$colors,
  cleanDat.target = cleanDat,
  matchCSV = "Synthetic_eigengene_exact_match_members.csv"
)

# 2) Compare template vs target (optionally with covariates)
pvals <- compareMEs(
  MEs.template = WGCNA::orderMEs(WGCNA::moduleEigengenes(t(cleanDatMEGA), netMEGA$colors)$eigengenes),
  MEs.target   = MEsSynth,
  template.name = "UPenn.vasc",
  target.name   = "ROSMAP400.vasc",
  Grouping.template = numericMeta.UPenn$Group,
  Grouping.target   = numericMeta.ROSMAP$Dx
)

# 3) Clean/repair module color assignments (optional)
res <- kMEcleanup(
  cleanDat = cleanDat,
  MEs = MEsSynth,
  synthMEmembersCSV = "Synthetic_eigengene_exact_match_members.csv",
  keepColors = TRUE
)
```

---

## Reproducibility Notes

* All three functions **coerce inputs to matrices/data.frames** and sanitize ME column naming (`"ME"` prefix handled; `"grey"` dropped where appropriate).
* Console messaging is designed to **trace decisions** and termination conditions.
* `kMEcleanup()` suppresses WGCNA chatter during eigengene recomputation to keep logs readable.

---

## License

MIT.

---

## Citation

If these functions contribute to your work, please cite **WGCNA** and link to this GitHub repository; this synthetic ME workflow has been used in the following publications:

1.  Dai J, Johnson ECB, Dammer EB, Duong DM, Gearing M, Lah JJ, Levey AI, Wingo TS, Seyfried NT: Effects of APOE Genotype on Brain Proteomic Network and Cell Type Changes in Alzheimer's Disease. Frontiers in Molecular Neuroscience 2018, Volume 11 - 2018.
2.  Dammer EB, Ping L, Duong DM, Modeste ES, Seyfried NT, Lah JJ, Levey AI, Johnson ECB: Multi-platform proteomic analysis of Alzheimer’s disease cerebrospinal fluid and plasma reveals network biomarkers associated with proteostasis and the matrisome. Alzheimer's Research & Therapy 2022, 14(1):174.
3.  Dammer EB, Shantaraman A, Ping L, Duong DM, Gerasimov ES, Ravindran SP, Gudmundsdottir V, Frick EA, Gomez GT, Walker KA et al: Proteomic analysis of Alzheimer's disease cerebrospinal fluid reveals alterations associated with APOE e4 and atomoxetine treatment. Science Translational Medicine 2024, 16(753):eadn3504.
4.  Gutierrez-Quiceno L, Dammer EB, Johnson AG, Webster JA, Shah R, Duong D, Yin L, Seyfried NT, Alvarez VE, Stein TD et al: A proteomic network approach resolves stage-specific molecular phenotypes in chronic traumatic encephalopathy. Molecular Neurodegeneration 2021, 16(1):40.
5.  Johnson ECB, Dammer EB, Duong DM, Ping L, Zhou M, Yin L, Higginbotham LA, Guajardo A, White B, Troncoso JC et al: Large-scale proteomic analysis of Alzheimer’s disease brain and cerebrospinal fluid reveals early changes in energy metabolism associated with microglia and astrocyte activation. Nature Medicine 2020, 26(5):769-780.
6.  Levites Y, Dammer EB, Ran Y, Tsering W, Duong D, Abreha M, Gadhavi J, Lolo K, Trejo-Lopez J, Phillips J et al: Integrative proteomics identifies a conserved Beta amyloid responsome, novel plaque proteins, and pathology modifiers in Alzheimer's disease. Cell Reports Medicine 2024, 5(8).
7.  Saloner R, Staffaroni AM, Dammer EB, Johnson ECB, Paolillo EW, Wise A, Heuer HW, Forsberg LK, Lario-Lago A, Webb JD et al: Large-scale network analysis of the cerebrospinal fluid proteome identifies molecular signatures of frontotemporal lobar degeneration. Nature Aging 2025, 5(6):1143-1158.
8.  Shantaraman A, Dammer EB, Ugochukwu O, Duong DM, Yin L, Carter EK, Gearing M, Chen-Plotkin A, Lee EB, Trojanowski JQ et al: Network proteomics of the Lewy body dementia brain reveals presynaptic signatures distinct from Alzheimer’s disease. Molecular Neurodegeneration 2024, 19(1):60.
9.  Trautwig AN, Fox EJ, Dammer EB, Shantaraman A, Ping L, Duong DM, Watson CM, Wu F, Asress S, Guo Q et al: Network analysis of the cerebrospinal fluid proteome reveals shared and unique differences between sporadic and familial forms of amyotrophic lateral sclerosis. Molecular Neurodegeneration 2025, 20(1):58.


