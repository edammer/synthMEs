#' Build synthetic module eigengenes (MEs) in a target cohort
#' using a template cohort's modules and kME ranking.
#'
#' Requires WGCNA loaded.
#'
#' @param minimumMEmembers integer (default 4)
#' @param topPercent numeric in (0,1] (default 0.20)
#' @param minKmeThreshToRescue numeric in [0,1] (default 0.70)
#' @param cleanDat.template matrix/data.frame of template expression (rows = features, cols = samples)
#' @param netColors.template vector of module colors (length = nrow(cleanDat.template))
#' @param cleanDat.target matrix/data.frame of target expression (rows = features, cols = samples)
#' @param matchCSV output CSV file name for matched markers (default
#'        "Synthetic_eigengene_exact_match_members.csv")
#' @return data.frame of synthetic MEs for the target cohort

getSynthMEs <- function(
  minimumMEmembers = 4,
  topPercent        = 0.20,
  minKmeThreshToRescue = 0.70,
  cleanDat.template,
  netColors.template,
  cleanDat.target,
  matchCSV = "Synthetic_eigengene_exact_match_members.csv"
) {
  # --- Basic checks -----------------------------------------------------------
  if (!requireNamespace("WGCNA", quietly = TRUE)) {
    stop("Package 'WGCNA' is required. Please install/load it.")
  }
  if (!is.matrix(cleanDat.template) && !is.data.frame(cleanDat.template)) {
    stop("cleanDat.template must be a matrix or data.frame with rows=features, cols=samples.")
  }
  if (!is.matrix(cleanDat.target) && !is.data.frame(cleanDat.target)) {
    stop("cleanDat.target must be a matrix or data.frame with rows=features, cols=samples.")
  }
  if (length(netColors.template) != nrow(cleanDat.template)) {
    stop("Length of netColors.template must equal nrow(cleanDat.template).")
  }
  if (topPercent <= 0 || topPercent > 1) {
    stop("topPercent must be in (0,1].")
  }
  if (minimumMEmembers < 1) {
    stop("minimumMEmembers must be >= 1.")
  }

  # Coerce to matrix to avoid data.frame pitfalls
  cleanDat.template <- as.matrix(cleanDat.template)
  cleanDat.target   <- as.matrix(cleanDat.target)

  # --- Template MEs (consensus cohort) ---------------------------------------
  MEList <- WGCNA::moduleEigengenes(t(cleanDat.template), colors = netColors.template)
  MEs.template <- WGCNA::orderMEs(MEList$eigengenes)                # has "ME" prefix
  rownames(MEs.template) <- colnames(cleanDat.template)              # samples as rows

  # Drop grey from template ME sets for downstream matching
  keepCols         <- setdiff(colnames(MEs.template), "MEgrey")
  MEs.template.nog <- MEs.template[, keepCols, drop = FALSE]

  # Prepare ME matrix (with/without "ME" prefix) for kME calculation
  tmpMEs <- MEs.template
  # ensure 'MEgrey' removed in the prefixed version as well
  tmpMEs <- tmpMEs[, keepCols, drop = FALSE]

  # --- kME in the template to rank members within each module -----------------
  kMEdat <- WGCNA::signedKME(t(cleanDat.template), tmpMEs, corFnc = "bicor")

  # --- Build overlap lists of target features per module ----------------------
  targetSymbols <- rownames(cleanDat.target)
  allColors     <- unique(netColors.template)
  colorlist     <- setdiff(allColors, "grey")

  # For each module color: rank members by kME, take top topPercent (or rescue by threshold),
  # then intersect with target features; keep only if >= minimumMEmembers.
  thismod <- thismodTop <- overlap <- vector("list", length(colorlist))
  names(thismod) <- names(thismodTop) <- names(overlap) <- colorlist

  for (eachColor in colorlist) {
    idx <- which(netColors.template == eachColor)
    # data.frame of (Symbol, kME) sorted by kME desc
    kmeCol <- paste0("kME", eachColor)
    df <- data.frame(
      Symbols = rownames(cleanDat.template)[idx],
      kME     = kMEdat[idx, kmeCol],
      row.names = NULL
    )
    df <- df[order(df$kME, decreasing = TRUE), ]
    thismod[[eachColor]] <- df

    # take top topPercent (at least 1 row if any exist)
    topN <- max(1, round(nrow(df) * topPercent, 0))
    topDF <- df[seq_len(topN), , drop = FALSE]

    # primary overlap
    overlap1 <- intersect(targetSymbols, topDF$Symbols)

    # if not enough, attempt "rescue" using kME threshold
    if (length(overlap1) < minimumMEmembers) {
      below <- which(df$kME < minKmeThreshToRescue)
      if (length(below) > 0) {
        take <- seq_len(max(1, below[1] - 1))
      } else {
        take <- seq_len(nrow(df))  # no entries below threshold -> take all
      }
      topDF <- df[take, , drop = FALSE]
      overlap1 <- intersect(targetSymbols, topDF$Symbols)
    }

    # keep only if we meet the minimum members criterion
    if (length(overlap1) >= minimumMEmembers) {
      overlap[[eachColor]]   <- overlap1
      thismodTop[[eachColor]] <- topDF
    } else {
      overlap[[eachColor]] <- NULL
      thismodTop[[eachColor]] <- NULL
    }
  }

  # Remove modules with no overlap
  overlap <- overlap[!vapply(overlap, is.null, logical(1))]
  if (length(overlap) == 0) {
    warning("No modules met the overlap criteria in the target cohort; returning empty MEsSynth.")
    # Still create a zero-column data.frame with rownames as target samples
    MEsSynth <- data.frame(row.names = colnames(cleanDat.target))
    # Print the requested dimension info
    message("dim(MEsSynth): ", paste(dim(MEsSynth), collapse = " x "))
    message("dim(MEs)     : ", paste(dim(MEs.template.nog), collapse = " x "))
    message("No modules generated syntheticMEs.")
    # Write an empty CSV with proper header
    utils::write.csv(data.frame(), file = matchCSV, row.names = FALSE)
    return(MEsSynth)
  }

  # --- Assemble marker expression matrix and parallel color vector ------------
  cleanDatMarkers <- matrix(nrow = 0, ncol = ncol(cleanDat.target))
  colnames(cleanDatMarkers) <- colnames(cleanDat.target)
  markerColors <- character(0)

  for (mod in names(overlap)) {
    rowsToBind <- cleanDat.target[match(overlap[[mod]], targetSymbols), , drop = FALSE]
    cleanDatMarkers <- rbind(cleanDatMarkers, rowsToBind)
    markerColors    <- c(markerColors, rep(mod, nrow(rowsToBind)))
  }

  # --- Synthetic MEs in the target cohort ------------------------------------
  MEsSynth <- WGCNA::moduleEigengenes(t(cleanDatMarkers), markerColors, impute = FALSE)$eigengenes
  rownames(MEsSynth) <- colnames(cleanDat.target)

  # Keep only columns that exist in the template ME set (sans grey), in the same order
  keepOrder <- intersect(colnames(MEs.template.nog), colnames(MEsSynth))
  MEsSynth  <- MEsSynth[, keepOrder, drop = FALSE]

  # --- Console outputs as requested ------------------------------------------
  message("dim(MEsSynth): ", paste(dim(MEsSynth), collapse = " x "))
  message("dim(MEs)     : ", paste(dim(MEs.template.nog), collapse = " x "))

  if (ncol(MEsSynth) < ncol(MEs.template.nog)) {
    message("Synthetic MEs missing for ", ncol(MEs.template.nog) - ncol(MEsSynth),
            " template module(s).")
  } else {
    message("All modules generated syntheticMEs.")
  }

  # --- Write CSV of exact-match members used per module -----------------------
  # Create a padded data.frame so all columns have equal length
  maxLen <- max(vapply(overlap, length, integer(1)))
  padCol <- function(x, n) { c(x, rep(NA_character_, n - length(x))) }
  memberDF <- stats::setNames(
    data.frame(lapply(overlap, padCol, n = maxLen), check.names = FALSE),
    names(overlap)
  )
  utils::write.csv(memberDF, file = matchCSV, row.names = FALSE)

  return(MEsSynth)
}



#' Compare template and target MEs with optional covariate adjustment and plots
#'
#' @param MEs.template   data.frame/matrix of template MEs (rows=samples, cols=ME*color)
#' @param MEs.target     data.frame/matrix of target  MEs (rows=samples, cols=ME*color)
#' @param template.name  character label for template (used in PDF filename & titles)
#' @param target.name    character label for target   (used in PDF filename & titles)
#' @param Grouping.template factor/character grouping (length = nrow(MEs.template))
#' @param Grouping.target   factor/character grouping (length = nrow(MEs.target))
#' @param sex.template   OPTIONAL factor/character/NULL (same length as Grouping.template)
#' @param PMI.template   OPTIONAL numeric/NULL (same length as Grouping.template)
#' @param sex.target     OPTIONAL factor/character/NULL (same length as Grouping.target)
#' @param PMI.target     OPTIONAL numeric/NULL (same length as Grouping.target)
#'
#' @return list with two named numeric vectors of p-values:
#'         list(template = p_template, target = p_target)
#'
compareMEs <- function(
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
) {
  # ---- packages ----
  if (!requireNamespace("beeswarm", quietly = TRUE)) {
    stop("Package 'beeswarm' is required for plotting. Install it or remove plotting.")
  }
  if (!requireNamespace("WGCNA", quietly = TRUE)) {
    stop("Package 'WGCNA' is required. Please install/load it.")
  }

  # ---- coerce & align ----
  MEs.template <- as.data.frame(MEs.template, check.names = FALSE)
  MEs.target   <- as.data.frame(MEs.target,   check.names = FALSE)

  if(!grepl('^ME',colnames(MEs.template)[1])) colnames(MEs.template)<-paste0("ME",colnames(MEs.template))
  if(!grepl('^ME',colnames(MEs.target)[1])) colnames(MEs.target)<-paste0("ME",colnames(MEs.target))

  # Drop MEgrey from both sets if present
  keepCols <- setdiff(colnames(MEs.template), "MEgrey")
  MEs.template <- MEs.template[, keepCols, drop = FALSE]
  keepCols <- setdiff(colnames(MEs.target), "MEgrey")
  MEs.target   <- MEs.target[, keepCols, drop = FALSE]

  # ---- helper: build design and run mlm for p-values ----
  .pvals_for_MEs <- function(MEs, Grouping, Sex = NULL, PMI = NULL) {
    # Build covariate data.frame
    df <- data.frame(Group = factor(Grouping), stringsAsFactors = FALSE)

    useSex <- !is.null(Sex)
    usePMI <- !is.null(PMI)

    if (useSex) df$Sex <- factor(Sex)
    if (usePMI) df$PMI <- as.numeric(PMI)

    # formula depending on available covariates
    rhs <- "Group"
    if (useSex) rhs <- paste(rhs, "Sex", sep = " + ")
    if (usePMI) rhs <- paste(rhs, "PMI", sep = " + ")
    fml <- as.formula(paste0("data.matrix(MEs) ~ ", rhs))

    fit <- lm(fml, data = df)                 # mlm
    s   <- summary(fit)

    pvec <- rep(NA_real_, ncol(MEs))
    names(pvec) <- colnames(MEs)

    # summary.mlm returns a list-like object with one summary per response
    for (i in seq_len(ncol(MEs))) {
      fi <- s[[i]]$fstatistic
      # Guard in case of degenerate models
      if (!is.null(fi) && all(is.finite(fi))) {
        pvec[i] <- pf(fi[1], fi[2], fi[3], lower.tail = FALSE)
      }
    }
    pvec
  }

  # ---- p-values for template & target ----
  p_template <- .pvals_for_MEs(MEs.template, Grouping.template, sex.template, PMI.template)
  p_target   <- .pvals_for_MEs(MEs.target,   Grouping.target,   sex.target,   PMI.target)

  # ---- ordered labels table (mirrors your approach, but robust) ----
  # Derive number of modules from target columns
  nModules <- length(colnames(MEs.target))
  orderedLabels <- cbind(
    paste0("M", seq_len(nModules)),
    WGCNA::labels2colors(seq_len(nModules))
  )

  # ---- group ordering with priority patterns ----
  patterns <- c("Control", "CT", "control", "WT", "wt")
  get_priority <- function(x) {
    for (i in seq_along(patterns)) {
      if (grepl(patterns[i], x, ignore.case = TRUE)) return(i)
    }
    Inf
  }

  template.groups <- names(table(Grouping.template))
  template.groups <- template.groups[order(sapply(template.groups, get_priority))]

  target.groups <- names(table(Grouping.target))
  target.groups <- target.groups[order(sapply(target.groups, get_priority))]

  # ---- color helpers ----
  col2hex_safe <- function(colname) {
    # Accept either an R color name or hex
    if (grepl("^#", colname)) {
      sub("^#", "", colname)
    } else {
      rgb <- try(grDevices::col2rgb(colname), silent = TRUE)
      if (inherits(rgb, "try-error")) {
        # Fallback to gray if unknown
        return("777777")
      }
      sprintf("%02X%02X%02X", rgb[1], rgb[2], rgb[3])
    }
  }

  alpha99 <- "99"  # ~60% opacity in hex

  # ---- plotting ----
  pdf(file = paste0(template.name, "(template)-", target.name, " Synthetic Eigenproteins.pdf"),
      width = 8.5, height = 12)
  on.exit(try(dev.off(), silent = TRUE), add = TRUE)

  op <- par(mfrow = c(4, 2), mar = c(5, 6, 4, 2))
  on.exit(try(par(op), silent = TRUE), add = TRUE)

  # To ensure consistent loop order, iterate over target ME columns
  for (i in colnames(MEs.target)) {
    mod <- substring(i, 3, 40)    # strip "ME"
    boxCol <- mod
    transCol <- paste0("#", col2hex_safe(mod), alpha99)

    # Special-case overrides retained from your snippet
    if (mod == "black")        { boxCol <- "#444444BB"; transCol <- "#44444499" }
    if (mod == "midnightblue") { boxCol <- "#191970BB"; transCol <- "#19197099" }
    if (mod == "darkmagenta")  { boxCol <- "#8B008BBB"; transCol <- "#8B008B99" }
    if (mod == "brown4")       { boxCol <- "#8B2323BB"; transCol <- "#8B232399" }
    if (mod == "magenta4")     { boxCol <- "#8B008BBB"; transCol <- "#8B008B99" }

    cexValTemplate <- if (nrow(MEs.template) < 40) 1.75 else 1.0
    cexValTarget   <- if (nrow(MEs.target)   < 40) 1.75 else 1.0

    # --- left panel: TEMPLATE (if this ME exists there) ---
    if (i %in% colnames(MEs.template)) {
      titCol <- ifelse(signif(p_template[i], 2) < 0.05, "red", "black")
      boxplot(MEs.template[, i] ~ factor(Grouping.template, levels = template.groups),
              col = boxCol, ylab = "Eigenprotein Value",
              main = paste0(
                {
                  idx <- match(mod, orderedLabels[, 2]); lab <- ifelse(is.na(idx), "", orderedLabels[idx, 1])
                  paste0(lab, " ", mod, ".", template.name)
                },
                "\nOne-way ANOVA p: ", signif(p_template[i], 2)
              ),
              xlab = NULL, las = 2, outline = FALSE, col.main = titCol)
      beeswarm::beeswarm(MEs.template[, i] ~ factor(Grouping.template, levels = template.groups),
                         method = "swarm", add = TRUE, corralWidth = 0.5, vertical = TRUE,
                         pch = 21, bg = transCol, col = "black", cex = cexValTemplate, corral = "gutter")
    } else {
      frame(); mtext(paste0(i, " not in MEs.template."), side = 3, adj = 0.5, padj = 0.5)
    }

    # --- right panel: TARGET ---
    titCol <- ifelse(signif(p_target[i], 2) < 0.05, "red", "black")
    boxplot(MEs.target[, i] ~ factor(Grouping.target, levels = target.groups),
            col = boxCol, ylab = "Eigenprotein Value",
            main = paste0(i, ".", target.name, " (Synthetic)",
                          "\nOne-way ANOVA p: ", signif(p_target[i], 2)),
            xlab = NULL, las = 2, outline = FALSE, col.main = titCol)
    beeswarm::beeswarm(MEs.target[, i] ~ factor(Grouping.target, levels = target.groups),
                       method = "swarm", add = TRUE, corralWidth = 0.5, vertical = TRUE,
                       pch = 21, bg = transCol, col = "black", cex = cexValTarget, corral = "gutter")
  }

  invisible(list(template = p_template, target = p_target))
}



#' Clean/repair module color assignments using kME rules and (optionally) a members CSV
#'
#' @param cleanDat           matrix/data.frame; rows = features (assays), cols = samples. [REQUIRED]
#' @param MEs                matrix/data.frame of module eigengenes (rows = samples, cols = ME*color). [REQUIRED]
#'                           Column names may be "ME<color>" or just "<color>"; "grey" column (if any) will be dropped.
#' @param NETcolors          optional named character vector of colors for each feature in rownames(cleanDat);
#'                           if not supplied, you must provide synthMEmembersCSV.
#' @param synthMEmembersCSV  optional path to CSV whose columns are module colors and whose column entries
#'                           are feature IDs (rownames of cleanDat) used to build each synthetic ME.
#'                           Used ONLY if NETcolors is not supplied.
#' @param kMEmaxDiff         numeric; reassign if (kME.max - kME.intramodule) >= this value. Default 0.10
#' @param reassignIfGT       numeric; reassign grey items to their kME.max color if kME.max > this value. Default 0.30
#' @param greyIfLT           numeric; set to grey if kME.intramodule < this value. Default 0.30
#' @param keepColors         logical; if TRUE (default), do NOT relabel/reorder colors by size and
#'                           protect each color's initial hub from reassignment.
#'
#' @return list(kMEtable = data.frame, NETcolors = named vector, MEs = data.frame)
#'
kMEcleanup <- function(
  cleanDat,
  MEs,
  NETcolors = NULL,
  synthMEmembersCSV = NULL,
  kMEmaxDiff = 0.10,
  reassignIfGT = 0.30,
  greyIfLT = 0.30,
  keepColors = TRUE
) {
  if (!requireNamespace("WGCNA", quietly = TRUE)) stop("Package 'WGCNA' is required.")
  if (missing(cleanDat) || is.null(cleanDat)) stop("'cleanDat' must be supplied.")
  if (missing(MEs)      || is.null(MEs))      stop("'MEs' must be supplied.")
  if (is.null(NETcolors) && is.null(synthMEmembersCSV)) {
    stop("If 'NETcolors' is not supplied, you must supply 'synthMEmembersCSV'.")
  }

  cleanDat <- as.matrix(cleanDat)
  MEs      <- as.data.frame(MEs, check.names = FALSE)

  # Normalize ME names, drop grey from stored MEs
  if (substr(colnames(MEs)[1], 1, 2) == "ME") colnames(MEs) <- sub("^ME", "", colnames(MEs))
  if ("grey" %in% colnames(MEs)) MEs[["grey"]] <- NULL
  tmpMEs <- MEs
  colnames(tmpMEs) <- paste0("ME", colnames(tmpMEs))
  tmpMEs <- tmpMEs[, !grepl("^MEgrey$", colnames(tmpMEs)), drop = FALSE]

  # Derive NETcolors from CSV if not provided
  if (is.null(NETcolors)) {
    if (!file.exists(synthMEmembersCSV)) stop(paste0("Cannot find file ", synthMEmembersCSV, "!"))
    synthMEmembers <- utils::read.csv(file = synthMEmembersCSV, header = TRUE, check.names = FALSE)
    synthMEcolors  <- colnames(synthMEmembers)

    if (!identical(sort(synthMEcolors), sort(colnames(MEs)))) {
      stop(paste0("Synthetic ME colors in ", synthMEmembersCSV, " do not match ME column colors."))
    }
    NETcolors <- rep("grey", nrow(cleanDat))
    names(NETcolors) <- rownames(cleanDat)
    for (thisCol in synthMEcolors) {
      members <- na.omit(synthMEmembers[[thisCol]])
      present <- intersect(members, rownames(cleanDat))
      NETcolors[present] <- thisCol
      missing <- setdiff(members, rownames(cleanDat))
      if (length(missing)) warning("Missing from cleanDat: ", paste(missing, collapse = ", "))
    }
  } else {
    if (is.null(names(NETcolors))) stop("'NETcolors' must be a NAMED vector matching rownames(cleanDat).")
    if (!all(rownames(cleanDat) %in% names(NETcolors))) {
      stop("All rownames(cleanDat) must be present in names(NETcolors).")
    }
    NETcolors <- NETcolors[rownames(cleanDat)]
  }

  # Initial kME (independent of NETcolors) for decisions
  kMEdat <- WGCNA::signedKME(t(cleanDat), tmpMEs, corFnc = "bicor")

  safe_grey_count <- function(cols) { tab <- table(cols); if (!"grey" %in% names(tab)) 0L else as.integer(tab[["grey"]]) }
  message(sprintf("Grey assay count, before cleanup: %s", safe_grey_count(NETcolors)))
  message(sprintf("  (%.2f%% grey)", 100 * safe_grey_count(NETcolors) / nrow(cleanDat)))

  compute_hubs <- function(NETcolors, kMEdat) {
    hubs <- integer(0)
    cols <- setdiff(unique(NETcolors), "grey")
    for (col in cols) {
      rows <- which(NETcolors == col)
      if (!length(rows)) next
      kcol <- paste0("kME", col)
      if (!kcol %in% colnames(kMEdat)) next
      hubRow <- rows[which.max(kMEdat[rows, kcol])]
      hubs <- c(hubs, hubRow)
    }
    names(hubs) <- setdiff(unique(NETcolors), "grey")[seq_along(hubs)]
    hubs
  }
  protectedIdx <- if (keepColors) compute_hubs(NETcolors, kMEdat) else integer(0)

  retry <- TRUE; iter <- 1L; maxIter <- 30L
  loop_end_reason <- NULL  # "clean", "max_iter"

  while (retry) {
    cat(sprintf("\nkME table Cleanup, processing iteration %d...", iter))

    if (keepColors && length(protectedIdx)) {
      for (col in names(protectedIdx)) {
        idx <- protectedIdx[[col]]
        if (is.finite(idx)) NETcolors[idx] <- col
      }
    }

    kMEintramoduleVector <- vapply(seq_len(nrow(kMEdat)), function(r) {
      col <- NETcolors[r]
      if (is.na(col) || col == "grey") return(NA_real_)
      kcol <- paste0("kME", col)
      if (kcol %in% colnames(kMEdat)) as.numeric(kMEdat[r, kcol]) else NA_real_
    }, numeric(1))

    toGrey <- which(!is.na(kMEintramoduleVector) & kMEintramoduleVector < greyIfLT)
    if (keepColors && length(toGrey) && length(protectedIdx)) toGrey <- setdiff(toGrey, unname(protectedIdx))
    if (length(toGrey)) NETcolors[toGrey] <- "grey"

    kMEmaxVec       <- apply(kMEdat, 1, function(x) max(as.numeric(x), na.rm = TRUE))
    kMEmaxIdx       <- apply(kMEdat, 1, which.max)
    kMEmaxColorsVec <- gsub("^kME", "", colnames(kMEdat)[kMEmaxIdx])

    kMEintramoduleVector <- vapply(seq_len(nrow(kMEdat)), function(r) {
      col <- NETcolors[r]
      if (is.na(col) || col == "grey") return(NA_real_)
      kcol <- paste0("kME", col)
      if (kcol %in% colnames(kMEdat)) as.numeric(kMEdat[r, kcol]) else NA_real_
    }, numeric(1))

    reassign_from_grey <- which(NETcolors == "grey" & is.finite(kMEmaxVec) & (kMEmaxVec > reassignIfGT))
    if (keepColors && length(reassign_from_grey) && length(protectedIdx)) {
      reassign_from_grey <- setdiff(reassign_from_grey, unname(protectedIdx))
    }

    diffs <- kMEmaxVec - ifelse(is.na(kMEintramoduleVector), 1, kMEintramoduleVector)
    diffTooBig <- which(is.finite(diffs) & (diffs >= kMEmaxDiff))
    if (keepColors && length(diffTooBig) && length(protectedIdx)) {
      diffTooBig <- setdiff(diffTooBig, unname(protectedIdx))
    }

    idx_reassign <- union(reassign_from_grey, diffTooBig)
    if (length(idx_reassign)) {
      NETcolors[idx_reassign] <- kMEmaxColorsVec[idx_reassign]
    }

    if (!keepColors) {
      tabNoGrey <- sort(table(NETcolors[NETcolors != "grey"]), decreasing = TRUE)
      if (length(tabNoGrey)) {
        oldColors <- names(tabNoGrey)
        newColors <- WGCNA::labels2colors(seq_along(oldColors))
        map <- stats::setNames(newColors, oldColors)
        NETcolors[NETcolors != "grey"] <- map[NETcolors[NETcolors != "grey"] ]
      }
    }

    # >>> silence eigengene recomputation chatter <<<
    MEList <- suppressMessages({
    tmp <- NULL
    capture.output({
      tmp <- WGCNA::moduleEigengenes(
        t(cleanDat),
        colors = NETcolors,
        verbose = 0
      )
    }, file = NULL)
    tmp
    })
    MEs <- WGCNA::orderMEs(MEList$eigengenes)
    if (substr(colnames(MEs)[1], 1, 2) == "ME") colnames(MEs) <- sub("^ME", "", colnames(MEs))
    if ("grey" %in% colnames(MEs)) MEs[["grey"]] <- NULL

    # tmpMEs = prefixed (no-grey) for signedKME only
    tmpMEs <- MEs
    colnames(tmpMEs) <- paste0("ME", colnames(tmpMEs))
    tmpMEs <- tmpMEs[, !grepl("^MEgrey$", colnames(tmpMEs)), drop = FALSE]
    kMEdat <- WGCNA::signedKME(t(cleanDat), tmpMEs, corFnc = "bicor")

    # rolling hub protection for NEXT iteration (names = colors)
    if (keepColors) {
      compute_hubs <- function(NETcolors, kMEdat) {
        cols <- setdiff(unique(NETcolors), "grey")
        out <- integer(0)
        for (col in cols) {
          rows <- which(NETcolors == col)
          kcol <- paste0("kME", col)
          if (length(rows) && kcol %in% colnames(kMEdat)) {
            out <- c(out, rows[which.max(kMEdat[rows, kcol])])
            names(out)[length(out)] <- col
          }
        }
        out
      }
      protectedIdx <- compute_hubs(NETcolors, kMEdat)
    }

    # ---- STOP CRITERIA (fresh, post-update; exclude hubs if keepColors) ----
    # current intramodule kME (numeric; NA for grey)
    kMEintramodule_current <- vapply(seq_len(nrow(kMEdat)), function(r) {
      col <- NETcolors[r]
      if (is.na(col) || col == "grey") return(NA_real_)
      kcol <- paste0("kME", col)
      if (kcol %in% colnames(kMEdat)) as.numeric(kMEdat[r, kcol]) else NA_real_
    }, numeric(1))

    # max kME per row
    kMEmax_current <- apply(kMEdat, 1, function(x) max(as.numeric(x), na.rm = TRUE))

    # treat GREY rows as satisfied by setting intramodule to 1 (so diff <= 0)
    intramodule_for_check <- ifelse(is.na(kMEintramodule_current), 1, kMEintramodule_current)
    diff_current <- kMEmax_current - intramodule_for_check

    # EXCLUDE protected hubs from convergence checks when keepColors=TRUE
    check_idx <- seq_len(nrow(kMEdat))
    if (keepColors && length(protectedIdx)) {
      check_idx <- setdiff(check_idx, unname(protectedIdx))
    }

    # if all non-protected rows meet both thresholds, we are clean
    if (length(check_idx) == 0 ||
        (min(intramodule_for_check[check_idx], na.rm = TRUE) >= greyIfLT &&
         max(diff_current[check_idx],           na.rm = TRUE) <= kMEmaxDiff)) {
      message(sprintf("\nkME table 'clean' in %d iterations.", iter))
      loop_end_reason <- "clean"
      retry <- FALSE
    }

    iter <- iter + 1L
    if (retry && iter > maxIter) {
      warning("Reached maximum iterations without satisfying cleanup criteria.")
      loop_end_reason <- "max_iter"
      break
    }
  } # while

  # ---- guarantee an exit message even if none printed in-loop ----
  if (is.null(loop_end_reason)) {
    warning("Cleanup loop terminated without recording an exit reason; treating as max-iter.")
    loop_end_reason <- "max_iter"
  }

  # ----- final reporting and table construction (unchanged) -----
  modSizes <- sort(table(NETcolors), decreasing = TRUE)
  nonGreyNames <- setdiff(names(modSizes), "grey")
  nModules <- length(nonGreyNames)

  if (keepColors) {
    orderedModules <- cbind(Mnum = paste0("M", seq_len(nModules)), Color = nonGreyNames)
  } else {
    orderedModules <- cbind(Mnum = paste0("M", seq_len(nModules)),
                            Color = WGCNA::labels2colors(seq_len(nModules)))
  }
  sizesAligned <- as.integer(modSizes[match(orderedModules[, "Color"], names(modSizes))])

  msg <- paste0(
    "\nModules after cleanup in Size Rank order:\n",
    paste(apply(cbind(orderedModules, Size = sizesAligned), 1, function(z) paste(z, collapse = " ")),
          collapse = "\n")
  )
  message(msg)
  message(sprintf("\nGrey assay count, after cleanup: %s", safe_grey_count(NETcolors)))
  message(sprintf("  (%.2f%% grey)", 100 * safe_grey_count(NETcolors) / nrow(cleanDat)))

  orderedModulesWithGrey <- rbind(c("M0", "grey"), orderedModules)
  sortKey <- vapply(seq_len(nrow(cleanDat)), function(r) {
    col <- NETcolors[r]
    if (!is.na(col) && col != "grey") {
      modPair <- orderedModulesWithGrey[match(col, orderedModulesWithGrey[, 2]), ]
      kcol <- paste0("kME", col)
      kval <- if (kcol %in% colnames(kMEdat)) round(as.numeric(kMEdat[r, kcol]), 4) else NA_real_
      paste(paste(modPair, collapse = " "), "|", kval)
    } else {
      paste0("grey|AllKmeAvg:", round(mean(as.numeric(kMEdat[r, ]), na.rm = TRUE), 4))
    }
  }, character(1))

  kMEtable <- cbind(
    Index   = seq_len(nrow(cleanDat)),
    AssayID = rownames(cleanDat),
    Color   = NETcolors,
    kMEdat,
    sortKey
  )
  kMEtable <- kMEtable[order(kMEtable[, "sortKey"], decreasing = TRUE), , drop = FALSE]

  return(list(kMEtable = kMEtable, NETcolors = NETcolors, MEs = tmpMEs))
}
