# User parameters ----

n_cores <- 11
n_cores <- if (Sys.info()[["sysname"]] == "Windows") {
  1
} else {
  min(n_cores, parallel::detectCores())
}

misc_dir <- "misc"
plot_files <- Sys.glob(file.path(misc_dir, "plot_*.bin"))

PlotDdaFile <- function(plot_file, with_feature) {
  prefix <- if (with_feature) "spec_w_feat_" else "spec_wo_feat_"
  plot_id <- substr(
    basename(substr(plot_file, 1, nchar(plot_file) - 4)),
    6,
    999
  )
  pdf_file <- paste0(prefix, plot_id, ".pdf")

  pdf(pdf_file, paper = "a4r", width = 0, height = 0)
  par(mar = c(4, 4, 3, 1))
  layout(matrix(c(1, 1, 2), 3, 1))

  print(pdf_file)
  data_file <- file(plot_file, "rb")

  while (TRUE) {
    spectrum_name <- readBin(data_file, character())
    if (length(spectrum_name) == 0) {
      break
    }

    feature_mz <- readBin(data_file, numeric(), size = 4, endian = "little")
    shape <- readBin(data_file, numeric(), size = 4, endian = "little")
    smooth <- readBin(data_file, numeric(), size = 4, endian = "little")
    peak_rt <- readBin(data_file, numeric(), size = 4, endian = "little")
    spectrum_mz <- readBin(data_file, numeric(), size = 4, endian = "little")
    spectrum_rt <- readBin(data_file, numeric(), size = 4, endian = "little")

    rt_left <- readBin(
      data_file,
      integer(),
      size = 2,
      signed = FALSE,
      endian = "little"
    )
    rt_right <- readBin(
      data_file,
      integer(),
      size = 2,
      signed = FALSE,
      endian = "little"
    )

    xic_length <- readBin(
      data_file,
      integer(),
      size = 2,
      signed = FALSE,
      endian = "little"
    )
    xic_matrix <- matrix(
      readBin(
        data_file,
        numeric(),
        size = 4,
        n = xic_length * 2,
        endian = "little"
      ),
      nrow = 2
    )
    rt_xic <- xic_matrix[1, ]
    intensity_xic <- xic_matrix[2, ]

    library_mass <- readBin(data_file, numeric(), size = 4, endian = "little")
    dot_product <- readBin(data_file, numeric(), size = 4, endian = "little")

    peak_length <- readBin(
      data_file,
      integer(),
      size = 1,
      signed = FALSE,
      endian = "little"
    )
    experimental_matrix <- matrix(
      readBin(
        data_file,
        numeric(),
        size = 4,
        n = peak_length * 2,
        endian = "little"
      ),
      nrow = 2
    )
    experimental_mz <- experimental_matrix[1, ]
    experimental_intensity <- experimental_matrix[2, ]

    collision_energy <- readBin(
      data_file,
      numeric(),
      size = 4,
      endian = "little"
    )

    peak_length <- readBin(
      data_file,
      integer(),
      size = 1,
      signed = FALSE,
      endian = "little"
    )
    library_matrix <- matrix(
      readBin(
        data_file,
        numeric(),
        size = 4,
        n = peak_length * 2,
        endian = "little"
      ),
      nrow = 2
    )
    library_mz <- library_matrix[1, ]
    library_intensity <- library_matrix[2, ]

    peak_length <- readBin(
      data_file,
      integer(),
      size = 1,
      signed = FALSE,
      endian = "little"
    )
    matched_library_matrix <- matrix(
      readBin(
        data_file,
        numeric(),
        size = 4,
        n = peak_length * 2,
        endian = "little"
      ),
      nrow = 2
    )
    matched_library_mz <- matched_library_matrix[1, ]
    matched_library_intensity <- matched_library_matrix[2, ]

    matched_experimental_matrix <- matrix(
      readBin(
        data_file,
        numeric(),
        size = 4,
        n = peak_length * 2,
        endian = "little"
      ),
      nrow = 2
    )
    matched_experimental_mz <- matched_experimental_matrix[1, ]
    matched_experimental_intensity <- matched_experimental_matrix[2, ]

    if ((with_feature && peak_rt == 0) || (!with_feature && peak_rt > 0)) {
      next
    }
    if (dot_product <= 0) {
      next
    }

    y_limit <- max(experimental_intensity)
    feature_or_spectrum_mz <- if (with_feature) feature_mz else spectrum_mz
    feature_title <- if (with_feature) {
      paste0(
        "MS feature: (",
        round(feature_mz, 3),
        "m/z, ",
        round(peak_rt, 2),
        "min), "
      )
    } else {
      ""
    }
    collision_energy_title <- if (collision_energy > 0) {
      paste0(", CE=", round(collision_energy, 2))
    } else {
      ""
    }
    plot_title <- paste0(
      feature_title,
      "MSMS: (",
      round(spectrum_mz, 3),
      "m/z, ",
      round(spectrum_rt, 2),
      "min)",
      collision_energy_title,
      ", mass_diff(lib-MS)=",
      format(
        round(library_mass - feature_or_spectrum_mz, digits = 4),
        scientific = FALSE
      ),
      ", score=",
      round(dot_product, 2),
      "\n",
      spectrum_name
    )

    plot(
      x = NA,
      type = "n",
      ylim = c(-y_limit, y_limit),
      xlim = c(0, feature_or_spectrum_mz),
      main = plot_title,
      xlab = "m/z",
      ylab = "normalised intensity"
    )

    abline(h = 0, col = gray(0, 0.2))

    max_library_intensity <- 0
    for (i in seq_along(library_mz)) {
      if (library_mz[i] < spectrum_mz - 0.1) {
        max_library_intensity <- library_intensity[i]
        break
      }
    }
    if (max_library_intensity <= 0) {
      max_library_intensity <- max(library_intensity)
    }

    max_library_intensity <- max_library_intensity / y_limit

    # Library spectrum, downward.
    library_intensity <- -library_intensity / max_library_intensity
    lines(x = library_mz, y = library_intensity, col = "blue", type = "h")
    text(
      library_mz[1:20],
      library_intensity[1:20],
      labels = round(library_mz[1:20], 2),
      cex = 0.7,
      col = gray(0, 0.2)
    )
    text(
      matched_library_mz,
      -matched_library_intensity / max_library_intensity,
      labels = round(matched_library_mz, 2),
      cex = 0.7,
      col = "blue"
    )

    # Experimental spectrum, upward.
    lines(
      x = experimental_mz,
      y = experimental_intensity,
      col = "red",
      type = "h"
    )
    text(
      experimental_mz[1:20],
      experimental_intensity[1:20],
      labels = round(experimental_mz[1:20], 2),
      cex = 0.7,
      col = gray(0, 0.2)
    )
    text(
      matched_experimental_mz,
      matched_experimental_intensity,
      labels = round(matched_experimental_mz, 2),
      cex = 0.7,
      col = "red"
    )

    legend(
      "topleft",
      legend = c("experimental", "library"),
      col = c("red", "blue"),
      lty = c(1, 1)
    )

    abline(v = spectrum_mz, col = gray(0, 0.1))

    # EIC plot.
    eic_title <- if (with_feature) {
      paste0("shape=", round(shape, 2), ", SN=", round(smooth, 2))
    } else {
      ""
    }
    plot(
      x = rt_xic,
      y = intensity_xic,
      pch = ".",
      cex = 3,
      xlab = "RT (minutes)",
      ylab = "intensity",
      main = eic_title
    )

    rect(
      rt_xic[rt_left],
      par("usr")[3],
      rt_xic[rt_right],
      par("usr")[4],
      col = gray(0, 0.1),
      border = NA
    )

    polygon(
      x = c(rt_xic[rt_left:rt_right], rt_xic[rt_right], rt_xic[rt_left]),
      y = c(intensity_xic[rt_left:rt_right], 0, 0),
      col = gray(0, 0.2),
      border = NA
    )

    abline(v = spectrum_rt, col = "red")
    abline(v = peak_rt, col = "red", lty = 2, lwd = 2)

    legend(
      "topleft",
      legend = c("MSMS", "apex"),
      lty = c(1, 2),
      lwd = c(1, 2),
      col = c("red", "red")
    )
  }

  dev.off()
  close(data_file)

  return(invisible(NULL))
}

results <- parallel::mclapply(
  plot_files,
  function(plot_file) {
    return(PlotDdaFile(plot_file, TRUE))
  },
  mc.cores = n_cores
)
# results <- parallel::mclapply(
#   plot_files,
#   function(plot_file) {
#     return(PlotDdaFile(plot_file, FALSE))
#   },
#   mc.cores = n_cores
# )
