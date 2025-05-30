#' Plot scib single methods - RNA
#' @description Plotting scib-metrics result of the integration for each single RNA task in an overview table.
#' Integration methods are ranked from top to bottom based on an 'Overall score', calculated as weighted
#' sum of 'Batch correction' (default w:0.4) and 'Bio conservation' (default w:0.6).
#'
#' @return \code{plotSingleTaskRNA} saves in `outdir` one overview table for each task in three formats (.pdf/.tiff/.png)
#' and a .csv file for each task containing the ranked summary table scores.
#'
#' @param csv_metrics_path path to a .csv file output of scib that contains the metrics calculated
#' across one or multiple RNA tasks.
#' @param outdir output directory where the plots and .csv will be saved.
#' @param weight_batch number in [0,1] to use as weight for the batch correction metrics. Weight for
#' bio conservation is calculated as 1-weight_batch
#' @examples
#' plotSingleTaskRNA(csv_metrics_path = "./data/metrics_RNA_allTasks.csv")
#' plotSingleTaskRNA(csv_metrics_path = "./data/metrics_RNA_immune_human.csv")
#'
#'

library(tibble)
library(RColorBrewer)
library(dynutils)
library(stringr)
library(Hmisc)
library(plyr)


library(dplyr)
library(scales)
library(ggimage)
library(cowplot)

scIB_knit_table <- function(
  data,
  column_info,
  row_info,
  palettes,
  usability = FALSE,
  atac = FALSE,
  atac_best = FALSE
) {
  # no point in making these into parameters
  row_height <- 1.1
  row_space <- .1
  row_bigspace <- .5
  col_width <- 1.1
  col_space <- .2
  col_bigspace <- .5
  segment_data <- NULL

  # DETERMINE ROW POSITIONS
  if (!"group" %in% colnames(row_info) || all(is.na(row_info$group))) {
    row_info$group <- ""
    row_groups <- tibble(group = "")
    plot_row_annotation <- FALSE
  } else {
    plot_row_annotation <- TRUE
  }

  row_pos <-
    row_info %>%
    group_by(group) %>%
    mutate(group_i = row_number()) %>%
    ungroup() %>%
    mutate(
      row_i = row_number(),
      colour_background = group_i %% 2 == 1,
      do_spacing = c(FALSE, diff(as.integer(factor(group))) != 0),
      ysep = ifelse(do_spacing, row_height + 2 * row_space, row_space),
      y = - (row_i * row_height + cumsum(ysep)),
      ymin = y - row_height / 2,
      ymax = y + row_height / 2
    )

  # DETERMINE COLUMN POSITIONS
  if (!"group" %in% colnames(column_info) || all(is.na(column_info$group))) {
    column_info$group <- ""
    plot_column_annotation <- FALSE
  } else {
    plot_column_annotation <- TRUE
  }

  column_info <-
    column_info %>%
    add_column_if_missing(width = col_width, overlay = FALSE)


  column_pos <-
    column_info %>%
    mutate(
      do_spacing = c(FALSE, diff(as.integer(factor(group))) != 0),
      xsep = case_when(
        overlay ~ c(0, -head(width, -1)),
        do_spacing ~ col_bigspace,
        TRUE ~ col_space
      ),
      xwidth = case_when(
        overlay & width < 0 ~ width - xsep,
        overlay ~ -xsep,
        TRUE ~ width
      ),
      xmax = cumsum(xwidth + xsep),
      xmin = xmax - xwidth,
      x = xmin + xwidth / 2
    )


  ##########################
  #### CREATE GEOM DATA ####
  ##########################

  # gather circle data
  ind_circle <- which(column_info$geom == "circle")
  if(length(ind_circle) > 0){
    dat_mat <- as.matrix(data[, ind_circle])
    col_palette <- data.frame(metric = colnames(dat_mat),
                              group = column_info[match(colnames(dat_mat), column_info$id), "group"])

    col_palette$name_palette <- lapply(col_palette$group, function(x) palettes[[as.character(x)]])

    circle_data <- data.frame(label = unlist(lapply(colnames(dat_mat),
                                                    function(x) rep(x, nrow(dat_mat)))),
                              x0 = unlist(lapply(column_pos$x[ind_circle],
                                                 function(x) rep(x, nrow(dat_mat)))),
                              y0 = rep(row_pos$y, ncol(dat_mat)),
                              r = row_height/2*as.vector(sqrt(dat_mat))
                              )
    for(l in unique(circle_data$label)){
      ind_l <- which(circle_data$label == l)
      circle_data[ind_l, "r"] <- rescale(circle_data[ind_l, "r"], to = c(0.05, 0.55), from = range(circle_data[ind_l, "r"], na.rm = T))
    }

    colors <- NULL


    for(i in 1:ncol(dat_mat)){
      palette <- colorRampPalette(rev(brewer.pal(9, col_palette$name_palette[[i]])))(nrow(data)-sum(is.na(dat_mat[,i])))
      colors <- c(colors, palette[rank(dat_mat[,i], ties.method = "average", na.last = "keep")])
    }

    circle_data$colors <- colors
  }


  # gather bar data
  ind_bar <- which(column_info$geom == "bar")
  dat_mat <- as.matrix(data[, ind_bar])

  col_palette <- data.frame(metric = colnames(dat_mat),
                            group = column_info[match(colnames(dat_mat), column_info$id), "group"])

  col_palette$name_palette <- lapply(col_palette$group, function(x) palettes[[as.character(x)]])


  rect_data <- data.frame(label = unlist(lapply(colnames(dat_mat),
                                                function(x) rep(x, nrow(dat_mat)))),
                          method = rep(row_info$id, ncol(dat_mat)),
                          value = as.vector(dat_mat),
                          xmin = unlist(lapply(column_pos[ind_bar, "xmin"],
                                               function(x) rep(x, nrow(dat_mat)))),
                          xmax = unlist(lapply(column_pos[ind_bar, "xmax"],
                                               function(x) rep(x, nrow(dat_mat)))),
                          ymin = rep(row_pos$ymin, ncol(dat_mat)),
                          ymax = rep(row_pos$ymax, ncol(dat_mat)),
                          xwidth = unlist(lapply(column_pos[ind_bar, "xwidth"],
                                                 function(x) rep(x, nrow(dat_mat))))
  )
  rect_data <- rect_data %>%
    add_column_if_missing(hjust = 0) %>%
    mutate(
      xmin = xmin + (1 - value) * xwidth * hjust,
      xmax = xmax - (1 - value) * xwidth * (1 - hjust)
    )

  colors <- NULL
  for(i in 1:ncol(dat_mat)){
    palette <- colorRampPalette(rev(brewer.pal(9, col_palette$name_palette[[i]])))(nrow(data)-sum(is.na(dat_mat[,i])))
    colors <- c(colors, palette[rank(dat_mat[,i], ties.method = "average", na.last = "keep")])
  }

  rect_data$colors <- colors



  # gather text data
  ind_text <- which(column_info$geom == "text")
  dat_mat <- as.matrix(data[, ind_text])

  if(atac){
    colnames(dat_mat)[1] <- "Method"
  }

  text_data <- data.frame(label_value = as.vector(dat_mat),
                          group = rep(colnames(dat_mat), each = nrow(dat_mat)),
                          xmin = unlist(lapply(column_pos[ind_text, "xmin"],
                                               function(x) rep(x, nrow(dat_mat)))),
                          xmax = unlist(lapply(column_pos[ind_text, "xmax"],
                                              function(x) rep(x, nrow(dat_mat)))),
                          ymin = rep(row_pos$ymin, ncol(dat_mat)),
                          ymax = rep(row_pos$ymax, ncol(dat_mat)),
                          size = 4, fontface = "plain", stringsAsFactors = F)

  text_data$colors <- "black"
  text_data[text_data$label_value == "HVG", "colors"] <- "darkgreen"
  text_data[text_data$label_value == "FULL", "colors"] <- "grey30"

  # replace scaled/unscaled with +/-
  text_data$label_value <- mapvalues(text_data$label_value, from = c("yes", "no"),
                                     to = c("+", "-"))

  text_data[text_data$label_value == "+" | text_data$label_value == "-", "size"] <- 5
  text_data[text_data$label_value == "+" | text_data$label_value == "-", "fontface"] <- "bold"

  text_data[text_data$label_value == "genes" | text_data$label_value == "peaks" | text_data$label_value == "windows", "size"] <- 3


  # ADD top3 ranking for each bar column
  if(usability || atac_best){
    cols_bar <- unique(rect_data$label)
    cols_bar <- as.character(cols_bar[!is.na(cols_bar)])
    for(c in cols_bar){
      rect_tmp <- rect_data[rect_data$label == c,]
      rect_tmp <- add_column(rect_tmp, "label_value" = as.character(rank(-rect_tmp$value, ties.method = "min")))
      rect_tmp <- rect_tmp[rect_tmp$label_value %in% c("1", "2", "3"), c("label_value", "xmin", "xmax", "ymin", "ymax")]
      rect_tmp <- add_column(rect_tmp, "size" = 2.5, .after = "ymax")
      rect_tmp <- add_column(rect_tmp, "colors" = "black", .after = "size")
      rect_tmp <- add_column(rect_tmp, "fontface" = "plain", .after = "colors")
      rect_tmp <- add_column(rect_tmp, "group" = "top3", .after = "fontface")
      text_data <- bind_rows(text_data, rect_tmp)
    }
  }



  # ADD COLUMN NAMES
  df <- column_pos %>% filter(id != "Method") %>% filter(id != "Ranking")

  if (nrow(df) > 0) {
    segment_data <- segment_data %>% bind_rows(
      df %>% transmute(x = x, xend = x, y = -.3, yend = -.1, size = .5)
    )
    text_data <-
      bind_rows(
        text_data,
        df %>% transmute(
          xmin = x, xmax = x, ymin = 0, ymax = -0.5,
          angle = 30, vjust = 0, hjust = 0,
          label_value = id,
          size = 3
        )
      )
  }


  # GENERATE ROW ANNOTATION
  if (plot_row_annotation) {
    row_annotation <-
      row_pos %>%
      select(group, ymin, ymax) %>%
      group_by(group) %>%
      summarise(
        ymin = min(ymin),
        ymax = max(ymax),
        y = (ymin + ymax) / 2
      ) %>%
      ungroup() %>%
      mutate(xmin = -.5, xmax = 5) %>%
      filter(!is.na(group), group != "")

    text_data <- text_data %>% bind_rows(
      row_annotation %>%
        transmute(xmin, xmax, ymin = ymax + row_space, label_value = group %>% gsub("\n", " ", .),
                  hjust = 0, vjust = .5, fontface = "bold", size = 4) %>%
        mutate(ymax = ymin + row_height)
    )
  }

  # gather image data
  ind_img <- which(column_info$geom == "image")
  if(length(ind_img) > 0){
    dat_mat <- as.matrix(data[, ind_img])

    image_data <- data.frame(x = unlist(lapply(column_pos$x[ind_img],
                                                function(x) rep(x, nrow(dat_mat)))),
                             y = rep(row_pos$y, ncol(dat_mat)),
                             image = mapvalues(dat_mat, from = c("graph", "embed", "gene"),
                                               to = c("./img/graph.png", "./img/embedding.png", "./img/matrix.png")),
                             stringsAsFactors = FALSE
                             )

  }

  suppressWarnings({
    minimum_x <- min(column_pos$xmin, segment_data$x, segment_data$xend,
                     text_data$xmin, na.rm = TRUE)
    maximum_x <- max(column_pos$xmax, segment_data$x, segment_data$xend,
                     text_data$xmax, na.rm = TRUE)
    minimum_y <- min(row_pos$ymin, segment_data$y, segment_data$yend,
                     text_data$ymin, na.rm = TRUE)
    maximum_y <- max(row_pos$ymax, segment_data$y, segment_data$yend,
                     text_data$ymax, na.rm = TRUE)
  })

  ####################################
  ###   CREATE HARDCODED LEGENDS   ###
  ####################################

  x_min_output <- minimum_x+0.5
  x_min_scaling <- minimum_x + 5.5
  x_min_ranking <- ifelse(atac, minimum_x + 5.5, minimum_x + 10.5)
  x_min_score <-  ifelse(atac, minimum_x + 11, minimum_x + 17)

  leg_max_y <- minimum_y - .5

  # Create legend for Output
  leg_min_x <- x_min_output
  output_title_data <- data.frame(xmin = leg_min_x,
                                xmax = leg_min_x+ 2,
                                ymin = leg_max_y - 1,
                                ymax = leg_max_y,
                                label_value = "Output",
                                hjust = 0, vjust = 0,
                                fontface = "bold",
                                size = 3)

  output_img <- data.frame(x = leg_min_x+0.5,
                           y = c(leg_max_y-2, leg_max_y-3.2,leg_max_y-4.4),
                           image = c("./img/matrix.png", "./img/embedding.png", "./img/graph.png")
                           )
  if(atac || atac_best){
    output_text <- data.frame(xmin = leg_min_x+1.5,
                              xmax = leg_min_x+3,
                              ymin = c(leg_max_y-2.2, leg_max_y-3.4,leg_max_y-4.6),
                              ymax = c(leg_max_y-2.2, leg_max_y-3.4,leg_max_y-4.6),
                              label_value = c("feature", "embed", "graph"),
                              hjust = 0, vjust = 0,
                              fontface = "plain",
                              size = 3)
  } else{
    output_text <- data.frame(xmin = leg_min_x+1.5,
                            xmax = leg_min_x+3,
                            ymin = c(leg_max_y-2.2, leg_max_y-3.4,leg_max_y-4.6),
                            ymax = c(leg_max_y-2.2, leg_max_y-3.4,leg_max_y-4.6),
                            label_value = c("gene", "embed", "graph"),
                            hjust = 0, vjust = 0,
                            fontface = "plain",
                            size = 3)
  }

  text_data <- bind_rows(text_data, output_text, output_title_data)
  image_data <- bind_rows(image_data, output_img)

  # Create legend for scaling
  if(!atac && !atac_best){
  leg_min_x <- x_min_scaling
  scaling_title_data <- data.frame(xmin = leg_min_x,
                                  xmax = leg_min_x+ 2,
                                  ymin = leg_max_y - 1,
                                  ymax = leg_max_y,
                                  label_value = "For spatial",
                                  hjust = 0, vjust = 0,
                                  fontface = "bold",
                                  size = 3)

  scaling_text <- data.frame(xmin = c(leg_min_x, leg_min_x+1),
                            xmax = c(leg_min_x+0.5, leg_min_x+3),
                            ymin = c(rep(leg_max_y-2,2), rep(leg_max_y-3,2)),
                            ymax = c(rep(leg_max_y-1,2), rep(leg_max_y-2,2)),
                            label_value = c("+", ": yes", "-", ": no"),
                            hjust = 0, vjust = 0,
                            fontface = c("bold","plain", "bold", "plain"),
                            size = c(5,3,5,3))

  text_data <- bind_rows(text_data, scaling_title_data, scaling_text)
  }

  # CREATE LEGEND for ranking colors
  leg_min_x <- x_min_ranking
  rank_groups <- as.character(column_info[column_info$geom == "bar", "group"])

  if(usability){
    rank_minimum_x <- list("RNA" = leg_min_x,
                           "Simulation" = leg_min_x+1,
                           "Usability" = leg_min_x+2,
                           "Scalability" = leg_min_x+3)
    leg_max_x <- leg_min_x+3
  } else if(atac_best){
    rank_minimum_x <- list("ATAC_windows" = leg_min_x,
                           "ATAC_peaks" = leg_min_x+1,
                           "ATAC_genes" = leg_min_x+2)
    leg_max_x <- leg_min_x+2
  } else{
    rank_minimum_x <- list("Score overall" = leg_min_x,
                           "Removal of batch effects" = leg_min_x+1,
                           "Cell type label variance" = leg_min_x+2)
    leg_max_x <- leg_min_x+2
  }

  rank_title_data <- data.frame(xmin = leg_min_x,
                                xmax = leg_min_x+ 2,
                                ymin = leg_max_y - 1,
                                ymax = leg_max_y,
                                label_value = "Ranking",
                                hjust = 0, vjust = 0,
                                fontface = "bold")

  for(rg in rank_groups){
    rank_palette <- colorRampPalette(rev(brewer.pal(9, palettes[[rg]])))(5)


    rank_data <- data.frame(xmin = rank_minimum_x[[rg]],
                            xmax = rank_minimum_x[[rg]] + .8,
                            ymin = seq(leg_max_y-4, leg_max_y - 2, by = .5),
                            ymax = seq(leg_max_y-3.5, leg_max_y -1.5, by = .5),
                            border = TRUE,
                            colors = rank_palette
    )
    rect_data <- bind_rows(rect_data, rank_data)

  }

  # create arrow for ranking
  arrow_data <- data.frame(x = leg_max_x + 1.5,
                           xend = leg_max_x +1.5,
                           y = leg_max_y-4,
                           yend = leg_max_y -1.5)


  # add text next to the arrow
  arrow_text <- data.frame(xmin = leg_max_x +2,
                           xmax = leg_max_x +2.5,
                           ymin = c(leg_max_y-2, leg_max_y-4),
                           ymax = c(leg_max_y-1.5, leg_max_y-3.5 ),
                           label_value = c("1", as.character(nrow(data))),
                           hjust = 0, vjust = 0, size = 2.5)


  text_data <- bind_rows(text_data, rank_title_data, arrow_text)

  # CREATE LEGEND for circle scores
  # circle legend
  if(!usability && !atac_best){
    cir_minimum_x <- x_min_score

    cir_legend_size <- 1
    cir_legend_space <- .1

    cir_legend_dat <-
      data.frame(
        value = seq(0, 1, by = .2),
        r = row_height/2*seq(0, 1, by = .2)
      )
    cir_legend_dat$r <- rescale(cir_legend_dat$r, to = c(0.05, 0.55), from = range(cir_legend_dat$r, na.rm = T))

    x0 <- vector("integer", nrow(cir_legend_dat))
    for(i in 1:length(x0)){
      if(i == 1){
        x0[i] <- cir_minimum_x + cir_legend_space + cir_legend_dat$r[i]
      }
      else {
        x0[i] <- x0[i-1] + cir_legend_dat$r[i-1] + cir_legend_space + cir_legend_dat$r[i]
      }
    }

    cir_legend_dat$x0 <- x0
    cir_legend_min_y <- leg_max_y-4
    cir_legend_dat$y0 <- cir_legend_min_y + 1 + cir_legend_dat$r

    cir_legend_dat$colors <- NULL
    cir_maximum_x <- max(cir_legend_dat$x0)

    cir_title_data <- data_frame(xmin = cir_minimum_x,
                                 xmax = cir_maximum_x,
                                 ymin = leg_max_y -1,
                                 ymax = leg_max_y,
                                 label_value = "Score",
                                 hjust = 0, vjust = 0, fontface = "bold")

    cir_value_data <- data.frame(xmin = cir_legend_dat$x0 - cir_legend_dat$r,
                                 xmax = cir_legend_dat$x0 + cir_legend_dat$r,
                                 ymin = cir_legend_min_y,
                                 ymax = cir_legend_min_y +3,
                                 hjust = .5, vjust = 0, size = 2.5,
                                 label_value = ifelse(cir_legend_dat$value %in% c(0, 1),
                                                      paste0(cir_legend_dat$value*100, "%"), ""))

    circle_data <- bind_rows(circle_data, cir_legend_dat)
    text_data <- bind_rows(text_data, cir_title_data, cir_value_data)



  }

  minimum_y <- min(minimum_y, min(text_data$ymin, na.rm = TRUE))

  ########################
  ##### COMPOSE PLOT #####
  ########################

  g <-
    ggplot() +
    coord_equal(expand = FALSE) +
    scale_alpha_identity() +
    scale_colour_identity() +
    scale_fill_identity() +
    scale_size_identity() +
    scale_linetype_identity() +
    cowplot::theme_nothing()

  # PLOT ROW BACKGROUNDS
  df <- row_pos %>% filter(colour_background)
  if (nrow(df) > 0) {
    g <- g + geom_rect(aes(xmin = min(column_pos$xmin)-.25, xmax = max(column_pos$xmax)+.25, ymin = ymin - (row_space / 2), ymax = ymax + (row_space / 2)), df, fill = "#DDDDDD")
  }



  # PLOT CIRCLES
  if (length(ind_circle) > 0) {
    g <- g + ggforce::geom_circle(aes(x0 = x0, y0 = y0, fill= colors, r = r), circle_data, size=.25)
  }


  # PLOT RECTANGLES
  if (nrow(rect_data) > 0) {
    # add defaults for optional values
    rect_data <- rect_data %>%
      add_column_if_missing(alpha = 1, border = TRUE, border_colour = "black") %>%
      mutate(border_colour = ifelse(border, border_colour, NA))

    g <- g + geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = colors, colour = border_colour, alpha = alpha), rect_data, size = .25)
  }


  # PLOT TEXT
  if (nrow(text_data) > 0) {
    # add defaults for optional values
    text_data <- text_data %>%
      add_column_if_missing(
        hjust = .5,
        vjust = .5,
        size = 3,
        fontface = "plain",
        colors = "black",
        lineheight = 1,
        angle = 0
      ) %>%
      mutate(
        angle2 = angle / 360 * 2 * pi,
        cosa = cos(angle2) %>% round(2),
        sina = sin(angle2) %>% round(2),
        alphax = ifelse(cosa < 0, 1 - hjust, hjust) * abs(cosa) + ifelse(sina > 0, 1 - vjust, vjust) * abs(sina),
        alphay = ifelse(sina < 0, 1 - hjust, hjust) * abs(sina) + ifelse(cosa < 0, 1 - vjust, vjust) * abs(cosa),
        x = (1 - alphax) * xmin + alphax * xmax,
        y = (1 - alphay) * ymin + alphay * ymax
      ) %>%
      filter(label_value != "")
    # Set fontface for legend bold
    text_data[text_data$label_value == "Ranking", "fontface"] <- "bold"
    # Set fontface for ranking numbers bold
    if(usability || atac_best){
    text_data[1:nrow(data), "fontface"] <- "bold"
    }
    # subset text_data to left-aligned rows
    text_data_left <- text_data[which(text_data$group == "Method" | text_data$group == "top3"), ]
    text_data <- text_data[-which(text_data$group == "Method" | text_data$group == "top3"), ]

    g <- g + geom_text(aes(x = x, y = y, label = label_value, colour = colors, hjust = hjust, vjust = vjust, size = size, fontface = fontface, angle = angle), data = text_data)

    text_data_left[text_data_left$group == "Method", "x"] <- text_data_left[text_data_left$group == "Method", "x"] - 3
    if(usability || atac_best){
    text_data_left[text_data_left$group == "top3", "x"] <- text_data_left[text_data_left$group == "top3", "xmin"] + .3
    text_data_left[text_data_left$group == "Method", "x"] <- text_data_left[text_data_left$group == "Method", "x"] + .5
    }
    g <- g + geom_text(aes(x = x, y = y, label = label_value, colour = colors, hjust = "left", vjust = vjust, size = size, fontface = fontface, angle = angle), data = text_data_left)
  }



  # PLOT SEGMENTS
  if (nrow(segment_data) > 0) {
    # add defaults for optional values
    segment_data <- segment_data %>% add_column_if_missing(size = .5, colour = "black", linetype = "solid")

    g <- g + geom_segment(aes(x = x, xend = xend, y = y, yend = yend, size = size, colour = colour, linetype = linetype), segment_data)
  }

  # PLOT ARROW RANKING
  if (nrow(arrow_data) > 0) {
    # add defaults for optional values
    arrow_data <- arrow_data %>% add_column_if_missing(size = .5, colour = "black", linetype = "solid")

    g <- g + geom_segment(aes(x = x, xend = xend, y = y, yend = yend, size = size, colour = colour, linetype = linetype), arrow_data, arrow = arrow(length = unit(0.1, "cm")), lineend = "round", linejoin = "bevel")
  }

  # PLOT IMAGES
  if(length(ind_img) > 0){
    for(r in 1:nrow(image_data)){
      g <- g + cowplot::draw_image(image = image_data$image[r], x = image_data[r, "x"]-.5, y = image_data[r, "y"]-.5)
    }

  }

  # ADD SIZE
  # reserve a bit more room for text that wants to go outside the frame
  minimum_x <- minimum_x - 2
  maximum_x <- maximum_x + 5
  minimum_y <- minimum_y - 2
  maximum_y <- maximum_y + 4

  g$width <- maximum_x - minimum_x
  g$height <- maximum_y - minimum_y

  g <- g + expand_limits(x = c(minimum_x, maximum_x), y = c(minimum_y, maximum_y))






  return(g)

}



add_column_if_missing <- function(df, ...) {
  column_values <- list(...)
  for (column_name in names(column_values)) {
    default_val <- rep(column_values[[column_name]], nrow(df))

    if (column_name %in% colnames(df)) {
      df[[column_name]] <- ifelse(is.na(df[[column_name]]), default_val, df[[column_name]])
    } else {
      df[[column_name]] <- default_val
    }
  }
  df
}


plotSingleTaskBatch <- function(csv_metrics_path, outdir = ".", weight_batch = 0.5){

  metrics_tab_lab <- read.csv(csv_metrics_path, sep = ",")

  # get metrics names from columns
  metrics <- colnames(metrics_tab_lab)[-1]
  metrics <- gsub("\\.", "/", metrics)
  metrics <- gsub("_", " ", metrics)
  metrics <- plyr::mapvalues(metrics, from = c("ASW label/batch", "graph conn", "iLISI"),
                             to = c("Batch ASW","graph connectivity", "graph iLISI"))

  print(metrics)
  # metrics names as they are supposed to be ordered
  group_batch <- c('Batch ASW spatial', 'graph iLISI spatial', 'graph connectivity spatial', 'kBET spatial')
  group_bio <- c("Batch ASW omics", "graph iLISI omics", "graph connectivity omics", "kBET omics")
  # set original values of number of metrics
  n_metrics_batch_original <- sum(group_batch %in% metrics)
  n_metrics_bio_original <- sum(group_bio %in% metrics)

  # order metrics present in the table
  matching.order <- match(c(group_batch, group_bio), metrics)
  metrics.ord <- metrics[matching.order[!is.na(matching.order)]]

  # get methods info from rownames
  methods_info_full  <- as.character(metrics_tab_lab[,1])

  # in case methods names start with /
  if(substring(methods_info_full[1], 1, 1) == "/"){
  methods_info_full <- sub("/", "", methods_info_full)
  }

  # Remove trvae full
  ind.trvae_full <- grep("trvae_full", methods_info_full)
  if(length(ind.trvae_full) > 0){
    methods_info_full <- methods_info_full[-ind.trvae_full]
    metrics_tab_lab <- metrics_tab_lab[-ind.trvae_full,]
  }

  # data scenarios to be saved in file name
  data.scenarios <- unique(unlist(sapply(str_split(methods_info_full, "/"), function(x) x[1])))




  ###### Plot one figure for each data task
  for (dt.sc in data.scenarios){
    ind.scen <- grep(paste0(dt.sc, "/"), methods_info_full)
    methods_info <- methods_info_full[ind.scen]
    metrics_tab_sub <- metrics_tab_lab[ind.scen, ]

    # info on scaling data
    scaling <- sapply(str_split(methods_info, "/"), function(x) x[3])


    # info on HVG selection
    hvg <- sapply(str_split(methods_info, "/"), function(x) x[4])
    hvg <- plyr::mapvalues(hvg, from = c("hvg", "full_feature"), to = c("HVG", "FULL"))

    methods <- sapply(str_split(methods_info, "/"), function(x) x[5])

    methods_name <- sapply(str_split(methods, "_"), function(x) x[1])
    # methods_name <- capitalize(methods_name)
    methods_name <- plyr::mapvalues(methods_name,
                                    from = c("Seurat", "Seuratrpca", "Mnn", "Bbknn", "Trvae", "scVI", "Liger", "Combat", "Saucie", "Fastmnn", "Desc", "Scanvi", "Scgen"),
                                    to = c("Seurat v3 CCA", "Seurat v3 RPCA", "MNN", "BBKNN", "trVAE", "scVI", "LIGER", "ComBat", "SAUCIE", "fastMNN", "DESC", "scANVI*", "scGen*"))


    method_groups <- sapply(str_split(methods, "_"), function(x) x[2])
    method_groups <- plyr::mapvalues(method_groups,
                                     from = c("knn", "embed", "full"),
                                     to = c("graph", "embed", "gene"))



    ##### Create dataframe
    metrics_tab <- as.data.frame(metrics_tab_sub[, -1])
    metrics_tab[metrics_tab == ""] <- NA
    colnames(metrics_tab) <- metrics

    #add Methods column
    metrics_tab <- add_column(metrics_tab, "Method" = methods_name, .before = 1)

    # reorder columns by metrics
    col.ordered <- c("Method", metrics.ord)
    metrics_tab <- metrics_tab[, col.ordered]

    ## Remove columns that are full NAs
    na.col <- apply(metrics_tab, 2, function(x) sum(is.na(x)) == nrow(metrics_tab))
    # redefine numbers of metrics per group
    if(sum(colnames(metrics_tab)[na.col] %in%  group_batch) > 0){
      n_metrics_batch <- n_metrics_batch_original - sum(colnames(metrics_tab)[na.col] %in%  group_batch)
    } else {
      n_metrics_batch <- n_metrics_batch_original
    }

    if(sum(colnames(metrics_tab)[na.col] %in% group_bio) > 0){
      n_metrics_bio <- n_metrics_bio_original - sum(colnames(metrics_tab)[na.col] %in% group_bio)
    } else{
      n_metrics_bio <- n_metrics_bio_original
    }

    metrics_tab <- metrics_tab[, !na.col]


    ## Scores should be already scaled [0,1] - however, we aim to compute the scores based on the min-max scaled metrics
    scaled_metrics_tab <- as.matrix(metrics_tab[, -1])
    scaled_metrics_tab <- apply(scaled_metrics_tab, 2, function(x) scale_minmax(x))

    # calculate average score by group and overall
    score_group_batch <- rowMeans(scaled_metrics_tab[, 1:n_metrics_batch], na.rm = T)
    score_group_bio <- rowMeans(scaled_metrics_tab[, (1+n_metrics_batch):ncol(scaled_metrics_tab)],
                             na.rm = T)

    score_all <- (weight_batch*score_group_batch + (1-weight_batch)*score_group_bio)

    metrics_tab <- add_column(metrics_tab, "Overall score" = score_all, .after = "Method")
    metrics_tab <- add_column(metrics_tab, "Spatial batch correction" = score_group_batch, .after = "Overall score")
    metrics_tab <- add_column(metrics_tab, "Omics batch correction" = score_group_bio, .after = 3+n_metrics_batch)

    metrics_tab <- add_column(metrics_tab, "Output" = method_groups, .after = "Method")
    metrics_tab <- add_column(metrics_tab, "Features" = hvg, .after = "Output")
    metrics_tab <- add_column(metrics_tab, "For spatial" = scaling, .after = "Features")

    # order methods by the Overall score
    metrics_tab <- metrics_tab[order(metrics_tab$`Overall score`,  decreasing = T), ]
    write.csv(metrics_tab, file = paste0(outdir, "/", dt.sc, "_summary_scores.csv"), quote = F)


    # Delete rows that are empty
    rowsNA <- which(is.na(metrics_tab$`Overall score`))
    if(length(rowsNA) >0){
      metrics_tab <- metrics_tab[-rowsNA, ]
    }


    # Defining column_info, row_info and palettes
    row_info <- data.frame(id = metrics_tab$Method)

    column_info <- data.frame(id = colnames(metrics_tab),
                              group = c("Text", "Image", "Text", "Text", "Score overall",
                                        rep("Removal of batch effects", (1 + n_metrics_batch)),
                                        rep("Cell type label variance", (1 + n_metrics_bio))),
                              geom = c("text", "image", "text", "text", "bar", "bar",
                                       rep("circle", n_metrics_batch), "bar", rep("circle", n_metrics_bio)),
                              width = c(8,2.5,2,1.5,2,2, rep(1,n_metrics_batch), 2, rep(1,n_metrics_bio)),
                              overlay = F)

    # defining colors palette
    palettes <- list("Score overall" = "YlGnBu",
                     "Removal of batch effects" = "BuPu",
                     "Cell type label variance" = "RdPu")


    g <- scIB_knit_table(data = metrics_tab, column_info = column_info, row_info = row_info, palettes = palettes, usability = F)
    now <- Sys.time()
    ggsave(paste0(outdir, "/", format(now, "%Y%m%d_%H%M%S_"), dt.sc, "_summary_metrics.pdf"), g, device = cairo_pdf, width = 297, height = 420, units = "mm")
    # ggsave(paste0(outdir, "/", format(now, "%Y%m%d_%H%M%S_"), dt.sc, "_summary_metrics.tiff"), g, device = "tiff", dpi = "retina", width = 297, height = 420, units = "mm")
    ggsave(paste0(outdir, "/", format(now, "%Y%m%d_%H%M%S_"), dt.sc, "_summary_metrics.png"), g, device = "png", dpi = "retina", width = 297, height = 420, units = "mm")
  }
}



plotSingleTaskCluster <- function(csv_metrics_path, outdir = ".", weight_batch = 0.5){

  metrics_tab_lab <- read.csv(csv_metrics_path, sep = ",")

  # get metrics names from columns
  metrics <- colnames(metrics_tab_lab)[-1]
  # metrics <- gsub("\\.", "/", metrics)
  metrics <- gsub("_", " ", metrics)
  # metrics <- plyr::mapvalues(metrics, from = c("ASW label/batch", "graph conn", "iLISI"),
  #                            to = c("Batch ASW","graph connectivity", "graph iLISI"))

  print(metrics)
  # metrics names as they are supposed to be ordered
  group_batch <- c("Simluation dataset ARI spatial" ,"Simluation dataset NMI spatial", "MERFISH dataset ARI spatial" ,   "MERFISH dataset NMI spatial"  ,"Xenium dataset ARI spatial","Xenium dataset NMI spatial")
  group_bio <- c("Simluation dataset ARI omics",  "Simluation dataset NMI omics"  , "MERFISH dataset ARI omics" , "MERFISH dataset NMI omics" , "Xenium dataset ARI omics"  , "Xenium dataset NMI omics"  )
  # group_batch <- c('nbr ari', 'nbr nmi')
  # group_bio <- c('gex ari', 'gex nmi')

  # set original values of number of metrics
  n_metrics_batch_original <- sum(group_batch %in% metrics)
  n_metrics_bio_original <- sum(group_bio %in% metrics)
  print(n_metrics_batch_original)
  print(n_metrics_bio_original)

  # order metrics present in the table
  matching.order <- match(c(group_batch, group_bio), metrics)
  metrics.ord <- metrics[matching.order[!is.na(matching.order)]]

  # get methods info from rownames
  methods_info_full  <- as.character(metrics_tab_lab[,1])

  # in case methods names start with /
  if(substring(methods_info_full[1], 1, 1) == "/"){
  methods_info_full <- sub("/", "", methods_info_full)
  }

  # Remove trvae full
  ind.trvae_full <- grep("trvae_full", methods_info_full)
  if(length(ind.trvae_full) > 0){
    methods_info_full <- methods_info_full[-ind.trvae_full]
    metrics_tab_lab <- metrics_tab_lab[-ind.trvae_full,]
  }

  # data scenarios to be saved in file name
  data.scenarios <- unique(unlist(sapply(str_split(methods_info_full, "/"), function(x) x[1])))


  ###### Plot one figure for each data task
  for (dt.sc in data.scenarios){
    ind.scen <- grep(paste0(dt.sc, "/"), methods_info_full)
    methods_info <- methods_info_full[ind.scen]
    metrics_tab_sub <- metrics_tab_lab[ind.scen, ]

    # info on scaling data
    scaling <- sapply(str_split(methods_info, "/"), function(x) x[3])


    # info on HVG selection
    hvg <- sapply(str_split(methods_info, "/"), function(x) x[4])
    hvg <- plyr::mapvalues(hvg, from = c("hvg", "full_feature"), to = c("HVG", "FULL"))

    methods <- sapply(str_split(methods_info, "/"), function(x) x[5])

    methods_name <- sapply(str_split(methods, "_"), function(x) x[1])
    # methods_name <- capitalize(methods_name)
    methods_name <- plyr::mapvalues(methods_name,
                                    from = c("Seurat", "Seuratrpca", "Mnn", "Bbknn", "Trvae", "scVI", "Liger", "Combat", "Saucie", "Fastmnn", "Desc", "Scanvi", "Scgen"),
                                    to = c("Seurat v3 CCA", "Seurat v3 RPCA", "MNN", "BBKNN", "trVAE", "scVI", "LIGER", "ComBat", "SAUCIE", "fastMNN", "DESC", "scANVI*", "scGen*"))


    method_groups <- sapply(str_split(methods, "_"), function(x) x[2])
    method_groups <- plyr::mapvalues(method_groups,
                                     from = c("knn", "embed", "full"),
                                     to = c("graph", "embed", "gene"))

    ##### Create dataframe
    metrics_tab <- as.data.frame(metrics_tab_sub[, -1])
    metrics_tab[metrics_tab == ""] <- NA
    colnames(metrics_tab) <- metrics

    #add Methods column
    metrics_tab <- add_column(metrics_tab, "Method" = methods_name, .before = 1)

    # reorder columns by metrics
    col.ordered <- c("Method", metrics.ord)
    metrics_tab <- metrics_tab[, col.ordered]

    ## Remove columns that are full NAs
    na.col <- apply(metrics_tab, 2, function(x) sum(is.na(x)) == nrow(metrics_tab))
    # redefine numbers of metrics per group
    if(sum(colnames(metrics_tab)[na.col] %in%  group_batch) > 0){
      n_metrics_batch <- n_metrics_batch_original - sum(colnames(metrics_tab)[na.col] %in%  group_batch)
    } else {
      n_metrics_batch <- n_metrics_batch_original
    }

    if(sum(colnames(metrics_tab)[na.col] %in% group_bio) > 0){
      n_metrics_bio <- n_metrics_bio_original - sum(colnames(metrics_tab)[na.col] %in% group_bio)
    } else{
      n_metrics_bio <- n_metrics_bio_original
    }

    metrics_tab <- metrics_tab[, !na.col]


    ## Scores should be already scaled [0,1] - however, we aim to compute the scores based on the min-max scaled metrics
    scaled_metrics_tab <- as.matrix(metrics_tab[, -1])
    scaled_metrics_tab <- apply(scaled_metrics_tab, 2, function(x) scale_minmax(x))

    # calculate average score by group and overall
    score_group_batch <- rowMeans(scaled_metrics_tab[, 1:n_metrics_batch], na.rm = T)
    score_group_bio <- rowMeans(scaled_metrics_tab[, (1+n_metrics_batch):ncol(scaled_metrics_tab)],
                             na.rm = T)

    score_all <- (weight_batch*score_group_batch + (1-weight_batch)*score_group_bio)

    metrics_tab <- add_column(metrics_tab, "Overall score" = score_all, .after = "Method")
    metrics_tab <- add_column(metrics_tab, "Spatial metrics" = score_group_batch, .after = "Overall score")
    metrics_tab <- add_column(metrics_tab, "Omics metrics" = score_group_bio, .after = 3+n_metrics_batch)

    metrics_tab <- add_column(metrics_tab, "Output" = method_groups, .after = "Method")
    metrics_tab <- add_column(metrics_tab, "Features" = hvg, .after = "Output")
    metrics_tab <- add_column(metrics_tab, "For spatial" = scaling, .after = "Features")

    # order methods by the Overall score
    metrics_tab <- metrics_tab[order(metrics_tab$`Overall score`,  decreasing = T), ]
    write.csv(metrics_tab, file = paste0(outdir, "/", dt.sc, "_summary_scores.csv"), quote = F)


    # Delete rows that are empty
    rowsNA <- which(is.na(metrics_tab$`Overall score`))
    if(length(rowsNA) >0){
      metrics_tab <- metrics_tab[-rowsNA, ]
    }


    # Defining column_info, row_info and palettes
    row_info <- data.frame(id = metrics_tab$Method)

    column_info <- data.frame(id = colnames(metrics_tab),
                              group = c("Text", "Image", "Text", "Text", "Score overall",
                                        rep("Removal of batch effects", (1 + n_metrics_batch)),
                                        rep("Cell type label variance", (1 + n_metrics_bio))),
                              geom = c("text", "image", "text", "text", "bar", "bar",
                                       rep("circle", n_metrics_batch), "bar", rep("circle", n_metrics_bio)),
                              width = c(8,2.5,2,1.5,2,2, rep(1,n_metrics_batch), 2, rep(1,n_metrics_bio)),
                              overlay = F)

    # defining colors palette
    palettes <- list("Score overall" = "YlGnBu",
                     "Removal of batch effects" = "BuPu",
                     "Cell type label variance" = "RdPu")


    g <- scIB_knit_table(data = metrics_tab, column_info = column_info, row_info = row_info, palettes = palettes, usability = F)
    now <- Sys.time()
    ggsave(paste0(outdir, "/", format(now, "%Y%m%d_%H%M%S_"), dt.sc, "_summary_metrics.pdf"), g, device = cairo_pdf, width = 297, height = 420, units = "mm")
    # ggsave(paste0(outdir, "/", format(now, "%Y%m%d_%H%M%S_"), dt.sc, "_summary_metrics.tiff"), g, device = "tiff", dpi = "retina", width = 297, height = 420, units = "mm")
    ggsave(paste0(outdir, "/", format(now, "%Y%m%d_%H%M%S_"), dt.sc, "_summary_metrics.png"), g, device = "png", dpi = "retina", width = 297, height = 420, units = "mm")

  }
}
