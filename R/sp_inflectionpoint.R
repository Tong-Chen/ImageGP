# 计算拐点, 代码取自ROSE
numPts_below_line <- function(myVector, slope, x) {
  yPt <- myVector[x]
  b <- yPt - (slope * x)
  xPts <- 1:length(myVector)
  return(sum(myVector <= (xPts * slope + b)))
}

sp_inflectionpoint <- function (data,
                                which_col,
                                slope = NULL,
                                color_plot = "red",
                                color_abline = 8,
                                color_point = 2,
                                type = "l",
                                lower = 1,
                                lty = 2,
                                saveplot = NULL,
                                keep_point = "greater",
                                ...)
{

  if (class(data) == "character") {
    enhancer = sp_readTable(data, header = F)
  } else {
    enhancer <- data
  }

  # head(enhancer)

  H3K27ac = sort(enhancer[, which_col])


  if (!is.null(saveplot)) {
    base_plot_save(saveplot, ...)
  }

  plot(H3K27ac, col = color_plot, type = type)


  inputVector <- H3K27ac
  #set those regions with more control than ranking equal to zero
  inputVector[inputVector < 0] <- 0

  # This is the slope of the line we want to slide. This is the diagonal.
  if (sp.is.null(slope)) {
    slope <- (max(inputVector) - min(inputVector)) / length(inputVector)
  }
  # Find the x-axis point where a line passing through that point has the minimum number
  # of points below it. (ie. tangent)。
  # 该点就是切点
  xPt <- floor(
    optimize(
      numPts_below_line,
      lower = 1,
      upper = length(inputVector),
      myVector = inputVector,
      slope = slope
    )[[1]]
  )

  y_cutoff <-
    inputVector[xPt] #The y-value at this x point. This is our cutoff.

  b <- y_cutoff - (slope * xPt)
  abline(v = xPt,
         h = y_cutoff,
         lty = lty,
         col = color_abline)
  points(xPt,
         y_cutoff,
         pch = 16,
         cex = 0.9,
         col = color_point)
  abline(coef = c(b, slope), col = 2)
  title(
    paste(
      "x=",
      xPt,
      "\ny=",
      signif(y_cutoff, 3),
      "\nFold over Median=",
      signif(y_cutoff / median(inputVector), 3),
      "x\nFold over Mean=",
      signif(y_cutoff / mean(inputVector), 3),
      "x",
      sep = ""
    )
  )

  #Number of regions with zero signal
  axis(4,
       sum(inputVector == 0),
       sum(inputVector == 0),
       col.axis = "blue",
       col = "blue")



  ## 超级增强子cluster
  if (keep_point == "greater") {
    keeppoint <- enhancer[enhancer[, which_col] >= y_cutoff, ]
  } else if (keep_point == "little"){
    keeppoint <- enhancer[enhancer[, which_col] < y_cutoff, ]
  }
  write.table(
    keeppoint,
    file = "keep_point.xls",
    sep = "\t",
    quote = F,
    row.names = F
  )

  if (!is.null(saveplot)) {
    dev.off()
  }
}
