

#' Generate waterfall plot using R package waterfalls.
#'
#' @param waterfallsinput A data.frame containing two columns,
#' one with the values, the other with the labels.
#' @param labels the labels corresponding to each vector, marked on the x-axis.
#' @param values a numeric vector making up the heights of the rectangles in the waterfall.
#' @param rect_text_labels (character) a character vector of the same length as values that are placed on the rectangles.
#' @param rect_text_size size of the text in the rectangles.
#' @param put_rect_text_outside_when_value_below (numeric) the text labels accompanying a rectangle of this height
#'  will be placed outside the box: below if it's negative; above if it's positive.
#' @param calc_total (logical, default: FALSE) should the final pool of the waterfall be calculated (and placed on the chart).
#' @param draw_lines (logical, default: TRUE) should lines be drawn between successive rectangles.
#' @param fill_colours Colours to be used to fill the rectangles, in order.
#'  Disregarded if fill_by_sign is TRUE (the default).
#' @param lines_anchors a character vector of length two specifying the horizontal placement of the drawn lines relative to the preceding and successive rectangles, respectively.
#' @param linetype the linetype for the draw_lines.
#' @param total_rect_color the color of the final rectangle.
#' @param total_rect_text (character) the text in the middle of the rectangle of the total rectangle.
#' @param total_rect_text_color the color of the final rectangle's label text.
#' @param total_axis_text (character) the text appearing on the axis underneath the total rectangle.
#' @param rect_width (numeric) the width of the rectangle, relative to the space between each label factor.
#' @param draw_axis_x (character) one of "none", "behind", "front" whether to draw an x.axis line and whether to draw it behind or in front of the rectangles, default is behind.
#' @param rect_border the border around each rectangle. Choose NA if no border is desired.
#' @param x_label The X axis name
#' @param y_label The Y axis name
#' @param ...
#'
#' @return pdf image
#' @export
#'
#' @examples
#'
#' waterfallsinput <- "test.file"
#' waterfalls_plot(waterfallsinput)
#'
waterfalls_plot <- function(waterfallsinput, sep="\t", row.names=NULL, header=T,
                            quote="", comment="", check.names=F,labels,values=NULL,
                            rect_text_labels='',rect_text_size=1,
                            put_rect_text_outside_when_value_below=1,
                            calc_total=FALSE, draw_lines=TRUE, fill_colours=NULL,
                            lines_anchors=c("right", "left"),
                            linetype="dashed", total_rect_color="black",
                            total_rect_text_color="white",
                            total_axis_text="Total", total_rect_text,
                            rect_width = 0.9,draw_axis.x="behind",
                            rect_border="white",scale_y_to_waterfall= TRUE,
                            fill_by_sign=FALSE, theme_text_family="",
                            x_label=NULL,y_label=NULL,...){


  data_m <- read.table(waterfallsinput,sep=sep, row.names=row.names, header=header,
                       quote=quote, comment=comment, check.names=check.names)

    rect_text_labels = paste(levels(data_m[,1]),'\n',data_m[,2])

  if (is.null(values)){
    p_waterfalls<-waterfall(.data = data_m,
                            rect_text_labels=rect_text_labels)
    } else {
      p_waterfalls<-waterfall(values=values,
                              rect_text_labels=rect_text_labels)
    }


  if (!is.null(x_label)){
    p_waterfalls <- p_waterfalls + xlab(x_label)
  }

  if (!is.null(y_label)){
    p_waterfalls <- p_waterfalls + ylab(y_label)
  }

  ggsave(paste0(waterfallsinput,"waterfalls.pdf"),plot=p_waterfalls)

}
