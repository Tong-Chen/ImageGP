

# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'
#   Generate DOC:              'Ctrl + Shift + Alt + r'


#' Flower plot could be treated as one kind of venn diagram but only showing common items
#' like OTUs or genes among all groups and total items (or total items excluding common
#' items) for each group.
#'
#' @param input Input data file (first line as header line, the first column is the name of
#' genes or OTUs or otehr things one wants to compare, the second column is the group name which genes belong to,
#' tab seperated)
#'
#' ```
#' Gene	Sample
#' g1	Set1
#' a1	Set3
#' b4	Set1
#' .
#' .
#' c1	Set3
#' ```
#' @param label_total_num_items Label total number of items in for each group (when True) or label number of items in each group
#' after substracting numbe of common items.
#'
#' @param group_color_alpha The transparency of each ellipse color. Default 0.6.
#'
#' @param common_col_alpha The transparency of common circle color. Default 0.6.
#' @inheritParams base_plot_save
#' @inheritParams flower_plot_inner
#'
#' @param ... Parameters givento \code{\link{base_plot_save}}
#'
#' @return An image
#' @export
#'
#' @examples
#'
#' flowerinput <- "test.file"
#' flower(flowerinput)
#'
flower_plot <- function(input, sep="\t", row.names=NULL, header=T,
                   quote="", comment="", check.names=F,
                   item_variable = NULL,
                   set_variable = NULL,
                   start=90, a=0.5, b=2, r=1,
                   group_color="Spectral",
				   group_color_alpha=0.6,
				   label_total_num_items=TRUE,
				   saveplot=NULL,
                   label="core",common_color="white",common_color_alpha=1,saveppt=FALSE, ...) {

  data_m <- read.table(input, sep=sep, row.names=row.names, header=header,
                       quote=quote, comment=comment, check.names=check.names)

  datam_nodup <- unique(data_m)

  flower_input <- table(datam_nodup[,set_variable])
  sample <- names(flower_input)
  num_of_groups <- length(flower_input)
  total_num <- as.vector(flower_input)

  time_gene <- as.vector(table(datam_nodup[,item_variable]))
  core_num <-length(which(time_gene==num_of_groups))

  if(! label_total_num_items){
  	total_num <- total_num - core_num
  }

  group_color = generate_color_list(group_color, num_of_groups, group_color_alpha)
  commom_color = generate_color_list(common_color, 1, common_color_alpha)


  if(!is.null(saveplot)) {
    base_plot_save(saveplot, ...)
  }

  flower_plot_inner(sample=sample,
                    total_num=total_num, core_num=core_num,
                    start=start, a=a, b=b, r=r,
                    group_color=group_color,
                    label=label,common_color=common_color)


  if(!is.null(saveplot)) {
    dev.off()
  }
  if (saveppt){
	  library(eoffice)
	  library(ggplotify)
	  plot.new()
      p  <- as.ggplot(as.grob(function ()
	  flower_plot_inner(sample=sample,
                      total_num=total_num, core_num=core_num,
                      start=start, a=a, b=b, r=r,
                      group_color=group_color,
                      label=label,common_color=common_color)))
    eoffice::topptx(p, filename = paste0(saveplot,".pptx"))
    dev.off()
  }

}



#' Flower plot could be treated as one kind of venn diagram but only showing common items
#' like OTUs or genes among all groups and total items (or total items excluding common
#' items) for each group.
#'
#' Modified from http://blog.sciencenet.cn/blog-3406804-1159241.html
#'
#' This function is not planned for public usages.
#'
#' @param sample A vector of sample names.
#'
#' Like
#' ```
#' c("Grp1", "Grp2", "Grp3")
#' ```
#'
#' @param total_num Number of total or specififc items for each group.
#'
#' like
#'
#' ```
#' c(20, 30, 40)
#' ```
#'
#' @param core_num Number of items common to all groups.
#' @param r Set the size of the center circle.
#'
#' @param group_color Set the color of the petal ellipse (each group), with input format，like：c('#6181BD4E','#F348004E','#64A10E4E'...) or
#' supply a RColorBrewer color set like "Set1", "Set2", "Set3", "YlOrRd"
#' (check http://www.sthda.com/english/wiki/colors-in-r for more).
#'
#'
#' @param common_color The color of the center circle. Default "white".
#'
#' @param label The name of the center circle.
#'
#' @param start Start position of first ellipse. Default 90 represents starting from 0 clock.
#'
#' @inheritParams plotrix::draw.ellipse
#'
#' @return A pdf image.
#' @export
#'
#' @examples

#' flower_plot_inner(sample = sample_id, total_num = total_num, core_num = core_num)
#'
flower_plot_inner <- function(sample, total_num, core_num, start=90, a=0.5, b=2, r=1,
                        group_color=rgb(135, 206, 235, 150, max = 255),
                        group_color_alpha=0.6,
                        common_color_alpha=0.6,
                        label="core",common_color="white",...) {
  par( bty = 'n', ann = F, xaxt = 'n', yaxt = 'n', mar = c(1,1,1,1))
  plot(c(0,10),c(0,10),type='n')
  n   <- length(sample)
  deg <- 360 / n
  if (length(group_color)==1 && group_color == rgb(135, 206, 235, 150, max = 255)){
    group_color =c(rep(group_color, n))
  }
  res <- lapply(1:n, function(t){
    plotrix::draw.ellipse(x = 5 + cos((start + deg * (t - 1)) * pi / 180),
                          y = 5 + sin((start + deg * (t - 1)) * pi / 180),
                          col = group_color[t],
                          border = group_color[t],
                          a = a, b = b, angle = deg * (t - 1))
    text(x = 5 + 2.5 * cos((start + deg * (t - 1)) * pi / 180),
         y = 5 + 2.5 * sin((start + deg * (t - 1)) * pi / 180),
         total_num[t])

    if (deg * (t - 1) < 180 && deg * (t - 1) > 0 ) {
      text(x = 5 + 3.3 * cos((start + deg * (t - 1)) * pi / 180),
           y = 5 + 3.3 * sin((start + deg * (t - 1)) * pi / 180),
           sample[t],
           srt = deg * (t - 1) - start,
           adj = 1,
           cex = 1
      )
    } else {
      text(x = 5 + 3.3 * cos((start + deg * (t - 1)) * pi / 180),
           y = 5 + 3.3 * sin((start + deg * (t - 1)) * pi / 180),
           sample[t],
           srt = deg * (t - 1) + start,
           adj = 0,
           cex = 1
      )
    }
  })
  plotrix::draw.circle(x = 5, y = 5, r = r, col = common_color, border = NA)
  text(x = 5, y = 5, label=paste(label,core_num))
}



