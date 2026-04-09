
if (FALSE){
	install.packages("ggplot2", repo="http://cran.us.r-project.org")
	install.packages("ggfortify", repo="http://cran.us.r-project.org")
	install.packages("data.table", repo="http://cran.us.r-project.org")
	install.packages("ggrepel", repo="http://cran.us.r-project.org")
}
library(plyr)
library(ggplot2)
library(grid)
library(data.table, quietly=T)
library(ggfortify)

if (FALSE) {
	library("ggrepel")
}



#data <- read.table(file="system_file/pca/pca.txt", sep="\t", quote="", comment="", header=T, row.names=1, check.names=F)

if (1 == 0){
        data <- read.table(file="system_file/pca/pca.txt", sep="\t", header=T, row.names=1,
                check.names=F, quote="", comment="")
} else if (1 == 1) {
        data <- read.table(file="system_file/pca/pca.txt", sep="\t", header=T, row.names=NULL,
                check.names=F, quote="", comment="")
        rownames_data <- make.names(data[,1],unique=T)
        data <- data[,-1,drop=F]
        rownames(data) <- rownames_data
}

data <- data[rowSums(abs(data))!=0, ]

data$mad <- apply(data, 1, mad) 

data <- data[data$mad>0.5, ]

data <- data[order(data$mad, decreasing=T), 1:(dim(data)[2]-1)]

dim_data <- dim(data)

data_row_num <- dim_data[1]

if (5000 != 0 & 5000 < data_row_num) {
	data <- data[1:5000, ]
}

data <- as.data.frame(t(data))

sampleL = rownames(data)

if ("system_file/pca/pcaphenodata_.txt" == "") {
	data_t_label <- data
	data_t_label$group = sampleL
	data_t_label$Row.names = sampleL
} else {
	grp_data <- read.table("system_file/pca/pcaphenodata_.txt", sep="\t", quote="", header=T, check.names=F, row.names=NULL, comment="")
	rownames(grp_data) <- grp_data[,1]
	data_t_label <- merge(data, grp_data, by=0, all.x=T)
	rownames(data_t_label) <- data_t_label$Row.names
	data_t_label <- data_t_label[match(sampleL, data_t_label$Row.names), ]
	#data_t_label <- grp_data[na.omit(sampleL, rownames(grp_data)), ]
}

shape_order <- c()

if ("Batch" != "c_t_c_t0304") {
if (length(shape_order) > 1) {
	data_t_label$Batch <- factor(data_t_label$Batch, levels=shape_order, ordered=T)
} else {
	data_t_label$Batch <- factor(data_t_label$Batch)
}
}

color_order <- c()

if (length(color_order) > 1) {
	data_t_label$Conditions <- factor(data_t_label$Conditions, levels=color_order, ordered=T)
}


if ("nolog" != "nolog"){
	data <- nolog(data + 0)
}



color_v <- c()

if ("Batch" != "c_t_c_t0304") {
	shape_level <- length(unique(data_t_label$Batch))
	shapes = (1:shape_level)%%30
}

pca <- prcomp(data, scale=TRUE)

rotation = pca$rotation

  outputprefix = "system_file/pca/pca.txt"
  if (outputprefix == "") {
  outputprefix = "system_file/pca/pca.txt"
}

write.table(rotation, file=paste0(outputprefix,".pca.weights.xls"), sep="\t", quote=F, row.names=T, col.names=T)

system("sed -i 's/^/ID\t/' system_file/pca/pca.txt.pca.weights.xls")

x = pca$x

write.table(x, file=paste0(outputprefix,".pca.pcs.xls"), sep="\t", quote=F, row.names=T,  col.names=T)

system("sed -i 's/^/ID\t/' system_file/pca/pca.txt.pca.pcs.xls")

# sdev: standard deviation of the principle components.
# Square to get variance
percentVar <- pca$sdev^2 / sum( pca$sdev^2)

percentVar2 <- as.data.frame(percentVar)
rownames(percentVar2) <- colnames(x)

write.table(percentVar2, file=paste0(outputprefix,".pca.pc_variance.xls"), sep="\t", quote=F, row.names=T)


if (2 == 2) {
	p <- autoplot(pca, data=data_t_label, alpha=0) + ggtitle("")+coord_fixed()
	
	if ("Diameters" != ""){
		p <- p + aes(size=Diameters)
	}
	if ("Conditions" != "c_t_c_t0304") {
		p <- p + aes(colour=Conditions)
		if (length(color_v) >= 2) {
			if (is.numeric(data_t_label$Conditions)){
				p <- p + scale_colour_gradient(low=color_v[1], high=color_v[2], name="Conditions")
			}else {
				p <- p + scale_color_manual(values=color_v)
			}
		}
	}

	if ("Batch" != "c_t_c_t0304") {
		p <- p + aes(shape=Batch)
		if ( shape_level > 6) {
			p <- p + scale_shape_manual(values=shapes)
		}
	}
	#if (("Diameters" != "") && ("Conditions" != "c_t_c_t0304") && ("Batch" != "c_t_c_t0304")) {
	#	p <- autoplot(pca, data=data_t_label, colour="Conditions", shape="Batch", size="Diameters", alpha=0)  
	#} else if (("Conditions" != "c_t_c_t0304") && ("Batch" != "c_t_c_t0304")) {
	#	p <- autoplot(pca, data=data_t_label, colour="Conditions", shape="Batch",alpha=0)  
	#} else if (("Diameters" != "") && ("Batch" != "c_t_c_t0304")) {
	#	p <- autoplot(pca, data=data_t_label, shape="Batch", size="Diameters", alpha=0)  
	#} else if (("Diameters" != "") && ("Conditions" != "c_t_c_t0304")) {
	#	p <- autoplot(pca, data=data_t_label, colour="Conditions", size="Diameters", alpha=0)  
	#} else if ("Diameters" != "") {
	#	p <- autoplot(pca, data=data_t_label, size="Diameters", alpha=0)  
	#} else if ("Conditions" != "c_t_c_t0304") {
	#	p <- autoplot(pca, data=data_t_label, colour="Conditions", alpha=0)  
	#} else if ("Batch" != "c_t_c_t0304") {
	#	p <- autoplot(pca, data=data_t_label, shape="Batch", alpha=0)  
	#} else {
	#	p <- autoplot(pca, data=data_t_label)  
	#}

	#if (("Conditions" != "c_t_c_t0304") && length(color_v) == 2) {
	#	p <- p + scale_colour_gradient(low=color_v[1], high=color_v[2], name="Conditions")
	#}

	#if (("Batch" != "c_t_c_t0304") && shape_level > 6) {
	#	p <- p + scale_shape_manual(values=shapes)
	#}


	if (FALSE) {
		#p <- p + geom_text(aes(label=Row.names), position="identity",
		#hjust=0, size=0, check_overlap=FALSE)
		if ("0" != "") {
			p <- p + geom_text_repel(aes(label=Row.names), show.legend=F, size=0)
		} else {
			p <- p + geom_text_repel(aes(label=Row.names), show.legend=F)
		}
	}

	p <- p + xlab(paste0("PC1 (", round(percentVar[1]*100), "% variance)")) + 
		ylab(paste0("PC2 (", round(percentVar[2]*100), "% variance)"))

	p <- p 
  


  ggsave(filename = paste0(outputprefix, '.pca.pdf'),p)
  
} else {
	library(scatterplot)	
	if ("Conditions" != "c_t_c_t0304") { 
		group = data_t_label$Conditions
		colorA <- rainbow(length(unique(group)))
		
		colors <- colorA[as.factor(group)]
		
		colorl <- colorA[as.factor(unique(group))]
	}

	if ("Batch" != "c_t_c_t0304") { 
		# 获得PCH symbol列表
		group <- data_t_label$Batch
		pch_l <- as.numeric(as.factor(unique(group)))
		# 产生每个样品的pch symbol
		pch <- pch_l[as.factor(group)]
	}

	pc <- as.data.frame(pca$x)
	outputprefix = "system_file/pca/pca.txt"
  if (outputprefix == "") {
  outputprefix = "system_file/pca/pca.txt"
}
	filename = paste0(outputprefix, '.pca.pdf')

	pdf(filename)
	scatterplot3d(x=pc$PC1, y=pc$PC2, z=pc$PC3, pch=pch, color=colors, xlab=paste0("PC1 (", round(percentVar[1]*100), "% variance)"), ylab=paste0("PC2 (", round(percentVar[2]*100), "% variance)"), zlab=paste0("PC3 (", round(percentVar[3]*100), "% variance)"))

	legend(-3,8, legend=levels(as.factor(Conditions)), col=colorl, pch=pch_l, xpd=T, horiz=F, ncol=6)
	dev.off()
}

