#绘制带物种数据的RDA环境因子图
library(vegan)

##读取数据
#读入物种数据，细菌门水平丰度表（OTU 水平数据量太大，后续的置换检验和变量选择过程很费时间，不方便做示例演示）
phylum <- read.delim('phylum_table.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
#读取环境数据
env <- read.delim('env_table.txt', row.names = 1, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)

#数据选择是否中心化
#直接使用原始数据，不做转化。对于群落物种组成数据来讲（因为通常包含很多 0 值），不是很推荐
rda_result <- rda(phylum~., env, scale = FALSE)

#选择数据的转化方法，如hellinger等
#物种数据 Hellinger 预转化（处理包含很多 0 值的群落物种数据时，推荐使用）
phylum_hel <- decostand(phylum, method = 'hellinger')

#使用全部的环境数据
rda_tb <- rda(phylum_hel~., env, scale = FALSE)

##R2 校正
#RsquareAdj() 提取 R2
r2 <- RsquareAdj(rda_tb)
rda_noadj <- r2$r.squared	#原始 R2
rda_adj <- r2$adj.r.squared	#校正后的 R2



##变量选择
#计算方差膨胀因子
vif.cca(rda_tb)

##可以看出有几个变量的VIF值特别大（超过10甚至20），所以有必要剔除一些变量
#多元回归变量筛选通常有三种模式：前向选择（forward selection）、后向选择（backward selection）以及逐步选择（stepwise selection，也称双向选择，forward-backward selection）。其中前向选择在RDA分析中最为常用，以下将以前向选择为例简介RDA模型中的变量选择方法。
#vegan 包 ordistep() 前向选择，基于 999 次置换检验


#vegan 包 ordistep() 前向选择，基于 999 次置换检验
rda_tb_forward_p <- ordistep(rda(phylum_hel~1, env, scale = FALSE), scope = formula(rda_tb), direction = 'forward', permutations = 999)

#vegan 包 ordiR2step() 前向选择，基于 999 次置换检验
rda_tb_forward_r <- ordiR2step(rda(phylum_hel~1, env, scale = FALSE), scope = formula(rda_tb), R2scope = rda_adj, direction = 'forward', permutations = 999)

#细节部分查看
summary(rda_tb_forward_r, scaling = 1)

#比较选择前后校正后 R2 的差异
RsquareAdj(rda_tb)$adj.r.squared
RsquareAdj(rda_tb_forward_r)$adj.r.squared





##ggplot2 作图，以前向选择后的简约模型 rda_tb_forward_r 为例，展示前两轴，II 型标尺，双序图，默认使用物种加权和计算的样方坐标
#提取样方和环境因子排序坐标，前两轴，II 型标尺
rda_tb_forward_r.scaling2 <- summary(rda_tb_forward_r, scaling = 2)
rda_tb_forward_r.site <- data.frame(rda_tb_forward_r.scaling2$sites)[1:2]
rda_tb_forward_r.env <- data.frame(rda_tb_forward_r.scaling2$biplot)[1:2]


#读取样本分组数据（附件“group.txt”）
group <- read.delim('group.txt', sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
#合并样本分组信息，构建 ggplot2 作图数据集
rda_tb_forward_r.site$sample <- rownames(rda_tb_forward_r.site)
rda_tb_forward_r.site <- merge(rda_tb_forward_r.site, group, by = 'sample')
rda_tb_forward_r.env$sample <- NA
rda_tb_forward_r.env$group <- rownames(rda_tb_forward_r.env)

library(ggplot2)
#ggplot2 作图
p <- ggplot(rda_tb_forward_r.site, aes(RDA1, RDA2)) +
  geom_point(aes(color = group)) +
  scale_color_manual(values = c('red', 'orange', 'green3')) +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), legend.title = element_blank(), legend.key = element_rect(fill = 'transparent')) + 
  labs(x = 'RDA1 (42.91%)', y = 'RDA2 (9.80%)') +
  geom_vline(xintercept = 0, color = 'gray', size = 0.5) + 
  geom_hline(yintercept = 0, color = 'gray', size = 0.5) +
  geom_segment(data = rda_tb_forward_r.env, aes(x = 0,y = 0, xend = RDA1,yend = RDA2), arrow = arrow(length = unit(0.1, 'cm')), size = 0.3, color = 'blue') +
  geom_text(data = rda_tb_forward_r.env, aes(RDA1 * 1.1, RDA2 * 1.1, label = group), color = 'blue', size = 3)
p

# 保存图形到PDF文件
ggsave("RDA_plot.pdf", plot = p, width = 8, height = 6)
