

# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'
#   Generate DOC:              'Ctrl + Shift + Alt + r'

# 1. 一般不在函数里面加载包，可放在外面运行测试时，如在 Plot.Rmd 中
# 2. 函数的描述需要为英文
# 3. 函数名前面增加 sp_ 前缀，以免与其它包的函数名重复，引起冲突
# 4. 变量名尽量通用，比如 otu_table 就比 phylum 好
# 5. 输入数据需要是常见格式，比如 OTU 表或物种表都是一行一个 OTU 或物种，需要转置时在程序内转置
# 6. 要对数据做格式效验和一致性判断
# 7. 尽量利用已有的函数和写法
# 8. 用到函数时要查看函数的帮助，看哪些参数是有意义的，提取出来作为参数，如下面的otu_scale_methods

#'  Generate RDA plot with species table and env table
#'
#' @param otu_table Normalized OTU/Species abundance data frame or data file (with header line, the first column will be treated as the row names, tab separated)
#' @param env_table Environment factor data frame or data file (with header line, the first column will be treated as the row names, tab separated)
#' @param sample_min_reads_count The minimum allowed sample reads count. Default 10000 meaning samples with total reads count less than 10000 would be filtered.
#' @param rare_otu_count Definite rare OTU. OTU with abundance in all samples less than given value would be filtered.
#' @param otu_scale_method Popular (and effective) standardization methods for community ecologists
#' like total, max, frequency, normalize, range, rank, rrank, standardize, hellinger, chi.square, rclr,
#' log, clr, alr.
#' @inheritParams sp_ggplot_layout
#' @inheritParams sp_pcoa
#' @param ... Parameters given to `sp_ggplot_layout`
#'
#' @return A ggplot2 object
#' @export
#'
#' @examples
#'
#' ## Not run:
#' input <- "KO-OE_all.txt"
#' volcano_plot(input,log2fc_var='Log2FoldChange',fdr_var='Padj',status_col_var='',
#'              title="sd",label="Label",log10_transform_fdr=TRUE,point_size=5)
#' ## End(Not run)


sp_rda <-
  function(
    otu_table,
    env_table,
    metadata,
    sample_min_reads_count = 10000,
    rare_otu_count = 5,
    otu_scale_method = "hellinger"
    ) {

  options(warn = -1)
  options(scipen = 999)

  if ("character" %in% class(otu_table)) {
    otu_table <- sp_readTable(otu_table, row.names = 1, stringsAsFactors = T, check.names = T)
  } else if (!"data.frame" %in% class(data)) {
    stop("Unknown input format for `data` parameter.")
  }

  # 样本测序量评估，根据测序量调整最小样本量，默认为10000
  # 按样本量筛选，低于最小样本量丢弃，100个样品只有76个大于10000条序列
  otu_table = otu_table[,colSums(otu_table) >= sample_min_reads_count]
  otu_table = as.data.frame(t(otu_table))

  # 过滤掉rare OTU，即所有样品累计小于5条reads的OTU
  otu_table = otu_table[, colSums(otu_table) > rare_otu_count]


  if ("character" %in% class(env_table)) {
    env_table <- sp_readTable(env_table, row.names = 1, stringsAsFactors = T, check.names = T)
  } else if (!"data.frame" %in% class(data)) {
    stop("Unknown input format for `data` parameter.")
  }

  # 去除缺失NA的数据，防止以下分析报错
  env_table = na.omit(env_table)
  # 将环境因子进行log+1自然对数转化，使得环境因子数据更均一
  env = log1p(env)

  # 选取 otu 表和 env表共有的样品，保证样品的顺序一致和名字匹配
  # matchedL = match_two_df(otu_table, env_table, way="row-row")
  # otu_table = matchedL$df1
  # env_table = matchedL$df2

  if ("character" %in% class(metadata)) {
    metadata <- sp_readTable(metadata, row.names = 1, stringsAsFactors = T, check.names = T)
  } else if (!"data.frame" %in% class(data)) {
    stop("Unknown input format for `data` parameter.")
  }

  # 选取 otu 表、env表、metdata共有的样品，保证样品的顺序一致和名字匹配
  otu_rowname = rownames(otu_table)
  env_rowname = rownames(env_table)
  metadata_rownames = rownames(metadata)

  common_rownameL <- intersect(otu_rowname, env_rowname)
  common_rownameL <- intersect(common_rownameL, metadata_rownames)

  otu_table <- otu_table[common_rownameL,,drop=F]
  env_table <- env_table[common_rownameL,,drop=F]
  metadata <- metadata[common_rownameL,,drop=F]

  ## Hellinger预转化（处理包含很多0值的群落物种数据时，推荐使用）
  if (!sp.is.null(otu_scale_method)){
    otu_table = decostand(otu_table, otu_scale_method)
  }


  # Adonis分析环境因子显著性
  # 见 env.Rmd 中的部分


  #数据选择是否中心化
  #直接使用原始数据，不做转化。对于群落物种组成数据来讲（因为通常包含很多 0 值），不是很推荐
  rda_result <- rda(otu_table ~ ., env_table, scale = FALSE)

  #选择数据的转化方法，如hellinger等
  #物种数据 Hellinger 预转化（处理包含很多 0 值的群落物种数据时，推荐使用）
  otu_table_hel <- decostand(otu_table, method = 'hellinger')

  #使用全部的环境数据
  rda_tb <- rda(otu_table_hel ~ ., env, scale = FALSE)

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
  rda_tb_forward_p <-
    ordistep(
      rda(otu_table_hel ~ 1, env, scale = FALSE),
      scope = formula(rda_tb),
      direction = 'forward',
      permutations = 999
    )

  #vegan 包 ordiR2step() 前向选择，基于 999 次置换检验
  rda_tb_forward_r <-
    ordiR2step(
      rda(otu_table_hel ~ 1, env, scale = FALSE),
      scope = formula(rda_tb),
      R2scope = rda_adj,
      direction = 'forward',
      permutations = 999
    )

  #细节部分查看
  summary(rda_tb_forward_r, scaling = 1)

  #比较选择前后校正后 R2 的差异
  RsquareAdj(rda_tb)$adj.r.squared
  RsquareAdj(rda_tb_forward_r)$adj.r.squared





  ##ggplot2 作图，以前向选择后的简约模型 rda_tb_forward_r 为例，展示前两轴，II 型标尺，双序图，默认使用物种加权和计算的样方坐标
  #提取样方和环境因子排序坐标，前两轴，II 型标尺
  rda_tb_forward_r.scaling2 <-
    summary(rda_tb_forward_r, scaling = 2)
  rda_tb_forward_r.site <-
    data.frame(rda_tb_forward_r.scaling2$sites)[1:2]
  rda_tb_forward_r.env <-
    data.frame(rda_tb_forward_r.scaling2$biplot)[1:2]


  #读取样本分组数据（附件“group.txt”）
  group <-
    read.delim(
      'group.txt',
      sep = '\t',
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  #合并样本分组信息，构建 ggplot2 作图数据集
  rda_tb_forward_r.site$sample <- rownames(rda_tb_forward_r.site)
  rda_tb_forward_r.site <-
    merge(rda_tb_forward_r.site, group, by = 'sample')
  rda_tb_forward_r.env$sample <- NA
  rda_tb_forward_r.env$group <- rownames(rda_tb_forward_r.env)

  library(ggplot2)
  #ggplot2 作图
  p <- ggplot(rda_tb_forward_r.site, aes(RDA1, RDA2)) +
    geom_point(aes(color = group)) +
    scale_color_manual(values = c('red', 'orange', 'green3')) +
    theme(
      panel.grid = element_blank(),
      panel.background = element_rect(color = 'black', fill = 'transparent'),
      legend.title = element_blank(),
      legend.key = element_rect(fill = 'transparent')
    ) +
    labs(x = 'RDA1 (42.91%)', y = 'RDA2 (9.80%)') +
    geom_vline(xintercept = 0,
               color = 'gray',
               size = 0.5) +
    geom_hline(yintercept = 0,
               color = 'gray',
               size = 0.5) +
    geom_segment(
      data = rda_tb_forward_r.env,
      aes(
        x = 0,
        y = 0,
        xend = RDA1,
        yend = RDA2
      ),
      arrow = arrow(length = unit(0.1, 'cm')),
      size = 0.3,
      color = 'blue'
    ) +
    geom_text(
      data = rda_tb_forward_r.env,
      aes(RDA1 * 1.1, RDA2 * 1.1, label = group),
      color = 'blue',
      size = 3
    )
  p

  # 保存图形到PDF文件
  ggsave("RDA_plot.pdf",
         plot = p,
         width = 8,
         height = 6)

}
