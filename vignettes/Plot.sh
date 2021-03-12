line.r -f line.data -x Pos -y value -l Variable -c Variable -m TRUE
line.r -f line.data -x Pos -y value -m TRUE
line.r -f line.data -x Pos -y value -l Variable -m TRUE -s auto -V "cyan,purple"
# 下行代码在R脚本中运行，图中文字带反斜线
line.r -f line.data -x Pos -y value -l Variable -m TRUE -s auto -V "cyan,purple" -b "1000,2000" -P "100
0,4000" -B '\-1 kb,1 kb' -U "\-5 kb,5 kb"

splot_line.sh -f line.data -x Pos -y value -l Variable -c Variable -m TRUE
splot_line.sh -f line.data -x Pos -y value -m TRUE
splot_line.sh -f line.data -x Pos -y value -l Variable -m TRUE -s auto -V "cyan,purple"
splot_line.sh -f line.data -x Pos -y value -l Variable -m TRUE -s auto -V "cyan,purple" -P "-5000,0,5000" -
b "-1000,1000" -B "-1 kb,1 kb" -U "-5 kb,TSS,5 kb"


flowerplot.r -f flower.txt -c green
flowerplot.r -f flower.txt -A 1 -B 1.5 -r 1 -n FALSE -c Set2


splot_flowerplot.sh -f flower.txt -c green
splot_flowerplot.sh -f flower.txt -A 1 -B 1.5 -r 1 -n FALSE -c Set2



sp_boxplot.r -f box.txt  -m T -x Gene -y Expr -l Group
sp_boxplot.r -f box.txt  -m T -x Gene -y Expr -l Group  -z T -V T
sp_boxplot.r -f box.txt  -m T -x Gene -y Expr -l Group -M Set3 -v T -F T
sp_boxplot.r -f box.txt  -m T -x Gene -y Expr -l Group -M 'green,yellow,red'  -z T -k T  -v T
sp_boxplot.r -f box.txt  -m T -x Gene -y Expr -l Group -M 'green,yellow,red'  -z T -k T  -v T
sp_boxplot.r -f box.txt  -m T -x Gene -y Expr -l Group -M Set3 -v T -e Group -N 3 -W 1


splot_boxplot.sh -f box.txt  -m T -x Gene -y Expr -l Group
splot_boxplot.sh -f box.txt  -m T -x Gene -y Expr -l Group  -z T -V T
splot_boxplot.sh -f box.txt  -m T -x Gene -y Expr -l Group -M Set3 -v T -F T
splot_boxplot.sh -f box.txt  -m T -x Gene -y Expr -l Group -M 'green,yellow,red'  -z T -k T  -v T
splot_boxplot.sh -f box.txt  -m T -x Gene -y Expr -l Group -M 'green,yellow,red'  -z T -k T  -v T
splot_boxplot.sh -f box.txt  -m T -x Gene -y Expr -l Group -M Set3 -v T -e Group -N 3 -W 1




sp_pheatmap.r -f exprTable.txt -A 90 -c TRUE -d TRUE
sp_pheatmap.r -f exprTable.txt -z Set2
sp_pheatmap.r -f exprTable.txt -z "green,red" -s row
sp_pheatmap.r -f exprTable.txt -c TRUE -l log2 -z "YlOrRd"

sp_pheatmap.r -f exprTable.txt -A 90 -c TRUE -d TRUE -u exprTable.annorow.txt -v exprTable.annocol.txt
sp_pheatmap.r -f exprTable.txt -A 90 -c TRUE -d TRUE -u exprTable.annorow.txt -v exprTable.annocol.txt -S 'Type=c(TF="red",Enzyme="green"), Count=c("grey","blue")'


sp_pheatmap.r -f exprTable.txt -b "quantile" -i 20
sp_pheatmap.r -f exprTable.txt -b "0,5,10,20,40"  -z "YlGnBu"


splot_pheatmap.sh -f exprTable.txt -A 90 -c TRUE -d TRUE
splot_pheatmap.sh -f exprTable.txt -z Set2
splot_pheatmap.sh -f exprTable.txt -z "green,red" -s row
splot_pheatmap.sh -f exprTable.txt -c TRUE -l log2 -z "YlOrRd"

splot_pheatmap.sh -f exprTable.txt -A 90 -c TRUE -d TRUE -u exprTable.annorow.txt -v exprTable.annocol.txt
splot_pheatmap.sh -f exprTable.txt -A 90 -c TRUE -d TRUE -u exprTable.annorow.txt -v exprTable.annocol.txt -S 'Type=c(TF=\"red\",Enzyme=\"green\"), Count=c(\"grey\",\"blue\")'


splot_pheatmap.sh -f exprTable.txt -b "quantile" -i 20
splot_pheatmap.sh -f exprTable.txt -b "0,5,10,20,40"  -z "YlGnBu"

# Plot.Rmd中用对 sp_volcano_plot 函数测试了 4 次，这里分别用 R 和 Bash 脚本基于相同
# 的数据和参数测试 4 次。

sp_volcano.r -f volcano.txt -l log2FoldChange -d padj -s level

sp_volcano.r -f volcano.txt -l log2FoldChange -d padj -S "0.05,1" -p 'red,blue,black'

sp_volcano.r -f volcano.txt -l log2FoldChange -d padj -S "0.05,1" -p 'red,blue,black' -c TRUE

sp_volcano.r -f volcano.txt -l log2FoldChange -d padj -s level -P Symbol

# 目前的测试是
# transferRscriptToBashScript.py -i sp_volcano.r >splot_volcano.sh 可以直接转换 R 脚本为所需的
# bash 脚本
# 下面的测试可以运行成功

splot_volcano.sh -f volcano.txt -l log2FoldChange -d padj -s level

splot_volcano.sh -f volcano.txt -l log2FoldChange -d padj -S "0.05,1" -p 'red,blue,black'

splot_volcano.sh -f volcano.txt -l log2FoldChange -d padj -S "0.05,1" -p 'red,blue,black' -c TRUE

splot_volcano.sh -f volcano.txt -l log2FoldChange -d padj -s level -P Symbol


splot_venn2.sh -f vennDiagram.data -a Gene -b Sample

splot_venn2.sh -f vennDiagram.data -a Gene -b Sample -c "Set1, Set2,Set3"

sp_enrichment.R -f enrichment.data -x "SampleGroup" -y "Description" -c "Qvalue" -l "Qvalue"  -s "Count"
sp_enrichment.R -f enrichment.data -x "GeneRatio" -y "Description" -c "Qvalue" -l "Qvalue" -q "Count" -r "SampleGroup"
sp_enrichment.R -f goeast.enrich.txt -x log_odds_ratio -y Term -c p -l p -q q -s q -r Ontology -p right -a 0 -C FALSE -i FALSE -v "Pastel1"  -w 25 -W 25

splot_enrichment.sh -f enrichment.data -x "SampleGroup" -y "Description" -c "Qvalue" -l "Qvalue"  -s "Count" -r "SampleGroup" -p "right"  -C FALSE -i 60 -v "#89E767, #E84921" -t "Enrichment plot for multiple groups DE genes" -X "Percentage of DE genes in enriched GO terms" -w 20 -W 12.36
splot_enrichment.sh -f enrichment.data -x "GeneRatio" -y "Description" -c "Qvalue" -l "Qvalue" -q "Count" -r "SampleGroup"
splot_enrichment.sh -f goeast.enrich.txt -x log_odds_ratio -y Term -c p -l p -q q -s q -r Ontology -p right -a 0 -C FALSE -i FALSE -v "Pastel1" -w 25 -W 25

sp_upsetview.R -f upset.wide.data -v 0
sp_upsetview.R -f vennDiagram.data -v 2
sp_upsetview.R -f upset.wide.data -v 0 -n 5 -r "degree" -d FALSE -q "Samp1,Samp3"
sp_upsetview.R -f upset.wide.data -v 0 -n 5 -r "degree" -d FALSE -q "Samp1,Samp3" -N "Samp1,Samp3"

splot_upsetview.sh -f upset.wide.data -v 0
splot_upsetview.sh -f vennDiagram.data -v 2
splot_upsetview.sh -f upset.wide.data -v 0 -n 5 -r "degree" -d FALSE -q "Samp1,Samp3" -N "Samp1,Samp3"



scatter.R -f scatter_demo1.txt -x Gene -y Cluster -c Expr -d Percent -l Expr -M "white,blue"
scatter.R -f scatter.txt -x "X_val" -y Y_val -c Color -s Shape -d Size -l Samp -m "1,3,2" -n "2,1,3" -C "grp2,grp1,grp3" -S "cluster2,cluster1" -D 2 -j TRUE -e Color -A "free_y"
scatter.R -f scatter_demo1.txt -x Gene -y Cluster -c Expr -d Percent -l Expr -M "white,blue"


splot_scatter.sh -f scatter_demo1.txt -x Gene -y Cluster -c Expr -d Percent -l Expr -M "white,blue"
splot_scatter.sh -f scatter.txt -x "X_val" -y Y_val -c Color -s Shape -d Size -l Samp -m "1,3,2" -n "2,1,3" -C "grp2,grp1,grp3" -S "cluster2,cluster1" -D 2 -j TRUE -e Color -A "free_y"
splot_scatter.sh -f scatter_demo1.txt -x Gene -y Cluster -c Expr -d Percent -l Expr -M "white,blue"


sp_barplot.R -f bar.data -x ID -y Exper -c Gene -A TRUE -b "scale_y_log10()" -B stack -s identity -v TRUE -T TRUE
sp_barplot.R -f exprTable.txt -m F -B fill -T TRUE


splot_barplot.sh -f bar.data -x ID -y Exper -c Gene -A TRUE -b "scale_y_log10()" -B stack -s identity -v TRUE -T TRUE
splot_barplot.sh -f bar.data -x ID -c Gene -y Exper -b log2 -a 0 -B stack -s identity -p "right" -u TRUE -v TRUE -A TRUE -T TRUE -P FALSE -r 0 -D 2  -w 25 -W 15
splot_barplot.sh -m FALSE -f exprTable.txt -b NULL -a 0 -p "top" -u TRUE -T TRUE -P FALSE -r 0 -w 25 -W 15
splot_barplot.sh -m FALSE -f exprTable.txt -b NULL -a 0 -p "top" -u TRUE -T TRUE -P FALSE -r 0 -w 25 -W 15 -i TRUE


sp_histogram.R -f histogram_demo2.txt -S density  -P line


splot_histogram.sh -f histogram_demo2.txt -S density  -P line -i identity
