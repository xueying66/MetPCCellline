
# Figure.2a ---------------------------------------------------------------

load( file = '/home/liuxueying/Proj_MetPCCelline/output/data/mutation.burden.df.RData')

ggplot(mutation.burden.df, aes(x = rank, y = mutation.burden, color = type2))+
  geom_point(size = 3)+ 
  geom_text_repel(data = subset(mutation.burden.df, type =='CCLE' & mutation.burden>9),
                  aes(label = id),
                  nudge_y = 0.1)+
  xlab('Rank') + ylab('Mutation burden')+
  geom_hline(yintercept = 9, color = "black")+
  scale_color_manual(values = type_colors) +  # 这里指定颜色
  theme_bw()+
  theme(
    
    plot.title = element_text(hjust = 0.5, size = 25),
    axis.title = element_text(size = 25, hjust = 0.5, color = 'black'),     # 坐标轴标题加粗
    axis.text = element_text(size = 25, color = 'black'),      # 坐标轴刻度加粗
    panel.border = element_rect(color = "black", linewidth = 1.5, fill = NA), # 加粗边框
    axis.ticks.length = unit(0.4, "cm"),         # 增加刻度线长度
    
  )
ggsave('/home/liuxueying/Proj_MetPCCelline/output/plot/three.dataset.mutation.burden.pdf', width = 12,height = 10)




# Figure.2b ---------------------------------------------------------------

load(file = '/home/liuxueying/Proj_MetPCCelline/output/MSI.Rdata')



p1 <- oncoPrint(
  mutation.profile,
  get_type = function(x) {ifelse(x == 1, "MUT", "background")},
  row_names_side = "left",          # 行名放置在左侧
  show_heatmap_legend = F,      # 隐藏图例
  alter_fun = alter_fun,            # 突变样式函数
  col = col,                        # 颜色定义
  show_column_names = FALSE,        # 不显示列名
  show_row_names = T, # 不显示行名
  top_annotation = col.annotation,  # 顶部注释
  right_annotation = NULL,          # 确保右侧无条形图
  show_pct = FALSE,                 # 禁用百分比显示
  row_names_gp = gpar(fontsize = 40, fontface = 'bold'),
  column_order = sample_order,
  row_order = gene.order
)


))


p2 <- Heatmap(mutation.frequency,
              show_column_dend = FALSE, 
              col=colorRamp2(c(0,1),c('blue','red')),
              row_names_gp = gpar(fontsize = 40,fontface='bold'),
              width = unit(4, "cm"),
              show_heatmap_legend = T,##显示图例
              top_annotation=col.ha,
              show_column_names = F,
              cluster_columns = F,
              cluster_rows = F,
              column_order = c("MSI.mutation.fequency", "MSI.no.mutation.fequency"), # 强制指定列顺序
              row_order = gene.order,
              name = 'Mutation \nfrequency',
              heatmap_legend_param = list(
                title_gp = gpar(fontsize = 20, fontface = "bold"), # 图例标题字体
                labels_gp = gpar(fontsize = 20),                  # 图例标签字体
                legend_height = unit(10, "cm"),                    # 图例高度
                legend_width = unit(5, "cm")                      # 图例宽度
              ))





pdf('/home/liuxueying/Proj_MetPCCelline/output/plot/mutation.profile.MMR.pdf', width = 30, height = 15)
p1+p2
dev.off()



# Fig.2c ------------------------------------------------------------------

load(file = '/home/liuxueying/Proj_MetPCCelline/output/data/de_res.RData')

ggplot(de.res, aes(x = log2FoldChange, y = -log10(padj),color=sig))+
  geom_point()+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black")+
  geom_vline(xintercept = c(-1,1), linetype = "dashed", color = "black")+
  geom_point()+
  scale_color_manual(values=c("#2f5688","#BBBBBB","#CC0000"))+
  scale_size_continuous(range = c(0,1))+
  theme(legend.title = element_blank() )+
  geom_label_repel(data = subset(de.res, gene %in% de_genes),
                   max.overlaps = getOption("ggrepel.max.overlaps", default = 20),
                   aes(label = gene),
                   size = 8, 
                   color = 'black',fontface = "bold.italic",
                   segment.color = "black",   
                   segment.linetype = "solid", 
                   segment.size = 0.5 ,
                   box.padding = 3,        
                   point.padding = 0.2 ) +
  xlab("log2FC")+
  ylab("-log10 p-adj")+
  theme_bw()+
  theme(
    legend.position = 'none',
    plot.title = element_text(hjust = 0.5, size = 30),
    axis.title = element_text(size = 30, hjust = 0.5, color = 'black'),     # 坐标轴标题加粗
    axis.text = element_text(size = 30, color = 'black'),      # 坐标轴刻度加粗
    panel.border = element_rect(color = "black", linewidth = 1.5, fill = NA), # 加粗边框
    axis.ticks.length = unit(0.4, "cm"),         # 增加刻度线长度
    
  )
ggsave('/home/liuxueying/Proj_MetPCCelline/output/plot/de.hyper.non.volcano.pdf', width = 10, height = 10)


# Fig.2d ------------------------------------------------------------------

load(file = '/home/liuxueying/Proj_MetPCCelline/output/data/de.enrich.result.RData')

ggplot(de.enrich.result, aes(x = gene_ratio, y = reorder(Term, count,decreasing = F))) +
  geom_point(aes(size = count, fill = P.value), shape =21) +
  scale_fill_gradient(low = "red", high = "blue") +
  
  labs(
    
    size = "Count",
    x = "Gene Ratio",
    y = "",
    title = "MSigDB Hallmark"
  ) +
  theme_bw()+
  theme(
    plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
    axis.title.y = element_text(size = 25, color = 'black'),
    axis.title.x = element_text(size = 25,  color = 'black'),
    axis.text.x = element_text(size = 22,  color = 'black'),
    axis.text.y = element_text(size = 25, color = 'black'),
    panel.border = element_rect(color = "black", linewidth = 1.2, fill = NA),
    
    axis.ticks.length = unit(0.2, "cm"),  # Increase tick mark length
    axis.ticks = element_line(size = 1)  # Increase tick mark thickness
  )

ggsave('/home/liuxueying/Proj_MetPCCelline/output/plot/DE.MSI.non.MSI.ORA.pdf', width = 16, height = 12)

