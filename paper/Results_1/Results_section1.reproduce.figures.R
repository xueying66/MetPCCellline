
# Figure.1a ---------------------------------------------------------------


load(file = '/home/liuxueying/Proj_MetPCCelline/output/data/mutation.profile.RData')
load(file = '/home/liuxueying/Proj_MetPCCelline/output/data/ccle.mutation.profile.RData')
load(file = '/home/liuxueying/Proj_MetPCCelline/output/data/mutation.frequency.RData')
p1 <- oncoPrint(
  mutation.profile,
  get_type = function(x) {ifelse(x == 1, "MUT", "background")},
  row_names_side = "left",          
  show_heatmap_legend = FALSE,     
  alter_fun = alter_fun,            
  col = col,                        
  show_column_names = FALSE,       
  show_row_names = F, # 不显示行名
  top_annotation = col.annotation, 
  right_annotation = NULL,         
  show_pct = FALSE,                 
  row_names_gp = gpar(fontsize = 40, fontface = 'bold'),
  column_order = sample_order,
  row_order = gene.order)


p2 <- oncoPrint(
  ccle.mutation.profile,
  get_type = function(x) {ifelse(x == 1, "MUT", "background")},
  #row_names_side = "right",          
  show_heatmap_legend = FALSE,      
  alter_fun = alter_fun,            
  col = col,                        
  # show_column_names = T,        
  top_annotation = col.annotation,  
  right_annotation = NULL,          
  show_pct = FALSE,                 
  row_names_gp = gpar(fontsize = 40, fontface = 'bold'),
  show_row_names = F,
  
  row_order = gene.order,
  column_order = sample_order,
)



p3 <- grid.grabExpr({
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(nrow = 1, ncol = 2, widths = unit(c(9, 2), "null"))))
  
  pushViewport(viewport(layout.pos.col = 1))
  draw(p1, newpage = FALSE)
  upViewport()
  
  pushViewport(viewport(layout.pos.col = 2))
  draw(p2, newpage = FALSE)
  upViewport(2)
})

pdf('aa.pdf', width = 90, height = 40)
grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow = 1, ncol = 2, widths = unit(c(9, 2), "null"))))

pushViewport(viewport(layout.pos.col = 1))
draw(p1, newpage = FALSE)
upViewport()

pushViewport(viewport(layout.pos.col = 2))
draw(p2, newpage = FALSE)
upViewport(2)
dev.off()

p4 <- Heatmap(mutation.frequency,
              show_column_dend = FALSE, 
              col=colorRamp2(c(0,1),c('blue','red')),
              row_names_gp = gpar(fontsize = 40,fontface='bold'),
              width = unit(1, "cm"),
              show_heatmap_legend = T,
              top_annotation=col.ha,
              show_column_names = F,
              cluster_columns = F,
              cluster_rows = F,
              column_order = c("CCLE.mutation.fequency"), 
              row_order = gene.order,
              name = 'Mutation \nfrequency',
              heatmap_legend_param = list(
                title_gp = gpar(fontsize = 20, fontface = "bold"),
                labels_gp = gpar(fontsize = 20),                  
                legend_height = unit(10, "cm"),                   
                legend_width = unit(5, "cm")                      
              ))


grob_p4 <- grid.grabExpr({
  draw(p4, newpage = FALSE)
})



pdf("//Fig.1a.pdf", width = 90, height = 40)


grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow = 1, ncol = 2, widths = unit(c(4, 0.2), "null"))))


pushViewport(viewport(layout.pos.col = 1))
grid.draw(p3)
upViewport()


pushViewport(viewport(layout.pos.col = 2))
grid.draw(grob_p4)
upViewport(2)
dev.off()


# Figure.1b ---------------------------------------------------------------

load(file = '/home/liuxueying/Proj_MetPCCelline/output/data/number.genes.RData')
ggplot(number.genes, aes(x = name, y = mutation_burden)) +
  geom_col(width = 0.6, fill = "#D8CDEB") +
  labs(y = "Number of mutated genes") +
  xlab('')+
  theme_bw()+
  theme(
    legend.position = 'none',
    plot.title = element_text(hjust = 0.5, size = 15),
    axis.title = element_text(size = 15, hjust = 0.5, color = 'black'),     # 坐标轴标题加粗
    axis.text = element_text(size = 15, color = 'black'),      # 坐标轴刻度加粗
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.border = element_rect(color = "black", linewidth = 1.5, fill = NA), # 加粗边框
    axis.ticks.length = unit(0.2, "cm"),         # 增加刻度线长度
    
  )
ggsave('/home/liuxueying/Proj_MetPCCelline/plot/Number.of.mutated.genes.pdf', width = 5, height = 5)


# Figure.1c  --------------------------------------------------------------

load(file = '/home/liuxueying/Proj_MetPCCelline/output/data/mutation.burden.RData')
ggplot(mutation.burden,aes(x= name, y = mutation.burden))+
  geom_col(width = 0.6, fill = "#D8CDEB")+
  xlab('Rank')+
  ylab('Mutation burden')+
  theme_bw()+
  theme(
    legend.position = 'none',
    plot.title = element_text(hjust = 0.5, size = 15),
    axis.title = element_text(size = 15, hjust = 0.5, color = 'black'),     # 坐标轴标题加粗
    axis.text = element_text(size = 15, color = 'black'),      # 坐标轴刻度加粗
    
    panel.border = element_rect(color = "black", linewidth = 1.5, fill = NA), # 加粗边框
    axis.ticks.length = unit(0.2, "cm"),         # 增加刻度线长度
    
  )

ggsave('/home/liuxueying/Proj_MetPCCelline/output/plot/Fig.1c.pdf', width = 5.5, height = 5)

# Figure.1d ---------------------------------------------------------------
load(file = '/home/liuxueying/Proj_MetPCCelline/output/data/ratio.RData')

ggplot(ratio,aes(x= name, y = aa)) + 
  geom_col(width = 0.6, fill = "#D8CDEB")+
  ylab('Ratio')+
  xlab('')+
  theme_bw()+
  theme(
    legend.position = 'none',
    plot.title = element_text(hjust = 0.5, size = 15),
    axis.title = element_text(size = 15, hjust = 0.5, color = 'black'),     # 坐标轴标题加粗
    axis.text = element_text(size = 15, color = 'black'),      # 坐标轴刻度加粗
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.border = element_rect(color = "black", linewidth = 1.5, fill = NA), # 加粗边框
    axis.ticks.length = unit(0.2, "cm"),         # 增加刻度线长度
    
  )

ggsave('/home/liuxueying/Proj_MetPCCelline/output/plot/Figure.1d.pdf', width = 5, height = 5)


# Figure.1e ---------------------------------------------------------------


load(file = '/home/liuxueying/Proj_MetPCCelline/output/data/supplementary.table1.RData')
hotspot.maf = read.maf(maf = supplementary.table1)


pdf('./pcdata/paper/plot/section1/AR.pdf', width = 15, height = 10)
lollipopPlot(
  maf = hotspot.maf,
  gene = 'AR',
  AACol = 'HGVSp_Short',
  showMutationRate = TRUE,
  labelPos =c(702,742, 875,878)
  
)
dev.off()


pdf('./pcdata/paper/plot/section1/TP53.pdf', width = 15, height = 10)
lollipopPlot(
  maf = hotspot.maf,
  gene = 'TP53',
  AACol = 'HGVSp_Short',
  showMutationRate = TRUE,
  labelPos = c(282,273,248,245,220,179,175)
  
)
dev.off()



# Figure.1f ---------------------------------------------------------------

load(file = '/home/liuxueying/Proj_MetPCCelline/output/data/hot.mutation.code.ccle.RData')


p1 <- oncoPrint(
  hot.mutation.code.ccle,
  get_type = function(x) {ifelse(x == 1, "MUT", "background")},
  row_names_side = "left",          # 行名放置在左侧
  show_heatmap_legend = F,      # 隐藏图例
  alter_fun = alter_fun,            # 突变样式函数
  col = col,                        # 颜色定义
  show_column_names = T,        # 显示列名
  top_annotation = NULL,        # 确保顶部无条形图
  right_annotation = NULL,          # 确保右侧无条形图
  show_pct = FALSE,                 # 禁用百分比显示
  row_names_gp = gpar(fontsize = 20),    #行名字体
  column_names_gp = gpar(fontsize = 25), # 列名字体
  column_order = sample.order, #行名顺序
  column_names_rot = 45)        # 列名旋转45
pdf('/home/liuxueying/Proj_MetPCCelline/output/plot/Figure.1f.pdf', width = 12, height = 9)
p1
dev.off()


# Figure.1g ---------------------------------------------------------------

load(file = '/home/liuxueying/Proj_MetPCCelline/output/data/dm.cnv.RData')
ggplot(dm.cnv, aes(x = delta, y = -log10(p.adj),color=sig))+
  geom_point()+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black")+
  geom_vline(xintercept = c(-0.3,0.3), linetype = "dashed", color = "black")+
  
  scale_color_manual(values=c("#2f5688","#BBBBBB","#CC0000"))+
  scale_size_continuous(range = c(0,1))+
  theme(legend.title = element_blank() )+
  geom_label_repel(data = subset(rs.df, gene %in% 'AR'),
                   max.overlaps = getOption("ggrepel.max.overlaps", default = 20),
                   aes(label = gene),
                   size = 5, 
                   color = 'black',fontface = "bold.italic",
                   segment.color = "black",   
                   segment.linetype = "solid", 
                   segment.size = 0.5 ,
                   box.padding = 3,        
                   point.padding = 0.2 ) +
  xlab("Difference")+
  ylab("-log10 p-adj")+
  xlim(-0.9, 0.9) +
  ggtitle('')+
  theme_bw()+
  theme(
    legend.position = 'none',
    plot.title = element_text(hjust = 0.5, size = 15),
    axis.title = element_text(size = 15, hjust = 0.5, color = 'black'),     # 坐标轴标题加粗
    axis.text = element_text(size = 15, color = 'black'),      # 坐标轴刻度加粗
    axis.text.x = element_text( hjust = 1, vjust = 1),
    panel.border = element_rect(color = "black", linewidth = 1.5, fill = NA), # 加粗边框
    axis.ticks.length = unit(0.2, "cm"),         # 增加刻度线长度
    
  )


ggsave('/home/liuxueying/Proj_MetPCCelline/output/plot/all.dm.CNV.pdf', width = 5, height = 5)




