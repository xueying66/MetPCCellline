

# Figure.3a ---------------------------------------------------------------


load(file = '/home/liuxueying/Proj_MetPCCelline/output/data/met500.df.RData')
ggplot(met500.df, aes(x = rank, y = correlation, color = prostate)) +
  geom_point(size = 3) +  # 绘制点，设置大小
  labs(title = "", x = "Rank", y = "Transcriptome similarity") +
  scale_color_manual(values = c("prostate" = "red", "others" = "black")) +  
  geom_point(data = subset(df, name == "PC3_PROSTATE"), color = "red", size = 4) +
  geom_text_repel(data = subset(df, name == "PC3_PROSTATE"| name == 'LNCAPCLONEFGC_PROSTATE' | name == 'VCAP_PROSTATE' | name == 'MDAPCA2B_PROSTATE'), 
                  aes(label = gsub(x= name, pattern = '_PROSTATE','')), 
                  nudge_y = 0.02,  # 标签在纵坐标上的偏移量
                  color = "black",
                  size = 4,
                  arrow = arrow(length = unit(0.02, "npc")))+
  theme(legend.position = 'none')

ggsave('/home/liuxueying/Proj_MetPCCelline/output/plot/MET500.adenocarcinoma.pdf', plot = p1,width = 15, height = 10)


# Figure.3b ---------------------------------------------------------------

load(file = '/home/liuxueying/Proj_MetPCCelline/output/data/mat.RData')

pdf(file  = '/home/liuxueying/Proj_MetPCCelline/output/plot/bulk.site.specific.pdf',width=20,height=20)
chart.Correlation(mat, histogram=FALSE, pch=19,method = 'spearman',cex.cor.tuning = 0.7,cor.range = c(-1,1))
dev.off()


# Figure.3c ---------------------------------------------------------------

load(file = '/home/liuxueying/Proj_MetPCCelline/output/data/scRNA_df.RData')

ggplot(scRNA_df, aes(x = rank, y = correlation, color = prostate)) +
  geom_point(size = 2) +  
  labs(title = "Adenocarcinoma scRNA-seq", x = "Rank", y = "Transcriptome similarity") +
  scale_color_manual(values = c("prostate" = "red", "others" = "black")) +  
  geom_point(data = subset(df, prostate == "prostate"), color = "red", size = 3)+
  ggplot.style +
  theme(legend.position = 'none')+
  
  
  geom_segment(data = df_highlight, 
               aes(x = rank - 50, xend = rank, y = correlation, yend = correlation),  
               arrow = arrow(length = unit(0.01, "npc")), 
               color = "black") +
  
  
  geom_text(data = df_highlight, 
            aes(x = rank - 50, y = correlation, label = gsub("_PROSTATE", "", name)),  
            hjust = 1,  # 右对齐标签
            color = "black", 
            size = 4)


ggsave('/home/liuxueying/Proj_MetPCCelline/output/plot/scRNA-seq.Adenocarcinoma.pdf',width = 20,height=15 )


# Figure.3d ---------------------------------------------------------------
load( file = '/home/liuxueying/Proj_MetPCCelline/output/bone.ln.df.RData')
ggplot(bone.ln.df, aes(x= ln.cor.median, y=bone.cor.median))+
  geom_point(size=8)+
  theme_bw()+
  xlab('Lymph node')+
  ylab('Bone')+
  ggtitle('')+
  xlim(c(0,0.3))+ylim(c(0,0.3))+geom_abline(slope = 1,intercept = 0)+
  annotate(
    'text',
    x = 0.2,            # x 轴位置（在两个柱中间）
    y = 0.15,           # y 轴位置（高于图顶端）
    label = paste0("Spearman-rank correlation = ", round(bone.ln.cor, digits = 2)), # p 值标签
    size = 5,           # 字体大小
    color = "black",      # 字体颜色
    hjust = 0           # 水平对齐方式
  )+theme_bw()+
  theme(
    legend.position = 'none',
    plot.title = element_text(hjust = 0.5, size = 25),
    axis.title = element_text(size = 25, hjust = 0.5, color = 'black'),     # 坐标轴标题加粗
    axis.text = element_text(size = 25, color = 'black'),      # 坐标轴刻度加粗
    panel.border = element_rect(color = "black", linewidth = 1.5, fill = NA), # 加粗边框
    axis.ticks.length = unit(0.4, "cm"),         # 增加刻度线长度
    
  )


ggsave('/home/liuxueying/Proj_MetPCCelline/output/plot/bone.ln.cor.pdf',width = 20,height=15 )






