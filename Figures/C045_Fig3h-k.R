rds<-readRDS('expr_matrix_q3_log2plus1.rds') #DSP data with q3 Normalization stored as Seurat rds
rds$Group<-paste0(rds$Location,'.',rds$Treatment)
rds$Group<-factor(rds$Group,levels=c('stroma.Pre-treatment','stroma.Post-treatment','tumor.Pre-treatment','tumor.Post-treatment'))
table(rds$Group)

#用治疗TherapeuticResponse分组
dataset<-read.csv('DSP_patients_Clinical.xls',header=T,sep='\t')
rds$TherapeuticResponse<-plyr::mapvalues(x =rds$Patient_ID,from = dataset$Patient_ID,to = dataset$TherapeuticResponse)
exprs <- data.frame(FetchData(object = rds, vars = 'SOX4'))
rds$SOX4<-exprs$SOX4
data<-rds@meta.data
write.table(data,file='DSP_SOX4-expression_patients_Clinical.xls',quote=F,sep='\t',row.names=F)



###----------------------------------------------------------------------------
#Fig3h
library(ggplot2)
dataset<-read.table('DSP_SOX4-expression_patients_Clinical.xls',header=T,sep='\t')
data1<-subset(dataset,Location=='tumor')[,c('roi','SOX4')]
colnames(data1)<-c('roi','SOX4.tumor')
data2<-subset(dataset,Location=='stroma')[,c('roi','SOX4')]
colnames(data2)<-c('roi','SOX4.stroma')
data<-merge(data1,data2,by='roi')
# 假设你的数据在 data 中
# 1. 计算 Pearson 相关系数和 p 值
cor_test <- cor.test(data$SOX4.tumor, data$SOX4.stroma, method = "pearson")
pcc_value <- cor_test$estimate   # 相关系数
p_value   <- cor_test$p.value    # p值

pdf("Fig3h.pdf", width = 6, height = 5)
ggplot(data, aes(x = SOX4.tumor, y = SOX4.stroma)) +
  geom_point(shape = 21,                     # 圆圈形状
             fill = "#ABD6EE",              # 填充颜色
             color = "black",               # 边框颜色
             size = 3,                      # 点大小
             stroke = 1,                    # 边框粗细
             alpha = 0.8) +                 # 透明度
  geom_smooth(method = "lm", 
              se = TRUE,                    # 显示置信区间阴影
              fill = "#ABD6EE",             # 阴影颜色
              alpha = 0.3,                  # 阴影透明度
              size = 0.8, 
              color = "#ABD6EE", 
              linetype = "dashed") +        # 拟合线虚线
  theme_classic(base_size = 14) +
  theme(
    axis.line = element_line(size = 0.6, colour = "black"),
    axis.ticks = element_line(size = 0.6, colour = "black"),
    axis.text = element_text(size = 12, colour = "black"),
    axis.title = element_text(size = 14, colour = "black", face = "bold")
  ) +
  annotate("text",
           x = 7.5,
           y = 9,
           label = paste0("PCC = ", round(pcc_value, 3),
                          "\nPvalue = ", signif(p_value, 3)),
           hjust = 1, vjust = 1, size = 4, color = "#ABD6EE", fontface = "bold")

dev.off()





###----------------------------------------------------------------------------
#Fig3i
data<-read.csv('Signature_gene_sets_gseva_group.csv',header=T,check.names=F)
comparisons<-list(c('SOX4-high','SOX4-low'))
#'stroma.Pre−treatment.high','stroma.Pre−treatment.low','stroma.Post−treatment.high','stroma.Post−treatment.low'
genes<-c('Angiogenesis',"TCell Terminal Differentiation","Bystander","TAM activation","Immunosuppression","DC maturation")
da<-subset(data,Group %in% c('stroma.Pre-treatment','stroma.Post-treatment'))
for (gene in genes){
	print(gene)
	plot_df<-da
	# 改列名方便 ggplot 调用
	plot_df$expression<-plot_df[,gene]	
	vmin<-min(plot_df[,gene])
	vmax<-max(plot_df[,gene])
	yvalue<-(vmin+vmax)/2.0
	ysize<-3
	# 绘图
	ptest1<- stat_compare_means(method = "wilcox.test",
						 comparisons = comparisons,
						 #label = "p.signif",
						 size = ysize) #,vjust = 1) )
	ptest2<- stat_compare_means(method='wilcox.test',
						 comparisons = comparisons,paired=FALSE,label = "p.format",
							   label.x.npc = "center",label.y.npc = "top",vjust = 1,size=ysize,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),symbols = c("***", "**", "*", "ns")),hide.ns = F)  # 根据 P 值映射颜色,不work aes(color = ifelse(..p.. < 0.05, "significant", "not significant"))
	p <- ggplot(plot_df, aes(x = SOX4.group, y = expression, fill = SOX4.group)) +
	  geom_violin(scale = "width", trim = FALSE) +
	  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.5) +
	 stat_summary (fun.data = function (x) data.frame (y=yvalue,label = paste (round(mean (x),2))), geom="text",size=ysize, vjust = -0.8, position = position_dodge(width = 0.9)) + ## 使用 position_dodge 分散标签
	  ptest1 +
	  facet_wrap(~Group, scales = "free_y", nrow = 1) +  # 一行展示所有细胞类型scale_fill_manual(values = group_colors) +
	  theme_bw(base_size = 10) +
	  theme(axis.text.x = element_text(color="black",size=10,angle=30, hjust = 1),
		   axis.text.y = element_text(color="black",size=8,angle=0),
			axis.title.x=element_text(size = 10),axis.title.y=element_text(size = 8),
			strip.text = element_text(size = 10),
			panel.grid = element_blank(),
			legend.position = "none") +
	  ylim(vmin, vmax+0.8) +
	  ylab(paste0('GSVA score')) +
	  xlab("") +
	  ggtitle(gene)
	# 保存 PDF
	pdf(paste0(gene, "_stroma-Pre-Post-SOX4-high-vs-low_VlnPlot_wilcox.test.pdf"), width = 5, height = 2.5)
	print(p)
	dev.off()
}




###----------------------------------------------------------------------------
#Fig3j poor vs good
rds$Group<-paste0(rds$Location,'.',rds$Treatment)
rds$Group<-factor(rds$Group,levels=c('stroma.Pre-treatment','stroma.Post-treatment','tumor.Pre-treatment','tumor.Post-treatment'))
table(rds$Group)
#用治疗TherapeuticResponse分组
dataset<-read.csv('DSP_patients_Clinical.xls',header=T,sep='\t')
rds$TherapeuticResponse<-plyr::mapvalues(x =rds$Patient_ID,from = dataset$Patient_ID,to = dataset$TherapeuticResponse)

tmp<-rds
gene<-'SOX4'
exprs <- data.frame(FetchData(object = rds, vars = gene))
tmp$SOX4<-exprs$SOX4
gene<-'SOX4'
comparisons<-list(c('poor','good'))
vmin<-min(tmp@meta.data[,gene])
vmax<-max(tmp@meta.data[,gene])
yvalue<-9


# 提取表达矩阵（注意 slot="data"）
expr <- FetchData(tmp, vars = gene, slot = "data")
expr$cell <- rownames(expr)

# 合并信息
plot_df <- cbind(expr, tmp@meta.data[, c("TherapeuticResponse", "Group")])
plot_df$Group<-factor(plot_df$Group,levels=c('stroma.Pre-treatment','stroma.Post-treatment','tumor.Pre-treatment','tumor.Post-treatment'))


# 改列名方便 ggplot 调用
colnames(plot_df)[1] <- "expression"
# 绘图
ptest1<- stat_compare_means(method = "wilcox.test",
					 comparisons = comparisons,
					 #label = "p.signif",
					 size = ysize,
					 vjust = 0.5) 
p <- ggplot(plot_df, aes(x = TherapeuticResponse, y = expression, fill = TherapeuticResponse)) +
  geom_violin(scale = "width", trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.5) +
 stat_summary (fun.data = function (x) data.frame (y=yvalue,label = paste (round(mean (x),2))), geom="text",size=ysize, vjust = -0.8, position = position_dodge(width = 0.9)) + ## 使用 position_dodge 分散标签
  ptest1 +
  facet_wrap(~Group, scales = "free_y", nrow = 1) +  # 一行展示所有细胞类型scale_fill_manual(values = group_colors) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(color="black",size=10,angle=30, hjust = 1),
	   axis.text.y = element_text(color="black",size=8,angle=0),
		axis.title.x=element_text(size = 10),axis.title.y=element_text(size = 8),
		strip.text = element_text(size = 10),
		panel.grid = element_blank(),
		legend.position = "none") +
  ylim(vmin, vmax*1.1) +
  ylab(paste0(gene, " expression")) +
  xlab("") +
  ggtitle(gene)


# 保存 PDF
#pdf(paste0(gene, "_Location-Treatment-TherapeuticResponse_VlnPlot_t.test-cmp1.pdf"), width = 10, height = 2.5)
pdf(paste0('Fig3j.pdf'), width = 10, height = 2.5)
print(p)
dev.off()




###----------------------------------------------------------------------------
#Fig3k Pre-treatment vs Post-treatment
rds$Group<-paste0(rds$Location,'.',rds$TherapeuticResponse)
rds$Group<-factor(rds$Group,levels=c('stroma.poor','stroma.good','tumor.poor','tumor.good'))
table(rds$Group)

tmp<-rds
exprs <- data.frame(FetchData(object = rds, vars = ge))
tmp$SOX4<-exprs$SOX4
gene<-'SOX4'
comparisons<-list(c('Pre-treatment','Post-treatment'))
vmin<-min(tmp@meta.data[,gene])
vmax<-max(tmp@meta.data[,gene])
yvalue<-9

# 提取表达矩阵（注意 slot="data"）
expr <- FetchData(tmp, vars = gene, slot = "data")
expr$cell <- rownames(expr)

# 合并信息
plot_df <- cbind(expr, tmp@meta.data[, c("Group", "Treatment")])
plot_df$Group<-factor(plot_df$Group,levels=c('stroma.poor','stroma.good','tumor.poor','tumor.good'))
plot_df$Treatment<-factor(plot_df$Treatment,levels=c('Pre-treatment','Post-treatment'))

# 改列名方便 ggplot 调用
colnames(plot_df)[1] <- "expression"

# 绘图
ptest1<- stat_compare_means(method = "wilcox.test",
					 comparisons = comparisons,
					 #label = "p.signif",
					 size = ysize,
					 vjust = 0.5) 
p <- ggplot(plot_df, aes(x = Treatment, y = expression, fill = Treatment)) +
  geom_violin(scale = "width", trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.5) +
 stat_summary (fun.data = function (x) data.frame (y=yvalue,label = paste (round(mean (x),2))), geom="text",size=ysize, vjust = -0.8, position = position_dodge(width = 0.9)) + ## 使用 position_dodge 分散标签
  ptest1 +
  facet_wrap(~Group, scales = "free_y", nrow = 1) +  # 一行展示所有细胞类型scale_fill_manual(values = group_colors) +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(color="black",size=10,angle=30, hjust = 1),
	   axis.text.y = element_text(color="black",size=8,angle=0),
		axis.title.x=element_text(size = 10),axis.title.y=element_text(size = 8),
		strip.text = element_text(size = 10),
		panel.grid = element_blank(),
		legend.position = "none") +
  ylim(vmin, vmax*1.1) +
  ylab(paste0(gene, " expression")) +
  xlab("") +
  ggtitle(gene)


# 保存 PDF
#pdf(paste0(gene, "_Location-Treatment-TherapeuticResponse_VlnPlot_t.test-cmp2.pdf"), width = 10, height = 2.5)
pdf(paste0('Fig3k.pdf'), width = 10, height = 2.5)
print(p)
dev.off()

