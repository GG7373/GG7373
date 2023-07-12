
###############################################################
getwd()
coont =read.csv("D:/R/DESer2-method/20230421Lafin/count.csv",row.names = 1)
dim(coont)

washcoont = coont[rowSums(coont) != 0,]
dim(washcoont)

write.csv(washcoont,file="D:/R/WASHCOUNT.csv",
          quote = F)
fenzu = read.csv("D:/R/DESer2-method/20230421Lafin/mymeta.csv",stringsAsFactors = T)

colnames(washcoont) == fenzu$id
View(fenzu)
library(DESeq2)

dds = DESeqDataSetFromMatrix(countData=washcoont, #清洗后数据
                             colData=fenzu, #分组文件
                             design= ~ dex) #分组数据集列名  形成矩阵
dds = DESeq(dds)

#####################################################################

one = results(dds, c("dex", "DW","MOD"))

#挑选数据进行对比

res1 = results(dds)#存储结果

head(res1)

class(one)#判定数据类型

washone = data.frame(one)#将结果转变为数据框

library(dplyr)
#调用新包 显示基因的上调或下调
washone %>% 
  mutate(group = case_when( 
    log2FoldChange >= 0.263 & pvalue <= 0.05 ~ "UP",
    log2FoldChange <= -0.322 & pvalue <= 0.05 ~ "DOWN",
    TRUE ~ "NOT_CHANGE"
  )) -> washedone
table(washedone$group)
write.csv(washedone,file="D:/R/666/MOD-DWEL.csv",
          quote = F)
View(washedone)
oneup = washedone[washedone$group=="UP",]
onedown = washedone[washedone$group=="DOWN",]

View(oneup)
class(oneup)

write.csv(oneup,file="D:/R/666/MOD-DWEL-UP.csv",
          quote = F)
write.csv(onedown,file="D:/R/666/MOD-DWEL-DOWN.csv",
          quote = F)
############

two = results(dds,c("dex", "MOD","NOR"))#挑选数据进行对比

res2 = results(dds)#存储结果

class(two)#判定数据类型
washtwo = data.frame(two)#将结果转变为数据框
library(dplyr)
#调用新包 显示基因的上调或下调
washtwo %>% 
  mutate(group = case_when( 
    log2FoldChange >= 0.263 & pvalue <= 0.05 ~ "UP",
    log2FoldChange <= -0.322 & pvalue <= 0.05 ~ "DOWN",
    TRUE ~ "NOT_CHANGE"
  )) -> washedtwo

table(washedtwo$group)
View(washedtwo)
write.csv(washedtwo,file="D:/R/666/MOD-NOR.csv",
          quote = F)
twoup = washedtwo[washedtwo$group=="UP",]
twodown = washedtwo[washedtwo$group=="DOWN",]
View(twoup)
write.csv(washedtwo,file="D:/R/666/MOD-NOR.csv",
          quote = F)
write.csv(twoup,file="D:/R/666/MOD-NOR-UP.csv",
          quote = F)
write.csv(twodown,file="D:/R/666/MOD-NOR-DOWN.csv",
          quote = F)


###########################

three = results(dds, c("dex", "WE", "MOD"))
View(three)
res3 = results(dds)#存储结果
head(res3)

class(three)#判定数据类型
washthree = data.frame(three)#将结果转变为数据框
library(dplyr)
#调用新包 显示基因的上调或下调
washthree %>% 
  mutate(group = case_when( 
    log2FoldChange >= 0.263 & pvalue <= 0.05 ~ "UP",
    log2FoldChange <= -0.352 & pvalue <= 0.05 ~ "DOWN",
    TRUE ~ "NOT_CHANGE"
  )) -> washedthree
table(washedthree$group)
write.csv(washedthree,file="D:/R/666/WE-MOD.csv",
          quote = F)

threeup = washedthree[washedthree$group=="UP",]
threedown = washedthree[washedthree$group=="DOWN",]

View(oneup)
class(oneup)

write.csv(threeup,file="D:/R/666/MOD-WEL-UP.csv",
          quote = F)
write.csv(threedown,file="D:/R/666/MOD-WEL-DOWN.csv",
          quote = F)
#############
four = results(dds, c("dex", "MOD", "WE"))

res4 = results(dds)#存储结果
head(res4)

class(four)#判定数据类型
washfour = data.frame(four)#将结果转变为数据框
library(dplyr)
#调用新包 显示基因的上调或下调
washfour %>% 
  mutate(group = case_when( 
    log2FoldChange >= 0.263 & pvalue <= 0.05 ~ "UP",
    log2FoldChange <= -0.322 & pvalue <= 0.05 ~ "DOWN",
    TRUE ~ "NOT_CHANGE"
  )) -> washedfour
table(washedfour$group)
fourup = washedfour[washedfour$group=="UP",]
fourdown = washedfour[washedfour$group=="DOWN",]
View(washedfour)
write.csv(washedfour,file="D:/R/MOD-WE.csv",
          quote = F)
write.csv(fourup,file="D:/R/MOD-WE-up.csv",
          quote = F)
write.csv(fourdown,file="D:/R/MOD-WE-DOWN.csv",
          quote = F)
