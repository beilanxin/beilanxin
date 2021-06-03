library(tidyverse)
library(dplyr)
options(scipen = 200)
#把测序cluster信息和misa位点信息整合起来，构建输出总表output0
multi <- read.table("multi_table.txt")
only <- read.table("1v1_table.txt")
colnames(only) <- c("SSR_ID", "Chromosome", "SSR_Order", "SSR_Type", "Unit_by_MISA", "Length", "MISA_Start", "MISA_End", "Cluster_Rep_Seq", "BLAST_Start", "Blast_End")
colnames(multi) <- c("SSR_ID", "Chromosome", "SSR_Order", "SSR_Type", "Unit_by_MISA", "Length", "MISA_Start", "MISA_End", "Cluster_Rep_Seq", "BLAST_Start", "Blast_End")
#为什么要merge以后再筛选呢，哈哈
only <- filter(only, !Cluster_Rep_Seq %in% only$Cluster_Rep_Seq[duplicated(only$Cluster_Rep_Seq)])#都是因为复杂微卫星的存在导致一个cluster的唯一的blast结果对应到两个SSR位点，删掉
multi <- filter(multi, !Cluster_Rep_Seq %in% multi$Cluster_Rep_Seq[duplicated(multi$Cluster_Rep_Seq)])#一对多原因：（1）一条cluster序列blast到基因组上两个位点；（2）复杂微卫星的存在导致一个cluster的唯一的blast结果对应到两个SSR位点。删掉
cluster <- read.table("cluster_rep_result.txt")
cluster <- cluster[ , c(1, 2, 4, 6)]
colnames(cluster) <- c("SSLP", "Cluster_ID", "Cluster_Rep_Seq", "Unit_by_Seq")
any(duplicated(cluster))
all(multi$Cluster_Rep_Seq %in% cluster$Cluster_Rep_Seq)
all(only$Cluster_Rep_Seq %in% cluster$Cluster_Rep_Seq)
multi_m <- merge(multi, cluster, by = "Cluster_Rep_Seq", all.x = T, all.y = F)
only_m <- merge(only, cluster, by = "Cluster_Rep_Seq", all.x = T, all.y = F)
dim(multi)
dim(multi_m)
dim(only)
dim(only_m)
any(duplicated(multi_m))
any(duplicated(only_m))
any(multi_m$id %in% only_m$id)
CLUSTER0 <- rbind(multi_m, only_m) %>% arrange(Cluster_ID)
#开始对CLUSTER0的重复数据进行筛选
#（1）对一个cluster中有两个以上SSR的类型下手（全部删掉）
CLUSTER1 <- CLUSTER0[!CLUSTER0$Cluster_Rep_Seq %in% CLUSTER0[grep("_2$", CLUSTER0$Cluster_ID), ]$Cluster_Rep_Seq, ]
#（2）对一个ssr对应多个cluster的情况下手
#暴力删掉聚成很多类的cluster，保留其中SSLP最大的
#选出SSLP最大的（有简单的方法为什么要写循环嘞）
CLUSTER1 <- CLUSTER1 %>% arrange(desc(SSLP)) %>% arrange(SSR_ID)
CLUSTER2 <- CLUSTER1[!duplicated(CLUSTER1$SSR_ID), ]
dim(CLUSTER2)
any(duplicated(CLUSTER2$Cluster_ID))
any(duplicated(CLUSTER2$Cluster_Rep_Seq))
any(duplicated(CLUSTER2$SSR_ID))
CLUSTER <- CLUSTER2
#CLUSTER的重复数据筛选完成
#把misa总的结果整合上
#misa就是misa软件的输出文件,川哥的是用提取出来的输入进misa那肯定很多复杂微卫星发现不了嘛
misa <- read.table("perfectssr_150cut1.txt")
colnames(misa) <- c("SSR_ID", "Chromosome", "SSR_Order", "SSR_Type", "Unit_by_MISA", "Length", "MISA_Start", "MISA_End")
CLUSTERpMISA <- merge(misa, CLUSTER[ , -c(3:9)], by = "SSR_ID", all = T)
dim(misa)
dim(CLUSTERpMISA)

#把引物信息整合
CLUSTERpMISA$Primer_ID <- paste(CLUSTERpMISA$Chromosome, ".1:", as.character(CLUSTERpMISA$MISA_Start-150), "-", as.character(CLUSTERpMISA$MISA_End+150), sep = "")
primers_only <- read.table("primers_only_1v1.txt")
primers_only$V1 <- substr(primers_only$V2, 1, nchar(primers_only$V2)-2)
primers_only <- primers_only[ , -c(2, 4)]
colnames(primers_only) <- c("Primer_ID", "F_Primer_Start", "R_Primer_End")
primers_seq <- read.table("WB_ssr_noc.primers.txt", fill = TRUE, na.strings = "NA")
primers_seq <- primers_seq[-1, ]
colnames(primers_seq) <- c("Primer_ID", "F_Primer_ID", "F_Primer_Seq", "R_Primer_ID", "R_Primer_Seq")
primer <- merge(primers_only, primers_seq, by = "Primer_ID", all.x = T, all.y = F)
CLUSTERpMISApPRIMER <- merge(CLUSTERpMISA, primer, by = "Primer_ID", all.x = T, all.y = F) %>% arrange(SSR_ID)
which((((as.numeric(CLUSTERpMISApPRIMER$MISA_Start) -20) < as.numeric(CLUSTERpMISApPRIMER$F_Primer_Start)) | (as.numeric(CLUSTERpMISApPRIMER$MISA_End + 20) > as.numeric(CLUSTERpMISApPRIMER$R_Primer_End))) == F)
which((((as.numeric(CLUSTERpMISApPRIMER$MISA_Start) -20) < as.numeric(CLUSTERpMISApPRIMER$F_Primer_Start)) | (as.numeric(CLUSTERpMISApPRIMER$MISA_End + 20) > as.numeric(CLUSTERpMISApPRIMER$R_Primer_End))) == T)
which((((as.numeric(CLUSTERpMISApPRIMER$MISA_Start) -20) < as.numeric(CLUSTERpMISApPRIMER$F_Primer_Start)) & (as.numeric(CLUSTERpMISApPRIMER$MISA_End + 20) > as.numeric(CLUSTERpMISApPRIMER$R_Primer_End))) == T)
CLUSTERpMISApPRIMER$Chromosome <- paste(CLUSTERpMISApPRIMER$Chromosome, ".1", sep = "")
filter(CLUSTERpMISApPRIMER, R_Primer_End > 0) %>% View()
filter(CLUSTERpMISApPRIMER, R_Primer_End > 0)[ , c(1, 3, 8, 9)] %>% write.table("genome_for_sep.txt", quote = F, na = "NA", col.names = TRUE, row.names = F)
filter(CLUSTERpMISApPRIMER, SSLP >= 2, R_Primer_End > 0) %>% dim()
filter(CLUSTERpMISApPRIMER, SSLP >= 2, R_Primer_End > 0)[ , c(1, 3, 8, 9)] %>% write.table("seq_for_sep.txt", quote = F, na = "NA", col.names = TRUE, row.names = F)
list1_sep <- read.table("genome_for_sep_sorted_interval_100kb.list", header = T)
colnames(list1_sep)[2] <- "Genome_Sep"
list2_sep <- read.table("seq_for_sep_sorted_interval_10kb.list", header = T)
colnames(list2_sep)[2] <- "Seq_Sep"
CLUSTERpMISApPRIMER <- merge(CLUSTERpMISApPRIMER, list1_sep, by = "Primer_ID", all.x = T, all.y = F) %>% arrange(SSR_ID)
CLUSTERpMISApPRIMER <- merge(CLUSTERpMISApPRIMER, list2_sep, by = "Primer_ID", all.x = T, all.y = F) %>% arrange(SSR_ID)
dim(CLUSTERpMISApPRIMER)
CLUSTERpMISApPRIMER$Chromosome_ID <- factor(as.numeric(substr(CLUSTERpMISApPRIMER$Chromosome, 7, 9))-695, levels = 1:37)
write.csv(filter(CLUSTERpMISApPRIMER,  R_Primer_End > 1), "1_genome.csv", quote = F, na = "NA", col.names = TRUE, row.names = F)
write.csv(filter(CLUSTERpMISApPRIMER, SSLP > 0, R_Primer_End > 1), "2_sequencing.csv", quote = F, na = "NA", col.names = TRUE, row.names = F)

table(CLUSTERpMISApPRIMER$In_1Mb_Sep_1)
table(CLUSTERpMISApPRIMER$In_1Mb_Sep_2)
CLUSTERpMISApPRIMER$Chromosome_ID <- gsub(37, "X", CLUSTERpMISApPRIMER$Chromosome_ID)
CLUSTERpMISApPRIMER$Chromosome_ID <- factor(CLUSTERpMISApPRIMER$Chromosome_ID, levels = c(1:36, "X"))
CLUSTERpMISApPRIMER$Unit <- str_extract(CLUSTERpMISApPRIMER$Unit_by_MISA,"(?<=\\().+?(?=\\))")
annotation <- read.table("perfectssr_150cut_annotated_simp.txt")
annotation$Primer_ID <- paste(annotation$V1, ":", as.character(annotation$V2-150), "-", as.character(annotation$V3+150), sep = "")
annotation$V2 == CLUSTERpMISApPRIMER$MISA_Start
annotation <- annotation[ , -c(1:3)]
colnames(annotation)[c(1, 2)] <- c("Annotationf", "Annotation")
CLUSTERpMISApPRIMER <- merge(CLUSTERpMISApPRIMER, annotation, by = "Primer_ID") %>% arrange(SSR_ID)
CLUSTERpMISApPRIMER[stringr::str_detect(CLUSTERpMISApPRIMER$Annotationf, "exon"), ]$SSLP %>% mean(na.rm = T)
CLUSTERpMISApPRIMER$SSLP %>% mean(na.rm = T)
filter(CLUSTERpMISApPRIMER, R_Primer_End > 0)$SSLP %>% mean(na.rm = T)
t.test(CLUSTERpMISApPRIMER[stringr::str_detect(CLUSTERpMISApPRIMER$Annotationf, "exon"), ]$SSLP, CLUSTERpMISApPRIMER$SSLP)


for(i in 1:nrow(CLUSTERpMISApPRIMER)){
  if(CLUSTERpMISApPRIMER$Unit[i] %in% c("AC", "CA", "GT", "TG")){
    CLUSTERpMISApPRIMER$Unitss[i] <- "AC"
  }else if(CLUSTERpMISApPRIMER$Unit[i] %in% c("AG", "GA", "CT", "TC")){
    CLUSTERpMISApPRIMER$Unitss[i] <- "AG"
    }else if(CLUSTERpMISApPRIMER$Unit[i] %in% c("AT", "TA")){
      CLUSTERpMISApPRIMER$Unitss[i] <- "AT"
    }else if(CLUSTERpMISApPRIMER$Unit[i] %in% c("CG", "Gc")){
      CLUSTERpMISApPRIMER$Unitss[i] <- "CG"
    }else{CLUSTERpMISApPRIMER$Unitss[i] <- "Others"}
}

CLUSTERpMISApPRIMER$With_Primer[is.na(CLUSTERpMISApPRIMER$F_Primer_Start)] <- "N"
CLUSTERpMISApPRIMER$With_Primer[CLUSTERpMISApPRIMER$F_Primer_Start > 0] <- "Y"
table(CLUSTERpMISApPRIMER$With_Prime)
CLUSTERpMISApPRIMER$Annotation[CLUSTERpMISApPRIMER$Annotation == "lncRNA"] <- "lnc_RNA"
which(CLUSTERpMISApPRIMER$Annotation == "t_RNA" )
CLUSTERpMISApPRIMER$SSLP[is.na(CLUSTERpMISApPRIMER$SSLP)] <- "Unknown"
write.csv(CLUSTERpMISApPRIMER, "all.csv", quote = F, na = "NA", col.names = TRUE, row.names = F)

