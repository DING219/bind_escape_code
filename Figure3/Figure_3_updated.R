

########################################Figure 2 clear version#######################################
rm(list=ls())
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(cowplot)
library(slider)
library(dplyr)
library(lubridate)
library(ggtrendline)
library(data.table)
library(patchwork)
##lodad data

Alpha <- read.csv("D:/all_project/Italy_project/manuscript/figures/updated_figure2/updated/Alpha_Italy_bind_esc_socre1.csv")
Beta <- read.csv("D:/all_project/Italy_project/manuscript/figures/updated_figure2/updated/Beta_Italy_bind_esc_socre1.csv")
Delta <- read.csv("D:/all_project/Italy_project/manuscript/figures/updated_figure2/updated/Delta_Italy_bind_esc_socre1.csv")
Gamma <- read.csv("D:/all_project/Italy_project/manuscript/figures/updated_figure2/updated/Gamma_Italy_bind_esc_socre1.csv")
Omicron <- read.csv("D:/all_project/Italy_project/manuscript/figures/updated_figure2/updated/Omicron_Italy_bind_esc_socre1.csv")


Alpha$VOC <- "Alpha"
Beta$VOC <- "Beta"
Delta$VOC<- "Delta"
Gamma$VOC <- "Gamma"
Omicron$VOC <- "Omicron"


library(purrr)
file_list <- list.files(path = "D:/all_project/Italy_project/manuscript/figures/updated_figure2/Omicron_accession_ID", 
                        pattern = "*.csv")  # 获取当前目录下所有CSV文件的文件名
path_list <- paste0("D:/all_project/Italy_project/manuscript/figures/updated_figure2/Omicron_accession_ID/",file_list)
data_list <- map(path_list, read.csv,header = F)   # 使用map()函数逐个读取CSV文件并存储为数据框列表

# 使用 sub 函数去掉 .csv 后缀
filenames_no_csv <- sub("\\.csv$", "", file_list)
names(data_list) <- filenames_no_csv

# 使用 lapply 将每个data.frame中的name添加为一列
my_list_with_name <- lapply(names(data_list), function(name) {
  df <- data_list[[name]]
  df$name <- name  # 添加name列
  names(df)[1] <- "Accession_ID"
  return(df)
})

# 将list转换为单个data.frame
omicron_id_lineage <- do.call(rbind, my_list_with_name)
omicron_id_lineage <- omicron_id_lineage %>%
  mutate(VOCs = paste0("Omicron(", gsub("^(BA\\.[12345])(\\..*)?$", "\\1*", omicron_id_lineage$name), ")" ) ) %>%
  select(-name)

Omicron1 <- Omicron %>% 
  left_join(omicron_id_lineage, by = "Accession_ID")  %>%
  select(-VOC) %>% 
  rename(VOC = VOCs) %>% #drop_na(VOCs) %>% 
  mutate(VOC = gsub("Omicron\\(BA\\.[345]\\*[,]?\\*?\\)", "Omicron(BA.3-5*)", VOC) ) %>% 
  mutate(VOC = replace_na(VOC, "Others"))

data.plot <- rbind(Alpha, Beta, Delta, Gamma, Omicron1)
data.plot$date <- as.Date(data.plot$date)

data.plot1_1 <- data.plot %>% 
  filter(binding_score > -4 & escape_score < 1) %>% 
  filter(VOC != "Beta" & VOC !="Gamma")

#table(data.plot1$VOC)
data.plot1.1 <- rbind(Alpha, Beta, Delta, Gamma, Omicron)
data.plot1.1$date <- as.Date(data.plot1.1$date)

data.plot1.1 <- data.plot1.1 %>% 
  filter(binding_score > -4 & escape_score < 1) %>% 
  filter(VOC != "Beta" & VOC !="Gamma")
#calculate weekly mean values of each voc binding and immune escape 
data.plot1.1$week <- floor_date(data.plot1.1$date, "week")
week_avg <- data.plot1.1 %>% 
  group_by(week, VOC) %>%
  summarise(binding_score_week_avg = mean(binding_score, na.rm = TRUE),
            escape_score_week_avg = mean(escape_score, na.rm = TRUE))

week_avg_voc_Alpha <- week_avg %>% filter(VOC == 'Alpha' ) %>% 
  filter(week >= ymd("2020-11-01") & week <= ymd("2021-09-01"))

week_avg_voc_Delta <- week_avg %>% filter(VOC == 'Delta' ) %>% 
  filter(week >= ymd("2021-01-01") & week <= ymd("2022-03-01"))

week_avg_voc_Omicorn <- week_avg %>% filter( VOC == "Omicron" )
week_avg_voc_total <- rbind(week_avg_voc_Alpha, week_avg_voc_Delta, week_avg_voc_Omicorn)


####################################################################################################################
##################################Figure 2A plot binding_score#######################################################################################
####################################################################################################################

startDate = ymd("2020-02-15")
endDate = ymd("2022-06-15")

data.plot1_1 <- data.plot1_1 %>% filter(VOC != "Others")

figure2A <- ggplot(data.plot1_1) +
  geom_point(aes(x = date, y = binding_score, shape = VOC, color = VOC), size = 2, alpha = 0.3) +
  scale_shape_discrete(solid = F)+
  geom_line(data = week_avg_voc_total, aes(x = week, y = binding_score_week_avg, group = VOC), linewidth = 0.9, color = "gray40")+
  scale_x_date(date_labels = "%Y-%m-%d", date_breaks = "2 months", 
               limits = c(startDate,endDate), expand = c(0.02,0.02)) +
  scale_y_continuous(limits=c(-4,1),breaks = seq(1,-4))+
  theme_bw() + ylab("Binding") + xlab("")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5,size=14),
        axis.text.y = element_text(size = 14)) +
  theme(axis.title.y = element_text(size = 16) ) +
  theme(legend.text=element_text(size=11), legend.position =c(0.2,0.35),legend.background = element_rect(fill = NA), 
        plot.title=element_text(size=16, face="bold",hjust = 0.5))+
  #ggtitle("Daily binding of VOCs in Italy") +
  theme( panel.grid = element_blank()) +
  #guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1)))
  scale_color_manual(values=c('Omicron(BA.2*)' = "#1565C0FF", 'Omicron(BA.1*)'= "#64B5F6FF", 
                             'Omicron(B.1.1.529)' =  "black", 'Omicron(BA.3-5*)'= "#984EA3",
                             #'Others' = "#787878FF",
                             'Gamma' = "#00BFC4",
                             'Alpha' = "#F8766D", 'Delta' = "#00BA38",
                             'Beta' = "#B79F00", 'Non-VOCs' = "#F564E3"))

####################################################################################################################
###########################################################Figure 2B plot immune escape#############################
####################################################################################################################
figure2B <- ggplot(data.plot1_1) +
  geom_point(aes(x = date, y = escape_score, shape = VOC, color = VOC), size = 2, alpha = 0.3) +
  scale_shape_discrete(solid = F)+
  geom_line(data = week_avg_voc_total, aes(x = week, y = escape_score_week_avg, group = VOC), linewidth = 0.9, color = "gray40") +
  scale_x_date(date_labels = "%Y-%m-%d", date_breaks = "2 months", 
               limits = c(startDate,endDate), expand = c(0.02,0.02)) +
  scale_y_continuous(limits=c(0,1),breaks = seq(0,1,0.2)) + 
  theme_bw() + ylab("Immune escape") +xlab("")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5,size=14),
        axis.text.y = element_text(size = 14)) +
  theme(axis.title.y = element_text(size = 16) ) +
  theme(legend.text=element_text(size=11), 
        plot.title=element_text(size=16, face="bold",hjust = 0.5))+
  #ggtitle("Daily Immune escape of VOCs in Italy") +
  theme( panel.grid = element_blank()) +
  scale_color_manual(values=c('Omicron(BA.2*)' = "#1565C0FF", 'Omicron(BA.1*)'= "#64B5F6FF", 
                             'Omicron(B.1.1.529)' =  "black", 'Omicron(BA.3-5*)'= "#984EA3",
                             'Others' = "#787878FF",
                             'Gamma' = "#00BFC4",
                             'Alpha' = "#F8766D", 'Delta' = "#00BA38",
                             'Beta' = "#B79F00", 'Non-VOCs' = "#F564E3"))
# 
# f1 <-figure2A
# f2 <-figure2B + theme(legend.position="none")
# #f3 <-figure2C + theme(legend.position="none")
# 
# 
# layout <- "
# AB
# "
# 
# figure_2A_2B <- f1+f2+
#   plot_layout(design = layout, heights = c(3,3, 1)) +
#   plot_annotation(tag_levels = 'A') & 
#   theme(plot.tag = element_text(size = 21, hjust = -0.1, vjust = 0.1,face = "bold"),plot.tag.position = c(0, 0.99))
# #ggsav
# ggsave(figure_2A_2B , file='D:/all_project/Italy_project/manuscript/figures/updated_figure2/figure_2A_2B.tiff', width=16, height=8.5) 

#######################################################################################
#############################Figure 2C plot binding and immune escape#############################
#######################################################################################

# figure2C <- ggtrendline(data.plot1$binding_score, data.plot1$escape_score, model = "line2P",linecolor = "grey40", linewidth = 1.05, 
#                         CI.fill = NA,CI.color = NA,eq.x = -0.7, eq.y = 0.9,rrp.x = -0.7, rrp.y= 0.8, eSize = 4) + 
#   geom_point(data = data.plot1, aes(binding_score, escape_score,shape = VOC, color = VOC), size = 2, alpha = 0.3) + 
#   scale_shape_discrete(solid = F)+
#   scale_y_continuous(limits=c(0,1),breaks = seq(0,1,0.2))+
#   theme_bw() +
#   ylab("Immune escape") +xlab("Binding")+
#   theme(axis.text.x = element_text(size=14, hjust = 0.5, vjust = .5),
#         axis.text.y = element_text(size = 14)) +
#   theme(axis.title.y = element_text(size = 16),axis.title.x = element_text(size = 16,vjust = 15)) +
#   theme(legend.text=element_text(size=11), 
#         plot.title=element_text(size=16, face="bold",hjust = 0.5))+
#   #ggtitle("Correlation for virus binding and immune escape") +
#   theme( panel.grid = element_blank()) +
#   guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1)))

# 定义逆逻辑函数
inverse_logistic <- function(x) {
  log(x / (1 - x))
}

# 对escape_score应用逆逻辑函数
data.plot1_1$transformed_escape_score <- inverse_logistic(data.plot1_1$escape_score)
data.plot2 <- data.plot1_1 %>% filter(transformed_escape_score != -Inf)

# 使用线性模型拟合
model <- lm(transformed_escape_score ~ binding_score, data = data.plot2)
AIC(model)
# 获取模型的摘要以及系数
model_summary <- summary(model)
coefficients <- model_summary$coefficients
intercept <- coefficients["(Intercept)", "Estimate"]
slope <- coefficients["binding_score", "Estimate"]

# 计算r_value（Pearson相关系数）
r_value <- cor(data.plot2$binding_score, data.plot2$escape_score)
#cor.test(data.plot2$binding_score, data.plot2$escape_score)

data.plot2_1 <- data.plot2 %>% filter( date < endDate & date > startDate)

figure2C <- ggplot(data.plot2_1, aes(x = binding_score, y = escape_score)) +
  geom_point(data = data.plot2_1, aes(binding_score, escape_score,shape = VOC, color = VOC), size = 2, alpha = 0.3) +
  scale_shape_discrete(solid = F)+
  stat_smooth(method = "glm", method.args = list(family = "binomial"), se = F, linewidth = 0.9, color = "gray40") +
  scale_y_continuous(limits=c(0,1),breaks = seq(0,1,0.2))+
  theme_bw() +
  ylab("Immune escape") +xlab("Binding")+
  scale_color_manual(values=c('Omicron(BA.2*)' = "#1565C0FF", 'Omicron(BA.1*)'= "#64B5F6FF", 
                              'Omicron(B.1.1.529)' =  "black", 'Omicron(BA.3-5*)'= "#984EA3",
                              'Others' = "#787878FF",
                              'Gamma' = "#00BFC4",
                              'Alpha' = "#F8766D", 'Delta' = "#00BA38",
                              'Beta' = "#B79F00", 'Non-VOCs' = "#F564E3"))+
  theme(axis.text.x = element_text(size=14, hjust = 0.5, vjust = .5),
        axis.text.y = element_text(size = 14)) +
  theme(axis.title.y = element_text(size = 16),axis.title.x = element_text(size = 16,vjust = 15)) +
  theme(legend.text=element_text(size=11),
        plot.title=element_text(size=16, face="bold",hjust = 0.5))+
  #ggtitle("Correlation for virus binding and immune escape") +
  theme( panel.grid = element_blank()) +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1)))+
  annotate("text", x = -0.7, y = 0.8, 
           label = paste("R-value:", round(r_value, 2), 
                         "\n y = 1/[(1+exp(", abs(round(intercept, 2)), 
                         "+", abs(round(slope, 2)), "x)]") ) 


################################################################################################
##############################Figure 2D voc time-varying prevalence############################
################################################################################################

Alpha_Italy_bind_esc_socre <- read.csv("D:/all_project/Italy_project/manuscript/figures/updated_figure2/updated/Alpha_Italy_bind_esc_socre1.csv")
Beta_Italy_bind_esc_socre <- read.csv("D:/all_project/Italy_project/manuscript/figures/updated_figure2/updated/Beta_Italy_bind_esc_socre1.csv")
Delta_Italy_bind_esc_socre <- read.csv("D:/all_project/Italy_project/manuscript/figures/updated_figure2/updated/Delta_Italy_bind_esc_socre1.csv")
Gamma_Italy_bind_esc_socre <- read.csv("D:/all_project/Italy_project/manuscript/figures/updated_figure2/updated/Gamma_Italy_bind_esc_socre1.csv")
Omicron_Italy_bind_esc_socre <- read.csv("D:/all_project/Italy_project/manuscript/figures/updated_figure2/updated/Omicron_Italy_bind_esc_socre1.csv")
Non_voc_Italy_bind_esc_socre <- read.csv("D:/all_project/Italy_project/manuscript/figures/updated_figure2/updated/Non_voc_Italy_bind_esc_socre1.csv")

Alpha_Italy_bind_esc_socre$VOCs <- "Alpha"
Beta_Italy_bind_esc_socre$VOCs <- "Beta"
Delta_Italy_bind_esc_socre$VOCs <- "Delta"
Gamma_Italy_bind_esc_socre$VOCs <- "Gamma"
Omicron_Italy_bind_esc_socre$VOCs <- "Omicron"

Non_voc_Italy_bind_esc_socre$VOCs <- ifelse( Non_voc_Italy_bind_esc_socre$date < as.Date("2020-09-01"),"Non-VOCs",
                                             ifelse(Non_voc_Italy_bind_esc_socre$date > as.Date("2020-09-01") & Non_voc_Italy_bind_esc_socre$binding_score < 0.01,"Non-VOCs","Alpha"  ) )
# count_variant <- rbind(Alpha_Italy_bind_esc_socre, Beta_Italy_bind_esc_socre, 
#                        Delta_Italy_bind_esc_socre,Gamma_Italy_bind_esc_socre,
#                        Omicron_Italy_bind_esc_socre, Non_voc_Italy_bind_esc_socre ) %>%
#   mutate(date = ymd(date)) %>% 
#   mutate(month_index = floor_date(date,"1 month")) %>% 
#   select(month_index, VOCs) %>% group_by(month_index, VOCs) %>% 
#   summarise(N = n())

Omicron_Italy_bind_esc_socre_new <- Omicron_Italy_bind_esc_socre %>% 
  left_join(omicron_id_lineage, by = "Accession_ID")  %>%
  select(-VOCs.x) %>% 
  rename(VOCs = VOCs.y) %>% #drop_na(VOCs) %>% 
  mutate(VOCs = gsub("Omicron\\(BA\\.[345]\\*[,]?\\*?\\)", "Omicron(BA.3-5*)", VOCs) ) %>% 
  mutate(VOCs = replace_na(VOCs, "Others"))

data.plot2 <- rbind(Alpha_Italy_bind_esc_socre, Beta_Italy_bind_esc_socre, 
                    Delta_Italy_bind_esc_socre,Gamma_Italy_bind_esc_socre,
                    Omicron_Italy_bind_esc_socre_new) %>% filter(date > as.Date("2021-02-01")) %>% 
  rbind(Non_voc_Italy_bind_esc_socre) %>% mutate(date = ymd(date))

data.plot2$date <- as.Date(data.plot2$date)

data_count_voc <- data.plot2 %>% 
  #filter(date <= as.Date("2022-10-01")) %>% 
  mutate(month_index = floor_date(date,"1 month")) %>%
  group_by(month_index,VOCs) %>%
  summarize(Count = n()) 

startDate = ymd("2020-03-01")
endDate = ymd("2022-06-01")

figure2D <- ggplot(data = data_count_voc, aes(x=month_index, y = Count, fill=VOCs)) + 
  #annotate(geom = "rect",xmin=startDate,xmax=endDate,ymin=0,ymax=Inf,fill="#F39B7F99",alpha=0.1)+
  #geom_col(position = position_stack(reverse = F)) + 
  geom_bar(stat = "identity", position = "fill")+
  labs(y = "% of VOCs") + xlab("")+
  scale_x_date(date_labels = "%Y-%m-%d", date_breaks = "2 months", 
               limits = as.Date(c(startDate,endDate)), expand = c(0,0))+
  # scale_fill_discrete(labels = c("Dose 1", "Dose 2"))+
  #scale_y_continuous(limits=c(0,3000), breaks = seq(0,3000,500))+
  scale_y_continuous(labels = scales::percent_format(),expand = c(0,0) )+
  # scale_fill_discrete(breaks=c('Non-VOCs', 'Gamma','Alpha', 'Delta', 'Beta', 'Omicron(B.1.1.529)', 'Omicron(BA.1*)', 
  #                              'Omicron(BA.2*)', 'Omicron(BA.3*/4*/5*)', 'Others')) +  
  scale_fill_manual(values=c('Omicron(BA.2*)' = "#1565C0FF", 'Omicron(BA.1*)'= "#64B5F6FF", 
                             'Omicron(B.1.1.529)' =  "black", 'Omicron(BA.3-5*)'= "#984EA3",
                             'Others' = "#787878FF",
                             'Gamma' = "#00BFC4",
                             'Alpha' = "#F8766D", 'Delta' = "#00BA38",
                             'Beta' = "#B79F00", 'Non-VOCs' = "#F564E3"),
                    breaks=c('Non-VOCs', 'Gamma','Alpha', 'Delta', 'Beta', 'Omicron(B.1.1.529)', 'Omicron(BA.1*)', 
                             'Omicron(BA.2*)', 'Omicron(BA.3-5*)', 'Others')) +
  #scale_fill_manual(values=c('dose 1' = "#E64B35FF","dose 2" = "#F39B7F99", "1 dose (with previous infection)" = "#4DBBD5FF",
  #"dose 3" = "#00A087FF", "dose 4" = "#3C5488FF"))+, margin = margin(0,1,0,0,'line')
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5,size=14),
        axis.text.y = element_text(size = 14)) +
  theme(axis.title.y = element_text(size = 16) ) +
  theme(legend.text=element_text(size=9), legend.position = "top",legend.key.size = unit(0.25, 'cm'),
        plot.title=element_text(size=14, face="bold",hjust = 0.5),panel.grid=element_blank())

####################################################################################################################
##########################Figure 2E binding and immune escape in RBD sites##########################################################
####################################################################################################################


RBD_deep_scan_data <- read.csv("D:/all_project/Italy_project/manuscript/figures/updated_figure2/bind_and_escape_data/RBD_deep_scan_data.csv")
ref_seq <- "MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHVSGTNGTKRFDNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLGVYYHKNNKSWMESEFRVYSSANNCTFEYVSQPFLMDLEGKQGNFKNLREFVFKNIDGYFKIYSKHTPINLVRDLPQGFSALEPLVDLPIGINITRFQTLLALHRSYLTPGDSSSGWTAGAAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQDVNCTEVPVAIHADQLTPTWRVYSTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICASYQTQTNSPRRARSVASQSIIAYTMSLGAENSVAYSNNSIAIPTNFTISVTTEILPVSMTKTSVDCTMYICGDSTECSNLLLQYGSFCTQLNRALTGIAVEQDKNTQEVFAQVKQIYKTPPIKDFGGFNFSQILPDPSKPSKRSFIEDLLFNKVTLADAGFIKQYGDCLGDIAARDLICAQKFNGLTVLPPLLTDEMIAQYTSALLAGTITSGWTFGAGAALQIPFAMQMAYRFNGIGVTQNVLYENQKLIANQFNSAIGKIQDSLSSTASALGKLQDVVNQNAQALNTLVKQLSSNFGAISSVLNDILSRLDKVEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRASANLAATKMSECVLGQSKRVDFCGKGYHLMSFPQSAPHGVVFLHVTYVPAQEKNFTTAPAICHDGKAHFPREGVFVSNGTHWFVTQRNFYEPQIITTDNTFVSGNCDVVIGIVNNTVYDPLQPELDSFKEELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNEVAKNLNESLIDLQELGKYEQYIKWPWYIWLGFIAGLIAIVMVTIMLCCMTSCCSCLKGCCSCGSCCKFDEDDSEPVLKGVKLHYT"
ref_seq_RBD <- paste0(unlist(strsplit(ref_seq, split = ""))[331:531], collapse = '')
dic_mut <- substr(RBD_deep_scan_data$mutation, start =2, stop = 6)
msa_results_copy <- readRDS( paste0(paste0("D:/all_project/Italy_project/manuscript/figures/updated_figure2/European_alg_ouput/", "Italy_alg"), '.rds') )

###immune escape calculation
imu_data <- read.csv('D:/all_project/Italy_project/manuscript/figures/updated_figure2/bind_and_escape_data/imm_eascape.csv')

imu_escape_data <- imu_data %>%
  tidyr::unite(col = "mut_type", c(wildtype, site, mutation), sep = "", remove = T) %>%
  group_by(mut_type) %>%
  summarise(
    mean_eas_score = mean(mut_escape)
  )
dic_eas <- substr(imu_escape_data$mut_type, start =2, stop = 6)

mut_index <- c(  "G339D", "R346K", "R346T", "A348S", "S371F", "S371L", "S373P", "S375F", "T376I", "T376A", "P384L", "D405N", "R408S", "Q414K", "K417N",
                 "T430S", "N439K", "N440K", "G446V", "G446S", "L452R", "L455F", "A475S", "S477N", "T478K", "P479S", "E484K", "E484Q", "E484A", "F486V",
                 "F490L", "Q493R", "G496S", "Q498R", "N501Y", "Y505H", "A520S", "A522P", "A522S")
data_mut <- data.frame(num = 1:length(mut_index), mutation = mut_index)

data_binding <- RBD_deep_scan_data[RBD_deep_scan_data$mutation %in% mut_index, ] %>% 
  select(mutation, bind_avg)
data_imm_esc <- imu_escape_data[imu_escape_data$mut_type %in% mut_index, ] %>% 
  rename(mutation = 'mut_type')

data_plot_RBD <- data_mut %>% left_join(data_binding, by = "mutation") %>%
  left_join(data_imm_esc, by = "mutation")
#write.csv(data_plot_RBD, "D:/all_project/Italy_project/manuscript/2024-08-17/data_plot_RBD.csv" ) 

#write.csv(data_plot, "D:/all_project/Italy_project/manuscript/codes/data_mut_bind_imm_point.csv")

# 将数据转换为长格式
data_long <- pivot_longer(data_plot_RBD, 
                          cols = c(bind_avg, mean_eas_score),
                          names_to = "Measurement", 
                          values_to = "value")
# 确保mutation为因子并按正确顺序排序
data_long$mutation <- factor(data_long$mutation, levels = unique(data_long$mutation))

# 使用ggplot2绘制图表，使用facet_grid按measurement分面
figure2E <- ggplot(data_long, aes(x = mutation, y = value, color = Measurement)) +
  geom_line(aes(group = Measurement )) +  # 添加这行来绘制折线
  geom_point(size = 3, shape = 16) +theme_bw()+
  facet_grid(Measurement~., scales = 'free_y',
             labeller=labeller(Measurement = c(bind_avg = "Binding affinity", 
                                               mean_eas_score = "Escape fraction")))+ 
  theme(strip.text.y = element_text(size = 14 ) )+
  #theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "", y = "") +
  theme(axis.text.x = element_blank(),axis.text.y = element_text(size = 14),legend.position="none",
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank())+
  geom_hline(data = data_long %>% filter(Measurement == "bind_avg"), aes(yintercept = 0), linetype = "dashed")

################################################################################################
################################Figure 2F Mutation prevalence########################################
################################################################################################


Alpha_data_long <- Alpha_Italy_bind_esc_socre %>%
  separate_rows(label_mut, sep = "-")

Beta_data_long <- Beta_Italy_bind_esc_socre %>%
  separate_rows(label_mut, sep = "-")

Delta_data_long <- Delta_Italy_bind_esc_socre %>%
  separate_rows(label_mut, sep = "-")

Gamma_data_long <- Gamma_Italy_bind_esc_socre %>%
  separate_rows(label_mut, sep = "-")

Omicron_data_long <- Omicron_Italy_bind_esc_socre %>%
  separate_rows(label_mut, sep = "-")

Omicron_data_long1 <- Omicron_Italy_bind_esc_socre_new %>%
  separate_rows(label_mut, sep = "-")

mut_data <- rbind(Alpha_data_long, Beta_data_long, Delta_data_long, Gamma_data_long, Omicron_data_long)


startDate = ymd("2020-10-01")
endDate = ymd("2022-10-01")

mut_data %>% 
  mutate(date = ymd(date)) %>% 
  filter(date >= startDate & date <= endDate) %>% 
  group_by(VOCs) %>%
  reframe(Count_mut = n())


Preval_mut_data <- mut_data %>% 
  mutate(date = ymd(date)) %>% 
  filter(date >= startDate & date <= endDate) %>% 
  mutate( Total_voc_mut = c(rep(25981, 25981),rep(377, 377),rep(74700, 74700), rep(6691, 6691), rep(436223, 436223)) ) %>% 
  group_by(VOCs,label_mut,Total_voc_mut) %>%
  reframe(Count_mut = n()) %>% 
  group_by(VOCs,label_mut ) %>%
  summarise(Preval_mut = Count_mut/Total_voc_mut ) %>% 
  drop_na() %>% slice(-1) %>% filter(Preval_mut > 0.0005) %>% 
  filter(VOCs == "Alpha"| VOCs =="Delta"| VOCs =="Omicron")

#write.csv(mut_data, "D:/all_project/Italy_project/manuscript/figures/updated_figure2/mut_data.csv")


# 计算每个月不同VOCs中的label_mut频率
frequency_data <- mut_data  %>%
  mutate(date = ymd(date)) %>%
  filter(date >= startDate & date <= endDate) %>%
  mutate(month = floor_date(date, "month")) %>%
  group_by(month, VOCs) %>%
  count(label_mut) %>%
  mutate(proportion = n / sum(n)) %>%
  group_by(VOCs, label_mut) %>%
  summarise(mean_preva = mean(proportion,na.rm = T)) %>%
  drop_na() %>% slice(-1) %>% filter(mean_preva > 0.001) %>%
  filter(VOCs == "Alpha"| VOCs =="Delta"| VOCs =="Omicron")


#ggtitle("S-gene mutations in lineages")
###########################################################################################################
###########################################################################################################

sample_count_month <- data_count_voc %>% group_by(month_index) %>% summarise(month_total_count = sum(Count))

#############################################updated###################################################
mut_data1 <- rbind(Alpha_data_long, Beta_data_long, Delta_data_long, Gamma_data_long, Omicron_data_long1)

frequency_data1 <- mut_data1  %>% 
  mutate(date = ymd(date)) %>% 
  filter(date <= endDate) %>%
  mutate(month_index = floor_date(date, "1 month")) %>% 
  group_by(month_index) %>%
  count(label_mut) %>% drop_na() %>% left_join(sample_count_month) %>% mutate(month_mut_freq = n/month_total_count)

#"VOCs  label_mut Preval_mut mutation_simplified"  

Alpha_period <- frequency_data1 %>% filter(month_index <= ymd("2021-05-01") & month_index >= ymd("2021-02-01")) %>% 
  group_by(month_index, label_mut) %>% summarise(month_mut_freq = mean(month_mut_freq)) %>% mutate(VOCs = "Alpha")
Delta_period <- frequency_data1 %>% filter(month_index == ymd("2021-08-01")) %>% mutate(VOCs = "Delta")
BA.1_period <- frequency_data1 %>% filter(month_index == ymd("2022-01-01")) %>% mutate(VOCs = "BA.1*")
BA.2_period <- frequency_data1 %>% filter(month_index == ymd("2022-04-01")) %>% mutate(VOCs = "BA.2*")

Preval_mut_data1 <- rbind(Alpha_period, Delta_period, BA.1_period, BA.2_period)
Preval_mut_data1 <- Preval_mut_data1[Preval_mut_data1$label_mut %in% unique(Preval_mut_data$label_mut), ]

mutation_simplified1 <- Preval_mut_data1 %>% pull(label_mut) %>% 
  unique()
mutation_simplified1 <- mutation_simplified1[order(as.numeric(gsub("\\D", "", mutation_simplified1)))]
mutation_simplified1

VOCs <- Preval_mut_data1 %>% pull(VOCs) %>% unique()
blank <- crossing(VOCs, mutation_simplified1)
blank$mutation_simplified1 <- factor(blank$mutation_simplified1, 
                                     levels = mutation_simplified1)
Preval_mut_data1$mutation_simplified1 <- factor(Preval_mut_data1$label_mut, 
                                                levels = mutation_simplified1)

Preval_mut_data1$VOCs <- factor(Preval_mut_data1$VOCs, 
                                levels = c( "Alpha", "Delta", "BA.1*", "BA.2*"))

MUTATIONPALETTE = c("#fff7f3",  "#fa9fb5","#f768a2",
                    "#f5428a")
lightBorders = F
borderColour = ifelse(lightBorders, "#FFFFFF", "#555555")

figure2F <- ggplot(Preval_mut_data1, aes(x =mutation_simplified1, y = VOCs, 
                                         fill = month_mut_freq)) + geom_tile(colour = borderColour, 
                                                                             fill = "#FFFFFF", data = blank) + geom_tile(colour = borderColour) + 
  theme_minimal()  + xlab("") + labs(fill = "Prevalence")+
  scale_fill_gradientn(colours = MUTATIONPALETTE, 
                       limits = c(0, 1), labels = scales::percent) + 
  #labs(caption = "Enabled by data from GISAID (https://gisaid.org/)") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 14), panel.grid = element_blank(), 
        legend.position = "bottom", axis.text.y = element_text(size = 14),axis.title.y = element_text(size = 16),
        plot.caption = element_text(size = 18))# + 


#############################################updated###################################################
# data_plot_mut_freq <- frequency_data1[frequency_data1$label_mut %in% mutation_simplified1,]
# ggplot(data = data_plot_mut_freq, aes(x = month_index, y = month_mut_freq, group= label_mut, color=label_mut)) + 
#   geom_line()
# 
# Alpha_data <- data_plot_mut_freq %>% filter(month_index == ymd("2021-04-01")) %>% mutate(VOCs = "Alpha")
# Delta_data <- data_plot_mut_freq %>% filter(month_index == ymd("2021-09-01")) %>% mutate(VOCs = "Delta")
# Omicron_data <- data_plot_mut_freq %>% filter(month_index == ymd("2022-03-01")) %>% mutate(VOCs = "Omicron")
# 
# Preval_mut_data1 <- rbind(Alpha_data, Delta_data, Omicron_data)
# Preval_mut_data1 <- Preval_mut_data1[Preval_mut_data1$label_mut %in% unique(Preval_mut_data$label_mut), ]
# 
# mutation_simplified1 <- Preval_mut_data1 %>% pull(label_mut) %>% 
#   unique()
# mutation_simplified1 <- mutation_simplified1[order(as.numeric(gsub("\\D", "", mutation_simplified1)))]
# mutation_simplified1
# 
# VOCs <- Preval_mut_data1 %>% pull(VOCs) %>% unique()
# blank <- crossing(VOCs, mutation_simplified1)
# blank$mutation_simplified1 <- factor(blank$mutation_simplified1, 
#                                      levels = mutation_simplified1)
# Preval_mut_data1$mutation_simplified1 <- factor(Preval_mut_data$label_mut, 
#                                                 levels = mutation_simplified1)
# 
# MUTATIONPALETTE = c("#fff7f3",  "#fa9fb5", 
#                     "#f768a1")
# lightBorders = F
# borderColour = ifelse(lightBorders, "#FFFFFF", "#555555")
# 
# figure2F <- ggplot(Preval_mut_data1, aes(x = mutation_simplified1, y = VOCs, 
#                                          fill = Preval_mut)) + geom_tile(colour = borderColour, 
#                                                                          fill = "#FFFFFF", data = blank) + geom_tile(colour = borderColour) + 
#   theme_minimal()  + xlab("") + labs(fill = "Prevalence")+
#   scale_fill_gradientn(colours = MUTATIONPALETTE, 
#                        limits = c(0, 1), labels = scales::percent) + 
#   #labs(caption = "Enabled by data from GISAID (https://gisaid.org/)") + 
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.9, hjust = 1, size = 14), panel.grid = element_blank(), 
#         legend.position = "bottom", axis.text.y = element_text(size = 14),axis.title.y = element_text(size = 16),
#         plot.caption = element_text(size = 18))# + 

##################################Combine all figure##################################

f1 <-figure2A
f2 <-figure2B + theme(legend.position="none")
f3 <-figure2C + theme(legend.position="none")
f4 <-figure2D
f5 <-figure2E
f6 <-figure2F

layout <- "
AB
CD
EE
ff
"

figure_2A_2B_2C_2D_2E_2F <- f4+f1+f2+f3+f5+f6+
  plot_layout(design = layout, heights = c(1,1,1.15,0.35)) +
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 21, hjust = -0.1, vjust = 0.1,face = "bold"),plot.tag.position = c(0, 0.99)) 
#ggsave(figure_2A_2B_2C_2D_2E_2F, file='D:/all_project/Italy_project/manuscript/figures/updated_figure2/figure_2A_2B_2C_2D_2E_2F11.pdf', width=16, height=13) 
ggsave(figure_2A_2B_2C_2D_2E_2F, file='D:/all_project/Italy_project/manuscript/figures/updated_figure2/figure_2A_2B_2C_2D_2E_2F_10_11.tiff', width=12, height=16) 

# 
# library(export)
# 
# graph2office(  x=figure_2A_2B_2C_2D_2E_2F,#需要输出的图形
#                file="D:/all_project/Italy_project/manuscript/figures/updated_figure2/figure_2A_2B_2C_2D_2E_2F11_updated12",#输出后图形的名字
#                type = c("PPT"),#输出图形的格式
#                width = 16,#图形在宽度
#                height = 13#图形在高度
# )














