library(ggplot2)
library(tidyverse)

cont_var <- read.csv("D:/all_project/Italy_project/manuscript/figures/suppl/count_variant.csv")


all_mut <- cont_var %>% separate_rows(label_mut, sep = "-") %>% select(-binding_score, -escape_score)


mut_index <- c(  "G339D", "R346K", "R346T", "A348S", "S371F", "S371L", "S373P", "S375F", "T376I", "T376A", "P384L", "D405N", "R408S", "Q414K", "K417N",
                 "T430S", "N439K", "N440K", "G446V", "G446S", "L452R", "L455F", "A475S", "S477N", "T478K", "P479S", "E484K", "E484Q", "E484A", "F486V",
                 "F490L", "Q493R", "G496S", "Q498R", "N501Y", "Y505H", "A520S", "A522P", "A522S")  

# 
# weekly_counts <- cont_var %>% na.omit() %>% filter(label_mut != "") %>%  
#   mutate(date = ymd(date)) %>% 
#   mutate(week_index = floor_date(date,"1 week")) %>%
#   group_by(week_index, label_mut) %>%
#   summarise(count = n())

weekly_counts <- all_mut %>% na.omit() %>% filter(label_mut != "") %>%  
  mutate(date = ymd(date)) %>% 
  mutate(week_index = floor_date(date,"1 week")) %>%
  group_by(week_index, label_mut) %>%
  summarise(count = n())

weekly_counts <- weekly_counts[weekly_counts$label_mut %in% mut_index,]

weekly_totals <- cont_var %>% na.omit() %>% filter(label_mut != "")  %>%  
  mutate(date = ymd(date)) %>% 
  mutate(week_index = floor_date(date,"1 week")) %>%
  group_by(week_index) %>%
  summarise(total_count = n())

# 合并每周突变次数和总样本数，计算频率
weekly_freq <- merge(weekly_counts, weekly_totals, by = "week_index")
weekly_freq <- weekly_freq %>%
  mutate(frequency = count / total_count)
weekly_freq$frequency[1:16] <- 0

# frequency_data1 <- mut_data  %>%  na.omit() %>% filter(label_mut != "") %>% 
#   mutate(date = ymd(date)) %>% 
#   filter(date <= endDate) %>%
#   mutate(week_index = floor_date(date, "1 week")) %>% 
#   group_by(week_index) %>%
#   count(label_mut) %>% drop_na() %>% left_join(sample_count_week) %>% mutate(week_mut_freq = n/week_total_count)
# 
# frequency_data1 <- frequency_data1[frequency_data1$label_mut %in% mut_index,]

ggplot(weekly_freq, aes(x = week_index, y = frequency, color = label_mut)) +
  geom_line() +
  scale_x_date(date_breaks = "2 month",expand = c(0,0), 
               limits = as.Date(c("2020-03-01","2022-06-15"))) +
  scale_y_continuous(limits=c(0,1), breaks = seq(0,1,0.2)) +   
  #geom_hline(yintercept = 100, lty=5, colour="grey", linewidth = 0.6)+ rosybrown1
  theme_bw() + ylab("Weekly mutation prevalence (%)") + xlab("")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5,size=12),
        axis.text.y = element_text(size = 12)) +
  theme(axis.title.y = element_text(size = 14) ) +
  theme(legend.text=element_text(size=11), legend.position = "none",
        plot.title=element_text(size=16, face="bold",hjust = 0.5))+
  theme( panel.grid = element_blank())+  
  geom_vline(xintercept = ymd("2021-02-03"), lty=5,colour="grey", linewidth = 0.5)+
  geom_vline(xintercept = ymd("2021-06-26"), lty=5,colour="grey", linewidth = 0.5)+
  geom_vline(xintercept = ymd("2021-09-28"), lty=5,colour="grey", linewidth = 0.5)+
  annotate("text", x=ymd("2020-08-10"), y=0.9, label= "Non-VOCs",col="grey", size=3.5)+
  annotate("text", x=ymd("2021-04-15"), y=0.9, label= "Alpha",col="grey", size=3.5)+
  annotate("text", x=ymd("2021-08-12"), y=0.9, label= "Delta",col="grey", size=3.5)+
  annotate("text", x=ymd("2022-01-01"), y=0.9, label= "Omicron",col="grey", size=3.5)+
  geom_vline(xintercept = ymd("2021-08-20"), lty=5,colour="rosybrown1", linewidth = 0.8)+
  geom_vline(xintercept = ymd("2021-11-06"), lty=5,colour="rosybrown1", linewidth = 0.8)+
  annotate("text", x=ymd("2021-09-28"), y=0.5, label= "Valley",col="rosybrown1", size=3.5)
# 
# weekly_freq1 <- weekly_freq
# 
# grep("E484A|S477N|G339D|Y505H|Q498R|S373P|S375F|K417N|Q493R", weekly_freq$label_mut, perl =T)
# weekly_freq %>% mutate(label_mut1 = grep("E484A|S477N|G339D|Y505H|Q498R|S373P|S375F|K417N|Q493R", label_mut, perl =T))

# set1 <- c("E484A","S477N","G339D","Y505H","Q498R", "S373P","S375F","K417N")
# set2 <- c("G446S","S371L","G496S")
# set3 <- c("S371F","T376A","D405N","R408S")

mut_index1 <- c( "R346K", "R346T", "A348S", "S371F", "T376I", "P384L", "Q414K", 
                 "T430S", "N439K", "N440K", "G446V", "G446S", "L452R", "L455F", 
                 "A475S", "T478K", "P479S", "E484K", "E484Q", "E484A", "F486V",
                 "F490L", "Q493R", "N501Y", "A520S", "A522P", "A522S")  

#weekly_counts <- weekly_counts[weekly_counts$label_mut %in% mut_index,]
weekly_freq <- weekly_freq[weekly_freq$label_mut %in% mut_index1,]

weekly_freq$label_mut[grep("E484A", weekly_freq$label_mut, perl =T)] <- "Set1"
weekly_freq$label_mut[grep("G446S", weekly_freq$label_mut, perl =T)] <- "Set2"
weekly_freq$label_mut[grep("S371F", weekly_freq$label_mut, perl =T)] <- "Set3"

weekly_freq_slice <- weekly_freq %>% group_by(label_mut) %>% arrange(week_index) %>% 
  group_by(label_mut) %>% slice_head(n = 1) %>% arrange(week_index)

weekly_freq$label_mut <- factor(weekly_freq$label_mut, levels = weekly_freq_slice$label_mut )


p1 <- ggplot(weekly_freq, aes(x = week_index, y = frequency)) +
  geom_line(color = "red", linewidth = 1) +
  scale_x_date(date_breaks = "2 month",expand = c(0.01,0), 
               limits = as.Date(c("2020-12-01","2022-05-10"))) + 
  scale_y_continuous(limits=c(0,1), breaks = seq(0,1,0.2)) +   
  #geom_hline(yintercept = 100, lty=5, colour="grey", linewidth = 0.6)+ rosybrown1
  theme_bw() + ylab("Mutation prevalence") + xlab("")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5,size=12),
        axis.text.y = element_text(size = 12)) +
  theme(axis.title.y = element_text(size = 14) ) +
  theme(legend.text=element_text(size=11), legend.position = "none",
        plot.title=element_text(size=16, face="bold",hjust = 0.5))+
  theme( panel.grid = element_blank())+  
  #geom_vline(xintercept = ymd("2020-12-15"), lty=5,colour="grey", linewidth = 0.5)+
  geom_vline(xintercept = ymd("2021-05-15"), lty=5,colour="grey", linewidth = 0.5)+
  geom_vline(xintercept = ymd("2021-12-10"), lty=5,colour="grey", linewidth = 0.5)+
  #annotate("text", x=ymd("2020-08-10"), y=0.9, label= "Non-VOCs",col="grey", size=3.5)+
  annotate("text", x=ymd("2021-03-01"), y=0.9, label= "Alpha",col="grey", size=3.5)+
  annotate("text", x=ymd("2021-08-20"), y=0.9, label= "Delta",col="grey", size=3.5)+
  annotate("text", x=ymd("2022-02-20"), y=0.9, label= "Omicron",col="grey", size=3.5)+
  annotate("rect", xmin=ymd("2021-05-01"), xmax=ymd("2021-06-27"), ymin=-Inf, ymax=Inf,fill='lightgoldenrod2', alpha = .3)+
  annotate("rect", xmin=ymd("2021-07-31"), xmax=ymd("2021-11-18"), ymin=-Inf, ymax=Inf,fill='#33CCCC', alpha = .3)+
  annotate("rect", xmin=ymd("2021-12-01"), xmax=ymd("2021-12-20"), ymin=-Inf, ymax=Inf,fill='seagreen3', alpha = .35)+
  annotate("rect", xmin=ymd("2021-12-27"), xmax=ymd("2022-02-01"), ymin=-Inf, ymax=Inf,fill='seagreen1', alpha = .3)+
  facet_wrap(~label_mut, ncol=4, scales = "free_y")

p1 
  # geom_vline(xintercept = ymd("2021-07-31"), lty=5,colour="green", linewidth = 0.3)+
  # geom_vline(xintercept = ymd("2021-11-18"), lty=5,colour="green", linewidth = 0.3)+
  # #annotate("text", x=ymd("2021-10-01"), y=0.5, label= "Valley1",col="rosybrown1", size=3.5)+
  # geom_vline(xintercept = ymd("2021-12-01"), lty=5,colour="green", linewidth = 0.3)+
  # geom_vline(xintercept = ymd("2021-12-15"), lty=5,colour="green", linewidth = 0.3)+
  # #annotate("text", x=ymd("2021-12-07"), y=0.5, label= "Valley2",col="rosybrown1", size=3.5)
  # geom_vline(xintercept = ymd("2021-12-23"), lty=5,colour="green", linewidth = 0.3)+
  # geom_vline(xintercept = ymd("2022-02-10"), lty=5,colour="green", linewidth = 0.3) + 
  # facet_wrap(~label_mut  , ncol=4, scales = "free_y")
  #  
# p1
# library(eoffice)
# topptx(p1,"D:/all_project/Italy_project/manuscript/figures/suppl/mutation_pre_by_wrap_face3.pptx",width = 18,height = 12)

 
#ggplot_build(p)$data
#write.csv(weekly_freq, "D:/all_project/Italy_project/manuscript/figures/updated_figure2/weekly_freq.csv")

weekly_freq_main <- weekly_freq %>% filter(label_mut == "E484K" | label_mut == "L452R" | label_mut == "T478K" | 
                         label_mut == "Set1" | label_mut == "Set2" | label_mut == "Set3")

# 将label_mut转换为因子，并按所需顺序排列
weekly_freq_main$label_mut <- factor(weekly_freq_main$label_mut, 
                                     levels = c("E484K", "L452R", "T478K", "Set1", "Set2", "Set3"))


p2 <- ggplot(weekly_freq_main, aes(x = week_index, y = frequency)) +
  geom_line(color = "red", linewidth = 1) +
  scale_x_date(date_breaks = "2 month",expand = c(0.01,0), 
               limits = as.Date(c("2020-12-01","2022-05-10"))) + 
  scale_y_continuous(limits=c(0,1), breaks = seq(0,1,0.2)) +   
  #geom_hline(yintercept = 100, lty=5, colour="grey", linewidth = 0.6)+ rosybrown1
  theme_bw() + ylab("") + xlab("")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5,size=12),
        axis.text.y = element_text(size = 12)) +
  theme(axis.title.y = element_text(size = 14) ) +
  theme(legend.text=element_text(size=11), legend.position = "none",
        plot.title=element_text(size=16, face="bold",hjust = 0.5))+
  theme( panel.grid = element_blank())+  
  #geom_vline(xintercept = ymd("2020-12-15"), lty=5,colour="grey", linewidth = 0.5)+
  geom_vline(xintercept = ymd("2021-05-15"), lty=5,colour="grey", linewidth = 0.5)+
  geom_vline(xintercept = ymd("2021-12-10"), lty=5,colour="grey", linewidth = 0.5)+
  #annotate("text", x=ymd("2020-08-10"), y=0.9, label= "Non-VOCs",col="grey", size=3.5)+
  annotate("text", x=ymd("2021-03-01"), y=0.9, label= "Alpha",col="grey", size=3.5)+
  annotate("text", x=ymd("2021-08-20"), y=0.9, label= "Delta",col="grey", size=3.5)+
  annotate("text", x=ymd("2022-02-20"), y=0.9, label= "Omicron",col="grey", size=3.5)+
  annotate("rect", xmin=ymd("2021-05-01"), xmax=ymd("2021-06-27"), ymin=-Inf, ymax=Inf,fill='lightgoldenrod2', alpha = .3)+
  annotate("rect", xmin=ymd("2021-07-31"), xmax=ymd("2021-11-18"), ymin=-Inf, ymax=Inf,fill='#33CCCC', alpha = .3)+
  annotate("rect", xmin=ymd("2021-12-01"), xmax=ymd("2021-12-20"), ymin=-Inf, ymax=Inf,fill='seagreen3', alpha = .35)+
  annotate("rect", xmin=ymd("2021-12-27"), xmax=ymd("2022-02-01"), ymin=-Inf, ymax=Inf,fill='seagreen1', alpha = .3)+
  facet_wrap(~label_mut  , ncol=3, scales = "free_y")

p2




