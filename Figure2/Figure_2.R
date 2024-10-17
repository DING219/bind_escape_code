

rm(list=ls())

library(tidyverse)
library(RColorBrewer)
library(cowplot)
library(ggrepel)
library(lubridate)
library(scales)
library(ggpubr)

source("D:/all_project/Italy_project/0_data_process/input/data_uesd.R")
###Vaccine data processing
startDate = ymd("2020-12-27")
endDate = ymd("2022-03-05")
index_date = ymd("2020-01-31")

vaccine_data_raw <- read.csv("C:/Users/USER01/Documents/Zhao-personal/important_literature/vacinee_imunity/The effect of COVID-19 vaccination in Italy and/covid19-opendata-vaccini-master/dati/somministrazioni-vaccini-latest.csv")

names(vaccine_data_raw)
# vaccine_data_Italy <- vaccine_data_raw %>%
#   mutate(data = as.Date(data)) %>%
#   select(-area,-N1,-N2,-ISTAT, -reg) %>%
#   filter(forn != "Novavax" & forn !="Sanofi") %>%
#   mutate(vacc_type = if_else(forn == "Janssen" | forn == "Vaxzevria (AstraZeneca)","Viral_vector_vacc" , "mRNA_vacc")) %>%
#   pivot_longer(cols = d1:db3, names_to = "Dose",values_to = "Values_for_each_dose") %>%
#   pivot_longer(cols = m:f, names_to = "Gender",values_to = "Values_for_each_gender")
# 
# 
# data_paper <- vaccine_data_raw %>%
#   mutate(data = as.Date(data)) %>%
#   filter(data >= startDate & data <= as.Date("2021-07-27"))%>%
#   select(data, d1:d2) %>%
#   group_by(data) %>% 
#   summarise(d1 = sum(d1), d2 = sum(d2))%>% 
#   pivot_longer(d1:d2,names_to = "Dose",values_to = "Count")
# 
# data_paper %>% ggplot(aes(x=data, y = Count/1e6, fill=Dose, inherit.aes = F)) + 
#   #annotate(geom = "rect",xmin=startDate,xmax=endDate,ymin=0,ymax=Inf,fill="#F39B7F99",alpha=0.1)+
#   geom_col(position = position_stack(reverse = TRUE) )

vaccine_data_Italy_doses <- vaccine_data_raw %>% 
  mutate(data = as.Date(data)) %>% 
  filter(forn != "Novavax" & forn !="Sanofi") %>% 
  mutate(vacc_type = if_else(forn == "Janssen" | forn == "Vaxzevria (AstraZeneca)","Viral_vector_vacc" , "mRNA_vacc")) %>% 
  group_by(data) %>% 
  summarise(d1 = sum(d1), d2 = sum(d2),
            dpi = sum(dpi),db1 = sum(db1),
            db2 = sum(db2),db3 = sum(db3))

vaccine_data_Italy_doses_tot <- rbind(data.frame(data = seq(as.Date("2020-02-21"), as.Date("2020-12-26"), by = 1),
                                                 d1 = rep(0, length(seq(as.Date("2020-02-21"), as.Date("2020-12-26"), by = 1))),
                                                 d2 = rep(0, length(seq(as.Date("2020-02-21"), as.Date("2020-12-26"), by = 1))),
                                                 dpi = rep(0, length(seq(as.Date("2020-02-21"), as.Date("2020-12-26"), by = 1))),
                                                 db1 = rep(0, length(seq(as.Date("2020-02-21"), as.Date("2020-12-26"), by = 1))),
                                                 db2 = rep(0, length(seq(as.Date("2020-02-21"), as.Date("2020-12-26"), by = 1))),
                                                 db3 = rep(NA, length(seq(as.Date("2020-02-21"), as.Date("2020-12-26"), by = 1))) ), vaccine_data_Italy_doses)

####B mobility index


mob_data <- data.frame(date = merged_mob$date, mobility_index = merged_mob_ave)
mob_data <- mob_data %>% filter(date <= ymd("2022-03-05") & date >= ymd("2020-02-21"))

figure1B <- ggplot(mob_data, aes(x = date)) + 
  geom_line( aes( y = mobility_index*100), color = "#6A8B21", linewidth = 0.7 ) + 
  scale_x_date(date_breaks = "2 month",expand = c(0.031,0.031)) +
  scale_y_continuous(limits=c(0,200), breaks = seq(0,200,40)) +   
  geom_hline(yintercept = 100, lty=5, colour="grey", linewidth = 0.6)+
  theme_bw() + ylab("Daily mobility\nindex (%)") + xlab("")+
  geom_vline(xintercept = ymd("2020-12-15"), lty=5,colour="grey", linewidth = 0.8)+
  geom_vline(xintercept = ymd("2021-05-15"), lty=5,colour="grey", linewidth = 0.8)+
  geom_vline(xintercept = ymd("2021-12-02"), lty=5,colour="grey", linewidth = 0.8)+
  annotate("text", x=ymd("2020-07-01"), y=188, label= "Non-VOCs",col="grey", size=5)+
  annotate("text", x=ymd("2021-03-01"), y=188, label= "Alpha",col="grey", size=5)+
  annotate("text", x=ymd("2021-08-20"), y=188, label= "Delta",col="grey", size=5)+
  annotate("text", x=ymd("2022-02-01"), y=188, label= "Omicron",col="grey", size=5)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5,size=14),
        axis.text.y = element_text(size = 14)) +
  theme(axis.title.y = element_text(size = 20) ) +
  theme(legend.text=element_text(size=13), legend.position = c(0.1,0.8),
        plot.title=element_text(size=18, face="bold",hjust = 0.5))+
  theme( panel.grid = element_blank())


Italy_reported_cases$new_cases[120] <- mean(Italy_reported_cases$new_cases[c(117:119,121:123)])
Italy_reported_cases$Italy_cu_sum_cases <- cumsum(Italy_reported_cases$new_cases)
#Italy_reported_cases
#write.csv(Italy_reported_cases, "D:/all_project/Italy_project/manuscript/Italy_reported_cases.csv")
####C daily reported cases
Italy_reported_cases$log_Italy_cu_sum_cases <- log10(Italy_reported_cases$Italy_cu_sum_cases)

figure1C <- ggplot(Italy_reported_cases) +
  #geom_point(aes(x = date, y = new_cases), size = 3, alpha = 0.5, color = "#6666CC99") +
  geom_line(aes(x = date, y = new_cases), linewidth = 0.4, color="#E64B35FF")+
  scale_x_date(date_breaks = "2 month",expand = c(0.031,0.031)) + 
  scale_y_continuous(transform = log10_trans(),
                     breaks = trans_breaks("log10", function(x) 10^x) (c( 10, 100, 1000,1e4,1e5,1e6, 1e7, 1e8)) ,
                     labels = trans_format("log10", math_format(10^.x) ), limits=c(10, 1e8) ,
                     sec.axis = sec_axis(~.,
                                         name = 'Cumulative Cases\n(log10 scale)', 
                                         breaks = trans_breaks("log10", function(x) 10^x)(c(10, 100, 1000, 1e4, 1e5, 1e6, 1e7, 1e8)),
                                         labels = trans_format("log10", math_format(10^.x)) )
  ) + 
  geom_line(aes(x = date, y = Italy_cu_sum_cases ),linewidth = 0.7, color = "#E8BF80") +
  theme_bw() +
  geom_vline(xintercept = ymd("2020-12-15"), lty=5,colour="grey", linewidth = 0.8)+
  geom_vline(xintercept = ymd("2021-05-15"), lty=5,colour="grey", linewidth = 0.8)+
  geom_vline(xintercept = ymd("2021-12-02"), lty=5,colour="grey", linewidth = 0.8)+
  labs(x = "", y = "Daily incidence\n(log10 scale)",width=1)+
  annotate("text", x=ymd("2020-07-01"), y=10^7.6, label= "Non-VOCs",col="grey", size=5)+
  annotate("text", x=ymd("2021-03-01"), y=10^7.6, label= "Alpha",col="grey", size=5)+
  annotate("text", x=ymd("2021-08-20"), y=10^7.6, label= "Delta",col="grey", size=5)+
  annotate("text", x=ymd("2022-02-01"), y=10^7.6, label= "Omicron",col="grey", size=5)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5,size=14),
        axis.text.y = element_text(size = 16, color="#E64B35FF")) +
  theme(axis.title.y = element_text(size = 20, color="#E64B35FF"),axis.text.y.right = element_text(color = "#E8BF80"),
        axis.title.y.right = element_text(size = 20,color = "#E8BF80")) +
  theme(legend.text=element_text(size=13), legend.position =c(0.2,0.3),
        plot.title=element_text(size=18, face="bold",hjust = 0.5))+
  theme( panel.grid = element_blank() ) 

####A vaccine coverage
vaccine_data_Italy_doses_tot <- data.frame(vaccine_data_Italy_doses_tot %>% filter(data <= ymd("2022-03-05") & data >= ymd("2020-02-21")))

figure1A <- ggplot(vaccine_data_Italy_doses_tot, aes(x = data)) + 
  geom_line( aes( y = 100*cumsum(d1)/59500579 , color = "d1"), linewidth = 0.7 ) + 
  geom_line( aes( y = 100*cumsum(d2)/59500579 , color = "d2") , linewidth = 0.7) + 
  geom_line( aes( y = 100*cumsum(db1)/59500579, color = "db1") , linewidth = 0.7) +
  scale_x_date(date_breaks = "2 month",expand = c(0.031,0.031), limits=c(as.Date("2020-02-21"),as.Date("2022-04-01"))) +
  scale_y_continuous(limits=c(0,100), breaks = seq(0,100,20)) +   theme_bw() + ylab("Cumulative vaccine\ncoverage (%)") + xlab("")+
  geom_vline(xintercept = ymd("2020-12-15"), lty=5,colour="grey", linewidth = 0.8)+
  geom_vline(xintercept = ymd("2021-05-15"), lty=5,colour="grey", linewidth = 0.8)+
  geom_vline(xintercept = ymd("2021-12-02"), lty=5,colour="grey", linewidth = 0.8)+
  annotate("text", x=ymd("2020-07-01"), y=94, label= "Non-VOCs",col="grey", size=5)+
  annotate("text", x=ymd("2021-03-01"), y=94, label= "Alpha",col="grey", size=5)+
  annotate("text", x=ymd("2021-08-20"), y=94, label= "Delta",col="grey", size=5)+
  annotate("text", x=ymd("2022-02-01"), y=94, label= "Omicron",col="grey", size=5)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5,size=14),
        axis.text.y = element_text(size = 14)) +
  theme(axis.title.y = element_text(size = 20) ) +
  theme(legend.text=element_text(size=9), legend.position = c(0.2,0.4), 
        legend.background = element_rect(fill = "transparent"),legend.spacing.y = unit(0.1, 'pt'),
        legend.title=element_text(size=10))+
  theme( panel.grid = element_blank()) +
  scale_colour_manual(name = "Dose", values = c("d1" = "#3C5487", "d2" =  "#00A088","db1"="#E64B34" ),
                      labels = c('1 dose', '2 doses', '3 doses') )


fig_cowplot <- plot_grid(figure1C, figure1A, figure1B, rel_heights = c(1, 1), 
                         align = 'hv', 
                         axis = "lr")
fig_cowplot 
ggsave(fig_cowplot, file='D:/all_project/Italy_project/manuscript/figures/updated_figure1/figure_2_10_10_supply.tiff', width=16, height=12) 


####supplementary XY plot for daily incidence and mobility

# data_mob_inc <- mob_data %>% merge(Italy_reported_cases, by = "date") %>% 
#   mutate(VOC_index = ifelse(date <= ymd("2020-12-15"), "Non-VOCs", 
#                             ifelse(date <= ymd("2021-05-15") & date > ymd("2020-12-15"), "Alpha", 
#                                    ifelse(date > ymd("2021-05-15") & date <= ymd("2021-12-10"), "Delta", "Omicron" ) )) ) %>% 
#   mutate(cases_log = log(new_cases)) %>% 
#   mutate(VOC_index = factor(VOC_index, levels = c("Non-VOCs",  "Alpha", "Delta","Omicron")) )
# 
# 
# figure_XY_mob_case <- ggplot(data_mob_inc, aes(x=mobility_index, y= cases_log, color=as.factor(VOC_index)  )) + 
#   geom_point(size=3, alpha = 0.4) + theme_bw()+
#   facet_wrap(~VOC_index , dir="h", scales = "free" ) +
#   scale_x_continuous( expand = c(0.15,0))+
#   #geom_smooth(method = "lm", se = T) +
#   stat_cor(label.y =c(3.4,6.5,5.7,7.5), r.digits = 3,p.accuracy = 0.01 )+
#   labs(x = "Mobility", y = "Daily incidence\n(log10 scale)", width=1)+
#   theme(legend.position="none",strip.text.x = element_text(size=13), 
#         axis.title.x = element_text(size=20),axis.text.x = element_text(size=13),axis.text.y = element_text(size = 13),
#         axis.title.y = element_text(size = 20),panel.grid = element_blank() )
# 
# library(patchwork)  
# 
# 
# layout <- "
# AB
# CD
# "
# 
# figure_2A_2B_2C <- figure1C + figure1B  + figure1A + plot_spacer() + 
#   plot_layout(design = layout, heights = c(2,2))
# 
# ggsave(figure_2A_2B_2C, file='D:/all_project/Italy_project/manuscript/figures/updated_figure1/figure_2_10_10_supply.tiff', width=17, height=12) 
# 
# 
# library(export)
# 
# graph2office(  x=figure_2A_2B_2C_2D,#需要输出的图形
#                file="D:/all_project/Italy_project/manuscript/figures/updated_figure2/updated/figure_2_08_27",#输出后图形的名字
#                type = c("PPT"),#输出图形的格式
#                width = 20,#图形在宽度
#                height = 12#图形在高度
# )
# 
