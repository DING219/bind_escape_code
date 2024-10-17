
rm(list=ls())

library(lubridate)
library(ggplot2)
library(patchwork)
library(cowplot)

model_fit_res <- read.csv("D:/all_project/Italy_project/manuscript/figures/updated_figure3/model_fit_res.csv")
model_fit_res$date <- ymd(model_fit_res$date)

# beta_v_res <- read.csv("D:/all_project/Italy_project/manuscript/figures/updated_figure3/beta_v_res.csv")
# beta_v_res$date <- ymd(beta_v_res$date)
beta_v_res_updated <- read.csv("D:/all_project/Italy_project/manuscript/figures/updated_figure3/beta_v_res_updated.csv")
beta_v_res_updated$date <- ymd(beta_v_res_updated$date)
#write.csv(beta_v_res, "D:/all_project/Italy_project/manuscript/figures/updated_figure4/output_07_12_2024/fitness_supp.csv")

effective_immunity_res <- read.csv("D:/all_project/Italy_project/manuscript/figures/updated_figure3/effective_immunity_res.csv")
effective_immunity_res$date <- ymd(effective_immunity_res$date)

beta_v_res$beta_v <- beta_v_res$beta_v*(1 - effective_immunity_res$eff_imm)
population_immunity <- read.csv("D:/all_project/Italy_project/manuscript/figures/updated_figure3/population_immunity_res_updated.csv")
population_immunity$date <- ymd(population_immunity$date)
population_immunity_res <- data.frame(date = population_immunity$date, pop_imm_without_rec = population_immunity$aju_3,  pop_imm = population_immunity$aju_4 )


##############################################base plot version###########################

library(lubridate)
library(ggplot2)
library(patchwork)
library(cowplot)


p1 <- ~{
  par(mar = c(0, 6, 3, 5))
  date_index <- which(seq(ymd("2020-02-21"), ymd('2022-03-05'), by='1 day')  %in% seq(ymd("2020-03-01"), ymd('2022-03-05'), by='2 month') == T)
  plot(model_fit_res$new_cases,  cex.axis = 1.2,cex.lab = 1.7,
       lty = 1, col=rgb(0.4,0.4,0.8,0.6), xlab = "", xaxt="n",pch = 19,cex=0.7,yaxt="n", ylim = c(0,300000),
       ylab = "Daily incidence\n(thousands)"#, main = "Model fitting for daily confirmed cases (Italy)"
  )
  lines(model_fit_res$fit_mean, col = "#E64B35FF", pch = 19, lwd=1.3)
  polygon(c(rev(1:744), 1:744), c(rev(model_fit_res$fit_0.025), model_fit_res$fit_0.975), 
          col = rgb(0.8,0.7,0.7,0.6) , border = NA)
  #axis(1, at = date_index,labels = model_fit_res$date[date_index], las = 2, cex.axis = 1.2)
  axis(2, at = seq(0,300000, by = 50000),labels = seq(0,300, by = 50), las = 2, cex.axis = 1.2)
  abline(v=date_index,lty=2,col="lightgray")
  abline(h=seq(0,300000, by = 50000),lty=2,col="lightgray")
  abline(v=which(seq(ymd("2020-02-21"), ymd('2022-03-05'), by='1 day') %in% ymd("2020-12-15") == T),lty=5,col="grey", lwd = 3)
  abline(v=which(seq(ymd("2020-02-21"), ymd('2022-03-05'), by='1 day') %in% ymd("2021-05-15") == T),lty=5,col="grey", lwd = 3)
  abline(v=which(seq(ymd("2020-02-21"), ymd('2022-03-05'), by='1 day') %in% ymd("2021-12-02") == T),lty=5,col="grey", lwd = 3)
  text(220, 280000,  "Non-VOCs",cex = 1.3, pos = 2, col = "grey")
  text(420, 280000,  "Alpha",cex = 1.3, pos = 2, col = "grey")
  text(590, 280000,  "Delta",cex = 1.3, pos = 2, col = "grey")
  text(782, 280000,  "Omicron",cex = 1.3, pos = 2, col = "grey")
  # abline(v=which(seq(ymd("2020-02-21"), ymd('2022-03-05'), by='1 day') %in% ymd("2020-11-27") == T),lty=6,col='#1E90FF', lwd = 2)
  # text(500,250000,  "1 and 2 doses",cex = 0.9, pos = 2, col = '#1E90FF')
  # abline(v=which(seq(ymd("2020-02-21"), ymd('2022-03-05'), by='1 day') %in% ymd("2021-09-13") == T),lty=6,col='#1E90FF', lwd = 2)
  # text(740,250000,  "3 doses",cex = 0.9, pos = 2, col = '#1E90FF')
}



p2 <- ~{
  par(mar = c(0, 6, 3, 5))
  date_index <- which(seq(ymd("2020-02-21"), ymd('2022-03-05'), by='1 day')  %in% seq(ymd("2020-03-01"), ymd('2022-03-05'), by='2 month') == T)
  plot(4*beta_v_res_updated$beta_v/mean(beta_v_res_updated$beta_v[1:28]),  cex.axis = 1.2,cex.lab = 1.7,
       ylim = c(0,12), las = 2,yaxt = "n",
       col = "#E64B35FF", xlab = "", xaxt="n", cex=0.7,lwd=3, lty = 1,
       ylab = "",type = "l" #, main = expression(paste("Fitness"," (",beta[v], ")")) 
  )
  lines(4*beta_v_res_updated$Re, col = rgb(0.7,0.7,0.7,0.6), pch = 19,lty = 2, lwd=3)
  #lines(silded_Rt_mean, col = 2, pch = 19, lwd=2)
  #polygon(c(rev(1:744), 1:744), c(rev(silded_beta_v_0.975), silded_beta_v_0.025), col = rgb(0.7,0.7,0.7,0.6) , border = NA)
  #axis(1, at = date_index,labels = model_fit_res$date[date_index], las = 2, cex.axis = 1.2)
  axis(2, at = seq(0,12, by = 2),labels = seq(0,12, by = 2), las = 2, cex.axis = 1.2, col = "#E64B35FF")
  axis(4, at = seq(0,12, by = 2),labels = seq(0,12, by = 2), las = 2, cex.axis = 1.2, col = rgb(0.7,0.7,0.7,0.6))
  mtext("ACE2-dependent Virus\nfitness",side=2,line=2.5, col = "#E64B35FF", cex = 1.7)  
  mtext("Effective reproduction number",side=4,line=2.5, col = rgb(0.7,0.7,0.7,0.6), cex = 1.7)  
  abline(v=date_index,lty=2,col="lightgray")
  abline(h=seq(0,12, by = 2),lty=2,col="lightgray")
  abline(v=which(seq(ymd("2020-02-21"), ymd('2022-03-05'), by='1 day') %in% ymd("2020-12-15") == T),lty=5,col="grey", lwd = 3)
  abline(v=which(seq(ymd("2020-02-21"), ymd('2022-03-05'), by='1 day') %in% ymd("2021-05-15") == T),lty=5,col="grey", lwd = 3)
  abline(v=which(seq(ymd("2020-02-21"), ymd('2022-03-05'), by='1 day') %in% ymd("2021-12-02") == T),lty=5,col="grey", lwd = 3)
  text(220, 11.2,  "Non-VOCs",cex = 1.3, pos = 2, col = "grey")
  text(420, 11.2,  "Alpha",cex = 1.3, pos = 2, col = "grey")
  text(590, 11.2,  "Delta",cex = 1.3, pos = 2, col = "grey")
  text(782, 11.2,  "Omicron",cex = 1.3, pos = 2, col = "grey")
}

p3 <- ~{
  par(mar = c(9, 6, 3, 5))
  date_index <- which(seq(ymd("2020-02-21"), ymd('2022-03-05'), by='1 day')  %in% seq(ymd("2020-03-01"), ymd('2022-03-05'), by='2 month') == T)
  plot(population_immunity_res$pop_imm,  cex.axis = 1.2,cex.lab = 1.7,ylim = c(0,1), las = 2,
       col = "#E64B35FF", xlab = "", xaxt="n", cex=0.7,lwd=3, lty = 1,yaxt = "n",
       ylab = "Immunity\n(against wild type (%))", type = "l"#, main = "Vaccine induced immunity"
  )
  lines(population_immunity_res$pop_imm_without_rec, col = "#00A087FF", pch = 19, lwd=3)
  #polygon(c(rev(1:744), 1:744), c(rev(silded_beta_v_0.975), silded_beta_v_0.025), col = rgb(0.7,0.7,0.7,0.6) , border = NA)
  legend(0,0.4, legend = c("Total", "Vaccine-induced"),
         col =  c("#E64B35FF", "#00A087FF" ), bty="n", lwd = 2)
  axis(1, at = date_index,labels = model_fit_res$date[date_index], las = 2, cex.axis = 1.2,cex.lab = 1.8)
  axis(2, at = seq(0,1, by = 0.2),labels = seq(0,100, by = 20), las = 2, cex.axis = 1.2)
  abline(v=date_index,lty=2,col="lightgray")
  abline(h=seq(0,1, by = 0.2),lty=2,col="lightgray")
  abline(v=which(seq(ymd("2020-02-21"), ymd('2022-03-05'), by='1 day') %in% ymd("2020-12-15") == T),lty=5,col="grey", lwd = 3)
  abline(v=which(seq(ymd("2020-02-21"), ymd('2022-03-05'), by='1 day') %in% ymd("2021-05-15") == T),lty=5,col="grey", lwd = 3)
  abline(v=which(seq(ymd("2020-02-21"), ymd('2022-03-05'), by='1 day') %in% ymd("2021-12-02") == T),lty=5,col="grey", lwd = 3)
  text(220, 0.93,  "Non-VOCs",cex = 1.3, pos = 2, col = "grey")
  text(420, 0.93,  "Alpha",cex = 1.3, pos = 2, col = "grey")
  text(590, 0.93,  "Delta",cex = 1.3, pos = 2, col = "grey")
  text(782, 0.93,  "Omicron",cex = 1.3, pos = 2, col = "grey")
  mtext("Date", side = 1, line = 7, cex = 1.8)
}


p4 <- ~{
  par(mar = c(9, 6, 3, 5))
  date_index <- which(seq(ymd("2020-02-21"), ymd('2022-03-05'), by='1 day')  %in% seq(ymd("2020-03-01"), ymd('2022-03-05'), by='2 month') == T)
  plot(effective_immunity_res$eff_imm,  cex.axis = 1.2,cex.lab = 1.7,ylim = c(0,1), las = 2,
       col = "#E64B35FF", xlab = "", xaxt="n", cex=0.7,lwd=3, lty = 1,yaxt = "n",
       ylab = "Effective immunity\n(against VOCs (%))", type = "l" #, main = "Effective immunity"
  )
  #lines(silded_Rt_mean, col = 2, pch = 19, lwd=2)
  #polygon(c(rev(1:744), 1:744), c(rev(silded_beta_v_0.975), silded_beta_v_0.025), col = rgb(0.7,0.7,0.7,0.6) , border = NA)
  axis(1, at = date_index,labels = model_fit_res$date[date_index], las = 2, cex.axis = 1.2)
  axis(2, at = seq(0,1, by = 0.2),labels = seq(0,100, by = 20), las = 2, cex.axis = 1.2)
  abline(v=date_index,lty=2,col="lightgray")
  abline(h=seq(0,1, by = 0.2),lty=2,col="lightgray")
  abline(v=which(seq(ymd("2020-02-21"), ymd('2022-03-05'), by='1 day') %in% ymd("2020-12-15") == T),lty=5,col="grey", lwd = 3)
  abline(v=which(seq(ymd("2020-02-21"), ymd('2022-03-05'), by='1 day') %in% ymd("2021-05-15") == T),lty=5,col="grey", lwd = 3)
  abline(v=which(seq(ymd("2020-02-21"), ymd('2022-03-05'), by='1 day') %in% ymd("2021-12-02") == T),lty=5,col="grey", lwd = 3)
  text(220, 0.93,  "Non-VOCs",cex = 1.3, pos = 2, col = "grey")
  text(420, 0.93,  "Alpha",cex = 1.3, pos = 2, col = "grey")
  text(590, 0.93,  "Delta",cex = 1.3, pos = 2, col = "grey")
  text(782, 0.93,  "Omicron",cex = 1.3, pos = 2, col = "grey")
  mtext("Date", side = 1, line = 7, cex = 1.8)
}

figure_3A_3B_3C_3D <- plot_grid(p1, p2, p3, p4, rel_heights = c(1, 1.4),scale = 1, 
                                labels = "AUTO",label_size = 20,  align = "hv", axis = "tb")
#ggsave(figure_3A_3B_3C_3D, file='D:/all_project/Italy_project/manuscript/figures/updated_figure3/figure_3A_3B_3C_3D.pdf', width=15, height=10.5) 
ggsave(figure_3A_3B_3C_3D, file='D:/all_project/Italy_project/manuscript/figures/updated_figure3/figure_3A_3B_3C_3D2_updated_10_10.tiff', width=16, height=11.5) 











