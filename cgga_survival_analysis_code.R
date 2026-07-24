#install.packages(c("dplyr","tidyr","survival","survminer","ggplot2","gridExtra"))

library(dplyr)
library(tidyr)
library(survival)
library(survminer)
library(ggplot2)
library(gridExtra)


clin <- clinical.file
clin$OS     <- as.numeric(clin$OS)
clin$Censor <- as.numeric(clin$Censor)

expr <- expr.file
gene.col <- names(expr)[1]

#functions
get_gene_data <- function(gene){
  
  expr %>%
    filter(.data[[gene.col]] == gene) %>%
    tidyr::pivot_longer(
      cols = -all_of(gene.col),
      names_to = "CGGA_ID",
      values_to = "Expression"
    ) %>%
    inner_join(clin, by = "CGGA_ID")
  
}

filter_group <- function(data, prs_type){
  
  data %>%
    filter(
      PRS_type == prs_type,
      Grade == "WHO IV"
    )
  
}


#### EGFR primary ####


egfr.dat <- get_gene_data("EGFR")
egfr.primary <- filter_group(egfr.dat, "Primary")

med.exp <- median(egfr.primary$Expression)

n.high <- sum(egfr.primary$Expression >= med.exp)
n.low  <- sum(egfr.primary$Expression < med.exp)

egfr.primary$status <- ifelse(
  egfr.primary$Expression >= med.exp,
  paste0("High (", n.high, ")"),
  paste0("Low (", n.low, ")")
)

high_lab <- paste0("High (", n.high, ")")
low_lab  <- paste0("Low (", n.low, ")")

egfr.primary$status <- factor(
  egfr.primary$status,
  levels = c(high_lab, low_lab)
)

fit <- survfit(
  Surv(OS, Censor) ~ status,
  data = egfr.primary
)

# Graph
egfr.plot.primary <- ggsurvplot(
  fit,
  data = egfr.primary,
  pval = TRUE,
  xlab = "Time (days)",
  pval.coord = c(3500, 0.98),
  ggtheme = theme_light() +
    theme(
      plot.title = element_text(hjust = 0.5)
    ),  
  surv.median.line = "hv",
  title = "EGFR (WHO IV Primary)",
  palette = c("#FF5C00", "#6900A8"),
  legend.title = "Expression",
  legend.labs = c(high_lab, low_lab)
)

#print(egfr.plot.primary)


#### EGFR Recurrent ####

egfr.recurrent <- filter_group(egfr.dat, "Recurrent")

med.exp <- median(egfr.recurrent$Expression)

n.high <- sum(egfr.recurrent$Expression >= med.exp)
n.low  <- sum(egfr.recurrent$Expression < med.exp)

egfr.recurrent$status <- ifelse(
  egfr.recurrent$Expression >= med.exp,
  paste0("High (", n.high, ")"),
  paste0("Low (", n.low, ")")
)

high_lab <- paste0("High (", n.high, ")")
low_lab  <- paste0("Low (", n.low, ")")

egfr.recurrent$status <- factor(
  egfr.recurrent$status,
  levels = c(high_lab, low_lab)
)

fit <- survfit(
  Surv(OS, Censor) ~ status,
  data = egfr.recurrent
)

# Graph 
egfr.plot.recurrent <- ggsurvplot(
  fit,
  data = egfr.recurrent,
  pval = TRUE,
  xlab = "Time (days)",
  pval.coord = c(450, 0.98),
  ggtheme = theme_light() +
    theme(
      plot.title = element_text(hjust = 0.5)
    ), 
  surv.median.line = "hv",
  title = "EGFR (WHO IV Recurrent)",
  palette = c("#FF5C00", "#6900A8"),
  legend.title = "Expression",
  legend.labs = c(high_lab, low_lab)
)

#print(egfr.plot.recurrent)


#ggsave(filename = "mcl1.plot.recurrent.png", plot = mcl1.plot.recurrent$plot, width = 6, height = 4, dpi = 300)

