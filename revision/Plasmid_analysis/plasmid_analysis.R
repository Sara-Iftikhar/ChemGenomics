library(tidyverse)
library(patchwork)
library(scales)

# ============================================================
# 1. SETTINGS AND FUNCTIONS
# ============================================================

base <- "/Users/daneshm/Documents/ChemicalGenomic/revision"
src <- file.path(base, "plasmid_LMM_full")
out <- file.path(src, "supplementary_material_two_models")
dir.create(out, recursive=TRUE, showWarnings=FALSE)

excluded <- "Col440II_1"

reps <- c("Col440I_1","ColRNAI_1","ColKP3_1","IncFIA(HI1)_1_HI1",
          "IncFIB(K)_1_Kpn3","IncFIB(Mar)_1_pNDM-Mar","IncFII_1_pKP91",
          "IncHI1B_1_pNDM-MAR","IncL/M(pOXA-48)_1_pOXA-48","IncR_1")

cargo7 <- c("blaOXA48","blaNDM","blaCTXM15","iuc","ybt","blaOXA232","blaKPC2")
cargo_labels <- c(blaOXA48="blaOXA-48",blaNDM="blaNDM",blaCTXM15="blaCTX-M-15",
                  iuc="iuc",ybt="ybt",blaOXA232="blaOXA-232",blaKPC2="blaKPC-2")

model_labels <- c("M1: ST-adjusted","M2: ST + cargo")
model_colors <- c("M1: ST-adjusted"="#2166AC","M2: ST + cargo"="#B2182B")

read_result <- function(...) read_csv(file.path(...),show_col_types=FALSE)
save_csv <- function(x,name) write_csv(x,file.path(out,paste0(name,".csv")))

library(ragg)

save_fig <- function(p, name, w, h) {
  
  ggsave(paste0(name, ".pdf"), p, width = w, height = h,
         
         device = cairo_pdf)
  
  ggsave(paste0(name, ".png"), p, width = w, height = h,
         
         device = ragg::agg_png, dpi = 300)
  
}

save_fig(p1, "Figure_S1_LMM", 8, 4.8)

theme_pub <- function() theme_bw(base_size=10,base_family="Arial") +
  theme(panel.grid=element_blank(),
        axis.text=element_text(size=9,color="black"),
        axis.title=element_text(size=10),
        strip.text=element_text(size=10,face="bold"),
        legend.title=element_text(size=9),
        legend.text=element_text(size=9),
        plot.title=element_text(size=11,face="bold"),
        plot.tag=element_text(size=14,face="bold"),
        plot.margin=margin(6,6,6,6))

theme_set(theme_pub())

rep_factor <- function(x) factor(x,levels=rev(reps))
model_factor <- function(x) factor(x,levels=c("M1","M2"),labels=model_labels)
safe_median <- function(x) if(any(is.finite(x))) median(x[is.finite(x)]) else NA_real_
safe_percentage <- function(x) if(any(!is.na(x))) 100*mean(x,na.rm=TRUE) else NA_real_
safe_cor <- function(x,y,method="pearson") {
  ok <- is.finite(x)&is.finite(y)
  if(sum(ok)<3 || sd(x[ok])==0 || sd(y[ok])==0) return(NA_real_)
  cor(x[ok],y[ok],method=method)
}
safe_reduction <- function(x,baseline)
  if(length(baseline)!=1 || !is.finite(baseline) || baseline<=0)
    rep(NA_real_,length(x)) else 100*(1-x/baseline)

# ============================================================
# 2. LOAD DATA AND EXCLUDE Col440II_1
# ============================================================

lmm <- read_result(src,"LMM_all_results.csv") %>%
  filter(Replicon!=excluded,Model %in% c("M1_ST","M3_Extended")) %>%
  mutate(Model=recode(Model,"M1_ST"="M1","M3_Extended"="M2"))

cooccur <- read_result(src,"replicon_cargo_cooccurrence.csv") %>%
  filter(Replicon!=excluded)

vif <- read_result(src,"replicon_VIF.csv") %>%
  filter(Replicon!=excluded,Model=="M3_Extended") %>%
  mutate(Model="M2")

matched <- read_result(src,"within_ST_matching","matched_fitness_results.csv") %>%
  filter(Replicon!=excluded)

diagnostics <- read_result(src,"within_ST_matching","matching_diagnostics.csv") %>%
  filter(Replicon!=excluded)

balance <- read_result(src,"within_ST_matching","cargo_balance.csv") %>%
  filter(Replicon!=excluded)

cn <- read_result(src,"copy_number","copy_number_data.csv") %>%
  filter(Replicon!=excluded)

cn_results <- read_result(src,"copy_number","copy_number_fitness_associations.csv") %>%
  filter(Replicon!=excluded)

stopifnot(
  all(c("Condition","Replicon","Model","Beta","SE","P","FDR","Status") %in% names(lmm)),
  setequal(unique(lmm$Model),c("M1","M2")),
  !any(duplicated(lmm[c("Condition","Replicon","Model")]))
)

save_csv(lmm,"LMM_all_results_two_models")

# ============================================================
# 3. OVERALL MODEL SUMMARY
# ============================================================

model_sum <- lmm %>%
  group_by(Model) %>%
  summarise(
    N_comparisons=n(),
    N_tested=sum(is.finite(P)),
    N_tested_FDR=sum(is.finite(FDR)),
    N_significant=sum(is.finite(FDR)&FDR<0.05),
    Percent_significant=safe_percentage(if_else(is.finite(FDR),FDR<0.05,NA)),
    Median_abs_beta=safe_median(abs(Beta)),
    Median_SE=safe_median(SE),
    .groups="drop"
  ) %>%
  arrange(Model) %>%
  mutate(
    Reduction_significant_pct=safe_reduction(N_significant,N_significant[Model=="M1"]),
    Reduction_median_abs_beta_pct=safe_reduction(Median_abs_beta,Median_abs_beta[Model=="M1"])
  )

save_csv(model_sum,"Table_LMM_overall_summary")

# ============================================================
# 4. FIGURE S1: SIGNIFICANT ASSOCIATIONS
# ============================================================

lmm_summary <- lmm %>%
  group_by(Replicon,Model) %>%
  summarise(
    N_tested=sum(is.finite(P)),
    N_significant=sum(is.finite(FDR)&FDR<0.05),
    Median_abs_beta=safe_median(abs(Beta)),
    Median_SE=safe_median(SE),
    .groups="drop"
  )

save_csv(lmm_summary,"Table_S1_LMM")

p1 <- lmm_summary %>%
  mutate(Replicon=rep_factor(Replicon),Model=model_factor(Model)) %>%
  ggplot(aes(N_significant,Replicon,fill=Model)) +
  geom_col(position=position_dodge(width=0.75),width=0.7) +
  scale_fill_manual(values=model_colors) +
  labs(x="FDR-significant fitness associations",y=NULL,fill=NULL) +
  theme(legend.position="top")

save_fig(p1,"Figure_S1_LMM",8,4.8)

# ============================================================
# 5. PAIRED MODEL COMPARISON
# ============================================================

comparison <- lmm %>%
  select(Condition,Replicon,Model,Beta,SE,FDR,Status) %>%
  pivot_wider(names_from=Model,values_from=c(Beta,SE,FDR,Status)) %>%
  mutate(
    Evaluable=is.finite(Beta_M1)&is.finite(Beta_M2),
    Direction_agreement=if_else(Evaluable,sign(Beta_M1)==sign(Beta_M2),NA),
    Significant_M1=if_else(is.finite(FDR_M1),FDR_M1<0.05,NA),
    Significant_M2=if_else(is.finite(FDR_M2),FDR_M2<0.05,NA),
    SE_ratio=if_else(Evaluable&is.finite(SE_M1)&is.finite(SE_M2)&SE_M1>0&SE_M2>0,
                     SE_M2/SE_M1,NA_real_),
    SE_change_pct=100*(SE_ratio-1),
    Absolute_beta_change=if_else(Evaluable,abs(Beta_M2)-abs(Beta_M1),NA_real_),
    Attenuation_pct=if_else(Evaluable&abs(Beta_M1)>1e-8,
                            100*(1-abs(Beta_M2)/abs(Beta_M1)),NA_real_)
  )

save_csv(comparison,"Table_LMM_paired_model_comparison")

effect_data <- comparison %>% filter(Evaluable)
if(nrow(effect_data)==0) stop("No paired finite model coefficients found.")

summarise_effects <- function(dat) dat %>%
  summarise(
    N_matched=n(),
    Pearson=safe_cor(Beta_M1,Beta_M2),
    Spearman=safe_cor(Beta_M1,Beta_M2,"spearman"),
    Direction_agreement_pct=safe_percentage(Direction_agreement),
    Median_abs_beta_M1=safe_median(abs(Beta_M1)),
    Median_abs_beta_M2=safe_median(abs(Beta_M2)),
    Median_attenuation_pct=safe_median(Attenuation_pct),
    N_evaluable_attenuation=sum(is.finite(Attenuation_pct)),
    Median_SE_M1=safe_median(SE_M1),
    Median_SE_M2=safe_median(SE_M2),
    Median_SE_ratio=safe_median(SE_ratio),
    Median_SE_change_pct=safe_median(SE_change_pct),
    Percent_SE_increased=safe_percentage(if_else(is.finite(SE_ratio),SE_ratio>1,NA)),
    N_evaluable_SE=sum(is.finite(SE_ratio)),
    N_significant_M1=sum(Significant_M1,na.rm=TRUE),
    N_significant_M2=sum(Significant_M2,na.rm=TRUE),
    N_significant_both=sum(Significant_M1&Significant_M2,na.rm=TRUE),
    N_significant_both_concordant=sum(Significant_M1&Significant_M2&Direction_agreement,na.rm=TRUE),
    .groups="drop"
  )

baseline_hits <- effect_data %>% filter(Significant_M1 %in% TRUE)
adjusted_hits <- effect_data %>% filter(Significant_M2 %in% TRUE)
shared_hits <- effect_data %>% filter(Significant_M1 %in% TRUE,Significant_M2 %in% TRUE)

effect_summary <- summarise_effects(effect_data)
baseline_hit_summary <- summarise_effects(baseline_hits)

publication_summary <- bind_rows(
  summarise_effects(effect_data) %>% mutate(Subset="All paired associations"),
  summarise_effects(baseline_hits) %>% mutate(Subset="M1 significant"),
  summarise_effects(adjusted_hits) %>% mutate(Subset="M2 significant"),
  summarise_effects(shared_hits) %>% mutate(Subset="Both models significant")
) %>% relocate(Subset)

save_csv(effect_summary,"Table_LMM_effect_concordance")
save_csv(baseline_hit_summary,"Table_LMM_baseline_significant_hits")
save_csv(publication_summary,"Table_LMM_publication_summary")

# ============================================================
# 6. SIGNIFICANCE OVERLAP
# ============================================================

significance_summary <- comparison %>%
  filter(!is.na(Significant_M1),!is.na(Significant_M2)) %>%
  summarise(
    N_comparable=n(),
    M1_only=sum(Significant_M1&!Significant_M2),
    M2_only=sum(!Significant_M1&Significant_M2),
    Both_significant=sum(Significant_M1&Significant_M2),
    Both_significant_concordant=sum(Significant_M1&Significant_M2&Direction_agreement,na.rm=TRUE),
    Neither_significant=sum(!Significant_M1&!Significant_M2),
    N_M1_significant=sum(Significant_M1),
    N_M2_significant=sum(Significant_M2),
    Retention_pct=if(N_M1_significant>0) 100*Both_significant/N_M1_significant else NA_real_
  )

rep_effect_summary <- effect_data %>% group_by(Replicon) %>% summarise_effects()
rep_baseline_summary <- baseline_hits %>% group_by(Replicon) %>% summarise_effects()
retained_hits <- comparison %>% filter(Significant_M2 %in% TRUE) %>% arrange(FDR_M2)

save_csv(significance_summary,"Table_LMM_significance_overlap")
save_csv(rep_effect_summary,"Table_LMM_replicon_effect_concordance")
save_csv(rep_baseline_summary,"Table_LMM_replicon_baseline_hits")
save_csv(retained_hits,"Table_LMM_cargo_adjusted_significant_hits")

# ============================================================
# 7. FIGURE S2: COEFFICIENT CONCORDANCE
# ============================================================

cor_label <- sprintf(
  "n = %s\nPearson r = %.3f\nSpearman rho = %.3f\nDirection agreement = %.1f%%",
  comma(effect_summary$N_matched),effect_summary$Pearson,
  effect_summary$Spearman,effect_summary$Direction_agreement_pct
)

effect_limits <- range(c(effect_data$Beta_M1,effect_data$Beta_M2),finite=TRUE)

p2 <- ggplot(effect_data,aes(Beta_M1,Beta_M2)) +
  geom_point(alpha=0.30,size=0.8,color="#2166AC") +
  geom_hline(yintercept=0,color="grey70",linewidth=0.3) +
  geom_vline(xintercept=0,color="grey70",linewidth=0.3) +
  geom_abline(slope=1,intercept=0,linetype="dashed",linewidth=0.4) +
  annotate("text",x=-Inf,y=Inf,label=cor_label,hjust=-0.08,vjust=1.2,size=3,family="Arial") +
  coord_equal(xlim=effect_limits,ylim=effect_limits) +
  labs(x="ST-adjusted coefficient (M1)",y="Cargo-adjusted coefficient (M2)")

save_fig(p2,"Figure_S2_Effect_Concordance",5.5,5)

# ============================================================
# 8. EFFECT ATTENUATION AND PRECISION
# ============================================================

p_attenuation <- baseline_hits %>%
  filter(is.finite(Attenuation_pct)) %>%
  mutate(Replicon=rep_factor(Replicon)) %>%
  ggplot(aes(Attenuation_pct,Replicon)) +
  geom_boxplot(width=0.6,fill="#D1E5F0",outlier.alpha=0.3) +
  geom_vline(xintercept=0,linetype="dashed") +
  labs(x="Absolute coefficient attenuation (%)",y=NULL)

p_precision <- baseline_hits %>%
  filter(is.finite(SE_ratio),SE_ratio>0) %>%
  mutate(Replicon=rep_factor(Replicon)) %>%
  ggplot(aes(SE_ratio,Replicon)) +
  geom_boxplot(width=0.6,fill="#FDDBC7",outlier.alpha=0.3) +
  geom_vline(xintercept=1,linetype="dashed") +
  scale_x_log10() +
  labs(x="Standard-error ratio (M2/M1)",y=NULL)

save_fig(
  (p_attenuation|p_precision)+plot_annotation(tag_levels="A"),
  "Figure_LMM_Attenuation_and_Precision",12,5
)

# ============================================================
# 9. FIGURE S3: ISOLATE-LEVEL CO-OCCURRENCE
# ============================================================

cooccur_plot <- cooccur %>%
  filter(Replicon %in% reps,Cargo %in% cargo7) %>%
  mutate(
    Replicon=rep_factor(Replicon),
    Cargo=factor(Cargo,levels=cargo7),
    Prevalence=P_cargo_given_replicon,
    Label=if_else(is.na(Prevalence),"NA",sprintf("%.0f",100*Prevalence))
  )

p3 <- ggplot(cooccur_plot,aes(Cargo,Replicon,fill=Prevalence)) +
  geom_tile(color="grey65",linewidth=0.3) +
  geom_text(aes(label=Label),size=2.5,family="Arial") +
  scale_fill_gradient(low="white",high="#2166AC",limits=c(0,1),
                      breaks=c(0,0.5,1),labels=c("0%","50%","100%"),na.value="grey90") +
  scale_x_discrete(labels=cargo_labels,expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  coord_fixed() +
  labs(x=NULL,y=NULL,fill="Prevalence") +
  theme(axis.text.x=element_text(angle=45,hjust=1),panel.border=element_blank())

save_fig(p3,"Figure_S3_Cooccurrence",7.5,4.8)

# ============================================================
# 10. FIGURE S4: VIF DIAGNOSTICS
# ============================================================

vif_cap <- max(c(35,vif$VIF[is.finite(vif$VIF)]))

vif_plot <- vif %>%
  mutate(
    Replicon=rep_factor(Replicon),
    VIF_plot=if_else(is.infinite(VIF),vif_cap,VIF),
    Label=case_when(is.infinite(VIF)~"Inf",is.na(VIF)~"NA",TRUE~sprintf("%.1f",VIF)),
    VIF_category=case_when(is.na(VIF)~NA_character_,VIF>=10~"High",
                           VIF>=5~"Moderate",TRUE~"Low")
  )

save_csv(vif,"Table_S4_VIF_seven_cargo")

p4 <- ggplot(vif_plot,aes(VIF_plot,Replicon,fill=VIF_category)) +
  geom_col(width=0.65,na.rm=TRUE) +
  geom_vline(xintercept=5,linetype="dashed",color="grey45",linewidth=0.4) +
  geom_vline(xintercept=10,linetype="dotted",color="grey45",linewidth=0.4) +
  geom_text(aes(label=Label),hjust=-0.2,size=3,na.rm=TRUE,family="Arial") +
  scale_fill_manual(values=c(Low="#2166AC",Moderate="#FDAE61",High="#B2182B")) +
  scale_x_continuous(expand=expansion(mult=c(0,0.15))) +
  labs(x="Variance inflation factor (VIF)",y=NULL,fill=NULL) +
  theme(legend.position="none")

save_fig(p4,"Figure_S4_VIF",7,4.8)

# ============================================================
# 11. FIGURE S5: ONT PHYSICAL LINKAGE
# ============================================================

ont_path <- file.path(base,"plasmid/plasmid_best_hits_with_Kleborate.tsv")
if(!file.exists(ont_path)) stop("ONT input file not found: ",ont_path)

ont_delim <- if(grepl("\t",readLines(ont_path,n=1))) "\t" else ","
ont <- read_delim(ont_path,delim=ont_delim,show_col_types=FALSE)

detect_gene <- function(x,pattern)
  as.integer(str_detect(replace_na(as.character(x),""),
                        regex(pattern,ignore_case=TRUE)))

ont <- ont %>%
  mutate(
    Plasmid=paste(file,SEQUENCE,sep="__"),
    blaOXA48=detect_gene(Bla_Carb_acquired,"OXA[-_]?48(?![0-9])"),
    blaNDM=detect_gene(Bla_Carb_acquired,"NDM"),
    blaCTXM15=detect_gene(Bla_ESBL_acquired,"CTX[-_]?M[-_]?15(?![0-9])"),
    blaOXA232=detect_gene(Bla_Carb_acquired,"OXA[-_]?232(?![0-9])"),
    blaKPC2=detect_gene(Bla_Carb_acquired,"KPC[-_]?2(?![0-9])"),
    iuc=as.integer(!is.na(Aerobactin)&!Aerobactin %in% c("","-","0")),
    ybt=as.integer(!is.na(Yersiniabactin)&!Yersiniabactin %in% c("","-","0"))
  )

ont_cargo <- ont %>%
  group_by(Plasmid) %>%
  summarise(across(all_of(cargo7),~as.integer(any(.x==1,na.rm=TRUE))),.groups="drop")

ont_reps <- ont %>%
  filter(GENE %in% reps) %>%
  distinct(Plasmid,Replicon=GENE)

ont_cooccur <- ont_reps %>%
  left_join(ont_cargo,by="Plasmid") %>%
  pivot_longer(all_of(cargo7),names_to="Cargo",values_to="Present") %>%
  group_by(Replicon,Cargo) %>%
  summarise(
    N_plasmids=n(),
    N_positive=sum(Present==1,na.rm=TRUE),
    Prevalence=N_positive/N_plasmids,
    .groups="drop"
  )

save_csv(ont_cooccur,"Table_S5_ONT_Physical_Linkage")

ont_plot <- ont_cooccur %>%
  complete(Replicon=reps,Cargo=cargo7) %>%
  mutate(
    Replicon=rep_factor(Replicon),
    Cargo=factor(Cargo,levels=cargo7),
    Label=if_else(is.na(Prevalence),"NA",sprintf("%.0f",100*Prevalence))
  )

p5 <- ggplot(ont_plot,aes(Cargo,Replicon,fill=Prevalence)) +
  geom_tile(color="grey65",linewidth=0.3) +
  geom_text(aes(label=Label),size=2.5,family="Arial") +
  scale_fill_gradient(low="white",high="#2166AC",limits=c(0,1),
                      breaks=c(0,0.5,1),labels=c("0%","50%","100%"),na.value="grey90") +
  scale_x_discrete(labels=cargo_labels,expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  coord_fixed() +
  labs(x=NULL,y=NULL,fill="Physical linkage") +
  theme(axis.text.x=element_text(angle=45,hjust=1),panel.border=element_blank())

save_fig(p5,"Figure_S5_ONT_Physical_Linkage",7.5,4.8)

# ============================================================
# 12. FIGURE S6: MATCHING DIAGNOSTICS
# ============================================================

match_pairs <- diagnostics %>%
  filter(Method %in% c("Exact","Propensity")) %>%
  mutate(Replicon=rep_factor(Replicon),
         ST=paste0("ST",str_remove(as.character(ST),"^ST")))

p6a <- ggplot(match_pairs,aes(Replicon,N_pairs,fill=Method)) +
  geom_col(position="dodge",width=0.7) +
  facet_wrap(~ST,nrow=1) +
  coord_flip() +
  labs(x=NULL,y="Matched pairs",fill=NULL) +
  theme(legend.position="top")

balance_plot <- balance %>%
  filter(is.finite(SMD)) %>%
  mutate(Abs_SMD=abs(SMD),Stage=factor(Stage,levels=c("Before","After")))

p6b <- ggplot(balance_plot,aes(Stage,Abs_SMD,fill=Stage)) +
  geom_boxplot(width=0.5,outlier.alpha=0.3) +
  geom_hline(yintercept=0.1,linetype="dashed",color="red") +
  facet_wrap(~Method) +
  labs(x=NULL,y="Absolute standardized mean difference") +
  theme(legend.position="none")

p6 <- (p6a/p6b)+plot_layout(heights=c(1.5,1))
save_fig(p6,"Figure_S6_Matching_Diagnostics",9,8)

# ============================================================
# 13. FIGURE S7: MATCHED FITNESS EFFECTS
# ============================================================

matched_plot <- matched %>%
  filter(is.finite(Beta),Method %in% c("Exact","Propensity")) %>%
  mutate(
    Replicon=rep_factor(Replicon),
    ST=paste0("ST",str_remove(as.character(ST),"^ST"))
  )

p7 <- ggplot(matched_plot,aes(Replicon,Beta,fill=Method)) +
  geom_boxplot(position=position_dodge(width=0.7),width=0.6,outlier.alpha=0.3) +
  geom_hline(yintercept=0,linetype="dashed") +
  facet_wrap(~ST,nrow=1) +
  coord_flip() +
  labs(x=NULL,y="Matched fitness difference",fill=NULL) +
  theme(legend.position="top")

save_fig(p7,"Figure_S7_Matched_Fitness",9,5)

# ============================================================
# 14. FIGURE S8: COPY-NUMBER DISTRIBUTION
# ============================================================

cn_iso <- cn %>%
  group_by(Isolate,Replicon) %>%
  summarise(CopyNumber=safe_median(CopyNumber),.groups="drop") %>%
  filter(is.finite(CopyNumber),CopyNumber>0)

cn_summary <- cn_iso %>%
  group_by(Replicon) %>%
  summarise(
    N_isolates=n_distinct(Isolate),
    Median_CN=median(CopyNumber),
    Q1_CN=quantile(CopyNumber,0.25),
    Q3_CN=quantile(CopyNumber,0.75),
    IQR_CN=IQR(CopyNumber),
    .groups="drop"
  )

save_csv(cn_summary,"Table_S8_Copy_Number")

p8 <- ggplot(cn_iso,aes(reorder(Replicon,CopyNumber,median),CopyNumber)) +
  geom_boxplot(width=0.55,fill="grey85",outlier.shape=NA) +
  geom_jitter(width=0.12,alpha=0.35,size=0.8,color="#2166AC") +
  scale_y_log10() +
  coord_flip() +
  labs(x=NULL,y="Chromosome-normalized copy-number proxy")

save_fig(p8,"Figure_S8_Copy_Number",8,4.8)

# ============================================================
# 15. FIGURE S9: COPY NUMBER–FITNESS ASSOCIATIONS
# ============================================================

cn_plot <- cn_results %>%
  filter(is.finite(rho)) %>%
  mutate(Replicon=rep_factor(Replicon))

cn_association_summary <- cn_results %>%
  group_by(Replicon) %>%
  summarise(
    N_evaluable=sum(is.finite(rho)),
    Median_rho=safe_median(rho),
    Median_abs_rho=safe_median(abs(rho)),
    .groups="drop"
  )

save_csv(cn_association_summary,"Table_Copy_Number_Fitness_Summary")

p9 <- ggplot(cn_plot,aes(Replicon,rho)) +
  geom_boxplot(width=0.55,fill="grey85",outlier.alpha=0.3) +
  geom_hline(yintercept=0,linetype="dashed") +
  coord_flip() +
  labs(x=NULL,y="Spearman correlation with fitness")

save_fig(p9,"Figure_S9_Copy_Number_Fitness",8,4.8)

# ============================================================
# 16. COMBINED PUBLICATION TABLE
# ============================================================

linkage_summary <- ont_cooccur %>%
  filter(is.finite(Prevalence),Prevalence>0) %>%
  mutate(
    Cargo_label=recode(Cargo,!!!cargo_labels),
    Linkage=sprintf("%s (%.0f%%)",Cargo_label,100*Prevalence)
  ) %>%
  arrange(Replicon,desc(Prevalence)) %>%
  group_by(Replicon) %>%
  summarise(Physically_linked_cargo=paste(Linkage,collapse="; "),.groups="drop")

ont_counts <- ont_reps %>%
  group_by(Replicon) %>%
  summarise(N_ONT_plasmids=n_distinct(Plasmid),.groups="drop")

publication_table <- tibble(Replicon=reps) %>%
  left_join(vif %>% select(Replicon,VIF),by="Replicon") %>%
  left_join(linkage_summary,by="Replicon") %>%
  left_join(ont_counts,by="Replicon") %>%
  left_join(cn_summary,by="Replicon") %>%
  mutate(
    Physically_linked_cargo=case_when(
      !is.na(Physically_linked_cargo)~Physically_linked_cargo,
      !is.na(N_ONT_plasmids)&N_ONT_plasmids>0~"None detected",
      TRUE~"Not evaluable"
    ),
    Copy_number_median_IQR=if_else(
      is.finite(Median_CN),
      sprintf("%.2f (%.2f-%.2f)",Median_CN,Q1_CN,Q3_CN),NA_character_
    ),
    Copy_number_median_IQR_width=if_else(
      is.finite(Median_CN),
      sprintf("%.2f (%.2f)",Median_CN,IQR_CN),NA_character_
    )
  ) %>%
  select(Replicon,VIF,N_ONT_plasmids,Physically_linked_cargo,
         Copy_number_median_IQR,Copy_number_median_IQR_width,
         N_isolates,Median_CN,Q1_CN,Q3_CN,IQR_CN)

save_csv(publication_table,"Table_S4_Replicon_cargo_copy_number")

# ============================================================
# 17. FINAL FOUR-PANEL FIGURE S10
# ============================================================

# A: Significant fitness associations
# B: Coefficient concordance
# C: ONT physical linkage
# D: Matched fitness differences

figure_S10 <- wrap_plots(
  p1,p2,
  p5,p7,
  ncol=2,
  widths=c(1,1),
  heights=c(1,1.2)
) +
  plot_annotation(tag_levels="A") &
  theme(
    text=element_text(family="Arial"),
    plot.tag=element_text(size=16,face="bold",family="Arial"),
    plot.margin=margin(8,8,8,8)
  )

save_fig(
  figure_S10,
  "Figure_S10_Replicon_fitness_and_cargo",
  16,12
)

# ============================================================
# 18. FINAL VALIDATION AND RESULTS
# ============================================================

stopifnot(
  !excluded %in% lmm$Replicon,
  !excluded %in% cooccur$Replicon,
  !excluded %in% vif$Replicon,
  !excluded %in% matched$Replicon,
  !excluded %in% diagnostics$Replicon,
  !excluded %in% balance$Replicon,
  !excluded %in% cn$Replicon,
  !excluded %in% cn_results$Replicon,
  !excluded %in% ont_cooccur$Replicon,
  !excluded %in% publication_table$Replicon
)

cat("\n========== ANALYSIS COMPLETE ==========\n")
cat("Excluded replicon:",excluded,"\n")
cat("Models:",paste(sort(unique(lmm$Model)),collapse=", "),"\n")
cat("LMM replicons:",n_distinct(lmm$Replicon),"\n")
cat("LMM conditions:",n_distinct(lmm$Condition),"\n")
cat("Unique comparisons:",nrow(distinct(lmm,Replicon,Condition)),"\n")
cat("Total model results:",nrow(lmm),"\n")
cat("Matched coefficients:",nrow(effect_data),"\n")
cat("ONT plasmids:",n_distinct(ont$Plasmid),"\n")
cat("ONT isolates:",n_distinct(ont$file),"\n")
cat("Matched ST-replicon combinations:",nrow(distinct(matched,ST,Replicon)),"\n")
cat("Copy-number isolates:",n_distinct(cn_iso$Isolate),"\n")
cat("Output directory:",out,"\n")

cat("\n========== MODEL SUMMARY ==========\n")
print(model_sum,width=Inf)

cat("\n========== EFFECT SUMMARY ==========\n")
print(publication_summary,n=Inf,width=Inf)

cat("\n========== SIGNIFICANCE OVERLAP ==========\n")
print(significance_summary,width=Inf)

cat("\n========== REPLICON-SPECIFIC COMPARISON ==========\n")
print(rep_effect_summary,n=Inf,width=Inf)

cat("\n========== BASELINE HITS ==========\n")
print(rep_baseline_summary,n=Inf,width=Inf)

cat("\n========== CARGO-ADJUSTED SIGNIFICANT HITS ==========\n")
print(retained_hits,n=Inf,width=Inf)

cat("\n========== VIF ==========\n")
print(vif,n=Inf,width=Inf)

cat("\n========== COPY-NUMBER SUMMARY ==========\n")
print(cn_summary,n=Inf,width=Inf)

cat("\n========== COPY-NUMBER FITNESS SUMMARY ==========\n")
print(cn_association_summary,n=Inf,width=Inf)

cat("\n========== PUBLICATION TABLE ==========\n")
print(publication_table,n=Inf,width=Inf)