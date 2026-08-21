design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

design

dge <- estimateDisp(dge, design)

fit <- glmQLFit(dge, design)

contrast_drought <- makeContrasts(
  Drought - Control,
  levels = design
)

contrast_heat <- makeContrasts(
  Heat - Control,
  levels = design
)

contrast_heat_vs_drought <- makeContrasts(
  Heat - Drought,
  levels = design
)

de_drought <- glmQLFTest(fit, contrast = contrast_drought)

de_heat <- glmQLFTest(fit, contrast = contrast_heat)

de_hd <- glmQLFTest(fit, contrast = contrast_heat_vs_drought)

tab_drought <- topTags(
  de_drought,
  n = Inf
)$table

tab_heat <- topTags(
  de_heat,
  n = Inf
)$table

tab_hd <- topTags(
  de_hd,
  n = Inf
)$table

head(tab_drought)

head(tab_heat)

sig_drought <- rownames(tab_drought[tab_drought$FDR < 0.05, ])
sig_heat <- rownames(tab_heat[tab_heat$FDR < 0.05, ])
sig_hd <- rownames(tab_hd[tab_hd$FDR < 0.05, ])

candidate_genes <- unique(c(
  sig_drought,
  sig_heat,
  sig_hd
))

length(candidate_genes)


