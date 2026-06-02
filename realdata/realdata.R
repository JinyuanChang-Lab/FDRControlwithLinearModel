library(knockoff)
library(glmnet)
library(scalreg)
library(doParallel)

source("gm.R")
source("Inverse_real.R")

set.seed(2024)

#Test the PI drug
drug_class = 'PI'

#Download the data and read it into data frames.

base_url = 'http://hivdb.stanford.edu/pages/published_analysis/genophenoPNAS2006'
#gene_url = paste(base_url, 'DATA', paste0(drug_class, '_DATA.txt'), sep='/')
#gene_df = read.delim(gene_url, na.string = c('NA', ''), stringsAsFactors = FALSE)

tsm_url = paste(base_url, 'MUTATIONLISTS', 'NP_TSM', drug_class, sep='/')
tsm_df = read.delim(tsm_url, header = FALSE, stringsAsFactors = FALSE)
names(tsm_df) = c('Position', 'Mutations')

#For the gene_df, read it manually
gene_df = read.delim("PI_DATA.txt", na.string = c('NA', ''), stringsAsFactors = FALSE)


#Remove error flags or nonstandard mutation codes.

grepl_rows <- function(pattern, df) {
  cell_matches = apply(df, c(1,2), function(x) grepl(pattern, x))
  apply(cell_matches, 1, all)
}

pos_start = which(names(gene_df) == 'P1')
pos_cols = seq.int(pos_start, ncol(gene_df))
valid_rows = grepl_rows('^(\\.|-|[A-Zid]+)$', gene_df[,pos_cols])
gene_df = gene_df[valid_rows,]

# Flatten a matrix to a vector with names from concatenating row/column names.
flatten_matrix <- function(M, sep='.') {
  x <- c(M)
  names(x) <- c(outer(rownames(M), colnames(M),
                      function(...) paste(..., sep=sep)))
  x
}

# Construct preliminary design matrix.
muts = c(LETTERS, 'i', 'd')
X = outer(muts, as.matrix(gene_df[,pos_cols]), Vectorize(grepl))
X = aperm(X, c(2,3,1))
dimnames(X)[[3]] <- muts
X = t(apply(X, 1, flatten_matrix))
mode(X) <- 'numeric'

# Remove any mutation/position pairs that never appear in the data.
X = X[,colSums(X) != 0]

# Extract response matrix.
Y = gene_df[,4:(pos_start-1)]


# Comparison 
Comparison <- function (X, y, q) {
  # Log-transform the drug resistance measurements.
  y = log(y)
  
  # Remove patients with missing measurements.
  missing = is.na(y)
  y = y[!missing]
  X = X[!missing,]
  
  # Remove predictors that appear less than 3 times.
  X = X[,colSums(X) >= 3]
  
  # Remove duplicate predictors.
  X = X[,colSums(abs(cor(X)-1) < 1e-4) == 1]
  
  # Original X and y
  X0 = X
  y0 = y
  
  # Centering X and y
  X = scale(X, center = TRUE, scale = FALSE)
  y = scale(y, center = TRUE, scale = FALSE)
  
  # Run Algorithm 1 and Algorithm 2 and MX knockoff
  p = ncol(X)
  n = nrow(X)
  
  # Transformation matrix
  ATrans = rbind(diag(p),diag(p))
  BTrans = rbind(diag(p),-diag(p))
  Trans = cbind(ATrans, BTrans)
  
  # Run the knockoff filter to achieve knockoff X
  KX = create.second_order(X, method = c("asdp"), shrink = F)
  
  # Augmented linear model 
  AX = cbind(X,KX)
  fit = cv.glmnet(AX, y, nlambda=100, alpha=1, standardize=T)
  Ab = as.vector(coef(fit, s=fit$lambda.1se))[-1]
  
  # Inverse directly
  G = t(AX)%*%AX/n
  Inv = InverseLinfty(G, n)
  
  # Debiased lasso
  Adb = as.numeric(Ab + (Inv%*%t(AX)%*%(y - AX %*% Ab))/n)
  w2 = Adb[1:p] - Adb[(p+1):(2*p)]
  w1 = Adb[1:p] + Adb[(p+1):(2*p)]
  
  # Scaled lasso
  Scaledlasso = scalreg(AX, y)
  hsigma = Scaledlasso$hsigma
  
  # Test statistic
  Cov = Trans%*%Inv%*%G%*%t(Inv)%*%Trans
  t1 = sqrt(n)*w1/(hsigma*sqrt(diag(Cov[1:p, 1:p])))
  t2 = sqrt(n)*w2/(hsigma*sqrt(diag(Cov[(p+1):(2*p), (p+1):(2*p)])))
  
  # p-value
  p1 = 2*(1 - pnorm(abs(t1)))
  p2 = 2*(1 - pnorm(abs(t2)))
  p21 = p2
  
  # adjusted p-value
  p2a = p.adjust(p2, method = "BH")
  
  # Knockoff-assisted BH
  index2 = seq(along = p2a)[p2a <= q]
  M1_selected = colnames(X[ ,index2])
  
  # Knockoff-assisted BBH
  p21[p1 > sqrt(q)] = 1
  p21a = p.adjust(p21, method = "BH")
  index21 = seq(along = p21a)[p21a <= sqrt(q)]
  M2_selected = colnames(X[ ,index21])
  
  # Run GM method
  GM = gm(y0, X0, q = q)
  GM_index = GM$gm_selected
  GM_selected = colnames(X0[ ,GM_index])
  
  # Run the knockoff filter.
  knock.gen = function(x) create.fixed(x, method='equi')
  Fixresult = knockoff.filter(X0, y0, fdr=q, knockoffs=knock.gen, statistic=stat.glmnet_lambdasmax)
  knockoff_selected = names(Fixresult$selected)
  
  list(Method1 = M1_selected, Method2 = M2_selected, GM = GM_selected, Knockoff = knockoff_selected)
}



# Run for all 7 responses
# Test fdr = 0.1 
results01 = lapply(Y, function(y) Comparison(X, y, 0.1))

# Calculate TP FP and FDP
get_position <- function(x)
  sapply(regmatches(x, regexpr('[[:digit:]]+', x)), as.numeric)

comparisons01 <- lapply(results01, function(drug) {
  lapply(drug, function(selected) {
    positions = unique(get_position(selected)) # remove possible duplicates
    discoveries = length(positions)
    false_discoveries = length(setdiff(positions, tsm_df$Position))
    list(true_discoveries = discoveries - false_discoveries,
         false_discoveries = false_discoveries,
         fdp = false_discoveries / max(1, discoveries))
  })
})


# Test fdr = 0.2
results02 = lapply(Y, function(y) Comparison(X, y, 0.2))
comparisons02 <- lapply(results02, function(drug) {
  lapply(drug, function(selected) {
    positions = unique(get_position(selected)) # remove possible duplicates
    discoveries = length(positions)
    false_discoveries = length(setdiff(positions, tsm_df$Position))
    list(true_discoveries = discoveries - false_discoveries,
         false_discoveries = false_discoveries,
         fdp = false_discoveries / max(1, discoveries))
  })
})


# Test fdr = 0.05
results005 = lapply(Y, function(y) Comparison(X, y, 0.05))
comparisons005 <- lapply(results005, function(drug) {
  lapply(drug, function(selected) {
    positions = unique(get_position(selected)) # remove possible duplicates
    discoveries = length(positions)
    false_discoveries = length(setdiff(positions, tsm_df$Position))
    list(true_discoveries = discoveries - false_discoveries,
         false_discoveries = false_discoveries,
         fdp = false_discoveries / max(1, discoveries))
  })
})


###################################################################

cnames = c("Algorithm1", "Algorithm2", "GM", "Knockoff", "B-BH")
rnames = c("true_discoveries", "false_discoveries")

#For BBH, we directly extract the data according to Sarkar and Tang (2022).


######### APV 0.05
apv005 = matrix(c(comparisons005$APV$Method1$true_discoveries, comparisons005$APV$Method1$false_discoveries,
                  comparisons005$APV$Method2$true_discoveries, comparisons005$APV$Method2$false_discoveries,
                  comparisons005$APV$GM$true_discoveries, comparisons005$APV$GM$false_discoveries,
                  comparisons005$APV$Knockoff$true_discoveries, comparisons005$APV$Knockoff$false_discoveries,
                  6, 0), nrow = 2, ncol = 5, byrow = FALSE, dimnames = list(rnames, cnames))

pdf("./realdata_plot/apv005.pdf", width = 10, height = 8)
barplot(apv005, main = "Resistance to APV", 
        col = c('navy', 'orange'), ylim = c(0, 40), cex.axis = 2, cex.main = 2, cex.names = 1.5)
dev.off()

######### ATV 0.05
atv005 = matrix(c(comparisons005$ATV$Method1$true_discoveries, comparisons005$ATV$Method1$false_discoveries,
                  comparisons005$ATV$Method2$true_discoveries, comparisons005$ATV$Method2$false_discoveries,
                  comparisons005$ATV$GM$true_discoveries, comparisons005$ATV$GM$false_discoveries,
                  comparisons005$ATV$Knockoff$true_discoveries, comparisons005$ATV$Knockoff$false_discoveries,
                  0, 0), nrow = 2, ncol = 5, byrow = FALSE, dimnames = list(rnames, cnames))

pdf("./realdata_plot/atv005.pdf", width = 10, height = 8)
barplot(atv005, main = "Resistance to ATV", 
        col = c('navy', 'orange'), ylim = c(0, 40), cex.axis = 2, cex.main = 2, cex.names = 1.5)
dev.off()

######### IDV 0.05
idv005 = matrix(c(comparisons005$IDV$Method1$true_discoveries, comparisons005$IDV$Method1$false_discoveries,
                  comparisons005$IDV$Method2$true_discoveries, comparisons005$IDV$Method2$false_discoveries,
                  comparisons005$IDV$GM$true_discoveries, comparisons005$IDV$GM$false_discoveries,
                  comparisons005$IDV$Knockoff$true_discoveries, comparisons005$IDV$Knockoff$false_discoveries,
                  12, 0), nrow = 2, ncol = 5, byrow = FALSE, dimnames = list(rnames, cnames))

pdf("./realdata_plot/idv005.pdf", width = 10, height = 8)
barplot(idv005, main = "Resistance to IDV", 
        col = c('navy', 'orange'), ylim = c(0, 40), cex.axis = 2, cex.main = 2, cex.names = 1.5)
dev.off()

######### LPV 0.05
lpv005 = matrix(c(comparisons005$LPV$Method1$true_discoveries, comparisons005$LPV$Method1$false_discoveries,
                  comparisons005$LPV$Method2$true_discoveries, comparisons005$LPV$Method2$false_discoveries,
                  comparisons005$LPV$GM$true_discoveries, comparisons005$LPV$GM$false_discoveries,
                  comparisons005$LPV$Knockoff$true_discoveries, comparisons005$LPV$Knockoff$false_discoveries,
                  5, 0), nrow = 2, ncol = 5, byrow = FALSE, dimnames = list(rnames, cnames))

pdf("./realdata_plot/lpv005.pdf", width = 10, height = 8)
barplot(lpv005, main = "Resistance to LPV", 
        col = c('navy', 'orange'), ylim = c(0, 40), cex.axis = 2, cex.main = 2, cex.names = 1.5)
dev.off()

######### NFV 0.05
nfv005 = matrix(c(comparisons005$NFV$Method1$true_discoveries, comparisons005$NFV$Method1$false_discoveries,
                  comparisons005$NFV$Method2$true_discoveries, comparisons005$NFV$Method2$false_discoveries,
                  comparisons005$NFV$GM$true_discoveries, comparisons005$NFV$GM$false_discoveries,
                  comparisons005$NFV$Knockoff$true_discoveries, comparisons005$NFV$Knockoff$false_discoveries,
                  10, 1), nrow = 2, ncol = 5, byrow = FALSE, dimnames = list(rnames, cnames))

pdf("./realdata_plot/nfv005.pdf", width = 10, height = 8)
barplot(nfv005, main = "Resistance to NFV", 
        col = c('navy', 'orange'), ylim = c(0, 40), cex.axis = 2, cex.main = 2, cex.names = 1.5)
dev.off()

######### RTV 0.05
rtv005 = matrix(c(comparisons005$RTV$Method1$true_discoveries, comparisons005$RTV$Method1$false_discoveries,
                  comparisons005$RTV$Method2$true_discoveries, comparisons005$RTV$Method2$false_discoveries,
                  comparisons005$RTV$GM$true_discoveries, comparisons005$RTV$GM$false_discoveries,
                  comparisons005$RTV$Knockoff$true_discoveries, comparisons005$RTV$Knockoff$false_discoveries,
                  7, 1), nrow = 2, ncol = 5, byrow = FALSE, dimnames = list(rnames, cnames))

pdf("./realdata_plot/rtv005.pdf", width = 10, height = 8)
barplot(rtv005, main = "Resistance to RTV", 
        col = c('navy', 'orange'), ylim = c(0, 40), cex.axis = 2, cex.main = 2, cex.names = 1.5)
dev.off()

######### APV 0.1
apv01 = matrix(c(comparisons01$APV$Method1$true_discoveries, comparisons01$APV$Method1$false_discoveries,
                 comparisons01$APV$Method2$true_discoveries, comparisons01$APV$Method2$false_discoveries,
                 comparisons01$APV$GM$true_discoveries, comparisons01$APV$GM$false_discoveries,
                 comparisons01$APV$Knockoff$true_discoveries, comparisons01$APV$Knockoff$false_discoveries,
                 7, 1), nrow = 2, ncol = 5, byrow = FALSE, dimnames = list(rnames, cnames))

pdf("./realdata_plot/apv01.pdf", width = 10, height = 8)
barplot(apv01, main = "Resistance to APV", 
        col = c('navy', 'orange'), ylim = c(0, 40), cex.axis = 2, cex.main = 2, cex.names = 1.5)
dev.off()

######### ATV 0.1
atv01 = matrix(c(comparisons01$ATV$Method1$true_discoveries, comparisons01$ATV$Method1$false_discoveries,
                 comparisons01$ATV$Method2$true_discoveries, comparisons01$ATV$Method2$false_discoveries,
                 comparisons01$ATV$GM$true_discoveries, comparisons01$ATV$GM$false_discoveries,
                 comparisons01$ATV$Knockoff$true_discoveries, comparisons01$ATV$Knockoff$false_discoveries,
                 0, 0), nrow = 2, ncol = 5, byrow = FALSE, dimnames = list(rnames, cnames))

pdf("./realdata_plot/atv01.pdf", width = 10, height = 8)
barplot(atv01, main = "Resistance to ATV", 
        col = c('navy', 'orange'), ylim = c(0, 40), cex.axis = 2, cex.main = 2, cex.names = 1.5)
dev.off()

######### IDV 0.1
idv01 = matrix(c(comparisons01$IDV$Method1$true_discoveries, comparisons01$IDV$Method1$false_discoveries,
                 comparisons01$IDV$Method2$true_discoveries, comparisons01$IDV$Method2$false_discoveries,
                 comparisons01$IDV$GM$true_discoveries, comparisons01$IDV$GM$false_discoveries,
                 comparisons01$IDV$Knockoff$true_discoveries, comparisons01$IDV$Knockoff$false_discoveries,
                 14, 2), nrow = 2, ncol = 5, byrow = FALSE, dimnames = list(rnames, cnames))

pdf("./realdata_plot/idv01.pdf", width = 10, height = 8)
barplot(idv01, main = "Resistance to IDV", 
        col = c('navy', 'orange'), ylim = c(0, 40), cex.axis = 2, cex.main = 2, cex.names = 1.5)
dev.off()

######### LPV 0.1
lpv01 = matrix(c(comparisons01$LPV$Method1$true_discoveries, comparisons01$LPV$Method1$false_discoveries,
                 comparisons01$LPV$Method2$true_discoveries, comparisons01$LPV$Method2$false_discoveries,
                 comparisons01$LPV$GM$true_discoveries, comparisons01$LPV$GM$false_discoveries,
                 comparisons01$LPV$Knockoff$true_discoveries, comparisons01$LPV$Knockoff$false_discoveries,
                 5, 0), nrow = 2, ncol = 5, byrow = FALSE, dimnames = list(rnames, cnames))

pdf("./realdata_plot/lpv01.pdf", width = 10, height = 8)
barplot(lpv01, main = "Resistance to LPV", 
        col = c('navy', 'orange'), ylim = c(0, 40), cex.axis = 2, cex.main = 2, cex.names = 1.5)
dev.off()

######### NFV 0.1
nfv01 = matrix(c(comparisons01$NFV$Method1$true_discoveries, comparisons01$NFV$Method1$false_discoveries,
                 comparisons01$NFV$Method2$true_discoveries, comparisons01$NFV$Method2$false_discoveries,
                 comparisons01$NFV$GM$true_discoveries, comparisons01$NFV$GM$false_discoveries,
                 comparisons01$NFV$Knockoff$true_discoveries, comparisons01$NFV$Knockoff$false_discoveries,
                 12, 3), nrow = 2, ncol = 5, byrow = FALSE, dimnames = list(rnames, cnames))

pdf("./realdata_plot/nfv01.pdf", width = 10, height = 8)
barplot(nfv01, main = "Resistance to NFV", 
        col = c('navy', 'orange'), ylim = c(0, 40), cex.axis = 2, cex.main = 2, cex.names = 1.5)
dev.off()

######### RTV 0.1
rtv01 = matrix(c(comparisons01$RTV$Method1$true_discoveries, comparisons01$RTV$Method1$false_discoveries,
                 comparisons01$RTV$Method2$true_discoveries, comparisons01$RTV$Method2$false_discoveries,
                 comparisons01$RTV$GM$true_discoveries, comparisons01$RTV$GM$false_discoveries,
                 comparisons01$RTV$Knockoff$true_discoveries, comparisons01$RTV$Knockoff$false_discoveries,
                 9, 1), nrow = 2, ncol = 5, byrow = FALSE, dimnames = list(rnames, cnames))

pdf("./realdata_plot/rtv01.pdf", width = 10, height = 8)
barplot(rtv01, main = "Resistance to RTV", 
        col = c('navy', 'orange'), ylim = c(0, 40), cex.axis = 2, cex.main = 2, cex.names = 1.5)
dev.off()

######### APV 0.2
apv02 = matrix(c(comparisons02$APV$Method1$true_discoveries, comparisons02$APV$Method1$false_discoveries,
                 comparisons02$APV$Method2$true_discoveries, comparisons02$APV$Method2$false_discoveries,
                 comparisons02$APV$GM$true_discoveries, comparisons02$APV$GM$false_discoveries,
                 comparisons02$APV$Knockoff$true_discoveries, comparisons02$APV$Knockoff$false_discoveries,
                 9, 2), nrow = 2, ncol = 5, byrow = FALSE, dimnames = list(rnames, cnames))

pdf("./realdata_plot/apv02.pdf", width = 10, height = 8)
barplot(apv02, main = "Resistance to APV", 
        col = c('navy', 'orange'), ylim = c(0, 40), cex.axis = 2, cex.main = 2, cex.names = 1.5)
dev.off()

######### ATV 0.2
atv02 = matrix(c(comparisons02$ATV$Method1$true_discoveries, comparisons02$ATV$Method1$false_discoveries,
                 comparisons02$ATV$Method2$true_discoveries, comparisons02$ATV$Method2$false_discoveries,
                 comparisons02$ATV$GM$true_discoveries, comparisons02$ATV$GM$false_discoveries,
                 comparisons02$ATV$Knockoff$true_discoveries, comparisons02$ATV$Knockoff$false_discoveries,
                 0, 0), nrow = 2, ncol = 5, byrow = FALSE, dimnames = list(rnames, cnames))

pdf("./realdata_plot/atv02.pdf", width = 10, height = 8)
barplot(atv02, main = "Resistance to ATV", 
        col = c('navy', 'orange'), ylim = c(0, 40), cex.axis = 2, cex.main = 2, cex.names = 1.5)
dev.off()

######### IDV 0.2
idv02 = matrix(c(comparisons02$IDV$Method1$true_discoveries, comparisons02$IDV$Method1$false_discoveries,
                 comparisons02$IDV$Method2$true_discoveries, comparisons02$IDV$Method2$false_discoveries,
                 comparisons02$IDV$GM$true_discoveries, comparisons02$IDV$GM$false_discoveries,
                 comparisons02$IDV$Knockoff$true_discoveries, comparisons02$IDV$Knockoff$false_discoveries,
                 15, 4), nrow = 2, ncol = 5, byrow = FALSE, dimnames = list(rnames, cnames))

pdf("./realdata_plot/idv02.pdf", width = 10, height = 8)
barplot(idv02, main = "Resistance to IDV", 
        col = c('navy', 'orange'), ylim = c(0, 40), cex.axis = 2, cex.main = 2, cex.names = 1.5)
dev.off()

######### LPV 0.2
lpv02 = matrix(c(comparisons02$LPV$Method1$true_discoveries, comparisons02$LPV$Method1$false_discoveries,
                 comparisons02$LPV$Method2$true_discoveries, comparisons02$LPV$Method2$false_discoveries,
                 comparisons02$LPV$GM$true_discoveries, comparisons02$LPV$GM$false_discoveries,
                 comparisons02$LPV$Knockoff$true_discoveries, comparisons02$LPV$Knockoff$false_discoveries,
                 5, 0), nrow = 2, ncol = 5, byrow = FALSE, dimnames = list(rnames, cnames))

pdf("./realdata_plot/lpv02.pdf", width = 10, height = 8)
barplot(lpv02, main = "Resistance to LPV", 
        col = c('navy', 'orange'), ylim = c(0, 40), cex.axis = 2, cex.main = 2, cex.names = 1.5)
dev.off()

######### NFV 0.2
nfv02 = matrix(c(comparisons02$NFV$Method1$true_discoveries, comparisons02$NFV$Method1$false_discoveries,
                 comparisons02$NFV$Method2$true_discoveries, comparisons02$NFV$Method2$false_discoveries,
                 comparisons02$NFV$GM$true_discoveries, comparisons02$NFV$GM$false_discoveries,
                 comparisons02$NFV$Knockoff$true_discoveries, comparisons02$NFV$Knockoff$false_discoveries,
                 12, 6), nrow = 2, ncol = 5, byrow = FALSE, dimnames = list(rnames, cnames))

pdf("./realdata_plot/nfv02.pdf", width = 10, height = 8)
barplot(nfv02, main = "Resistance to NFV", 
        col = c('navy', 'orange'), ylim = c(0, 40), cex.axis = 2, cex.main = 2, cex.names = 1.5)
dev.off()

######### RTV 0.2
rtv02 = matrix(c(comparisons02$RTV$Method1$true_discoveries, comparisons02$RTV$Method1$false_discoveries,
                 comparisons02$RTV$Method2$true_discoveries, comparisons02$RTV$Method2$false_discoveries,
                 comparisons02$RTV$GM$true_discoveries, comparisons02$RTV$GM$false_discoveries,
                 comparisons02$RTV$Knockoff$true_discoveries, comparisons02$RTV$Knockoff$false_discoveries,
                 11, 2), nrow = 2, ncol = 5, byrow = FALSE, dimnames = list(rnames, cnames))

pdf("./realdata_plot/rtv02.pdf", width = 10, height = 8)
barplot(rtv02, main = "Resistance to RTV", 
        col = c('navy', 'orange'), ylim = c(0, 40), cex.axis = 2, cex.main = 2, cex.names = 1.5)
dev.off()
