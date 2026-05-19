## SES-Related health Inequality Indices
### Complete Monte Carlo Simulation Study

rm(list = ls())

library(data.table)
library(DescTools)

set.seed(12345)
pop = fread(file.choose())

###Function for proposed indices
proposed_ind = function(y, s) {
  n = length(y)
  y_min = min(y)
  y_bar = mean(y)
  x = y - y_min          # UD health
  z = x / y_bar         # RUD health
  t = s / mean(s)       # relative SES
  z_bar = mean(z)
  t_bar = mean(t)
  model = lm(z ~ t)
  b = coef(model)[2]
  IT = sum(z^2) / n
  IV = sum((z - z_bar)^2) / sum(z^2)
  I_SES_V = (b^2 * sum((t - t_bar)^2)) / sum((z - z_bar)^2)
  I_SES_T = (b^2 * sum((t - t_bar)^2)) / sum(z^2)
  result = data.frame(
    IT = IT,
    IV = IV,
    I_SES_V = I_SES_V,
    I_SES_T = I_SES_T
  )
  return(result)
}

###Function for existing indices
existing_ind = function(y, s) {
  n = length(y)
  rank_s = rank(s, ties.method = "average")
  r = rank_s / n
  mu = mean(y)
  CI = (2 / (n * mu)) * sum(y * r) - 1
  ymin = min(y)
  ymax = max(y)
  Erreygers_CI = 4 * mu * CI / (ymax - ymin)
  Gini = DescTools::Gini(y, unbiased = FALSE, na.rm = TRUE)
  y_ratio = y / mean(y)
  theil_part = ifelse(y_ratio == 0,
                      0,
                      y_ratio * log(y_ratio))
  Theil = mean(theil_part)
  Atkinson = DescTools::Atkinson(y, parameter = 0.5, na.rm = TRUE)
  result = data.frame(
    CI = CI,
    Erreygers_CI = Erreygers_CI,
    Gini = Gini,
    Theil = Theil,
    Atkinson = Atkinson
  )
  return(result)
}


###Population true values: proposed indices
res_s1 = proposed_ind(pop$health_s1, pop$ses_s1)
res_s1$Scenario = "S1 Low inequality + weak relation"

res_s2 = proposed_ind(pop$health_s2, pop$ses_s2)
res_s2$Scenario = "S2 Moderate inequality + moderate relation"

res_s3 = proposed_ind(pop$health_s3, pop$ses_s3)
res_s3$Scenario = "S3 High tail inequality + strong relation"

res_s4 = proposed_ind(pop$health_s4, pop$ses_s4)
res_s4$Scenario = "S4 Bimodal society"

res_s5 = proposed_ind(pop$health_s5, pop$ses_s5)
res_s5$Scenario = "S5 Nonlinear logarithmic relation"

res_s6 = proposed_ind(pop$health_s6, pop$ses_s6)
res_s6$Scenario = "S6 Threshold relation"

res_s7 = proposed_ind(pop$health_s7, pop$ses_s7)
res_s7$Scenario = "S7 Bounded proportion health"

res_s8 = proposed_ind(pop$health_s8, pop$ses_s8)
res_s8$Scenario = "S8 Binary health"

proposed_true = rbind(res_s1, res_s2, res_s3, res_s4,
                      res_s5, res_s6, res_s7, res_s8)

rownames(proposed_true) = NULL
proposed_true = proposed_true[, c("Scenario",
                                  "IT",
                                  "IV",
                                  "I_SES_V",
                                  "I_SES_T")]

### Population true values: existing indices
ex_s1 = existing_ind(pop$health_s1, pop$ses_s1)
ex_s1$Scenario = "S1 Low inequality + weak relation"

ex_s2 = existing_ind(pop$health_s2, pop$ses_s2)
ex_s2$Scenario = "S2 Moderate inequality + moderate relation"

ex_s3 = existing_ind(pop$health_s3, pop$ses_s3)
ex_s3$Scenario = "S3 High tail inequality + strong relation"

ex_s4 = existing_ind(pop$health_s4, pop$ses_s4)
ex_s4$Scenario = "S4 Bimodal society"

ex_s5 = existing_ind(pop$health_s5, pop$ses_s5)
ex_s5$Scenario = "S5 Nonlinear logarithmic relation"

ex_s6 = existing_ind(pop$health_s6, pop$ses_s6)
ex_s6$Scenario = "S6 Threshold relation"

ex_s7 = existing_ind(pop$health_s7, pop$ses_s7)
ex_s7$Scenario = "S7 Bounded proportion health"

ex_s8 = existing_ind(pop$health_s8, pop$ses_s8)
ex_s8$Scenario = "S8 Binary health"

existing_true = rbind(ex_s1, ex_s2, ex_s3, ex_s4,
                      ex_s5, ex_s6, ex_s7, ex_s8)

rownames(existing_true) = NULL
existing_true = existing_true[, c("Scenario",
                                  "CI",
                                  "Erreygers_CI",
                                  "Gini",
                                  "Theil",
                                  "Atkinson")]

###Combine true population results
true_results = merge(proposed_true,
                     existing_true,
                     by = "Scenario")

print(true_results)


### Monte Carlo settings

B = 5000
sample_sizes = c(200, 500, 1000)
mc_all = data.frame()


### Monte Carlo simulation
for (n_smp in sample_sizes) {
  cat("Running sample size:", n_smp, "\n")
  for (bb in 1:B) { 
    smp_id = sample(1:nrow(pop),
                    size = n_smp,
                    replace = FALSE)
    smp = pop[smp_id, ]
    p1 = proposed_ind(smp$health_s1, smp$ses_s1)
    e1 = existing_ind(smp$health_s1, smp$ses_s1)
    r1 = cbind(p1, e1)
    r1$Scenario = "S1 Low inequality + weak relation"
    p2 = proposed_ind(smp$health_s2, smp$ses_s2)
    e2 = existing_ind(smp$health_s2, smp$ses_s2)
    r2 = cbind(p2, e2)
    r2$Scenario = "S2 Moderate inequality + moderate relation"
    p3 = proposed_ind(smp$health_s3, smp$ses_s3)
    e3 = existing_ind(smp$health_s3, smp$ses_s3)
    r3 = cbind(p3, e3)
    r3$Scenario = "S3 High tail inequality + strong relation"
    p4 = proposed_ind(smp$health_s4, smp$ses_s4)
    e4 = existing_ind(smp$health_s4, smp$ses_s4)
    r4 = cbind(p4, e4)
    r4$Scenario = "S4 Bimodal society"
    p5 = proposed_ind(smp$health_s5, smp$ses_s5)
    e5 = existing_ind(smp$health_s5, smp$ses_s5)
    r5 = cbind(p5, e5)
    r5$Scenario = "S5 Nonlinear logarithmic relation"
    p6 = proposed_ind(smp$health_s6, smp$ses_s6)
    e6 = existing_ind(smp$health_s6, smp$ses_s6)
    r6 = cbind(p6, e6)
    r6$Scenario = "S6 Threshold relation"
    p7 = proposed_ind(smp$health_s7, smp$ses_s7)
    e7 = existing_ind(smp$health_s7, smp$ses_s7)
    r7 = cbind(p7, e7)
    r7$Scenario = "S7 Bounded proportion health"
    p8 = proposed_ind(smp$health_s8, smp$ses_s8)
    e8 = existing_ind(smp$health_s8, smp$ses_s8)
    r8 = cbind(p8, e8)
    r8$Scenario = "S8 Binary health"
    temp = rbind(r1, r2, r3, r4,
                 r5, r6, r7, r8)
    temp$Sample_Size = n_smp
    temp$Replication = bb
    mc_all = rbind(mc_all, temp)
  }
}

rownames(mc_all) = NULL


###Save raw Monte Carlo results
save_raw = "D:/Lab/SES-Related Health Inequality/final SES paper/R files/MC_raw_results.csv"
fwrite(mc_all, save_raw)
cat("Raw MC results saved.\n")


### Convert MC results to long format
mc_dt = as.data.table(mc_all)
true_dt = as.data.table(true_results)
measure_cols = c("IT", "IV", "I_SES_V", "I_SES_T",
                 "CI", "Erreygers_CI", "Gini",
                 "Theil", "Atkinson")
mc_long = melt(mc_dt,
               id.vars = c("Scenario", "Sample_Size", "Replication"),
               measure.vars = measure_cols,
               variable.name = "Measure",
               value.name = "Estimate")
true_long = melt(true_dt,
                 id.vars = "Scenario",
                 measure.vars = measure_cols,
                 variable.name = "Measure",
                 value.name = "True_Value")


### Merge estimates with true values
mc_long = merge(mc_long,
                true_long,
                by = c("Scenario", "Measure"))

mc_long$Error = mc_long$Estimate - mc_long$True_Value


###Monte Carlo performance summary

mc_summary = mc_long[, .(
  
  True_Value = mean(True_Value),
  
  Mean_Estimate = mean(Estimate, na.rm = TRUE),
  
  Bias = mean(Error, na.rm = TRUE),
  
  Abs_Bias = abs(mean(Error, na.rm = TRUE)),
  
  SD = sd(Estimate, na.rm = TRUE),
  
  RMSE = sqrt(mean(Error^2, na.rm = TRUE))
  
), by = .(Scenario, Sample_Size, Measure)]


### Save summary results

save_summary = "D:/Lab/SES-Related Health Inequality/final SES paper/R files/MC_summary_results.csv"

fwrite(mc_summary, save_summary)

cat("MC summary results saved.\n")


### Print results
print(true_results)
print(mc_summary)


### Final Tables for Paper

### Table 1: Decomposition percentage table

table1_decomposition = copy(proposed_true)

table1_decomposition$Percent_IV = 100 * table1_decomposition$IV
table1_decomposition$Percent_I_SES_V = 100 * table1_decomposition$I_SES_V
table1_decomposition$Percent_I_SES_T = 100 * table1_decomposition$I_SES_T

table1_decomposition$Check_identity = table1_decomposition$IV *
  table1_decomposition$I_SES_V

table1_decomposition$Difference = table1_decomposition$I_SES_T -
  table1_decomposition$Check_identity

print(table1_decomposition)

fwrite(table1_decomposition,
       "D:/Lab/SES-Related Health Inequality/final SES paper/R files/Table1_Decomposition_Percentage.csv")



### Table 2: Comparison with existing measures


table2_comparison = merge(proposed_true[, c("Scenario", "I_SES_T")],
                          existing_true,
                          by = "Scenario")

print(table2_comparison)

fwrite(table2_comparison,
       "D:/Lab/SES-Related Health Inequality/final SES paper/R files/Table2_Comparison_Existing_Proposed.csv")



### Table 3: Monte Carlo stability table
### Mean estimate and SD by sample size


mc_stability = mc_summary[, c("Scenario",
                              "Sample_Size",
                              "Measure",
                              "True_Value",
                              "Mean_Estimate",
                              "SD")]

table3_stability = dcast(as.data.table(mc_stability),
                         Scenario + Measure + True_Value ~ Sample_Size,
                         value.var = c("Mean_Estimate", "SD"))

print(table3_stability)

fwrite(table3_stability,
       "D:/Lab/SES-Related Health Inequality/final SES paper/R files/Table3_MC_Stability.csv")



### Table 4: Rank-preserving SES transformation sensitivity


rank_sensitivity = data.frame()

scenario_names = c(
  "S1 Low inequality + weak relation",
  "S2 Moderate inequality + moderate relation",
  "S3 High tail inequality + strong relation",
  "S4 Bimodal society",
  "S5 Nonlinear logarithmic relation",
  "S6 Threshold relation",
  "S7 Bounded proportion health",
  "S8 Binary health"
)

health_list = list(pop$health_s1, pop$health_s2, pop$health_s3, pop$health_s4,
                   pop$health_s5, pop$health_s6, pop$health_s7, pop$health_s8)

ses_list = list(pop$ses_s1, pop$ses_s2, pop$ses_s3, pop$ses_s4,
                pop$ses_s5, pop$ses_s6, pop$ses_s7, pop$ses_s8)

for (i in 1:8) {
  
  y = health_list[[i]]
  s = ses_list[[i]]
  
  ses_forms = c("Original SES", "Log SES", "Squared SES")
  
  s1 = s
  s2 = log(s)
  s3 = s^2
  
  ses_transformed = list(s1, s2, s3)
  
  for (j in 1:3) {
    
    ss = ses_transformed[[j]]
    
    p_temp = proposed_ind(y, ss)
    e_temp = existing_ind(y, ss)
    
    temp = data.frame(
      Scenario = scenario_names[i],
      SES_Form = ses_forms[j],
      CI = e_temp$CI,
      I_SES_V = p_temp$I_SES_V,
      I_SES_T = p_temp$I_SES_T
    )
    
    rank_sensitivity = rbind(rank_sensitivity, temp)
  }
}

print(rank_sensitivity)

fwrite(rank_sensitivity,
       "D:/Lab/SES-Related Health Inequality/final SES paper/R files/Table4_Rank_Sensitivity.csv")









