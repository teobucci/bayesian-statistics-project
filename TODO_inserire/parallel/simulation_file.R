file_name = paste0("./simulation_",idx,".Rdat")
cat('\n inizio rep = ',idx,'\n')

set.seed(123*idx)

# generate data
data = rnorm(n=10)

# run sampler
fit = lm(data ~ 1)

# save results
result = list("coef" = coef(fit), "residuals" = fit$residuals)
save(result, file = file_name)


