# Gamma dist params
shapes_rate_u <- c(30, 100, 1000)
shapes_rate_r <- c(1, 3, 10)
mean_r <- 1e-8
mean_u <- 1e-8

# Geom dist params
avg_rec_spans <- c(1e+3, 1e+4, 1e+5)
avg_mut_spans <- c(1e+3, 1e+4, 1e+5)
exon_lengths <- 1000
num_exons <- c(600, 1500, 3000)
L <- 3e+7
mult_r_u <- c(0.05, 0.1, 0.25)

models <- setDT(crossing(L, num_exons, exon_lengths, mult_r_u,
                         mean_r, mean_u, shapes_rate_u, shapes_rate_r,
                         avg_rec_spans, avg_mut_spans))
models$model <- 1:nrow(models)

fwrite(models, "models.csv")
