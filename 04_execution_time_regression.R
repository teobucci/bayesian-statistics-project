p = c(5, 10, 15, 20, 25, 30, 35, 40, 45, 50)

execution_time = c(
    49.63541 / 60,
    59.66805 / 60,
    1.126382,
    1.29735,
    1.573691,
    1.951381,
    2.314018,
    2.611061,
    3.177836,
    3.701382
)

degree = 2
model_poly <- lm(execution_time ~ poly(p, degree = degree))

# choose how many extra point to the right to predict
extra_p = 100

x_grid = seq(range(p)[1], range(p)[2] + extra_p, length.out = 100)
newdata = data.frame('p' = x_grid)
preds = predict(model_poly, newdata = newdata, se = T)
se.bands = cbind(preds$fit + 2 * preds$se.fit, preds$fit - 2 * preds$se.fit)

output_path = file.path("output")
dir.create(output_path, showWarnings = FALSE, recursive = TRUE)

pdf(file = file.path(output_path, "execution_time.pdf"))

plot(
    p,
    execution_time,
    xlim = c(range(p)[1], range(p)[2] + extra_p), 
    ylim = range(preds$fit),
    cex = 1,
    col = "black",
    main = "Execution time against number of nodes p",
    pch = 16
)
lines(x_grid, preds$fit, lwd = 2, col = "blue")
matlines(x_grid,
         se.bands,
         lwd = 1,
         col = "blue",
         lty = 3)
grid()

dev.off()
