p = c(5,10,15,20,25,30,35,40,45,50)

execution_time = c(49.63541/60,59.66805/60,1.126382,1.29735,1.573691,1.951381,2.314018,2.611061,3.177836,3.701382)


degree = 2
model_poly <- lm(execution_time ~ poly(p, degree = degree))

x_grid = seq(range(p)[1], range(p)[2], length.out = 100)
newdata = data.frame('p' = x_grid)
preds = predict(model_poly, newdata = newdata, se = T)
se.bands = cbind(preds$fit + 2 * preds$se.fit, preds$fit - 2 * preds$se.fit)
plot(p, execution_time, xlim = range(p), cex = 1, col = "black", main = "Execution time against number of nodes", pch=16)
lines(x_grid, preds$fit, lwd = 2, col = "blue")
matlines(x_grid, se.bands, lwd = 1, col = "blue", lty = 3)
grid()


