# plotPCA <- function(pca, pc1, pc2, col, covariates, selectedCov, bty="o"){
#     require(plotly)
#     
#     pc1pc2 <- cbind(pca[,as.numeric(pc1)], pca[,as.numeric(pc2)])
# 
#     ax_x <- list(
#         title = paste("PC", as.numeric(pc1)),
#         zeroline = FALSE,
#         showline = TRUE,
#         showticklabels = TRUE,
#         showgrid = TRUE
#     )
# 
#     ax_y <- list(
#         title = paste("PC", as.numeric(pc2)),
#         zeroline = FALSE,
#         showline = TRUE,
#         showticklabels = TRUE,
#         showgrid = TRUE
#     )
# 
# 
#     plot(as.data.frame(pc1pc2),x=~V1, y=~V2, type = "scatter", mode="markers", text="", color = unique(covariates[,match(selectedCov, colnames(covariates))]), size = 2) %>%
#         layout(xaxis = ax_x, yaxis = ax_y)
#     # xMin <- min(pca[,as.numeric(pc1)])
#     # xMax <- max(pca[,as.numeric(pc1)])
#     # xRange <- xMax - xMin
#     # xlim <- c(xMin-0.05*xRange, xMax+0.20*xRange)
#     # xlab <- paste("PC",as.numeric(pc1), " scores", sep="")
#     # ylab <- paste("PC",as.numeric(pc2), " scores", sep="")
#     # plot(pca[,as.numeric(pc1)], pca[,as.numeric(pc2)],
#     #      col = col, pch = 19, cex = 2, xlab = xlab,
#     #      ylab = ylab, xlim = xlim,
#     #      main = "Principal component analysis (PCA)",
#     #     cex.main = 1.5, cex.lab = 1.5, bty = bty)
#     # uColor <- unique(col)
#     # uCov   <- unique(covariates[,match(selectedCov, colnames(covariates))])
#     # 
#     # #legend("bottomright", legend = uCov, pch = 19, col = uColor,
#     #   #     cex = 1.5, title = selectedCov, bty = "n")
#     # grid()
# }
# 
plotPCA <- function(pca, pc1, pc2, col, covariates, selectedCov, bty="o"){
    xMin <- min(pca[,as.numeric(pc1)])
    xMax <- max(pca[,as.numeric(pc1)])
    xRange <- xMax - xMin
    xlim <- c(xMin-0.05*xRange, xMax+0.20*xRange)
    xlab <- paste("PC",as.numeric(pc1), " scores", sep="")
    ylab <- paste("PC",as.numeric(pc2), " scores", sep="")
    
    plot(pca[,as.numeric(pc1)], pca[,as.numeric(pc2)],
         col = col, pch = 19, cex = 2, xlab = xlab,
         ylab = ylab, xlim = xlim,
         main = "Principal component analysis (PCA)",
         cex.main = 1.5, cex.lab = 1.5, bty = bty)
    uColor <- unique(col)
    uCov   <- unique(covariates[,match(selectedCov, colnames(covariates))])
    
    
    grid()
}
