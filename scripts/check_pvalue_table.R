ppm15_ptm0 <- read.delim("ppm5_ptm0.txt", header=FALSE)
ppm15_ptm1 <- read.delim("ppm5_ptm1.txt", header=FALSE)
ppm15_ptm2 <- read.delim("ppm5_ptm2.txt", header=FALSE)

peak <- unique(ppm15_ptm0$V1)

for (i in 1:length(peak)) {
  ptm0 <- ppm15_ptm0[which(ppm15_ptm0$V1 == peak[i]),]
  ptm1 <- ppm15_ptm1[which(ppm15_ptm1$V1 == peak[i]),]
  ptm2 <- ppm15_ptm2[which(ppm15_ptm2$V1 == peak[i]),]
  cat("peak = ", peak[i] , "\n")
  for (j in 1:length(ptm0$V3)) {
    if (ptm0$V3[j] != 1) {
      if (!(ptm0$V3[j] < ptm1$V3[j] && ptm1$V3[j] < ptm2$V3[j])) {
        cat(ptm0$V3[j], ptm1$V3[j], ptm2$V3[j], "\n")
      }
    } else {
      if (ptm1$V3[j] != 1) {
        cat(ptm1$V3[j], ptm2$V3[j], "\n")
      } else if (ptm2$V3[j] != 1) {
        cat(ptm2$V3[j], "\n")
      }
    }
  }
}
