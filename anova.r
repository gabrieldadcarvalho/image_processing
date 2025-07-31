img <- matrix(
  c(
    10, 100, 10,
    50, 10, 50,
    10, 100, 10,
    50, 10, 50
  ),
  nrow = 4,
  byrow = TRUE
)
m <- nrow(img)
n <- ncol(img)

media_linhas <- rowMeans(img)
media_colunas <- colMeans(img)
media_total <- mean(img)

sse <- matrix(0, nrow = m, ncol = n)
ssa <- media_linhas - media_total
sst <- media_colunas - media_total

for (i in 1:m) {
  for (j in 1:n) {
    sse[i, j] <- img[i, j] - media_linhas[i] - media_colunas[j] + media_total
  }
}

sse <- sum(sse^2)
SSa <- sum((ssa)^2) * ncol(img) # fator linha
SSt <- sum((sst)^2) * nrow(img) # fator coluna

print(paste("SSe:", sse))
print(paste("SSa:", SSa))
print(paste("SSt:", SSt))

df1 <- ncol(img) - 1 #
df2 <- (nrow(img) * ncol(img)) - nrow(img) - ncol(img) + 1 

F <- (SSt / df1) / (sse / df2)
Lim <- qf(0.975, df1, df2)
Liminf <- qf(0.025, df1, df2)

print(paste("F-value:", F))
print(paste("Critical value (F(0.975)):", Lim))
print(paste("Critical value (F(0.025)):", Liminf))

print("Média Linhas:")
print(media_linhas)
print(paste("Média Total:", media_total))
