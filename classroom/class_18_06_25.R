# ~~~~~~~~~~~~~~~~~~~~~~~~
# IMPLEMENTAÇÂO MANUAL
# ~~~~~~~~~~~~~~~~~~~~~~~~

# ..> Tranformada de Fourier 2D
manual_dft2d <- function(mat) {
  # mat: Imagem
  M <- nrow(mat)
  N <- ncol(mat)
  F <- matrix(complex(real = 0, imaginary = 0), nrow = M, ncol = N)

  for (u in 0:(M - 1)) {
    for (v in 0:(N - 1)) {
      sum_val <- 0 + 0i
      for (x in 0:(M - 1)) {
        for (y in 0:(N - 1)) {
          angle <- -2 * pi * ((u * x) / M + (v * y) / N)
          sum_val <- sum_val + mat[x + 1, y + 1] * exp(1i * angle)
        }
      }
      F[u + 1, v + 1] <- sum_val
    }
  }

  return(F)
}

# ..> Tranformada Inversa de Fourier 2D

manual_idft2d <- function(F) {
  # F: Matrix com transformada de fourier
  M <- nrow(F)
  N <- ncol(F)
  f <- matrix(0, nrow = M, ncol = N)

  for (x in 0:(M - 1)) {
    for (y in 0:(N - 1)) {
      sum_val <- 0 + 0i
      for (u in 0:(M - 1)) {
        for (v in 0:(N - 1)) {
          angle <- 2 * pi * ((u * x) / M + (v * y) / N)
          sum_val <- sum_val + F[u + 1, v + 1] * exp(1i * angle)
        }
      }
      f[x + 1, y + 1] <- Re(sum_val) / (M * N)
    }
  }

  return(f)
}

fftshift <- function(mat) {
  nr <- nrow(mat)
  nc <- ncol(mat)
  mat[c((nr / 2 + 1):nr, 1:(nr / 2)), c((nc / 2 + 1):nc, 1:(nc / 2))]
}

# ...> Aplicando funções em dados iid

# [[Criando uma matriz/imagem N x N sem estrutura]]
set.seed(1)
N <- 30
mat <- matrix(rgamma(N * N, 2, 2), nrow = N)

hist(mat, freq = FALSE)
curve(dgamma(x, 2, 2), lwd = 2, add = T)
points(density(mat), col = "green", t = "l", lwd = 3, add = T)

# [[Empregando DFT]]
F <- manual_dft2d(mat)
F[1:4, 1:4]
mean(F)

summary(as.numeric(Re(F)))
summary(as.numeric(Im(F)))
summary(as.numeric(Mod(F)))

par(mfrow = c(2, 2))
hist(as.numeric(Re(F)), freq = FALSE)
points(density(Re(F)), col = "green", t = "l", lwd = 3, add = T)
hist(as.numeric(Im(F)), freq = FALSE)
points(density(Im(F)), col = "green", t = "l", lwd = 3, add = T)
hist(as.numeric(Mod(F)), freq = FALSE)
points(density(Mod(F)), col = "green", t = "l", lwd = 3, add = T)
hist(as.numeric(Arg(F)), freq = FALSE)
points(density(Arg(F)), col = "green", t = "l", lwd = 3, add = T)


# [[Empregando IDFT]]
reconstructed <- manual_idft2d(F)

# [[Avaliando recostrução por erro absoluto médio]]
# Esperado: Deve ser muito pequeno
max(abs(mat - reconstructed))

# Visualização

par(mfrow = c(2, 2))
image(
  # t(img_gray[nrow(img_gray):1, ])
  mat,
  col = gray((0:255) / 255),
  main = "Original Image",
  axes = FALSE
)
image(
  # t(img_gray[nrow(img_gray):1, ])
  reconstructed,
  col = gray((0:255) / 255),
  main = "Original Image",
  axes = FALSE
)
image(
  # t(img_gray[nrow(img_gray):1, ])
  Mod(F),
  col = gray((0:255) / 255),
  main = "Original Image",
  axes = FALSE
)
image(
  # t(img_gray[nrow(img_gray):1, ])
  Arg(F),
  col = gray((0:255) / 255),
  main = "Original Image",
  axes = FALSE
)


mag <- Mod(F)
mag_log <- log(1 + mag)
mag_shifted <- fftshift(mag)
image(
  mag_shifted,
  col = gray((0:255) / 255),
  main = "Original Image",
  axes = FALSE
)


# ...> Aplicando funções em dados espaciais

# [[Função Spatial AR Normal]]
rspatialARnorm <-
  function(N, M, alpha, phi10, phi01, phi11, sigma) {
    R <- matrix(
      rnorm((N + 1) * (M + 1), 0, sigma),
      N + 1,
      M + 1
    )

    for (i in 2:(M + 1)) {
      for (j in 2:(N + 1)) {
        # i=2; j=3;
        R[j, i] <- rnorm(
          1,
          alpha +
            phi10 * R[j - 1, i] +
            phi01 * R[j, i - 1] +
            phi11 * R[j - 1, i - 1],
          sigma
        )
      }
    }

    R[-1, -1]
  }

# [[Geração e visualização]]
N <- M <- 50
sigma <- 10
alpha <- 1
phi10 <- 0.9
phi01 <- 0.
phi11 <- 0.
Xobs <- rspatialARnorm(N, M, alpha, phi10, phi01, phi11, sigma)
par(mar = c(0, 0, 0, 0), mfrow = c(1, 1))
image(
  Xobs,
  col = gray((0:255) / 255),
  main = "Original Image",
  axes = FALSE
)


# [[Empregando DFT]]
F <- manual_dft2d(Xobs)
F[1:4, 1:4]
mean(F)

summary(as.numeric(Re(F)))
summary(as.numeric(Im(F)))
summary(as.numeric(Mod(F)))

par(mfrow = c(2, 2))
hist(as.numeric(Re(F)), freq = FALSE)
points(density(Re(F)), col = "green", t = "l", lwd = 3, add = T)
hist(as.numeric(Im(F)), freq = FALSE)
points(density(Im(F)), col = "green", t = "l", lwd = 3, add = T)
hist(as.numeric(Mod(F)), freq = FALSE)
points(density(Mod(F)), col = "green", t = "l", lwd = 3, add = T)
hist(as.numeric(Arg(F)), freq = FALSE)
points(density(Arg(F)), col = "green", t = "l", lwd = 3, add = T)

# [[Empregando IDFT]]
reconstructed <- manual_idft2d(F)

# [[Avaliando recostrução por erro absoluto médio]]
# Esperado: Deve ser muito pequeno
max(abs(Xobs - reconstructed))

# Visualização
par(mfrow = c(2, 2))
image(
  # t(img_gray[nrow(img_gray):1, ])
  Xobs,
  col = gray((0:255) / 255),
  main = "Original Image",
  axes = FALSE
)
image(
  # t(img_gray[nrow(img_gray):1, ])
  reconstructed,
  col = gray((0:255) / 255),
  main = "Original Image",
  axes = FALSE
)
image(
  # t(img_gray[nrow(img_gray):1, ])
  Mod(F),
  col = gray((0:255) / 255),
  main = "Original Image",
  axes = FALSE
)
image(
  # t(img_gray[nrow(img_gray):1, ])
  Arg(F),
  col = gray((0:255) / 255),
  main = "Original Image",
  axes = FALSE
)


mag <- Mod(F)
mag_log <- log(1 + mag)
mag_shifted <- fftshift(mag_log)
par(mfrow = c(1, 1), mar = c(0, 0, 0, 0))
image(
  mag_shifted,
  col = gray((0:255) / 255),
  main = "Original Image",
  axes = FALSE
)

# ~~~~~~~~~~~~~~~~~~~~~~~~

# ~~~~~~~~~~~~~~
# Etapas (RECONSTRUÇÃO)
# ~~~~~~~~~~~~~~
# Leitura da imagem	& readJPEG() & Matriz da imagem
# FFT 2D	          & fftw2d()	 & Transformada de Fourier
# Visualização	    &Mod(), log(), image()	& Espectro de magnitude
# Reconstrução	    & fftw2d(..., inverse=TRUE)	& Imagem reconstruída
# ~~~~~~~~~~~~~~

# 1. Install and Load Required Packages

install.packages(c("jpeg", "fftwtools"))
library(jpeg)
library(fftwtools)

# 2. Load a Grayscale Image

img <- readJPEG(
  "/home/gabrieldadcarvalho/github/image_processing/classroom/cat.jpg"
)
# img <- readImage("E:/Work Activities/8. Digital Library/RANDOM_FIELDS/Books/Material/CODE/TEST1.jpg")

# If image is color, convert to grayscale (simple average)
if (length(dim(img)) == 3) {
  img_gray <- apply(img, c(1, 2), mean)
} else {
  img_gray <- img
}

# 3. Apply 2D Fourier Transform

fft_img <- fftw2d(img_gray)

# Get magnitude (for visualization)
mag <- Mod(fft_img)
mag_log <- log(1 + mag)
# log scale for better viewing

# 4. Center Low Frequencies

fftshift <- function(mat) {
  nr <- nrow(mat)
  nc <- ncol(mat)
  mat[c((nr / 2 + 1):nr, 1:(nr / 2)), c((nc / 2 + 1):nc, 1:(nc / 2))]
}

mag_shifted <- fftshift(mag_log)

# 5. Visualize Original and Frequency Domain

par(mfrow = c(1, 2))

# Original grayscale image
image(
  t(img_gray[nrow(img_gray):1, ]),
  col = gray((0:255) / 255),
  main = "Original Image",
  axes = FALSE
)

# Log-scaled, centered magnitude spectrum
image(
  t(mag_shifted[nrow(mag_shifted):1, ]),
  col = gray((0:255) / 255),
  main = "Fourier Magnitude Spectrum",
  axes = FALSE
)

# 6. (Optional) Inverse FFT to Reconstruct Image

reconstructed <- Re(fftw2d(fft_img, inverse = TRUE)) / length(fft_img)

image(
  t(reconstructed[nrow(reconstructed):1, ]),
  col = gray((0:255) / 255),
  main = "Reconstructed Image",
  axes = FALSE
)


# ~~~~~~~~~~~~~~~~~~~~~~
# Aplicação:
# Filtro Passa-Baixa com FFT para Remoção de Ruído
# ~~~~~~~~~~~~~~~~~~~~~~~
# 1. Aplicar a FFT 2D à imagem.
# 2. Zerar (ou atenuar) as altas frequências
# (bordas, ruído).
# 3. Fazer a FFT inversa para recuperar a
# imagem suavizada.
# ~~~~~~~~~~~~~~~~~~~~~~

fftshift <- function(mat) {
  nr <- nrow(mat)
  nc <- ncol(mat)
  r_idx <- c((nr / 2 + 1):nr, 1:(nr / 2))
  c_idx <- c((nc / 2 + 1):nc, 1:(nc / 2))
  return(mat[r_idx, c_idx])
}


# Filtro1
low_pass_filter <-
  function(fft_matrix, radius) {
    dim_x <- nrow(fft_matrix)
    dim_y <- ncol(fft_matrix)
    center_x <- floor(dim_x / 2) + 1
    center_y <- floor(dim_y / 2) + 1
    mask <- matrix(0, nrow = dim_x, ncol = dim_y)

    for (i in 1:dim_x) {
      for (j in 1:dim_y) {
        distance <-
          sqrt((i - center_x)^2 + (j - center_y)^2)
        if (distance < radius) {
          mask[i, j] <- 1
        }
      }
    }

    # Transformar a máscara em complexa, com parte imaginária 0
    mask <- as.complex(mask)

    return(mask)
  }

# Reading image
img <- readJPEG(
  "/home/gabrieldadcarvalho/github/image_processing/classroom/cat.jpg"
)
# If image is color, convert to grayscale (simple average)
if (length(dim(img)) == 3) {
  img_gray <- apply(img, c(1, 2), mean)
} else {
  img_gray <- img
}

# plotting
N <- dim(img_gray1 <- img_gray)
noise <- runif(N[1] * N[2])
dim(noise) <- N
img_gray <- img_gray + noise
par(mfrow = c(1, 2))
image(
  t(img_gray1[nrow(img_gray1):1, ]),
  col = gray((0:255) / 255),
  main = "Imagem Original",
  axes = FALSE
)
image(
  t(img_gray[nrow(img_gray):1, ]),
  col = gray((0:255) / 255),
  main = "Imagem Original",
  axes = FALSE
)


# Reaplica FFT com fftw2d
fft_img <- fftw2d(img_gray)

# Centraliza
fft_centered <- fftshift(fft_img)

# Gera a máscara compatível com a matriz complexa
# Filtro1
mask <-
  low_pass_filter(fft_centered, radius = 20)


# Aplica a máscara
# Filtro1
fft_filtered <- fft_centered * mask

# Desfaz centralização
# Filtro1
fft_filtered_uncentered <- fftshift(fft_filtered)

# Inversa e reconstrução
# Filtro1
img_filtered <-
  Re(fftw2d(fft_filtered_uncentered, inverse = TRUE)) /
    length(fft_filtered_uncentered)
img_filtered <-
  (img_filtered - min(img_filtered)) /
    (max(img_filtered) - min(img_filtered)) # normaliza

# Plota
par(mfrow = c(1, 2))
image(
  t(img_gray[nrow(img_gray):1, ]),
  col = gray((0:255) / 255),
  main = "Imagem Original Cotaminada",
  axes = FALSE
)
image(
  t(img_filtered[nrow(img_filtered):1, ]),
  col = gray((0:255) / 255),
  main = "Imagem Suavizada1",
  axes = FALSE
)
