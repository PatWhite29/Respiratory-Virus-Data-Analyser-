library("seqinr")
library("ggplot2")
library("DECIPHER")
library(ape)

#variantes
covid <- read.fasta("covid.fasta")
alemania <- read.fasta("alemania.fasta")
argentina <- read.fasta("argentina.fasta")
australia <- read.fasta("australia.fasta")
israel <- read.fasta("israel.fasta")
italia <- read.fasta("italia.fasta")
polonia <- read.fasta("polonia.fasta")
uk <- read.fasta("uk.fasta")
usa <- read.fasta("usa.fasta")
vietnam <- read.fasta("vietnam.fasta")
zimbawbue <- read.fasta("zimbawbue.fasta")
uruguay <- read.fasta("uruguay.fasta")
japan <- read.fasta("japan.fasta")
gales <- read.fasta("gales.fasta")
brazil <- read.fasta("brazil.fasta")
luxenburgo <- read.fasta("luxenburgo.fasta")
ghana <- read.fasta("ghana.fasta")
espana <- read.fasta("espana.fasta")
qatar <- read.fasta("qatar.fasta")
russia <- read.fasta("russia.fasta")

#longitud
lcovid <- length(covid[[1]])
lalemania <- length(alemania[[1]])
largentina <- length(argentina[[1]])
laustralia <- length(australia[[1]])
lisrael <- length(israel[[1]])
litalia <- length(italia[[1]])
lpolonia <- length(polonia[[1]])
luk <- length(uk[[1]])
lusa <- length(usa[[1]])
lvietnam <- length(vietnam[[1]])
lzimbawbue <- length(zimbawbue[[1]])
luruguay <- length(uruguay[[1]])
ljapan <- length(japan[[1]])
lgales <- length(gales[[1]])
lbrazil <- length(brazil[[1]])
lluxenburgo <- length(luxenburgo[[1]])
lghana <- length(ghana[[1]])
lespana <- length(espana[[1]])
lqatar <- length(qatar[[1]])
lrussia <- length(russia[[1]])

#Grafica comparando % bases de variantes
contarbases <- function(secuencia) {
  secuencia <- toupper(paste(secuencia, collapse = ""))
  secuencia <- unlist(strsplit(secuencia, "", fixed = TRUE))
  bases <- c('A', 'T', 'G', 'C')
  cuentas <- sapply(bases, function(base) sum(secuencia == base))
  data.frame(base = bases, cantidad = cuentas)
}
variantes_nombres <- c("COVID-19", "Alemania", "Argentina", "Australia", "Israel", "Italia",
                       "Polonia", "UK", "USA", "Vietnam", "Zimbabwe", "Uruguay", "Japón",
                       "Gales", "Brasil", "Luxemburgo", "Ghana", "España", "Qatar", "Rusia")

variantes_archivos <- c("covid.fasta", "alemania.fasta", "argentina.fasta", "australia.fasta",
                        "israel.fasta", "italia.fasta", "polonia.fasta", "uk.fasta",
                        "usa.fasta", "vietnam.fasta", "zimbawbue.fasta", "uruguay.fasta",
                        "japan.fasta", "gales.fasta", "brazil.fasta", "luxenburgo.fasta",
                        "ghana.fasta", "espana.fasta", "qatar.fasta", "russia.fasta")

datos_variantes <- lapply(variantes_archivos, function(archivo) {
  secuencia <- read.fasta(archivo)[[1]]
  contarbases(sapply(secuencia, as.character))
})
for (i in seq_along(datos_variantes)) {
  datos_variantes[[i]]$virus <- variantes_nombres[i]
}
dataframe <- do.call(rbind, datos_variantes)
ggplot(dataframe, aes(x = base, y = cantidad, fill = virus)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Comparacion de bases entre variantes de virus por pas",
       x = "Base",
       y = "Cantidad",
       fill = "País")

#arbol
secuencias <- lapply(list("alemania.fasta", "argentina.fasta", "australia.fasta", "israel.fasta",
                          "italia.fasta", "polonia.fasta", "uk.fasta", "usa.fasta",
                          "vietnam.fasta", "zimbawbue.fasta", "uruguay.fasta", "japan.fasta",
                          "gales.fasta", "brazil.fasta", "luxenburgo.fasta", "ghana.fasta",
                          "espana.fasta", "qatar.fasta", "russia.fasta"), readDNAStringSet)

secuencia_no_alineada <- do.call(c, secuencias)

#alinear las secuencias
secuencia_alineada <- DECIPHER::AlignSeqs(secuencia_no_alineada)
secuencia_alineada_dnabin <- as.DNAbin(secuencia_alineada)
matriz_distancia <- dist.dna(secuencia_alineada_dnabin, model = "JC", pairwise.deletion = TRUE)
arbol_filogenetico <- njs(matriz_distancia)
arbol_filogenetico <- ladderize(arbol_filogenetico)

#graficar
plot(arbol_filogenetico, no.margin=TRUE, edge.width=2)
title("Árbol Filogenético de Variantes de COVID-19 por País") 