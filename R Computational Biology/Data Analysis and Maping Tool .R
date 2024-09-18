#Patricio Blanco Rafols A01642057
library(seqinr)
library(ggplot2)

#Variantes
covid19 <-read.fasta("covid.fasta")
alpha <-read.fasta("Alpha.fasta")
beta <-read.fasta("beta2.fasta")
gamma <-read.fasta("gamma2.fasta")
delta <- read.fasta("delta2.fasta")
omnicron <-read.fasta("omincron2.fasta")
usa <- read.fasta("usa.fasta")
zimbabwe <- read.fasta("zimbawbue.fasta")
vietnam <- read.fasta("vietnam.fasta")
uk <- read.fasta("uk.fasta")
polonia <- read.fasta("polonia.fasta")
italia <- read.fasta("italia.fasta")
israel <- read.fasta("israel.fasta")
australia <- read.fasta("australia.fasta")
argentina <- read.fasta("argentina.fasta")
alemania <- read.fasta("alemania.fasta")

#longitud Variantes
lcovid <- length(covid19[[1]])
lalpha <- length(alpha[[1]])
lbeta <- length(beta[[1]])
lgamma <- length(gamma[[1]])
ldelta <- length(delta[[1]])
lomicron <- length(omnicron[[1]])
lusa <- length(usa[[1]])
lzimbabwe <- length(zimbabwe[[1]])
lvietnam <- length(vietnam[[1]])
luk <- length(uk[[1]])
lpolonia <- length(polonia[[1]])
litalia <- length(italia[[1]])
lisrael <- length(israel[[1]])
laustralia <- length(australia[[1]])
largentina <- length(argentina[[1]])
lalemania <- length(alemania[[1]])
print(lcovid)
print(lalpha)
print(lbeta)
print(lgamma)
print(ldelta)
print(lomicron)
print(lusa)
print(lzimbabwe)
print(lvietnam)
print(luk)
print(lpolonia)
print(litalia)
print(lisrael)
print(laustralia)
print(largentina)
print(lalemania)

#% de GC
GC(covid19[[1]])*100
GC(alpha[[1]])*100
GC(beta[[1]])*100
GC(gamma[[1]])*100
GC(delta[[1]])*100
GC(omnicron[[1]])*100
GC(usa[[1]])*100
GC(zimbabwe[[1]])*100
GC(vietnam[[1]])*100
GC(uk[[1]])*100
GC(polonia[[1]])*100
GC(italia[[1]])*100
GC(israel[[1]])*100
GC(australia[[1]])*100
GC(argentina[[1]])*100
GC(alemania[[1]])*100

#Secuencias complementarias de las variantes 
compcovid <- comp(covid19[[1]])
compalpha <- comp(alpha[[1]])
compbeta <- comp(beta[[1]])
compgamma <- comp(gamma[[1]])
compdelta <- comp(delta[[1]])
compomicron <- comp(omnicron[[1]])
compusa <- comp(usa[[1]])
compzimbabwe <- comp(zimbabwe[[1]])
compvietnam <- comp(vietnam[[1]])
compuk <- comp(uk[[1]])
comppolonia <- comp(polonia[[1]])
compitalia <- comp(italia[[1]])
compisrael <- comp(israel[[1]])
compaustralia <- comp(australia[[1]])
compargentina <- comp(argentina[[1]])
compalemania <- comp(alemania[[1]])
print(compcovid)
print(compalpha)
print(compbeta)
print(compgamma)
print(compdelta)
print(compomicron)
print(compusa)
print(compzimbabwe)
print(compvietnam)
print(compuk)
print(comppolonia)
print(compitalia)
print(compisrael)
print(compaustralia)
print(compargentina)
print(compalemania)

#Grafica
#funcion que cuenta bases
contarbases <- function(secuencia) {
  secuencia <- toupper(paste(secuencia, collapse = ""))
  secuencia <- unlist(strsplit(secuencia, "", fixed = TRUE))
  bases <- c('A', 'T', 'G', 'C')
  cuentas <- sapply(bases, function(base) sum(secuencia == base))
  data.frame(base = bases, cantidad = cuentas)
}

variantes_nombres <- c("COVID-19", "Alpha", "Beta", "Gamma", "Delta", "Omicron",
                       "USA", "Zimbabwe", "Vietnam", "UK", "Polonia", "Italia",
                       "Israel", "Australia", "Argentina", "Alemania")

variantes_archivos <- c("covid.fasta", "Alpha.fasta", "beta2.fasta", "gamma2.fasta",
                        "delta2.fasta", "omincron2.fasta", "usa.fasta", "zimbawbue.fasta",
                        "vietnam.fasta", "uk.fasta", "polonia.fasta", "italia.fasta",
                        "israel.fasta", "australia.fasta", "argentina.fasta", "alemania.fasta")

datos_variantes <- lapply(variantes_archivos, function(archivo) {
  secuencia <- read.fasta(archivo)[[1]]
  contarbases(sapply(secuencia, as.character))
})

#agregar nombres a variantes 
for (i in seq_along(datos_variantes)) {
  datos_variantes[[i]]$virus <- variantes_nombres[i]
}

#unir en un solo dataframe
dataframe <- do.call(rbind, datos_variantes)

#grafico
graph <- ggplot(dataframe, aes(x = base, y = cantidad, fill = virus)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Comparacion de bases entre variantes",
       x = "Bases",
       y = "Cantidad",
       fill = "Variante del virus")
print(graph)