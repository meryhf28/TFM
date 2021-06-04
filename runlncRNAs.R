# Diagrama de Venn con el solapamiento del programa FEELnc, CPC y CPAT y obtención de la lista de genes lncRNAs finales
# Se carga la librería openxlsx
library(openxlsx)

# Se cargan los archivos correspondientes con la función read.xlsx
cpc <- read.xlsx("/Users/maria hernandez/Desktop/cpc_names.xlsx")
cpat <- read.xlsx("/Users/maria hernandez/Desktop/cpat_names.xlsx")
feelnc <- read.xlsx("/Users/maria hernandez/Desktop/lncRNA_gene.xlsx")

# Se crea una lista con la lista de lncRNAs de los tres programas
lists <- list(feelnc, cpc, cpat) 

# Se carga la librería gplots
library(gplots)
# Se usa la función venn que acepta una lista de conjuntos como argumento, indicando para cada elemento la pertenencia a cada conjunto
venn <- venn(lists)
# Se convierte a data.frame y se obtienen las ocho primeras filas
venn_num <- as.data.frame(head(venn,8))
# Se elimina la primera fila de la variable venn_num
venn_num <- venn_num[-1,]
# Se ordena el data.frame de la manera que el posterior gráfico requiere
venn_num <- rbind.data.frame(venn_num[6,],venn_num[3,],venn_num[2,],
                        venn_num[4,],venn_num[5,],venn_num[1,],venn_num[7,])

# Se carga la librería colorfulVennPlot
library(colorfulVennPlot)
# Se grafica el diagrama de Venn
plotVenn3d(as.vector(venn_num[,1]),labels=c('CPAT','CPC','FEELnc'),Colors=c("chartreuse4","darkorange","salmon","lightskyblue","plum3","lightgoldenrod1","firebrick"))

# Se almacenan las intersecciones en una nueva variable llamada inters
inters <- attr(venn,"intersections")
# Se exporta un fichero xlsx con los nombres de los genes lncRNAs que solapan en los tres programas
write.xlsx(inters$`A:B:C`,file="lncRNAs_finales.xlsx")
