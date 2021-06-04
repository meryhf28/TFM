## Análisis de expresión diferencial

# Se carga la librería edgeR
library(edgeR)

# Se carga la matriz de recuento conjunta
data = read.table("/Users/maria hernandez/Desktop/htseq_count.txt", header=T, row.names=1, com='')

# Se eliminan los genes no diferencialmente expresados
rnaseqMatrix = data[rowSums(cpm(data) > 1) >= 2,]

# Número de genes eliminados mediante filtrado
nrow(data)-nrow(rnaseqMatrix)

# Se crea el factor conditions
conditions = factor(c(rep("F", 2), rep("M", 2)))

# Se crea un objeto DGEList
exp_study = DGEList(counts=rnaseqMatrix, group=conditions)

## Transformación de los datos 

# Se carga la librería ggplot2 y ggpubr
library(ggplot2)
library(ggpubr)

# Se transforman los datos mediante logaritmo en base 2
pseudocounts <- log2(exp_study$counts + 1)

# Se representan los datos filtrados sin transformar mediante ggplot
notransformado <- ggplot(data.frame(exp_study$counts), aes(x=M1M_1.1)) + 
  xlab("counts") + ylab("Frecuencia") +
  geom_histogram(fill="gray40", binwidth=1000) 
  
# Se representan los datos filtrados y transformados mediante ggplot
transformado <- ggplot(data.frame(pseudocounts), aes(x=M1M_1.1)) +
  xlab(expression(log[2](counts + 1))) + ylab("Frecuencia") +
  geom_histogram(colour="white", fill="gray40", binwidth=0.6)
  
# Se representan ambos gráficos en una misma ventana
ggarrange(notransformado, transformado, ncol=2)

## Representaciones gráficas
### Boxplot
# Se carga la librería reshape
library(reshape)

# Se crea un nuevo conjunto de datos para poder graficar el siguiente diagrama de cajas 
df <- melt(data.frame(pseudocounts), variable_name="Samples")
# Se añade la columna Sexp al nuevo data frame creado
df <- data.frame(df, Sexo=substr(df$Samples, 1, nchar(as.character(df$Samples))-6))
# Se representan los pseudocounts para cada una de las muestras mediante un boxplot
ggplot(df, aes(x=Samples, y=value, fill=Sexo)) + 
    geom_boxplot() + xlab("") + ylab(expression(log[2](count + 1)))

### MA-plot entre muestras
# Se crea la función MAplot para representar un MA-plot entre dos muestras con el fin de simplificar el código
MAplot <- function(sample1, sample2, title) {
x <- pseudocounts[, sample1]
y <- pseudocounts[, sample2]
# Valores de M
M <- x - y
# Valores de A
A <- (x + y)/2
# Se crea un data frame con los valores de M y de A
df <- data.frame(A, M)
# Se representa dichos valores mediante el MA-plot
ggplot(df, aes(x = A, y = M)) + geom_point(size = 1, alpha = 1/5) + 
  ylim(-6,5) + geom_hline(yintercept = 0, color = "blue3") + 
  stat_smooth(se = FALSE, method = "gam", color = "red3") + ggtitle(title) +
  theme(plot.title = element_text(hjust = 0.5))
# Se aplica la función creada anteriormente para representar MA-plots entre las muestras de los individuos femeninos y masculinos
fema <- MAplot(1, 2, "F1 vs F2")
male <- MAplot(3, 4, "M1 vs M2")
# Se representan ambos gráficos en una misma ventana
ggarrange(fema,male)

### Heatmap de la matriz de distancias entre muestras
# Se carga la librería pheatmap
library(pheatmap)
# Se calcula la matriz de distancia euclídea
mat.dist = as.matrix(dist(t(pseudocounts)))
mat.dist = mat.dist/max(mat.dist)
# Se realiza el heatmap
pheatmap(mat.dist,display_numbers=F, cluster_rows = T, cluster_cols = T, cutree_rows=1, clustering_method="ward.D2", clustering_distance_rows="euclidean")

### Gráfico de componentes principales
# Se cargan la librería ggrepel
library(ggrepel)

# Ajustes del gráfico
data <- prcomp(t(pseudocounts),scale=F)
dataDf <- data.frame(data$x)
loads <- round(data$sdev^2/sum(data$sdev^2)*100,1)
# Gráfico entre la PC1 y la PC2
p1 <- ggplot(dataDf,aes(x=PC1, y=PC2)) +
  theme_classic() +
  geom_hline(yintercept = 0, color = "gray70") +
  geom_vline(xintercept = 0, color = "gray70") +
  geom_point(aes(color = as.factor(c(rep("F",2), rep("M",2)))), 
             alpha = 0.55, size = 3) +
  coord_cartesian(xlim = c(min(data$x[,1])-5,max(data$x[,1])+5)) +
  scale_fill_discrete(name = "Group")
# Evita la superposición de etiquetas
p1 + geom_text_repel(aes(y = PC2 + 0.25, 
                         label = rownames(dataDf)),segment.size = 0.25, size = 3) + 
  labs(x = c(paste("PC1",loads[1],"%")),y=c(paste("PC2",loads[2],"%"), title="Gráfico de PCA"),
       colour="Sexo") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values=c("indianred1","dodgerblue1"))

## Normalización de los datos
exp_study = calcNormFactors(exp_study)

## Identificación de genes diferencialmente expresados
### Estimación de los parámetros de dispersión
exp_study = estimateCommonDisp(exp_study)
exp_study = estimateTagwiseDisp(exp_study)

### Representación de las estimaciones de dispersión por gen junto con la relación de dispersión media ajustada
plotBCV(exp_study)
exp_study$common.dispersion

## Análisis de expresión diferencial
# Prueba exacta
et = exactTest(exp_study, pair=c("F", "M"))
# Extracción de los resultados
tTags = topTags(et,n=NULL)
# Print number of up/down significant genes at FDR = 0.05  significance level
summary(de <- decideTestsDGE(et, p=.05))
detags <- rownames(exp_study)[as.logical(de)]

result_table <- tTags$table

# Se exporta un excel con los resultados
write.xlsx(result_table, file='resultados_DEG.xlsx', row.names=T)

# Se carga un archivo con los resultados del análisis de expresión diferencial y un archivo con los nombres de la lista de los lncRNAs obtenidos después del solapamiento de FEELnc, CPC y CPAT
DEG <- read.xlsx("/Users/maria hernandez/Desktop/resultados_DEG.xlsx")
lncRNAs <- read.xlsx("/Users/maria hernandez/Desktop/lncRNAs_finales.xlsx")

# Se carga la librería dplyr
library(dplyr)
# Se filtra de manera que se devuelve el archivo de los DEG con las coincidencias que han habido en el archivo de lncRNAs. Se compara la primera columna (X1), por ello se especifica "by("X1")".
coincidencia <- semi_join(DEG, lncRNAs)
# Se cambia el nombres de las filas de la variable coincidencia por la primera columna de la variable coincidencia
rownames(coincidencia) <- coincidencia[,1]
# Se elimina la primera columna de la variable coincidencia
coincidencia[,1] <- NULL
# Se ordena la variable coincidencia en orden descendente
coincidencia <- arrange(coincidencia,desc(abs(logFC)))
# Se cambia el nombres de las filas de la variable DEG por la primera columna de la variable DEG
rownames(DEG) <- DEG[,1]
# Se elimina la primera columna de la variable DEG
DEG[,1] <- NULL
# Se ordena la variable DEG en orden descendente
DEG <- arrange(DEG,desc(abs(logFC)))
  
# Se crea una función que permite calcular el gráfico MA-plot con ggplot
DEG_MA_plot <- function(datos,x,y,n){
datos1 <- cbind.data.frame(datos,x=x,y=y)
ggplot(datos1, aes(x = x, y = y)) +
  geom_point(cex=2, col=ifelse(datos1$FDR<=0.05,"royalblue1","red"), pch=16, alpha= 0.8)   +
  geom_hline(yintercept = 2,lty="dashed",alpha = 0.8) +
  geom_hline(yintercept = -2,lty="dashed", alpha = 0.8) +
  geom_label_repel(size = 5, data= datos1[1:n,], aes(label=rownames(datos1)[1:n]),
                   max.overlaps = 30, box.padding = 0.5, max.iter = 10000000,
                   colour="slateblue4",segment.colour="royalblue4") + 
  theme(axis.title=element_text(size=14), 
        plot.title = element_text(hjust = 0.5, size=16)) +
  labs(x="logCounts",y="logFC", title="MA-plot")
}

# MA-plot de los genes y los genes lncRNAs diferencialmente expresados sin aplicar ningún criterio
DEG_MA_plot(DEG,DEG$logCPM,DEG$logFC, 5)
DEG_MA_plot(coincidencia,coincidencia$logCPM,coincidencia$logFC, 5)

# Se crea una función que permite calcular el gráfico Volcano plot con ggplot
DEG_Volcano_plot <- function(datos,x,y,n){
datos1 <- cbind.data.frame(datos,x=x,y=y)
ggplot(datos1, aes(x = x, y = y)) +
  geom_point(cex=2, col=ifelse(datos1$FDR<=0.05,"royalblue1","red"), 
             pch=16, alpha= 0.8)  +
  geom_vline(xintercept = 2,lty="dashed",alpha = 0.8) +
  geom_vline(xintercept = -2, lty="dashed", alpha = 0.8) +
  geom_label_repel(size = 5, data= datos1[1:n,], aes(label=rownames(datos1)[1:n]),
                   max.overlaps = 20, max.iter = 10000000, box.padding = 0.3,
                   colour="slateblue4",segment.colour="royalblue4") + 
  theme(axis.title=element_text(size=14), 
        plot.title = element_text(hjust = 0.5, size=16)) +
  labs(x="logFC",y="-1*log10(FDR)", title="Volcano plot")
}

# Volcano plot de los genes y los genes lncRNAs diferencialmente expresados sin aplicar ningún criterio
DEG_Volcano_plot(DEG,DEG$logFC,
                 -1*log10(DEG$FDR), 5)
DEG_Volcano_plot(coincidencia,coincidencia$logFC,
                 -1*log10(coincidencia$FDR), 5)

# Se filtra en base a los siguientes criterios: FDR < 0.5 y |logFC|>2 para obtener el listado de lncRNAs diferencialmente expresados
coincidencia <- subset(coincidencia,coincidencia$FDR < 0.5 & abs(coincidencia$logFC)>2)

# Se exporta un excel con los resultados
write.xlsx(coincidencia, file='lncRNAs_diferencialmente_expresados.xlsx', row.names=T)

## Heatmap de los genes seleccionados entre muestras
# Se transforman los datos
pseudocountsnorm <- cpm(exp_study, log=T, prior.count = 1)
# Se seleccionan los valores de los pseudocounts normalizados de los genes seleccionados
norm_counts <- pseudocountsnorm[row.names(coincidencia),]
  
# Se crea una paleta de color
my_palette <- colorRampPalette(c("indianred1","dodgerblue1"))(n=4)
# Representa el mapa de calor con los pseudocounts normalizados de los genes seleccionados
pheatmap(norm_counts, display_numbers=F, cluster_rows=T, cluster_cols=T, cutree_cols=2, clustering_method="ward.D2",clustering_distance_rows="euclidean",show_rownames=F,cex=1, col=my_palette)
