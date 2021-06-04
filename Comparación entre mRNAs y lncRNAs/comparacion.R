# Comparación entre mRNAs y lncRNAs

# Se carga la librería openxlxs
library(openxlsx)
# Se cargan los archivos que contienen la longitud y el porcentaje GC de lncRNAs y mRNAs
lncRNA <- read.xlsx("/Users/maria hernandez/Desktop/lncRNA_GC_length.xlsx")
mRNA <- read.xlsx("/Users/maria hernandez/Desktop/mRNA_GC_length.xlsx")
# Se carga el archivo que contiene el número de exones de lncRNAs y mRNAs
exones <- read.xlsx("/Users/maria hernandez/Desktop/exones.xlsx")

# Se convierte en numérico las variables GC de las dos archivos cargados anteriormente
lncRNA$GC <- as.numeric(lncRNA$GC)
mRNA$GC <- as.numeric(mRNA$GC)

# Se unen los dos data frames en uno llamado genes
genes <- merge(lncRNA, mRNA, all=T)

# Se carga la librería ggplot2
library(ggplot2)
# Se representa un diagrama de cajas de la longitud de secuencia de mRNAs y lncRNAs
ggplot(genes, aes(x=Type, y=Length, fill=Type)) + 
    geom_boxplot() + xlab("") + ylab("Longitud de secuencia") +  theme(legend.title=element_blank())

# Se representa un diagrama de cajas del contenido GC de mRNAs y lncRNAs
ggplot(genes, aes(x=Type, y=GC, fill=Type)) + 
    geom_boxplot() + xlab("") + ylab("Contenido GC") +
   theme(legend.title=element_blank())

# Se divide la variable GC de genes en varios intervalos
sort_genes <- subset(genes,genes$Type==sort(genes$Type))
genes_cut_GC <- cut(genes$GC, breaks = seq(20, 80, by = 10))
genes_type <- c(rep("lncRNA",296),rep("mRNA",4577))

# Se unen las dos variables creadas anteriormente en un data.frame
genes_GC <- data.frame(genes_cut_GC, genes_type)
# Se cambia el nombre de las columnas
colnames(genes_GC) <- c("GC","Type")

# Se representa un gráfico de barras del contenido GC de mRNAs y lncRNAs
ggplot(genes_GC, aes(x=GC, fill=Type)) + 
  geom_bar(position="dodge")  +
  theme(legend.title=element_blank()) +
  labs(y="", x="Contenido GC")
  
# Se representa un gráfico de barras del número de exones de mRNAs y lncRNAs
ggplot(exones, aes(x=Exones, fill=Type)) + 
  geom_bar(position="dodge")  +
  theme(legend.title=element_blank()) +
  labs(y="", x="Número de exones")

# Se calcula el mínimo y el máximo de la longitud de los lncRNAs
min(lncRNA$Length)
max(lncRNA$Length)
# Se calcula la media y la mediana de la longitud de los lncRNAs
mean(lncRNA$Length)
median(lncRNA$Length)

# Se calcula el mínimo y el máximo de la longitud de los mRNAs
min(mRNA$Length)
max(mRNA$Length)
# Se calcula la media y la mediana de la longitud de los mRNAs
mean(mRNA$Length)
median(mRNA$Length)

# Se calcula la media y la mediana del contenido GC de los lncRNAs
mean(lncRNA$GC)
median(lncRNA$GC)

# Se calcula la media y la mediana del contenido GC de los mRNAs
mean(mRNA$GC)
median(mRNA$GC)
