## Conversión de gff3 a gtf

# Se carga la librería rtrackplayer
library("rtracklayer")

# Se importa el archivo solea_v4.1.gff3
dataset <- import("solea_v4.1.gff3")

# Se exporta el fichero solea_v4.1.gtf en formato gtf
export(dataset, "solea_v4.1.gtf", "gtf")
