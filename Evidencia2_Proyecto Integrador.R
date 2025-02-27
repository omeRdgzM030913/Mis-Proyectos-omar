#Evidencia 2

#Jocelyn Ileana Balderas Sánchez - A01798528
#Miguel Angel Galicia Sánchez - A01750744
#Omar Rodríguez Montiel - A0175836
#José Eduardo Rosas Ponciano - A01784461
#Hernán Gael Romero García - A01769342

#Link Vídeo: https://youtu.be/nMvIc9noiz4

library(readr)

# Leer archivos y guardar en vectores separados
Bovino_V <- read_lines("Bovino.txt")
Felino_V <- read_lines("Felino.txt")
Ganso_V <- read_lines("ganso.txt")
Murcielago_H_V <- read_lines("Murcielago_Herradura.txt")
Murcielago_R_V <- read_lines("Murcielago_Rat.txt")
Pangolin_V <- read_lines("Pangolin.txt")
Murcielago_P_V <- read_lines("murcielago_p.txt")
Porcino_V <- read_lines("Porcino.txt")
Roedor_V <- read_lines("roedor.txt")
Wuhan_V <- read_lines("wuhan.txt")

# convirtiendo vectores en secuencias
Bovino_seq <- paste(Bovino_V, collapse="")
Felino_seq <- paste(Felino_V, collapse="")
Ganso_seq <- paste(Ganso_V, collapse="")
Murcielago_H_seq <- paste(Murcielago_H_V, collapse="")
Murcielago_R_seq <- paste(Murcielago_R_V, collapse="")
Pangolin_seq <- paste(Pangolin_V, collapse="")
Murcielago_P_seq <- paste(Murcielago_P_V, collapse="")
Porcino_seq <- paste(Porcino_V, collapse="")
Roedor_seq <- paste(Roedor_V, collapse="")
Wuhan_seq <- paste(Wuhan_V, collapse="")

# función para graficar la distribución de bases
plot_ADN <- function(seq, filename) {
  total <- nchar(seq) #cálculo de la longitud de las secuencias
  Na <- nchar(gsub("[^A]", "", seq))
  Nt <- nchar(gsub("[^T]", "", seq))
  Ng <- nchar(gsub("[^G]", "", seq))
  Nc <- nchar(gsub("[^C]", "", seq))
  NgC <- Ng + Nc
  percentages <- c(Na, Nt, Ng, Nc) / total * 100 #porcentaje de las secuencias
  
  #imprimir porcentajes de cada secuencia de todas las variantes
  print(paste("Variante:", filename))
  print(paste("A:", round(percentages[1],2), "%"))
  print(paste("T:", round(percentages[2],2), "%"))
  print(paste("G:", round(percentages[3],2), "%"))
  print(paste("C:", round(percentages[4],2), "%"))
  print(paste("G+C:", round(percentages[4],2) + round(percentages[3],2), "%")) #imprimir la suma de G y C
  
  return(percentages)
}

# matriz de porcentajes para todas las variantes
percentages_mat <- matrix(nrow = 10, ncol = 4)

#ciclo for para ejecutar la función plot_ADN para cada una de las 5 secuencias de ADN
for (i in 1:10) {
  seq <- switch(i,
                Bovino_seq,
                Felino_seq,
                Ganso_seq,
                Murcielago_H_seq,
                Murcielago_R_seq,
                Pangolin_seq,
                Murcielago_P_seq,
                Porcino_seq,
                Roedor_seq,
                Wuhan_seq)
  percentages_mat[i, ] <- plot_ADN(seq, paste(c("Bovino", "Felino", "Ganso", "Murcielago_H", "Murcielago_R", "Pangolin", "Murcielago_P", "Porcino", "Roedor", "Wuhan")[i], ".txt", sep=""))
  # Llamado a la función plot_ADN para cada secuencia
}

#Se requiere expandir la grafica para visualizarla por completo
# Graficar la matriz de porcentajes
barplot(t(percentages_mat), col = c("chartreuse2", "deepskyblue", "red3", "darkorange"),
        xlab = "Base de ADN", ylab = "Porcentaje", 
        main = "Número de bases de ADN que\n componen a todos los\n coronavirus investigados",
        legend.text = c("A", "T", "G", "C"), 
        names.arg = c("Bovino", "Felino", "Ganso", "Murcielago_H", "Murcielago_R", "Pangolin", "Murcielago_P", "Porcino", "Roedor", "Wuhan"),
        beside = TRUE,
        args.legend = list(x = "topleft", y = "topleft"))


#Creamos el vector con las secuencias

secuencias <- c(Bovino_seq,Felino_seq,Ganso_seq,Murcielago_H_seq,Murcielago_R_seq,Pangolin_seq,Porcino_seq,Roedor_seq,Wuhan_seq)
names(secuencias) <- c("Bovino",'Felino','Ganso','Murcielago_H','Murcielago_R','Pangolin','Porcino','Roedor','Wuhan')

library(msa)

SecAA <- AAStringSet(secuencias)

#Creamos el alineamiento multiple

Alin <- msa(SecAA,"ClustalW")

#Matriz con la distancia

library(seqinr)
Al2<-msaConvert(Alin, "seqinr::alignment")

#Distancia entre los alineamientos

disM <-  dist.alignment(Al2, "identity")
clust <-  hclust(disM)

#Importamos la libreria 'ape' y convertimos el clust en un objeto Phylo

library(ape)
arbol = as.phylo(clust)

#Graficamos el arbol filogenetico y le ponemos color

plot(arbol, type="fan",tip.color = c("purple", "purple", "pink", "gray",
                                     
                                     "green", "blue", "red", "red",
                                     
                                     "red"))