#María Fernanda Bravo
#Librerías
library(igraph)
library(igraphdata)
library(Biostrings)
# Cargar
data("yeast")
#Red
plot(yeast)

#1.Encuentre las diez proteínas más conectadas
ordenado=sort(degree(yeast),decreasing = TRUE)
print("Estas son las 10 proteínas con el mayor número de conexiones");print(ordenado[1:10])

#2.El diámetro de la red
print(paste("El diámetro de la red es ", diameter(yeast)))

#3.El coeficiente de clusterización
##### Para cada nodo
for(i in 1:length(V(yeast))){if (V(yeast)[i]) {print(paste("Coeficiente de clusterización del nodo ", i));print(transitivity(yeast, type= "local", vids=i))}}
##### Para toda la red
print(paste("Coeficiente de clusterización para la red es ", transitivity(yeast)))

#4.El porcentaje de conexiones respecto al total
##El número de nodos de yeast
totalnodos=length(V(yeast))
##El número de conexiones de yeast
totaledges=length(E(yeast))
##Máximo de conexiones posibles para la red (apliqué fórmula)
maxedges=(totalnodos*(totalnodos-1))/2
print(paste("Este es el máximo de conexiones ", maxedges))
##Porcentaje de conexiones respecto al total (Es una regla de 3)
print(paste("El procentaje de conexiones respecto al total es ",(totaledges*100)/maxedges))

#5.El promedio de conectividades
print(paste("El degree promedio es ", mean(degree.distribution(yeast))))

#6.Encuentre los caminos entre una pareja de las 10 proteínas más conectadas 
shortest_paths(yeast,V(yeast)["YPR110C"],V(yeast)["YBR251W"])

#7.Seleccione 10 nodos al azar de la red y
#calcule el promedio de las distancias  Hacer esto 10 veces más y compare este número.
##Tam son los nodos que va a eliminar
tam=sample(yeast,10)
##i son las repeticiones 
for(i in 1:10) {tam=(tam+(i*10))-((i-1)*10);yeast2=delete.vertices(yeast, tam);print(paste("Distancia media", i, mean_distance(yeast2)))}

### CORRECIÓN
yeast2= delete.vertices (yeast, sample (1:2617,10))
for( i in 1:10) {
  print(mean_distance(yeast2))
  yeast2= delete.vertices(yeast2, sample(1:length(V(yeast2)),10))
}


#8. Discute dentro del programa en R este último resultado, 
#esta discusión debe incluir propiedades de la red, 
#pero también suposiciones de por qué una red biológica está "construida" de esta manera. 
### La distancia se refiere a el número mínimo de conexiones que se requieren para llegar de un nodo a otro. 
### En este caso, la distancia media de la red original (yeast) es de 5.095629. 
### Cada que se eliminan 10 nodos al azar, la distancia media se mantiene en 5.09, variando
### solo milésimas, lo cual indica que es una red de mundo pequeño. Este tipo de redes son
### Características de redes biológicas, lo cual es una ventaja respecto a que si es eliminado un nodo,
### no afecta drasticamente las demás relaciones que se mantienen. Por lo tanto, la red mantien su comportamiento