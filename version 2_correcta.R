# Load libraries ----------------------------------------------------------
# install.packages('pacman')
install.packages("spgwr")
require(pacman)
pacman::p_load(
  openxlsx, tidyverse, readxl, gtools, glue, KernSmooth, spatstat, crayon, ggspatial, 
  hrbrthemes, RColorBrewer, geodata, fs, sf, classInt,ggpubr,sp,gstat,spgwr,spatstat
)

g <- gc(reset = T)
rm(list = ls())
options(scipen = 999, warn = -1)

# Load data ---------------------------------------------------------------
file <- './tble/Datos_Tesis2.xls'
shts <- excel_sheets(file)
tble <- read_excel(file, sheet = shts)

# Definir las coordenadas en grados decimales (latitud, longitud)
elementos<-data.frame(tble$lon, tble$lat,tble$As)


names(elementos)[1]="longitud"
names(elementos)[2]="latitud"
names(elementos)[3]="As"
datos_arsenico<-elementos


## Crear un objeto sf con las coordenadas de latitud y longitud
utm_crs <- st_crs("+proj=utm +zone=19 +datum=WGS84")
sf_datos_arsenico <- st_as_sf(datos_arsenico, coords = c("longitud", "latitud"), crs = 4326)
datos_sp <- as(sf_datos_arsenico, "Spatial")


# Definir la proyección UTM de destino (por ejemplo, Zona 19)
utm_crs <- CRS("+proj=utm +zone=19 +datum=WGS84")

# Transformar los datos a UTM
datos_sp_utm <- spTransform(datos_sp, utm_crs)
# Convertir datos_arsenico_utm a un data.frame
datos_arsenico_df <- as.data.frame(datos_sp_utm )

#transformar los valores de arsenico con log

datos_arsenico_df$As<-log(datos_arsenico_df$As)


# Load data ---------------------------------------------------------------
proj <- '+proj=utm +zone=19 +datum=WGS84'
#proj <- '+proj=utm +zone=19 +south +ellps=GRS80 +units=m +no_defs +type=crs'

# Polygons 
zone <- st_read('./shpf/raw/LIMITE_URBANO_CENSAL_C17.shp') %>% 
  filter(NOM_CATEG == 'CIUDAD') %>% 
  st_transform(crs = proj)

# Points 
pnt.as <- st_read('./gpkg/pnt_as.gpkg')

# Points as a table 
tbl.as <- pnt.as %>% 
  st_coordinates() %>% 
  as_tibble() %>% 
  setNames(c('x', 'y')) %>% 
  mutate(
    As = pnt.as$value
  )


#tengo una tabla con los valores de la longitu, latitud y arsenico con los datos ar transformado y latitud y longitu en utm corectas

# Functions to use --------------------------------------------------------
F_ballon_h1 <- function(s_i, percentil = 0.20){ #retorna el ancho de banda h_1
  cat(s_i, '\t')
  Data_ballon<-data.frame()
  for(j in 1:397){
    if(s_i!=j){
      valor=abs(tbl.as$As[s_i]-tbl.as$As[j])
      Data_ballon[j,1]<-valor
    }
  }
  Data_ballon <- na.exclude(Data_ballon)
  h_1 = quantile(Data_ballon$V1,percentil)
  return(h_1)
}

Interior<- function(s_i){ # retorna los valores que ingresan al kernel en cada punto si
  cat(s_i, '\t')
  valor_1= tbl.as$As[s_i]
  Data_k<-data.frame()
  for (j in 1:397){
    valor_2= tbl.as$As[j]
    valor=(valor_1 - valor_2)
    Data_k[j,1]<-valor
  }
  return(Data_k)
}


#Funcion kernel Gaussiano

kernel_gg<- function(dis,h_1) {
  exp_val <- -(dis^2) / (2 * h_1^2)
  const <- 1 / (sqrt(2 * pi) * h_1)
  return(const * exp(exp_val))
}

#kernel evaluado en cada punto
Valor_kernel<- function(interior){
  data<-data.frame()
  for (j in 1:397){
    v=interior$V1[j]
    h_1=tbl.as$prcn_20[j]
    valor=kernel_gg(v,h_1)
    data[j,1]<-valor
  }
  return(data)
}








F_ballon_h2 <- function(s_i, percentil){
  
  cat(s_i, '\t')
  gid <- tbl.as[s_i,]
  dfm <- data.frame()
  
  for(j in 1:nrow(tbl.as)){
    cat('Row: ', j, '\n')
    if(s_i != j){
      valor = abs(2.617 - tbl.as$As[j]) 
      dfm[j,1] <- valor
    } else { 
      print('==')    
    }
  }
  
  vls <- dfm$V1
  vls <- na.exclude(vls)
  h_2 <- quantile(vls, percentil)
  return(h_2)
  
}


#k_2 cursiva

K_2<- function(x,s_j,h_2){ #me entrega el valor de K_2 cursiva
  Z=tbl.as$As[s_j]
  up= (x-Z)/h_2
  k<-integrate(function(y) (0.5)*exp(-(y^2)/0.5), lower = Inf, upper = up)
  return(k$value)}


# H1 ----------------------------------------------------------------------
## Cálculo de la diferencia de los valores del arsenico -------------------
prcn.20 <- map(1:nrow(tbl.as), function(i){F_ballon_h1(s_i = i, percentil = 0.20)})
prcn.20 <- as.numeric(unlist(prcn.20))
tbl.as <- mutate(tbl.as, prcn_20 = prcn.20)

# Would lead to take h2 as the m2-percentile of the positive values |x − Z(sj)|, for all j and some m2 ∈ (0, 1).





# x = 2.617  17.4 segun tume luego log(17.4)

# H2 ----------------------------------------------------------------------
vles.h2 <- map(1:nrow(tbl.as), function(i)F_ballon_h2(s_i = i, percentil = 0.10))
vles.h2 <- as.numeric(unlist(vles.h2))
tbl.as <- mutate(tbl.as, h2 = vles.h2)

#los indices 352 al 355 tienen malos los h_2



#Estimar ecuacion 6----------------------------------------------------

F_6<-function(s_i,h_1,h_2,x){
  numerador=0
  interior = Interior(s_i)
  V_kernel=Valor_kernel(interior)
  for (j in 1:397){
    v_kernel=V_kernel$V1[j]
    v_k2=K_2(x,j,h_2)
    valor= v_kernel * v_k2
    numerador=numerador+ valor
  }
  valor=numerador/ (sum(V_kernel))
  return(valor)
}
  


vles.6 <- map(1:nrow(tbl.as), function(i)F_6(s_i = i, h_1=tbl.as$prcn_20[i],h_2=tbl.as$h2[i], x=2.617))
vles.6 <- as.numeric(unlist(vles.6))
tbl.as <- mutate(tbl.as, F_6 = vles.6) # añade la estimacion a la tabla de F_6

# Calcular el kernel por cada ancho de banda que tenemos, diferente para cada punto --------------------
tbl.as
dst <- pnt.as %>% terra::vect() %>% terra::distance(., .) %>% as.data.frame() %>% mutate(gid = paste0('V', 1:nrow(pnt.as))) %>% gather(id, value, -gid) %>% as_tibble() %>% filter(value != 0) 
dst <- dst[1:(nrow(dst)/2),]
dst.min <- min(dst$value) # Distancia minima: 42.2
dst.avg <- mean(dst$value) # Distancia promedio: 3035.5
res <- dst.min * 9 # 380 metros

# To create the grid ------------------------------------------------------
grid <- terra::rast(xmin = ext(zone)[1], xmax = ext(zone)[2], ymin = ext(zone)[3], ymax = ext(zone)[4], 
                    crs = crs(pnt.as),
                    resolution = res
)
values(grid) <- 1
grid <- terra::crop(grid, vect(zone))
grid <- terra::mask(grid, vect(zone))
grid.pnts <- terra::as.data.frame(grid, xy = T) %>% as_tibble()

gg.grid <- ggplot() +  
  geom_sf(data = zone, fill = NA, col = 'grey30') + 
  geom_point(data = grid.pnts, aes(x = x, y = y), col = 'grey20', size = 0.2) +
  geom_sf(data = pnt.as, col = 'darkred', size = 0.7, shape = 2) +
  coord_sf(crs = proj, datum = proj) + 
  labs(x = 'UTM X', y = 'UTM Y') +
  theme_light() + 
  theme(axis.text.y = element_text(angle = 90, hjust = 0.5, size = 5), 
        axis.title = element_text(size = 7, face = 'bold'),
        axis.text.x = element_text(size = 5))

gg.grid
ggsave(plot = gg.grid, filename = './png/maps/grid_v1.png', units = 'in', width = 7, height = 10, dpi = 300)

# Crear un objeto SpatialPointsDataFrame para la grilla
 
cord_grid <- cbind(grid.pnts$x, grid.pnts$y) # Crear una matriz de coordenadas
proj4string = CRS("+proj=utm +zone=19 +datum=WGS84")
# Crear el SpatialPointsDataFrame
S<- SpatialPointsDataFrame(coords = cord_grid, 
                                           data = grid.pnts,
                                           proj4string = CRS("+proj=utm +zone=19 +datum=WGS84"))





# Crear un objeto SpatialPointsDataFrame para los datos

cord<- cbind(tbl.as$x, tbl.as$y)
# Crear el SpatialPointsDataFrame
As_utm<- SpatialPointsDataFrame(coords = cord, 
                                  data = tbl.as,
                                  proj4string = CRS("+proj=utm +zone=19 +datum=WGS84"))


corde_observaciones=coordinates(As_utm)




#----podemos calcular la distancia entre un par de coordenadas y todos los puntos de un objeto espacial. 
H<- function(S,cordenada,percentil){
  # Calcular la distancia entre las coordenadas y todos los puntos en spdf
  distancias <- spDistsN1(S, cordenada)
  #calcular el percentil
  h=quantile(distancias, percentil)
  return(h)
  
  }


vles.H <- map(1:nrow(tbl.as), function(i)H(S, cordenada = coordinates(As_utm)[i,],percentil = 0.15))
vles.H <- as.numeric(unlist(vles.H))
tbl.as <- mutate(tbl.as, H = vles.H) # añade la estimacion a la tabla de F_6

#Calcular la distania entre dos puntos espaciales

Distancia<- function(s_j,s_i){
  x1=s_j[1]
  x2=s_i[1]
  y1=s_j[2]
  y2=s_i[2]
  distance <- sqrt((x2 - x1)^2 + (y2 - y1)^2)
  return(distance)
}

vls_s <- map(1:nrow(tbl.as),function(i)Distancia(coordinates(S)[1,], coordinates(As_utm)[i,]))
vls_s <- as.numeric(unlist(vls_s))


#KERNEL EN DOS DIMENSIONES, GAUSSIANO------

# Función para calcular el kernel gaussiano
gaussian_kernel <- function(distance, sigma) {
  kernel <- (1 / (2 * pi * sigma^2)) * exp(-0.5 * (distance^2) / sigma^2)
  return(kernel)
}
#Evaluando el kernel para un cordnada fija de S
kernel_g <-map(1:nrow(tbl.as),function(i)gaussian_kernel(vls_s[i], tbl.as$H[i]))
kernel_g<- as.numeric(unlist(kernel_g))




estimacion<-function(kernel, estimado){ # kernel es un valor numerico de 397 y estimado valor numerico de 397
  valor=0
  for (i in 1:397){
    v=kernel[i]* estimado[i]
    valor=valor+v
  }
  e=valor/sum(kernel)
  return(e)
  
}
data<-data.frame()
for ( j in 1:nrow(S)){
  S_cordenada= coordinates(S)[j,]
  vls_s <- map(1:nrow(tbl.as),function(i)Distancia(S_cordenada, coordinates(As_utm)[i,]))
  vls_s <- as.numeric(unlist(vls_s))
  #Evaluando el kernel para un cordnada fija de S
  kernel_g <-map(1:nrow(tbl.as),function(i)gaussian_kernel(vls_s[i], tbl.as$H[i]))
  kernel_g<- as.numeric(unlist(kernel_g))
  estimado <- estimacion(kernel_g, tbl.as$F_6)
  data[j,1]<-estimado
}




grid.pnts <- mutate(grid.pnts, F_s = data$V1 )
grid.pnts <- mutate(grid.pnts, plot = 1-data$V1 )
df <- subset(grid.pnts, select = -lyr.1)




