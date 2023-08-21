```{r Import packages}
library(tidyverse)
library(tmap)
library(geojsonio)
library(plotly)
library(rgdal)
library(broom)
library(mapview)
library(crosstalk)
library(sf)
library(sp)
library(spdep)
library(car)
library(fs)
library(janitor)
library(tidypredict)
library(spatialreg)
```

```{r import data}
#import community point data
Sh_community2023<-fs::dir_info(here::here("new"))%>%
  #$ means exact match
  dplyr::filter(str_detect(path, 
                           "Tessellation_sh_2023_new.shp$"))%>%
  dplyr::select(path)%>%
  dplyr::pull()%>%
  #read in the file in
  sf::st_read()

selected_hexagons <- Sh_community2023 %>% 
  filter(MEDIAN_UPI != 0)
qtm(Sh_community2023)
qtm(selected_hexagons)
```

```{r check the original dataset}
#let's map our dependent variable to see if the join has worked:
tmap_mode("plot")
custom_breaks <- c(0, 500, 2000, 5000, 20000, 80000)

map_plot <- tm_shape(selected_hexagons) +
  tm_fill(col = "MEDIAN_UPI", palette = "Blues", title = "Median UPI", breaks = custom_breaks) +
  tm_borders(lwd = 0.5) +
  tm_layout(legend.title.size = 1,
            legend.text.size = 0.8,
            legend.bg.color = "white")
map_plot

q <- qplot(x = `MEDIAN_uni`, 
           y = `MEDIAN_UPI`, 
           data=selected_hexagons)
#plot with a regression line - note, I've added some jitter here as the x-scale is rounded
q + stat_smooth(method="lm", se=FALSE, size=1) + 
  geom_jitter()
```

```{r try the basic OLS regression model}
#run the linear regression model and store its outputs in an object called model1
Regressiondata<- selected_hexagons%>%
  dplyr::select(MEDIAN_UPI,MEDIAN_uni)

#now model
model1 <- Regressiondata %>%
  lm(MEDIAN_UPI ~ MEDIAN_uni,
     data=.)

#show the summary of those outputs
summary(model1)


library(tidymodels)
# set the model
lm_mod <- linear_reg()
# fit the model
lm_fit <- 
  lm_mod %>% 
  fit(MEDIAN_UPI ~ MEDIAN_uni,
     data=Regressiondata)
# we cover tidy and glance in a minute...
tidy(lm_fit)
glance(lm_fit)
```

```{r test regression assumption 1}
#let's check the distribution of these variables first

ggplot(selected_hexagons, aes(x=`MEDIAN_UPI`)) + 
  geom_histogram(aes(y = ..density..),
                 binwidth = 5) + 
  geom_density(colour="red", 
               size=1, 
               adjust=1)
ggplot(selected_hexagons, aes(x=`MEAN_UPI`)) + 
  geom_histogram(aes(y = ..density..),
                 binwidth = 5) + 
  geom_density(colour="red", 
               size=1, 
               adjust=1)

library(ggplot2)
ggplot(selected_hexagons, aes(x=MEDIAN_uni)) + 
  geom_histogram()
ggplot(selected_hexagons, aes(x=MEAN_unit_)) + 
  geom_histogram()
```

```{r transform variables}
ggplot(selected_hexagons, aes(x=log(MEDIAN_uni))) + 
  geom_histogram()

symbox(~MEDIAN_uni, 
       selected_hexagons, 
       na.rm=T,
       powers=seq(-3,3,by=.5))
#log is the most suitable

symbox(~MEAN_unit_, 
       selected_hexagons, 
       na.rm=T,
       powers=seq(-3,3,by=.5))
#log is the most suitable

ggplot(selected_hexagons, aes(x=log(MEDIAN_uni))) + 
  geom_histogram()
ggplot(selected_hexagons, aes(x=log(MEAN_unit_))) + 
  geom_histogram()

# check the distribution again after transforming the variable
qplot(x = log(MEDIAN_uni), 
      y = MEDIAN_UPI,
      data=selected_hexagons)

```

```{r regression model2}
Regressiondata2<- selected_hexagons%>%
  dplyr::select(MEDIAN_UPI,
         MEDIAN_uni)
model2 <- lm(MEDIAN_UPI ~ log(MEDIAN_uni), data = Regressiondata2)

#show the summary of those outputs
summary(model2)

```


```{r test regression assumption 2-5}
#save the residuals into your dataframe

model2_data <- model2 %>%
  augment(., Regressiondata2)

#plot residuals
model2_data%>%
dplyr::select(.resid)%>%
  pull()%>%
  qplot()+ 
  geom_histogram() 

library(performance)
check_model(model2, check="all")

#run durbin-watson test for autocorrelation
DW <- durbinWatsonTest(model2)
tidy(DW)
# the DW statistics for our model is 1.56, so some indication of autocorrelation, which should check for spatial autocorrelation
```

```{r spatial autocorrelation test}
# also add them to the shapelayer
selected_hexagons <- selected_hexagons %>%
  mutate(model2resids = residuals(model2))

#now plot the residuals
library(RColorBrewer)
tmap_mode("view")
custom_breaks <- c(-2000,-500,0,500,2000,5000,10000,25000)
map_plot <- tm_shape(selected_hexagons) +
  tm_polygons("model2resids",
              palette = "OrRd", breaks = custom_breaks) +
  tm_borders(lwd = 0.5) +
  tm_layout(legend.title.size = 1,
            legend.text.size = 0.8,
            legend.bg.color = "white"
            )
map_plot
#the result suggests that there could well be some spatial autocorrelation biasing our model

st_write(selected_hexagons, "new/selected_hexagons.shp")

```


```{r test for spatial autocorrelation more systematically}
#calculate the centroids of all hexagons
coordsW <- selected_hexagons%>%
  st_centroid()%>%
  st_geometry()
plot(coordsW)

#Now we need to generate a spatial weights matrix 
#(remember from the lecture a couple of weeks ago). 
#We'll start with a simple binary matrix of queen's case neighbours
sh_nb <- selected_hexagons %>%
  poly2nb(., queen=T)

#or nearest neighbours
knn_wards <-coordsW %>%
  knearneigh(., k=4)

sh_knn <- knn_wards %>%
  knn2nb()

#plot them
plot(sh_nb, st_geometry(coordsW), col="red")
plot(sh_knn, st_geometry(coordsW), col="blue")

#create a spatial weights matrix object from these weights
#sh.queens_weight <- sh_nb %>%
#  nb2listw(., style="W")
sh.knn_4_weight <- sh_knn %>%
  nb2listw(., style="W")

Nearest_neighbour <- selected_hexagons %>%
  st_drop_geometry()%>%
  dplyr::select(model2resids)%>%
  pull()%>%
  moran.test(., sh.knn_4_weight)%>%
  tidy()
Nearest_neighbour
#Moran's I based on Knn spatial weight matrix is 0.2150085, therefore, some weak to moderate spatial autocorrelation in our residuals, presence of some spatial autocorrelation could be leading to biased estimates of our parameters and significance values

```

```{r SPatial Lag Model (Queen case lag) }
#Original Model
model2 <- lm(MEDIAN_UPI ~ log(MEDIAN_uni), data = selected_hexagons)
summary(model2)

# Queen case lag


```


```{r SPatial Lag Model (Knn case lag) }
#Original Model
model2 <- lm(MEDIAN_UPI ~ log(MEDIAN_uni), data = selected_hexagons)
summary(model2)
glance(model2)

#Knn case lag model
#run a spatially-lagged regression model
model2_knn4 <- lagsarlm(MEDIAN_UPI ~ log(MEDIAN_uni), 
               data = selected_hexagons, 
               nb2listw(sh_knn, 
                        style="C"), 
               method = "eigen")
tidy(model2_knn4)
glance(model2_knn4)

#compare the original model and Knn case lag model
library(lmtest)
lrtest(model2_knn4, model2)

#write out the residuals, and found that there is almost no longer exhibiting spatial autocorrelation (Moran's I is about -0.01125208)
selected_hexagons <- selected_hexagons %>%
  mutate(model2_knn4_resids = residuals(model2_knn4))
KNN4Moran <- selected_hexagons %>%
  st_drop_geometry()%>%
  dplyr::select(model2_knn4_resids)%>%
  pull()%>%
  moran.test(., sh.knn_4_weight)%>%
  tidy()
KNN4Moran

## comparing the AICs of tow models (model2 and Knn case lag model, AIC of model2 is 22214.49, and AIC of knn case lag model is 22110.97)

## interpret the result like: Using the 4 nearest neighbours instead of just considering all adjacent zones in the spatial weights matrix, the size and significance of the spatially lagged term changes quite dramatically. In the 4 nearest neighbour model it is both quite large, positive and statistically significant (<0.05), conversely the effects of unauthorised absence and (log (median house price)) are reduced (CASA0005 contents)
```

```{r spatial error model}
sem_model1 <- errorsarlm(MEDIAN_UPI ~ log(MEDIAN_uni), 
               data = selected_hexagons,
               nb2listw(sh_knn, style="C"), 
               method = "eigen")
tidy(sem_model1)
summary(sem_model1)
```

```{r Lagrange multiplier test}
library(spdep)

sh.knn_weight_ROW <- sh_knn %>%
  nb2listw(., style="C")

lm.LMtests(model2, sh.knn_weight_ROW, test = c("LMerr","LMlag","RLMerr","RLMlag","SARMA"))
```

```{r GWR}
library(spgwr)
coordsW2 <- st_coordinates(coordsW)
selected_hexagons2 <- cbind(selected_hexagons,coordsW2)
GWRbandwidth <- gwr.sel(MEDIAN_UPI ~ log(MEDIAN_uni), 
                  data = selected_hexagons2,                         coords=cbind(selected_hexagons2$X, selected_hexagons2$Y),
                  adapt=T)
# adapt=T here means to automatically find the proportion of observations for the weighting using k nearest neighbours (an adaptive bandwidth)

GWRbandwidth

#run the gwr model
gwr.model = gwr(MEDIAN_UPI ~ log(MEDIAN_uni),
                  data = selected_hexagons2, 
                coords=cbind(selected_hexagons2$X, selected_hexagons2$Y), 
                adapt=GWRbandwidth,
                #matrix output
                hatmatrix=TRUE,
                #standard error
                se.fit=TRUE)

#print the results of the model
gwr.model
```

```{r GWR results}
results <- as.data.frame(gwr.model$SDF)
names(results)

#attach coefficients to original SF
selected_hexagons2 <- selected_hexagons %>%
  mutate(coef_medianprice = results$log.MEDIAN_uni.)
custom_breaks <- c(-15000,-5000,-2000,-500,0,500,2000,5000,12000)
tm_shape(selected_hexagons2) +
  tm_polygons(col = "coef_medianprice", 
              palette = "OrRd", 
              alpha = 0.5, breaks = custom_breaks)

```
```{r GWR-run the significance test}
#run the significance test
sigTest = abs(gwr.model$SDF$"log(MEDIAN_uni)")-2 * gwr.model$SDF$"log(MEDIAN_uni)_se"
#store significance results
selected_hexagons2 <- selected_hexagons2 %>%
  mutate(GWRmedianPriceSig = sigTest)
custom_breaks <- c(-15000,-5000,-2000,-500,0,500,2000,5000,15000)
tm_shape(selected_hexagons2) +
  tm_polygons(col = "GWRmedianPriceSig", 
               palette = "OrRd", breaks = custom_breaks)
```

```{r}
st_write(selected_hexagons2, "new/selected_hexagons2.shp")
```

