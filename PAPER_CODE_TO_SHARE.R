library(tidyverse)
library(randomForest)

### datasets for this script can be provided by GitHub links or downloaded from respective Zenodo repository

### script to calculate C loss by local combustion factors and NBR severity levels is available upon request

### script to calculate C stock by forest attributes (to train RF model) is available upon request

### script to process raster data to make spatial predictions is available upon request

### GEE scripts to download imagery and these raster data are available upon request

######### train C model on data of 2016 forest inventory of CEZ and then validate it on independent dataset ###

### carbon data:
carbon.data <- read_csv('https://raw.githubusercontent.com/Janzeero-PhD/C_loss_Chornobyl_2022/main/SHARE_carbon.csv') 
# data set already contains calculated C stock per biomass compartment

carbon.all.alt <- carbon.data %>%
  #select(-12, -17) %>%
  mutate(C = rowSums(.[,c(9:11, 13:18)])) %>%
  rename(CONTROL = ID) %>%
  dplyr::select(CONTROL, C)

# extract data
sample.Sentinel.2016.alt <- read_csv('https://raw.githubusercontent.com/Janzeero-PhD/C_loss_Chornobyl_2022/main/SHARE_Sentinel_2016.csv')

sample.Sentinel.2016.alt$C <- carbon.all.alt$C[match(sample.Sentinel.2016.alt$CONTROL,
                                                     carbon.all.alt$CONTROL)]
summary(sample.Sentinel.2016.alt$C)

sample.Sentinel.2016.alt <- sample.Sentinel.2016.alt %>%
  drop_na()

# test prediction RF:
sample.alt.RF <- randomForest(
  x = sample.Sentinel.2016.alt[, c(3:5, 9, 11:15)],
  y = sample.Sentinel.2016.alt$C,
  importance = T, na.rm = T)
sample.alt.RF
varImpPlot(sample.alt.RF)

### validation
# extract data
validated.Sentinel.2016.alt <- read_csv('https://raw.githubusercontent.com/Janzeero-PhD/C_loss_Chornobyl_2022/main/SHARE_Sentinel_2016_validation.csv')

validated.Sentinel.2016.alt$C <- carbon.live.live$C[match(validated.Sentinel.2016.alt$CONTROL,
                                                          carbon.all.alt$CONTROL)]
validated.Sentinel.2016.alt <- validated.Sentinel.2016.alt %>%
  drop_na()

# test RF model:
validated.Sentinel.2016.alt$predicted <- predict(sample.alt.RF, 
                                                 validated.Sentinel.2016.alt[, c(3:5, 9, 11:15)
                                                                             ])

gmfr.alt.RF.val <- gmfr.data(validated.Sentinel.2016.alt$C, validated.Sentinel.2016.alt$predicted)
gmfr.alt.RF.val

mean(validated.Sentinel.2016.alt$C)
22.6/54.5*100 # 37.2% of the mean

theme_set(theme_bw())
gmfr_graph_RF_val_alt <- gmfr.alt.RF.val[[1]] + 
  scale_x_continuous(limits = c(0, 150)) +
  scale_y_continuous(limits = c(0, 150)) +
  annotate("text", label = c(expression(paste("RMSE = 22.6 ", 'Mg C ', ha^-1))), x = 25, y = 150, size = 2) +
  annotate('text', label = c(expression(paste('AC = ', 0.45))), x = 25, y = 140, size = 2, col = 'blue') +
  annotate('text', label = c(expression(paste('AC sys = ', 0.84))), x = 25, y = 130, size = 2, col = 'blue') +
  annotate('text', label = c(expression(paste('AC uns = ', 0.60))), x = 25, y = 120, size = 2, col = 'blue') +
  labs(x = c(expression(paste("Validation ground truth values, ", MgC))),
       y = c(expression(paste("Predicted values, ", MgC))))
gmfr_graph_RF_val_alt
ggsave('test_RF_val_C_alt.jpg', dpi = 600, units = 'in', width = 3.5, height = 3.5)




########################### C delta (%) model based on SAR backscatter predictors #############################

# VH dif (2016-2020)
VH.dif.data <- read_csv('https://raw.githubusercontent.com/Janzeero-PhD/C_loss_Chornobyl_2022/main/SHARE_VH_dif.csv')

# data to compare
layer <- read_csv('https://raw.githubusercontent.com/Janzeero-PhD/C_loss_Chornobyl_2022/main/SHARE_layer.csv')

layer.comparison <- layer %>%
  dplyr::select(ID, level, C_2020, C_2016, C_dif)

layer.comparison$VH_dif <- VH.dif.data$MEDIAN[match(layer.comparison$ID,
                                                    VH.dif.data$VALUE)]

ggplot(layer.comparison, aes(C_dif*-1, VH_dif*-1)) +
  geom_point() +
  coord_fixed(xlim = c(0,30),
              ylim = c(0,30))

ggplot(layer.comparison, aes(level, VH_dif*-1)) +
  geom_boxplot(aes(group = level))

layer.comparison <- na.omit(layer.comparison)

gmfr.C.VH <- gmfr.data(layer.comparison$C_dif*-1, layer.comparison$VH_dif*-1)
gmfr.C.VH

mean(layer.comparison$C_dif*-1) # 9.3%
sd(layer.comparison$C_dif*-1)


theme_set(theme_bw())
gmfr.C.VH.graph <- gmfr.C.VH[[1]] + 
  scale_x_continuous(limits = c(0, 30)) +
  scale_y_continuous(limits = c(0, 30)) +
  annotate("text", label = c(expression(paste("RMSE = 5.5%"))), x = 5, y = 29, size = 2) +
  annotate('text', label = c(expression(paste('AC = ', 0.18))), x = 5, y = 28, size = 2, col = 'blue') +
  annotate('text', label = c(expression(paste('AC sys = ', 0.90))), x = 5, y = 27, size = 2, col = 'blue') +
  annotate('text', label = c(expression(paste('AC uns = ', 0.27))), x = 5, y = 26, size = 2, col = 'blue') +
  labs(x = c(expression(paste("C loss (inventory data), %"))),
       y = c(expression(paste("SAR backscatter loss, %"))))
gmfr.C.VH.graph
ggsave('gmfr_C_VH_graph.jpg', dpi = 600, units = 'in', width = 3.5, height = 3.5)


layer.comparison$level <- as.factor(layer.comparison$level)

VH.prior.data <- read_csv('https://raw.githubusercontent.com/Janzeero-PhD/C_loss_Chornobyl_2022/main/SHARE_VH_prior.csv')
layer.comparison$VH_prior <- VH.prior.data$MEDIAN[match(layer.comparison$ID,
                                                        VH.prior.data$VALUE)]

VV.prior.data <- read_csv('https://raw.githubusercontent.com/Janzeero-PhD/C_loss_Chornobyl_2022/main/SHARE_VV_prior.csv')
layer.comparison$VV_prior <- VV.prior.data$MEDIAN[match(layer.comparison$ID,
                                                        VV.prior.data$VALUE)]

VV.dif.data <- read_csv('https://raw.githubusercontent.com/Janzeero-PhD/C_loss_Chornobyl_2022/main/SHARE_VV_dif.csv')
layer.comparison$VV_dif <- VV.dif.data$MEDIAN[match(layer.comparison$ID,
                                                    VV.dif.data$VALUE)]

### model:

VH.test.RF <- randomForest(x = layer.comparison[,c(6, 8, 10, 11)],
                           y = layer.comparison$C_dif,
                           importance = T)
VH.test.RF

layer.comparison$C_dif_predict_RF <- predict(VH.test.RF, data = layer.comparison)

ggplot(layer.comparison, aes(C_dif*-1, C_dif_predict_RF*-1)) +
  geom_point() +
  coord_fixed(xlim = c(0,30),
              ylim = c(0,30))

gmfr.C.dif.RF <- gmfr.data(layer.comparison$C_dif*-1, layer.comparison$C_dif_predict_RF*-1)
gmfr.C.dif.RF

theme_set(theme_bw())
gmfr_C_dif_RF <- gmfr.C.dif.RF[[1]] + 
  scale_x_continuous(limits = c(0, 30)) +
  scale_y_continuous(limits = c(0, 30)) +
  annotate("text", label = c(expression(paste("RMSE = 3.9 %"))), x = 5, y = 29, size = 2) +
  annotate('text', label = c(expression(paste('AC = ', -0.31))), x = 5, y = 28, size = 2, col = 'blue') +
  annotate('text', label = c(expression(paste('AC sys = ', 0.73))), x = 5, y = 27, size = 2, col = 'blue') +
  annotate('text', label = c(expression(paste('AC uns = ', -0.05))), x = 5, y = 26, size = 2, col = 'blue') +
  labs(x = c(expression(paste("C loss (inventory data), %"))),
       y = c(expression(paste("C loss (based on SAR backscatter), %"))))
gmfr_C_dif_RF
ggsave('gmfr_C_RF_full_graph.jpg', dpi = 600, units = 'in', width = 3.5, height = 3.5)
