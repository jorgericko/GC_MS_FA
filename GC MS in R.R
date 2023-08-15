

#Read GC files
setwd("~/Projects/BETO/ISR/FA_Analysis_R")
gcms_data=read.csv("~/Projects/BETO/ISR/GC-MS/gcms_data_3.csv",header = FALSE)
gcms.df <- data.frame(gcms_data)
colnames(gcms.df) <- c("sample_id", "peak", "rt","area","date","isr","rep","day_cat","batch")



##### CALIBRATION CURVES
#######



## Define how many calibration curves you want to make


###
## Define the injection times for standards in these calibration curves / Batches


#### Get conceentrations for standards, now this time trying retention timpes all standards combined.
### Look at your chromatograms and see if you can do that!

STD <-subset(gcms.df, grepl("STD_[0-9]*|BLANK*|MIX*", sample_id)) # standards second run

#Select reference samples use chromatograms and retention times to choose.

#lowestreference = VFA_STD_10_B2_1
# highest reference VFA_STD_10_B2_2

### Assign Fatty acid IDs based on retention times, for now (in the future will use machine learning incorporating scan values)
#Clean data 
STD_clean <- STD[STD$rt >= 9.1 & STD$rt <= 17.1,] #from STD_10_low_rt and STD_10_high_rt data... need to modify it in the future.

#Here we have to manually view the STD_10 chromatograms for first and second injection within a batch 
#and clean data based on that. We will first remove peaks from standards at early or late retention times


STD_low <- subset(STD_clean, grepl("VFA_STD_10_B2_1", sample_id))
STD_high <- subset(STD_clean, grepl("VFA_STD_10_B2_2", sample_id))

FA=c("C2","C1","C3","ISOC4","C4","ISOC5","C5","ISOC6","C6","C7")

#create a new data frame for acid names and retention time range
FA_rt <- data.frame(FA,STD_low$rt, STD_high$rt)

# Add value1 and value2 columns to my_df, with +1 added to each value

library(dplyr)

FA_rt_range <- FA_rt %>% 
  mutate(STD_low.rt = FA_rt$STD_low.rt -0.1,
         STD_high.rt = FA_rt$STD_high.rt + 0.1)


colnames(FA_rt_range) <- c("acid", "rt.min", "rt.max")

## assign acid names in a new column in the STD clean data set


i = 0
STD_clean$acid <- NA 

for (rt in STD_clean$rt) {
  i = i + 1
  acid_indices <- which(rt >= FA_rt_range$rt.min & rt < FA_rt_range$rt.max)
  #print(acid_indices)
  if (length(acid_indices) == 0) {
    STD_clean$acid[i] <- NA
  } else {
    STD_clean$acid[i] <- FA_rt_range$acid[acid_indices][1]
    #print(FA_rt_range$acid[acid_indices][1])
    
  }
}


## Now print Sample IDs with NAs and inspect those chromatograms if they are noise clean again
print(STD_clean[is.na(STD_clean$acid), ])

### Are they noise? then clean

STD_clean_2 <- STD_clean %>% filter(!is.na(acid))


# create vector of valid values

valid_values <- c("C2", "C3", "C4")

# filter rows where column2 value is not valid
STD_clean_3 <- STD_clean_2[grepl("*MIX*", STD_clean_2$sample_id) & STD_clean_2$acid %in% valid_values, ]

# Create a vector of valid values for MIX 
valid_values <- c("C2", "C3", "C4")

# Filter rows
STD_clean_3 <- subset(STD_clean_2, !grepl("MIX", sample_id) | acid %in% valid_values)

library(stringr)



STD_ref <- STD_clean_3 %>%
  mutate(mM = if_else(str_detect(sample_id, "BLANK"), 0, as.numeric(str_extract(sample_id, "(?<=STD_)\\d+|(?<=MIX_)\\d+"))))


#Filter data set by batch to generate calibration curves

STD_B1 <- subset(STD_ref, !grepl("B[23]", sample_id) & grepl("STD_[0-9]+| BLANK|MIX", sample_id)) #standards first run
STD_B2 <- subset(STD_ref, !grepl("B[3]", sample_id) & grepl("STD_[0-9]*_B2+| BLANK*_B2|MIX_[0-9]*_B2", sample_id)) #standards first run
STD_B3 <- subset(STD_ref, !grepl("B[2]", sample_id) & grepl("STD_[0-9]*_B3+| BLANK*_B2|MIX_[0-9]*_B3", sample_id)) #standards first run



#For the future all standards should have B1,B2,B3  ... Samples probably too as well, that will facilitate instead of messing with injection times


STD_B1$acid <- as.factor(STD_B1$acid)
STD_B2$acid <- as.factor(STD_B2$acid)
STD_B3$acid <- as.factor(STD_B3$acid)


library(patchwork)
library(ggplot2)

lm_list_B1 <- list()
plots_B1 <- list()  # initialize empty list to store plots
results_table_B1 <- data.frame(Acid = character(),
                            RSQ = numeric(),
                            Intercept = numeric(),
                            Slope = numeric(),
                            stringsAsFactors = FALSE)


for (fa in levels(STD_B1$acid)) {
  
  subset_fa <- subset(STD_B1, acid == fa)
  lm_list_B1[[fa]] <- lm(mM ~ area, data = subset_fa)
  
  plot_fa <- ggplot(data = subset_fa, aes(x = area, y = mM)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE, formula = y ~ x) +
    labs(title = fa, x = "Area", y = "mM") 
  
  plots_B1[[fa]] <- plot_fa  # add plot to list
  
  # Get model summary
  model_summary <- summary(lm_list_B1[[fa]])
  
  # Extract R-squared, intercept, and slope from model summary
  rsq <- round(summary(lm_list_B1[[fa]])$r.squared, 3)
  intercept <- round(model_summary$coefficients[1], 3)
  slope <- model_summary$coefficients[2]
  
  # Add results to table
  results_table_B1 <- rbind(results_table_B1, data.frame(Acid = fa,
                                                   RSQ = rsq,
                                                   Intercept = intercept,
                                                   Slope = slope,
                                                   stringsAsFactors = FALSE))
  
  
}

combined_plot_B1 <- wrap_plots(plots_B1)  # combine all plots using patchwork
print(combined_plot_B1)  # display the combined plot
print(results_table_B1) # Print results table





lm_list_B2 <- list()
plots_B2 <- list()  # initialize empty list to store plots
results_table_B2 <- data.frame(Acid = character(),
                               RSQ = numeric(),
                               Intercept = numeric(),
                               Slope = numeric(),
                               stringsAsFactors = FALSE)

for (fa in levels(STD_B2$acid)) {
  
  subset_fa <- subset(STD_B2, acid == fa)
  lm_list_B2[[fa]] <- lm(mM ~ area, data = subset_fa)
  
  plot_fa <- ggplot(data = subset_fa, aes(x = area, y = mM)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE, formula = y ~ x) +
    labs(title = fa, x = "Area", y = "mM") 
  
  plots_B2[[fa]] <- plot_fa  # add plot to list
  
  # Get model summary
  model_summary <- summary(lm_list_B2[[fa]])
  
  # Extract R-squared, intercept, and slope from model summary
  rsq <- round(summary(lm_list_B2[[fa]])$r.squared, 3)
  intercept <- round(model_summary$coefficients[1], 3)
  slope <- model_summary$coefficients[2]
  
  # Add results to table
  results_table_B2 <- rbind(results_table_B2, data.frame(Acid = fa,
                                                         RSQ = rsq,
                                                         Intercept = intercept,
                                                         Slope = slope,
                                                         stringsAsFactors = FALSE))
  
  
}

combined_plot_B2 <- wrap_plots(plots_B2)  # combine all plots using patchwork
print(combined_plot_B2)  # display the combined plot
print(results_table_B2) # Print results table






lm_list_B3 <- list()
plots_B3 <- list()  # initialize empty list to store plots
results_table_B3 <- data.frame(Acid = character(),
                               RSQ = numeric(),
                               Intercept = numeric(),
                               Slope = numeric(),
                               stringsAsFactors = FALSE)

for (fa in levels(STD_B3$acid)) {
  
  subset_fa <- subset(STD_B3, acid == fa)
  lm_list_B3[[fa]] <- lm(mM ~ area, data = subset_fa)
  
  plot_fa <- ggplot(data = subset_fa, aes(x = area, y = mM)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE, formula = y ~ x) +
    labs(title = fa, x = "Area", y = "mM") 
  
  plots_B3[[fa]] <- plot_fa  # add plot to list
  
  # Get model summary
  model_summary <- summary(lm_list_B3[[fa]])
  
  # Extract R-squared, intercept, and slope from model summary
  rsq <- round(summary(lm_list_B3[[fa]])$r.squared, 3)
  intercept <- round(model_summary$coefficients[1], 3)
  slope <- model_summary$coefficients[2]
  
  # Add results to table
  results_table_B3 <- rbind(results_table_B3, data.frame(Acid = fa,
                                                         RSQ = rsq,
                                                         Intercept = intercept,
                                                         Slope = slope,
                                                         stringsAsFactors = FALSE))
  
  
}

combined_plot_B3 <- wrap_plots(plots_B3)  # combine all plots using patchwork
print(combined_plot_B3)  # display the combined plot
print(results_table_B3) # Print results table


#assign acid to gc-ms sample data 
######


#clean df, remove std

gcms.df.nostd <- subset(gcms.df, !grepl("STD|MIX|BLANK", sample_id)) #standards first run
samples <- gcms.df.nostd[gcms.df.nostd$rt >= 9.1 & gcms.df.nostd$rt <= 17.1,] #from STD_10_low_rt and STD_10_high_rt data... need to modify it in the future.

######
######


colnames(FA_rt_range) <- c("acid", "rt.min", "rt.max")

## assign acid names in a new column in the STD clean data set

j = 0
samples$acid <- NA 
for (i in 1:length(samples$rt)) {
  rt <- samples$rt[i]
  acid_indices <- which(rt >= FA_rt_range$rt.min & rt < FA_rt_range$rt.max)
  if (length(acid_indices) == 0) {
    samples$acid[i] <- NA
  } else {
    samples$acid[i] <- FA_rt_range$acid[acid_indices][1]
    # check for repeated rt values and assign the same acid
    repeat_rt_indices <- which(samples$rt == rt & !is.na(samples$acid))
    if (length(repeat_rt_indices) > 1) {
      samples$acid[repeat_rt_indices] <- samples$acid[i]
    }
  }
}

### CLEAN

## Now print Sample IDs with NAs and inspect those chromatograms if they are noise clean again
not_assigned=samples[is.na(samples$acid), ]
samples.clean <- samples %>% filter(!is.na(acid))


#Modify or ISR column
library(stringr)
samples.clean$isr <- paste0("ISR.", gsub("^\\D*(\\d+).*", "\\1", samples.clean$sample_id))
samples.clean$isr <- ifelse(grepl("_D", samples.clean$sample_id) & samples.clean$rep == "D", paste0(samples.clean$isr, ".control"), samples.clean$isr)

#clean C6 data from top

head(samples.clean)
samples.clean_2=samples.clean
samples.clean_2$area <- ifelse(!samples.clean_2$isr %in% c("ISR.3", "ISR.01", "ISR.0", "ISR.05") & samples.clean_2$acid == "C6", 0, samples.clean_2$area)

write.csv(samples.clean_2, file = "data_acids_toclean.csv")

samples.clean_3<-read.csv("data_acids_clean_C6.csv",header = TRUE)



####
## ESTIMATE CONCENTRATIONS
#####

mw <- c(60.052, 46.03, 74.08, 88.11, 88.11, 102.13, 102.13, 116.160, 116.160, 130.1849)
names(mw) <- c("C2", "C1", "C3", "ISOC4", "C4", "ISOC5", "C5", "ISOC6", "C6", "C7")


###
### BATCH 1
####


B1.df <- subset(samples.clean_3, batch =="B1") #standards first run

# create empty mM column
B1.df$mM <- NA
B1.df$gL <- NA

# loop through each row
for (i in seq_along(B1.df$acid)) {
  acid <- B1.df$acid[i]
  lm_model <- lm_list_B1[[acid]]
  new_data <- B1.df[i, ]
  B1.df$mM[i] <- predict(lm_model, newdata = new_data)
  B1.df$gL[i] <- B1.df$mM[i] * mw[acid] / 1000}

# convert mM column to numeric format and round to 2 decimal places
B1.df$mM <- as.numeric(round(B1.df$mM, 2))
B1.df$gL <- as.numeric(round(B1.df$gL, 4))





###
### BATCH 2
####



B2.df <- subset(samples.clean_3, batch =="B2") #standards first run

# create empty mM column
B2.df$mM <- NA
B2.df$gL <- NA

# loop through each row
for (i in seq_along(B2.df$acid)) {
  acid <- B2.df$acid[i]
  lm_model <- lm_list_B2[[acid]]
  new_data <- B2.df[i, ]
  B2.df$mM[i] <- predict(lm_model, newdata = new_data)
  B2.df$gL[i] <- B2.df$mM[i] * mw[acid] / 1000}

# convert mM column to numeric format and round to 2 decimal places
B2.df$mM <- as.numeric(round(B2.df$mM, 2))
B2.df$gL <- as.numeric(round(B2.df$gL, 4))




###
### BATCH 2
####



B3.df <- subset(samples.clean_3, batch =="B3") #standards first run

# create empty mM column
B3.df$mM <- NA
B3.df$gL <- NA

# loop through each row
for (i in seq_along(B3.df$acid)) {
  acid <- B3.df$acid[i]
  lm_model <- lm_list_B3[[acid]]
  new_data <- B3.df[i, ]
  B3.df$mM[i] <- predict(lm_model, newdata = new_data)
  B3.df$gL[i] <- B3.df$mM[i] * mw[acid] / 1000}

# convert mM column to numeric format and round to 2 decimal places
B3.df$mM <- as.numeric(round(B3.df$mM, 2))
B3.df$gL <- as.numeric(round(B3.df$gL, 4))



#Combine data and remove negative values

all.data=rbind(B1.df,B2.df,B3.df)

data_acids_0 <- all.data[all.data$mM > 0, ]

write.csv(data_acids_0, file = "data_acids.csv")


### ************************************
###  LINE PLOTS ----
### ************************************

source(file="libraries.R")
source(file="functions.R")



library(dplyr)
detach("package:plyr", unload = TRUE)

duplicate_rows <- data_acids_0 %>%
  group_by(sample_id,acid) %>%
  summarise(n = n()) %>%
  filter(n > 1)

duplicate_rows


data_acids_clean <- data_acids_0 %>%
  group_by(sample_id,acid) %>%
  filter(gL == max(gL)) %>%
  distinct()


duplicate_rows_check <- data_acids_clean %>%
  group_by(sample_id,acid) %>%
  summarise(n = n()) %>%
  filter(n > 1)

duplicate_rows_check #check there are no duplicate rows

#clean data before plots
data_acids_clean_2 <- subset(data_acids_clean, select = -c(peak, rt,area,date,batch,mM))

# Create the day_label vector
day_label <- setNames(c(0.00, 0.44, 1.52, 3.46, 6.59, 11.63, 19.67), c("D0", "D1", "D2", "D3", "D4", "D5", "D6"))

# Create a new column with the day categories
data_acids_clean_2$day <- day_label[match(as.character(data_acids_clean_2$day_cat), names(day_label))]


# Create a new column with the day categories
data_acids_clean_2$day <- day_label[match(as.character(data_acids_clean_2$day_cat), names(day_label))]


library(dplyr)

data_acids_clean_3 <- data_acids_clean_2 %>% 
  filter(!grepl("^\\s*05_C", sample_id))



#Reformat data frame for plotting
data_acids <- pivot_wider(data_acids_clean_3, names_from = acid, values_from = gL)
data_acids[is.na(data_acids)] = 0
data_acids$total <- rowSums(data_acids[,6:15])
data_acids=as.data.frame(data_acids)

data_acids$isr <- as.factor(data_acids$isr)

#plot


C1=ggline(data_acids, x = "day", y = "C1",numeric.x.axis = TRUE, xlab = "day", ylab = "Formic Acid [g/L]",
          error.plot = "errorbar",add = c("mean_se","jitter"),add.params = list(size = 0.5),color = "isr", shape = "isr",point.size = 1, palette = c("black","red","#CD534CFF","#8F7700FF", "#EFC000FF","blue","#0073C2FF","forestgreen","green","purple","pink" ))+ theme_pubr()+labs_pubr()


C2=ggline(data_acids, x = "day", y = "C2",numeric.x.axis = TRUE, xlab = "day", ylab = "Acetic Acid [g/L]",
          error.plot = "errorbar",add = c("mean_se","jitter"),add.params = list(size = 0.5),color = "isr", shape = "isr",point.size = 1, palette = c("black","red","#CD534CFF","#8F7700FF", "#EFC000FF","blue","#0073C2FF","forestgreen","green","purple","pink" ))+ theme_pubr()+labs_pubr()

C3=ggline(data_acids, x = "day", y = "C3",numeric.x.axis = TRUE, xlab = "day", ylab = "Propionic Acid [g/L]",
          error.plot = "errorbar",add = c("mean_se","jitter"),add.params = list(size = 0.5),color = "isr", shape = "isr",point.size = 1, palette = c("black","red","#CD534CFF","#8F7700FF", "#EFC000FF","blue","#0073C2FF","forestgreen","green","purple","pink" ))+ theme_pubr()+labs_pubr()


C4=ggline(data_acids, x = "day", y = "C4",numeric.x.axis = TRUE, xlab = "day", ylab = "Butyric Acid [g/L]",
          error.plot = "errorbar",add = c("mean_se","jitter"),add.params = list(size = 0.5),color = "isr", shape = "isr",point.size = 1, palette = c("black","red","#CD534CFF","#8F7700FF", "#EFC000FF","blue","#0073C2FF","forestgreen","green","purple","pink" ))+ theme_pubr()+labs_pubr()


ISOC4=ggline(data_acids, x = "day", y = "ISOC4",numeric.x.axis = TRUE, xlab = "day", ylab = "Iso-Butyric Acid [g/L]",
          error.plot = "errorbar",add = c("mean_se","jitter"),add.params = list(size = 0.5),color = "isr", shape = "isr",point.size = 1, palette = c("black","red","#CD534CFF","#8F7700FF", "#EFC000FF","blue","#0073C2FF","forestgreen","green","purple","pink" ))+ theme_pubr()+labs_pubr()



ISOC5=ggline(data_acids, x = "day", y = "ISOC5",numeric.x.axis = TRUE, xlab = "day", ylab = "Iso-Valeric Acid [g/L]",
             error.plot = "errorbar",add = c("mean_se","jitter"),add.params = list(size = 0.5),color = "isr", shape = "isr",point.size = 1, palette = c("black","red","#CD534CFF","#8F7700FF", "#EFC000FF","blue","#0073C2FF","forestgreen","green","purple","pink" ))+ theme_pubr()+labs_pubr()


C5=ggline(data_acids, x = "day", y = "C5",numeric.x.axis = TRUE, xlab = "day", ylab = "Valeric Acid [g/L]",
             error.plot = "errorbar",add = c("mean_se","jitter"),add.params = list(size = 0.5),color = "isr", shape = "isr",point.size = 1, palette = c("black","red","#CD534CFF","#8F7700FF", "#EFC000FF","blue","#0073C2FF","forestgreen","green","purple","pink" ))+ theme_pubr()+labs_pubr()


ISOC6=ggline(data_acids, x = "day", y = "ISOC6",numeric.x.axis = TRUE, xlab = "day", ylab = "Iso-Caproic Acid [g/L]",
          error.plot = "errorbar",add = c("mean_se","jitter"),add.params = list(size = 0.5),color = "isr", shape = "isr",point.size = 1, palette = c("black","red","#CD534CFF","#8F7700FF", "#EFC000FF","blue","#0073C2FF","forestgreen","green","purple","pink" ))+ theme_pubr()+labs_pubr()

C6=ggline(data_acids, x = "day", y = "C6",numeric.x.axis = TRUE, xlab = "day", ylab = "Caproic Acid [g/L]",
             error.plot = "errorbar",add = c("mean_se","jitter"),add.params = list(size = 0.5),color = "isr", shape = "isr",point.size = 1, palette = c("black","red","#CD534CFF","#8F7700FF", "#EFC000FF","blue","#0073C2FF","forestgreen","green","purple","pink" ))+ theme_pubr()+labs_pubr()



C7=ggline(data_acids, x = "day", y = "C7",numeric.x.axis = TRUE, xlab = "day", ylab = "Heptanoic Acid [g/L]",
          error.plot = "errorbar",add = c("mean_se","jitter"),add.params = list(size = 0.5),color = "isr", shape = "isr",point.size = 1, palette = c("black","red","#CD534CFF","#8F7700FF", "#EFC000FF","blue","#0073C2FF","forestgreen","green","purple","pink" ))+ theme_pubr()+labs_pubr()


total=ggline(data_acids, x = "day", y = "total",numeric.x.axis = TRUE, xlab = "day", ylab = "Total [g/L]",
          error.plot = "errorbar",add = c("mean_se","jitter"),add.params = list(size = 0.5),color = "isr", shape = "isr",point.size = 1, palette = c("black","red","#CD534CFF","#8F7700FF", "#EFC000FF","blue","#0073C2FF","forestgreen","green","purple","pink" ))+ theme_pubr()+labs_pubr()


tf=20

C2=ggpar(C2,xlim=c(0,tf))
C3=ggpar(C3,xlim=c(0,tf))
C4=ggpar(C4,xlim=c(0,tf))
C5=ggpar(C5,xlim=c(0,tf))
C6=ggpar(C6,xlim=c(0,tf))
C7=ggpar(C7,xlim=c(0,tf))
ISOC4=ggpar(ISOC4,xlim=c(0,tf))
ISOC5=ggpar(ISOC5,xlim=c(0,tf))
ISOC6=ggpar(ISOC6,xlim=c(0,tf))
total=ggpar(total,xlim=c(0,tf))



ISOC4
ISOC6
ISOC5



plot_list = list(C2,C3,C4,C5,C6,total) 

#plot_list = list(ISOC4,ISOC5,ISOC6) 


ggarrange(plotlist=plot_list)


#### GAS DATA

gas=read.csv("gas.csv",header = TRUE)

gas_2 <- gas %>% 
  filter(!grepl("^\\s*05_C", sample_id))


CH4=ggbarplot(gas_2, x = "isr", y = "CH4_mL",numeric.x.axis = FALSE, xlab = "ISR", ylab = "CH4 [ml]",
             error.plot = "errorbar",add = c("mean_se","jitter"),add.params = list(size = 0.5),color = "isr", palette = c("purple","pink","forestgreen","green","blue","#0073C2FF","#8F7700FF", "#EFC000FF","red","#CD534CFF", "black" ))+ theme_pubr()+labs_pubr()




CO2=ggbarplot(gas_2, x = "isr", y = "CO2_mL",numeric.x.axis = FALSE, xlab = "ISR", ylab = "CO2 [ml]",
              error.plot = "errorbar",add = c("mean_se","jitter"),add.params = list(size = 0.5),color = "isr", palette = c("purple","pink","forestgreen","green","blue","#0073C2FF","#8F7700FF", "#EFC000FF","red","#CD534CFF", "black" ))+ theme_pubr()+labs_pubr()




plot_list = list(C2,C3,C4,C5,C6,total) 

#plot_list = list(ISOC4,ISOC5,ISOC6) 

plot_list = list(C4,total,CH4) 

ggarrange(plotlist=plot_list)



#### then running the libraries.R (source(file="libraries.R")) is redundant but doesn't hurt
source(file="libraries.R")
#### The following loads some functions we use to plot and look at our data (description is here and in the "functions.R" file)
source(file="functions.R")




detach("package:dplyr", unload=TRUE)


data_acids_stacked<-(data_acids %>% pivot_longer(
  cols = C2:total,
  names_to = c("acid"),
  values_to = "conc_gL"))


all_data_acids_sum <- data_summary(data_acids_stacked, varname="conc_gL",groupnames=c("isr","day", "acid"))


write.csv(all_data_acids_sum, "data_acids_summary.csv", row.names=FALSE)






### Last day plots

data_acids_lastday=subset(data_acids, day == '6.59')



data_acids_lastday

total_last=ggbarplot(data_acids_lastday, x = "isr", y = "total",numeric.x.axis = FALSE, xlab = "ISR", ylab = "g/L",
              error.plot = "errorbar",add = c("mean_se","jitter"),add.params = list(size = 0.5),color = "isr", palette = c("purple","pink","forestgreen","green","blue","#0073C2FF","#8F7700FF", "#EFC000FF","red","#CD534CFF", "black" ))+ theme_pubr()+labs_pubr()

total_last
