library(ggplot2);library(dplyr);library(xlsx)
data <- read.xlsx2('C:/Users/cma86/Desktop/plasma_film.xlsx',
                  sheetName = 'Processed',
                  colClasses = NA)

ggplot(data,aes(Wafer,Weight,color=Donor))+geom_point()
data %>% group_by(Wafer) %>% summarize(n=n(),weight=mean(Weight),SE=sd(Weight)/sqrt(n))

block <-aov(Weight~Donor+Wafer,data=data);summary(block)

avg_weight <- mean(data$Weight)
SE <- sd(data$Weight)/sqrt(nrow(data))
avg_weight;SE*qt(1-0.025,nrow(data))
