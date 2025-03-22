#*************     Notes     ************
#** 1.Need to prepare the chosen points excel file *
#**   for generating regular or independent time plots *
#** 2.Set the circuit number at line 34.*
#** 3.Set regular or independent condition at line 35.*

library(data.table)
library(ggplot2)
library(deSolve)
library(readxl)

# load the equations file
source("Circuit 1-4 equations - steps.R")
source("helper-funcs.R")

# Regular or independent,choose one ,another to a commnet sign
# Regular chosen points
parameters_data <- read_excel("resdf-circuit1 regular sweep phase plots_chosen points.xlsx")
# Independent chosen points
# parameters_data <- read_excel("resdf-circuit1 independent sweep phase plots_chosen points.xlsx")

# set a circuit number & circuit_independent

# Loop through each row of parameters_data
for (i in 1:nrow(parameters_data)) {
  # Extract K_H, K_TF, and HNOorder for the current row
  K_H_this <- parameters_data$K_H[i]
  K_TF_this <- parameters_data$K_TF[i]
  HNOorder_this <- parameters_data$HNOorder[i]
  
  parameters = parameters_default
  parameters["K_TF"] = K_TF_this
  parameters["K_H"] = K_H_this
  parameters["circuit_number"] = 1
  parameters["circuit_independent"]=FALSE
  
  # plot the time simulation
  times <-seq(-15,30,by=.05)
  Hstepout = as.data.table(ode(y = state, parms = parameters, 
                               func = H_step, 
                               times = times))
  
  # add the column for HIF1A
  Hstepout[, HIF1A := if (time<parameters["steptime"]) 0 else {if (time<parameters["osctime"]) 1 else periodic(time)}, by=time]
  Hstepoutdt_melt = melt.data.table(data = Hstepout, id.vars = "time")
  Hstepoutdt_melt[, Condition:= if (time<parameters["steptime"]) "Normoxia" else {if (time<parameters["osctime"]) "Hypoxia" else "Oscillation"}, by=time]
  Hstepoutdt_melt$Condition = factor(Hstepoutdt_melt$Condition, levels = c("Normoxia", "Hypoxia", "Oscillation"))
  
  # re-level the factor Hstepoutdt_melt$variable to have HIF1A first
  Hstepoutdt_melt$variable = factor(Hstepoutdt_melt$variable, levels = c("HIF1A", "TF", "Target"))
  
  # Calculate means for time > osctime
  means <- Hstepoutdt_melt[time > 21, .(mean_value = mean(value)), by = variable]
  
  # Clean the HNOorder string for the file name
  clean_HNOorder <- gsub("[<>]", "_", HNOorder_this)
  
  # Save the plot as a TIFF file
  if (parameters["circuit_independent"] == FALSE)
  {file_name <- paste0("Circuit ", parameters["circuit_number"]," K_H(",K_H_this,")", " K_TF(",K_TF_this,") ", clean_HNOorder, " regular time plot.tiff")
  } else { file_name <- paste0("Circuit ", parameters["circuit_number"]," K_H(",K_H_this,")", " K_TF(",K_TF_this,") ", clean_HNOorder, " independent time plot.tiff") }
  tiff(file_name, width = 4.5, height = 2.5, units = "in", res = 300, compression = "lzw+p")
  
  p<-ggplot(data = Hstepoutdt_melt[((time > -9)&(time <0))|
                                     ((time >6)&(time <15))|
                                     ((time >21))], 
            mapping = aes(x=time, y = value)) +
    geom_line(mapping = aes(color=Condition),linewidth=0.6) +
    facet_grid( variable ~ ., scales = "free_y") +
    theme_classic() +
    theme(axis.text = element_text(color="black"),  
          axis.title = element_text(size=14),                
          plot.title = element_text(size=12,hjust = 0.3)) +  
    
    # Add horizontal lines for means in time > osctime
    geom_segment(data = means, aes(x = 21, xend = max(Hstepoutdt_melt$time), 
                                   y = mean_value, yend = mean_value), linetype = "dashed", color = "black") +
    # Add text annotations for means
    geom_text(data = means, aes(x = max(Hstepoutdt_melt$time), y = mean_value, label = round(mean_value, 2)), 
              vjust = -0.15, hjust = 1.1, color = "black", size = 4) +
    
    # Add title
    ggtitle(if (parameters["circuit_independent"] == FALSE) {
    paste0("Circuit ",parameters["circuit_number"]," regular / K_H=",round(K_H_this,3), " K_TF=",round(K_TF_this,3)," / ",HNOorder_this)
    } else {paste0("Circuit ",parameters["circuit_number"]," independent / K_H=",round(K_H_this,3), " K_TF=",round(K_TF_this,3)," / ",HNOorder_this,systemtime())})
  
  print(p)
  dev.off()
}