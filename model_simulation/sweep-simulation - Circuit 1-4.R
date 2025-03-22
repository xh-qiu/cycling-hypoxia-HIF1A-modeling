#*************     Notes     ************
#** 1.Set the circuit number at line 19.*
#** 2.Set regular or independent condition at line 20.*

library(data.table)
library(ggplot2)
library(deSolve)
library(officer)
library(rvg)


# load the equations file
source("Circuit 1-4 equations - steps.R")
# load the helper functions
source("helper-funcs.R")

# set a circuit number & circuit_independent logical
parameters = parameters_default
parameters["circuit_number"] = 1
parameters["circuit_independent"]=TRUE

resdf = data.frame(K_TF = numeric(), K_H = numeric(), H_mean_Target = numeric(), N_mean_Target = numeric(), O_mean_Target = numeric())
plot_data = data.frame(K_TF = numeric(), K_H = numeric(), Value = numeric(), Condition = character())

for (K_TF_this in seq(from = 0.01, to = 2, length.out=20))  {
  for (K_H_this in seq(from =0.01, to = 2, length.out=20)) {

    # change the parameters
    
    parameters["K_TF"] = K_TF_this
    parameters["K_H"] = K_H_this
    
    # simulate the ODE model
    times <-seq(-15,30,by=.05)
    Hstepout = as.data.table(ode(y = state, parms = parameters, 
                             func = H_step, 
                             times))
    # add the column for HIF1A
    Hstepout[, HIF1A := if (time<parameters["steptime"]) 0 else {if (time<parameters["osctime"]) 1 else periodic(time)}, by=time]
    
    # calculate average Target values for H, N, and O
    N_mean_Target = mean(Hstepout[(time< 0 ) & (time > - 9), Target])
    H_mean_Target = mean(Hstepout[(time > 6 ) & (time < 15 ), Target])
    O_mean_Target = mean(Hstepout[(time > 21) & (time < 30 ), Target])
    
    # append to resdf
    resdf = rbind(resdf, data.frame(K_TF = K_TF_this, K_H = K_H_this, H_mean_Target = H_mean_Target, N_mean_Target = N_mean_Target, O_mean_Target = O_mean_Target))
    plot_data = rbind(plot_data, data.frame(K_TF = K_TF_this, K_H = K_H_this, Value = H_mean_Target, Condition = "H"))
    plot_data = rbind(plot_data, data.frame(K_TF = K_TF_this, K_H = K_H_this, Value = N_mean_Target, Condition = "N"))
    plot_data = rbind(plot_data, data.frame(K_TF = K_TF_this, K_H = K_H_this, Value = O_mean_Target, Condition = "O"))
  }
}

resdf = as.data.table(resdf)

resdf[, HNOorder := getHNOorder(H_mean_Target, N_mean_Target, O_mean_Target)]

# Create the plot
# Set the file names based on the circuit number and condition
if (parameters["circuit_independent"] == FALSE)
  {file_name01 <- paste0("Circuit ", parameters["circuit_number"] ," regular sweep phase plots.tiff")
  file_name02 <- paste0("Circuit ", parameters["circuit_number"] ," regular sweep phase plots.pptx")
  file_name <- paste0("resdf-circuit", parameters["circuit_number"], " regular sweep phase plots.xlsx")
} else {file_name01 <- paste0("Circuit ", parameters["circuit_number"] ," independent sweep phase plots.tiff")
        file_name02 <- paste0("Circuit ", parameters["circuit_number"] ," independent sweep phase plots.pptx")
        file_name <- paste0("resdf-circuit", parameters["circuit_number"], " independent sweep phase plots.xlsx")}

# Save the plot as a TIFF file
tiff(file_name01, width = 7, height = 6, units = "in", res = 300, compression = "lzw+p")

# Tile plot with HNOorder
p_HNO <- ggplot(resdf, aes(x = K_H, y = K_TF, fill = HNOorder)) +
  geom_tile(  ) + theme_minimal() + scale_fill_brewer(palette = "Paired") + theme(axis.text = element_text(colour = "black"))

# Print the plot to the tiff device
print(p_HNO)

# source plot to ppt function
source("add_plot_to_ppt_function.R")

# Create PPT object
ppt <- read_pptx()

# Add the third plot to a slide
ppt <- add_plot_to_ppt(ppt, p_HNO)

print(ppt, target = file_name02)
openxlsx::write.xlsx(resdf, file_name)

dev.off() 