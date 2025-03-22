# Circuit 1-4 equations - steps

# Function to calculate periodic value
periodic <- function(t, high_time = 1, low_time = 0.5) {
  period <- high_time + low_time
  t <- t %% period
  if (t < high_time) {
    return(1)
  } else if (t < period) {
    return(0)
  }
}

# Parameters
parameters_default <- c(alpha_TF = 1, beta_TF=1, K_H = 0.35, 
                        alpha_Target = 1, K_TF = 0.75, beta_Target=1, 
                        N=4, A_Target = 1,
                        steptime=0,
                        osctime=15,
                        circuit_number = 1,
                        circuit_independent = FALSE
)

state <- c(TF = 0, Target = 0)

H_step <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    # Set H =0 before steptime, H = 1 between steptime and osctime, and H = periodic(t, 1, .5) after osctime
    if (t<steptime) H=0
    else {
      if (t<osctime)
        H=1
      else
        H = periodic(t)
    }
    # Calculate dTF and dTarget based on circuit_number and circuit_independent logic value and related equations
    if (circuit_number == 1) {
      dTF = (H/K_H)^N/((H/K_H)^N+1)* beta_TF - alpha_TF * TF
      if (circuit_independent == FALSE) {
        dTarget = beta_Target * ((TF/K_TF)^N/((TF/K_TF)^N+1) + A_Target )/(1 + (H/K_H)^N) - alpha_Target * Target 
      } else {dTarget = beta_Target * ((TF/K_TF)^N/((TF/K_TF)^N+1) + A_Target /(1 + (H/K_H)^N)) - alpha_Target * Target }
      
    } else if (circuit_number == 2) {
      dTF = (beta_TF)/((H/K_H)^N+1) - alpha_TF * TF
      if (circuit_independent == FALSE) {
        dTarget = beta_Target * (((H/K_H)^N/((H/K_H)^N+1))*(TF/K_TF)^N/((TF/K_TF)^N+1) + A_Target ) - alpha_Target * Target
      } else {dTarget = beta_Target * (((H/K_H)^N*A_Target/((H/K_H)^N+1))+(TF/K_TF)^N/((TF/K_TF)^N+1)) - alpha_Target * Target}  
    } else if (circuit_number == 3) {
      dTF = (H/K_H)^N/((H/K_H)^N+1)* beta_TF - alpha_TF * TF
      if (circuit_independent == FALSE) {
      dTarget = beta_Target * ((H/K_H)^N/((H/K_H)^N+1) + A_Target )/(1 + (TF/K_TF)^N) - alpha_Target * Target
      } else {dTarget = beta_Target * ((H/K_H)^N*A_Target/((H/K_H)^N+1) + 1/((TF/K_TF)^N+1) ) - alpha_Target * Target}
    } else if (circuit_number == 4) {
      dTF = (beta_TF)/((H/K_H)^N+1) - alpha_TF * TF
      if (circuit_independent == FALSE) {
      dTarget = beta_Target*(1/(((H/K_H)^N+1)*((TF/K_TF)^N+1))+A_Target) - alpha_Target * Target
      } else {dTarget = beta_Target*(A_Target/((H/K_H)^N+1)+1/((TF/K_TF)^N+1)) - alpha_Target * Target}
    } 
    
    list(c(dTF, dTarget))
  })
}
