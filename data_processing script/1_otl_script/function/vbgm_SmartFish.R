# PROJECT: FWO PHD - WP1
# Tuan Anh Bui (10/05/2021)

# Function to fit von Bertalanffy growth model (VBGM) to ILVO SmartFish data

vbgm_smartfish <- function(data) {
  
  # Fit Von Bertalanffy Growth Model ----
  # Model formulation====
  # Variables
  L = data$SPE.Length
  A = data$SPA.Age
  
  # Parameters====
  Linf = max(L) 
  K = 0.2  # reasonably low, e.g. K=0.2
  t0 = 0   # Start at age 0
  
  print("fitting model")
  
  # Model formulation====
  VBGM <- L ~ Linf * (1 - exp(-K*(A-t0))) # Von Bertalanffy Equation
  fit_vbgm <- nls(VBGM, start = list (Linf = Linf, K = K, t0 = t0),
                  control = nls.control(maxiter = 1000)) # use the predefined starting values for the VB parameters
  
  # Model output====
  summary(fit_vbgm)
  # Linf: The theoretical maximum size is estimated at 61.2 cm
  # K: the growth coefficient (or Brody coefficient) describes how quickly the Linf is reached
  # t0: the extrapolated and theoretical age at which length of the fish is zero.
  #     t0 is the embryonic age, which should be close to zero for oviparous species with small eggs
  coef(fit_vbgm)
  
  data <- data       %>% mutate(Linf       = summary(fit_vbgm)$parameters[1,1],
                                #Linf_se    = summary(fit_vbgm)$parameters[1,2],
                                #Linf_upper = Linf + 2*Linf_se,
                                #Linf_lower = Linf - 2*Linf_se,
                                
                                K       = summary(fit_vbgm)$parameters[2,1],
                                #K_se    = summary(fit_vbgm)$parameters[2,2],
                                #K_upper = K + 2*K_se,
                                #K_lower = K - 2*K_se,
                                
                                t0       = summary(fit_vbgm)$parameters[3,1],
                                #t0_se    = summary(fit_vbgm)$parameters[3,2],
                                #t0_upper = t0 + 2*t0_se,
                                #t0_lower = t0 - 2*t0_se
                                )
  
  # prediction
  # L ~ Linf * (1 - exp(-K*(A-t0)))
  
  data <- data       %>% mutate(pred_vbgm = Linf*(1 - exp(-K*(SPA.Age - t0))),
                                #pred_vbgm_upper = Linf_upper*(1 - exp(-K_upper*(SPA.Age - t0_upper))),
                                #pred_vbgm_lower = Linf_upper*(1 - exp(-K_lower*(SPA.Age - t0_lower)))
                                )
  
  print("plotting")
  
  # Plot model outputs====
  # parameters
  para_vbgm <- paste("Linf =", round(unique(data$Linf),2), "|",
                     "K =", round(unique(data$K),2), "|",
                     "t0 =",round(unique(data$t0),2), sep = " ")
  
  p <- ggplot() + 
    geom_line(data = data, aes(x = SPA.Age, y = pred_vbgm)) +
    #geom_ribbon(data = data, aes(x = SPA.Age, ymin = pred_vbgm_lower, ymax = pred_vbgm_upper), alpha = 0.3) +
    geom_point(data = data, 
               aes(x = SPA.Age, y = SPE.Length), alpha = 0.1) +
    annotate("text", x = 10, y = 20, label = para_vbgm) +
    labs(title = paste("SOL", unique(data$HAU.IcesAreaGroup), sep = " "), 
         x = "age", 
         y= "length (mm)") +
    theme_bw() # remove grey background
  
  return(p)
  
}
