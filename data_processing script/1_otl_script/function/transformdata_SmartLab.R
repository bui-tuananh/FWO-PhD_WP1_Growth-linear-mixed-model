# ********************************
# OTOLITH SMARTLAB DATA PROCESS **
# ********************************

# INFO 

# PROJECT: FWO PHD - WP1
# Tuan Anh Bui (06/05/2022)

# Function to transform SmartLab/SmartDots data into AnnulusDiameter and AnnulusDiameterIncrement

# Input: data from SmartLab/SmartDots 
# Output: AnnulusDiameter, AnnulusDiameterIncrement

# Note: The XY in SmartDots seems to be from the top-left, while the XY in R is bottom-left
# this only influences how we plot the points, not the calculation of the OtolithWidth

SmartDots_to_OtolithIncrement <- function(otl_data) {
  
  # Find X, Y of projected dots on the linear line of 2 edge dots
  # formula: (Xp, Yp) = [(X,Y)*v]*v 
  # Xp, Yp - X, Y projected dots
  # X, Y   - X,Y original dots
  # v      - transforming vector - v = (1/sqrt(1 + Intercept^2))*[1; Intercept] 
  
  # more explanation of the method see: https://math.stackexchange.com/questions/3020773/determine-position-of-projected-points-onto-a-line
  
  # workflow
  # 1. find linear line Y = intercept + slope*X of 2 edge dots (dots with largest DotIndex)
  # 2. calculate transforming vector
  # 3. calculate the projected Xp, Yp
  # 4. calculate the distance - AnnulusDiameter - on the linear line of 2 edge dots
  
  # list of Sample_ID - each otolith (sample) has an unique Sample_ID
  list_otl <- sort(unique(otl_data$Sample_ID)) 
  
  # create an output empty dataframe
  otl_increment <- otl_data %>% select(Sample_ID, Outcome_LabTechnician, DotIndex)
  otl_increment <- otl_increment[FALSE,]
  
  for (i in 1:length(list_otl)) {
    
    print(paste0("Processing sample ", list_otl[i]))
    
    # 1. find linear line Y = intercept + slope*X of 2 edge dots (dots with largest DotIndex)
    
    # change otl_temp to otb_data_sub
    otb_data_sub <- otl_data %>% filter(Sample_ID == list_otl[i]) %>% select(Sample_ID, Outcome_LabTechnician, Outcome_Comment, X1:DotIndex)
    
    # process for each lab technician
    list_LabTechnician <- unique(otb_data_sub$Outcome_LabTechnician)
    
    for (l in 1:length(list_LabTechnician)) {
      otl_temp <- otb_data_sub %>% filter(Outcome_LabTechnician == list_LabTechnician[l])
      
      
      # create otl_maxwidth - width of otolith and to calculate line 
      otl_maxwidth <- otl_temp %>% select(Outcome_Comment, X2, Y2)                
      
      if (any(is.na(otl_maxwidth$Outcome_Comment))  == T) { #if the otolith lacks Outcome_Comment, stop the loop
        
        print(paste("SampleID", list_otl[i], "lacks Outcome_Comment"), sep = " ")
        break # exit loop
      
      }
      
      otl_maxwidth <- otl_maxwidth[FALSE,]
      
      list_edge <- unique(otl_temp$Outcome_Comment)
      
      for (e in 1:length(list_edge)) {
        otl_edge <- otl_temp %>% filter(Outcome_Comment == list_edge[e])
        
        # select the last row which has information of the last line in SmartDots: start point (X1,Y1), end point (X2,Y2)
        # we want the endpoint, which is the furthest point in the otolith edge
        otl_edge <- otl_edge[nrow(otl_edge),] %>% select(Outcome_Comment, X2, Y2) 
        
        otl_maxwidth <- rbind(otl_maxwidth, otl_edge) 
        
      }
      
      lm <- lm(Y2 ~ X2, data = otl_maxwidth)
      
      intercept <- lm$coefficients[1]   # intercept
      slope     <- lm$coefficients[2]   # slope
      
      # 2. calculate transforming vector v
      v0 <- c(1, slope)
      v <- (1/sqrt(1 + slope^2))*v0
      
      # remove X1-Y2 to keep coordinates of annotation (X,Y) only
      otl_temp <- otl_temp %>% select(Sample_ID, Outcome_LabTechnician, Outcome_Comment, X, Y, DotIndex)
      
      # 3. calculate the projected Xp, Yp
      otl_temp <- otl_temp %>% mutate(vX = v[1],
                                      vY = v[2],
                                      XYv = X*vX + (Y - intercept)*vY, # (X,Y)*v
                                      Xp = XYv*vX,                     # [(X,Y)*v]*v
                                      Yp = XYv*vY + intercept)         # [(X,Y)*v]*v
      
      ggplot() + 
        geom_point(data = otl_temp, aes(x = Xp, y = Yp)) +
        geom_point(data = otl_temp, aes(x = X, y = Y, color = as.factor(DotIndex)), alpha = 0.5) +
        geom_point(data = otl_maxwidth, aes(x = X2, y = Y2), shape = 2, color = "red") + #points at 2 edges
        xlim(0, max(otl_temp$X) + 100) +
        ylim(0, max(otl_temp$Y) + 100) +
        theme_bw()
      
      # 4. calculate the distance - AnnulusDiameter - on the linear line of 2 edge dots
      # using the law of sine
      # AnnulusDiameter*sin(90deg) = (Xp2-Xp1)*(sin(90deg - i_degree)
      
      i_degree = abs(atan(slope)/(pi/180))
      
      # create an empty dataframe and list of DotIndex
      otl_diameter <- otl_temp[FALSE,]
      
      list_DotIndex <- unique(otl_temp$DotIndex)
      
      for (s in 1:length(list_DotIndex)) {
        otl_index <- otl_temp %>% filter(DotIndex == list_DotIndex[s])
        otl_index <- unique(otl_index)
        
        otl_diameter_ori <- abs(otl_index$Xp[1] - otl_index$Xp[2])
        otl_diameter_proj <- otl_diameter_ori/sin((90-i_degree)*pi/180)
        
        otl_index$AnnulusDiameter.pixel <- otl_diameter_proj
        
        otl_diameter <- rbind(otl_diameter, otl_index)
        
      }
      
      otl_diameter <- otl_diameter %>% select(Sample_ID, Outcome_LabTechnician, DotIndex, AnnulusDiameter.pixel)
      otl_diameter <- unique(otl_diameter)
      
      # Calculate otolith increment - using diff function
      otl_diameter$AnnulusDiameterIncrement.pixel <- c(otl_diameter$AnnulusDiameter.pixel[1], diff(otl_diameter$AnnulusDiameter.pixel))
      
      # Calculate otolith width using the edges
      otl_width_ori <- abs(otl_maxwidth$X2[1] - otl_maxwidth$X2[2])
      otl_width_proj <- otl_width_ori/sin((90-i_degree)*pi/180)
      
      otl_diameter$OtolithWidth.pixel <- otl_width_proj
      
      # binding otl_width into otl
      otl_increment <- rbind(otl_increment, otl_diameter)
    }
    
    
  }
  
  # join otl_increment to the otolith data from SmartLab
  otl  <- otl_data %>% select(-(X1:Y),
                              -Outcome_ID, 
                              -Outcome_Comment, 
                              -Outcome_DateCreation)
  otl <- unique(otl)
  
  otl <- left_join(otl, otl_increment, by = c("Sample_ID", "Outcome_LabTechnician", "DotIndex"))
  
  # convert pixel to mm (File_Scale: X pixel/1000 um)
  otl <- otl %>% mutate(File_Scale                  = as.numeric(File_Scale),
                        AnnulusDiameter.um          = AnnulusDiameter.pixel/File_Scale*1000,
                        AnnulusDiameterIncrement.um = AnnulusDiameterIncrement.pixel/File_Scale*1000,
                        OtolithWidth.um             = OtolithWidth.pixel/File_Scale*1000) 
  
  return(otl)
}
