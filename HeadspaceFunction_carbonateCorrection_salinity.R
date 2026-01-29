#####################################################################
# Rheadspace_GHG.R
# with carbonate correction and salinity correction
#

# R function that calculates dissolved pCO2, pCH4, and pN2O from raw GC data from headspace
# equilibration.  For pCO2, the complete headspace output accounts for carbonate
# equilbrium according to Rheadspace.R by Marce, Kim, and Prairie as in 
# Koschorreck, M., Y.T. Prairie, J. Kim, and R. Marcé. 2020. Technical note: CO2 is not like CH4 – limits of the headspace method to analyse pCO2 in water. Biogeosciences, 2021
# 
# This function modifies the Rheadspace function by including pCH4 and pN2O processing using 
# CH4 and N2O Henry's Law constants according to Sander 2015 for freshwater option and
# Wiesenburg and Guinasso 1979 and Weiss and Price 1980 (for estuarine and marine options).
#

#####################################################################

Rheadspace <-  function(...){
  arguments <- list(...)
  
  # test arguments and initialize variables
  if (is.data.frame(arguments[[1]])) {
    input.table=arguments[[1]]
    if (dim(input.table)[2]!=16){
      stop("You should input a data frame with 16 columns. See the readme file or comments in the function", call.=FALSE)
    }else{
      Site = as.character(input.table$site)
      Timestamp = input.table$datetime.EST
      mCH4_headspace = input.table$HS.mCH4.before #the mCH4 (ppmv) of the headspace "before" equilibration
      mCO2_headspace = input.table$HS.mCO2.before #the mCO2 (ppmv) of the headspace "before" equilibration
      mN2O_headspace = input.table$HS.mN2O.before #the mN2O (ppmv) of the headspace "before" equilibration
      mCH4_eq = input.table$HS.mCH4.after #the measured mCH4 (ppmv) of the headspace "after" equilibration
      mCO2_eq = input.table$HS.mCO2.after #the measured mCO2 (ppmv) of the headspace "after" equilibration
      mN2O_eq = input.table$HS.mN2O.after #the measured mN2O (ppmv) of the headspace "after" equilibration
      temp_insitu = input.table$Temp.insitu #in situ water temperature in degrees celsius
      temp_eq = input.table$Temp.equil #the water temperature after equilibration in degree celsius
      alk = input.table$Alkalinity.measured #Total alkalinity (micro eq/L) of the water sample
      vol_gas = input.table$Volume.gas #Volume of gas in the headspace vessel (mL)
      vol_water = input.table$Volume.water #Volume of water in the headspace vessel (mL)   
      Bar.pressure = input.table$Bar.pressure #Barometric pressure at field conditions in kPa. 101.325 kPa = 1 atm   
      c_constants = input.table$Constants #Constants for carbonate equilibrium (1=Freshwater; 2=Estuarine; 3=Marine) 
      Salinity = input.table$Salinity #Salinity in PSU. 
    } 
  } else if (length(arguments)==16) {
    Site = as.character(arguments[[1]])
    Timestamp = arguments[[2]]
    mCH4_headspace = arguments[[3]] #the mCH4 (ppmv) of the headspace "before" equilibration
    mCO2_headspace = arguments[[4]] #the mCO2 (ppmv) of the headspace "before" equilibration
    mN2O_headspace = arguments[[5]] #the mN2O (ppmv) of the headspace "before" equilibration
    mCH4_eq = arguments[[6]] #the measured mCH4 (ppmv) of the headspace "after" equilibration
    mCO2_eq = arguments[[7]] #the measured mCO2 (ppmv) of the headspace "after" equilibration
    mN2O_eq = arguments[[8]] #the measured mN2O (ppmv) of the headspace "after" equilibration
    temp_insitu = arguments[[9]] #in situ water temperature in degrees celsius
    temp_eq = arguments[[10]] #the water temperature after equilibration in degree celsius
    alk = arguments[[11]] #Total alkalinity (micro eq/L) of the water sample
    vol_gas = arguments[[12]] #Volume of gas in the headspace vessel (mL)
    vol_water = arguments[[13]] #Volume of water in the headspace vessel (mL)   
    Bar.pressure = arguments[[14]] #Barometric pressure at field conditions in kPa. 101.325 kPa = 1 atm   
    c_constants = arguments[[15]] #Constants for carbonate equilibrium (1=Freshwater; 2=Estuarine; 3=Marine) 
    Salinity = arguments[[16]] #Salinity in PSU. Set to zero if Constants = 1
  } else {
    stop("You should input either a data frame or a vector of 16 values. See the readme file or comments in the function", call.=FALSE)
  }
  
  
  
  #initialization of variables
  pGHG_orig <- data.frame(matrix(NA,length(mCO2_headspace),17))
  names(pGHG_orig) <- c("Site",
                        "Timestamp",
                        "mCO2 complete headspace (ppmv)",
                        "pCO2 complete headspace (micro-atm)", 
                        "CO2 concentration complete headspace (micro-mol/L)", 
                        "pH", 
                        "mCO2 simple headspace (ppmv)", 
                        "pCO2 simple headspace (micro-atm)", 
                        "CO2 concentration simple headspace (micro-mol/L)", 
                        "% error",
                        "mCH4 simple headspace (ppmv)", 
                        "pCH4 simple headspace (micro-atm)", 
                        "CH4 concentration simple headspace (micro-mol/L)",
                        "mN2O simple headspace (ppmv)", 
                        "pN2O simple headspace (micro-atm)", 
                        "N2O concentration simple headspace (micro-mol/L)",
                        "constants")
  
  R <- 0.082057338 #L atm K-1 mol-1
  
  #the function uniroot cannot handle vectors, so we need a loop
  for (i in 1:length(mCO2_headspace)){ 
    
    AT = alk[i]*(1e-6) #conversion to mol/L
    
    #Constants of the carbonate ewuilibrium
    # Kw = the dissociation constant of H2O into H+ and OH-
    # Kh = the solubility of CO2 in water - equilibration conditions
    # Kh2 = the solubility of CO2 in water - in situ field conditions
    # K1 = the equilibrium constant between CO2 and HCO3-
    # K2 = the equilibrium constant between HCO3- and CO3 2-
    
    # Solubility coefficients from Weiss (1974) with Sal=0 for freshwater option
    # Dissociation of water from Dickson and Riley (1979)
    
    if (c_constants[i] == 1) {
      
      #Millero, F. (1979). The thermodynamics of the carbonate system in seawater
      #Geochimica et Cosmochimica Acta 43(10), 1651 1661.  
      K1=10^-(-126.34048+6320.813/(temp_eq[i]+273.15)+19.568224*log(temp_eq[i]+273.15))
      K2=10^-(-90.18333+5143.692/(temp_eq[i]+273.15)+14.613358*log(temp_eq[i]+273.15))
      Kw = exp(148.9652-13847.26/(temp_eq[i]+273.15)-23.6521*log(273.15+temp_eq[i]))
      
      #Weiss contants for mol/kg·atm instead of mol/L·atm were used in the first version of the code by mistake. Correcting now to the right values in mol/L·atm
      #Kh = 10^((-60.2409+93.4517*(100/(273.15+temp_eq[i]))+23.3585*log((273.15+temp_eq[i])/100))/log(10)) # mol/L/atm equilibration conditions
      #Kh2 = 10^((-60.2409+93.4517*(100/(273.15+temp_insitu[i]))+23.3585*log((273.15+temp_insitu[i])/100))/log(10)) # mol/L/atm original conditions
      Kh.co2 = 10^((-58.0931+90.5069*(100/(273.15+temp_eq[i]))+22.2940*log((273.15+temp_eq[i])/100))/log(10)) # mol/L/atm equilibration conditions
      Kh2.co2 = 10^((-58.0931+90.5069*(100/(273.15+temp_insitu[i]))+22.2949*log((273.15+temp_insitu[i])/100))/log(10)) # mol/L/atm original conditions
      Kh.ch4 = 0.000014*exp(1900*(1/(273.15+temp_eq[i])-1/298.15))*101325/1000# mol/L/atm equilibration conditions
      Kh2.ch4 = 0.000014*exp(1900*(1/(273.15+temp_insitu[i])-1/298.15))*101325/1000 # mol/L/atm original conditions
      Kh.n2o = 0.00024*exp(2700*(1/(273.15+temp_eq[i])-1/298.15))*101325/1000 # mol/L/atm equilibration conditions
      Kh2.n2o = 0.00024*exp(2700*(1/(273.15+temp_insitu[i])-1/298.15))*101325/1000 # mol/L/atm original conditions
      
      
    } else if (c_constants[i] == 2) {
      
      #Millero, F. (2010). Carbonate constants for estuarine waters Marine and Freshwater
      #Research 61(2), 139. As amended by Orr et al. 2015.
      pK10=(-126.34048+6320.813/(temp_eq[i]+273.15)+19.568224*log(temp_eq[i]+273.15))
      A1 = 13.4038*Salinity[i]^0.5 + 0.03206*Salinity[i] - 5.242e-5*Salinity[i]^2
      B1 = -530.659*Salinity[i]^0.5 - 5.8210*Salinity[i]
      C1 = -2.0664*Salinity[i]^0.5
      pK1 = pK10 + A1 + B1/(temp_eq[i]+273.15) + C1*log(temp_eq[i]+273.15)
      K1 = 10^-pK1;
      pK20=(-90.18333+5143.692/(temp_eq[i]+273.15)+14.613358*log(temp_eq[i]+273.15))
      A2 = 21.3728*Salinity[i]^0.5 + 0.1218*Salinity[i] - 3.688e-4*Salinity[i]^2
      B2 = -788.289*Salinity[i]^0.5 - 19.189*Salinity[i]
      C2 = -3.374*Salinity[i]^0.5
      pK2 = pK20 + A2 + B2/(temp_eq[i]+273.15) + C2*log(temp_eq[i]+273.15)
      K2 = 10^-pK2;
      
      Kw=exp(148.9652-13847.26/(temp_eq[i]+273.15)-23.6521*log(273.15+temp_eq[i])+sqrt(Salinity[i])*(118.67/(temp_eq[i]+273.15)-5.977+1.0495*log(273.15+temp_eq[i]))-0.01615*Salinity[i])
      
      #Weiss contants for mol/kg·atm instead of mol/L·atm were used in the first version of the code by mistake. Correcting now to the right values in mol/L·atm
      #Kh = exp(-60.2409 + 93.4517 * (100 / (273.15 + temp_eq[i])) + 23.3585 * log( (273.15 + temp_eq[i]) / 100 )+Salinity[i]*(0.023517-0.023656*(273.15+temp_eq[i])/100+0.0047036*((273.15+temp_eq[i])/100)^2))
      #Kh2 = exp(-60.2409 + 93.4517 * (100 / (273.15 + temp_insitu[i])) + 23.3585 * log( (273.15 + temp_insitu[i]) / 100 )+Salinity[i]*(0.023517-0.023656*(273.15+temp_insitu[i])/100+0.0047036*((273.15+temp_insitu[i])/100)^2))
      Kh.co2 = exp(-58.0931 + 90.5069 * (100 / (273.15 + temp_eq[i])) + 22.2940 * log( (273.15 + temp_eq[i]) / 100 )+Salinity[i]*(0.027766-0.025888*(273.15+temp_eq[i])/100+0.0050578*((273.15+temp_eq[i])/100)^2))
      Kh2.co2 = exp(-58.0931 + 90.5069 * (100 / (273.15 + temp_insitu[i])) + 22.2940 * log( (273.15 + temp_insitu[i]) / 100 )+Salinity[i]*(0.027766-0.025888*(273.15+temp_insitu[i])/100+0.0050578*((273.15+temp_insitu[i])/100)^2))
      
      #Constants for Bunsen coefficient CH4 from Wiesenburg and Guinasso (1979) see Wanninkhof 2014
      ch4A1 = -68.8862
      ch4A2 = 101.4956
      ch4A3 = 28.7314
      ch4B1 = -0.076146
      ch4B2 = 0.04397
      ch4B3 = -0.006872
      
      #calculate  Bunsen coefficient (B.ch4), density and K0 (Kh.ch4) based on temp and salinity
      B.ch4 = exp(ch4A1+ch4A2*(100/(273.15+temp_eq[i]))+ch4A3*log((273.15+temp_eq[i])/100)+Salinity[i]*(ch4B1+ch4B2*((273.15+temp_eq[i])/100)+ch4B3*((273.15+temp_eq[i])/100)^2))# mol/L/atm equilibration conditions
      B.ch4_density = ((999.842594+0.06793952*temp_eq[i]-0.00909529*temp_eq[i]^2+0.0001001685*temp_eq[i]^3-0.000001120083*temp_eq[i]^4+0.000000006536332*temp_eq[i]^5)+(0.824493-0.0040899*temp_eq[i]+0.000076438*temp_eq[i]^2-0.00000082467*temp_eq[i]^3+0.0000000053875*temp_eq[i]^4)*Salinity[i]+(-0.00572466+0.00010227*temp_eq[i]-0.0000016546*temp_eq[i]^2)*Salinity[i]^1.5+(0.00048314)*Salinity[i]^2)/1000
      Kh.ch4 = B.ch4/(B.ch4_density*22.413)
      
      B2.ch4 = exp(ch4A1+ch4A2*(100/(273.15+temp_insitu[i]))+ch4A3*log((273.15+temp_insitu[i])/100)+Salinity[i]*(ch4B1+ch4B2*((273.15+temp_insitu[i])/100)+ch4B3*((273.15+temp_insitu[i])/100)^2)) # mol/L/atm original conditions
      B2.ch4_density = ((999.842594+0.06793952*temp_insitu[i]-0.00909529*temp_insitu[i]^2+0.0001001685*temp_insitu[i]^3-0.000001120083*temp_insitu[i]^4+0.000000006536332*temp_insitu[i]^5)+(0.824493-0.0040899*temp_insitu[i]+0.000076438*temp_insitu[i]^2-0.00000082467*temp_insitu[i]^3+0.0000000053875*temp_insitu[i]^4)*Salinity[i]+(-0.00572466+0.00010227*temp_insitu[i]-0.0000016546*temp_insitu[i]^2)*Salinity[i]^1.5+(0.00048314)*Salinity[i]^2)/1000
      Kh2.ch4 = B2.ch4/(B2.ch4_density*22.413)
      
      # Constants for solubility K0 for N2O from Weiss and Price (1980) see Wanninkhof 2014
      n2oA1 = -62.7062
      n2oA2 = 97.3066
      n2oA3 = 24.1406
      n2oB1 = -0.058420
      n2oB2 = 0.033193
      n2oB3 = -0.0051313
      
      Kh.n2o = exp(n2oA1+n2oA2*(100/(273.15+temp_eq[i]))+n2oA3*log((273.15+temp_eq[i])/100)+Salinity[i]*(n2oB1+n2oB2*((273.15+temp_eq[i])/100)+n2oB3*((273.15+temp_eq[i])/100)^2))# mol/L/atm equilibration conditions
      Kh2.n2o = exp(n2oA1+n2oA2*(100/(273.15+temp_insitu[i]))+n2oA3*log((273.15+temp_insitu[i])/100)+Salinity[i]*(n2oB1+n2oB2*((273.15+temp_insitu[i])/100)+n2oB3*((273.15+temp_insitu[i])/100)^2)) # mol/L/atm original conditions
      
    } else if (c_constants[i] == 3) {
      
      #Dickson, A. G., Sabine, C. L., and Christian, J. R. (2007): Guide to best practices for
      #ocean CO2 measurements, PICES Special Publication 3, 191 pp.
      K1=10^(-3633.86/ (temp_eq[i] + 273.15)+61.2172-9.67770*log(temp_eq[i]+273.15)+0.011555*Salinity[i]-0.0001152*Salinity[i]^2)
      K2 = 10^(-417.78/ (temp_eq[i] + 273.15) - 25.9290 + 3.16967*log(temp_eq[i]+273.15)+0.01781*Salinity[i]-0.0001112*Salinity[i]^2)
      
      Kw=exp(148.9652-13847.26/(temp_eq[i]+273.15)-23.6521*log(273.15+temp_eq[i])+sqrt(Salinity[i])*(118.67/(temp_eq[i]+273.15)-5.977+1.0495*log(273.15+temp_eq[i]))-0.01615*Salinity[i])
      
      #Weiss contants for mol/kg·atm instead of mol/L·atm were used in the first version of the code by mistake. Correcting now to the right values in mol/L·atm
      #Kh = exp(-60.2409 + 93.4517 * (100 / (273.15 + temp_eq[i])) + 23.3585 * log( (273.15 + temp_eq[i]) / 100 )+Salinity[i]*(0.023517-0.023656*(273.15+temp_eq[i])/100+0.0047036*((273.15+temp_eq[i])/100)^2))
      #Kh2 = exp(-60.2409 + 93.4517 * (100 / (273.15 + temp_insitu[i])) + 23.3585 * log( (273.15 + temp_insitu[i]) / 100 )+Salinity[i]*(0.023517-0.023656*(273.15+temp_insitu[i])/100+0.0047036*((273.15+temp_insitu[i])/100)^2))
      Kh.co2 = exp(-58.0931 + 90.5069 * (100 / (273.15 + temp_eq[i])) + 22.2940 * log( (273.15 + temp_eq[i]) / 100 )+Salinity[i]*(0.027766-0.025888*(273.15+temp_eq[i])/100+0.0050578*((273.15+temp_eq[i])/100)^2))
      Kh2.co2 = exp(-58.0931 + 90.5069 * (100 / (273.15 + temp_insitu[i])) + 22.2940 * log( (273.15 + temp_insitu[i]) / 100 )+Salinity[i]*(0.027766-0.025888*(273.15+temp_insitu[i])/100+0.0050578*((273.15+temp_insitu[i])/100)^2))
      
      #Constants for Bunsen coefficient CH4 from Wiesenburg and Guinasso (1979) see Wanninkhof 2014
      ch4A1 = -68.8862
      ch4A2 = 101.4956
      ch4A3 = 28.7314
      ch4B1 = -0.076146
      ch4B2 = 0.04397
      ch4B3 = -0.006872
      
      #calculate  Bunsen coefficient (B.ch4), density and K0 (Kh.ch4) based on temp and salinity
      B.ch4 = exp(ch4A1+ch4A2*(100/(273.15+temp_eq[i]))+ch4A3*log((273.15+temp_eq[i])/100)+Salinity[i]*(ch4B1+ch4B2*((273.15+temp_eq[i])/100)+ch4B3*((273.15+temp_eq[i])/100)^2))# mol/L/atm equilibration conditions
      B.ch4_density = ((999.842594+0.06793952*temp_eq[i]-0.00909529*temp_eq[i]^2+0.0001001685*temp_eq[i]^3-0.000001120083*temp_eq[i]^4+0.000000006536332*temp_eq[i]^5)+(0.824493-0.0040899*temp_eq[i]+0.000076438*temp_eq[i]^2-0.00000082467*temp_eq[i]^3+0.0000000053875*temp_eq[i]^4)*Salinity[i]+(-0.00572466+0.00010227*temp_eq[i]-0.0000016546*temp_eq[i]^2)*Salinity[i]^1.5+(0.00048314)*Salinity[i]^2)/1000
      Kh.ch4 = B.ch4/(B.ch4_density*22.413)
      
      B2.ch4 = exp(ch4A1+ch4A2*(100/(273.15+temp_insitu[i]))+ch4A3*log((273.15+temp_insitu[i])/100)+Salinity[i]*(ch4B1+ch4B2*((273.15+temp_insitu[i])/100)+ch4B3*((273.15+temp_insitu[i])/100)^2)) # mol/L/atm original conditions
      B2.ch4_density = ((999.842594+0.06793952*temp_insitu[i]-0.00909529*temp_insitu[i]^2+0.0001001685*temp_insitu[i]^3-0.000001120083*temp_insitu[i]^4+0.000000006536332*temp_insitu[i]^5)+(0.824493-0.0040899*temp_insitu[i]+0.000076438*temp_insitu[i]^2-0.00000082467*temp_insitu[i]^3+0.0000000053875*temp_insitu[i]^4)*Salinity[i]+(-0.00572466+0.00010227*temp_insitu[i]-0.0000016546*temp_insitu[i]^2)*Salinity[i]^1.5+(0.00048314)*Salinity[i]^2)/1000
      Kh2.ch4 = B2.ch4/(B2.ch4_density*22.413)
      
      # Constants for solubility K0 for N2O from Weiss and Price (1980) see Wanninkhof 2014
      n2oA1 = -62.7062
      n2oA2 = 97.3066
      n2oA3 = 24.1406
      n2oB1 = -0.058420
      n2oB2 = 0.033193
      n2oB3 = -0.0051313
      
      Kh.n2o = exp(n2oA1+n2oA2*(100/(273.15+temp_eq[i]))+n2oA3*log((273.15+temp_eq[i])/100)+Salinity[i]*(n2oB1+n2oB2*((273.15+temp_eq[i])/100)+n2oB3*((273.15+temp_eq[i])/100)^2))# mol/L/atm equilibration conditions
      Kh2.n2o = exp(n2oA1+n2oA2*(100/(273.15+temp_insitu[i]))+n2oA3*log((273.15+temp_insitu[i])/100)+Salinity[i]*(n2oB1+n2oB2*((273.15+temp_insitu[i])/100)+n2oB3*((273.15+temp_insitu[i])/100)^2)) # mol/L/atm original conditions
      
      
    } else {
      print(i)
      stop("Option for carbonate equilibrium constants should be a number between 1 and 3", call.=FALSE)
      
    }
    
    HS.ratio <- vol_gas[i]/vol_water[i] #Headspace ratio (=vol of gas/vol of water)
    
    #The following calculations assume 1 atm, this is corrected later for measured pressure in the field.
    
    #DIC at equilibrium
    co2 <- Kh.co2 * mCO2_eq[i]/1000000
    h_all <- polyroot(c(-(2*K1*K2*co2),-(co2*K1+Kw),AT,1))
    real<-Re(h_all)
    h <-real[which(real>0)]
    DIC_eq <- co2 * (1 + K1/h + K1 * K2/(h * h))
    
    #DIC in the original sample
    DIC_ori <- DIC_eq + (mCO2_eq[i] - mCO2_headspace[i])/1000000/(R*(temp_eq[i]+273.15))*HS.ratio
    
    #pCO2 in the original sample
    h_all <- polyroot(c(-(K1*K2*Kw),K1*K2*AT-K1*Kw-2*DIC_ori*K1*K2,AT*K1-Kw+K1*K2-DIC_ori*K1,AT+K1,1))
    real<-Re(h_all)
    h <-real[which(real>0)]
    
    co2 <- h* (DIC_ori * h * K1/(h * h + K1 * h + K1 * K2)) / K1
    
    pGHG_orig[i,1] <- as.character(Site[i])
    pGHG_orig[i,2] <- Timestamp[i]
    pGHG_orig[i,3] <- co2/Kh2.co2*1000000
    pGHG_orig[i,4] <- co2/Kh2.co2*1000000*Bar.pressure[i]/101.325
    pGHG_orig[i,5] <- co2*1000000
    pGHG_orig[i,6] <- -log10( h )
    
    
    #Calculation not accounting for alkalinity effects and associated error
    
    #concentration and total mass in the water sample assuming ideal gas from the pCO2 measured at the headspace
    CO2_solution <- mCO2_eq[i]/1000000*Kh.co2 #mol/L
    CO2_solution_mass <- CO2_solution * vol_water[i]/1000 #mol
    
    #mass of CO2 in the measured headspace
    final_C_headspace_mass <- mCO2_eq[i]/1000000*(vol_gas[i]/1000) / (R * (temp_eq[i]+273.15)) #mol
    
    mols_headspace <- mCO2_headspace[i]/1000000*(vol_gas[i]/1000)/(R * (temp_eq[i]+273.15)) #mol PV / RT = n
    
    #implication: mass, concentration, and partial pressure of CO2 in the original sample (aount in sample and headspace after equilibration minus original mass in the headspace)
    Sample_CO2_mass <- CO2_solution_mass + final_C_headspace_mass - mols_headspace #mol
    Sample_CO2_conc <- Sample_CO2_mass/(vol_water[i]/1000) #mol/L
    pGHG_orig[i,7] <- Sample_CO2_conc/Kh2.co2*1000000 #ppmv
    pGHG_orig[i,8] <- Sample_CO2_conc/Kh2.co2*1000000*Bar.pressure[i]/101.325 # micro-atm
    pGHG_orig[i,9] <- Sample_CO2_conc*1000000 # micro-mol/L
    #calculation of the error
    pGHG_orig[i,10] <- (pGHG_orig[i,6]-pGHG_orig[i,2])/pGHG_orig[i,2] *100  #%
    
    #concentration and total mass in the water sample assuming ideal gas from the pCH4 measured at the headspace
    CH4_solution <- mCH4_eq[i]/1000000*Kh.ch4 #mol/L
    CH4_solution_mass <- CH4_solution * vol_water[i]/1000 #mol
    
    #mass of CH4 in the measured headspace
    final_C_headspace_mass.ch4 <- mCH4_eq[i]/1000000*(vol_gas[i]/1000) / (R * (temp_eq[i]+273.15)) #mol
    
    mols_headspace.ch4 <- mCH4_headspace[i]/1000000*(vol_gas[i]/1000)/(R * (temp_eq[i]+273.15)) #mol PV / RT = n
    
    #implication: mass, concentration, and partial pressure of CH4 in the original sample (aount in sample and headspace after equilibration minus original mass in the headspace)
    Sample_CH4_mass <- CH4_solution_mass + final_C_headspace_mass.ch4 - mols_headspace.ch4 #mol
    Sample_CH4_conc <- Sample_CH4_mass/(vol_water[i]/1000) #mol/L
    pGHG_orig[i,11] <- Sample_CH4_conc/Kh2.ch4*1000000 #ppmv
    pGHG_orig[i,12] <- Sample_CH4_conc/Kh2.ch4*1000000*Bar.pressure[i]/101.325 # micro-atm
    pGHG_orig[i,13] <- Sample_CH4_conc*1000000 # micro-mol/L
    
    #concentration and total mass in the water sample assuming ideal gas from the pN2O measured at the headspace
    N2O_solution <- mN2O_eq[i]/1000000*Kh.n2o #mol/L
    N2O_solution_mass <- N2O_solution * vol_water[i]/1000 #mol
    
    #mass of N2O in the measured headspace
    final_N_headspace_mass <- mN2O_eq[i]/1000000*(vol_gas[i]/1000) / (R * (temp_eq[i]+273.15)) #mol
    
    mols_headspace.n2o <- mN2O_headspace[i]/1000000*(vol_gas[i]/1000)/(R * (temp_eq[i]+273.15)) #mol PV / RT = n
    
    #implication: mass, concentration, and partial pressure of N2O in the original sample (aount in sample and headspace after equilibration minus original mass in the headspace)
    Sample_N2O_mass <- N2O_solution_mass + final_N_headspace_mass - mols_headspace.n2o #mol
    Sample_N2O_conc <- Sample_N2O_mass/(vol_water[i]/1000) #mol/L
    pGHG_orig[i,14] <- Sample_N2O_conc/Kh2.n2o*1000000 #ppmv
    pGHG_orig[i,15] <- Sample_N2O_conc/Kh2.n2o*1000000*Bar.pressure[i]/101.325 # micro-atm
    pGHG_orig[i,16] <- Sample_N2O_conc*1000000 # micro-mol/L
    pGHG_orig[i,17] <- c_constants[i]
  }
  
  
  return(pGHG_orig) #Output data frame
  
}

g <- Rheadspace(h)

