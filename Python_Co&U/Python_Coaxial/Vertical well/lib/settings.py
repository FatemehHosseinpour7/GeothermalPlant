#NAME = 'Temparossa'
#NAME = 'Gela'
#NAME = 'Villafortuna1'
#NAME = 'Cinzano'
NAME = 'Benigno'
#NAME = 'Trecate4'
# NAME = 'Casteggio'
# NAME = 'Cammarata'

#1
WELL_VILLAFORTUNA1 = {
    'DEPTH': 6202,
    'DZ': 10,
    
    # Coaxial
    'EXT Radius': 15e-2/2, #m
    'ISO Radius': 13.9e-2/2, #m
    'INT Radius': 7.7e-2/2, #m
    
    # Utube
    'WELL Radius' : 150e-3/2, #m
    'TUBE Radius' : 43e-3/2, #m
    'W Distance' : 80e-3, #m 
    
    # Time Analysis
    'START': 0.5, #year
    'STEP': 0.5, #year
    'STOP': 10, #year
    
    # Flow
    'PORTATA': 3, #kg/s
    # Temperature of water entering the well
    'T WATER IMMESSA': 15, #C°
    # Steps relating to the selected years of analysis
    'STEP T': [10], #year
    # Soil temperature profile in relation to depth
    'A': 0.025, 
    'B': 15,
    # Parameters of the heat transfer fluid (water) 100°C and 2 bar
    'CPW' : 4215.45, #j/kg/k
    'COND WATER' : 0.679146, #W/m/K
    'RHOW' : 958.395, #kg/m^3
    'VISC WATER' : 2.8177e-4, #Pa s
    'COND ISO' : 0.025, #W/m/K
    
    'Utube CPW': 4186, #j/kg/k
    'Utube COND WATER' : 0.679146, #W/m/K
    'Utube RHOW' : 958.395, #kg/m^3
    'Utube VISC WATER' : 2.8177e-4, #Pa s
    'Utube COND ISO' : 0.025, #W/m/K
    'COND GROUT' : 2, #W/m/K
    
    'VISC WATER Utube' : 8.4e-4 #Pa s

}

WELL_CINZANO = {
    'DEPTH': 1378.65,
    'DZ': 10,

    # Lateral
    'L' : 1000, #m
    
    # Coaxial
    'EXT Radius': 15e-2/2, #m
    'ISO Radius': 13.9e-2/2, #m
    'INT Radius': 7.7e-2/2, #m
    
    
    # Time Analysis
    'START': 0.5, #year
    'STEP': 0.5, #year
    'STOP': 10, #year
    
    # Flow
    'PORTATA': 3, #kg/s
    # Temperature of water entering the well
    'T WATER IMMESSA': 15, #C°
    # Steps relating to the selected years of analysis
    'STEP T': [8, 10], #year
    # Soil temperature profile in relation to depth
    'A': 0.015, 
    'B': 15,
    # Parameters of the heat transfer fluid (water) 100°C and 2 bar
    'CPW' : 4215.45, #j/kg/k
    'COND WATER' : 0.679146, #W/m/K
    'RHOW' : 958.395, #kg/m^3
    'VISC WATER' : 2.8177e-4, #Pa s
    'COND ISO' : 0.025, #W/m/K
      
}

#3
WELL_BENIGNO = {
    'DEPTH': 2701,
    'DZ': 10,

    # Lateral
    'L' : 1000, #m
    
    # Coaxial
    'EXT Radius': 15e-2/2, #m
    'ISO Radius': 13.9e-2/2, #m
    'INT Radius': 7.7e-2/2, #m
    
    
    # Time Analysis
    'START': 0.5, #year
    'STEP': 0.5, #year
    'STOP': 10, #year
    
    # Flow
    'PORTATA': 3, #kg/s
    # Temperature of water entering the well
    'T WATER IMMESSA': 15, #C°
    # Steps relating to the selected years of analysis
    'STEP T': [10], #year
    # Soil temperature profile in relation to depth
    'A': 0.015, 
    'B': 15,
    # Parameters of the heat transfer fluid (water) 100°C and 2 bar
    'CPW' : 4215.45, #j/kg/k
    'COND WATER' : 0.679146, #W/m/K
    'RHOW' : 958.395, #kg/m^3
    'VISC WATER' : 2.8177e-4, #Pa s
    'COND ISO' : 0.025, #W/m/K
      
}

#2
WELL_TRECATE4 = {
    'DEPTH': 6282,
    'DZ': 10,
    
    # Coaxial
    'EXT Radius': 15e-2/2, #m
    'ISO Radius': 13.9e-2/2, #m
    'INT Radius': 7.7e-2/2, #m
    
    # Utube
    
    'WELL Radius' : 150e-3/2, #m
    'TUBE Radius' : 43e-3/2, #m
    'W Distance' : 80e-3, #m 
    
    # Time Analysis
    'START': 0.5, #year
    'STEP': 0.5, #year
    'STOP': 10, #year
    
    # Flow
    'PORTATA': 3, #kg/s
    # Temperature of water entering the well
    'T WATER IMMESSA': 50, #C°
    # Steps relating to the selected years of analysis
    'STEP T': [3, 5], #year
    # Soil temperature profile in relation to depth
    'A': 0.025, 
    'B': 15,
    
    # Parameters of the heat transfer fluid (water) 100°C and 2 bar
    'CPW' : 4215.45, #j/kg/k
    'COND WATER' : 0.679146, #W/m/K
    'RHOW' : 958.395, #kg/m^3
    'VISC WATER' : 2.8177e-4, #Pa s
    'COND ISO' : 0.025, #W/m/K
        
    'Utube CPW': 4186, #j/kg/k
    'Utube COND WATER' : 0.679146, #W/m/K
    'Utube RHOW' : 958.395, #kg/m^3
    'Utube VISC WATER' : 2.8177e-4, #Pa s
    'Utube COND ISO' : 0.025, #W/m/K
    'COND GROUT' : 2, #W/m/K
    
    'VISC WATER Utube' : 8.4e-4 #Pa s
    
}

#3
WELL_TEMPAROSSA = {
    'DEPTH': 5042,
    'DZ': 10,
    
    # Coaxial
    'EXT Radius': 15e-2/2, #m
    'ISO Radius': 13.9e-2/2, #m
    'INT Radius': 7.7e-2/2, #m

    # Utube 

    'WELL Radius' : 150e-3/2, #m
    'TUBE Radius' : 43e-3/2, #m
    'W Distance' : 80e-3, #m 
    
    # Time Analysis
    'START': 0.5, #year
    'STEP': 0.5, #year
    'STOP': 10, #year
    
    # Flow
    'PORTATA': 3, #kg/s
    # Temperature of water entering the well
    'T WATER IMMESSA': 50, #C°
    # Steps relating to the selected years of analysis
    'STEP T': [3,5], #year
    
    # Soil temperature profile in relation to depth
    'A': 0.015, 
    'B': 15,
    
    # Parameters of the heat transfer fluid (water) 100°C and 2 bar
    'CPW' : 4215.45, #j/kg/k
    'COND WATER' : 0.679146, #W/m/K
    'RHOW' : 958.395, #kg/m^3
    'VISC WATER' : 2.8177e-4, #Pa s
    'COND ISO' : 0.025, #W/m/K
    
    'Utube CPW': 4186, #j/kg/k
    'Utube COND WATER' : 0.679146, #W/m/K
    'Utube RHOW' : 958.395, #kg/m^3
    'Utube VISC WATER' : 2.8177e-4, #Pa s
    'Utube COND ISO' : 0.025, #W/m/K
    'COND GROUT' : 2, #W/m/K
    

    'VISC WATER Utube' : 8.4e-4 #Pa s

}

#4
WELL_GELA = {
    'DEPTH': 3446,
    'DZ': 10,

    # Coaxial
    'EXT Radius': 15e-2/2, #m
    'ISO Radius': 13.9e-2/2, #m
    'INT Radius': 7.7e-2/2, #m

    # Utube

    'WELL Radius' : 150e-3/2, #m
    'TUBE Radius' : 43e-3/2, #m
    'W Distance' : 80e-3, #m
    
    # Time Analysis
    'START': 0.5, #year
    'STEP': 0.5, #year
    'STOP': 10, #year
    
    # Flow
    'PORTATA': 3, #kg/s
    # Temperature of water entering the well
    'T WATER IMMESSA': 50, #C°
    # Steps relating to the selected years of analysis
    'STEP T': [3, 5], #year
    # Soil temperature profile in relation to depth
    'A': 0.025, 
    'B': 15,
    
    # Parameters of the heat transfer fluid (water) 100°C and 2 bar
    'CPW' : 4215.45, #j/kg/k
    'COND WATER' : 0.679146, #W/m/K
    'RHOW' : 958.395, #kg/m^3
    'VISC WATER' : 2.8177e-4, #Pa s
    'COND ISO' : 0.025, #W/m/K
   
    'Utube CPW': 4186, #j/kg/k
    'Utube COND WATER' : 0.679146, #W/m/K
    'Utube RHOW' : 958.395, #kg/m^3
    'Utube VISC WATER' : 2.8177e-4, #Pa s
    'Utube COND ISO' : 0.025, #W/m/K
    'COND GROUT' : 2, #W/m/K

    'VISC WATER Utube' : 8.4e-4 #Pa s

}

#5
WELL_CASTEGGIO = {
    'DEPTH': 2501,
    'DZ': 10,

    # Coaxial
    'EXT Radius': 15e-2/2, #m
    'ISO Radius': 13.9e-2/2, #m
    'INT Radius': 7.7e-2/2, #m

     # Utube 

    'WELL Radius' : 150e-3/2, #m
    'TUBE Radius' : 43e-3/2, #m
    'W Distance' : 80e-3, #m
    
    # Time Analysis
    'START': 0.5, #year
    'STEP': 0.5, #year
    'STOP': 10, #year
    # Flow
    'PORTATA': 3, #kg/s
    # Temperature of water entering the well
    'T WATER IMMESSA': 50, #C°
    # Steps relating to the selected years of analysis
    'STEP T': [3, 5], #year
    
    # Soil temperature profile in relation to depth
    'A': 0.020, 
    'B': 15,
    # Parameters of the heat transfer fluid (water) 100°C and 2 bar
    'CPW' : 4215.45, #j/kg/k
    'COND WATER' : 0.679146, #W/m/K
    'RHOW' : 958.395, #kg/m^3
    'VISC WATER' : 2.8177e-4, #Pa s
    'COND ISO' : 0.025, #W/m/K
   
    'Utube CPW': 4186, #j/kg/k
    'Utube COND WATER' : 0.679146, #W/m/K
    'Utube RHOW' : 958.395, #kg/m^3
    'Utube VISC WATER' : 2.8177e-4, #Pa s
    'Utube COND ISO' : 0.025, #W/m/K
    'COND GROUT' : 2, #W/m/K

    'VISC WATER Utube' : 8.4e-4 #Pa s

}

#6
WELL_CAMMARATA = {
    'DEPTH': 3567.1,
    'DZ': 10,

    # Coaxial
    'EXT Radius': 15e-2/2, #m
    'ISO Radius': 13.9e-2/2, #m
    'INT Radius': 7.7e-2/2, #m

    # Utube 

    'WELL Radius' : 150e-3/2, #m
    'TUBE Radius' : 43e-3/2, #m
    'W Distance' : 80e-3, #m
    
    # Time Analysis
    'START': 0.5, #year
    'STEP': 0.5, #year
    'STOP': 10, #year
    # Flow
    'PORTATA': 3, #kg/s
    # Temperature of water entering the well
    'T WATER IMMESSA': 50, #C°
    # Steps relating to the selected years of analysis
    'STEP T': [3, 5], #year
    
    # Soil temperature profile in relation to depth
    'A': 0.030 , 
    'B': 15,
    
    # Parameters of the heat transfer fluid (water) 100°C and 2 bar
    'CPW' : 4215.45, #j/kg/k
    'COND WATER' : 0.679146, #W/m/K
    'RHOW' : 958.395, #kg/m^3
    'VISC WATER' : 2.8177e-4, #Pa s
    'COND ISO' : 0.025, #W/m/K
   
    'Utube CPW': 4186, #j/kg/k
    'Utube COND WATER' : 0.679146, #W/m/K
    'Utube RHOW' : 958.395, #kg/m^3
    'Utube VISC WATER' : 2.8177e-4, #Pa s
    'Utube COND ISO' : 0.025, #W/m/K
    'COND GROUT' : 2, #W/m/K

    'VISC WATER Utube' : 8.4e-4 #Pa s

}

SENSIBILITY = {
    # CPROCK
    'CPROCK_START': 650, #j/kg/k
    'CPROCK_STEP': 50, #j/kg/k
    'CPROCK_STOP': 1050, #j/kg/k
    'CP_CONDROCK': 2.5, #W/m/K
    'CP_RHOROCK': 2125, #kg/m^3
    # CONDROCK
    'CONDROCK_START': 1.5, #W/m/K
    'CONDROCK_STEP': 0.25, #W/m/K
    'CONDROCK_STOP': 3.5, #W/m/K
    'COND_CPROCK': 820, #j/kg/k
    'COND_RHOROCK': 2125, #kg/m^3
    # RHOROCK
    'RHOROCK_START': 1000, #kg/m^3
    'RHOROCK_STEP': 200, #kg/m^3
    'RHOROCK_STOP': 2500, #kg/m^3
    'RHO_CPROCK': 820, #j/kg/k
    'RHO_CONDROCK': 2.5, #W/m/K
        #ALPHA SENSIBILITY TRECATE4
    'A_CPROCK_START': 1775, #j/kg/k
    'A_CPROCK_STOP': 3500, #j/kg/k
    'A_CONDROCK_START': 2.5, #j/kg/k
    'A_CONDROCK_STOP': 5, #j/kg/k
    'A_RHOROCK_START': 2125, #j/kg/k
    'A_RHOROCK_STOP': 4300, #j/kg/k
    'VALUE_STEP': 10
    }
