#NAME = 'Villafortuna1'
#NAME = 'Cinzano'
NAME = 'Benigno'

#1
WELL_VILLAFORTUNA1 = {
    'DEPTH': 6202,
    'DZ': 50,

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
    'T WATER IMMESSA': 154.3544 , #C°
    # Steps relating to the selected years of analysis
    'STEP T': [3], #year
    # Soil temperature profile in relation to depth
    'A': 0.025, 
    'B': 15,
    # Parameters of the heat transfer fluid (water) 100°C and 2 bar
    'CPW' : 4215.45, #j/kg/k
    'COND WATER' : 0.679146, #W/m/K
    'RHOW' : 958.395, #kg/m^3
    'VISC WATER' : 2.8177e-4, #Pa s
    'COND ISO' : 0.025, #W/m/K
      
}

#2
WELL_CINZANO = {
    'DEPTH': 1378.65,
    'DZ': 50,

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
    'T WATER IMMESSA': 31.7446, #C°
    # Steps relating to the selected years of analysis
    'STEP T': [3], #year
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
    'DEPTH': 2700,
    'DZ': 50,

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
    'T WATER IMMESSA': 55.3163, #C°
    # Steps relating to the selected years of analysis
    'STEP T': [3], #year
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
