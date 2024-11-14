import numpy as np
import matplotlib.pyplot as plt
from lib.settings import NAME, WELL_TEMPAROSSA, WELL_VILLAFORTUNA1, WELL_TRECATE4, WELL_GELA, WELL_CASTEGGIO, WELL_CAMMARATA, WELL_CINZANO, WELL_BENIGNO

#Definition of Model Function

def Rs_line(Rs, Rs_max, alpha, tt, condrock, rext):
    rs_function = lambda tt: (1/2/condrock)*np.log(2 * (np.sqrt(alpha * tt))/rext)
    print(rs_function(tt))
    Rs.append(rs_function(tt))
    Rs_max.append(max(rs_function(tt)))
    print(Rs_max)
    return Rs, Rs_max


def Rs_matrix(zz, tt, rext):
    Rs, Rs_max = [], []
    #Definition of alpha and Rs
    for ii in range(1,len(zz)):
        condrock=condrockfunz(zz[ii]) #W/m/K
        rhorock=rhorockfunz(zz[ii]) #kg/m^3
        cprock=cprockfunz(zz[ii]) #j/kg/k
        alpha = condrock/rhorock/cprock
        Rs, Rs_max = Rs_line(Rs, Rs_max, alpha, tt, condrock, rext)
    return Rs, Rs_max


def find_max(Rs, Rs_max):
    value = [0,0]
    for x in range(len(Rs)):
        i = np.where(Rs[x] == max(Rs_max))
        if len(i[0]) > 0:
            if value[0] < x and value[1] <= i[0][0]:
                value = [x, i[0][0]]
    return value


def twat_down_evaluation(zz, step, rext, vel, dhydr, portata, cpw, tground, twat_value):
    twat_down = []
    for ii in range(len(zz)):
        condrock=condrockfunz(zz[ii]) #W/m/K
        rhorock=rhorockfunz(zz[ii]) #kg/m^3
        cprock=cprockfunz(zz[ii]) #j/kg/k
        alpha = condrock/rhorock/cprock
        
        #The total coefficient is calculated every time because it changes with depth
        #The factor that changes is the resistance of the soil as its parameters changes with depth
        ktot = np.power(rext / condrock * np.log((2 * (np.sqrt(alpha * step))/rext)) + 1 / funzh(vel, dhydr),-1) #%W/m^2/K #far variare t
        
        #Thermal exchange between two depths at a distance dz
        num1 = dz * 2 * np.pi * rext * ktot
        den1 = portata * cpw
        aa = (num1 / (2 * den1)) + 1 # st_dz
        
        #Water temperature after heat exchange (temperature vector)
        twat_value = (num1 / den1 * funztground(zz[ii] - dz / 2) - (num1 / (2 * den1)-1) * twat_value) / aa #C°
        twat_down.append(twat_value)
    return twat_down


def twat_up_evaluation(zz, twat, step, rint, riso, condiso, portata, cpw):
    
    #I invert the vectors because now we start from the bottom and go up and so when We want the first value up we have to refer to the last one
    #Falling temperature    
    zz2 = np.flip(zz) #m
    twatdalbasso = np.flip(twat)
    vel2 = portata / (np.pi * np.power(rint,2)) / rhow #m/s
    hint = funzh(vel2, 2*rint) #W/m^2/K
    hexttubo = funzh(vel, 2*riso) #W/m^2/K
    
    #Heat transfer coefficient including the two convective resistances and the conductive one due to the insulator
    ktotinner = np.power(rint / riso * 1/ hexttubo + rint * np.log(riso / rint) * 1 / condiso + 1 / hint,-1) #W/m^2/K
    aa1 = (dz * 2 * np.pi * rint * ktotinner / (2 * portata * cpw)) + 1
    
    #Loop for the water that returns to the surface, the first value equal to the last of the calculation for the water that reaches the bottom
    twat_up_value = twat[-1] #C°
    twat_up = []
    for jj in range(0, len(zz2)):
        twat_up_value = (dz * 2 * np.pi * rint * ktotinner / (portata * cpw) * ((twatdalbasso[jj-1] + twatdalbasso[jj]) / 2) - (dz * 2 * np.pi * rint * ktotinner / (2 * portata * cpw) - 1) * twat_up_value) / aa1 #C°
        twat_up.append(twat_up_value)
    return twat_up

#Definition of Graphics

def plot_1(tt, value):
    f = plt.figure()
    f.set_figwidth(20)
    f.set_figheight(10)
    plt.plot(tt/3600/24/365, value, 'g')
    plt.xlabel('Time [years]')
    plt.ylabel('Thermal Resistance - Rs [m^2K/W]')
    plt.grid('both', linestyle='--')
    plt.xticks([1,2,3,4,5,6,7,8,9,10])
    
    plt.savefig('ThermalResistance_'+str(NAME)+'.png', dpi=1200)
    plt.show()


def plot_2(twat, twat2, tground, zz, step):
    # fig = plt.figure(figsize=(19.20,10.80))
    # plt.subplots(figsize=(15,15))
    fig, ax = plt.subplots()
    ax.plot(twat, zz, 'b', label='T downward water')
    ax.plot(twat2, np.flip(zz), 'r', label='T upward water')
    ax.plot(tground, zz, 'g', label='T ground')
    ax.set_ylim(zz[-1]+1000, 0)  # decreasing time
    ax.set_xlabel('Temperature [°C]')
    ax.set_ylabel('Depth [m]')
    ax.set_title('Evaluation {} years'.format(step/365/86400))
    ax.grid('both', linestyle='--')
    ax.legend()
    plt.savefig('Coaxial_profile_'+str(step/365/86400)+'y_'+str(NAME)+'.png', dpi=1200)
    plt.show()
    

#Parameters related to the geological formations
if NAME == 'Temparossa':
    condrockfunz = lambda bb: 3.0*(bb<2912)+2.34*(bb>=2912)*(bb<=depth) #W/m/K
    rhorockfunz = lambda bb: 2330*(bb<2912)+1917*(bb>=2912)*(bb<=depth) #kg/m^3
    cprockfunz = lambda bb: 808.6*(bb<2912)+829.4*(bb>=2912)*(bb<=depth) #j/kg/k
    WELL = WELL_TEMPAROSSA
elif NAME == 'Villafortuna1':
    condrockfunz = lambda bb: 0.3*(bb<1258)+1.61*(bb>=1258)*(bb<1405)+3.16*(bb>=1405)*(bb<5493)+3.5*(bb>=5493)*(bb<=depth) #W/m/K
    rhorockfunz = lambda bb: 1700*(bb<1258)+1890*(bb>=1258)*(bb<1405)+2359*(bb>=1405)*(bb<5493)+2480*(bb>=5493)*(bb<=depth) #kg/m^3
    cprockfunz = lambda bb: 800*(bb<1258)+1693*(bb>=1258)*(bb<1405)+821.111*(bb>=1405)*(bb<5493)+810.484*(bb>=5493)*(bb<=depth) #j/kg/k
    WELL = WELL_VILLAFORTUNA1
elif NAME == 'Trecate4':
    condrockfunz = lambda bb: 1.61*(bb<1632)+2.17*(bb>=1632)*(bb<5106)+3.50*(bb>=5106)*(bb<6189)+3.00*(bb>=6189)*(bb<=depth) #W/m/K
    rhorockfunz = lambda bb: 1890*(bb<1632)+1801*(bb>=1632)*(bb<5106)+2480*(bb>=5106)*(bb<6189)+2330*(bb>=6189)*(bb<=depth) #kg/m^3
    cprockfunz = lambda bb: 1696*(bb<1632)+830.09*(bb>=1632)*(bb<5106)+810.48*(bb>=5106)*(bb<6189)+821.11*(bb>=6189)*(bb<=depth) #j/kg/k
    WELL = WELL_TRECATE4
elif NAME == 'Gela':
    condrockfunz = lambda bb: 3.16*(bb<2117)+2.17*(bb>=2117)*(bb<2556)+3.12*(bb>=2556)*(bb<2860)+2.17*(bb>=2860)*(bb<=depth) #W/m/K
    rhorockfunz = lambda bb: 2359*(bb<2117)+1801*(bb>=2117)*(bb<2556)+2480*(bb>=2556)*(bb<2860)+1801*(bb>=2860)*(bb<=depth) #kg/m^3
    cprockfunz = lambda bb: 821.11*(bb<2117)+830.09*(bb>=2117)*(bb<2556)+810.48*(bb>=2556)*(bb<2860)+830.09*(bb>=2860)*(bb<=depth) #j/kg/k
    WELL = WELL_GELA
elif NAME == 'Casteggio':
    condrockfunz = lambda bb: 2.45*(bb<660)+3.28*(bb>=660)*(bb<2212)+2.77*(bb>=2212)*(bb<=depth) #W/m/K
    rhorockfunz = lambda bb: 1757*(bb<660)+2161*(bb>=660)*(bb<2212)+1787*(bb>=2212)*(bb<=depth) #kg/m^3
    cprockfunz = lambda bb: 1459*(bb<660)+1125*(bb>=660)*(bb<2212)+733.072*(bb>=2212)*(bb<=depth) #j/kg/k
    WELL = WELL_CASTEGGIO
elif NAME == 'Cammarata':
    condrockfunz = lambda bb: 3.16*(bb<2362.1)+2.77*(bb>=2362.1)*(bb<2567.1)+3.6*(bb>=2567.1)*(bb<=depth) #W/m/K
    rhorockfunz = lambda bb: 2359*(bb<2362.1)+2278*(bb>=2362.1)*(bb<2567.1)+2520*(bb>=2567.1)*(bb<=depth) #kg/m^3
    cprockfunz = lambda bb: 821.11*(bb<2362.1)+793.68*(bb>=2362.1)*(bb<2567.1)+807.94*(bb>=2567.1)*(bb<=depth) #j/kg/k
    WELL = WELL_CAMMARATA
elif NAME == 'Cinzano':
    condrockfunz = lambda bb: 4.44*(bb<100)+5.03*(bb>=100)*(bb<728)+4.44*(bb>=728)*(bb<1280)+2.77*(bb>=1280)*(bb<=depth) #W/m/K
    rhorockfunz = lambda bb: 1785*(bb<100)+1916*(bb>=100)*(bb<728)+1785*(bb>=728)*(bb<1280)+2278*(bb>=1280)*(bb<=depth) #kg/m^3
    cprockfunz = lambda bb: 1175*(bb<100)+1270*(bb>=100)*(bb<728)+1175*(bb>=728)*(bb<1280)+1808*(bb>=1280)*(bb<=depth) #j/kg/k
    WELL = WELL_CINZANO
elif NAME == 'Benigno':
    condrockfunz = lambda bb: 4.44*(bb<260)+2.45*(bb>=260)*(bb<1100)+2.77*(bb>=1100)*(bb<2300)+2.77*(bb>=2300)*(bb<=depth) #W/m/K
    rhorockfunz = lambda bb: 1785*(bb<260)+1757*(bb>=260)*(bb<1100)+2278*(bb>=1100)*(bb<2300)+2278*(bb>=2300)*(bb<=depth) #kg/m^3
    cprockfunz = lambda bb: 1175*(bb<260)+1459*(bb>=260)*(bb<1100)+1808*(bb>=1100)*(bb<2300)+1808*(bb>=2300)*(bb<=depth) #j/kg/k
    WELL = WELL_BENIGNO

#Geothermal gradient function
funztground = lambda z : a*z + b #C°

#Function of the calculation of the coeff. convection in the pipes. I find Nusselt from Dittus-Boelter see
funzh = lambda v, drif : np.power((rhow * v * drif / viscw), 0.8) * np.power(pr, 0.4) * condw * 0.023 / (drif) #W/m^2/K

#Depth
depth = WELL['DEPTH']
#Tubes parameters
rext = WELL['EXT Radius']
riso = WELL['ISO Radius']
rint = WELL['INT Radius']

#Operating time
T = 31536000
tt = np.arange(WELL['START']*T, WELL['STOP']*T, WELL['STEP']*T) #s

#discretization of depth
dz = WELL['DZ'] #m
zz = np.arange(0,depth,dz) #m

####################### Exchanger analysis
#Flow rate entered
portata = WELL['PORTATA']
#Set the temperature of the water entering the well
twat_value = WELL['T WATER IMMESSA']
#Steps relating to the years of analysis
t = []
for x in WELL['STEP T']:
    t.append(x*T) #S is a parameter that can vary
    
#Soil temperature profile in relation to depth
a, b = WELL['A'], WELL['B']
#Water parameters at 100C� and 2 bar
cpw = WELL['CPW'] #j/kg/k
condw = WELL['COND WATER'] #W/m/K
rhow = WELL['RHOW'] #kg/m^3
viscw = WELL['VISC WATER'] #Pa s
condiso = WELL['COND ISO'] #W/m/K

pr = cpw * viscw / condw

########### Block 1 ##########################
#Definition of the Rs function with depth over time 
Rs, Rs_max = Rs_matrix(zz, tt, rext)

#Find the depth of the maximum Rs
value = find_max(Rs, Rs_max)


#Graphic creation
plot_1(tt, Rs[value[0]])


####################### Exchanger analysis
#Coeff. of total heat exchange including soil and convective part in the hydraulic ring 
#Diameter for the calculation of the coeff. dimensionless of the ring 
dhydr = (4 * np.pi * (np.power(rext, 2) - np.power(riso, 2))) / (2 * np.pi * (rext + riso)) #m 
areac = np.pi * (np.power(rext, 2) - np.power(riso, 2)) #%m^2
vel = portata / areac / rhow #m/s

#Coefficient of convective heat transfer
hext = funzh(vel, dhydr) #W/m^2/K
#Calculate the ground temperature with the depth
tground1 = funztground(zz) #C°

#Loop for ring temp. profile
for step in t: 
    twat_down = twat_down_evaluation(zz, step, rext, vel, dhydr, portata, cpw, tground1, twat_value)
    twat_up = twat_up_evaluation(zz, twat_down, step, rint, riso, condiso, portata, cpw)
    plot_2(twat_down, twat_up, tground1, zz, step)

print(twat_up)
