import numpy as np
import matplotlib.pyplot as plt
from lib.settings import NAME, WELL_VILLAFORTUNA1, WELL_CINZANO, WELL_BENIGNO

#Definition of Model Function

def Rs_line(Rs, Rs_max, alpha, tt, condrock, rext):
    rs_function = lambda tt: (1/2/condrock)*np.log(2 * (np.sqrt(alpha * tt))/rext)
    #print(rs_function(tt))
    Rs.append(rs_function(tt))
    Rs_max.append(max(rs_function(tt)))
    #print(Rs_max)
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

def twat_evaluation(zz, step, time, rext, vel, dhydr, portata, cpw, tground, twat_in, L):
    twat_down_set = []
    twat_horizontal_set = []
    twat_up_set = []
    for ii in range(len(zz)):
        condrock=condrockfunz(zz[ii]) #W/m/K
        rhorock=rhorockfunz(zz[ii]) #kg/m^3
        cprock=cprockfunz(zz[ii]) #j/kg/k
        alpha = condrock/rhorock/cprock

        #The total coefficient is calculated every time because it changes with depth
        #The factor that changes is the resistance of the soil as its parameters changes with depth
        ktot = np.power(rext / condrock * np.log((2 * (np.sqrt(alpha * step))/rext)) + 1 / funzh(vel, dhydr),-1) #%W/m^2/K #far variare t

        #Thermal exchange between two depths at a distance dz
        num1 = zz[ii] * 2 * np.pi * rext * ktot
        den1 = portata * cpw
        aa = (num1 / (2 * den1)) + 1 # st_dz

        twat_down = (num1 / den1 * funztground(zz[ii]) - (num1 / (2 * den1)-1) * twat_in) / aa #C°
        twat_horizontal = (tground[ii] - ((tground[ii] - twat_down)/(np.exp((2 * np.pi * condrock * L)/(portata * cpw * np.log((np.pi * (alpha * time)**0.5)/(2 * rext))))))) #C°
        twat_out = (num1 / den1 * funztground(zz[ii]) - twat_horizontal * aa) / (num1 / (2 * den1)-1)  #C°
        twat_down_set.append(twat_down)
        twat_horizontal_set.append(twat_horizontal)
        twat_up_set.append(twat_out)
        print(f"L: {L}, Depth: {zz[ii]}, twat_in: {twat_in}")
        print(f"twat_down: {twat_down}")
        print(f"twat_horizontal: {twat_horizontal}")
        print(f"twat_out: {twat_out}")
        print()
        #twat_up = (dz * 2 * np.pi * rint * ktot / (portata * cpw) * ((twat_horizontal[jj-1] + twat_horizontal[jj]) / 2) - (dz * 2 * np.pi * rint * ktot / (2 * portata * cpw) - 1) * twat_horizontal) / aa1 #C°
    
    return twat_down_set, twat_horizontal_set, twat_up_set
 

#def twat_lateral_time_evaluation(tground, twat_down, Ktot, L, portata, cpw, Alpha, tt2, riso):
    twat_out2 = []
    for i in tt2:
        twat_out2.append(tground - ((tground - twat_down)/(np.exp((2 * np.pi * Ktot * L)/(portata * cpw * np.log((np.pi * (Alpha * i)**0.5)/(2 * riso))))))) #C°
    
    return twat_out2


#Definition of Graphics

#def plot_1(tt, value):
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


def plot_2(twat_in, twat_down_set, twat_up_set, twat_horizontal_set, tground, zz, step):
    #fig = plt.figure(figsize=(19.20,10.80))
    #plt.subplots(figsize=(15,15))
    fig, ax = plt.subplots()
    ax.set_ylim(zz[-1]+1000, 0)
    ax.set_ylabel('Depth [m]')
    ax.set_xlabel('Temperature [°C]')
    ax.set_title('Evaluation {} years'.format(step/365/86400))
    ax.grid('both', linestyle='--')
    ax.plot(tground, zz,'g', label='T ground')
    for i in range(len(zz)):
        ax.plot([twat_in, twat_down_set[i]], [0, zz[i]], 'b', label='T downward water')
        ax.plot([twat_down_set[i], twat_horizontal_set[i]], [zz[i], zz[i]],'y', label='T lateral water')
        ax.plot([twat_horizontal_set[i], twat_up_set[i]], [zz[i], 0],'r', label='T upward water')
        if i == 0:
            ax.legend(loc='lower left')
    plt.savefig('Output_Temperature_'+str(step/365/86400)+'y_'+str(NAME)+'.png', dpi=1200)
    plt.show()

def plot_3(twat_out2, tt):
    fig, ax = plt.subplots()
    ax.plot( tt/365/86400, twat_out2, 'y', label='T lateral water') 
    ax.set_xlabel('Time [year]')
    ax.set_ylabel('Temperature [°C]')
    ax.set_title('Evaluation at 5402 m depth')
    ax.legend()
    plt.show()
    


#Parameters related to the geological formations
if NAME == 'Villafortuna1':
    condrockfunz = lambda bb: 0.3*(bb<1258)+1.61*(bb>=1258)*(bb<1405)+3.16*(bb>=1405)*(bb<5493)+3.5*(bb>=5493)*(bb<=depth) #W/m/K
    rhorockfunz = lambda bb: 1700*(bb<1258)+1890*(bb>=1258)*(bb<1405)+2359*(bb>=1405)*(bb<5493)+2480*(bb>=5493)*(bb<=depth) #kg/m^3
    cprockfunz = lambda bb: 800*(bb<1258)+1693*(bb>=1258)*(bb<1405)+821.111*(bb>=1405)*(bb<5493)+810.484*(bb>=5493)*(bb<=depth) #j/kg/k
    WELL = WELL_VILLAFORTUNA1
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
#riso = WELL['ISO Radius']
#rint = WELL['INT Radius']

#Operating time
T = 31536000
tt = np.arange(WELL['START']*T, WELL['STOP']*T, WELL['STEP']*T) #s

#discretization of depth
dz = WELL['DZ'] #m
zz = np.arange(0,depth,dz) #m

####################### Exchanger analysis

#Length of the lateral portion
L = WELL['L']
#Flow rate entered
portata = WELL['PORTATA']
#Set the temperature of the water entering the well
twat_in = WELL['T WATER IMMESSA']
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
#print(value,Rs[130])

#Graphic creation
#plot_1(tt, Rs[value[0]])


####################### Exchanger analysis
#Coeff. of total heat exchange including soil and convective part in the hydraulic ring 
#Diameter for the calculation of the coeff. dimensionless of the ring 
dhydr = (4 * np.pi * np.power(rext, 2) ) / (2 * np.pi * rext ) #m 
areac = np.pi * np.power(rext, 2) #%m^2
vel = portata / areac / rhow #m/s
time = L / vel

#Coefficient of convective heat transfer
hext = funzh(vel, dhydr) #W/m^2/K
#Calculate the ground temperature with the depth
tground = funztground(zz) #C°

#Temperature variation in lateral part
#for step in t:
    #twat_down_set = twat_evaluation(zz, step, rext, vel, dhydr, portata, cpw, tground, twat_in, L)
#twat_down = twat_down_set[1]
#twat_out2 = twat_lateral_time_evaluation(tground, twat_down, Kr, L, portata, cpw, Alpha, tt, riso)
#plot_3(twat_out2, tt)

#Loop for ring temp. profile
for step in t: 
    twat_down_set, twat_horizontal_set, twat_up_set = twat_evaluation(zz, step, time, rext, vel, dhydr, portata, cpw, tground, twat_in, L)
    plot_2(twat_in, twat_down_set, twat_up_set, twat_horizontal_set, tground, zz, step)


