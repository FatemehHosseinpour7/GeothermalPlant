import numpy as np
import matplotlib.pyplot as plt
from lib.settings import NAME, WELL_VILLAFORTUNA1, WELL_TRECATE4

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

def twat_lateral_depth_evaluation(Tr, twat_down, L, portata, cpw, step, riso):
    twat_out = []
    for i in range(len(twat_down)):
        condrock=condrockfunz(zz[i]) #W/m/K
        rhorock=rhorockfunz(zz[i]) #kg/m^3
        cprock=cprockfunz(zz[i]) #j/kg/k
        alpha = condrock/rhorock/cprock
        
        twat_out.append(Tr - ((Tr - twat_down[i])/(np.exp((2 * np.pi * condrock * L)/(portata * cpw * np.log((np.pi * (alpha * step)**0.5)/(2 * riso))))))) #C°
    
    return twat_out 

def twat_lateral_time_evaluation(Tr, twat_down, Kr, L, portata, cpw, Alpha, tt2, riso):
    twat_out2 = []
    for i in tt2:
        twat_out2.append(Tr - ((Tr - twat_down)/(np.exp((2 * np.pi * Kr * L)/(portata * cpw * np.log((np.pi * (Alpha * i)**0.5)/(2 * riso))))))) #C°
    
    return twat_out2

def twat_up_evaluation(zz, twat_out, step, rint, riso, condiso, portata, cpw):
    
    #I invert the vectors because now we start from the bottom and go up and so when We want the first value up we have to refer to the last one
    #Falling temperature    
    zz2 = np.flip(zz) #m
    twatdalbasso = np.flip(twat_out)
    vel2 = portata / (np.pi * np.power(rint,2)) / rhow #m/s
    hint = funzh(vel2, 2*rint) #W/m^2/K
    hexttubo = funzh(vel, 2*riso) #W/m^2/K
    
    #Heat transfer coefficient including the two convective resistances and the conductive one due to the insulator
    ktotinner = np.power(rint / riso * 1/ hexttubo + rint * np.log(riso / rint) * 1 / condiso + 1 / hint,-1) #W/m^2/K
    aa1 = (dz * 2 * np.pi * rint * ktotinner / (2 * portata * cpw)) + 1
    
    #Loop for the water that returns to the surface, the first value equal to the last of the calculation for the water that reaches the bottom
    twat_up_value = twat_out[-1] #C°
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


def plot_2(twat_down, twat_up, tground, zz, step):
    # fig = plt.figure(figsize=(19.20,10.80))
    # plt.subplots(figsize=(15,15))
    fig, ax = plt.subplots()
    ax.plot(twat_down, zz, 'b', label='T downward water')
    ax.plot(twat_out, zz, 'y', label='T lateral water')
    ax.plot(twat_up, np.flip(zz), 'r', label='T upward water')
    ax.plot(tground, zz, 'g', label='T ground')
    ax.set_ylim(zz[-1]+1000, 0)  # decreasing time
    ax.set_xlabel('Temperature [°C]')
    ax.set_ylabel('Depth [m]')
    ax.set_title('Evaluation {} years'.format(step/365/86400))
    ax.grid('both', linestyle='--')
    ax.legend()
    plt.savefig('Coaxial_profile_'+str(step/365/86400)+'y_'+str(NAME)+'.png', dpi=1200)
    plt.show()

def plot_3(twat_out2, tt):
    fig, ax = plt.subplots()
    ax.plot( tt/365/86400, twat_out2, 'y', label='T lateral water') 
    ax.set_xlabel('Time [year]')
    ax.set_ylabel('Temperature [°C]')
    ax.set_title('Evaluation at {} m depth'.format(depth))
    ax.legend()
    plt.savefig("LateralTemperatureVariation"+NAME+".png")
    plt.show()
    

#Parameters related to the geological formations
if  NAME == 'Villafortuna1':
    condrockfunz = lambda bb: 0.3*(bb<1258)+1.61*(bb>=1258)*(bb<1405)+3.16*(bb>=1405)*(bb<5493)+3.5*(bb>=5493)*(bb<=depth) #W/m/K
    rhorockfunz = lambda bb: 1700*(bb<1258)+1890*(bb>=1258)*(bb<1405)+2359*(bb>=1405)*(bb<5493)+2480*(bb>=5493)*(bb<=depth) #kg/m^3
    cprockfunz = lambda bb: 800*(bb<1258)+1693*(bb>=1258)*(bb<1405)+821.111*(bb>=1405)*(bb<5493)+810.484*(bb>=5493)*(bb<=depth) #j/kg/k
    WELL = WELL_VILLAFORTUNA1
elif NAME == 'Trecate4':
    condrockfunz = lambda bb: 1.61*(bb<1632)+2.17*(bb>=1632)*(bb<5106)+3.50*(bb>=5106)*(bb<6189)+3.00*(bb>=6189)*(bb<=depth) #W/m/K
    rhorockfunz = lambda bb: 1890*(bb<1632)+1801*(bb>=1632)*(bb<5106)+2480*(bb>=5106)*(bb<6189)+2330*(bb>=6189)*(bb<=depth) #kg/m^3
    cprockfunz = lambda bb: 1696*(bb<1632)+830.09*(bb>=1632)*(bb<5106)+810.48*(bb>=5106)*(bb<6189)+821.11*(bb>=6189)*(bb<=depth) #j/kg/k
    WELL = WELL_TRECATE4


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
#Formation temperature at lateral part
Tr = WELL['Tr']
#Rock conductivity at lateral part
Kr = WELL['Kr']
#Length of the lateral portion
L = WELL['L']
#Thermal diffiusivity at lateral portion
Alpha = WELL['Alpha']
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
condiso = WELL['VISC WATER'] #W/m/K

pr = cpw * viscw / condw

########### Block 1 ##########################
#Definition of the Rs function with depth over time 
Rs, Rs_max = Rs_matrix(zz, tt, rext)

#Find the depth of the maximum Rs
value = find_max(Rs, Rs_max)
print(value,Rs[130])

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

#Temperature variation in lateral part
for step in t:
    twat_down = twat_down_evaluation(zz, step, rext, vel, dhydr, portata, cpw, tground1, twat_value)
twat_down = twat_down[1]
twat_out2 = twat_lateral_time_evaluation(Tr, twat_down, Kr, L, portata, cpw, Alpha, tt, riso)
plot_3(twat_out2, tt)

#Loop for ring temp. profile
for step in t: 
    twat_down = twat_down_evaluation(zz, step, rext, vel, dhydr, portata, cpw, tground1, twat_value)
    twat_out = twat_lateral_depth_evaluation(Tr, twat_down, L, portata, cpw, step, riso)
    twat_up = twat_up_evaluation(zz, twat_out, step, rint, riso, condiso, portata, cpw)
    plot_2(twat_down, twat_up, tground1, zz, step)


