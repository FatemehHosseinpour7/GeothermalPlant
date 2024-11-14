import numpy as np
import matplotlib.pyplot as plt
from lib.settings import NAME, WELL_VILLAFORTUNA1, WELL_CINZANO, WELL_BENIGNO, SENSIBILITY


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


valuestep = SENSIBILITY['VALUE_STEP']
#print(valuestep)
#print(SENSIBILITY['A_CONDROCK_START'])
########### Sensibilità ##########################
cprockvect = SENSIBILITY['CPROCK_START']
condrockvect = np.arange(SENSIBILITY['CONDROCK_START'], SENSIBILITY['CONDROCK_STOP'], (SENSIBILITY['CONDROCK_STOP']-SENSIBILITY['CONDROCK_START'])/valuestep)
rhorockvect = SENSIBILITY['RHOROCK_START']

#Definition of Model Function

def Rs_matrix(condrockvect, rhorockvect, cprockvect):
    Rs, Rs_max = [], []
    for x in range(len(condrockvect)):
        alpha = condrockvect[x]/rhorockvect/cprockvect
        # calcolo di come cambia la resistenza del terreno in funzione del tempo
        RS = (1/2/condrockvect[x])*np.log((2 * np.sqrt(alpha * tt))/rext)
        Rs.append(RS)
        Rs_max.append(max(RS))
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
    twat_down_set_depth1 = []
    twat_horizontal_set_depth1 = []
    twat_up_set_depth1 = []
    twat_down_set_depth2 = []
    twat_horizontal_set_depth2 = []
    twat_up_set_depth2 = []
    twat_down_set_depth3 = []
    twat_horizontal_set_depth3 = []
    twat_up_set_depth3 = []
    den1 = portata * cpw
    for ii in range(1, len(zz)):
        for x in range(len(condrockvect)):
            alpha = condrockvect[x]/rhorockvect/cprockvect
            condrock=condrockvect[x] #W/m/K
        #The total coefficient is calculated every time because it changes with depth
        #The factor that changes is the resistance of the soil as its parameters changes with depth
            ktot = np.power(rext / condrock * np.log((2 * (np.sqrt(alpha * step))/rext)) + 1 / funzh(vel, dhydr),-1) #%W/m^2/K #far variare t
        #Thermal exchange between two depths at a distance dz
            num1 = zz[ii] * 2 * np.pi * rext * ktot
            aa = (num1 / (2 * den1)) + 1 # st_dz

            twat_down = (num1 / den1 * funztground(zz[ii]) - (num1 / (2 * den1)-1) * twat_in) / aa #C°
            twat_horizontal = (tground[ii] - ((tground[ii] - twat_down)/(np.exp((2 * np.pi * condrock * L)/(portata * cpw * np.log((np.pi * (alpha * time)**0.5)/(2 * rext))))))) #C°
            twat_out = (num1 / den1 * funztground(zz[ii]) - twat_horizontal * aa) / (num1 / (2 * den1)-1)  #C°
            print(zz[ii])
            print(twat_down)
            print(twat_horizontal)
            print(twat_out)
            if ii == 1:
                twat_down_set_depth1.append(twat_down)
                twat_horizontal_set_depth1.append(twat_horizontal)
                twat_up_set_depth1.append(twat_out)
            elif ii == 2:
                twat_down_set_depth2.append(twat_down)
                twat_horizontal_set_depth2.append(twat_horizontal)
                twat_up_set_depth2.append(twat_out)
            else:
                twat_down_set_depth3.append(twat_down)
                twat_horizontal_set_depth3.append(twat_horizontal)
                twat_up_set_depth3.append(twat_out)
    
    return twat_down_set_depth1, twat_horizontal_set_depth1, twat_up_set_depth1,  twat_down_set_depth2, twat_horizontal_set_depth2, twat_up_set_depth2,  twat_down_set_depth3, twat_horizontal_set_depth3, twat_up_set_depth3
 

#Definition of Graphic
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
            ax.legend(loc='upper right')
    plt.savefig('Output_Temperature_'+str(step/365/86400)+'y_'+str(NAME)+'.png', dpi=1200)
    plt.show()


def plot_3(twat_in, twat_down_set_depth1, twat_horizontal_set_depth1, twat_up_set_depth1, twat_down_set_depth2, twat_horizontal_set_depth2, twat_up_set_depth2, twat_down_set_depth3, twat_horizontal_set_depth3, twat_up_set_depth3, tground, zz, step):
    #fig = plt.figure(figsize=(19.20,10.80))
    #plt.subplots(figsize=(15,15))
    fig, ax = plt.subplots()
    ax.set_ylim(zz[-1]+1000, 0)
    ax.set_ylabel('Depth [m]')
    ax.set_xlabel('Temperature [°C]')
    ax.set_title('Evaluation {} years'.format(step/365/86400))
    ax.grid('both', linestyle='--')
    ax.plot(tground, zz,'g', label='T ground')
    for i in range(1, len(zz)):
        if i == 1:
            for j in range(len(twat_down_set_depth1)):
                if j == 0:
                    ax.plot([twat_in, twat_down_set_depth1[j]], [0, zz[i]], 'purple', label='T downward water (min)')
                    ax.plot([twat_down_set_depth1[j], twat_horizontal_set_depth1[j]], [zz[i], zz[i]],'y', label='T lateral water')
                    ax.plot([twat_horizontal_set_depth1[j], twat_up_set_depth1[j]], [zz[i], 0],'purple', label='T upward water (min)')
                elif j == (len(twat_down_set_depth1) - 1):
                    ax.plot([twat_in, twat_down_set_depth1[j]], [0, zz[i]], 'black', label='T downward water (max)')
                    ax.plot([twat_down_set_depth1[j], twat_horizontal_set_depth1[j]], [zz[i], zz[i]],'y')
                    ax.plot([twat_horizontal_set_depth1[j], twat_up_set_depth1[j]], [zz[i], 0],'black', label='T upward water (max)')
                    ax.legend(loc='lower left', prop={'size': 6})
                else:
                    if j == 1:
                        ax.plot([twat_in, twat_down_set_depth1[j]], [0, zz[i]], 'b', label='T downward water')
                        ax.plot([twat_horizontal_set_depth1[j], twat_up_set_depth1[j]], [zz[i], 0],'r', label='T upward water')
                    else:
                        ax.plot([twat_in, twat_down_set_depth1[j]], [0, zz[i]], 'b')
                        ax.plot([twat_horizontal_set_depth1[j], twat_up_set_depth1[j]], [zz[i], 0],'r')
                    ax.plot([twat_down_set_depth1[j], twat_horizontal_set_depth1[j]], [zz[i], zz[i]],'y') 

        elif i == 2:
            for j in range(len(twat_down_set_depth2)):
                if j == 0:
                    ax.plot([twat_in, twat_down_set_depth2[j]], [0, zz[i]], 'purple', label='T downward water')
                    ax.plot([twat_down_set_depth2[j], twat_horizontal_set_depth2[j]], [zz[i], zz[i]],'y', label='T lateral water')
                    ax.plot([twat_horizontal_set_depth2[j], twat_up_set_depth2[j]], [zz[i], 0],'purple', label='T upward water')
                elif j == (len(twat_down_set_depth2) - 1):
                    ax.plot([twat_in, twat_down_set_depth2[j]], [0, zz[i]], 'black', label='T downward water')
                    ax.plot([twat_down_set_depth2[j], twat_horizontal_set_depth2[j]], [zz[i], zz[i]],'y', label='T lateral water')
                    ax.plot([twat_horizontal_set_depth2[j], twat_up_set_depth2[j]], [zz[i], 0],'black', label='T upward water')
                else:
                    ax.plot([twat_in, twat_down_set_depth2[j]], [0, zz[i]], 'b', label='T downward water')
                    ax.plot([twat_down_set_depth2[j], twat_horizontal_set_depth2[j]], [zz[i], zz[i]],'y', label='T lateral water')
                    ax.plot([twat_horizontal_set_depth2[j], twat_up_set_depth2[j]], [zz[i], 0],'r', label='T upward water')
        else:
            for j in range(len(twat_down_set_depth3)):
                if j == 0:
                    ax.plot([twat_in, twat_down_set_depth3[j]], [0, zz[i]], 'purple', label='T downward water')
                    ax.plot([twat_down_set_depth3[j], twat_horizontal_set_depth3[j]], [zz[i], zz[i]],'y', label='T lateral water')
                    ax.plot([twat_horizontal_set_depth3[j], twat_up_set_depth3[j]], [zz[i], 0],'purple', label='T upward water')
                elif j == (len(twat_down_set_depth3) - 1):
                    ax.plot([twat_in, twat_down_set_depth3[j]], [0, zz[i]], 'black', label='T downward water')
                    ax.plot([twat_down_set_depth3[j], twat_horizontal_set_depth3[j]], [zz[i], zz[i]],'y', label='T lateral water')
                    ax.plot([twat_horizontal_set_depth3[j], twat_up_set_depth3[j]], [zz[i], 0],'black', label='T upward water')
                else:
                    ax.plot([twat_in, twat_down_set_depth3[j]], [0, zz[i]], 'b', label='T downward water')
                    ax.plot([twat_down_set_depth3[j], twat_horizontal_set_depth3[j]], [zz[i], zz[i]],'y', label='T lateral water')
                    ax.plot([twat_horizontal_set_depth3[j], twat_up_set_depth3[j]], [zz[i], 0],'r', label='T upward water')

    plt.savefig('Output_Temperature_'+str(step/365/86400)+'y_'+str(NAME)+'.png', dpi=1200)
    plt.show()


#Geothermal gradient function
funztground = lambda z : a*z + b #C°

#Function of the calculation of the coeff. convection in the pipes. I find Nusselt from Dittus-Boelter see
funzh = lambda v, drif : np.power((rhow * v * drif / viscw), 0.8) * np.power(pr, 0.4) * condw * 0.023 / (drif) #W/m^2/K

#Depth
depth = WELL['DEPTH']
#Tubes parameters
rext = WELL['EXT Radius']

#Operating time
T = 31536000
tt = np.arange(WELL['START']*T, WELL['STOP']*T, WELL['STEP']*T) #s

#discretization of depth
dz = WELL['DZ'] #m
zz = np.arange(0,depth,dz) #m

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
condiso = WELL['VISC WATER'] #W/m/K

pr = cpw * viscw / condw


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
#Definition of the Rs function with depth over time 
print("hhh")
print(condrockvect)
Rs, Rs_max = Rs_matrix(condrockvect, rhorockvect, cprockvect)
value = find_max(Rs, Rs_max)
print(Rs)
print(Rs_max)

cp_Rs = []
##DEFINIZIONE GRAFICO
#for x in range(len(cprockvect)):
    #alpha = condrockvect[x]/rhorockvect[x]/cprockvect[x]
    #Rs = (1/2/condrockvect[x])*np.log((2 * np.sqrt(alpha * tt))/rext)
   # plt.plot(tt/3600/24/365, Rs, label = "α {}".format(round(alpha,8)))
    #plt.xlabel('Time [years]')    
    #plt.ylabel('Thermal Resistance - $R_s$ [$m^2K/W$]')
    #plt.grid('both', linestyle='--')
   # cp_Rs.append(Rs)
#plt.legend(fontsize='x-small', title='Thermal diffusivity [$m^2/s$]', title_fontsize='x-small')
#plt.savefig("alpha_"+NAME+".png", dpi=1200)
#plt.savefig("alpha_"+NAME+".svg", dpi=1200)
#plt.show()




#Loop for ring temp. profile
for step in t: 
    twat_down_set_depth1, twat_horizontal_set_depth1, twat_up_set_depth1, twat_down_set_depth2, twat_horizontal_set_depth2, twat_up_set_depth2, twat_down_set_depth3, twat_horizontal_set_depth3, twat_up_set_depth3 = twat_evaluation(zz, step, time, rext, vel, dhydr, portata, cpw, tground, twat_in, L)
    plot_3(twat_in, twat_down_set_depth1, twat_horizontal_set_depth1, twat_up_set_depth1, twat_down_set_depth2, twat_horizontal_set_depth2, twat_up_set_depth2, twat_down_set_depth3, twat_horizontal_set_depth3, twat_up_set_depth3, tground, zz, step)


