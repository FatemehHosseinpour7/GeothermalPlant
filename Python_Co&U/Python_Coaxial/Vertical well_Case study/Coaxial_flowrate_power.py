import numpy as np
import matplotlib.pyplot as plt
from lib.settings import NAME, WELL_VILLAFORTUNA1, WELL_TRECATE4

#DEFINIZIONE FUNZIONI MODELLO

def Rs_line(Rs, Rs_max, alpha, tt, condrock, rext):
    rs_function = lambda tt: (1/2/condrock)*np.log((2 * np.sqrt(alpha * tt))/rext)
    Rs.append(rs_function(tt))
    Rs_max.append(max(rs_function(tt)))
    return Rs, Rs_max


def Rs_matrix(zz, tt, rext):
    Rs, Rs_max = [], []
    # Definisco alpha and RS
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
        alpha = condrock / rhorock / cprock
        
        #Il coeff.totale viene calcolato ogni volta perch cambia con la profondità. Il fattore che cambia  la resistenza del
        #terreno poichè i suoi parametri cambiano con la profondit
        ktot = np.power(rext / condrock * np.log((2 * (np.sqrt(alpha * step))/rext)) + 1 / funzh(vel, dhydr),-1) #%W/m^2/K #far variare t 
        
        #Scambio termico tra due profondità a distanza dz
        num1 = dz * 2 * np.pi * rext * ktot
        den1 = portata * cpw
        aa = (num1 / (2 * den1)) + 1 # st_dz
        
        #Temperatura acqua dopo scambio termico (vettore di temperatura)
        twat_value = (num1 / den1 * funztground(zz[ii] - dz / 2) - (num1 / (2 * den1)-1) * twat_value) / aa #C°
        twat_down.append(twat_value)
    return twat_down


def twat_up_evaluation(zz, twat, step, rint, riso, condiso, portata, cpw):
    # inverto i vettori perch ora partiamo dal fondo e risaliamo e quindi quando 
    # noi vogliamo il primo valore di risalita ci dobbiamo riferire all'ultima
    # temperatura di discesa
    zz2 = np.flip(zz) #m
    twatdalbasso = np.flip(twat)
    vel2 = portata / (np.pi * np.power(rint,2)) / rhow #m/s
    hint = funzh(vel2, 2*rint) #W/m^2/K
    hexttubo = funzh(vel, 2*riso) #W/m^2/K 
    # coeff. di scambio termico comprensivo delle due resistenze convetive e di
    # quella conduttiva dovuta all'isolante
    ktotinner = np.power(rint / riso * 1/ hexttubo + rint * np.log(riso / rint) * 1 / condiso + 1 / hint,-1) #W/m^2/K
    aa1 = (dz * 2 * np.pi * rint * ktotinner / (2 * portata * cpw)) + 1
    # loop per l'acqua che torna in superficie il primo valore � uguale
    # all'ultimo del calcolo per l'acqua che arriva al fondo
    twat_up_value = twat[-1] #C°
    twat_up = []
    for jj in range(0, len(zz2)):
        twat_up_value = (dz * 2 * np.pi * rint * ktotinner / (portata * cpw) * ((twatdalbasso[jj-1] + twatdalbasso[jj]) / 2) - (dz * 2 * np.pi * rint * ktotinner / (2 * portata * cpw) - 1) * twat_up_value) / aa1 #C°
        twat_up.append(twat_up_value)
    return twat_up


def plot_1(tt, value):
    f = plt.figure()
    f.set_figwidth(20)
    f.set_figheight(10)
    plt.plot(tt/3600/24/365, value, 'g')
    plt.xlabel('Time [years]')
    plt.ylabel('Thermal Resistance - $R_s$ [$m^2K/W$]')
    plt.grid('both', linestyle='--')
    plt.xticks([1,2,3,4,5,6,7,8,9,10])
    plt.show()


def plot_2(twat, twat2, tground, zz, step):
    fig, ax = plt.subplots()
    ax.plot(twat, zz, 'b', label='t downward water')
    ax.plot(twat2, np.flip(zz), 'r', label='t upward water')
    ax.plot(tground, zz, 'y', label='t ground')
    ax.set_ylim(zz[-1]+1000, 0)  # decreasing time
    ax.set_xlabel('Temperature [°C]')
    ax.set_ylabel('Depth [m]')
    ax.set_title('Evaluation {} years'.format(step/365/86400))
    ax.grid('both', linestyle='--')
    ax.legend()
    plt.show()

def plot_3(portata, value, name_xvalue):
    # f = plt.figure()
    # f.set_figwidth(20)
    # f.set_figheight(10)
    plt.plot(portata, value, 'r')
    plt.ylabel(name_xvalue)
    # plt.yticks(value, None, fontsize=7)
    plt.xticks(fontsize=7)
    plt.ylim((min(twat_up), max(twat_up)))
    plt.xlim(0.5,(portata_step[-1]))
    plt.xlabel('Flowrate [$kgs^{-1}$]')
    plt.grid('both', linestyle='--')
    # plt.yticks(portata)
    plt.yticks(np.arange(int(min(value)-3), max(value)+8, step=10))
    # plt.set_xlim(1.5,2.5)
    plt.savefig("Coaxial_flowrate_"+NAME+".png", dpi=1800)
    plt.show()
    
def plot_4(portata, value, name_xvalue):
    # f = plt.figure()
    # f.set_figwidth(20)
    # f.set_figheight(10)
    plt.plot(portata, value, 'r')
    plt.ylabel(name_xvalue)
    # plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0)) #scientific notation
    # plt.xlabel(round(name_xvalue))
    plt.yticks(value, None, fontsize=7)
    #limitazione asse verticale
    plt.ylim((min(value), max(value)-2))
    #plt.ylim((78, max(value)))
    plt.xticks(fontsize=7)
    plt.xlim(0.5,(portata_step[-1]))
    #plt.xlim(1,(portata_step[-4]))
    plt.xlabel('Flowrate [$kgs^{-1}$]')
    plt.grid('both', linestyle='--')
    #plt.yticks(portata)
    #plt.set_xlim(1.5,2.5)
    plt.savefig("Coaxial_power_"+NAME+".png", dpi=1800)
    #plt.savefig("Coaxial_power_"+NAME+".svg", dpi=1200)
    plt.show()
    
# def plot_3(portata, value, name_xvalue):
#     f = plt.figure()
#     f.set_figwidth(20)
#     f.set_figheight(10)
#     plt.plot(value, portata, 'r')
#     plt.xlabel(name_xvalue)
#     plt.xticks(value, None, rotation='vertical')
#     print(min(value))
#     # plt.xlim((min(t_out), max(t_out)))
#     plt.ylabel('Portata [kgs-1]')
#     plt.grid('both', linestyle='--')
#     plt.xticks(np.arange(int(min(value)-1), max(value), step=5))
#     #plt.xticks(t_out)
#     plt.xlim((36, max(value)))
#     #plt.set_xlim(1.5,2.5)
#     plt.savefig("Coaxial_flowrate.png", dpi=300)
#     plt.show()

# def plot_4(portata, value, name_xvalue):
#     f = plt.figure()
#     f.set_figwidth(20)
#     f.set_figheight(10)
#     plt.plot(value, portata, 'r')
#     plt.xlabel(name_xvalue)
#     plt.xticks(value, None, rotation='vertical')
#     print(min(value))
#     # plt.xlim((min(value), max(value)))
#     plt.ylabel('Portata [kgs-1]')
#     plt.grid('both', linestyle='--')
#     plt.xticks(value)
#     plt.xlim((10732.558205503105, 438835.6740207797))
#     #plt.set_xlim(1.5,2.5)
#     plt.savefig("Coaxial_power.png", dpi=300)
#     plt.show()
    
#Parametri relativi alle formazioni geologiche
if NAME == 'Villafortuna1':
    condrockfunz = lambda bb: 0.3*(bb<1258)+1.61*(bb>=1258)*(bb<1405)+3.16*(bb>=1405)*(bb<5493)+3.5*(bb>=5493)*(bb<=depth) #W/m/K
    rhorockfunz = lambda bb: 1700*(bb<1258)+1890*(bb>=1258)*(bb<1405)+2359*(bb>=1405)*(bb<5493)+2480*(bb>=5493)*(bb<=depth) #kg/m^3
    cprockfunz = lambda bb: 800*(bb<1258)+1693*(bb>=1258)*(bb<1405)+821.111*(bb>=1405)*(bb<5493)+810.484*(bb>=5493)*(bb<=depth) #j/kg/k
    WELL = WELL_VILLAFORTUNA1
elif NAME == 'Trecate4':
    condrockfunz = lambda bb: 1.61*(bb<1632)+2.17*(bb>=1632)*(bb<5106)+3.50*(bb>=5106)*(bb<6189)+3.00*(bb>=6189)*(bb<=depth) #W/m/K
    rhorockfunz = lambda bb: 1890*(bb<1632)+1801*(bb>=1632)*(bb<5106)+2480*(bb>=5106)*(bb<6189)+2330*(bb>=6189)*(bb<=depth) #kg/m^3
    cprockfunz = lambda bb: 1696*(bb<1632)+830.09*(bb>=1632)*(bb<5106)+810.48*(bb>=5106)*(bb<6189)+821.11*(bb>=6189)*(bb<=depth) #j/kg/k
    WELL = WELL_TRECATE4


#Funzione gradiente geotermico
funztground = lambda z : a*z + b #C°

#Funzione del calcolo del coeff. convettivo nei tubi. Trovo Nusselt da dittus-boelter 
funzh = lambda v, drif : np.power((rhow * v * drif / viscw), 0.8) * np.power(pr, 0.4) * condw * 0.023 / (drif) #W/m^2/K

#profondità
depth = WELL['DEPTH']
#tubes parameters
rext = WELL['EXT Radius']
riso = WELL['ISO Radius']
rint = WELL['INT Radius']
#tempo di funzionamento
T = 31536000
tt = np.arange(WELL['START']*T, WELL['STOP']*T, WELL['STEP']*T) #s
#discretizzazione della profondit�
dz = WELL['DZ'] #m
zz = np.arange(0,depth,dz) #m 

####################### Analisi scambiatore
#portata immessa
portata = WELL['PORTATA']
portata_step = np.arange(0,portata+0.2,0.2) #m 
# imposto la temperatura di entrata dell'acqua nel pozzo
twat_value = WELL['T WATER IMMESSA']
#step relativi agli anni si analisi
t = []
for x in WELL['STEP T']:
    t.append(x*T) #s parametro che io posso far variare
# profilo temperatura del terreno in relazione alla profondit�
a, b = WELL['A'], WELL['B']
# water parameters at 100C� and 2 bar
cpw = WELL['CPW'] #j/kg/k
condw = WELL['COND WATER'] #W/m/K
rhow = WELL['RHOW'] #kg/m^3
viscw = WELL['VISC WATER'] #Pa s
condiso = WELL['VISC WATER'] #W/m/K

pr = cpw * viscw / condw

########### Blocco 1 ##########################
# Definizione della funzione Rs con la profonditànel tempo
Rs, Rs_max = Rs_matrix(zz, tt, rext)
# #individuazione della profondit� della massima Rs
value = find_max(Rs, Rs_max)
# creazione grafico
plot_1(tt, Rs[value[0]])


#####################à## Analisi scambiatore
# coeff. di scambio termico totale comprendendo terreno e parte convettiva nell'anello 
# diametro idraulico per il calcolo dei coeff. adimensionali dell'anello
dhydr = (4 * np.pi * (np.power(rext, 2) - np.power(riso, 2))) / (2 * np.pi * (rext + riso)) #m 
areac = np.pi * (np.power(rext, 2) - np.power(riso, 2)) #%m^2


tout, P = [], [] #creo lista valori
for x in range (len(portata_step)):
    portata = portata_step[x]
    vel = portata / areac / rhow #m/s

    # coefficiente di scambio termico convettivo
    hext = funzh(vel, dhydr) #W/m^2/K
    # calcolo la temperatura del suolo con la profondit�
    tground1 = funztground(zz) #C�

    # # loop per profilo temp anello;
    # for step in t: 
    twat_down = twat_down_evaluation(zz, t[0], rext, vel, dhydr, portata, cpw, tground1, twat_value)
    twat_up = twat_up_evaluation(zz, twat_down, t[0], rint, riso, condiso, portata, cpw)
    
    print(portata, twat_up[-1])
    tout.append(twat_up[-1])
    P.append((portata*cpw*(twat_up[-1]- twat_down[1])/1000))

print(P)
print(len(tout), len(P))
plot_3(portata_step[1:], tout[1:], 'Temperature [°C]')
# plot_3(tout[1:], portata_step[1:], 'Temperature [°C]')
plot_4(portata_step[1:], P[1:], 'Heat Power [kW]')

