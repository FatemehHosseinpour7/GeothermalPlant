import numpy as np
import matplotlib.pyplot as plt
from lib.settings import NAME, WELL_TEMPAROSSA, WELL_VILLAFORTUNA1, WELL_TRECATE4, WELL_GELA, WELL_CASTEGGIO, WELL_CAMMARATA, WELL_BENIGNO, WELL_CINZANO

#DEFINIZIONE FUNZIONI MODELLO

def Rs_line(Rs, Rs_max, alpha, tt, condrock, rext):
    rs_function = lambda tt: (1/2/condrock)*np.log((2 * np.sqrt(alpha * tt))/rext)
    Rs.append(rs_function(tt))
    Rs_max.append(max(rs_function(tt)))
    return Rs, Rs_max


def Rs_matrix(zz, tt, rext):
    Rs, Rs_max = [], []
    # Definisco alpha and Rx
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


def twat_Udown_evaluation(zz, step, rext, vel, dhydr, portata, cpw, tground, twat_value, resconv, resgrout):
    twat_down = []
    for ii in range(len(zz)):
        condrock=condrockfunz(zz[ii]) #W/m/K
        rhorock=rhorockfunz(zz[ii]) #kg/m^3
        cprock=cprockfunz(zz[ii]) #j/kg/k
        alpha = condrock / rhorock / cprock
        #in questo caso il coeff.totale viene calcolato ogni volta perchè cambia con la profondità. Il fattore che cambia la resistenza del
        #terreno poichè i suoi parametri cambiano con la profondità
        resterreno = rext / condrock * np.log(2 * np.power(alpha * step, 0.5)/rext) #m^2*K/W
        ktot = np.power((resconv + resgrout + resterreno), -1) #(m^2*K/W)^-1
        # print(ktot)
        # scambio termico tra due profondit� a distanza dz
        num1 = dz * np.pi * rext * ktot
        den1 = portata * cpw
        aa = (num1 / (2 * den1)) + 1 # st_dz
        # temperatura acqua dopo scambio termico (vettore di temperatura)
        twat_value = (num1 / den1 * funztground(zz[ii] - dz / 2) - (num1 / (2 * den1)-1) * twat_value) / aa #C°
        twat_down.append(twat_value)
    return twat_down


def twat_Uup_evaluation(zz, twat, step, rext, vel, dhydr, portata, cpw, tground, twat_value, resconv, resgrout):
    # inverto i vettori perchè ora partiamo dal fondo e risaliamo e quindi quando
    # noi vogliamo il primo valore di risalita ci dobbiamo riferire all'ultima
    # temperatura di discesa
    zz = np.flip(zz) #m
    twatdalbasso = np.flip(twat)
    twat_value = twat[-1] #C�
    twat_up = []
    for ii in range(0, len(zz)):
        condrock=condrockfunz(zz[ii]) #W/m/K
        rhorock=rhorockfunz(zz[ii]) #kg/m^3
        cprock=cprockfunz(zz[ii]) #j/kg/k
        alpha = condrock / rhorock / cprock
        
        #Il coeff.totale viene calcolato ogni volta perchè cambia con la profondità. Il fattore che cambia la resistenza del
        #terreno poichè i suoi parametri cambiano con la profondit�
        resterreno = rext / condrock * np.log((2 * (np.sqrt(alpha * step))/rext)) #m^2*K/W
        ktot = np.power((resconv + resgrout + resterreno), -1) #(m^2*K/W)^-1

        #Scambio termico tra due profondità a distanza dz
        num1 = dz * np.pi * rext * ktot
        den1 = portata * cpw
        aa = (num1 / (2 * den1)) + 1 # st_dz
        
        #Temperatura acqua dopo scambio termico (vettore di temperatura)
        twat_value = (num1 / den1 * funztground(zz[ii] - dz / 2) - (num1 / (2 * den1)-1) * twat_value) / aa #C�
        twat_up.append(twat_value)
    return twat_up

#DEFINIZIONE GRAFICI

def plot_1(tt, value):
    # f = plt.figure()
    # f.set_figwidth(20)
    # f.set_figheight(10)
    plt.plot(tt/3600/24/365, value, 'g')
    plt.xlabel('Time [years]')
    plt.ylabel('Thermal Resistance - Rs [m^2K/W]')
    plt.grid('both', linestyle='--')
    plt.xticks([1,2,3,4,5,6,7,8,9,10])
    #plt.savefig('ThermalResistance_'+str(NAME)+'.png', dpi=1200)
    plt.show()


def plot_2(twat, twat2, tground, zz, step):
    fig, ax = plt.subplots()
    ax.plot(twat, zz, 'b', label='T downward water')
    ax.plot(twat2, np.flip(zz), 'r', label='T upward water')
    ax.plot(tground, zz, 'g', label='T ground')
    ax.set_ylim(zz[-1]+1000, 0) # decreasing time
    ax.set_xlabel('Temperature [°C]')
    ax.set_ylabel('Depth [m]')
    ax.set_title('Evaluation {} years'.format(step/365/86400))
    ax.grid('both', linestyle='--')
    ax.legend()
    plt.savefig('Utube_profile_'+str(step/365/86400)+'y_'+str(NAME)+'.png', dpi=1200)
    plt.show()


#Parametri delle formazioni geologiche
if NAME == 'Temparossa':
    condrockfunz = lambda bb: 3.0*(bb<2912)+2.34*(bb>=2912)*(bb<=depth) #W/m/K
    rhorockfunz = lambda bb: 2330*(bb<2912)+1917*(bb>=2912)*(bb<=depth) #kg/m^3
    cprockfunz = lambda bb: 808.6*(bb<2912)+829.4*(bb>=2912)*(bb<=depth) #j/kg/k
    WELL = WELL_TEMPAROSSA
elif NAME == 'Villafortuna1':
    condrockfunz = lambda bb: 0.3*(bb<1258)+1.61*(bb>=1258)*(bb<1405)+3.16*(bb>=1405)*(bb<5493)+3.5*(bb>=5493)*(bb<=depth) #W/m/K
    rhorockfunz = lambda bb: 1700*(bb<1258)+1890*(bb>=1258)*(bb<1405)+2359*(bb>=1405)*(bb<5493)+2480*(bb>=5493)*(bb<=depth) #kg/m^3
    cprockfunz = lambda bb: 800*(bb<1258)+1696*(bb>=1258)*(bb<1405)+821.111*(bb>=1405)*(bb<5493)+810.484*(bb>=5493)*(bb<=depth) #j/kg/k
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
    

#Funzione radiente geotermico
funztground = lambda z : a*z + b #C°

#Funzione del calcolo del coeff. convettivo nei tubi. Trovo Nusselt da dittus-boelter
funzh = lambda v, drif : np.power((rhow * v * drif / viscw), 0.8) * np.power(pr, 0.4) * condw * 0.023 / (drif) #W/m^2/K
funznu = lambda v, drif : np.power((rhow * v * drif / viscw), 0.8) * np.power(pr, 0.4) * condw * 0.023


# profondità
depth = WELL['DEPTH']
#tubes parameters
rwell = WELL['WELL Radius']
rtube = WELL['TUBE Radius']
w = WELL['W Distance']
#tempo di funzionamento
T = 31536000
tt = np.arange(WELL['START']*T, WELL['STOP']*T, WELL['STEP']*T) #s
#discretizzazione della profondità
dz = WELL['DZ'] #m
zz = np.arange(0,depth,dz) #m

####################### Analisi scambiatore
# portata immessa
portata = WELL['PORTATA']
# imposto la temperatura di entrata dell'acqua nel pozzo
twat_value = WELL['T WATER IMMESSA']

# step relativi agli anni si analisi
t = []
for x in WELL['STEP T']:
    t.append(x*T) #s parametro che io posso far variare
# profilo temperatura del terreno in relazione alla profondità
a, b = WELL['A'], WELL['B']
# water parameters at 100C° and 2 bar
cpw = WELL['Utube CPW'] #j/kg/k
condw = WELL['Utube COND WATER'] #W/m/K
rhow = WELL['Utube RHOW'] #kg/m^3
viscw = WELL['Utube VISC WATER'] #Pa s
condgrout = WELL['COND GROUT'] #W/m/K
pr = cpw * viscw / condw

# area di resistenza tot (modello U-tube)
areac = np.pi * np.power(rtube,2)
vel = portata / areac / 1000 #m/s

# area di scambio
areascambio = np.pi * rwell * dz
deq = 2 * rtube * np.power((4 * w)/(np.pi * rtube) +1, 0.5)

# coeff. di scambio termico totale comprendendo terreno e parte convettiva nell'anello
dhydr = 2 * rtube 
hext = funzh(vel, dhydr)
hu_ext = funznu(vel, dhydr)
resgrout = np.log(2 * rwell / deq) / (np.pi * condgrout * dz) * areascambio #K*m^2/W
resconv = np.power((np.pi * dz * hu_ext * condw),-1) * areascambio #K*m^2/W

# calcolo la temperatura del suolo con la profondità (richiamo funzione)
tground1 = funztground(zz) #C°

########### Blocco 1 ##########################
# Definizione della funzione RS con la profondità nel tempo
Rs, Rs_max = Rs_matrix(zz, tt, rwell)
# individuazione della profondità della massima Rs
value = find_max(Rs, Rs_max)
# creazione grafico
plot_1(tt, Rs[value[0]])

for step in t:
    twat_down = twat_Udown_evaluation(zz, step, rwell, vel, dhydr, portata, cpw, tground1, twat_value, resconv, resgrout)
    # print(twat_down)
    twat_up = twat_Uup_evaluation(zz, twat_down, step, rwell, vel, dhydr, portata, cpw, tground1, twat_value, resconv, resgrout)
    plot_2(twat_down, twat_up, tground1, zz, step)

print(twat_up)