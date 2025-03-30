import matplotlib.pyplot as plt
from numpy import *
import math
import numpy as np
import matplotlib.animation as animation
from scipy.integrate import odeint
pi = math.pi

global G, M_Terre, R_Terre, mu

#constante
G = 6.674*10**(-11)
M_Terre = 5.972*10**24
R_Terre = 6370*10**3
mu = G*M_Terre
rho_0 = 1.292
r = 287.058
rho_kerosene=800

#CI
angle = 5*pi/12
x = R_Terre
y= 0
v_x = sqrt(mu/(R_Terre+400000))*2
v_y = 0#sqrt(mu/(R_Terre))*2
alt_max=400000*10

fusee = [x,v_x,y,v_y]

R_0 = [x,v_x,y,v_y]

#precision
dtheta = 1000
dt = 120*8
t_f = 2

t = np.linspace(0,t_f,dt)
theta = np.linspace(0,2*pi,dtheta)

R_x = np.linspace(0,2*pi,dtheta)
R_y = np.linspace(0,2*pi,dtheta)
Rh_x = np.linspace(0,2*pi,dtheta)
Rh_y = np.linspace(0,2*pi,dtheta)

for i in range(len(theta)):
    R_x[i] = R_Terre*math.cos(theta[i])
    R_y[i] = R_Terre*math.sin(theta[i])
    Rh_x[i] = (R_Terre+alt_max)*math.cos(theta[i])
    Rh_y[i] = (R_Terre+alt_max)*math.sin(theta[i])

dic_atmosphere = {"altitude":[R_Terre],"rho":[rho_0]}

def nul(Y,dY,t):
    return 0

def atmosphere(z):
    if z>=0 and z<=11e3:
        return rho_0*(1-70*10**(-3)/11/288*z)**(9.81/(70*10**(-3)/11*r))
    elif z>11e3 and z<=20e3:
        return atmosphere(11e3)*exp(-9.81/(r*218)*(z-11e3))
    elif z>20e3 and z<=32e3:
        return atmosphere(20e3)*((9*10**(-3)*z+218)/(9*10**(-3)*20e3+218))**(-9.91/(9*10**(-3)*r))
    elif z>32e3 and z<=48.5e3:
        return atmosphere(32e3)*(((42/16.5)*10**(-3)*z+230)/((42/16.5)*10**(-3)*32e3+230))**(-9.91/((42/16.5)*10**(-3)*r))
    elif z>48.5e3 and z<=53e3:
        return atmosphere(48.5e3)*exp(-9.81/(r*270)*(z-48.5e3))
    else :
        return 0

def collision(c_1,c_2):
    # d = sqrt((c_1["x"]-c_2["x"])*(c_1["x"]-c_2["x"])+(c_1["y"]-c_2["y"])*(c_1["y"]-c_2["y"]))
    return sqrt((c_1["x"][-1]-c_2["x"][-1])*(c_1["x"][-1]-c_2["x"][-1])+(c_1["y"][-1]-c_2["y"][-1])*(c_1["y"][-1]-c_2["y"][-1]))>=c_1["rayon"]+c_2["rayon"]

def gravitation(Y,dY,t,m=M_Terre):

    # a_x = -mu*Y[0]/sqrt(Y[0]*Y[0]+Y[1]*Y[1])**3
    # a_y = -mu*Y[1]/sqrt(Y[0]*Y[0]+Y[1]*Y[1])**3

    # return [dY[0],a_x,dY[1],a_y]
    return [dY[0],-G*m*Y[0]/sqrt(Y[0]*Y[0]+Y[1]*Y[1])**3,dY[1],-G*m*Y[1]/sqrt(Y[0]*Y[0]+Y[1]*Y[1])**3]

def gravi_fusee(Y,dY,t):
    dx=0
    dvx = 0
    dy=0
    dvy = 0

    for i in corps:
        if collision(corps[i],objet["fusee"])==True:
            # x = (Y[0]-corps[i]["x"][int(t)])
            # y = (Y[1]-corps[i]["y"][int(t)])
            # print(i,"x",corps[i]["x"][int(t/dt)],"y",corps[i]["y"][int(t/dt)],"t",int(t/dt),t)
            acc = gravitation([Y[0]-corps[i]["x"][int(t/dt)],Y[1]-corps[i]["y"][int(t/dt)]],dY,t,corps[i]["masse"])
        else :
            acc = [0,0,0,0]
            print("crash")
        dx+= acc[0]
        dvx+=acc[1] 
        dy+=acc[2]
        dvy+=acc[3]
    
    return [dx,dvx,dy,dvy]

def frottement(Y,dY,t):

    # ax = 0.5*atmosphere(sqrt(Y[0]*Y[0]+Y[1]*Y[1])-R_Terre)*0.4*1866*(dY[0]*dY[0]+dY[1]*dY[1])*Y[0]/sqrt(Y[0]*Y[0]+Y[1]*Y[1])
    # ay = 0.5*atmosphere(sqrt(Y[0]*Y[0]+Y[1]*Y[1])-R_Terre)*0.4*1866*(dY[0]*dY[0]+dY[1]*dY[1])*Y[1]/sqrt(Y[0]*Y[0]+Y[1]*Y[1])

    return [dY[0],-0.5*atmosphere(sqrt(Y[0]*Y[0]+Y[1]*Y[1])-R_Terre)*0.4*(5.375/2)**2*pi*(dY[0]*dY[0]+dY[1]*dY[1])*Y[0]/sqrt(Y[0]*Y[0]+Y[1]*Y[1]),\
            dY[1],-0.5*atmosphere(sqrt(Y[0]*Y[0]+Y[1]*Y[1])-R_Terre)*0.4*(5.375/2)**2*pi*(dY[0]*dY[0]+dY[1]*dY[1])*Y[1]/sqrt(Y[0]*Y[0]+Y[1]*Y[1])]

def Poussee(ve):
    return ve*ve*rho_kerosene*(2.15/2)*(2.15/2)*pi

def Force(Y,dY,t,liste_force=[gravitation]):
    somme = [0]*(len(Y)+len(dY))
    for fct in liste_force:
        somme += fct(Y,dY,t)
    return somme
    # somme = liste_force[0](Y,dY,t)
    # return somme

def acceleration(Y,dY,t,ve):
    l_acc = gravi_fusee(Y,dY,t)+Poussee(ve)/objet["fusee"]["masse"][-1]
    return []

def changement_orbite(Y):
    if (Y[0]*Y[3]-Y[1]*Y[2])/(sqrt(Y[0]*Y[0]+Y[2]*Y[2])*sqrt(Y[1]*Y[1]+Y[3]*Y[3]))>=0:
        signe = 1
    else :
        signe =-1
    # ctheta = Y[0]/sqrt(Y[0]*Y[0]+Y[2]*Y[2])
    # stheta = Y[2]/sqrt(Y[0]*Y[0]+Y[2]*Y[2])
    # v = sqrt(mu/sqrt(Y[0]*Y[0]+Y[2]*Y[2]))*signe*sqrt(1.3)
    # return [stheta,ctheta],v
    return [Y[2]/sqrt(Y[0]*Y[0]+Y[2]*Y[2]),Y[0]/sqrt(Y[0]*Y[0]+Y[2]*Y[2])],sqrt(mu/sqrt(Y[0]*Y[0]+Y[2]*Y[2]))*signe*sqrt(1.3)

def RK4(fct,Y,t_f,h,a=1,dt=1):
    temps = [0]
    # Yx = [Y[0]]
    # Yy = [Y[2]]
    # Yvx = [Y[1]]
    # Yvy = [Y[3]]
    # y1 = [Y[0],Y[2]]
    # d_y1 = [Y[1],Y[3]]
    # Yx = [Y["x"][0]]
    # Yy = [Y["y"][0]]
    # Yvx = [Y["vx"][0]]
    # Yvy = [Y["vy"][0]]
    y1 = [Y["x"][0],Y["y"][0]]
    d_y1 = [Y["vx"][0],Y["vy"][0]]
    y2 = [0,0]
    d_y2 = [0,0]
    y3 = [0,0]
    d_y3 = [0,0]
    y4 = [0,0]
    d_y4 = [0,0]
    t = 0
    # a = 1
    while t < t_f:
        if sqrt(y1[0]*y1[0]+y1[1]*y1[1])>= R_Terre and sqrt(y1[0]*y1[0]+y1[1]*y1[1])<= h:
	
            p1 = fct(y1,d_y1,t)
            for i in range(len(y1)):
                y2[i] = y1[i] + dt/2*d_y1[i]
                d_y2[i] = d_y1[i] + dt/2*p1[2*i+1]
            
            p2 = fct(y2,d_y2,t)
            for i in range(len(y1)):
                y3[i] = y1[i] + dt/2*d_y1[i] + dt**2/4*p1[2*i+1]
                d_y3[i] = d_y1[i] + dt/2*p2[2*i+1]
            
            p3 = fct(y3,d_y3,t)
            for i in range(len(y1)):
                y4[i] = y1[i] +dt*d_y1[i] + dt**2/2*p2[2*i+1]
                d_y4[i] = d_y1[i] +dt*p3[2*i+1]
            
            p4 = fct(y4,d_y4,t)
            for i in range(len(y1)):
                y1[i] = y1[i] + dt*d_y1[i] +dt**2*(p1[2*i+1]+p2[2*i+1]+p3[2*i+1])/6
                d_y1[i] = d_y1[i] + dt*(p1[2*i+1]+2*p2[2*i+1]+2*p3[2*i+1]+p4[2*i+1])/6
        
        elif sqrt(y1[0]*y1[0]+y1[1]*y1[1])>= R_Terre+h and a == 1:
            angle,v = changement_orbite([Y["x"][-1],Y["vx"][-1],Y["y"][-1],Y["vy"][-1]])
            d_y1[0]=v*angle[0]
            d_y1[1]=v*angle[1]
            y1[0] = Y["x"][-1]+d_y1[0]*dt
            y1[1] = Y["y"][-1]+d_y1[1]*dt
            a = 0

        elif sqrt(y1[0]*y1[0]+y1[1]*y1[1])>= R_Terre :
            p1 = fct(y1,d_y1,t)
            for i in range(len(y1)):
                y2[i] = y1[i] + dt/2*d_y1[i]
                d_y2[i] = d_y1[i] + dt/2*p1[2*i+1]
            
            p2 = fct(y2,d_y2,t)
            for i in range(len(y1)):
                y3[i] = y1[i] + dt/2*d_y1[i] + dt**2/4*p1[2*i+1]
                d_y3[i] = d_y1[i] + dt/2*p2[2*i+1]
            
            p3 = fct(y3,d_y3,t)
            for i in range(len(y1)):
                y4[i] = y1[i] +dt*d_y1[i] + dt**2/2*p2[2*i+1]
                d_y4[i] = d_y1[i] +dt*p3[2*i+1]
            
            p4 = fct(y4,d_y4,t)
            for i in range(len(y1)):
                y1[i] = y1[i] + dt*d_y1[i] +dt**2*(p1[2*i+1]+p2[2*i+1]+p3[2*i+1])/6
                d_y1[i] = d_y1[i] + dt*(p1[2*i+1]+2*p2[2*i+1]+2*p3[2*i+1]+p4[2*i+1])/6

        t = t + dt
        
        temps.append(t)
        
        # Yx.append(y1[0])
        # Yy.append(y1[1])
        # Yvx.append(d_y1[0])
        # Yvy.append(d_y1[1])
        Y["x"].append(y1[0])
        Y["y"].append(y1[1])
        Y["vx"].append(d_y1[0])
        Y["vy"].append(d_y1[1])
        corps["Terre"]["x"].append(0)
        corps["Terre"]["y"].append(0)

    # return [Yx,Yy],temps
    return Y,temps

def Euler(fct,Y,t_f,h,a=1):
    temps = [0]

while dic_atmosphere["altitude"][-1]< R_Terre+53e3:
    dic_atmosphere["altitude"].append(dic_atmosphere["altitude"][-1]+100)
    dic_atmosphere["rho"].append(atmosphere(dic_atmosphere["altitude"][-1]-R_Terre))


# sol = odeint(orbite,R_0,t)
# solution = RK4(gravitation,fusee,10000,400000)
objet = {"fusee":{"masse":[10000],"x":[R_Terre],"y":[0],"vx":[v_x],"vy":[v_y],"rayon":0,"Force":[gravitation],"temps":[0]}}
corps = {"Terre":{"masse":M_Terre,"x":[0],"y":[0],"vx":[0],"vy":[0],"rayon":R_Terre,"Force":[nul]},\
         "Lune":{"masse":7.347*10**22,"x":[406300000],"y":[0],"vx":[0],"vy":[-995],"rayon":1734400,"Force":[nul]}}


t_max = 86400*60*4
Y,temps=RK4(gravitation,corps["Lune"],t_max,400000,2,dt)
RK4(gravi_fusee,objet["fusee"],t_max,alt_max*90,1,dt)

def affichage():
    plt.figure(3)
    plt.plot(corps["Lune"]["x"],corps["Lune"]["y"])
    plt.plot(R_x,R_y)
    plt.title("orbite de la Lune")

    plt.figure(4)
    plt.plot(dic_atmosphere["altitude"],dic_atmosphere["rho"])
    plt.title("atmosphere")
    

# affichage
x = objet["fusee"]["x"]
y = objet["fusee"]["y"]

plt.figure(1)
plt.plot(objet["fusee"]["x"],objet["fusee"]["y"],label="fusee")
plt.plot(R_x,R_y,label="Terre")
plt.plot(Rh_x,Rh_y,label="Terre+400km")
plt.plot(corps["Lune"]["x"],corps["Lune"]["y"],label="Lune")
plt.legend()


def animate(k):
    i = min(k, len(x))
    line.set_data(x[:i], y[:i])
    point.set_data(x[i], y[i])
    Lune_p.set_data(corps["Lune"]["x"][i],corps["Lune"]["y"][i])
    Lune_l.set_data(corps["Lune"]["x"][:i],corps["Lune"]["y"][:i])
    return line, point, Lune_p,Lune_l

## Création de la figure et de l'axe
fig, ax = plt.subplots()

## Création de la ligne qui sera mise à jour au fur et à mesure
planete, = ax.plot(R_x,R_y, color='blue')
altitude, = ax.plot(Rh_x,Rh_y, color='red')
line, = ax.plot([],[], color='green')
point, = ax.plot([], [], ls="none", marker="o")
Lune_p,=ax.plot([], [], ls="none", marker="o")
Lune_l, = ax.plot([],[], color='orange')

#Gestion des limites de la fenêtre
lim_m = max([abs(np.min([np.min(x),min(y)])),np.max([max(y),max(x)])])
ax.set_xlim([-lim_m,lim_m])
ax.set_ylim([-lim_m,lim_m])
plt.title('mouvement Lune et fusée')


ani = animation.FuncAnimation(fig=fig, func=animate, frames=range(len(x)), interval=1, blit=True)

plt.show()
