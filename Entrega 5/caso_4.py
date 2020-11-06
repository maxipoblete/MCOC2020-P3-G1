from matplotlib.pylab import *
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import math
import sys







"""
----------------- CASO GENERAL ----------------------

TIEMPO DE SIMULACION:   ~1 DIA 

SIMULACION CADA:        30 MINUTOS


#----- CONDICIONES INICIAL  -------

u_k[:,:] = T                        # Temperatura Inicial

#----- CONDICIONES DE BORDE -------

T: Temperatura
G: Gradiente

u_k[0 ,:]  = T                      # Borde Izquierdo
u_k[-1,:]  = T                      # Borde Derecho 

u_k[:, 0]  = T                      # Borde Inferior
u_k[:,-1]  = T                      # Borde Superior


u_k[-1,:]   =  u_k[-2,:] + G*dx     # Gradiente Derecho
u_k[0,:]    =  u_k[1,:]  - G*dx     # Gradiente Izquierdo

u_k[:,-1] =  u_k[:,-2] + G*dy       # Gradiente Superior
u_k[:,0]  =  u_k[:,1]  - G*dy       # Gradiente Inferior

"""





#----- ELEMENTOS INICIALES PARA DEFINIR EL PROBLEMA

a      = 1 
b      = 0.5

Nx     = 15
Ny     = 30
dx     = b/Nx
dy     = a/Ny
h      = dx
coords = lambda i,j:(dx*i,dy*j)
x,y    = coords(2,2)
u_k    = zeros((Nx+1,Ny+1),dtype=double)
u_km1  = zeros((Nx+1,Ny+1),dtype=double)






#----- CONDICION INICIAL
u_k[:,:] = 10 



#----- ALGUNOS PARAMETROS
dt        = 0.01 
K         = 79.5
c         = 450
p         = 7800
alfa      = K*dt/(c*p*dx**2)
minuto    = 60
hora      = 3600
dia       = 86400
dt        = 1      * minuto
d_next    = 0.5    * hora
next_t    = 0
framenum  = 0 
T         = 1      * dia
Days      = 1.005  * T           # ---> DIAS A SIMULAR




#------------------- Vectores de Temperatura en puntos de interes.

U_A2_B2   =  zeros(int32(Days/dt))
U_A2_3B4  =  zeros(int32(Days/dt))
U_3A4_3B4 =  zeros(int32(Days/dt))





#--------------------------------------------- FUNCION 1 
def truncate(n,decimals=0):
    multiplier = 10 ** decimals
    return int(n*multiplier) /multiplier



#--------------------------------------------- FUNCION 2

def imshowbien(u):
    cmapchoice = cm.Spectral_r
    imshow(u.T[Ny::-1,:],cmap=cmapchoice,interpolation="Gaussian")
    cbar           = colorbar(extend="both",cmap=cmapchoice)
    ticks          = arange(0,35,5)
    ticks_Text     = ["{}º".format(deg) for deg in ticks]
    cbar.set_ticks(ticks)
    cbar.set_ticklabels(ticks_Text)
    clim(0,30)
    xlabel("b")
    ylabel("a")
    xTicks_N       = arange(0,Nx+1,3)
    yTicks_N       = arange(0,Ny+1,3)    
    xTicks         = [coords(i,0)[0] for i in xTicks_N]
    yTicks         = [coords(0,i)[1] for i in yTicks_N]
    xTicks_Text    = ["{0:.2f}".format(tick) for tick in xTicks]
    yTicks_Text    = ["{0:.2f}".format(tick) for tick in yTicks]    
    yTicks_Text = yTicks_Text[::-1]
    xticks(xTicks_N,xTicks_Text,rotation="vertical")
    yticks(yTicks_N,yTicks_Text)
    margins(0.2)
    subplots_adjust(bottom=0.15)
    subplots_adjust(right=0.75)
    
    
    
    



#--------------------------------------------- ALGORITMO

#----- Loop Tiempo -----

for k in range(int32(Days/dt)):
    t       = dt*(k+1)
    dias    = truncate(t/dia,0)
    horas   = truncate((t-dias*dia)/hora,0)
    minutos = truncate((t-dias*dia-horas*hora)/minuto ,0)
    titulo  = "k = {0:05.0f}".format(k) + " t={0:02.0f}d {1:02.0f}h {2:02.0f}m ".format(dias,horas,minutos)
    print (titulo)
    
    
    
    #-------   CONDICIONES DE BORDE
    u_k[0,:]   =  20                   #borde izq
    u_k[:,0]   =  20                   #borde inferior
    u_k[-1,:]  =  u_k[-2,:] + -10*dx   #gradiente derecho
    u_k[:,-1]  =  0                    #borde superior  
    
    

    #----- Loop Espacio -----
    for i in range(1,Nx):
        for j in range(1,Ny):
            #LAPLACE
            nabla_u_k  = (u_k[i-1,j] + u_k[i+1,j] + u_k[i,j-1] + u_k[i,j+1] - 4*u_k[i,j]) / h**2
            #FORWARD EULER
            u_km1[i,j] = u_k[i,j] + (alfa * nabla_u_k)
    u_k = u_km1
    
    
    #-------   CONDICIONES DE BORDE
    u_k[0,:]   =  20                   #borde izq
    u_k[:,0]   =  20                   #borde inferior
    u_k[-1,:]  =  u_k[-2,:] + -10*dx   #gradiente derecho
    u_k[:,-1]  =  0                    #borde superior  
         
    

    
    #------  VECTORES DE TEMPERATURA EN PUNTOS DE INTERES       #  a      b

    U_A2_B2[k]   =  u_k[ int(Nx/2)   , int(Ny/2)      ]         #  a/2    b/2
    U_A2_3B4[k]  =  u_k[ int(Nx/2)   , int(3*Ny/4)    ]         #  a/2   3b/4
    U_3A4_3B4[k] =  u_k[ int(3*Nx/4) , int(3*Ny/4)    ]         # 3a/4   3b/4


 

    #GRAFICAR SIGUIENTE
    if t > next_t:
        figure(1)
        imshowbien(u_k)
        titulo  = "k = {0:05.0f}".format(k) + " t={0:02.0f}d {1:02.0f}h {2:02.0f}m ".format(dias,horas,minutos-1)
        title(titulo)
        savefig("Ejemplo_caso_4/frame_{0:04.0f}.png".format(framenum))
        framenum  +=1
        next_t    += d_next
        close(1)




plt.figure(2)
plt.plot(array(range(int32(Days/dt)))/60, U_A2_B2    ,color="salmon"       , linewidth=2, label="Punto (a/2, b/2) ")
plt.plot(array(range(int32(Days/dt)))/60, U_A2_3B4   ,color="yellowgreen"  , linewidth=2, label="Punto (a/2, 3b/4)")
plt.plot(array(range(int32(Days/dt)))/60, U_3A4_3B4  ,color="turquoise"    , linewidth=2, label="Punto (3a/4, 3b/4)")



plt.title ("Evolución de la Temperatura - CASO 4 ")
plt.ylabel(" Temperatura [ºC] ")
plt.xlabel(" Tiempo [hrs] ")
plt.legend(loc=1)
plt.grid()
plt.show()




