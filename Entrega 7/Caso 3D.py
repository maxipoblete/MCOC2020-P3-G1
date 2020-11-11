from matplotlib.pylab import *
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import math
import sys
from calor_de_hidratacion import Calor_de_hidratacion 
from graficar import sensor 




#----- ELEMENTOS INICIALES PARA DEFINIR EL PROBLEMA



fname   = 'caso_2_intemperie.csv'
sensores_predichos  = [u1,u2,u3,u4,u5,u6,u7,u8,u9,u10,u11,u12,u13,u14,u15]
sensores_observados = [ sensor(fname)[1][0],
                        sensor(fname)[1][1],
                        sensor(fname)[1][2],
                        sensor(fname)[1][3],
                        sensor(fname)[1][4],
                        sensor(fname)[1][5],
                        sensor(fname)[1][6],
                        sensor(fname)[1][7],
                        sensor(fname)[1][8],
                        sensor(fname)[1][9],
                        sensor(fname)[1][10],
                        sensor(fname)[1][11],
                        sensor(fname)[1][12],
                        sensor(fname)[1][13],
                        sensor(fname)[1][14]]



tiempo_sensores = sensor(fname)[0]
sensores = sensor(fname)[1]





















x1 = 1.04
y1 = 0.54
z1 = 0.5


Nx     = int(100*1.04)
Ny     = int(100*0.54)
Nz     = int(100*0.5)


dx     = x1/Nx
dy     = y1/Ny
dz     = z1/Nz

h      = dx

coords = lambda i,j,k:(dx*i,dy*j,dz*k)
x,y,z    = coords(2,2,2)

u_k    = zeros(  ( Nx+1, Ny+1, Nz+1 )   ,dtype=double)
u_km1  = zeros(  ( Nx+1, Ny+1, Nz+1 )   ,dtype=double)





#----- CONDICION INICIAL
u_k[:,:,:] = 20



#----- ALGUNOS PARAMETROS
dt1       = 0.01 
K         = 1.495
c         = 1023
p         = 2476
alfa      = K*dt1/(c*p*dx**2)
minuto    = 60
hora      = 3600
dia       = 86400
dt        = 30      * minuto
d_next    = 0.5    * hora           # ---> INTERVALO DE TIEMPO A SIMULAR
next_t    = 0
framenum  = 0 
T         = 1      * dia
Days      = 0.05     * T           # ---> DIAS A SIMULAR





#----------- Vectores de Temperatura en puntos de interes.

u1     =  zeros(int32(Days/dt))
u2     =  zeros(int32(Days/dt))
u3     =  zeros(int32(Days/dt))
u4     =  zeros(int32(Days/dt))
u5     =  zeros(int32(Days/dt))
u6     =  zeros(int32(Days/dt))
u7     =  zeros(int32(Days/dt))
u8     =  zeros(int32(Days/dt))
u9     =  zeros(int32(Days/dt))
u10    =  zeros(int32(Days/dt))
u11    =  zeros(int32(Days/dt))
u12    =  zeros(int32(Days/dt))
u13    =  zeros(int32(Days/dt))
u14    =  zeros(int32(Days/dt))
u15    =  zeros(int32(Days/dt))







#--------------------------------------------- FUNCION 1 
def truncate(n,decimals=0):
    multiplier = 10 ** decimals
    return int(n*multiplier) /multiplier



#--------------------------------------------- FUNCION 2

def imshowbien(u):
    cmapchoice = cm.Spectral_r
    imshow(u.T[Ny::-1,:],cmap=cmapchoice,interpolation="None")
    cbar           = colorbar(extend="both",cmap=cmapchoice)
    ticks          = arange(0,75,10)
    ticks_Text     = ["{}º".format(deg) for deg in ticks]
    cbar.set_ticks(ticks)
    cbar.set_ticklabels(ticks_Text)
    clim(0,70)
    xlabel("Eje X")
    ylabel("Eje Y")
    xTicks_N       = arange(0,Nx+1,6.5)
    yTicks_N       = arange(0,Ny+1,6)    
    xTicks         = [coords(i,0,0)[0] for i in xTicks_N]
    yTicks         = [coords(0,i,0)[1] for i in yTicks_N]
    xTicks_Text    = ["{0:.2f}".format(tick) for tick in xTicks]
    yTicks_Text    = ["{0:.2f}".format(tick) for tick in yTicks]    
    yTicks_Text = yTicks_Text[::-1]
    xticks(xTicks_N,xTicks_Text,rotation="vertical")
    yticks(yTicks_N,yTicks_Text)
    margins(0.2)
    subplots_adjust(bottom=0.15)
    grid(True)
    
    
    



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
    u_ambiente = 20 + 10*sin((2*np.pi/T)*t)
    # u_k[ : , 0 ,: ]   =  20                   # BORDE FRONTAL 
    # u_k[ : ,-1 ,: ]   =  0                    # BORDE POSTERIOR 
    
    # u_k[: , : , 0  ]  = 50                    # BORDE INFERIOR
    u_k[: , : , -1 ]  = u_ambiente              # BORDE SUPERIOR 
    
    # u_k[ 0 , :  ,: ]  =  50                    # BORDE IZQUIEROD 
    # u_k[-1 , :  ,: ]  =  50                    # BORDE DERECHO 


    u_k[ 0 , :  ,: ]   =  u_k[1,:,:]  - 0*dx        # Gradiente Izquierdo
    u_k[-1 , :  ,: ]   =  u_k[-2,:,:] + 0*dx        # Gradiente Derecho

    u_k[ : , 0 ,: ]    =  u_k[:,1,:]  - 0*dy        # Gradiente Frontal
    u_k[ : ,-1 ,: ]    =  u_k[:,-2,:] + 0*dy        # Gradiente Posterior

    u_k[ : , : ,0 ]    =  u_k[:,:,1]    - 0*dz      # Gradiente Inferior
    # u_k[ : , : ,-1 ]   =  u_k[:,:,-2]  + 0*dz     # Gradiente Superior
     
    
    
    

    #----- Loop Espacio -----
    for i in range(1,Nx):
        for j in range(1,Ny):
            for w in range(1,Nz):
                #LAPLACE
                nabla_u_k  = (  u_k[i-1,j,w] + u_k[i+1,j,w] + u_k[i,j-1,w] + u_k[i,j+1,w]  + u_k[i,j,w-1] + u_k[i,j,w+1]  - 6*u_k[i,j,w] ) 
                #FORWARD EULER
                u_km1[i,j,w] =  dt1*Calor_de_hidratacion(t) + u_k[i,j,w] + ( alfa * nabla_u_k ) 
        
    u_k = u_km1
    
    

    #-------   CONDICIONES DE BORDE
    u_ambiente = 20 + 10*sin((2*np.pi/T)*t)
    # u_k[ : , 0 ,: ]   =  20                   # BORDE FRONTAL 
    # u_k[ : ,-1 ,: ]   =  0                    # BORDE POSTERIOR 
    
    # u_k[: , : , 0  ]  = 50                    # BORDE INFERIOR
    u_k[: , : , -1 ]  = u_ambiente              # BORDE SUPERIOR 
    
    # u_k[ 0 , :  ,: ]  =  50                    # BORDE IZQUIEROD 
    # u_k[-1 , :  ,: ]  =  50                    # BORDE DERECHO 


    u_k[ 0 , :  ,: ]   =  u_k[1,:,:]  - 0*dx        # Gradiente Izquierdo
    u_k[-1 , :  ,: ]   =  u_k[-2,:,:] + 0*dx        # Gradiente Derecho

    u_k[ : , 0 ,: ]    =  u_k[:,1,:]  - 0*dy        # Gradiente Frontal
    u_k[ : ,-1 ,: ]    =  u_k[:,-2,:] + 0*dy        # Gradiente Posterior

    u_k[ : , : ,0 ]    =  u_k[:,:,1]    - 0*dz      # Gradiente Inferior
    # u_k[ : , : ,-1 ]   =  u_k[:,:,-2]  + 0*dz     # Gradiente Superior
    


 


    #------  VECTORES DE TEMPERATURA EN PUNTOS DE INTERES


    u1[k] =  u_k[ 77, 27,25 ]
    u2[k] =  u_k[ 52, 27,25 ]
    u3[k] =  u_k[ 52, 30,25 ]  
    u4[k] =  u_k[ 52, 27,47 ]    
    u5[k] =  u_k[ 30, 27,25 ]     
    u6[k] =  u_k[ 52, 21,25 ]    
    u7[k] =  u_k[ 52, 27,36 ]    
    u8[k] =  u_k[ 52, 27,14 ]  
    u9[k] =  u_k[ 52, 39,25 ]   
    u10[k] = u_k[ 101,27,25 ]
    u11[k] = u_k[ 52, 27,30 ]
    u12[k] = u_k[ 28, 27,25 ]
    u13[k] = u_k[ 52, 27,50 ] # Externo
    u14[k] = u_k[ 52,15,25  ]
    u15[k] = u_k[ 52, 27,50 ] # Externo




    #GRAFICAR SIGUIENTE
    if t > next_t:
        
        figure(1)
        imshowbien( u_k[:,:,25] )   # ******
        titulo  = "k = {0:05.0f}".format(k) + " t={0:02.0f}d {1:02.0f}h {2:02.0f}m ".format(dias,horas,minutos-1)
        title(titulo)
        savefig("Ejemplo/frame_{0:04.0f}.png".format(framenum))
        framenum  +=1
        next_t    += d_next
        close(1)




fname   = 'caso_2_intemperie.csv'
sensores_predichos  = [u1,u2,u3,u4,u5,u6,u7,u8,u9,u10,u11,u12,u13,u14,u15]
sensores_observados = [ sensor(fname)[1][0],
                        sensor(fname)[1][1],
                        sensor(fname)[1][2],
                        sensor(fname)[1][3],
                        sensor(fname)[1][4],
                        sensor(fname)[1][5],
                        sensor(fname)[1][6],
                        sensor(fname)[1][7],
                        sensor(fname)[1][8],
                        sensor(fname)[1][9],
                        sensor(fname)[1][10],
                        sensor(fname)[1][11],
                        sensor(fname)[1][12],
                        sensor(fname)[1][13],
                        sensor(fname)[1][14]]



tiempo_sensores = sensor(fname)[0]

sensores = sensor(fname)[1]



for i in range(0,15):
    figure (i+2)
    plot( array(tiempo_sensores)/(3600*24)      , sensores_observados[i] ,label=f"Sensor Observado {i+1}" )
    plot( array(range(int32(Days/dt)))/(60*24)  , sensores_predichos[i]  ,label=f"Sensor Predicho {i+1}" )
    plt.title (f"Evolución de la Temperatura - Sensor {i+1}")
    plt.ylabel("Temperatura [ºC] ")
    plt.xlabel("Tiempo     [dias]    ")
    plt.legend()
    plt.grid()
    plt.show()
    
    








