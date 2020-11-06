from matplotlib.pylab import *
import matplotlib.pyplot as plt
import numpy as np
import math
L = 1.04
n = 100
dx = L / n


x = [0.104, 0.208, 0.416] 


K = 1.495    # W /m C
c = 1023     # J / kg C
ρ = 2476     # kg / m3

α = K/(c*ρ)




t_final = 100800   #28 horas







for xs in x:
    
    plt.figure()
    
    
    # FOURIER
    tf = [[],[]]
    for t in range(0,t_final,100):
    
        suma = 0
        for n in range(1,500): 
            fourier =  ((40*( 1 - (-1)**n )) / (n*np.pi)) * np.sin( (n*np.pi*xs)/(L)) * np.exp(-α * t * (((n*np.pi)/(L))**2) )
            suma += fourier
            
        tf[0].append(t)
        tf[1].append(suma)
    
    plt.plot(array(tf[0])/3600,tf[1],"--k",label="Serie de Fourier")
    
    
    
    
    
    # ECUACION DE CALOR CON dt = 1, 5, 10, 50 y 100
    
    
    steps = [1,5,10,50,100]
   
    # for s in steps:
        
        
    #     # PARAMETROS INICIALES 
    #     t_final = 100800
    #     L = 1.04
       
    #     # CANTIDAD DE COLUMNAS
    #     n = 1000                 
       
    #     dx = L / n 
        
    #     α = (K/(c*ρ)) * (s/(dx**2))


    #     tiempos = []
    #     for t in range(0,t_final,s):
    #         tiempos.append(t)
    #     #     print(t)
    #     # print (tiempos)
        
    #     Nt = (len(tiempos))


    #     u = zeros(( int(t_final/s) ,n+1))
        
        
    #     # CONDICIONES
    #     u[:,0]    = 0
    #     u[:,-1]   = 0
    #     u[0,1:n]   = 20
        
    #     print (u)
        
        
        
    """
    
    IMPORTANTE: ESTA PARTE DE AQUI NO NOS SALIA PORQUE TENIAMOS UN ERROR DE
    QUE EL RESULTADO DE LOS CALCULOS ERA "NOT A NUMBER" , SIN EMBARGO EL ALGORITMO ESTABA BUENO
    Y NO PUDIMOS ENCONTRAR EL ERROR, LO DEJAMOS COMO ESTA AHORA LUEGO DE VARIAS PRUEBAS
        
    """
        
    #     for c,t in enumerate(tiempos): 
            
    #         for i in range(1,n):
                

    #             # print (c,i)
    #             if c < 1007:
                    
                    
    #                 u[c+1, i] = round(( u[c, i] +  ( α * ( (u[c, i-1]) - (2*u[c, i]) + (u[c, i+1])) )),1)
            
    #     print (u)
        
    #     plt.plot(array(tiempos)/3600,u[:,104])
                
        

            
     
       
            
          
          
            
    # plt.ylim(0,25)
    plt.legend()
    plt.grid()
    # plt.savefig("Imagen",dpi=300)






























