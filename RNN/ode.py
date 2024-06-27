# -*- coding: utf-8 -*-
"""pr4
"""

"""## **Ejercicio 2)** Integrador de **Euler**
"""

def Euler(f, x0, t0, h, p): 
    return x0 + h*f(x0, t0, p)


def integrador_ode(m,f,xa,a,b,k,p):  #, c=lambda x,t,p:x):
    assert k>0
    n = len(xa)
    h = (b-a)/k
    w = np.zeros((n,k+1)) # Produce un array con forma y tipo especificada con los parametros, 
                          # lleno de ceros. la forma puede ser espcificada con un entero o tupla (n,k+1)    
    t = np.zeros(k+1)
    w[:,0] = xa           # actualiza la posicion inicial (columna de indice 0) de las variables con los valores 
                          # de las condiciones iniciales
    t[0] = a              # actualiza la posicion cero con el valor del tiempo inicial
    
    for j in range(k):    #Aca se produce la iteración en j 
        
        t[j+1] = t[j] + h                # iteracion tiempo 
        w[:,j+1] = m(f,w[:,j],t[j],h,p)  # iteracion de w 
        #w[:,j+1] = c(w[:,j+1],t[j+1],p)  # condicion sobre w

    return t,w


"""## **Ejercicio 7)** Integrador de **RK4**
"""

#1er paso metodo Runge-Kutta 4"
def RK4(f,t0,y0,h,p):
    k1 = f(t0, y0, p)
    k2 = f(t0+h/2, y0+h/2*k1, p)
    k3 = f(t0+h/2, y0+h/2*k2, p)
    k4 = f(t0+h, y0+h*k3, p)
    return y0 + (h/6)*(k1 + 2*k2 + 2*k3 + k4)      ### y0 + h/6(k1 + 2*k2 + 2*k3 + k4) No anda

def iteracion_ODE_multidimencional(Metodo,f,y0,a,b,N,p):
    t = np.zeros(N+1)
    w = np.zeros((len(y0),N+1))
    h = (b-a)/(N)
    t[0] = a
    w[:,0] = y0
    for i in range(1,N+1):
        t[i] = t[i-1]+h
        w[:,i] = Metodo(f,t[i-1],w[:,i-1],h,p)
    return t[:],w[:,:]


