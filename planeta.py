#!/usr/bin/env python
# -*- coding: utf-8 -*-
G=1
M=1
m=1

class Planeta(object):
    '''
    Complete el docstring.
    '''

    def __init__(self, condicion_inicial, alpha=0):
        '''
        __init__ es un método especial que se usa para inicializar las
        instancias de una clase.

        Ej. de uso:
        >> mercurio = Planeta([x0, y0, vx0, vy0])
        >> print(mercurio.alpha)
        >> 0.
        '''
        self.y_actual=condicion_inicial
        self.t_actual=0.
        self.alpha=alpha

    def ecuacion_de_movimiento(self):
        '''
        Implementa la ecuación de movimiento, como sistema de ecuaciónes de
        primer orden.
        '''
        x,y,vx,vy=self.y_actual
        fx=(2*self.alpha*G*M*x)/((x**2+y**2)**2)-(G*M*x)/((x**2+y**2)**(3/2.))
        fy=(2*self.alpha*G*M*y)/((x**2+y**2)**2)-(G*M*y)/((x**2+y**2)**(3/2.))
        return [vx, vy, fx, fy]

    def avanza_euler(self, dt):
        '''
        Toma la condición actual del planeta y avanza su posicion y velocidad
        en un intervalo de tiempo dt usando el método de Euler explícito. El
        método no retorna nada, pero re-setea los valores de self.y_actual.
        '''
        vx,vy,fx,fy=self.ecuacion_de_movimiento() #valores que voy a ocupar

        #método de Euler
        x=vx+dt*self.y_actual[0]
        y=vy+dt*self.y_actual[1]
        v_x=fx+dt*self.y_actual[2]
        v_y=fy+dt*self.y_actual[3]

        now=[x,y,v_x,v_y] #arrejunto los nuevos valores en un arreglín
        
        #redefino la posición en el nuevo tiempo:
        self.y_actual=now
        self.t_actual+=dt
        pass

    def avanza_rk4(self, dt):
        '''
        Similar a avanza_euler, pero usando Runge-Kutta 4.
        '''
        x0,y0,vx0,vy0=self.y_actual

        #definición de k's
        k1=self.ecuacion_de_movimiento() #este va así tal cual
        self.y_actual=[x0+dt*k1[0]/2.,y0+dt*k1[1]/2.,vx0+dt*k1[2]/2.,vy0+dt*k1[3]/2.] #redefino para k2
        k2=self.ecuacion_de_movimiento()
        self.y_actual=[x0+dt*k2[0]/2.,y0+dt*k2[1]/2.,vx0+dt*k2[2]/2.,vy0+dt*k2[3]/2.] #redefino para k3
        k3=self.ecuacion_de_movimiento()
        self.y_actual=[x0+dt*k3[0],y0+dt*k3[1],vx0+dt*k3[2],vy0+dt*k3[3]] #redefino para k4
        k4=self.ecuacion_de_movimiento()

        #paso de runge-kutta
        x=x0+(dt/6.)*(k1[0]+2*k2[0]+2*k3[0]+k4[0])
        y=y0+(dt/6.)*(k1[1]+2*k2[1]+2*k3[1]+k4[1])
        vx=vx0+(dt/6.)*(k1[2]+2*k2[2]+2*k3[2]+k4[2])
        vy=vy0+(dt/6.)*(k1[3]+2*k2[3]+2*k3[3]+k4[3])

        now=[x,y,vx,vy] #arrejunto los nuevos valores en un arreglín

        #redefino la posición en el nuevo tiempo:
        self.y_actual=now
        self.t_actual+=dt

        pass

    def avanza_verlet(self, dt):
        '''
        Similar a avanza_euler, pero usando Verlet.
        '''
        pass

    def energia_total(self):
        '''
        Calcula la energía total del sistema en las condiciones actuales.
        '''
        pass
