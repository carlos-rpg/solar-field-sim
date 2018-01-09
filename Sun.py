# -*- coding: utf-8 -*-
"""
Created on Wed Jul 30 13:29:28 2014

@author: Carlos R. Pascual
"""

############################### LIBRERÍAS ###############################

import numpy as np
from scipy import integrate, constants


############################### CONSTANTES ###############################

h = constants.h # Constante de Planck en J·s
c = constants.c # Velocidad de la luz en m/s
k = constants.k # Constante de Boltzmann
sigma = constants.sigma # Constante de Stefan-Boltzmann
T_sol = 5778 # Temperatura del sol en K (Cuerpo negro)


############################### FUNCIONES ###############################

def logNormal(x, mu, sigma):
    '''
    Devuelve la probabilidad de que ocurra la variable x con una 
    distribución log-normal de parámetros de forma mu y sigma.
	
    x: Variable aleatoria (array, float)
    mu: Factor de forma (float).
    sigma: Factor de forma (float).
    '''
    f1 = np.exp(-(np.log(x) - mu)**2 / (2*sigma**2))
    f2 = x*sigma*(2*np.pi)**0.5
    return f1/f2

def arrangePhotons(photons, numberOfBins=100):
    '''
    Ordena el array photons en una lista de numberOfBins objetos.
    
    photons: longitudes de onda de los fotones (array, float)
    numberOfBins: número de intervalos en que se van a dividir los elementos
                  de photons (integer).
    '''
    bins = []
    limits = np.arange(0, numberOfBins + 1) / (numberOfBins/10.0)
    for i in xrange(len(limits) - 1):
        b1n = photons[(photons > limits[i]) & (photons < limits[i+1])]
        bins.append(b1n)
    return (limits, bins)

def averageEnergy(wavelen_1, wavelen_2):
    '''
    Calcula el valor medio de la energía de un fotón en un intervalo
    de longitudes de onda (J). 
    
    wavelen_1: Longitud de onda, cota inferior del intervalo (float, um).
    wavelen_2: Longitud de onda, cota superior del intervalo (float, um).
    '''
    if wavelen_1 == 0.0: # Evita error por dividir entre cero.
        wavelen_1 = 0.0001
    wavelen_1 *= 1e-6
    wavelen_2 *= 1e-6
    return (h*c / (wavelen_2-wavelen_1)) * np.log(wavelen_2 / wavelen_1)
    
def solarUnitVector(solarAngles):
    '''
    Vector unitario que señala la posición del sol en el cielo.
    
    solarAngles: Array compuesto por [gamma, psi] en radianes.
    '''
    solarHeight = solarAngles[0]
    solarAzimut = solarAngles[1]
    x = -1 * np.cos(solarHeight) * np.sin(solarAzimut)
    y = -1 * np.cos(solarHeight) * np.cos(solarAzimut)
    z = np.sin(solarHeight)
    return np.array([x, y, z])
    
def deltaAngles(n):
    '''
    Genera n incertidumbres sobre la altura y el acimut calculados por
    el método angles, una por cada fotón a modelar. 
    
    n: Número de fotones a modelar (integer).
    '''
    radiuses = np.random.triangular(0, 0, 4.65e-3, n)
    angles = np.random.random(n) * 2*np.pi
    deltaPsi = radiuses * np.cos(angles)
    deltaGamma = radiuses * np.sin(angles)
    return np.array([deltaGamma, deltaPsi])


############################### CLASES ###############################

class Sun(object):
    
    def __init__(self, B, latitude, day, solarTime=12):
        '''
        B: Irradiancia normal directa (W/m^2).
        n: Número de fotones a modelar.
        latitude: Latitud del emplazamiento (grados).
        day: Día del año como entero entre 1 y 365.
        solarTime: Hora solar como número entre 0.0 y 23.0.
        '''
        self.B = float(B)
        self.latitude = latitude * np.pi/180
        self.day = int(day)
        self.solarTime = solarTime
        self.epsilon = float(B) / (sigma*T_sol**4)
        self.N = integrate.quad(self.N_l, 0, np.inf)[0]
        self.alphaFactors = []
        self.solarAngles = self.init_solarAngles()
        self.solarVector = solarUnitVector(self.solarAngles)
        
        if self.solarAngles[0] < 0.0:
            raise ValueError(str(self.solarTime) + ' hours is night time')
    
    def init_solarAngles(self):
        '''
        Altura y acimut solares [rad]. 
        '''
        # Calcula la declinación y el ángulo solar.
        delta = 0.4093*np.sin((360*(284 + self.day)/365.0) * np.pi/180)
        omega = 15*(self.solarTime - 12) * np.pi/180    
        # Cálculo de la altura solar
        f1 = np.cos(self.latitude) * np.cos(delta) * np.cos(omega)
        f2 = np.sin(self.latitude) * np.sin(delta)
        gamma = np.arcsin(f1 + f2) 
        # Cálculo del acimut solar
        f1 = np.sin(gamma) * np.sin(self.latitude) - np.sin(delta)
        f2 = np.cos(gamma) * np.cos(self.latitude)
        psi = np.arccos(f1 / f2)     
        # Criterio de signos para el acimut
        if self.solarTime < 12:
            return np.array([[gamma], [-psi]])
        return np.array([[gamma], [psi]])
    
    def getSolarAngles(self):
        '''
        Devuelve los ángulos solares gamma y psi.
        '''
        return self.solarAngles
        
    def getSolarVector(self):
        '''
        Devuelve el vector unitario que señala la posición del sol en el
        cielo.
        '''
        return self.solarVector
        
    def N_l(self, wavelen):
        '''
        Nº de fotones por ud. de tiempo, área y longitud de onda para una
        irradiancia B (fotones/s·m^2·um)
        
        wavelen: longitud de onda (micrómetros)
        '''
        wavelen *= 1e-6
        f1 = 2*np.pi*c/wavelen**4
        f2 = np.exp(h*c/k/wavelen/T_sol) - 1
        return self.epsilon*f1/f2
        
    def wavelenPDF(self, wavelen):
        '''
        Función densidad de probabilidad de cada longitud de onda emitida por 
        el sol.
        
        wavelen: longitud de onda (micrómetros)
        '''
        return self.N_l(wavelen) / self.N
                
    def calculateAlphaFactors(self, photons):
        '''
        Calcula los factores de proporcionalidad entre la cantidad de fotones
        real emitida por el sol en un intervalo de longitudes de onda y la
        cantidad de fotones a modelar.
        
        photons: Longitudes de onda de los fotones a simular (array, floats).
        '''
        alphaFactors = []
        limits, bins = arrangePhotons(photons)
        for i in xrange(len(bins) - 1):
            # N_i: Número de fotones en el bin
            N_i = len(bins[i]) 
            if N_i == 0:
                # alpha=0 si el bin está vacío
                alphaFactors.append(0.0)
                continue
            # N: Cantidad teórica de fotones en el bin.
            N = integrate.quad(self.N_l, limits[i], limits[i+1])[0]
            alphaFactors.append(N / N_i)
        return alphaFactors
        
    def generatePhotons(self, numberOfPhotons):
        '''
        Genera un array con la información de cada fotón creado.
        fila 0 -> Longitudes de onda (float, um).
        fila 1 -> Coordenada X del vector de dirección.
        fila 2 -> Coordenada Y del vector de dirección.
        fila 3 -> Coordenada Z del vector de dirección.
        
        numberOfPhotons: nº de fotones a modelar (integer).
        '''
        # Garantiza que n es un número entero.
        n = int(numberOfPhotons)
        # Genera n fotones siguiendo la distribución lognormal
        photonWavelengths = np.random.lognormal(-0.04243, 0.6398, n)
        # Factores alpha que corresponden a los fotones generados
        self.alphaFactors = self.calculateAlphaFactors(photonWavelengths)
        # Dirección de orígen para cada fotón generado.
        solarAngles = self.solarAngles + deltaAngles(n)
        return np.vstack([photonWavelengths, solarUnitVector(solarAngles)])
 
    def photonsToEnergy(self, photons):
        '''
        Dado un array de longitudes de onda, determina la irrandiancia que 
        le corresponde a dicho array en W/m^2.
        
        photons: Longitudes de onda (array, floats)
        '''
        limits, bins = arrangePhotons(photons[0])
        B = 0
        for i in xrange(len(bins) - 1):
            N_i = len(bins[i])
            if N_i == 0:
                continue
            factor = self.alphaFactors[i]
            E = averageEnergy(limits[i], limits[i+1])
            B += N_i * E * factor
        return B * 1e-6 # Conversión de unidades
