# -*- coding: utf-8 -*-
"""
Created on Fri Aug 08 20:54:30 2014

@author: Carlos R. Pascual
"""

############################### LIBRERÍAS ###############################

import numpy as np


############################### FUNCIONES ###############################

def glassRefractionIndex(wavelen):
    '''
    Devuelve el índice de refracción del vídrio en función de la longitud
    de onda incidente sobre el.
    
    wavelen: Longitud de onda (micrómetros).
    '''
    return 1.5
    
def glassReflectivity(incidenceAngle, wavelen, airToGlass):
    '''
    Determina mediante la ley de Fresnel la reflectividad del vídrio en 
    función del ángulo de incidencia y la longitud de onda del fotón, cuando 
    uno de los medios es aire atmosférico y el otro vídrio.
    
    incidenceAngle: Angulo de incidencia (float, rad).
    wavelen: Longitud de onda incidente (float, um).
    airToGlass: True si el fotón pasa de aire a cristal, falso si pasa de
                cristal a aire (boolean).
    '''
    # Determinar el sentido de propagación
    if airToGlass:
        n1 = 1.0003
        n2 = glassRefractionIndex(wavelen)
    else:
        n1 = glassRefractionIndex(wavelen)
        n2 = 1.0003
    # Calcular la reflectividad
    cosTheta1 = np.cos(incidenceAngle)
    cosTheta2 = np.cos(np.arcsin((n1/n2) * np.sin(incidenceAngle)))    
    f1 = ((n1*cosTheta1 - n2*cosTheta2) / (n1*cosTheta1 + n2*cosTheta2))**2
    f2 = ((n1*cosTheta2 - n2*cosTheta1) / (n1*cosTheta2 + n2*cosTheta1))**2
    return (f1 + f2) * 0.5
    
def coatReflectivity(wavelen):
    '''
    Devuelve la reflectividad del recubrimiento interno del heliostato en
    función de la longitud de onda incidente.
    
    wavelen: Longitud de onda incidente (micrómetros).
    '''
    return 0.92
    
def airAbsortivity(wavelen, distance):
    '''
    Devuelve la absortividad de una capa de aire de espesor conocido.
    
    wavelen: Longitud de onda (micrómetros)
    distance: Distancia viajada por la radiación (metros)
    '''
    k = 1.005e-5
    return 1 - np.exp(-k*distance)    
    
def collectorAbsortivity(wavelen):
    '''
    Devuelve la reflectividad de la superficie colectora de la caldera solar
    en función de la longitud de onda incidente.
    
    wavelen: Longitud de onda incidente (micrómetros).
    '''
    return 0.9


############################### CLASES ###############################

class Heliostat(object):
    
    def __init__(self, coordinates, targetPoint, solarVector, position,
                 height=5, wide=6.3, large=6.3):
        '''
        coordinates: Coordenadas XY del heliostato (tupla, metros).
        targetPoint: Punto donde apunta el heliostato (array, metros).
        solarVector: Dirección principal de los rayos (array).
        position: Fila y lugar en la fila (tupla).
        height: altura de la estructura del heliostato (float, metros).
        wide: Ancho de la superficie reflectante del heliostato (float, metros).
        large: Alto de la superficie reflectante del heliostato (float, metros).
        '''
        self.coordinates = np.array([coordinates[0], coordinates[1], height], 
                                    float)
        self.wide = float(wide)
        self.large = float(large)
        self.solarVector = solarVector
        self.position = position
        self.targetDistance = np.linalg.norm(targetPoint - self.coordinates)
        modulus = np.linalg.norm(targetPoint - self.coordinates)
        self.targetVector = (targetPoint - self.coordinates) / modulus
        self.incidenceAngle = self.init_incidenceAngle()
        self.orientationAngles = self.init_orientationAngles()
        
    def init_incidenceAngle(self):
        '''
        Determina en ángulo de incidencia de la radiación solar con un
        heliostato orientado hacia la caldera solar.
        '''
        modTargetVec = np.linalg.norm(self.targetVector)
        modSunVec = np.linalg.norm(self.solarVector)
        dotProduct = np.dot(self.solarVector, self.targetVector)
        return 0.5 * np.arccos(dotProduct / (modTargetVec*modSunVec))
        
    def init_orientationAngles(self):
        '''
        Determina el ángulo de inclinación del heliostato respecto de la 
        horizontal (suelo) y el ángulo de orientación respecto a la dirección
        sur.
        '''
        # Cálculo de alfa, inclinación del heliostato respecto al suelo.
        v1 = np.cross(-1*self.solarVector, self.targetVector)
        v2 = self.targetVector + -1*self.solarVector
        helioVector = np.cross(v1, v2) / np.linalg.norm(np.cross(v1, v2))
        groundVector = np.array([0, 0, 1])
        alpha = np.arccos(np.dot(helioVector, groundVector))
        # Cálculo de beta, orientación del heliostato respecto al sur.
#        v1, modV1 = self.targetVector[:2], np.linalg.norm(v1)
#        v2, modV2 = np.array([0, -1]), np.linalg.norm(v2)
#        beta = 0.5 * np.arccos(np.dot(v1, v2) / (modV1*modV2))
        return (alpha,)
        
    def getCoordinates(self):
        '''
        Devuelve las coordenadas XYZ de la superficie reflectante del 
        heliostato.
        '''
        return self.coordinates     
        
    def getSurface(self):
        '''
        Devuele la superficie reflectante que tiene el heliostato.
        '''
        return self.wide * self.large
        
    def getPosition(self):
        '''
        Devuelve la posición del heliostato en el campo. (a, b):
        a: Fila en la que se encuentra.
        b: Posición en la fila.
        '''
        return self.position
        
    def getOrientationAngles(self):
        '''
        Devuelve los ángulos de horientación respecto del suelo y la dirección
        sur.
        '''
        return self.orientationAngles
        
    def cosineEffect(self):
        '''
        Aplica el efecto coseno retirando una parte proporcional de los
        fotones dirigidos al heliostato.
        
        photons: Array 4 x n generado en la clase 'Sun'.
        '''
        return np.cos(self.incidenceAngle)

    def mirrorSurvivors(self, photons):
        '''
        Determina los elementos del array de fotones que sobreviven a las 
        interacciones con el heliostato.
        
        photons: array de fotones creado en la clase 'Sun' por el método
                 'generatePhotons'.
        '''
        s = []
        for i in xrange(5):
            odds = np.random.random(len(photons[0]))
            if i == 0:
                rhoGlass1 = glassReflectivity(self.incidenceAngle, photons[0], 
                                              True)
                s1 = photons[:, odds <= rhoGlass1]
                s.append(s1)
                photons = photons[:, odds < (1 - rhoGlass1)]
            elif i == 2 or i == 4:
                rhoGlass2 = glassReflectivity(self.incidenceAngle, photons[0], 
                                              False)
                s2 = photons[:, odds <= (1 - rhoGlass2)]
                s.append(s2)
                photons = photons[:, odds < rhoGlass2]
            else:
                rhoCoat = coatReflectivity(photons[0])
                photons = photons[:, odds < rhoCoat]
        return np.hstack(s)
        
    def airSurvivors(self, photons):
        '''
        Determina los elementos del array de fotones que sobreviven a las 
        interacciones con el aire atmosférico.
        
        photons: array de fotones creado en la clase 'Sun' por el método
                 'generatePhotons'.
        '''
        odds = np.random.random(len(photons[0]))
        alphaAtmosphere = airAbsortivity(photons[0], self.targetDistance)
        photons = photons[:, odds < (1 - alphaAtmosphere)]
        return photons
        
    def collectorSurvivors(self,photons):
        '''
        Determina los elementos del array de fotones que sobreviven a las 
        interacciones con la superficie colectora de la caldera solar.
        
        photons: array de fotones creado en la clase 'Sun' por el método
                 'generatePhotons'.
        '''
        odds = np.random.random(len(photons[0]))
        alphaCollector = collectorAbsortivity(photons[0])
        photons = photons[:, odds < alphaCollector]
        return photons
