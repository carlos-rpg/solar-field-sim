# -*- coding: utf-8 -*-
"""
Created on Thu Aug 21 17:39:06 2014

@author: Carlos
"""

############################# LIBRERÍAS ##############################

import numpy as np
import Sun as sn
import Field as fd

#np.seterr(over='ignore') # Ignora una advertencia inócua de overflow.


############################# INPUTS ##################################

B = 800 # Irradiancia de diseño (800 W/m^2)
latitude = 40 # latitud geográfica (40º Norte)
day = 172 # Día del año (172 solsticio de verano)
baseHeight = 10 # Base de la caldera solar (10 metros)
apertureAngle = 90 # Ángulo de apertura (90 grados)
firstRowDist = 40 # Distancia a la primera fila (80 metros)
numberOfPhotons = 10000 # Número de fotones para modelar (10.000)
boilerPower = 74e3 # Densidad de potencia del evaporador (74e3 W/m2)
overheaterPower = 52e3 # Densidad de potencia del sobrecalentador (52e3 W/m2)
surfaceSection = 5*5 # Superficie de una sección objetivo (m2)


############################ FUNCIONES ################################
    
def designField(targets, powerRequired):
    '''
    Diseña un campo solar para que se cumplan el requisito de potencia de la
    caldera solar para cada una de las divisiones de su pared colectora.
    
    targets: array con los objetivos a enfocar.
    powerRequired: Potencia necesaria en una sección de 5x5 m de la pared.
    '''
    # Para cada sección objetivo de la pared colectora:
    for target in targets:
        powerGathered = 0
        # Mientras no se cumpla el requisito de potencia:
        while powerGathered < powerRequired:
            field.placeHeliostat(target, solarVector)
            newHeliostat = field.getHeliostats()[-1]
            coordinates = newHeliostat.getCoordinates()
            # Si el heliostato está bajo la sombra de la caldera:
            if field.isTowerShading(coordinates, solarAngles):
                surfacesInShadow.append(newHeliostat.getSurface())
                continue # Siguiente iteración     
            # Interacciones con el heliostato:
            photonsLeft = newHeliostat.mirrorSurvivors(photons)
            irradianceAfterMirror = sun.photonsToEnergy(photonsLeft)
            opticPerformances.append(irradianceAfterMirror / B)
            # Intercacciones con el aire atmosférico:
            photonsLeft = newHeliostat.airSurvivors(photonsLeft)
            irradianceAfterAir = sun.photonsToEnergy(photonsLeft)
            airPerformances.append(irradianceAfterAir / irradianceAfterMirror)
            # Interacciones con la superficie colectora:
            photonsLeft = newHeliostat.collectorSurvivors(photonsLeft)
            irradianceAfterCollector = sun.photonsToEnergy(photonsLeft)
            wallPerformances.append(irradianceAfterCollector / irradianceAfterAir)        
            # Efecto coseno:
            cosine = newHeliostat.cosineEffect()
            cosineEffects.append(1 - cosine)
            # Área efectiva de reflexión:
            area = newHeliostat.getSurface()
            shadow = field.surfaceShaded(newHeliostat, solarVector)
            surfacesInShadow.append(shadow)
            # Balance final:
            powerGathered += irradianceAfterCollector * (area*cosine - shadow) 
   
def printResults():
    '''
    Muestra en pantalla información del campo solar diseñado.
    '''
    # Cálculos previos:
    lastHeliostat = field.getHeliostats()[-1]
    x, y = lastHeliostat.getCoordinates()[:2]
    mirrorSurface = field.getHeliostatNumber() * lastHeliostat.getSurface()
    radius = (x**2 + y**2)**0.5
    landSurface = (np.pi * (radius**2 - firstRowDist**2)) * apertureAngle/360
    absShadedSurface = sum(surfacesInShadow)
    relShadedSurface = (absShadedSurface / mirrorSurface) * 100
    occupation = (mirrorSurface / landSurface) * 100
    averageOptical = np.average(opticPerformances) *100
    averageAir = np.average(airPerformances) *100
    averageWall = np.average(wallPerformances) *100
    averageCosine = np.average(cosineEffects) *100
    # Resultados:
    field.draw()
    print
    print '------------ RESULTS ------------'
    print 'Heliostats:', str(field.getHeliostatNumber())
    print 'Land surface:', str(round(landSurface * 1e-4, 2)), 'hec'
    print 'land occupation:', str(round(occupation, 2)), '%'
    print 'Total mirror surface in shade:', str(round(absShadedSurface, 2)), 'm2'
    print 'Relative mirror surface in shade:', str(round(relShadedSurface, 2)), '%'
    print 'Cosine effect:', str(round(averageCosine, 2)), '%'
    print
    print 'Optical performance:', str(round(averageOptical, 2)), '%'
    print 'Air performance:', str(round(averageAir, 2)), '%'
    print 'Wall performance:', str(round(averageWall, 2)), '%'
                   

########################### SIMULACIÓN ################################

# Objetos previos a la simulación:    
sun = sn.Sun(B, latitude, day)
field = fd.Field(baseHeight, apertureAngle, firstRowDist)
photons = sun.generatePhotons(numberOfPhotons)
solarAngles = sun.getSolarAngles()
solarVector = np.transpose(sun.getSolarVector())[0]

# Contadores:
surfacesInShadow = []
opticPerformances = []
airPerformances = []
wallPerformances = []
cosineEffects = []

# Alimentar al evaporador:
boilerTargets = field.getTargetCoords()[:8]
boilerPowerRequired = boilerPower * surfaceSection
designField(boilerTargets, boilerPowerRequired)
   
# Alimentar al sobrecalentador:
overheaterTargets = field.getTargetCoords()[8:]
overveaterPowerRequired = overheaterPower * surfaceSection
designField(overheaterTargets, overveaterPowerRequired)

# Resultados finales de diseño:
printResults()
