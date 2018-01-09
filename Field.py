# -*- coding: utf-8 -*-
"""
Created on Tue Aug 05 19:40:41 2014

@author: Carlos R. Pascual
"""

############################### LIBRERÍAS ###############################

import numpy as np
import matplotlib.pyplot as plt
import Heliostat as ht


############################### FUNCIONES ###############################

def targetCoords(baseHeight):
    '''
    Para una superficie colectora dividida en 12 cuadros de 5 x 5 m, determina
    las coordenadas de los centros de cada cuadro en función de la altura de
    la base de la caldera solar.
    
    baseHeight: Altura desde el suelo hasta la base de la superficie colectora 
                (metros).
    '''
    x = np.ones(12) * 2.5
    for i in xrange(len(x)):
        if i % 2 == 1:
            x[i] *= -1
    y = np.zeros(12)
    z = np.array([2.5, 2.5, 7.5, 7.5, 12.5, 12.5, 17.5, 17.5, 
                  22.5, 22.5, 27.5, 27.5]) + baseHeight
    return np.transpose(np.vstack([x, y, z]))


############################### CLASES ###############################

class Field(object):
    
    def __init__(self, baseHeight, apertureAngle, firstRowDist,
                 heightHelio=5, wideHelio=6.3, largeHelio=6.3):
        '''
        baseHeight: Altura desde el suelo hasta la base de la superficie
                    colectora (float, metros).
        apertureAngle: Ángulo de apertura del campo solar (float, grados).
        firstRowDist: Distancia desde el orígen hasta la primera fila de
                      heliostatos (float, metros).
        heightHelio: Altura de la estructura de sujeción del heliostato 
                     (float, metros).
        wideHelio: Ancho de la superficie reflectante del heliostato 
                   (float, metros).
        largeHelio: Alto de la superficie reflectante del heliostato 
                    (float, metros).
        '''
        self.baseHeight = float(baseHeight)
        self.apertureAngle = apertureAngle * np.pi/180
        self.heightHelio = float(heightHelio)
        self.wideHelio = float(wideHelio)
        self.medDiameter = (float(wideHelio)**2 + float(largeHelio)**2)**0.5
        self.firstRowDist = float(firstRowDist)
        self.heliostats = []
        self.targetCoords = targetCoords(baseHeight)
               
    def getHeliostatNumber(self):
        '''
        Devuelve el nº de heliostatos que tiene el campo solar.
        '''
        return len(self.heliostats)
        
    def getTargetCoords(self):
        '''
        Devuelve las coordenadas del objetivo al que enfoca el heliostato.
        '''
        return self.targetCoords
        
    def getHeliostats(self):
        '''
        Devuelve la lista con los heliostatos del campo.
        '''
        return self.heliostats
        
    def deltaRadDist(self, radialDist, targetHeight, factor=0.25):
        '''
        Distancia radial entre la fila n y la n+2. 

        radialDist: Distancia de la fila de heliostatos a la caldera solar 
                    (metros).
        targetHeight: Altura a la que se encuentra el objetivo a enfocar 
                      (metros).
        factor: factor de corrección (float).
        '''
        theta = np.arctan((targetHeight - self.heightHelio*0.5) / radialDist)
        deltaRadDist = self.medDiameter * ((1.009/theta) - 0.063 + 0.04803*theta)
        return  deltaRadDist * factor
        
    def deltaAngDist(self, radialDist, targetHeight, factor=0.45):
        '''
        Distancia angular entre heliostatos de la misma fila. 

        radialDist: Distancia de la fila de heliostatos a la caldera solar 
                    (metros).
        targetHeight: Altura a la que se encuentra el objetivo a enfocar 
                      (metros).
        '''
        theta = np.arctan((targetHeight - self.heightHelio*0.5) / radialDist)
        deltaAngDist = self.medDiameter * (2.170 - 0.6589*theta + 1.247*theta**2)
        return deltaAngDist * factor
        
    def maxHeliosInRow(self, radialDist, targetHeight):
        '''
        Determina cuantos heliostatos caben en una fila situada a una 
        distancia radialDist de la caldera. 

        radialDist: Distancia de la fila a la caldera (metros).
        targetHeight: Altura a la que se encuentra el objetivo a enfocar 
                      (metros).
        '''
        arcLength = self.apertureAngle * radialDist
        angularSeparation = self.deltaAngDist(radialDist, targetHeight)
        return int(arcLength // angularSeparation)
        
    def assignCartesianCoords(self, rowDist, targetHeight, position, 
                              maxHeliosInRow):
        '''
        Asigna coordenadas cartesianas a un heliostato. 
        
        radialDist: Distancia a la fila donde se encuentra el heliostato (metros).
        targetHeight: Altura a la que se encuentra el objetivo a enfocar 
                      (metros).
        maxHeliosInRow: nº máximo de heliostatos que cabe en la fila.
        '''
        heliosDist = self.deltaAngDist(rowDist, targetHeight)
        # Si el numero max de heliostatos en la fila es impar:
        if maxHeliosInRow % 2 == 1:
            if position[1] % 2 == 1:
                n = (position[1] - 1) * 0.5
                angle = (heliosDist / rowDist) * n
            else:
                n = position[1] * 0.5
                angle = -(heliosDist / rowDist) * n
        # Si el número máximo no es impar:
        else:
            initAngle = heliosDist * 0.5 / rowDist
            if position[1] % 2 == 1:
                n = (position[1] - 1) * 0.5
                angle = initAngle + (heliosDist / rowDist) * n
            else:
                n = position[1] * 0.5
                angle = -initAngle - (heliosDist / rowDist) * (n - 1)
        # Coordenadas polares a cartesianas:
        x = rowDist * np.sin(angle)
        y = rowDist * np.cos(angle)
        return (x, y)

    def placeHeliostat(self, targetPoint, solarVector):
        '''
        Asigna una nueva posición a un heliostato y calcula sus coordenadas 
        cartesianas.
        
        targetPoint: Coordenadas XYZ del objetivo a enfocar (metros).
        solarVector: Vector de dirección solar.
        '''
        targetHeight = targetPoint[2]
        # Si se trata del primer heliostato del campo:
        if len(self.heliostats) == 0:
            rowDistance = self.firstRowDist
            maxHeliosInRow = self.maxHeliosInRow(rowDistance, targetHeight)
            position = (1, 1)
        # Si ya hay al menos un heliostato:
        else:
            xLast, yLast = self.heliostats[-1].getCoordinates()[:2]
            rowDistance = (xLast**2 + yLast**2)**0.5
            maxHeliosInRow = self.maxHeliosInRow(rowDistance, targetHeight)
            lastPosition = self.heliostats[-1].getPosition()
            # Si la fila no está llena aun:
            if lastPosition[1] < maxHeliosInRow:
                position = (lastPosition[0], lastPosition[1] + 1)
            # Si no caben más heliostatos en la fila:
            else:
                position = (lastPosition[0] + 1, 1)
                rowDistance += 0.5 * self.deltaRadDist(rowDistance, 
                                                       targetHeight)
                maxHeliosInRow = self.maxHeliosInRow(rowDistance, 
                                                     targetHeight)
        # Crear y añadir un nuevo heliostato a la lista:
        xyCoord = self.assignCartesianCoords(rowDistance, targetHeight, 
                                             position, maxHeliosInRow)
        newHelio = ht.Heliostat(xyCoord, targetPoint, solarVector, position)
        self.heliostats.append(newHelio)
        
    def draw(self):
        '''
        Representa gráficamente el campo solar.
        '''
        X, Y = [], []
        for helio in self.heliostats:
            x, y = helio.getCoordinates()[:2]
            X.append(x)
            Y.append(y)
        plt.xlabel('Coordenada X (m)')
        plt.ylabel('Coordenada Y (m)')
        plt.scatter(X, Y, s=5)
        plt.show()
        
    def isTowerShading(self, helioCoord, solarAngles):
        '''
        Determina si un heliostato está sombreado por la torre en el mediodía
        solar.
        
        helioCoord: Coordenadas cartesianas del heliostato (lista, metros)
        solarAngles: Ángulos solares.
        '''
        towerHeight = self.baseHeight + 30
        solarHeight = solarAngles[0]
        if solarHeight >= np.pi/2:
            return False
        xShadow = 5
        yShadow = towerHeight / np.tan(solarHeight)
        xCoord, yCoord = helioCoord[0], helioCoord[1]
        if abs(xCoord) <= xShadow and yCoord <= yShadow:
            return True
        return False
        
    def surfaceShaded(self, heliostat2, solarVector):
        '''
        Calcula la superficie que proyecta un heliostato sobre el que tiene
        detrás.
        
        heliostat2: Heliostato sombreado, 1 es para el que proyecta sombra.
        solarvector: Vector de posición solar.
        '''
        position2 = heliostat2.getPosition()
        coordinates2 = heliostat2.getCoordinates()
        heliostat1 = None
        # Buscar el heliostato que sombrea a heliostat2:
        for heliostat in self.heliostats:
            position1 = heliostat.getPosition()
            coordinates1 = heliostat.getCoordinates()
            # Si heliostat está en la fila anterior a heliostat2:
            if position1[0] == position2[0] - 1:
                xDistance = abs(coordinates2[0] - coordinates1[0])
                # Si existe solapamiento entre ellos:
                if xDistance < self.wideHelio:
                    horizontalShadow = self.wideHelio - xDistance
                    heliostat1 = heliostat
                    break
        # Si no recibe sombra de ninguno:
        if heliostat1 == None:
            return 0.0
        # Cálculo de la sombra vertical
        alpha1 = heliostat1.getOrientationAngles()[0]
        alpha2 = heliostat2.getOrientationAngles()[0]
        gamma = np.arcsin(solarVector[2])
        Y1 = coordinates1[1] + 0.5*self.heightHelio * np.cos(alpha1)
        Z1 = coordinates1[2] + 0.5*self.heightHelio * np.sin(alpha1)
        Y2 = coordinates2[1] - 0.5*self.heightHelio * np.cos(alpha2)
        Z2 = coordinates2[2] - 0.5*self.heightHelio * np.sin(alpha2)
        f1 = Y1*np.tan(gamma) + Y2*np.tan(alpha2) + (Z1 - Z2)
        f2 = np.tan(alpha2) + np.tan(gamma)
        Ys = f1 / f2
        # Si la sombra vertical es mas baja que la parte baja del heliostato:
        if Ys <= Y2:
            return 0.0
        Zs = np.tan(alpha2)*(Ys - Y2) + Z2
        verticalShadow = ((Ys - Y2)**2 + (Zs - Z2)**2)**0.5
        return horizontalShadow * verticalShadow
