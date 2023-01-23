#pip install pygmo

import matplotlib.pyplot as plt
from scipy.spatial import Delaunay
import numpy as np
import random as rd
import math
from scipy import integrate
import plotly.figure_factory as ff
import random as rand
import pygmo as pg

'''
Se define la función para la cual queremos crear la malla, incluyendo los 
límites en los que se definen x e y, y sus derivadas parciales.
'''

f = lambda y,x: (x*y*(x**2-y**2))/(x**2+y**2+1)
a = 0
b = 3
at = 0
bt = 3
fx = lambda y,x: (y*(x**4+x**2*(4*y**2+3)-y**2*(y**2+1)))/(x**2+y**2+1)**2
fy = lambda y,x: x*(-y**4-4*y**2*x**2-3*y**2+x**4+x**2)/(x**2+y**2+1)**2

def getArea(triangulos):
  '''
  Calcula el área de los triangulos del grid a partir de los pares ordenados 
  de sus vértices.
  '''

  areas = []
  for i in triangulos:
      x1 = i[0][0]
      x2 = i[1][0]
      x3 = i[2][0]
      y1 = i[0][1]
      y2 = i[1][1]
      y3 = i[2][1]
      area = (1/2)*abs(x1*(y2 - y3) + x2*(y3 - y1) + x3*(y1 - y2))
      areas.append(area)
  return np.array(areas)

def getAreaF(triangulos3d):
  '''
  Calcula el área de los triangulos de la malla a partir de los puntos 
  coordenados de sus vértices.
  '''

  areas = []
  for i in triangulos3d:
      x1 = i[0][0]
      x2 = i[1][0]
      x3 = i[2][0]
      y1 = i[0][1]
      y2 = i[1][1]
      y3 = i[2][1]
      z1 = i[0][2]
      z2 = i[1][2]
      z3 = i[2][2]
      u = np.array([x2-x1,y2-y1,z2-z1])
      v = np.array([x3-x1,y3-y1,z3-z1])
      area = np.linalg.norm(np.cross(u,v)/2)
      areas.append(area)
  return areas

def desv_est(areas):
  '''
  Calcula la desviación estándar del área de los triangulos
  '''

  m = len(areas)
  promArea = np.mean(areas)
  
  desv = 0
  for i in areas:
    desv = desv + abs(i-promArea)**2/m
  return desv

def error(triangulos,fy,fx):
  '''
  Calcula el error absoluto del área aproximada de la malla. Se calcula el área 
  real mediante una integral de superficie.
  '''

  areas = getAreaF(triangulos)
  g = lambda y,x: math.sqrt(1+fx(y,x)**2+fy(y,x)**2)
  areaReal = integrate.dblquad(g,0,3, lambda x: 0, lambda x: 3)[0]
  areaAprox = sum(areas)
  error = abs(areaReal-areaAprox)
  return error

def errorgrid(triangulos, areacuad):
  '''
  Calcula el error absoluto del área aproximada del grid. 
  '''

  areas = getArea(triangulos)
  areaAprox = sum(areas)
  errorgrid = abs(areacuad-areaAprox)
  return errorgrid

def getGenotipos(pop,grid):
  '''
  Definimos el genotipo de la población, es una traducción de los cromosomas.
  Cada genotipo es un conjunto de pares ordenados que representan los nodos en dicho grid.
  '''

  genotipos = []
  for p in pop:
    genotipo = []
    for j in range(len(p)):
      if p[j] == 1:
        genotipo.append(grid[j])
      else:
        continue
    genotipos.append(np.array(genotipo))
  return genotipos

def getMutation(crom,grid,esq1,esq2,esq3,esq4):
    j = rand.randint(0,len(crom)-1)
    i = rand.random()
    if i>= 0.5:
      if (crom[j] == 0 and not((grid[j] == esq1).all() or (grid[j] == esq2).all()) or ((grid[j] == esq3).all() or (grid[j] == esq4).all())):
        crom[j] = 1
    else:
      if (crom[j] == 1 and not((grid[j] == esq1).all() or (grid[j] == esq2).all()) or ((grid[j] == esq3).all() or (grid[j] == esq4).all())):
        crom[j] = 0
    return crom

def getCrossover(crom1, crom2,pc):
  length = len(crom1)
  r = rand.random()
  if length < 2 or r > pc:
      return crom1, crom2

  p = rand.randint(1, length - 1)
  return np.concatenate((crom1[0:p],crom2[p:])), np.concatenate((crom2[0:p],crom1[p:]))

def getDistancias(desviaciones, errores):
  distancia1 = []
  distancia2 = []
  desvEnum = []
  errEnum = []
  for i in range(len(desviaciones)):
    errEnum.append([errores[i],i])
    desvEnum.append([desviaciones[i],i])

  desvEnum = sorted(desvEnum, key=lambda x: x[0])
  errEnum = sorted(errEnum, key=lambda x: x[0])

  distancia1.append([float('Inf'),desvEnum[0][1]])
  distancia2.append([float('Inf'),errEnum[0][1]])
  
  for i in range(1,len(errores)-1):
    distancia1.append([(desvEnum[i+1][0]-desvEnum[i-1][0])/(float(max(desvEnum[:][0]))-float(min(desvEnum[:][0]))),desvEnum[i][1]])
    distancia2.append([(errEnum[i+1][0]-errEnum[i-1][0])/(float(max(errEnum[:][0]))-float(min(errEnum[:][0]))), errEnum[i][1]])
  distancia1.append([float('Inf'),desvEnum[-1][1]])
  distancia2.append([float('Inf'),errEnum[-1][1]])

  distancia1 = sorted(distancia1, key=lambda x: x[1])
  distancia2 = sorted(distancia2, key=lambda x: x[1])

  distancias = []
  for i in range(len(distancia1)):
    distancias.append((distancia1[i][0]+distancia2[i][0])/2)
  return distancias

def getSelection(ndr, distancias, pop, npop):
  fittest = []
  while len(fittest)<npop:
    r1 = rand.randint(0,npop-1)
    r2 = rand.randint(0,npop-1)
    if ndr[r1]<ndr[r2]:
      fittest.append(pop[r1])
    if ndr[r2]<ndr[r1]:
      fittest.append(pop[r2])
    elif ndr[r2]==ndr[r1]:
      if distancias[r1]>distancias[r2]:
        fittest.append(pop[r1])
      elif distancias[r2]>distancias[r1]:
        fittest.append(pop[r2])
      else: 
        continue
  return fittest

def NSGA2(grid, npop, niter, f,pc):  
    
  #Definimos los cromosomas de la población inicial, siempre incluimos las esquinas.
  pop0 = []
  esq1 = [0.,0.]
  esq2 = [0., 3.]
  esq3 = [3., 0.]
  esq4 = [3., 3.]
  
  for i in range(npop):
    crom = np.zeros(len(grid))
    for j in range(len(crom)):
      crom[j] = rand.randint(0,1)
      if ((grid[j] == esq1).all() or (grid[j] == esq2).all()) or ((grid[j] == esq3).all() or (grid[j] == esq4).all()):
        crom[j] = 1
      else:
        continue
    pop0.append(crom)

  #Para cada grid en genotipos hacemos una triangulación de Delaunay y calculamos
  #su error y las desviaciones estándar de las áreas de sus triangulos.
  genotipos = getGenotipos(pop0,grid)
  errores = []
  desviaciones = []
  for g in genotipos:
    tri = Delaunay(g)
    triangulos = g[tri.simplices]

    eval = []
    for x in g:
      eval.append(f(x[1],x[0]))
    eval = np.array(eval)
    
    triangulos3d = np.zeros([len(triangulos),len(triangulos[0]),3]) 
    for i in range(len(triangulos)):
      for j in range(len(triangulos[i])):
        triangulos3d[i][j] = np.append(triangulos[i][j],f(triangulos[i][j][1],triangulos[i][j][0]))

    areas3d = getAreaF(triangulos3d)
    error3d = error(triangulos3d, fy, fx)
    de = desv_est(areas3d)
    errores.append(error3d)
    desviaciones.append(de)
  plt.scatter(desviaciones,errores)
  plt.show()  

  papas = pop0.copy()
  
  #Se itera para niter generaciones
  for i in range(niter):
    
    #Creamos hijos de la poblacion y lo unimos a la población inicial
    hijos = []
    while len(hijos)<npop:
      r1 = rand.randint(0,npop-1)
      r2 = rand.randint(0,npop-1)
      if r1 == r2:
        continue
      else:
        a, b = getCrossover(papas[r1],papas[r2],pc)
        hijos.append(a)
        hijos.append(b)
    
    Rt = np.array(papas+hijos)

  
    #Mutación
    for i in range(2*npop):
      p = 1/len(crom)
      r = rand.random()
      if r<=p:
        Rt[i] = getMutation(Rt[i],grid,esq1,esq2,esq3,esq4)
    
    
    #Para cada grid en el genotipo hacemos una triangulación de Delaunay y calculamos
    #su error y las desviaciones estándar de las áreas de sus triangulos.
    genotipos = getGenotipos(Rt,grid)
    errores = []
    desviaciones = []
    for g in genotipos:
      tri = Delaunay(g)
      eval = []

      triangulos = g[tri.simplices]
      for x in g:
        eval.append(f(x[1],x[0]))
      eval = np.array(eval)
      
      triangulos3d = np.zeros([len(triangulos),len(triangulos[0]),3]) 
      for i in range(len(triangulos)):
        for j in range(len(triangulos[i])):
          triangulos3d[i][j] = np.append(triangulos[i][j],f(triangulos[i][j][1],triangulos[i][j][0]))

      areas3d = getAreaF(triangulos3d)
      error3d = error(triangulos3d, fy, fx)
      de = desv_est(areas3d)
      errores.append(error3d)
      desviaciones.append(de)
    
    
    #Ordenamiento no dominado, ndf es la lista de los frentes no dominados y ndr 
    #es el ranking de no dominancia.
    criterios = []
    for i in range(len(errores)):
      criterios.append([errores[i],desviaciones[i]])
    
    ndf, dl, dc, ndr = pg.fast_non_dominated_sorting(points = criterios)
    distancias = getDistancias(desviaciones, errores)
    papas = getSelection(ndr,distancias,Rt,npop)
    
    
  
  return papas

def main(a,b,at,bt,M,N, npop, niter, pc):
  '''
  Función principal.
  
  Parametros
  ----------
  
  [a,b]x[at,bt] es el intervalo en el que está definida la función.
  M y N son la cantidad de nodos en los ejes x e y respectivamente.
  npop es el número de individuos en la población, niter es el número 
  de iteraciones del algoritmo NSGA2 y pc es la probabilidad de cruce.
  '''  

  #Se define un grid uniforme de MxN nodos.
  gridl = []
  for x in np.linspace(a,b,M):
    for y in np.linspace(at,bt,N):
      gridl.append([x,y])
  grid = np.array(gridl)
  tri = Delaunay(grid)
  plt.triplot(grid[:,0], grid[:,1], tri.simplices)
  plt.plot(grid[:,0], grid[:,1], 'o')
  plt.show()
  triangulos1 = grid[tri.simplices]
  areas2d = getArea(triangulos1)
  error2d = errorgrid(triangulos1,9)

  triangulos3d1 = np.zeros([len(triangulos1),len(triangulos1[0]),3]) 
  for i in range(len(triangulos1)):
    for j in range(len(triangulos1[i])):
      triangulos3d1[i][j] = np.append(triangulos1[i][j],f(triangulos1[i][j][1],triangulos1[i][j][0]))
  areas3d = getAreaF(triangulos3d1)
  error3d = error(triangulos3d1, fy, fx)
  des = desv_est(areas3d)


  #Se llama al algoritmo genético
  cromosomas = NSGA2(grid,npop,niter, f,pc)

  genotipos = getGenotipos(cromosomas,grid)

  
  #Quitamos los elementos repetidos
  m = {}
  for elem in genotipos:
    m[str(elem)] = elem
  gridoptimos = []
  for value in m.values():
    gridoptimos.append(value)

  errores = []
  desviaciones = []
  for g in gridoptimos:
    tri = Delaunay(g)
    triangulos = g[tri.simplices]

    eval = []
    for x in g:
      eval.append(f(x[1],x[0]))
    eval = np.array(eval)
    
    triangulos3d = np.zeros([len(triangulos),len(triangulos[0]),3]) 
    for i in range(len(triangulos)):
      for j in range(len(triangulos[i])):
        triangulos3d[i][j] = np.append(triangulos[i][j],f(triangulos[i][j][1],triangulos[i][j][0]))

    areas3d = getAreaF(triangulos3d)
    error3d = error(triangulos3d, fy, fx)
    de = desv_est(areas3d)
    errores.append(error3d)
    desviaciones.append(de)

  plt.scatter(desviaciones,errores)
  plt.show()  

  for i in gridoptimos:
    print(len(i))  
    eval = []

    tri2 = Delaunay(i)
    plt.triplot(i[:,0], i[:,1], tri2.simplices)
    plt.plot(i[:,0], i[:,1], 'o', 'b')
    plt.show()
    triangulos2 = grid[tri2.simplices]
    
    #Se evalúa cada nodo en la función
    for x in i:
        eval.append(f(x[1],x[0]))
    eval = np.array(eval)
    

    triangulos3d = np.zeros([len(triangulos2),len(triangulos2[0]),3]) 
    for k in range(len(triangulos2)):
      for j in range(len(triangulos2[k])):
          triangulos3d[k][j] = np.append(triangulos2[k][j],f(triangulos2[k][j][1],triangulos2[k][j][0]))

    
    #Se grafica la malla
    fig = ff.create_trisurf(i[:,0],i[:,1], eval, colormap='Viridis', simplices = tri2.simplices)
    fig.show()
    
  

main(0,3,0,3,20,20,100,10,1)
