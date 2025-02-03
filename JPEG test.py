#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 22 21:49:59 2025

@author: Alexandre

Reste à faire
-Étape de compression décompréssion avec RLE (fait) et Huffman
-Mesure du score de qualité et de compression (fait)
-Ajout de sécurité bord et dimension
-Division d'une image en petit pixel et recomposition 
-Optimisation de l'algo de calcul de la DCT
-Ajouter des matrices de quantizations
-Couleurs
"""
#Import des packages nécessaire

from PIL import Image
import os
import numpy as np
import random
import matplotlib.pyplot as plt


D=8 #Dimension de l'espace 8 par défaut comme en JPEG

#Choix du Current working directory
path = ''
os.chdir(path)

            

class ImageD_D: #Classe sur un block  D par D
    def __init__(self, nom, matrice):  
        self.nom = nom #purement esthétique
        self.contenue = np.array(matrice) #Matrice de la DCT
        self.valide=True
    def check_matrice(self): #Vérifie la validité des donnés
        if self.contenue.shape != (D,D) :
            self.Valide = False
            print("error")
        for i in range(0,D):
             for j in range(0,D):
                if self.contenue[i,j]> 255 or self.contenue[i,j]<0:
                    self.valide = False
                    print("error")
                    break
    def show(self, trans): #Représente l'image trans permet de choisir si on veut l'image comprimer ou pas
        self.check_matrice()
        if self.valide:
            if trans==True: 
                M=self.comp
            else:
                M=self.contenue
            plt.figure(figsize=(D, D))
            plt.imshow(M, cmap="gray", vmin=0, vmax=255)
            plt.colorbar(label="Intensity (0-255)")
            plt.title("Visualisation image " +self.nom)
            plt.axis("off")  # Turn off axis labels for a cleaner look
            plt.show()
            return self
        else:
            print("Matrice sous un mauvais format")
    def DCT(self): # Calcule la DCT pixel par pixel
        M_DCT=[]
        for i in range(0,D):
            M_DCT.append([])
            for j in range(0,D):
                M_DCT[i].append(DCT(self.contenue,i,j))
        self.DCT=np.array(M_DCT)
        return self
    def quantization(self,Q): #Division par la matrice Q pour supprimer l'information la moins importante
        self.DCT=np.round(self.DCT/Q).astype(int)
        return self
    def RLE(self): #Première étape de compression
        self.line = [] 
        self.RLE=[]
        reverse = True  
        for diag_start in range(D): 
#cette étape lit la matrice diagonale par diagonal de sorte que les plus gros coeeficient soit au début généré à l'aide de l'IA
            i, j = 0, diag_start
            diagonal = []
            while j >= 0 and i < D:
                diagonal.append(self.DCT[i,j])
                i += 1
                j -= 1
            if reverse:
                diagonal.reverse()
            self.line.extend(diagonal)
            reverse = not reverse
        for diag_start in range(1, D):
            i, j = diag_start, D - 1
            diagonal = []
            while i < D and j >= 0:
                diagonal.append(self.DCT[i,j])
                i += 1
                j -= 1
            if reverse:
                diagonal.reverse()
            self.line.extend(diagonal)
            reverse = not reverse
#Cette seconde étape remplace les suite de N 0 par le caracthère spéciale [0,N]
        k=-1
        cpt=0
        for i in self.line: 
            if i==0:
                if k!=0:
                    self.RLE.append([0,1])
                    cpt+=1
                else:
                    self.RLE[cpt-1]=[0,self.RLE[cpt-1][1]+1]
            else:
                self.RLE.append(i)
                cpt+=1
            k=i
        return self
    def decompression():  #Algo inverse 
        pass
    def dequantization(self,Q): #inverse la quantisation
        self.DCT=self.DCT*Q
        return self
    def inv_DCT(self): #Calcule pixel par pixel la DCT inverse
        comp=[]
        for x in range(0,D):
            comp.append([])
            for y in range(0,D):
                comp[x].append(inv_DCT(self.DCT,x,y))
        self.comp=np.array(comp)
        self.comp=np.round(self.comp).astype(int)
        return self
    def score(self):
        self.erreur=(np.sum(abs(self.comp-self.contenue))/64)
        self.compression=(len(self.RLE)/64)*100
        print("L'erreur moyenne des pixels pour l'image",self.nom,"est :",self.erreur," px\nSon taux de compression est :",self.compression, "%")
        return self
    def compress(self): #Fait les étapes une par une jusqu'à la compression
        self.show(False).DCT().quantization(Q).RLE()
        return self
    def decompress(self): #Fait les étapes de décompression
        self.dequantization(Q).inv_DCT().score().show(True)

    
def a(n): #Calcul des coefficient de normalisation
    if n==0:
        return np.sqrt(1/D)
    else:
        return np.sqrt(2/D)

def DCT(matrice,i,j): #Fonction de calcule de la DCT en 1 point note: algo sub optimale
    DCT=0
    for x in range(0,D):
        for y in range(0,D):
            coef1=((2*x+1)*i*np.pi)/(2*D)
            coef2=((2*y+1)*j*np.pi)/(2*D)
            DCT=DCT+matrice[x,y]*np.cos(coef1)*np.cos(coef2)
    DCT= a(i)*a(j)*DCT
    return DCT
    
def inv_DCT(DCT,x,y): #Calcule la transformation inverse de la DCT
    i_DCT=0
    for u in range(0,D):
        for v in range(0,D):
            coef1=((2*x+1)*u*np.pi)/(2*D)
            coef2=((2*y+1)*v*np.pi)/(2*D)
            i_DCT=i_DCT+a(u)*a(v)*DCT[u,v]*np.cos(coef1)*np.cos(coef2)
    return i_DCT


#Quantization Matrix from https://hal.science/hal-03858141/file/qtable.pdf
Q=[1,2,2,3,5,8,10,12,2,2,3,4,5,12,12,11,3,3,3,5,8,11,14,11,3,3,4,6,10,17,16,12,4,4,7,11,14,22,21,15,5,7,11,13,16,21,23,18,10,13,16,17,21,24,24,20,14,18,19,20,22,20,21,20]       
Q=np.array(Q)
Q=Q.reshape((8,8))

#Génération d'images test entiérement noire, en damier noire blanc, avec un gradient fort avec un gradient plus faible, au hasard
B=[[0]*D]*D 

C=[]
for i in range(1,9):
     C.append([])
     for j in range(1,9):
        C[i-1].append(((i+j)%2)*255)
G=[]
for i in range(1,9):
     G.append([])
     for j in range(1,9):
        G[i-1].append(4*i*j-1)
S=[]
for i in range(1,9):
     S.append([])
     for j in range(1,9):
        S[i-1].append(255-(i+j)*5)     
R=[]
for i in range(1,9):
     R.append([])
     for j in range(1,9):
        R[i-1].append(random.randint(0,255))


#Test des images test
ImageD_D("Standard",S).compress().decompress()
ImageD_D("Black",B).compress().decompress()
ImageD_D("Chess",C).compress().decompress()
ImageD_D("Gradient",G).compress().decompress()
ImageD_D("Random",R).compress().decompress()
