# import matplotlib.pyplot as plt
import numpy as np

class Component:
    def __init__(self,info):
        # Definindo os campos
        self.Name = ''
        self.CASNumber = ''
        self.MolarMass = 0
        self.Tc = 0
        self.Pc = 0
        self.Vc = 0
        self.Zc = 0
        self.omega = 0
        # Definindo dados experimentais de Psat
        self.dados = np.zeros(1)
        # Definindo parâmetros da sigma-MTC
        self.vcosmo = 0
        self.sigma = [0]
        self.acosmo = 0
        self.z = [0]
        self.V_i = 0
        self.q = 0
        self.r = 0
        self.alphaprime = 8287227.862
        self.cNE = 0
        self.cNE_T = 0
        self.cHB = 43069954.705586
        self.v = 0
        self.sigmaHB = 0.008400
        self.rnorm = 0

        with open('Dados/'+str(info)+'.txt',"r") as arq:
            cont = arq.readlines()
            for i in range(0,len(cont)):
                exec('self.'+cont[i])


    def Rackett(self,T):
        vol = self.Vc*self.Zc**((1-T/self.Tc)**(2/7))
        return vol

