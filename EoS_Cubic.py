import numpy as np
import matplotlib.pyplot as plt

class Cubic_eos:
    # Definição da classe de equações de estado cúbicas
    def __init__(self,nome):
        # Definição do nome
        self.Name = nome
        # Entrada dos parâmetros
        if nome=='vdW':
            sig = 0
            e = 0
            om = 1 / 8
            psi = 27 / 64
            zc = 3 / 8
            alp = lambda Tr, omega: 1
        elif nome == "RK":
            sig = 1
            e = 0
            om = .08664
            psi = .42748
            zc = 1 / 3
            alp = lambda Tr, omega: 1/np.sqrt(Tr)
        elif nome == "SRK":
            sig = 1
            e = 0
            om = .08664
            psi = .42748
            zc = 1 / 3
            alp = lambda Tr,omega: (1+(.480+1.574*omega-.176*omega**2)*(1-Tr**.5))**2
        elif nome == "PR":
            sig=1 + np.sqrt(2)
            e= 1 - np.sqrt(2)
            om=.07780
            psi=.45724
            zc=.30740
            alp = lambda Tr, omega: (1 + (.37464+1.54226*omega-.26992*omega**2) * (1 - Tr ** .5)) ** 2
        self.sigma = sig;
        self.ep = e;
        self.Omega = om;
        self.Psi = psi;
        self.Zc = zc;
        self.alpha = alp
        self.R = 83.14

    def CalcP(eos,comp,T,V):
        # Cálculo da pressão pela equação de estado cúbica

        # Cálculo dos parâmetros a e b
        a = eos.Psi * eos.alpha(T / comp.Tc, comp.omega)* eos.R**2. * comp.Tc**2. / comp.Pc
        b = eos.Omega * eos.R * comp.Tc / comp.Pc

        P = eos.R * T / (V - b) - a / (V + eos.ep * b) / (V + eos.sigma * b)
        return P

    def CalcPhi(eos,comp,T,V):
        # Calculo do coeficiente de atividade de compostos puros

        # Cálculo dos parâmetros a e b
        a = eos.Psi*eos.alpha(T/comp.Tc,comp.omega)*eos.R**2.*comp.Tc**2./comp.Pc
        b = eos.Omega*eos.R*comp.Tc/comp.Pc
        # Cálculando a pressão
        P = eos.R * T / (V - b) - a / (V + eos.ep * b) / (V + eos.sigma * b)
        # Definindo o fator de compressibilidade
        Z = P*V/eos.R/T
        # Cálculo dos parâmetros da equação
        beta = b * P / eos.R / T
        if eos.Name=='vdW':
            I = beta/Z
        else:
            I = 1/(eos.sigma-eos.ep) * np.log((Z+eos.sigma*beta)/(Z+eos.ep*beta))
        q_ =  eos.Psi*eos.alpha(T/comp.Tc,comp.omega)/eos.Omega/(T/comp.Tc)
        # Coeficiente de fugacidade
        lnphi = (Z-1) - np.log(Z-beta) - q_*I
        return lnphi,P

    def EoSRoots(eos,comp,T,P):
        # Cálculo dos volumes da equação de estado
        alp = eos.alpha(T/comp.Tc,comp.omega)
        A= eos.Psi*alp*eos.R**2*comp.Tc**2/comp.Pc
        B= eos.Omega*eos.R*comp.Tc/comp.Pc
        CA = A*P/eos.R**2/T**2
        CB = P*B/eos.R/T
        if eos.Name=="vdW":
            Vec = [1, -(1 + CB), CA, -CA*CB]
        elif eos.Name == "PR":
            Vec = [1, -(1 - CB),(CA -3*CB**2 - 2*CB), -(CA*CB - CB**2 - CB**3)]
        else:
            Vec = [1,-1, (CA-CB -CB**2), -CA*CB]
        Z = np.roots(Vec)
        V = Z * eos.R * T / P
        return V

    def ELVPure(eos,comp,T):
        # Cálculo do equilíbrio de fases

        # Definindo P0
        P0 = comp.Pc * np.exp(5.4 * (comp.omega + 1) * (1 - comp.Tc / T))
        # Encontrando o volume
        Vol = eos.EoSRoots(comp, T, P0)
        # Definindo o volume da fase líquida e da vapor
        Vliq = np.min(Vol)
        Vvap = np.max(Vol)
        # Cálculo do coeficiente de fugacidade
        lnphi_liq, P = eos.CalcPhi(comp, T, Vliq)
        lnphi_vap, P = eos.CalcPhi(comp, T, Vvap)
        lnphi = np.array([lnphi_liq, lnphi_vap])
        # Definindo o contador
        cont = 0
        # Calculando a nova pressão
        P1 = P0
        P0 = (P1 * np.exp(lnphi[0]) / np.exp(lnphi[1]))
        # Entrando no loop
        while abs(np.exp(lnphi[0]) / np.exp(lnphi[1]) - 1) > 1e-4:
            # Encontrando o volume
            Vol = eos.EoSRoots(comp, T, P0)
            # Definindo o volume da fase líquida e da vapor
            Vliq = min(Vol)
            Vvap = max(Vol)
            # Cálculo do coeficiente de fugacidade
            lnphi_liq, P = eos.CalcPhi(comp, T, Vliq)
            lnphi_vap, P = eos.CalcPhi(comp, T, Vvap)
            lnphi = np.array([lnphi_liq, lnphi_vap])
            cont += 1
            # Calculando a nova pressão
            P1 = P0
            P0 = np.real(P1 * np.exp(lnphi[0]) / np.exp(lnphi[1]))

        return P0,[Vliq,Vvap]

    def PlotGraf_PV(eos,comp,T1=0):
        # Definindo vetores
        T = np.linspace(comp.Ttriple,comp.Tc,101)
        P = np.zeros(len(T))
        V = np.zeros([len(T),2])
        # Calculando os dados de equilíbrio
        for i in range(len(T)):
            [P[i], V[i]] = eos.ELVPure(comp,T[i])
        # Definindo alguma temperatura para construção da isoterma
        if T1 == 0:
            T1 = comp.Tc
        # Definindo vetores
        Vr = np.logspace(np.log10(np.min(V)),np.log10(np.max(V)),1001)
        Pr = np.zeros(len(Vr))
        # Calculando o ponto crítico
        Pc,Vc = eos.ELVPure(comp, comp.Tc)
        # Calculando P em função de V na isoterma
        for i in range(len(Vr)):
            # Verificando se está no envelope ou fora
            if T1 <= comp.Tc:
                # Excluindo a região metaestável/instável
                Peq, Veq = eos.ELVPure(comp, T1)
                if Vr[i]<np.min(Veq) or Vr[i]>np.max(Veq):
                    Pr[i] = eos.CalcP(comp,T1,Vr[i])
                else:
                    Pr[i] = Peq
            else:
                Pr[i] = eos.CalcP(comp, T1, Vr[i])
        # Plotando
        p = plt.loglog(V,P,'k',Vc[0],Pc,'or',Vr,Pr,markersize=4)
        plt.xlabel(r'Volume molar (cm$^3$/mol)')
        plt.ylabel('Pressão (bar)')
        plt.legend([p[0],p[3],p[2]],['Curva de equilíbrio','Isoterma à '+str(np.round(T1))+'K','Ponto crítico'])
        plt.show()

    def Residual_Enthalpy(eos,V,T,comp):
        # Calculando o fator de compressibilidade
        Z = V/eos.R/T*eos.CalcP(comp,T,V)
        # Calculando a pressão
        P = Z * eos.R * T / V
        # Calculando os parâmetros da equação
        a = eos.Psi * (eos.R * comp.Tc) ** 2 / comp.Pc * eos.alpha(T / comp.Tc, comp.omega)
        b = eos.Omega * eos.R * comp.Tc / comp.Pc
        B = b * P / eos.R / T
        # Para cada equação tem-se diferentes fórmulas
        if eos.Name == 'PR':
            kappa = (.37464 + 1.54226 * comp.omega - .62992 * comp.omega ** 2)
            dadT = -a * kappa / np.sqrt(eos.alpha(T / comp.Tc, comp.omega) * T * comp.Tc)
            Hr = Z - 1 + (T * dadT - a) / (eos.R * T * (eos.sigma - eos.ep) * b) * np.log(
                (Z + eos.sigma * B) / (Z + eos.ep * B))
        elif eos.Name == 'SRK':
            kappa = (0.480 + 1.574 * comp.omega - .176 * comp.omega ** 2)
            dadT = -a * kappa / np.sqrt(eos.alpha(T / comp.Tc, comp.omega) * T * comp.Tc)
            Hr = Z - 1 + (T * dadT - a) / (eos.R * T * (eos.sigma - eos.ep) * b) * np.log(
                (Z + eos.sigma * B) / (Z + eos.ep * B))
        elif eos.Name == 'RK':
            q = eos.Psi / eos.Omega * eos.alpha(T / comp.Tc, comp.omega) / (T / comp.Tc)
            Hr = Z - 1 - 1.5 * q * np.log((Z + B) / (Z))
        elif eos.Name == 'vdW':
            Hr = Z - 1 - a / (eos.R * T * V)

        return Hr

    def Residual_Entropy(eos,V,T,comp):
        # Calculando o fator de compressibilidade
        Z = V / eos.R / T * eos.CalcP(comp, T, V)
        # Calculando a pressão
        P = Z * eos.R * T / V
        # Calculando os parâmetros
        a = eos.Psi * (eos.R * comp.Tc) ** 2 / comp.Pc * eos.alpha(T / comp.Tc, comp.omega)
        b = eos.Omega * eos.R * comp.Tc / comp.Pc
        B = b * P / eos.R / T
        # Cada equação tem uma fórmula
        if eos.Name == 'PR':
            kappa = (.37464 + 1.54226 * comp.omega - .62992 * comp.omega ** 2)
            dadT = -a * kappa / np.sqrt(eos.alpha(T / comp.Tc, comp.omega) * T * comp.Tc)
            Sr = np.log(Z-B) + dadT/(eos.R*(eos.sigma-eos.ep)*b)*np.log((Z+eos.sigma*B)/(Z+eos.ep*B))
        elif eos.Name == 'SRK':
            kappa = (0.480 + 1.574 * comp.omega - .176 * comp.omega ** 2)
            dadT = -a * kappa / np.sqrt(eos.alpha(T / comp.Tc, comp.omega) * T * comp.Tc)
            Sr = np.log(Z-B) + dadT/(eos.R*(eos.sigma-eos.ep)*b)*np.log((Z+eos.sigma*B)/(Z+eos.ep*B))
        elif eos.Name == 'RK':
            q = eos.Psi / eos.Omega * eos.alpha(T / comp.Tc, comp.omega) / (T / comp.Tc)
            Sr = np.log(Z-B) -0.5*q*np.log((Z+B)/(Z))
        elif eos.Name == 'vdW':
            Sr = np.log(Z)+np.log((V-b)/V)

        return Sr


# eos = Cubic_eos('PR')

# comp = Component('Water',1.801530e+01,6.471300e+02,2.205500e+02,5.594780e+01,2.290000e-01,3.448610e-01,273.16,0.006)

# lnphi,P = eos.CalcPhi(comp,300,30)
# print(lnphi,P)

# V = eos.EoSRoots(comp,300,1)
# print(V)
# P = eos.ELVPure(comp,373.15)
# print(P)
# eos.PlotGraf_PV(comp,680)

# Ti = 500
# Pi = 50
# Vi = 6.402029253201094e+02
# Hr = eos.Residual_Enthalpy(Vi,Ti,comp)*Ti*8.314
# print(Hr)
# Sr = eos.Residual_Entropy(Vi,Ti,comp)*8.314
# print(Sr)
