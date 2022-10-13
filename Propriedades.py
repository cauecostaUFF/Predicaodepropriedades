import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import EoS_Cubic as EoS
import COMPONENT as Comp

# from matplotlib.widgets import Button
# from matplotlib.widgets import TextBox

st.session_state.visibility = "collapsed"
st.session_state.disabled = True
st.set_page_config(
    page_title="Predição de deltaH de vaporização",
    # page_icon="🧊",
    layout="centered",
    initial_sidebar_state="expanded"
)

# Definindo modelos
modelos = {0: r"van der Waals", 1: r"Redlich-Kwong", 2: r"Soave-Redlich-Kwong", 3: r"Peng-Robinson"}
def format_func(option):
    return modelos[option]
nome = ['vdW','RK','SRK','PR']
# Definindo substâncias e suas propriedades
subs = {0:"Água",
        1:"Metano",
        2:"Octano",
        3:"Ciclohexano",
        4:"Benzeno",
        5:"Metanol",
        6:"Etanol"}
index = [1076,1,23,99,242,477,478]
def format_subs(option):
    return subs[option]
main_col1,main_col2,main_col3 = st.columns([10,10,10])

with main_col2:
    st.write(r"Selecione a substância:")
    option_subs = st.selectbox(
        "Água",
        options=list(subs.keys()), format_func=format_subs,
        label_visibility=st.session_state.visibility,
        disabled=False,
    )

with main_col1:
    # Selecionando modelo
    st.write(r"Selecione a equação de estado:")
    # st.write(r"Modelo de Margules: G$^{E}$ = n$\beta$RTx$_A$x$_B$")
    option_eq = st.selectbox(
        "van der Waals",
        options=list(modelos.keys()),format_func=format_func,
        label_visibility=st.session_state.visibility,
        disabled=False,
    )
with main_col3:
    Graph = st.radio("Escolha o gráfico",
        ('Pressão de vapor', 'Entalpia de vaporização', 'Segundo coeficiente do virial'))

eos = EoS.Cubic_eos(nome[option_eq])
comp = Comp.Component(index[option_subs])
if Graph == 'Pressão de vapor':
    Texp = comp.dados[0]
    Pexp = comp.dados[1]
    Tspan = np.linspace(min(Texp), max(Texp), 200)
    Pvap = np.zeros(len(Tspan))
    for i in range(len(Tspan)):
        T = Tspan[i]
        [Pvap[i], Vc] = eos.ELVPure(comp, T)
    p = plt.semilogy(Texp, Pexp, 'ok',Tspan, Pvap,  markersize=4)
    plt.xlabel(r'Temperatura (K)')
    plt.ylabel(r'Pressão (bar)')
elif Graph == "Entalpia de vaporização":
    Texp = comp.H_vap[0]
    Hexp = comp.H_vap[1]
    R = 8.314 #J/molK
    Tspan = np.linspace(min(Texp),max(Texp),200)
    dHvap = np.zeros(len(Tspan))

    for i in range(len(Tspan)):
        T = Tspan[i]
        [Pc, Vc] = eos.ELVPure(comp, T)
        Hr_liq = eos.Residual_Enthalpy(Vc[0], T, comp)
        Hr_vap = eos.Residual_Enthalpy(Vc[1], T, comp)
        dHvap[i] = (Hr_vap - Hr_liq) * T * R / 1000

    p = plt.plot(Texp,Hexp,'ok',Tspan,dHvap,markersize=4)
    plt.xlabel(r'Temperatura (K)')
    plt.ylabel(r'$\Delta H^{vap}$ (kJ/mol)')
elif Graph == "Segundo coeficiente do virial":
    Texp = comp.Second_Virial_Coefficient[0]
    Bexp = comp.Second_Virial_Coefficient[1]
    Tspan = np.linspace(min(Texp), max(Texp), 200)
    B = np.zeros(len(Tspan))

    for i in range(len(Tspan)):
        T = Tspan[i]
        a = eos.Psi * eos.alpha(T / comp.Tc, comp.omega) * eos.R ** 2. * comp.Tc ** 2. / comp.Pc
        b = eos.Omega * eos.R * comp.Tc / comp.Pc
        B[i] = b - a / eos.R / T

    p = plt.plot( Texp, Bexp, 'ok',Tspan, B, markersize=4)
    plt.xlabel(r'Temperatura (K)')
    plt.ylabel(r'$B$ (cm$^3$/mol)')
plt.legend([p[0],p[1]],["Dados experimentais","Predição com a equação "+eos.Name])
plt.show()
st.pyplot(plt)
col_end1,col_end2,col_end3 = st.columns(3)
with col_end3:
    st.write("Desenvolvido por: Cauê Costa\nE-mail: cauecosta@id.uff.br")
with col_end2:
    st.image("UFF.png", width=150)
with col_end1:
    st.image("molmod.png", width=150)