#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Config file for the test detector using ntd technology.
Basically set up the simulation of the detector.

@author: misiak
"""

import sympy as sy

import ethem as eth

#==============================================================================
# SYSTEM
#==============================================================================
syst = eth.System()

### Defining time and frequency variables
time, freq = syst.time, syst.freq

### Defining the thermal system
### cryostat
cryo = eth.Thermostat(syst, 'b')
### absorber thermal bath
abso = eth.ThermalBath(syst, 'a')
### ntd phonon bath
phntd = eth.ThermalBath(syst, 'p')
### ntd thermal bath
thntd = eth.ThermalBath(syst, 'ntd')
### thermal leak between ntd and cryo
leak = eth.ThermalLink(cryo, phntd, 'leak')
### glue between absorber and ntd
glue = eth.ThermalLink(abso, phntd, 'glue')
### ep coupling
epcoup = eth.ThermalLink(phntd, thntd, 'ep')

### Chassis ground
ground = eth.Voltstat(syst, 'ground')
ground.voltage = 0
### Wire capacitance
capa = eth.Capacitor(syst, 'f')
### NTD resistance
elntd = eth.Resistor(capa, ground, 'ntd')

#==============================================================================
# PHYSICAL RELATIONS AND ADDITIONNAL SYMBOLS
#==============================================================================
# temperature of NTD resistor is temperature of the NTD electron bath
elntd.temperature = thntd.temperature

#bias current for AC electronic, as main_flux on capa
capa.ibias = sy.symbols('I_bias')
capa.current = capa.ibias

# NTD characteristics
R0, T0 = sy.symbols('R0, T0')
elntd.R0, elntd.T0 = R0, T0
elntd.resistivity = eth.ntd_char(R0, T0, thntd.temperature)

# Joule Power from NTD resistor to NTD electron bath
thntd.power = eth.joule_power(capa.voltage, elntd.resistivity)

# Volume of absorber from its mass
abso.density, abso.mass = sy.symbols('D_Ge, M_a')
abso.volume = abso.mass / abso.density

# Volume of NTD from its dimensions, same volume for phntd, thntd and elntd
elntd.height, elntd.length, elntd.width = sy.symbols('H_ntd, L_ntd, W_ntd')
elntd.volume = elntd.height * elntd.length * elntd.width
phntd.volume = thntd.volume = elntd.volume

# Thermal Capacity expression in germanium
thntd.th_capa_coeff, phntd.th_capa_coeff = sy.symbols('ce_Ge, cp_Ge')
abso.th_capa_coeff = thntd.th_capa_coeff
abso.th_capacity = abso.volume * abso.th_capa_coeff * abso.temperature**3
phntd.th_capacity = phntd.volume * thntd.th_capa_coeff * phntd.temperature**3
thntd.th_capacity = thntd.volume * phntd.th_capa_coeff * thntd.temperature


# Power expression in gold link
leak.surface, leak.cond_alpha, leak.cond_expo = sy.symbols('S_Au, g_Au, n_Au')
leak.power = eth.kapitsa_power(leak.surface*leak.cond_alpha,
                               leak.cond_expo,
                               leak.from_bath.temperature,
                               leak.to_bath.temperature)

# Power expression in glue link
glue.cond_alpha, glue.cond_expo = sy.symbols('g_glue, n_glue')
glue.power = eth.kapitsa_power(glue.cond_alpha,
                               glue.cond_expo,
                               glue.from_bath.temperature,
                               glue.to_bath.temperature)

# Power expression in epcoup link
epcoup.cond_alpha, epcoup.cond_expo = sy.symbols('g_ep, n_ep')
epcoup.power = eth.kapitsa_power(phntd.volume*epcoup.cond_alpha,
                                 epcoup.cond_expo,
                                 epcoup.from_bath.temperature,
                                 epcoup.to_bath.temperature)

#==============================================================================
# NOISE POWER
#==============================================================================
# TFN noise for each link
for link in [glue, leak, epcoup]:
    tfn = eth.tfn_noise(link.conductance,
                        link.from_bath.temperature,
                        link.to_bath.temperature)
    tfn = tfn**0.5 # to obtain the LPSD
    link.noise_flux['TFN '+link.label] = tfn

# Johnson noise for each
for resi in [elntd,]:
    john_voltage = eth.johnson_noise(resi.resistivity, resi.temperature)
    john_voltage = john_voltage**0.5 # to obtain the LPSD
    john_current = john_voltage/resi.resistivity # to obtain the noise current
    resi.noise_flux['Johnson '+resi.label] = john_current

# amplifier current noise (impact the system, and so the observer)
i_a1, i_a2, i_a3 = sy.symbols('i_a1, i_a2, i_a3')
capa.i_a1 = i_a1
capa.i_a2 = i_a2
capa.i_a3 = i_a3
noise_current = (i_a1**2 + i_a2**2 *freq + i_a3**2 *freq**2)**0.5
capa.noise_sys['Ampli. Current'] = noise_current

# dac current noise (impact system and so observer)
e_dac = sy.symbols('e_dac')
capa.e_dac = e_dac
dac_current = capa.e_dac * sy.I * 2 * sy.pi * freq * 10e-12
capa.noise_sys['DAC Current'] = dac_current

# amplifier voltage noise (impact the observer only)
e_a1, e_a2, e_a3 = sy.symbols('e_a1, e_a2, e_a3')
capa.e_a1 = e_a1
capa.e_a2 = e_a2
capa.e_a3 = e_a3
noise_voltage = (e_a1**2 + e_a2**2 /freq + e_a3**2 /freq**2)**0.5
capa.noise_obs['Ampli. voltage'] = noise_voltage

#==============================================================================
# UPDATING THE SYSTEM
#==============================================================================
syst.build_sym()

#==============================================================================
# EVENT PERTURBATION
#==============================================================================
energy, tau_therm, eps = sy.symbols('E, tau_th, eps', positive=True)

per = eth.Perturbation(syst, 
                       energy,
                       [1-eps, 0., eps, 0.],
                       [tau_therm, tau_therm, tau_therm, tau_therm])

#==============================================================================
# EVALUATION DICT
#==============================================================================
evad_const = {syst.kB : 1.3806485e-23, #J.K-1
              abso.density : 5.32, # g.cm-3
              thntd.th_capa_coeff : 1.03e-6, # J.K-2.cm-3
              phntd.th_capa_coeff : 2.66e-6, # J.K-4.cm-3
              }


evad_sys = {
            elntd.temperature : 20e-3, #K
            elntd.resistivity : 3e6, #Ohms
            capa.capacity : 70e-12, # F

}

evad_sys = {leak.surface : 0.40, # cm2
            glue.cond_alpha : 1.66e-7, # W.K-1
            glue.cond_expo : 1., # 1
            epcoup.cond_alpha : 216., # W.K-6.cm-3
            epcoup.cond_expo : 6., # 1
            leak.cond_alpha : 5e-3, # W.K-4.cm-2
            leak.cond_expo : 4., # 1
            capa.capacity : 70e-12, # F
            abso.mass : 32, # g
            elntd.length :0.2, # cm
            elntd.width :0.2, # cm
            elntd.height :0.045, # cm
            R0 : 5.79, # Ohms
            T0 : 3.17, # K
            cryo.temperature : 18e-3, # K
            capa.ibias: 1e-9, #A
#            cryo.temperature : 19e-3, # K
#            bias.voltage : 2., #V
}

evad_per = {
            tau_therm : 5e-3, # s
            energy : 1e3 * 1.6e-19, # J
            eps : 0.1, #fraction
}

evad_noise = {e_a1: 1.62e-9,
              e_a2: 9.62e-9,
              e_a3: 3.09e-8,
              i_a1: 8e-17,
              i_a2: 4.87e-17,
              i_a3: 7.5e-19,
              e_dac: 2.05e-8}

evad = dict()
evad.update(evad_const)
evad.update(evad_sys)
evad.update(evad_per)
evad.update(evad_noise)

def get_eval_dict():
    return evad

def get_system():
    return syst

### checking the completeness of the evaluation dictionnary
# free symbols without evaluation
free_set = set(syst.phi_vect)|{time,freq}

# checking the electro-thermal equations
ete_free = syst.eteq.subs(evad).free_symbols
assert ete_free.issubset(free_set)

## checking the event perturbation
#per_free = per.matrix.subs(evad).free_symbols
#assert per_free.issubset(free_set)

# checking the noise power
for e in syst.elements_list:

    if isinstance(e, eth.RealBath):

        for v in e.noise_obs.values():
            noise_free = v.subs(evad).free_symbols
            assert noise_free.issubset(free_set)

        for v in e.noise_sys.values():
            noise_free = v.subs(evad).free_symbols
            assert noise_free.issubset(free_set)

    if isinstance(e, eth.Link):

        for v in e.noise_flux.values():
            noise_free = v.subs(evad).free_symbols
            assert noise_free.issubset(free_set)

if __name__ == '__main__':
    eth.sys_scheme(syst)
