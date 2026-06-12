# MKWS Project: Stability Limits of Methane Combustion with EGR

**Author:** Jan Piasecki (Student ID: 333486)  
**Course:** Computer Methods in Combustion (MKWS)  
**University:** Warsaw University of Technology (WUT)  

## Project Overview
This repository contains a numerical analysis of the impact of Exhaust Gas Recirculation (EGR) on the propagation and stability of a one-dimensional, laminar methane flame. The study utilizes the **Cantera** software package and the **GRI-Mech 3.0** kinetic mechanism to investigate flame behavior across varying air excess ratios ($\lambda \in [0.8, 1.3]$) and $CO_2$ dilution levels.

The primary objective was to determine the critical flame quenching limits by applying a physical minimum velocity criterion ($S_L < 3.5$ cm/s) coupled with an advanced numerical bisection method, enabling the construction of an operational stability map.

## Repository Structure
* `MKWS.py` - The main Python script utilizing Cantera 1D FreeFlame to compute laminar burning velocities, peak temperatures, and detect quenching limits via the bisection method.
* `main.tex` - The LaTeX source code for the final scientific report.
* `Tabela_Wynikow.txt` - Raw numerical output data from the general parameter sweep.
* `Tabela_Granic_Gaszenia.txt` - Precise EGR quenching limits computed via the bisection method.
* `Wykres_1_Predkosc.png` - Plot showing the suppression of laminar burning velocity.
* `Wykres_2_Mapa.png` - The computed Stability Limits Map (flammability vs. quenching regions).
* `Wykres_3_Temperatura.jpg` - Plot showing the reduction of peak adiabatic flame temperature.

## Dependencies and Setup
To run the Python script locally, you need a working Python environment with the following packages installed:
* `cantera`
* `numpy`
* `matplotlib`
