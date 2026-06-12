# MKWS Project: Stability Limits of Methane Combustion with EGR

**Author:** Jan Piasecki (Student ID: 333486)  
**Course:** Computer Methods in Combustion (MKWS)  
**University:** Warsaw University of Technology (WUT)  

## Project Overview
This repository contains a numerical analysis of the impact of Exhaust Gas Recirculation (EGR) on the propagation and stability of a one-dimensional, laminar methane flame. The study utilizes the **Cantera** software package and the **GRI-Mech 3.0** kinetic mechanism to investigate flame behavior across varying air excess ratios ($\lambda \in [0.8, 1.3]$) and $CO_2$ dilution levels.

The primary objective was to determine the critical flame quenching limits by applying a physical minimum velocity criterion ($S_L < 3.5$ cm/s) coupled with an advanced numerical bisection method, enabling the construction of an operational stability map.

### Source Code & Documentation
* `MKWS.py` – The main Python script. 
* `main.tex` – The complete LaTeX source code for the scientific report.
* `MKWS_Report.pdf` – The final compiled project report.
* `README.md` – This project documentation.
* `.gitignore` – Specifies untracked files to ignore.

### `Results/` Directory (Outputs)
* `Wykres_1_Predkosc.png` – Plot illustrating the suppression of laminar burning velocity.
* `Wykres_2_Mapa.png` – The physical Stability Limits Map (flammability vs. quenching regions).
* `Wykres_3_Temperatura.png` – Plot demonstrating the peak adiabatic flame temperature reduction.
* `Tabela_Wynikow.txt` – Raw text output containing the general parameter sweep results.
* `Tabela_Granic_Gaszenia.txt` – Precise EGR quenching limits determined by the bisection method.

## Dependencies and Setup
To run the Python script locally, you need a working Python environment with the following packages installed:
* `cantera`
* `numpy`
* `matplotlib`
