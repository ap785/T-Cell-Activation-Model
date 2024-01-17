# T Cell Activation Mathematical Model
ODE-based Mathematical Model of T cell activation following antigen stimulation over time

This repository contains a mathematical model for T cell activation, implemented in Python using NumPy, SciPy, and Matplotlib. The model explores different scenarios, including variations with and without costimulatory signals.

## Table of Contents

- [Overview](#overview)
- [Usage](#usage)
- [Models](#models)
- [Parameters](#parameters)
- [Results](#results)
- [License](#license)

## Overview

The mathematical model simulates T cell activation under different conditions, considering the impact of costimulatory signals such as CD2, ICAM, CD28, CD27, and 4-1BB. The implemented models are based on ordinary differential equations (ODEs) to describe the dynamics of key variables over time.

![Figure Description]()


## Usage

To use the model, follow these steps:

1. Clone the repository to your local machine:

   ```bash
   git clone https://github.com/ap785/TCellActivationModel.git

## Parameters

$k_{on} = $ TCR-pMHC association rate ($\mu m^{2}s^{-1}$)

$k_{off} = $ TCR-pMHC dissociation rate ($s^{-1}$)

$k_{syn} = $ basal TCR upregulation ($s^{-1}\mu m^{-2}$)

$k_{basal} = $ basal (ligand-independent) TCR down-regulation ($s^{-1}$)

$k_{down} = $ ligand-induced TCR down-regulation ($s^{-1}$)

$k_{act} = $ production of switch / signalling amplification ($s^{-1}$)

$k_{prod} = $ cytokine production rate ($s^{-1}$)

$k_{deg} = $ Degredation of switch ($s^{-1}$)

$k_{cdeg} = $ Degredation of output ($s^{-1}$)

$k_{on_2} = $ costimulatory receptor-ligand association rate ($\mu m^{2}s^{-1}$)

$k_{off_2} = $ costimulatory receptor-ligand dissociation rate ($s^{-1}$)

$k_{syn_2} = $ basal costimulatory receptor upregulation ($s^{-1}\mu m^{-2}$)

$k_{basal_2} = $ basal (ligand-independent) costimulatory receptor down-regulation ($s^{-1}$)

$k_{down_2} = $ ligand-induced costimulatory receptor down-regulation ($s^{-1}$)

$k_{cd2a} = $ CD2 costimulation (proximal)

$k_{cd2b} = $ CD2 costimulation (distal)

$k_{icama} = $ LFA-1 costimulation (proximal)

$k_{icamb} = $ LFA-1 costimulation (distal)

$k_{28} = $ CD28 costimulation (distal)

Cstar

Cstar2

time_delay - 41BB upregulation delay after T cell activation




