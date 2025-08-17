Reconstructing Scalar Dark Energy- A Collider-Based Simulation Study

Overview

This repository contains the code and data for my senior thesis project, "Reconstructing Scalar Dark Energy: A Collider-Based Simulation Study" (UCI, 2025).

The project explores potential dark energy signatures under the quintessence model using simulated proton–proton collisions. Events were generated with MadGraph5, Pythia, and Delphes, and analyzed using Python and ROOT.

The final thesis paper (thesis.pdf) provides the complete research context, methodology, and results.

Repository Structure

/code

Reconstruction.py — Main analysis: reconstructing collision events

Histograms.py — Secondary script: generates histograms of results

/simulation_settings

run_card.dat — Simulation run parameters (event counts, cuts, etc.)

param_card.dat — Physics parameters for the simulated collisions (masses, widths, couplings, etc.)

Reconstructing Scalar Dark Energy- A Collider-Based Simulation Study - Edward Sahyoun.pdf — Full thesis paper

Note: Original .root simulation outputs (several GB) are not included due to size constraints. See thesis pdf for full results and analysis.

Usage

Collision events were originally simulated with MadGraph5, Pythia, and Delphes. The /data directory provides the input cards used to configure these simulations.

Run the main reconstruction script: python code/Reconstruction.py

Generate histograms from the processed data: python code/Histograms.py

(Visualizations and final analysis are included in thesis pdf.)
