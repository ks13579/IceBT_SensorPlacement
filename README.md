# Ice borehole thermometry: Sensor placement using greedy optimal sampling

This repository provides the data and code supporting the research article **"Ice Borehole Thermometry: Sensor Placement Using Greedy Optimal Sampling"** by K. Shaju, T. Laepple, N. Hrisch, and P. Zaspel.



## Description

This repository contains the following directories and files:

- `Forward_Model/` – For data generation.

- `Sensor placements/` – Code, data, and results.

- `Plots.ipynb`  – Steps to reproduce article plots.

- `requirements.txt` – Required Python libraries.

- `README.md` – Project overview and instructions (that's me!).


### Data generation

The data used for sensor placement methods consist of borehole temperature–depth profiles simulated using a heat transfer model based on various possible surface temperature time series. The code for the heat transfer model, along with the surface temperature time series data, is organized in the `Forward_Model/` directory.

We consider two sites to evaluate our results: EPICA-EDML and GRIP. The `Forward_Model/` directory includes the following subdirectories and files:

* `data/` contains the possible surface temperature time series for both sites.

* `src/` holds the code and relevant parameters for the heat transfer model.

* `EDML_borehole_simulation.ipynb` demonstrates how to simulate the borehole temperature profile for EPICA-EDML.

* `GRIP_borehole_simulation.ipynb` demonstrates the borehole profile simulation for GRIP.

The borehole temperature–depth profiles thus simulated are stored under `Sensor placements/data_borehole_simulations/`, as described in the Code and Data section below.

### Code, data and results

The data generated and the code used to produce the results presented in the article are organized in the `Sensor placements/` directory. This directory includes the following subdirectories and files:

* `code/` contains Python scripts (.py) used to generate the results (i.e., the plots shown in the article). These scripts utilize the data in `data_borehole_simulations/` and store output files in the `output/` directory.

* `data_borehole_simulations/` holds the necessary data, including borehole temperature–depth profiles and other relevant files.

* `output/` contains the generated results, which are used to create the figures in the article.

* `Sensor placement using greedy optimal sampling.ipynb` provides a step-by-step explanation of the greedy optimal sensor placement method, which is the core focus of the research article.

### Plots

The procedures for reproducing the article's plots and locating the corresponding results are detailed in `Plots.ipynb`.

### Requirements

The `requirements.txt` lists the Python libraries required to run the code.

### Instructions

`README.md` - That’s me 🙂 Your guide to the repository, with an overview and usage instructions.



## Software dependencies and installation

This project requires Python 3 (tested with Python 3.10) and the packages listed in `requirements.txt`.

To be able to execute the notebooks, install the Jupyter notebook for example:
```bash
`pip install notebook==6.5.2` 

