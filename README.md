# AI3: Protein-Ligand Binding Affinity Dataset

## Overview

The **AI3 (AWS-IIITH-Intel-Insilico) Dataset** project aims to create the world's largest protein-ligand binding affinity dataset through a two-stage computational process. 

### Stage 1

In the initial stage, **20,000 bound protein-ligand complexes (PLCs)** have been processed. These PLCs were meticulously downloaded from **RCSB PDB (RCSB-PDB, 2023)**, ensuring each complex includes the protein, ligand, cofactor(s), and crystal water molecules. To accommodate multimers, only subunits within 8Å of the ligand molecule were considered. The missing amino acid residues were modeled using **UCSF Chimera's modeler package (UCSF-Chimera, 2023)**, with protonation states determined via the **H++ webserver (H++, 2023)** at a biological pH of 7.4. The protonation states of ligands and cofactors were ascertained through the **Proteins Plus webserver (ZBH, 2023)**.

### Stage 2

The subsequent stage will incorporate **200,000 unbound or partially bound PLCs**, derived from the Stage 1 structures. For each original PLC, 10 additional structures will be generated to represent higher energy states.

## Methodology

The **Amber14 SB force field** defined protein and peptide atom parameters, with ligand and cofactor atoms detailed by **GAFF2**. Ligand and cofactor atoms were charged using **AM1BCC**. Systems were solvated using Amber's **tleap module**, with a solvation buffer of 10Å from the protein surface, and neutralized with Na+ or Cl- ions.

Energy minimization was conducted over 2000 steps using the steepest descent algorithm in **GROMACS**. This was followed by heating to 300K over 400ps in an NVT ensemble, preparing five systems with unique initial velocity distributions to mitigate uncertainty in binding affinity predictions. The production MD run was then executed under NPT conditions at 300K for 6ns, saving frames every 100ps. Binding free energy calculations utilized the last 4ns of simulation data.

The **molecular mechanics Poisson-Boltzmann surface area (MMPBSA) approach** estimated the binding free energy, treating the solvent as a dielectric continuum. This method, chosen for its efficiency, calculates energetics via molecular mechanics force fields, with solvation components derived from the Poisson-Boltzmann equation and changes in solvent accessible surface area. A single trajectory protocol provided robust and accurate binding affinity estimations from five independent runs.

## Repository Structure

- **`/IIITH`**: Codes for input generation for MD simulation.
- **`/AWS`**: Codes for MD simulation execution and protein-ligand binding free energy calculation using MMPBSA.
- **`/InsilicoMedicine`**: Codes for protein-ligand binding free energy calculation using Alchemical free energy calculation.

## References

- **RCSB-PDB, 2023**
- **UCSF-Chimera, 2023**
- **H++, 2023**
- **ZBH, 2023**

## Acknowledgments

This project is a collaborative effort between **AWS, IIITH, Intel, and Insilico Medicine**, aiming to advance the understanding of protein-ligand interactions through computational simulations.
