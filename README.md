# Opemm_Turtorial

This tutorial is used to simulate a protein-ligand system by openmm.  
I've seen a few examples of simulating protein-ligand complexes by openmm before this.  
The problem is the ligand force field is not included in openmm, which needs to be created manually.  
So in this tutorial, I use the gromacs topol file(`.itp`) to generate the ligand xml force field file (which is supported by openmm).    
After docking, you can generate the ligand top by acpype or sobtop. And then use `itp2xml.itp` to generate xml file.  

The itp2xml directory contains two python script  
- `itp2xml.py`
- `simulation.py`

## `itp2xml.py`
`itp2xml.py` is used to convert gromacs ligand itp file to xml file (openmm supported force field format).   
This script now only support one ligand.  
Usage:
```python
python itp2xml.py xxx.itp xxx.xml
```
- `xxx.itp`: gromacs format itp
- `xxx.xml`: set an output xml filename

## `simulation.py`
`simulation.py` is used to run simulation with openmm (you should ensure that openmm is installed in your workstation).   
Before running, you should open the script to change the protein and ligand filename.     
Also, simulation paramters can be changed at the beginning of the file as you wish.

protein pdbfile, ligand pdbfile, generated ligand xml file (by `itp2xml.py`) are nedded. 
Usage:
```python
python simulation.py
```
