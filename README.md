# Code to reproduce results for "paper title"

## ParseFreqCalc.py
Driver code which allows harmonic frequencies, normal mode eigenvectors, mode energies, etc to be parsed from a valid CRYSTAL frequency output file. 
Needed data from user:
1. Frequency output file
```python
import os
from ParseFreqCalc import ParseFreq
path_to_freq_file = "some_path/"
frequency_file = "some_file.out"
freq_calc = ParseFreq(os.join(path_to_freq_file, frequency_file), )

#get modes from the file
modes = freq_calc.modes
#returns a dictionary of modes with integer keys, starting at 1 as convention in crystal
info_for_mode_10 = modes[10]
#this returns a dictionary with the following keys 'EigenValue', 'WaveNumber', 'THz', 'ActiveIR', 'Intensity', 'ActiveRaman', 'EigenVector'
energy_in_wavnum = info_for_mode_10['WaveNumber']
mode_10_eigenvector = info_for_mode_10['EigenVector'] #NumPy array with shape (Atoms in Unit Cell X 3)

```

## ParseCalcVSCF.py 
Easily obtains coupled and uncoupled 3rd and 4th order force constants from a CRYSTAL VSCF output file. Can also write VIBPOT.DAT files from the parsed force constants. 
Needed data from user:
1. VSCF output file

## ScanEnergiesWithANI.py
Given an eigenvector for a normal mode, will generate a scan of energies as the equilibrium geometry of a crystal is perturbed by scaled multiples of that normal mode. Also works for pairs of normal modes. 
Needed data from user:
1. CIF of equilibirium crystal
2. Frequency output file
3. Path to folder to dump CIF files created intermediately (removed as they are generated)

## GenerateInputsForCRYSTAL.py 
Generates the input files for single point energies, frequency calculations, anharmonic force constant calculations, and VSCF calculations. 
Needed data from user:
1. CIF file for crystal including fractional coordinates, symmetry, cell parameters, etc
2. Basis set
3. Keywords for frequency / anharmonic calculations
4. Optional, shrink values, tolerance values 

## CalcForceConstants.py 
Uses the output for a frequency calculation from CRYSTAL to generate higher order force constants. 
Needed data from user:
1. Frequency output file
2. Path to folder to dump CIF files created intermediately (removed as they are generated)

# Calculating force constants
We use this paper "DOI for paper" following "Scheme 1". Here, the $i$ in a coupled force constant such as $\eta_{iij}$ is associated with the eigenvector $\nu$ and the $j$ corresponds to the eigenvector $\mu$.  (make sure we reference that we are quoting ourselves to not seem plagerizing lol)

The scaling factors $s_i$ or $s_j$ are calcated using the following formula from their energy in $cm^{-1}$.  $m_e = 1822.8848$ (units) $H_{eV} = 27.211386024367243 \frac{eV}{H}$ $\omega_{eV} = 8065.73 \frac{cm^{-1}}{eV}$ $B_A =0.529177 \frac{Bohr}{Angstrom} $  



### Uncoupled 3rd and 4th order force constants 

$\eta_{iii}=\frac{1}{2s_i^3}(-E_{-2}+2E_{-1}-2E_1+E_2)$

$\eta_{iiii}=\frac{1}{s_i^4}(6E_0+E_{-2}-4E_{-1}-4E_1+E_2)$

### Coupled 3rd and 4th order force constants 

$\eta_{iij}=\frac{1}{2s_i^2s_j}(2E_{0,-1}-2E_{0,1}-E_{-1,-1}+E_{-1,1}+E_{1,1})$

$\eta_{iiij}=\frac{1}{4s_i^3s_j}(E_{-2,-1}-E_{-2,1}-2E_{-1,-1}+2E_{-1,1}+2E_{1,-1}-2E_{1,1}-E_{2,-1}+E_{2,1})$

$\eta_{iijj}=\frac{1}{s_i^2s_j^2}(4E_0-2E_{0,-1}-2E_{0,1}-2E_{-1,0}+E_{-1,-1}-E_{-1,1}-2E_{1,0}+E_{1,-1}+E_{1,1})$


# Drawbacks and future directions
Presently, most of these programs require at least a convereged frequency calculation, which for large systems may incur an extensive time cost. However, with the present scheme used in calcualting the higher order force constants, to characterize the coupled constants between all possible modes in a system with $M$ vibrational modes a total of $N$ single point energy calculations is needed where $N = 1 + 4M + 12 \binom{M}{2}$. As the degrees of freedom increase in a system, such as a cyrstal for a material used in something like an organic semiconductor, the total number of points increases exponentially, and the time needed to compute even a fraction of the total possible coupled force constants becomes intractable. While using ANI to calculate these single point energies does not escape the need for their calculation, it still greatly accelerates their computation. Unfortunately, we also require that the force constants as computed by ANI be converted to a VIBPOT.DAT file to allow for the VSCF calculations to be computed. 

However, in future work, we aim to show that using the normal modes as computed by some machine learning potential such as the one from ANI can still produce reasonable energies in the application of VSCF. Likewise, this would allow people to use any real geometry, such as a molecule, rather than only a solid-state system, and compare against the modes computed from any software, including Gaussian. Moreover, we also intend to implement a python-based version of the VSCF code found in CRYSTAL to allow for more facile interaction between the force constants computed by ANI and the transition energies obtained from them. Due to differences in how CRYSTAL computes the eigenvectors and eigenvalues for a given matrix and how these are solved for in LAPACK, the software which powers these routined in NumPy and PyTorch, the transition energies calculated via an equivalent algorithm in python do not line up with ones obtained via CRYSTAL. 

Presently we just aim to show a method to bypass the computational cost incurred with just the CRYSTAL software system, currently the gold-standard for solid-sate DFT, to compute the higher-order force constants. 
