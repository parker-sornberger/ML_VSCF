# Code to reproduce results for "paper title"

## ParseFreqCalc.py
Driver code which allows harmonic frequencies, normal mode eigenvectors, mode energies, etc to be parsed from a valid CRYSTAL frequency output file. 
Needed data from user:
1. Frequency output file

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
