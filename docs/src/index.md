```@meta
CurrentModule = Equil
```

# Equil
Equil is a program for the calculation of equilibrium composition. There are different methods for calculating the equilibrium composition of a reacting mixture. For a system that consists of $N$ species and $K$ elements ($N > K$), the elements are conserved, and the classical approach involves the solution of $N+K$ non-linear equations. This approach calculates the moles of different species ($N$ numbers) at equilibrium and $K$ number of Lagrangian multipliers.  An alternate method solves $N-K$ non-linear equations involving equilibrium constants and $K$ equations describing component activities. The latter approach is adopted in this program.


Documentation for [Equil](https://github.com/vinodjanardhanan/Equil.jl).

## Installation
To install the package, use the following commands in the julia REPL
```julia
julia> using Pkg
julia> Pkg.add("Equil")
```

## General interfaces
```@index
```

```@autodocs
Modules = [Equil]
```

## Executing the code
On the Julia REPL 
```julia
julia>using Equil
julia>equil("equil.xml","../lib/")
```

## Input file
The method takes *file\_path* as the argument. The file_path points to the input XML file. The structure of the XML input file is shown below.

```
<?xml version="1.0" encoding="ISO-8859-1"?>
<equil>
    <gasphase>CH4 H2 CO CO2 H2O O2 N2 </gasphase>
    <molefractions>CH4=0.6, H2O=0.1, CO2=0.2, N2=0.1</molefractions>
    <T>1073.15</T>
    <p>1e5</p>
</equil>
```
The meaning of different tags is specified below.

- <equil> : The root XML tag for equilibrate
- <gasphase> : list of gas-phase species. The species names must be separated by white spaces or tab
- <molefractions> : Initial mole fraction of the gas-phase species. Instead of mole fractions, mass fractions may also be specified. In that case, the tag must be <massfractions>. You must ensure that the sum of mass or mole fractions specified is unity. There are no internal checks to ascertain this. 
- <T>: operating temperature in K
- <p>: initial pressure in Pa

## Output
The code generates only screen output. An example output is shown below.

```
Inititial condition:

Species 	 moles 		 molefraction
       CH4 	 6.7244e+00 	 6.0000e-01 
        H2 	 0.0000e+00 	 0.0000e+00 
        CO 	 0.0000e+00 	 0.0000e+00 
       CO2 	 2.2415e+00 	 2.0000e-01 
       H2O 	 1.1207e+00 	 1.0000e-01 
        O2 	 0.0000e+00 	 0.0000e+00 
        N2 	 1.1207e+00 	 1.0000e-01

Equilibrium composition @ T= 1073.15 K and p=100000.0 Pa

Species          moles           molefraction
       CH4       3.3912e+00      1.8973e-01
       CO2       1.2649e-02      7.0766e-04
       H2O       1.6324e-02      9.1329e-04
        CO       5.5621e+00      3.1118e-01
        O2       3.2984e-23      1.8454e-24
        H2       7.7709e+00      4.3476e-01
        N2       1.1207e+00      6.2703e-02
```        

In certain cases, you may want to calculate the equilibrum composition from within another program.
For such cases Equil package provides another interface *equilibrate*

```julia
julia>using Equil
julia>using IdealGas
julia>gasphase = ["CH4","H2", "CO", "CO2", "H2O", "O2" ,"N2"]
julia>mole_fracs = [0.6, 0.0, 0.0, 0.2, 0.1, 0.0, 0.1]
julia>thermo_file = get_path(lib_dir, "therm.dat")
julia>thermo_obj = create_thermo(gasphase, thermo_file)        
julia>T = 1073.15
julia>p = 1e5
julia>gasphase_in, moles_all, eq_molefracs  = equilibrate(T,p, thermo_obj, mole_fracs, gasphase)        
```
Instead of passing the species list and the mole fractions separately, you may pass a Dictionary of species names and the molefractions a described below


```julia
julia>using IdealGas
julia>using Equil
julia>species_comp = Dict("CH4"=>0.6,"H2"=>0.0, "CO"=>0.0, "CO2"=>0.2, "H2O"=>0.1, "O2"=>0.0 ,"N2"=>0.1)        
julia>thermo_file = get_path(lib_dir, "therm.dat")
julia>thermo_obj = create_thermo(collect(keys(species_comp)), thermo_file)        
julia>T = 1073.15
julia>p = 1e5
julia>gasphase_in, moles_all, eq_molefracs  = equilibrate(T,p, thermo_obj, species_comp)                
```

## Input file download
The xml input file and the *lib* directory containig other required input files may be downloaded from [here](https://github.com/vinodjanardhanan/Equil.jl/tree/main/test).
