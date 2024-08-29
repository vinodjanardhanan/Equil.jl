var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = Equil","category":"page"},{"location":"#Equil","page":"Home","title":"Equil","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Equil is a program for the calculation of equilibrium composition. There are different methods for calculating the equilibrium composition of a reacting mixture. For a system that consists of N species and K elements (N  K), the elements are conserved, and the classical approach involves the solution of N+K non-linear equations. This approach calculates the moles of different species (N numbers) at equilibrium and K number of Lagrangian multipliers.  An alternate method solves N-K non-linear equations involving equilibrium constants and K equations describing component activities. The latter approach is adopted in this program.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Documentation for Equil.","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"To install the package, use the following commands in the julia REPL","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> using Pkg\njulia> Pkg.add(\"Equil\")","category":"page"},{"location":"#General-interfaces","page":"Home","title":"General interfaces","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [Equil]","category":"page"},{"location":"#Equil.equilibrate-Tuple{AbstractString, AbstractString}","page":"Home","title":"Equil.equilibrate","text":"This function calculates equilibrium composition based on initial conditions     read from an input filter\n\nUsage\n\nequilibrate(input_file::AbstractString, lib_dir::AbstractString)\n\ninput_file ; input xml file, which specifies T, p, species list, molefracs\nlib_dir : the folder in which therm.dat file exists\n\nAlternatively the same be calculated based on user input, wherein the species composition is passed as a Dictionary\n\nUsage\n\nequilibrate(T::Float64, p::Float64, species_comp::Dict{String, Float64}, thermo_obj)\n\nT : Temperature in K \np : Pressure in Pa \nthermoobj : Thermodynamic object created using IdealGas.createthermo \nspecies_comp : Dictionary of species composition {String, Float64}\n\nAlternatively the same be calculated based on user input where in the species list and molefractions are passed in as separate arrays \n\nUsage\n\nequilibrate(T, p, thermo_obj, mole_fracs, gasphase)\n\nT : Temperature in K \np : Pressure in Pa \nthermoobj : Thermodynamic object created using IdealGas.createthermo \nmole_fracs : species mole fractions \ngasphase : list of gasphase species \n\n\n\n\n\n","category":"method"},{"location":"#Executing-the-code","page":"Home","title":"Executing the code","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"In order to calculate the equilibrium composition based on input file, do the following. On the Julia REPL ","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia>using Equil\njulia>equilibrate(\"equil.xml\",\"../lib/\")","category":"page"},{"location":"#Input-file","page":"Home","title":"Input file","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The method takes file_path as the argument. The file_path points to the input XML file. The structure of the XML input file is shown below.","category":"page"},{"location":"","page":"Home","title":"Home","text":"<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n<equil>\n    <gasphase>CH4 H2 CO CO2 H2O O2 N2 </gasphase>\n    <molefractions>CH4=0.6, H2O=0.1, CO2=0.2, N2=0.1</molefractions>\n    <T>1073.15</T>\n    <p>1e5</p>\n</equil>","category":"page"},{"location":"","page":"Home","title":"Home","text":"The meaning of different tags is specified below.","category":"page"},{"location":"","page":"Home","title":"Home","text":"<equil> : The root XML tag for equilibrate\n<gasphase> : list of gas-phase species. The species names must be separated by white spaces or tab\n<molefractions> : Initial mole fraction of the gas-phase species. Instead of mole fractions, mass fractions may also be specified. In that case, the tag must be <massfractions>. You must ensure that the sum of mass or mole fractions specified is unity. There are no internal checks to ascertain this. \n<T>: operating temperature in K\n<p>: initial pressure in Pa","category":"page"},{"location":"#Output","page":"Home","title":"Output","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The code generates screen output as well as a CSV file (chequil.csv). The *chequil.csv* file will be located in the same directory where the input xml file is located. An example output is shown below.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Inititial condition:\n\nSpecies \t moles \t\t molefraction\n       CH4 \t 6.7244e+00 \t 6.0000e-01 \n        H2 \t 0.0000e+00 \t 0.0000e+00 \n        CO \t 0.0000e+00 \t 0.0000e+00 \n       CO2 \t 2.2415e+00 \t 2.0000e-01 \n       H2O \t 1.1207e+00 \t 1.0000e-01 \n        O2 \t 0.0000e+00 \t 0.0000e+00 \n        N2 \t 1.1207e+00 \t 1.0000e-01\n\nEquilibrium composition @ T= 1073.15 K and p=100000.0 Pa\n\nSpecies          moles           molefraction\n       CH4       3.3912e+00      1.8973e-01\n       CO2       1.2649e-02      7.0766e-04\n       H2O       1.6324e-02      9.1329e-04\n        CO       5.5621e+00      3.1118e-01\n        O2       3.2984e-23      1.8454e-24\n        H2       7.7709e+00      4.3476e-01\n        N2       1.1207e+00      6.2703e-02","category":"page"},{"location":"","page":"Home","title":"Home","text":"In certain cases, you may want to calculate the equilibrum composition from within another program. For such cases Equil package provides another interface equilibrate","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia>using Equil\njulia>using IdealGas\njulia>gasphase = [\"CH4\",\"H2\", \"CO\", \"CO2\", \"H2O\", \"O2\" ,\"N2\"]\njulia>mole_fracs = [0.6, 0.0, 0.0, 0.2, 0.1, 0.0, 0.1]\njulia>thermo_file = get_path(lib_dir, \"therm.dat\")\njulia>thermo_obj = create_thermo(gasphase, thermo_file)        \njulia>T = 1073.15\njulia>p = 1e5\njulia>gasphase_in, moles_all, eq_molefracs  = equilibrate(T,p, thermo_obj, mole_fracs, gasphase)        ","category":"page"},{"location":"","page":"Home","title":"Home","text":"Instead of passing the species list and the mole fractions separately, you may pass a Dictionary of species names and the molefractions a described below.","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia>using IdealGas\njulia>using Equil\njulia>species_comp = Dict(\"CH4\"=>0.6,\"H2\"=>0.0, \"CO\"=>0.0, \"CO2\"=>0.2, \"H2O\"=>0.1, \"O2\"=>0.0 ,\"N2\"=>0.1)        \njulia>thermo_file = get_path(lib_dir, \"therm.dat\")\njulia>thermo_obj = create_thermo(collect(keys(species_comp)), thermo_file)        \njulia>T = 1073.15\njulia>p = 1e5\njulia>gasphase_in, moles_all, eq_molefracs  = equilibrate(T,p, thermo_obj, species_comp)                ","category":"page"},{"location":"","page":"Home","title":"Home","text":"The ch_equil.csv file will not be generated for call from the REPL or for an interface call.","category":"page"},{"location":"#Input-file-download","page":"Home","title":"Input file download","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The xml input file and the lib directory containig other required input files may be downloaded from here.","category":"page"}]
}
