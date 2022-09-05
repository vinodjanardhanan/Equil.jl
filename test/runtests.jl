using Equil
using Test
using IdealGas
using RxnHelperUtils

@testset "Equil.jl" begin
    if Sys.isapple() || Sys.islinux()
        lib_dir = "lib/"
    elseif Sys.iswindows()
        lib_dir = "lib\\"
    end
    @testset "Testing ch4 reforming" begin         
        retcode = equil("ch4/equil.xml", lib_dir)
        @test retcode == Symbol("Success")
    end
    @testset "Testing RWGS" begin         
        retcode = equil("rwgs/equil.xml", lib_dir)
        @test retcode == Symbol("Success")
    end
    @testset "Interface call-1" begin
        gasphase = ["CH4","H2", "CO", "CO2", "H2O", "O2" ,"N2"]
        mole_fracs = [0.6, 0.0, 0.0, 0.2, 0.1, 0.0, 0.1]
        thermo_file = get_path(lib_dir, "therm.dat")
        thermo_obj = create_thermo(gasphase, thermo_file)        
        T = 1073.15
        p = 1e5
        gasphase_in, moles_all, eq_molefracs  = equilibrate(T,p, thermo_obj, mole_fracs, gasphase)        
        @test eq_molefracs != mole_fracs
    end
    @testset "Interface call-2" begin
        species_comp = Dict("CH4"=>0.6,"H2"=>0.0, "CO"=>0.0, "CO2"=>0.2, "H2O"=>0.1, "O2"=>0.0 ,"N2"=>0.1)        
        thermo_file = get_path(lib_dir, "therm.dat")
        thermo_obj = create_thermo(collect(keys(species_comp)), thermo_file)        
        T = 1073.15
        p = 1e5
        gasphase_in, moles_all, eq_molefracs  = equilibrate(T,p, thermo_obj, species_comp)                
        @test eq_molefracs != collect(values(species_comp))
    end
end
