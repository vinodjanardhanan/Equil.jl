module Equil

using LightXML, Printf
using LinearAlgebra
using RxnHelperUtils, IdealGas, NonlinearSolve, NLsolve

include("Constants.jl")

export equilibrate

"""
This function calculates equilibrium composition based on initial conditions
    read from an input filter
#   Usage
    equilibrate(input_file::AbstractString, lib_dir::AbstractString)
-   input_file ; input xml file, which specifies T, p, species list, molefracs
-   lib_dir : the folder in which therm.dat file exists

Alternatively the same be calculated based on user input, wherein the species composition is
passed as a Dictionary
#   Usage
    equilibrate(T::Float64, p::Float64, species_comp::Dict{String, Float64}, thermo_obj)
-   T : Temperature in K 
-   p : Pressure in Pa 
-   thermo_obj : Thermodynamic object created using IdealGas.create_thermo 
-   species_comp : Dictionary of species composition {String, Float64}

Alternatively the same be calculated based on user input where in the species list and molefractions are
passed in as separate arrays 
#   Usage
    equilibrate(T, p, thermo_obj, mole_fracs, gasphase)
-   T : Temperature in K 
-   p : Pressure in Pa 
-   thermo_obj : Thermodynamic object created using IdealGas.create_thermo 
-   mole_fracs : species mole fractions 
-   gasphase : list of gasphase species 

"""
function equilibrate(input_file::AbstractString, lib_dir::AbstractString)
    #lib directory    
    xmldoc = parse_file(input_file)
    xmlroot = root(xmldoc)
    

    gasphase = get_collection_from_xml(xmlroot,"gasphase")
    thermo_file = get_path(lib_dir, "therm.dat")
    
    thermo_obj = IdealGas.create_thermo(gasphase, thermo_file)
    
    mole_fracs  = get_molefraction_from_xml(xmlroot,thermo_obj.molwt,gasphase)
    local T = get_value_from_xml(xmlroot,"T")
    local p = get_value_from_xml(xmlroot,"p")


    gasphase, moles, n_equil, mole_frac_final = equilibrate(T,p,thermo_obj, mole_fracs, gasphase)
    # println("species moles = ", n_equil)
    # println("total moles = ", sum(n_equil))
    # println("mole fractions = ", mole_frac_final)


    println("\nInitial condition:\n")
    println("Species \t moles \t\t molefraction")
    for k in eachindex(gasphase)
        @printf("%10s \t %.4e \t %.4e \n", gasphase[k], moles[k], mole_fracs[k])
    end
    eq_stream = open(output_file(input_file,"ch_equil.csv"),"w")
    write_csv(eq_stream, ["Species", "moles", "molefracs"])
    println("\nEquilibrium composition @ T= $T K and p=$p Pa\n")
    println("Species \t moles \t\t molefraction")
    for k in eachindex(gasphase)
        @printf("%10s \t %.4e \t %.4e\n", gasphase[k], n_equil[k], mole_frac_final[k])
        write_csv(eq_stream, [gasphase[k], n_equil[k], mole_frac_final[k]])
    end
    close(eq_stream)

    return Symbol("Success")
    
end

#=
This function calculates equilibrium compostion from user input
#   Usage
    equilibrate(T::Float64, p::Float64, species_comp::Dict{String, Float64}, thermo_obj)
-   T : Temperature in K 
-   p : Pressure in Pa 
-   thermo_obj : Thermodynamic object created using IdealGas.create_thermo 
-   species_comp : Dictionary of species composition {String, Float64}
=#
function equilibrate(T::Float64, p::Float64, thermo_obj, species_comp::Dict{String, Float64})
    gasphase = collect(keys(species_comp))
    molefracs = collect(values(species_comp))
    gasphase, moles, n_equil, mole_frac_final = equilibrate(T,p,thermo_obj, molefracs, gasphase)
    return gasphase, moles, n_equil, mole_frac_final
end

#=
This function calculates equilibrium compostion from user input
#   Usage
    equilibrate(T, p, thermo_obj, mole_fracs, gasphase)
-   T : Temperature in K 
-   p : Pressure in Pa 
-   thermo_obj : Thermodynamic object created using IdealGas.create_thermo 
-   mole_fracs : species mole fractions 
-   gasphase : list of gasphase species 
=#
function equilibrate(T, p, thermo_obj, mole_fracs, gasphase)

    n_total = 1.0
    moles = n_total .* mole_fracs
    
    H_all = IdealGas.H_all(thermo_obj,T)
    S_all = IdealGas.S_all(thermo_obj,T)
    G_all = H_all - T*S_all        

    G_all_by_RT = G_all ./ (R*T)

    

    # Find the elements in the system 
    elements = Array{String,1}()
    for i in thermo_obj.thermo_all
       append!(elements,collect(keys(i.composition)))
    end
    unique!(elements)

    elements_vector = zeros(length(elements))
    for i in thermo_obj.thermo_all
        for (k,v) in i.composition            
            j = get_index(String(k),elements)            
            elements_vector[j] += v* moles[get_index(i.name,gasphase)]
        end
    end

    q_species = (p_std/p) * exp.(-G_all_by_RT)
    
    #Formula coefficient matrix 
    ns = length(gasphase)
    ne = length(elements)
    B = zeros(Float64, ns, ne)
    for i in eachindex(thermo_obj.thermo_all)
        for (k,v) in thermo_obj.thermo_all[i].composition
            j = get_index(String(k),elements)
            B[i,j] = v
        end
    end
    
    # Species moles 
    λ0 = zeros(length(elements))
    # spelpt =  ones(length(thermo_obj.thermo_all))
    
    params = (
        B = B,
        b = elements_vector,
        q_species = q_species,
        thermo_obj = thermo_obj,
        elements = elements,
        gasphase = gasphase
    )

    function residual!(du, λ, p)
        B = p[:B]
        b = p[:b]
        thermo_obj = p[:thermo_obj]
        q_species = p[:q_species]
        elements = p[:elements]
        gasphase = p[:gasphase]

        # spelpt: same type as q_species (so it can hold Duals during AD)
        spelpt = similar(q_species, eltype(λ))


        # compute n_i = q_i * exp(sum_j a_ij * λ_j)
        for sp in thermo_obj.thermo_all
            eλ = zero(eltype(λ))   # AD-safe accumulator
            for (k, v) in sp.composition
                idx = get_index(String(k), elements)
                eλ += v * λ[idx]
            end
            sp_index = get_index(sp.name, gasphase)
            spelpt[sp_index] = q_species[sp_index] * exp(eλ)
        end

        # element residuals: B' * n - b  (B is species × elements)
        du .= B' * spelpt .- b

        return nothing

    end

    prob = NonlinearProblem(residual!, λ0, params)
    sol = solve(prob, NewtonRaphson())

    λ_final = Array(sol.u)   # convert to concrete Float64 array
    # println("Lagrange multipliers = ", λ_final)
    n_equil = [ q_species[i] * exp(sum(B[i,j] * λ_final[j] for j in 1:ne)) for i in 1:ns ]
    mole_frac_final = n_equil ./ sum(n_equil)
    return gasphase, moles, n_equil, mole_frac_final
end


end
