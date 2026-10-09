module Equil

using LightXML
using Printf
using LinearAlgebra
using RxnHelperUtils
using IdealGas
using NonlinearSolve
using NLsolve

include("Constants.jl")

export equilibrate


# ============================================================
# equilibrate from XML input file
# ============================================================

"""
    equilibrate(input_file::AbstractString, thermo_file::AbstractString)

Calculate equilibrium composition from an XML input file.

The XML file specifies:
    - T       : temperature [K]
    - p       : pressure [Pa]
    - gasphase species
    - initial mole fractions

`thermo_file` is the path to the thermodynamic database file.
"""
function equilibrate(input_file::AbstractString, thermo_file::AbstractString)

    # --------------------------------------------------------
    # Read XML file
    # --------------------------------------------------------

    xmldoc = parse_file(input_file)
    xmlroot = root(xmldoc)

    gasphase = get_collection_from_xml(xmlroot, "gasphase")


    thermo_obj = IdealGas.create_thermo( gasphase, thermo_file )

    mole_fracs = get_molefraction_from_xml( xmlroot, thermo_obj.molwt, gasphase )

    local T = get_value_from_xml(xmlroot, "T")
    local p = get_value_from_xml(xmlroot, "p")


    # --------------------------------------------------------
    # Equilibrium calculation
    # --------------------------------------------------------

    gasphase,
    moles,
    n_equil,
    mole_frac_final = equilibrate( T, p, thermo_obj, mole_fracs, gasphase )


    # --------------------------------------------------------
    # Print initial condition
    # --------------------------------------------------------

    println("\nInitial condition:\n")

    println("Species \t moles \t\t molefraction")

    for k in eachindex(gasphase)

        @printf( "%10s \t %.4e \t %.4e \n", gasphase[k], moles[k], mole_fracs[k] )

    end


    # --------------------------------------------------------
    # Write equilibrium results
    # --------------------------------------------------------
    eq_stream = open( output_file(input_file, "ch_equil.csv"), "w")

    write_csv( eq_stream, ["Species", "moles", "molefracs"] )
    println( "\nEquilibrium composition @ T= $T K and p=$p Pa\n" )
    println("Species \t moles \t\t molefraction")

    for k in eachindex(gasphase)
        @printf( "%10s \t %.4e \t %.4e\n", gasphase[k], n_equil[k], mole_frac_final[k] )
        write_csv(eq_stream, [gasphase[k], n_equil[k], mole_frac_final[k] ] )
    end

    close(eq_stream)
    return Symbol("Success")

end


# ============================================================
# equilibrate from Dictionary
# ============================================================

"""
    equilibrate(
        T::Float64,
        p::Float64,
        thermo_obj,
        species_comp::Dict{String,Float64}
    )

Calculate equilibrium composition from a dictionary:

    Dict(
        "CH4" => 0.5,
        "O2"  => 0.5
    )
"""
function equilibrate( T::Float64, p::Float64,  thermo_obj, species_comp::Dict{String,Float64})

    gasphase = collect(keys(species_comp))
    molefracs = collect(values(species_comp))
    gasphase, moles, n_equil, mole_frac_final = equilibrate( T, p,thermo_obj, molefracs, gasphase )
    return ( gasphase, moles, n_equil, mole_frac_final )

end


# ============================================================
# Main equilibrium calculation
# ============================================================

"""
    equilibrate(
        T,
        p,
        thermo_obj,
        mole_fracs,
        gasphase
    )

Calculate equilibrium composition using Gibbs free-energy
minimization subject to elemental conservation.

The calculation uses the Lagrange multiplier formulation:

    n_i = q_i exp(sum_j B[i,j] λ_j)

where

    B[i,j] = number of atoms of element j in species i

and the nonlinear equations are

    B' * n - b = 0

where `b` contains the initial elemental inventories.
"""
function equilibrate( T,p, thermo_obj,  mole_fracs, gasphase )

    # ========================================================
    # Basic information
    # ========================================================

    n_total = 1.0
    ns = length(gasphase)
    moles = n_total .* mole_fracs


    # ========================================================
    # IMPORTANT:
    #
    # Construct all thermodynamic information in GASPHASE
    # ORDER.
    #
    # Do NOT use thermo_obj.thermo_all indices directly for
    # the equilibrium vectors/matrices.
    # ========================================================

    thermo_names = [sp.name for sp in thermo_obj.thermo_all]


    # --------------------------------------------------------
    # Find thermodynamic object for each gasphase species
    # --------------------------------------------------------

    thermo_species = Vector{Any}(undef, ns)

    for i in 1:ns
        idx = get_index( gasphase[i], thermo_names )
        thermo_species[i] = thermo_obj.thermo_all[idx]
    end


    # ========================================================
    # Diagnostic species mapping
    # ========================================================

    # Uncomment this if you want to verify the mapping.
    #
    # println("\nSpecies mapping:")
    #
    # for i in 1:ns
    #     println(
    #         "gasphase[$i] = ",
    #         gasphase[i],
    #         "   thermo species = ",
    #         thermo_species[i].name
    #     )
    # end


    # ========================================================
    # Thermodynamic properties
    # ========================================================

    H_all = IdealGas.H_all( thermo_obj, T )
    S_all = IdealGas.S_all( thermo_obj, T )
    G_all = H_all .- T .* S_all


    # --------------------------------------------------------
    # Gibbs energy in GASPHASE ORDER
    # --------------------------------------------------------

    G_gas = zeros(Float64, ns)
    for i in 1:ns
        idx = get_index( gasphase[i], thermo_names )
        G_gas[i] = G_all[idx]
    end


    # Gibbs free energy divided by RT

    G_gas_by_RT = G_gas ./ (R * T)


    # ========================================================
    # Collect elements
    # ========================================================

    elements = String[]
    for sp in thermo_species
        append!( elements, String.(collect(keys(sp.composition))) )
    end

    unique!(elements)
    ne = length(elements)


    # ========================================================
    # Element inventory of initial mixture
    #
    # b[j] = total number of atoms of element j
    # ========================================================

    b = zeros(Float64, ne)
    for i in 1:ns
        sp = thermo_species[i]
        for (element, amount) in sp.composition
            j = get_index( String(element), elements )
            b[j] += amount * moles[i]
        end
    end


    # ========================================================
    # Formula coefficient matrix
    #
    # B[i,j] =
    # number of atoms of element j in species i
    #
    # IMPORTANT:
    # rows are in GASPHASE ORDER
    # ========================================================

    B = zeros(Float64, ns, ne)

    for i in 1:ns
        sp = thermo_species[i]
        for (element, amount) in sp.composition
            j = get_index( String(element), elements )
            B[i, j] = amount
        end

    end


    # ========================================================
    # Calculate q_i
    #
    # For ideal gases:
    #
    # n_i = q_i exp(sum_j B[i,j] λ_j)
    #
    # q_i = n_total * (p0/p) * exp(-G_i/RT)
    #
    # ========================================================

    p0 = 101325.0
    q_species = zeros(Float64, ns)
    for i in 1:ns
        q_species[i] = n_total * (p0 / p) * exp(-G_gas_by_RT[i])
    end


    # ========================================================
    # Diagnostic checks BEFORE solving
    # ========================================================

    # Initial elemental balance

    initial_element_balance = B' * moles
    balance_error = initial_element_balance .- b
    # This should be approximately zero.
    if maximum(abs.(balance_error)) > 1e-10
        println( "\nWARNING: Initial elemental balance is not zero." )
        println( "B' * moles = ", initial_element_balance )
        println( "b          = ", b )
        println( "error      = ", balance_error )
    end


    # ========================================================
    # Nonlinear system in Lagrange multipliers
    #
    # B' * n - b = 0
    #
    # n_i = q_i exp(B[i,:]' λ)
    # ========================================================

    λ0 = zeros(Float64, ne)


    params = ( B = B, b = b, q_species = q_species )


    # ========================================================
    # Residual function
    # ========================================================

    function residual!( du, λ, p )
        B = p.B
        b = p.b
        q_species = p.q_species

        # ----------------------------------------------------
        # Species mole numbers
        # ----------------------------------------------------

        n = similar( q_species, eltype(λ) )
        for i in 1:length(q_species)
            eλ = zero(eltype(λ))
            for j in 1:size(B, 2)
                eλ += B[i, j] * λ[j]
            end
            n[i] = q_species[i] * exp(eλ)
        end


        # ----------------------------------------------------
        # Element balance residual
        # ----------------------------------------------------

        du .= B' * n .- b
        return nothing

    end


    # ========================================================
    # Nonlinear problem
    # ========================================================

    prob = NonlinearProblem( residual!, λ0, params )

    # ========================================================
    # Solve
    # ========================================================

    sol = solve( prob, TrustRegion(); abstol = 1e-10, reltol = 1e-10, maxiters = 1000 )


    # ========================================================
    # Check solver status
    # ========================================================

    if !SciMLBase.successful_retcode(sol.retcode)
        println( "\nWARNING: Equilibrium solver did not converge." )
        println( "retcode = ", sol.retcode )
    end


    # ========================================================
    # Final Lagrange multipliers
    # ========================================================

    λ_final = Array(sol.u)


    # ========================================================
    # Back-calculate equilibrium mole numbers
    # ========================================================

    n_equil = zeros(Float64, ns)
    for i in 1:ns
        exponent = zero(Float64)
        for j in 1:ne
            exponent += B[i, j] * λ_final[j]
        end
        n_equil[i] = q_species[i] * exp(exponent)

    end


    # ========================================================
    # Equilibrium mole fractions
    # ========================================================

    mole_frac_final = n_equil ./ sum(n_equil)


    # ========================================================
    # Final elemental balance check
    # ========================================================

    final_element_balance =
        B' * n_equil

    final_balance_error =
        final_element_balance .- b

    if (abs.(final_balance_error) .> 1e-10) |> any
        println( "\nWARNING: Final elemental balance is not zero." )
        println( "B' * n_equil = ", final_element_balance )
        println( "b            = ", b )
        println( "error        = ", final_balance_error )
    end
    

    # Uncomment for detailed diagnostics:
    #
    # println("Elements = ", elements)
    # println("B = ")
    # display(B)
    #
    # println("Initial elemental balance = ",
    #         B' * moles)
    #
    # println("Required elemental balance = ",
    #         b)
    #
    # println("Final elemental balance = ",
    #         B' * n_equil)
    #
    # println("Equilibrium residual = ",
    #         final_balance_error)


    # ========================================================
    # Return
    # ========================================================

    return ( gasphase, moles, n_equil, mole_frac_final )

end


end