module Equil

using LightXML, Printf
using LinearAlgebra, Roots, RowEchelon
using RxnHelperUtils, IdealGas

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

    moles = (p/R/T)*mole_fracs
    
    println("\nInitial condition:\n")
    println("Species \t moles \t\t molefraction")
    for k in eachindex(gasphase)
        @printf("%10s \t %.4e \t %.4e \n", gasphase[k], moles[k], mole_fracs[k])
    end

    gasphase_in, moles_all, molefracs  =  equilibrate(T, p, thermo_obj, mole_fracs, gasphase)


    println("\nEquilibrium composition @ T= $T K and p=$p Pa\n")
    println("Species \t moles \t\t molefraction")
    for k in eachindex(gasphase_in)
        @printf("%10s \t %.4e \t %.4e\n", gasphase_in[k], moles_all[k], molefracs[k])        
    end

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
    equilibrate(T,p,thermo_obj, molefracs, gasphase)
    
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
    all_species_thermo = thermo_obj.thermo_all
    #Total number of moles under the given conditions
    n_total = p/R/T
    # #Number of moles of each component    
    moles = n_total*mole_fracs
    #Order the species according to inlet moles
    order_species(gasphase,moles)
    #Save the ordered list for later use
    gasphase_in = copy(gasphase)
    moles_in = copy(moles)
    

    #reorder  the thermo obj based on the ordered list
    reorder_thermo!(all_species_thermo, gasphase)


    #create elements from the species composition
    elements = Array{String,1}()
    for sp in all_species_thermo        
        append!(elements,collect(keys(sp.composition)))
    end
    unique!(elements)    


    #Create the stoichiometric coefficients for the indepdendent reactions
    stc = stoichiometric_coeffient_matrix(gasphase,moles,elements,all_species_thermo)   
    stc = round.(stc,digits=2)
    #The above method removes the inerts from the list. Therefore store the inlet moles of inerts for final calculation 
    
    inert_species = setdiff(gasphase_in,gasphase)
    inert_moles = Array{Float64,1}()    
    for i in inert_species
        append!(inert_moles,moles_in[get_index(i,gasphase_in)])
    end
    
    # remove the inert species from thermo object 
    for inert in inert_species
        remove_inert_from_thermo!(inert, all_species_thermo)
    end
        

    #Uncomment this for printing the independent reactions
    # print_reactions(stc,gasphase)

    #get the Gibb's free energy of all species to calculate Kp
    H_all = IdealGas.H_all(thermo_obj,T)
    S_all = IdealGas.S_all(thermo_obj,T)
    G_all = H_all - T*S_all        
    
    # #Calculate the Del G for each reaction 
    DelG = Array{Float64,1}()
    for i in 1:size(stc,1)        
        push!(DelG,sum(stc[i,:] .* G_all))
    end    
    molefracs = moles_in/sum(moles_in)
    

    Kp = exp.(-DelG/R/T)
    # Minimize the Gibb's energy of the mixture
    minimize_G!(moles, stc, thermo_obj, Kp, T, p)
    #Account for the inert moles 
    moles_all = copy(moles)
    append!(moles_all,inert_moles)
    molefracs = moles_all/sum(moles_all)
    gasphase_in = copy(gasphase)
    append!(gasphase_in,inert_species)
    
    return gasphase_in, moles_all, molefracs  
    
end


function remove_inert_from_thermo!(inert, thermo_all)
    inert_pos = 0
    for i in eachindex(thermo_all)
        if inert == thermo_all[i].name 
            inert_pos = i
            break
        end
    end
    if inert_pos != 0
        deleteat!(thermo_all, inert_pos)
    end
end

function reorder_thermo!(thermo_all, gasphase)    
    for k in eachindex(gasphase)
        for i in eachindex(thermo_all)
            if gasphase[k] == thermo_all[i].name && k != i 
                thermo_all[k], thermo_all[i] = thermo_all[i], thermo_all[k]
            end
        end
    end
end


#Bubble sort the species with the highest to lower order of moles    
function order_species(gasphase::Array{String}, moles::Array{Float64})
    for i in 1:length(gasphase)-1
        for j in i+1:length(gasphase)
            if moles[i] < moles[j]
                moles[i], moles[j] = moles[j], moles[i]
                gasphase[i], gasphase[j] = gasphase[j], gasphase[i]
            end
        end
    end
end


#=This function creates the stoichiometric coefficients for 
linearly indepdendent reactions based on null space
=#
function stoichiometric_coeffient_matrix(gasphase, moles, elements, all_species_thermo)
    #constuct the B matrix
    B = zeros(length(gasphase),length(elements))
    for i in eachindex(all_species_thermo)
        for (k,v) in all_species_thermo[i].composition
            j = get_index(String(k),elements)
            B[i,j] = v
        end
    end
    B = B'
    # Eliminate the inert element from B
    nr = Array{Int64,1}()
    for i in 1:size(B,1)
        if count(e->e>0, B[i,:]) < 2
            append!(nr,i)
        end
    end    
    ncount = 0    
    for i in nr
        B = B[1:end .!=i-ncount,:]
        ncount += 1
    end    
    # Eliminate the element constituting the inert from the elements list 
    deleteat!(elements,nr)    
    # Eliminate the empty columns (delete the inert)
    nc = Array{Int64,1}()        
    for i in 1:size(B,2)        
        if count(e->e>0,B[:,i]) < 1
            append!(nc,i)
        end
    end 
    ncount = 0
    for i in nc
        B = B[:,1:end .!= i-ncount]
        ncount += 1
    end    
    #The inert needs to be removed from the gasphase list and the moles as well
    deleteat!(gasphase,nc)
    deleteat!(moles,nc)
    
    if size(B,1) != size(B,2)
        n_null_space = length(gasphase)-rank(B)    
        A = zeros(n_null_space, length(gasphase)-n_null_space)
        AI = Matrix(I,n_null_space, n_null_space)
        B = vcat(B,hcat(A,AI))
        C = inv(B)
        nsv = C[:,end-n_null_space+1:end]        
        b = minimum(abs.(filter(e->e!=0,nsv)))
        nsv /= b
        return nsv'
    else
        B = rref(A)
        #check for rows containing only zero elements and augment with identity matrix 
        #to the right with size equal to the number of rows containing only zero as elements
        n_null_space = 0
        for i in 1:size(B,1)
            if count(e->e!=0, B(i,:)) == 0 
                B[i,i] = 1
                n_null_space += 1
            end
        end
        C = inv(B)
        nsv = C[:,end-n_null_space+1:end]        
        b = minimum(abs.(filter(e->e!=0,nsv)))
        nsv /= b
        return nsv'
    end

end

function print_reactions(stc::Matrix{Float64}, gasphase::Array{String})
    println("\nLinearly indepdendent reactions:\n")
    #Negative coefficients are assumed to be reactants     
    for i in 1:size(stc,1)
        rcts = ""
        prdts = ""
        for k in 1:length(gasphase)
            if stc[i,k] < 0 && stc[i,k] != 0
                rcts *= string(abs(stc[i,k]))*" "*gasphase[k]*" + "
            elseif stc[i,k] > 0 && stc[i,k] != 0
                prdts *= string(stc[i,k])*" "*gasphase[k]*" + "
            end
        end
        rcts = (strip(rcts))[1:end-1]
        prdts = (strip(prdts))[1:end-1]
        println(rcts * " > " * prdts)
    end
    
end

function minimize_G!(moles::Array{Float64},stc::Matrix, thermo_obj, Kp::Array{Float64} , T::Float64, p::Float64)        
    minimized = false
    nr = size(stc,1)
    Gmin = true 
    while Gmin
        mole_fracs = moles/sum(moles)                  
        G0 = IdealGas.Gmix(thermo_obj,T,p,mole_fracs)           
        for i in 1:nr            
            if !in(0.0,moles[findall(<(0),stc[i,:])])                
                equilibrate_reaction!(moles,stc[i,:],Kp[i],p)
            end
        end
        mole_fracs = moles/sum(moles)                  
        G = IdealGas.Gmix(thermo_obj,T,p,mole_fracs) 
        # println(G, "\t", G0, "\t", G-G0)     
        #check if the Gibbs energy has decreased             
        Gmin = abs(G-G0) < 1e-15 ? false : true
    end
    # mole_fracs = moles/sum(moles)                  
    # println(mole_fracs)
    # readline()    
end

    
function equilibrate_reaction!(moles::Array{Float64}, rxn_stc::Array{Float64}, Kp::Float64, p::Float64 )
    # println("Equilibrating rxn: ", rxn )
    dn_min = -1e15
    dn_max = 1e15   
    
    #The maximum possible extent of a reaction depends on the limiting concentration that is present 
    for i in eachindex(rxn_stc)        
        if rxn_stc[i] < 0
            dn_max = min(dn_max,-moles[i]/rxn_stc[i])
        end
        if rxn_stc[i] > 0
            dn_min = max(dn_min,-moles[i]/rxn_stc[i])
        end
    end    
    
    local moles_updated = copy(moles)    

    function root_func(x)                
        n1 = (moles_updated .+ rxn_stc*x)                 
        mole_fracs = n1/sum(n1)
        n_pp = mole_fracs * (p/p_std)
        prdt_p = 1.0
        prdt_r = 1.0
        for k in eachindex(rxn_stc)
            if rxn_stc[k] < 0
                prdt_r *= n_pp[k]^abs(rxn_stc[k])
            elseif rxn_stc[k] > 0                
                prdt_p *= n_pp[k]^rxn_stc[k]
            end            
        end
        return (Kp*prdt_r - prdt_p)                
    end    
    dn = find_zero(root_func,(dn_max,dn_min),Bisection())
    #moles = moles_updated .+ rxn*dn
    moles .+= rxn_stc*dn
    
end

end
