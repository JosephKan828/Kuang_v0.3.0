# ====================================================
# This script is for analyzing the growth rate and 
# phase speed of this system 
# by performing eigenvalue analysis.
# ====================================================

# ====================================================
# using package
# ====================================================
# Official packages
using HDF5
using LinearAlgebra

# Individual modules
include(joinpath(@__DIR__, "..", "src", "API.jl"))
using .API

# ====================================================
# Helper function
# ====================================================

function Operator_Matrix(
    kn :: Float64,
    p :: ModelParam,
    model_type :: Symbol,
    rad_type :: Symbol
)
    return Calc_Dynamics(kn, p, model_type) .+ Calc_Radiation(p, model_type, rad_type)
end

# ====================================================
# Main function
# ====================================================

function main(
    model_type :: Symbol,
    rad_type :: Symbol,
    output_path :: String
)

    # ------------------------------------------------
    # Unpack parameters
    # ------------------------------------------------

    ## File path for TOML file
    ModelPath  = joinpath(@__DIR__, "..", "Config", "ModelParams.toml")
    DomainPath = joinpath(@__DIR__, "..", "Config", "Domain.toml")
    
    ## Load model parameters
    @unpack_ModelParam model = Load(ModelParam, ModelPath)

    ## Load domain parameters
    @unpack_DomainParam domain = Load(DomainParam, DomainPath)
    
    # ------------------------------------------------
    # Calculate operator of model at each wavenumber
    # ------------------------------------------------

    # Setup wavenumber
    k :: Vector = collect(kmin:dk:kmax)

    # Array for saving operators
    Operators = Array{Matrix{Float64}, 1}(undef, length(k))

    for (i, kn) in enumerate(k)
        Operators[i] = Operator_Matrix(kn, model, model_type, rad_type)
    end

    # ------------------------------------------------
    # Calculate eigenvalues
    # ------------------------------------------------

    # Array for saving eigenvalues and eigenvectors
    GrowthRates = Array{Float64, 2}(undef, 6, length(k))
    PhaseSpeeds = Array{Float64, 2}(undef, 6, length(k))

    for (i, kn) in enumerate(k)
        op = Operators[i]
        evals = eigen(op).values

        # sort eigenvalues by their real part in descending order
        evals_sorted = sort(evals, by=real, rev=true)

        # Save growth rate and phase speed
        GrowthRates[:, i] .= real(evals_sorted)
        PhaseSpeeds[:, i] .= -(imag(evals_sorted) ./ kn) .* (4320000/86400)
    end

    # ------------------------------------------------
    # Save results
    # ------------------------------------------------

    h5open(output_path*"EigenAnalysis.h5", "w") do file
        write(file, "k", k)
        write(file, "GrowthRates", GrowthRates)
        write(file, "PhaseSpeeds", PhaseSpeeds)
    end

end

# ====================================================
# execute main function
# ====================================================
if abspath(PROGRAM_FILE) == @__FILE__

    # Parse command-line arguments for model_type and rad_type
    model_type_arg = ARGS[1]
    rad_type_arg = ARGS[2]
    output_path_arg = ARGS[3]

    # Convert strings to symbols
    model_type = Symbol(model_type_arg)
    rad_type = Symbol(rad_type_arg)
    output_path = String(output_path_arg)

    main(model_type, rad_type, output_path)
end