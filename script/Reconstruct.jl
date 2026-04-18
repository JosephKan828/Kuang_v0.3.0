# ====================================================
# This file is to reconstruct state vector back to 
# Fourier and Galerkin space
# ====================================================

# ====================================================
# Using package
# ====================================================
# Official package
using HDF5
using LinearAlgebra

# Individual module
include(joinpath(@__DIR__, "..", "src", "API.jl"))
using .API

# ====================================================
# Helper function
# ====================================================

"""Function for loading state"""
function _load_state(
        filepath :: String
    )

    data = Dict{String, Array}()

    h5open(filepath, "r") do f
        for key in keys(f)
            data[key] = read(f, key)
        end
    end

    return data
end

"""Function for get density profile"""
function _get_density(
    T0 :: Float64,
    Γ  :: Float64,
    p0 :: Float64,
    z  :: Vector{Float64},
    R  :: Float64,
    g  :: Float64
)
    # Get temperature profile
    T :: Vector{Float64} = T0 .- Γ .* z

    # Get pressure profile by hypsometric relation
    p :: Vector{Float64} = p0 .* (1 .- Γ.*z/T0) .^ (g/(R*Γ))

    # Get density profile by ideal gas law
    ρ :: Vector{Float64} = p ./ (R .* T)

    return ρ

end

"""Function for generate Fourier basis"""
function _get_Fourier(
        kn :: Float64,
        x  :: Vector{Float64}
    )

    # scale physical distance
    x_scale :: Vector{Float64} = x ./ 4.32e6

    return exp.((im * kn) .* x_scale)
end

"""Function for reconstructing state vector to Fourier space"""
function _reconstruct_Fourier(
    state_data :: Matrix{ComplexF64},
    FourierBasis :: Vector{ComplexF64}
)
    # get shape of inpute 
    nens, nt = size(state_data)
    nx = length(FourierBasis)

    # Reshape data and basis
    state_reshaped = reshape(state_data, nens, nt, 1) # shape: (nens, ntime, 1)
    state_flattened = reshape(state_reshaped, nens * nt, 1) # shape: (nens*ntime, 1)

    basis_reshaped = reshape(FourierBasis, 1, nx) # shape: (1, nx)

    # reconstruction
    FourierState = Array{ComplexF64}(undef, nens*nt, nx)

    mul!(FourierState, state_flattened, basis_reshaped) # shape: (nens*ntime, nx)
    
    FourierState_reshaped = reshape(FourierState, nens, nt, nx) # shape: (nens, ntime, nx)

    return FourierState_reshaped
end

"""Function for generate Galerkin basis"""
function _get_Galerkin(
    z  :: Vector{Float64},
    Γ  :: Float64,
    Γd :: Float64
)
    # Calculate structure for Galerkin basis 1 and 2
    G1 :: Vector{Float64} = π/2 .* sin.(π .* z ./ maximum(z))
    G2 :: Vector{Float64} = π/2 .* sin.(2*π .* z ./ maximum(z))

    # Calculate structure with stratification for Galerkin basis 1 and 2
    G1_strat :: Vector{Float64} = G1 .* (Γd - Γ)
    G2_strat :: Vector{Float64} = G2 .* (Γd - Γ)

    return G1, G2, G1_strat, G2_strat
end

"""Function for projecting Fourier state onto Galerkin basis"""
function _reconstruct_Galerkin(
    state::Array{Float64, 3},
    basis::Vector{Float64}
)
    # get shape of inpute 
    nens, nt, nx = size(state)
    nz = length(basis)

    # Reshape data and basis
    state_reshaped = reshape(state, nens, nt, nx, 1) # shape: (nens*ntime, nx)
    state_flattened = reshape(state_reshaped, nens * nt * nx, 1) # shape: (nens*ntime*nx, 1)

    basis_reshaped = reshape(basis, 1, nz) # shape: (1, nz)

    # reconstruction
    GalerkinState = Array{Float64}(undef, nens*nt*nx, nz)

    mul!(GalerkinState, state_flattened, basis_reshaped) # shape: (nens*ntime*nx, nz)
    GalerkinState_reshaped = reshape(GalerkinState, nens, nt, nx, nz) # shape: (nens, ntime, nx, nz)

    return GalerkinState_reshaped
    
end

# ====================================================
# Main function
# ====================================================
function main(
    input_path :: String
    )

    # ------------------------------------------------
    # Load data
    # ------------------------------------------------
    # unpack domain data to get information for basis
    DomainPath = joinpath(@__DIR__, "..", "Config", "Domain.toml")
    ConstantPath = joinpath(@__DIR__, "..", "Config", "Constant.toml")
    BackGroundPath = joinpath(@__DIR__, "..", "Config", "Background.toml")

    @unpack_DomainParam domain = Load(DomainParam, DomainPath)
    @unpack_ConstantParam constant = Load(ConstantParam, ConstantPath)
    @unpack_BackGroundParam background = Load(BackGroundParam, BackGroundPath)

    # Load state evolution
    # state_data = _load_state(input_path*"State.h5")
    # state_rad  = _load_state(input_path*"Radiation.h5")

    # create x and k series
    x :: Vector{Float64} = collect(-domain.X:domain.dx:domain.X)
    z :: Vector{Float64} = collect(0:domain.dz:domain.Z)
    k :: Vector{Float64} = collect(domain.kmin:domain.dk:domain.kmax)
    
    nx :: Int = length(x)
    nz :: Int = length(z)
    nk :: Int = length(k)

    # Get density profile for later use
    ρ = _get_density(
        background.T0,
        background.GammaT,
        background.p0,
        z,
        constant.Rd,
        constant.g
    )

    ρ_inv = 1.0 ./ reshape(ρ, 1, 1, 1, length(ρ)) # reshape to (1, 1, nz) for later broadcasting

    # Defining variable mapping
    state_vars = ["w1", "w2", "T1", "T2", "J1", "J2"]
    rad_vars   = ["qLW1", "qLW2", "qSW1", "qSW2", "tLW1", "tLW2", "tSW1", "tSW2", "wLW1", "wLW2", "wSW1", "wSW2"]

    # ------------------------------------------------
    # Generate basis
    # ------------------------------------------------
    # create Fourier basis for different wavenumber
    FourierBasis = Matrix{ComplexF64}(undef, nx, nk)

    # generate Fourier basis
    for (i, kn) in enumerate(k)
        FourierBasis[:, i] .= _get_Fourier(kn, x)
    end

    # generate Galerkin basis
    G1, G2, G1_strat, G2_strat = _get_Galerkin(z, background.GammaT, constant.g / constant.Cp)

    # ------------------------------------------------
    # Iteratively save state
    # ------------------------------------------------

    FourierPath = joinpath(input_path, "FourierState.h5")
    GalerkinPath = joinpath(input_path, "GalerkinState.h5")

    # Open file
    h5open(FourierPath, "w") do h5f
        h5open(GalerkinPath, "w") do h5g
            # Write dimension
            write(h5f, "x", x)
            write(h5g, "x", x)
            write(h5g, "z", z)

            # Process ALL variables one by one
            all_vars = ["w1", "w2", "T1", "T2", "J1", "J2", "qLW1", "qLW2", "qSW1", "qSW2", "tLW1", "tLW2", "tSW1", "tSW2", "wLW1", "wLW2", "wSW1", "wSW2"]
            
            @inbounds for var in all_vars
                # Select source file
                src_file = var in state_vars ? "State.h5" : "Radiation.h5"

                # Load data
                var_data = h5read(joinpath(input_path, src_file), var)

                nens, nt, _ = size(var_data)

                # pre-allocate: create empty dataset
                dset_f = create_dataset(h5f, var, datatype(ComplexF64), dataspace(nens, nt, nx, nk))
                dset_g = create_dataset(h5g, var, datatype(ComplexF64), dataspace(nens, nt, nx, nz, nk))

                # Select appropriate basis
                basis = if  var == "w1" G1
                elseif var == "w2" G2
                elseif var in ["T1", "J1", "qLW1", "qSW1", "tLW1", "tSW1", "wLW1", "wSW1"] G1_strat
                else G2_strat
                end

                # Compute and write
                Threads.@threads for j in 1:nk
                    ## Fourier reconstruction
                    f_slice = _reconstruct_Fourier(var_data[:, :, j], FourierBasis[:, j])

                    dset_f[:, :, :, j] = f_slice

                    ## Galerkin reconstruction
                    g_slice = _reconstruct_Galerkin(real.(f_slice), basis)

                    dset_g[:, :, :, :, j] = g_slice .* ρ_inv
                end

                var_data = nothing
                GC.gc()
            end

        end
    end

end

# ====================================================
# Execute main function
# ====================================================


if abspath(PROGRAM_FILE) == @__FILE__

    # Parse command-line arguments
    input_path_arg = ARGS[1]

    # Convert strings to symbols
    input_path = String(input_path_arg)

    main(input_path)
end
