# ====================================================
# This script is for unifying the execution of all
# scripts in this project.
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
# Helper functions
# ====================================================

# struct for radiative heating
    mutable struct RadHeating{
        T <: Number,
        A1 <: AbstractArray{Complex{T}, 3},
        A2 <: AbstractArray{Complex{T}, 3},
        A3 <: AbstractArray{Complex{T}, 3},
        A4 <: AbstractArray{Complex{T}, 3},
        A5 <: AbstractArray{Complex{T}, 3},
        A6 <: AbstractArray{Complex{T}, 3},
        A7 <: AbstractArray{Complex{T}, 3},
        A8 <: AbstractArray{Complex{T}, 3},
        A9 <: AbstractArray{Complex{T}, 3},
        A10 <: AbstractArray{Complex{T}, 3},
        A11 <: AbstractArray{Complex{T}, 3},
        A12 <: AbstractArray{Complex{T}, 3}
    }
        qLW1 :: A1
        qLW2 :: A2
        qSW1 :: A3
        qSW2 :: A4
        tLW1 :: A5
        tLW2 :: A6
        tSW1 :: A7
        tSW2 :: A8
        wLW1 :: A9
        wLW2 :: A10
        wSW1 :: A11
        wSW2 :: A12
    end

function DerivePhy_Radiation(
    p :: ModelParam, # Model parameter struct
    state :: SimResult # State variable dictionary
)
    # Unpack parameters
    @unpack_ModelParam p

    # Calculate moisture-radiative feedback
    ## Calculate parameter for moisture-radiative feedback
    q1LW = Rq1LW*RadqScale; q2LW = Rq2LW*RadqScale
    q1SW = Rq1SW*RadqScale; q2SW = Rq2SW*RadqScale
    
    qLW1 :: Array{ComplexF64, 3} = q1LW .* state.q
    qLW2 :: Array{ComplexF64, 3} = q2LW .* state.q
    qSW1 :: Array{ComplexF64, 3} = q1SW .* state.q
    qSW2 :: Array{ComplexF64, 3} = q2SW .* state.q
    
    # Calculate temperature-radiative feedback
    ## Calculate parameter for temperature-radiative feedback
    t11LW = Rt11LW*RadtScale; t12LW = Rt12LW*RadtScale
    t11SW = Rt11SW*RadtScale; t12SW = Rt12SW*RadtScale
    t21LW = Rt21LW*RadtScale; t22LW = Rt22LW*RadtScale
    t21SW = Rt21SW*RadtScale; t22SW = Rt22SW*RadtScale
    
    tLW1 :: Array{ComplexF64, 3} = t11LW .* state.T1 .+ t21LW .* state.T2
    tLW2 :: Array{ComplexF64, 3} = t12LW .* state.T1 .+ t22LW .* state.T2
    tSW1 :: Array{ComplexF64, 3} = t11SW .* state.T1 .+ t21SW .* state.T2
    tSW2 :: Array{ComplexF64, 3} = t12SW .* state.T1 .+ t22SW .* state.T2
    
    # Calculate cloud-radiative feedback
    ## Calculate parameter for cloud-radiative feedback
    w11LW = Rw11LW*RadwScale; w12LW = Rw12LW*RadwScale
    w11SW = Rw11SW*RadwScale; w12SW = Rw12SW*RadwScale
    w21LW = Rw21LW*RadwScale; w22LW = Rw22LW*RadwScale
    w21SW = Rw21SW*RadwScale; w22SW = Rw22SW*RadwScale
    
    wLW1 :: Array{ComplexF64, 3} = w11LW .* state.w1 .+ w21LW .* state.w2
    wLW2 :: Array{ComplexF64, 3} = w12LW .* state.w1 .+ w22LW .* state.w2
    wSW1 :: Array{ComplexF64, 3} = w11SW .* state.w1 .+ w21SW .* state.w2
    wSW2 :: Array{ComplexF64, 3} = w12SW .* state.w1 .+ w22SW .* state.w2
    
    return RadHeating(
        view(qLW1, :, :, :), view(qLW2, :, :, :), view(qSW1, :, :, :), view(qSW2, :, :, :),
        view(tLW1, :, :, :), view(tLW2, :, :, :), view(tSW1, :, :, :), view(tSW2, :, :, :),
        view(wLW1, :, :, :), view(wLW2, :, :, :), view(wSW1, :, :, :), view(wSW2, :, :, :)
    )
end

# ====================================================
# Main execution
# ====================================================
function main(
    model_type :: Symbol,
    rad_type :: Symbol,
    output_path :: String
)

    # ------------------------------------------------
    # Unpack the parameters
    # ------------------------------------------------

    ## File path for TOML file
    ModelPath  = joinpath(@__DIR__, "..", "Config", "ModelParams.toml")
    DomainPath = joinpath(@__DIR__, "..", "Config", "Domain.toml")
    InitPath   = joinpath(@__DIR__, "..", "Config", "Init.toml")

    ## Load model parameters
    @unpack_ModelParam model = Load(ModelParam, ModelPath)

    ## Load domain parameters
    @unpack_DomainParam domain = Load(DomainParam, DomainPath)
    
    ## Load Initial scaling
    @unpack_InitScale init = Load(InitScale, InitPath)

    # ------------------------------------------------
    # Setup Domain
    # ------------------------------------------------

    # Setup wavenumber domain
    k :: Vector = collect(kmin:dk:kmax)
    
    # Generate initial conditions
    ## Forming scale factors for each variable
    scale_factors = collect([
        init.w1scale, init.w2scale, init.T1scale, init.T2scale, init.qscale, init.Lscale
        ]) # Shape: (nstate,)

    ## Generate random initial conditions
    InitState = Array{ComplexF64, 3}(undef, 6, domain.Ens, length(k)) # Shape: (nstate, nens, nk)

    for (i, kn) in enumerate(k)

        # Calculate structure of most unstable mode of eigenvector
        operator = Calc_Dynamics(kn, model, model_type) .+ Calc_Radiation(model, model_type, rad_type)

        # Solving for eigenvector and eigenvalue
        evecs = eigen(operator).vectors
        evals = eigen(operator).values

        # most unstable mode
        max_idx = argmax(real.(evals))
        most_unstable_mode = evecs[:, max_idx]

        for ens in 1:domain.Ens
            random_scaling = randn(ComplexF64, 6) .* scale_factors

            InitState[:, ens, i] = most_unstable_mode .* random_scaling
        end

    end

    # InitState = reshape(scale_factors, 6, 1, 1) .* randn(ComplexF64, 6, domain.Ens, length(k)) # Shape: (nstate, nens, nk)

    # ------------------------------------------------
    # Run the time stepper
    # ------------------------------------------------
    sim_result = TimeStep(k, InitState, model, domain, model_type, rad_type)
    
    # ------------------------------------------------
    # Convert L to J1 and J2
    # ------------------------------------------------
    Calc_Convection!(sim_result, model.r0, model.rq)

    # ------------------------------------------------
    # Calculate radiative heating for each time step
    # ------------------------------------------------
    rad_heating = DerivePhy_Radiation(model, sim_result)

    # ------------------------------------------------
    # Save results
    # ------------------------------------------------
    # Save metadata and initial conditions for this simulation
    h5open(output_path*"Metadata.h5", "w") do file
        # save radiative heating scaling
        write(file, "RadqScale", model.RadqScale)
        write(file, "RadtScale", model.RadtScale)
        write(file, "RadwScale", model.RadwScale)

        # save number of ensemble
        write(file, "nEnsemble", Ens)

        # save initial conditions
        write(file, "InitState", InitState)
    end
    
    # Save state evolution
    h5open(output_path*"State.h5", "w") do file
        write(file, "w1", Array(sim_result.w1))
        write(file, "w2", Array(sim_result.w2))
        write(file, "T1", Array(sim_result.T1))
        write(file, "T2", Array(sim_result.T2))
        write(file, "q", Array(sim_result.q))
        write(file, "L", Array(sim_result.L))
        write(file, "J1", Array(sim_result.J1))
        write(file, "J2", Array(sim_result.J2))
    end

    # Save radiative heating
    h5open(output_path*"Radiation.h5", "w") do file
        write(file, "qLW1", Array(rad_heating.qLW1))
        write(file, "qLW2", Array(rad_heating.qLW2))
        write(file, "qSW1", Array(rad_heating.qSW1))
        write(file, "qSW2", Array(rad_heating.qSW2))
        write(file, "tLW1", Array(rad_heating.tLW1))
        write(file, "tLW2", Array(rad_heating.tLW2))
        write(file, "tSW1", Array(rad_heating.tSW1))
        write(file, "tSW2", Array(rad_heating.tSW2))
        write(file, "wLW1", Array(rad_heating.wLW1))
        write(file, "wLW2", Array(rad_heating.wLW2))
        write(file, "wSW1", Array(rad_heating.wSW1))
        write(file, "wSW2", Array(rad_heating.wSW2))
    end
end

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