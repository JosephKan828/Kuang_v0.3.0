# ====================================================
# Time stepper module
# ====================================================

module TimeStepper

    using LinearAlgebra
    using Base.Threads

    using ..Params: ModelParam, DomainParam, var"@unpack_ModelParam", var"@unpack_DomainParam"
    using ..Radiation: Calc_Radiation
    using ..Dynamics: Calc_Dynamics


    export TimeStep, SimResult

    # Define struct for output
    mutable struct SimResult{
        T <: Number,
        A1 <: AbstractArray{Complex{T}, 3},
        A2 <: AbstractArray{Complex{T}, 3},
        A3 <: AbstractArray{Complex{T}, 3},
        A4 <: AbstractArray{Complex{T}, 3},
        A5 <: AbstractArray{Complex{T}, 3},
        A6 <: AbstractArray{Complex{T}, 3},
        A7 <: AbstractArray{Complex{T}, 3},
        A8 <: AbstractArray{Complex{T}, 3}
    }
        w1 :: A1
        w2 :: A2
        T1 :: A3
        T2 :: A4
        q  :: A5
        L  :: A6
        J1 :: A7
        J2 :: A8
    end

    # time stepper function for single wavenumber
    function _single_stepper!(
        history_slice::AbstractArray{ComplexF64, 3}, # Pass in a view to avoid copies
        InitState::Matrix{ComplexF64}, 
        p::ModelParam,
        d::DomainParam,
        kn::Float64,
        model_type::Symbol,
        rad_type::Symbol
    )
        # @unpack_ModelParam p
        
        # Matrix Exponential
        Propagator = exp((Calc_Dynamics(kn, p, model_type) .+ Calc_Radiation(p, model_type, rad_type)) .* d.dt)

        nt = size(history_slice, 3)
        history_slice[:, :, 1] .= InitState
        
        # Use a reusable buffer to avoid allocations inside the loop
        current_buffer = copy(InitState)

        @inbounds for n in 2:nt
            # Directly update the history_slice (which is a view of the main Output array)
            next_step = view(history_slice, :, :, n)
            mul!(next_step, Propagator, current_buffer)
            current_buffer .= next_step
        end
    end

    # Public interface for time stepping
    function TimeStep(
        k :: Vector{Float64},
        InitState :: Array{ComplexF64, 3},
        p :: ModelParam,
        d :: DomainParam,
        model_type :: Symbol = :full,
        rad_type :: Symbol = :all
    )

        # Calculate sizes of dimensions
        dt, T = d.dt, d.T
        nt = Int(round(T / dt)) # Number of time steps
        nstate, nens, nk = size(InitState) # Number of states, ensembles, and wavenumbers

        # Arrange output array
        output_data = Array{ComplexF64, 4}(undef, nstate, nens, nt, nk)

        # Apply simulations using views to fill the main array directly
        @inbounds Threads.@threads for i in 1:nk
            # Passing a view of the i-th wavenumber slice
            _single_stepper!(
                view(output_data, :, :, :, i), 
                InitState[:, :, i], 
                p, d, k[i], model_type, rad_type
            )
        end

        return SimResult(
            view(output_data, 1, :, :, :),
            view(output_data, 2, :, :, :),
            view(output_data, 3, :, :, :),
            view(output_data, 4, :, :, :),
            view(output_data, 5, :, :, :),
            view(output_data, 6, :, :, :),
            similar(view(output_data, 6, :, :, :)),
            similar(view(output_data, 6, :, :, :))
        )
    end

end