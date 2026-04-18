# ====================================================
# Convective parameterization
# ====================================================

module Convection

    using LinearAlgebra

    using ..TimeStepper

    export Calc_Convection!

    function Calc_Convection!(
        state:: SimResult,
        r0:: Float64, rq :: Float64,
        )
        
        U = @. r0 * state.L + rq * (state.q - 1.5 * state.T1)

        state.J1 .= state.L .+ U
        state.J2 .= state.L .- U

        return nothing
    end

end