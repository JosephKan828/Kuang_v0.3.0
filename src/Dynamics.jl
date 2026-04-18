# ====================================================
# This module is for linear model kernel
# ====================================================

module Dynamics

    # Using modules
    ## Official packages
    using LinearAlgebra
    using Parameters

    ## Individual
    include(joinpath(@__DIR__, "Params.jl"))
    using ..Params: ModelParam, var"@unpack_ModelParam"

    # export variables
    export Calc_Dynamics

    # Full dynamics
    function _full_dynamics(kn::Float64, p::ModelParam)
        # This macro is now exported from Params and correctly brought into scope
        @unpack_ModelParam p

        A::Float64 = 1 - 2*f + (b2-b1)/F
        B::Float64 = 1+(b2+b1)/F - A*r0
        
        return Float64[
            -epsilon    0.0         (kn*c1)^2    0.0          0.0                  0.0;
             0.0       -epsilon      0.0         (kn*c2)^2    0.0                  0.0;
            -1.0        0.0         -1.5*rq       0.0          rq                  (1+r0);
             0.0       -1.0          1.5*rq       0.0         -rq                  (1-r0);
             a1         a2           1.5*rq*(d1-d2) 0.0        rq*(d2-d1)          -(d1*(1+r0)+d2*(1-r0));
             f/(B*tauL) (1-f)/(B*tauL) -1.5*A*rq/(B*tauL) 0.0  A*rq/(B*tauL)       -1/tauL
        ]
    end

    # TODO: 
    # derive one-way dynamical core
    function _oneway_dynamics()
        println("TBD")
    end

    # export switcher for the two modes of dcore
    function Calc_Dynamics(
        kn  ::Float64,
        p   ::ModelParam,
        mode::Symbol # :full or :oneway
        )
        if mode === :full
            return _full_dynamics(kn, p)
        elseif mode === :oneway
            return _oneway_dynamics()
        else
            error("Invalid mode. Choose either :full or :oneway.")
        end
    end


end