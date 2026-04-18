module Radiation

    # Using package
    ## Using official one
    using LinearAlgebra

    ## Using individual one
    include(joinpath(@__DIR__, "Params.jl"))
    using ..Params: ModelParam, var"@unpack_ModelParam"

    # Export public interface
    export Calc_Radiation

    # Calculate LW radiative heating in full model
    function _full_LW(p::ModelParam)
        
        # unpack parameters
        @unpack_ModelParam p

        # Calculate radiative heating
        Rw11 = Rw11LW*RadwScale
        Rw21 = Rw21LW*RadwScale
        Rw12 = Rw12LW*RadwScale
        Rw22 = Rw22LW*RadwScale
        Rt11 = Rt11LW*RadtScale
        Rt21 = Rt21LW*RadtScale
        Rt12 = Rt12LW*RadtScale
        Rt22 = Rt22LW*RadtScale
        Rq1  = Rq1LW*RadqScale
        Rq2  = Rq2LW*RadqScale

        # Calculate element at the matrix
        m51 = -(d1*Rw11 + d2*Rw12)*0.0
        m52 = -(d1*Rw21 + d2*Rw22)*0.0
        m53 = -(d1*Rt11 + d2*Rt12)*0.0
        m54 = -(d1*Rt21 + d2*Rt22)*0.0
        m55 = -(d1*Rq1  + d2*Rq2)*0.0

        return Float64[
            0.0   0.0   0.0   0.0   0.0   0.0;
            0.0   0.0   0.0   0.0   0.0   0.0;
            Rw11  Rw21  Rt11  Rt21  Rq1   0.0;
            Rw12  Rw22  Rt12  Rt22  Rq2   0.0;
            m51   m52   m53   m54   m55   0.0;
            0.0   0.0   0.0   0.0   0.0   0.0
        ]
    end

    # Calculate radiative heating in one-way model
    function _oneway_LW(p::ModelParam)
        println("TBD")
        return zeros(Float64, 6, 6)
    end

    # Calculate SW radiative heating in full model
    function _full_SW(p::ModelParam)
        
        # unpack parameters
        @unpack_ModelParam p

        # Calculate radiative heating
        Rw11 = Rw11SW*RadwScale
        Rw21 = Rw21SW*RadwScale
        Rw12 = Rw12SW*RadwScale
        Rw22 = Rw22SW*RadwScale
        Rt11 = Rt11SW*RadtScale
        Rt21 = Rt21SW*RadtScale
        Rt12 = Rt12SW*RadtScale
        Rt22 = Rt22SW*RadtScale
        Rq1  = Rq1SW*RadqScale
        Rq2  = Rq2SW*RadqScale

        # Calculate element at the matrix
        m51 = -(d1*Rw11 + d2*Rw12)
        m52 = -(d1*Rw21 + d2*Rw22)
        m53 = -(d1*Rt11 + d2*Rt12)
        m54 = -(d1*Rt21 + d2*Rt22)
        m55 = -(d1*Rq1  + d2*Rq2)

        return Float64[
            0.0   0.0   0.0   0.0   0.0   0.0;
            0.0   0.0   0.0   0.0   0.0   0.0;
            Rw11  Rw21  Rt11  Rt21  Rq1   0.0;
            Rw12  Rw22  Rt12  Rt22  Rq2   0.0;
            m51   m52   m53   m54   m55   0.0;
            0.0   0.0   0.0   0.0   0.0   0.0
        ]

    end

    # Calculate radiative heating in one-way model
    function _oneway_SW(p::ModelParam)
        println("TBD")
        return zeros(Float64, 6, 6)
    end

    # Public interface for radiation
    function Calc_Radiation(
        p          :: ModelParam, # Model parameter struct
        model_type :: Symbol,     # :oneway or :full
        rad_type   :: Symbol      # :LW or :SW or :all
    )

        # Pre-allocate matrix for output radiative heating
        RadMatrix :: Matrix{ComplexF64} = zeros(ComplexF64, 6, 6)

        # LW options
        if rad_type === :all || rad_type === :LW
            if model_type === :full
                RadMatrix .+= _full_LW(p)

            elseif model_type === :oneway
                RadMatrix .+= _oneway_LW(p)

            else
                error("Invalid model type. Use :full or :oneway.")
            end
        end

        if rad_type === :all || rad_type === :SW
            if model_type === :full
                RadMatrix .+= _full_SW(p)

            elseif model_type === :oneway
                RadMatrix .+= _oneway_SW(p)

            else
                error("Invalid model type. Use :full or :oneway.")
            end
        end

        return RadMatrix
    end

end