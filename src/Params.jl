""" Module for read and assigning model parameters"""

module Params

    # Using Module
    using TOML
    using Parameters: @with_kw, @unpack

    # export variables
    export ModelParam, DomainParam, BackGroundParam, InitScale, ConstantParam, Load
    export var"@unpack_ModelParam", var"@unpack_DomainParam", var"@unpack_BackGroundParam", var"@unpack_InitScale", var"@unpack_ConstantParam"

    # parameter structure
    ## Scaling on initial conditions
    @with_kw struct InitScale
        w1scale::Float64; w2scale::Float64
        T1scale::Float64; T2scale::Float64
        qscale ::Float64; Lscale ::Float64
    end

    ## domain setting information    
    @with_kw struct DomainParam
        T   ::Float64; X   ::Float64; Z  ::Float64; Ens :: Int
        dt  ::Float64; dx  ::Float64; dz ::Float64
        kmin::Float64; kmax::Float64; dk ::Float64
    end

    ## information in background state
    @with_kw struct BackGroundParam
        GammaT::Float64; T0::Float64; p0::Float64
    end

    ## Constant struct
    @with_kw struct ConstantParam
        g::Float64; Rd::Float64; Cp::Float64
    end

    ## physical parameters
    @with_kw struct ModelParam
        c1::Float64; c2::Float64; epsilon::Float64
        a1::Float64; a2::Float64
        d1::Float64; d2::Float64; m1::Float64; m2::Float64
        F::Float64; b1::Float64; b2::Float64; f::Float64; tauL::Float64
        g0::Float64; gq::Float64; r0::Float64; rq::Float64
        RadqScale::Float64; RadtScale::Float64; RadwScale::Float64
        Rq1LW::Float64; Rq2LW::Float64; Rq1SW::Float64; Rq2SW::Float64
        Rt11LW::Float64; Rt12LW::Float64; Rt21LW::Float64; Rt22LW::Float64
        Rt11SW::Float64; Rt12SW::Float64; Rt21SW::Float64; Rt22SW::Float64
        Rw11LW::Float64; Rw12LW::Float64; Rw21LW::Float64; Rw22LW::Float64
        Rw11SW::Float64; Rw12SW::Float64; Rw21SW::Float64; Rw22SW::Float64
    end

    # """Data type for union different struct types"""
    const StructParam = Union{InitScale, DomainParam, BackGroundParam, ModelParam, ConstantParam}

    function _flatten_dict!(dest::Dict{Symbol, Any}, src::Dict{String, Any})::Dict{Symbol, Any}
        for (key, value) in src
            if value isa Dict{String, Any}
                _flatten_dict!(dest, value)
            else
                dest[Symbol(key)] = value
            end
        end
        return dest
    end

    # """Load TOML file and save to different struct types"""
    function Load(T::Type{<:StructParam}, FilePath::String)::StructParam
        # parse TOML file
        data = TOML.parsefile(FilePath)

        # Flatten dictionary into a single layer
        flat_data = _flatten_dict!(Dict{Symbol, Any}(), data)

        # Find fields allowed in the struct
        allowed_fields = fieldnames(T)
        
        # dictionary for allowed field only
        clean_dict = Dict{Symbol, Any}(k => v for (k, v) in flat_data if k in allowed_fields)

        # return into struct
        return T(; clean_dict...)
    end

end
