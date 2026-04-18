# ====================================================
# Calculating the profile of background state
# ====================================================

module Background

    # """using Modules"""
    using LinearAlgebra
    using Parameters

    using ..Params: BackGroundParam, ConstantParam, Load

    # """export variables"""
    export Profiles

    # """functions"""
    function Profiles(
        BackgroundPath :: String,
        ConstantPath   :: String,
        Z              :: Vector{Float64}
    )
        # Load background parameters
        BackgroundParameterList = Load(BackGroundParam, BackgroundPath)
        ConstantParameterList   = Load(ConstantParam, ConstantPath)
        
        # Unpack parameters
        
        (; GammaT, T0, p0) = BackgroundParameterList
        (; g, Rd, Cp) = ConstantParameterList

        # Calculate temperature profile
        T :: Vector{Float64} = T0 .- GammaT .* Z

        # Calculate pressure profile
        p :: Vector{Float64} = p0 .* (T ./ T0) .^ (g / (GammaT * Rd))        

        # Calculate density profile
        rho :: Vector{Float64} = p ./ (T .* Rd)

        return T, p, rho
    end

end