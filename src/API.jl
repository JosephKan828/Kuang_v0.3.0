module API

include("Params.jl")
include("Dynamics.jl")
include("Radiation.jl")
include("Background.jl")
include("TimeStepper.jl")
include("Convection.jl")

using .Params
using .Dynamics
using .Radiation
using .Convection
using .Background
using .TimeStepper

# API of Params.jl
export ModelParam, DomainParam, BackGroundParam, InitScale, ConstantParam, Load
export var"@unpack_ModelParam", var"@unpack_DomainParam", var"@unpack_BackGroundParam", var"@unpack_InitScale", var"@unpack_ConstantParam"

# API of Dynamics.jl, Radiation.jl, Convection.jl
export Calc_Dynamics, Calc_Radiation, Calc_Convection!
export DerivePhy_Radiation, RadHeating

# API of Background.jl
export Profiles

# API of TimeStepper.jl
export TimeStep, SimResult

end