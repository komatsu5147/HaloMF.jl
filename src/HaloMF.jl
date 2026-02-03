module HaloMF
using Dierckx, ForwardDiff
export tinker08MF, tinker10MF
export bocquetMFhy, bocquetMFdm
export psMF, stMF, jenkinsMF
export stBias, stBias1, stBias2, stBias3
export tinker10Bias
export pbsBias, pbsBias1, pbsBias2
include("tinkerMF.jl")
include("bocquetMF.jl")
include("classicMF.jl")
include("stBias.jl")
include("tinkerBias.jl")
include("pbsBias.jl")
end # module
