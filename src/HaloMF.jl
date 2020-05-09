module HaloMF
using Dierckx
export tinker08MF, tinker10MF
export bocquetMFhy, bocquetMFdm
export psMF, stMF, jenkinsMF
include("tinkerMF.jl")
include("bocquetMF.jl")
include("classicMF.jl")
end # module
