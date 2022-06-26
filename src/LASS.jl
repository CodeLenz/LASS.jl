module LASS

using Statistics

include("bins.jl")
include("flass.jl")
include("dlass.jl")

export NBin_data
export Generate_bins
export Lass, dLass

end # module LASS
