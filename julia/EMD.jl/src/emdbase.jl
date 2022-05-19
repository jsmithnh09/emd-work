module Emdbase

using ..Util, ..Config

function emd(x::AbstractVector; kwargs...)
    
    # create the configuration object.
    cfg = EMDConfig(x, kwargs...)


end # emd


end # module