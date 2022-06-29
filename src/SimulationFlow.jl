#-------------------------------------------------------------------------------

struct SpeckleSim{U<:SpeckleReadout,V<:Correlation}
    dt::DateTime
    id::UUID
    params::SpeckleParams
    bs::Beamsplitter
    elapsed::Float64
    readout::Vector{U}
    corr::Vector{V}
end

function Base.length(sim::SpeckleSim)
    return length(sim.readout)
end

#-------------------------------------------------------------------------------
# this is what enters the simulation for each "run"
struct SpeckleInstance
    source::LightSource
    detect::Detector
    bs::Beamsplitter
    γint::γIntensity
    params::SpeckleParams
end

function SpeckleInstance(params::SpeckleParams)
    source = LightSource(
                         params.n,
                         params.Em,
                         params.νm,
                         params.σ,
                         params.fγ
                        )
    detect = Detector(
                      params.deadtime,
                      params.resolution,
                      params.jitter,
                      params.efficiency,
                      params.darkcounts
                     )
    bs = Beamsplitter(1,1)
    nb = nbar(params.duration,source)*bs.t
    γint = γIntensity(nb,params.duration,detect.resolution,source)

    return SpeckleInstance(source,  detect, bs, γint, params)
end

export SpeckleInstance
#-------------------------------------------------------------------------------
