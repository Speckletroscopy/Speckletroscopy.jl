################################################################################
# Counts are generated at the detector, so photon count generating functions
# should live here.
################################################################################
"""
	function beCount(nbar::Real)

Returns Bose-Einstein distributed counts for average count rate nbar
"""
function beCount(nbar::Real)
	p = 1/(nbar+1)
	f = p*nbar
	s = p
	r = rand()
	count = 0
	while r>s
		count += 1
		p *= f
		s += p
	end
	return count
end

################################################################################
# Container for detector properties
################################################################################
struct Detector
    deadtime::Float64
    resolution::Float64
    jitter::Float64
    efficiency::Float64
    darkcounts::Float64
    function Detector(deadtime::Float64,resolution::Float64,jitter::Float64,efficiency::Float64,darkcounts::Float64)
        @assert 0 ≤ efficiency ≤ 1 "Efficiency must be between 0 and 1"
        new(
            deadtime,
            resolution,
            jitter,
            efficiency,
            darkcounts
        )
    end
end

"""
    Detector(deadtime::Number,resolution::Number,jitter::Number,efficiency::Number,darkcounts::Number)

Returns a Detector object. Also accepts a dictionary with keywords matching the inputs below.

Inputs:
    deadtime   : dead time in nanoseconds
    resolution : time resolution in nanoseconds
    jitter     : timing jitter in nanoseconds
    efficiency : quantum efficiency ∈ [0,1]
    darkcounts : dark count rate in GHz
"""
function Detector(deadtime::Number,resolution::Number,jitter::Number,efficiency::Number,darkcounts::Number)
    return Detector(
        convert(Float64,deadtime),
        convert(Float64,resolution),
        convert(Float64,jitter),
        convert(Float64,efficiency),
        convert(Float64,darkcounts)
    )
end

function Detector(params::Dict)
    return Detector(
            params[:deadtime],
            params[:resolution],
            params[:jitter],
            params[:efficiency],
            params[:darkcounts]
        )
end

export Detector

################################################################################
# Data structures for storing the readout of the detector
################################################################################

# define a type for both dense and sparse sim results
abstract type SpeckleReadout end
export SpeckleReadout

#-------------------------------------------------------------------------------

"""
    DenseReadout(beamData::IndexedTable, corrData::IndexedTable)

Stores beam intensity from Speckle simulation in a dense format.
"""
struct DenseReadout <: SpeckleReadout
    beam1::Vector
    beam2::Vector
    function DenseReadout(beam1::Vector,beam2::Vector)
        @assert length(beam1) == length(beam2) "Beam readout vectors must have same length"
        new(beam1,beam2)
    end
end

export DenseReadout
#-------------------------------------------------------------------------------
"""
    SparseReadout(beamData::NDSparse, corrData::NDSparse)

Stores beam intensity from Speckle simulation in a sparse format.
"""
struct SparseReadout <: SpeckleReadout
    beam1::SparseVector
    beam2::SparseVector
    function SparseReadout(beam1::SparseVector,beam2::SparseVector)
        @assert length(beam1) == length(beam2) "Beam readout vectors must have same length"
        new(beam1,beam2)
    end
end

export SparseReadout

################################################################################
# Functions to generate counts from the detector
################################################################################
"""
    denseReadout(t::Number,source::LightSource,bs,::Beamsplitter,detect::Detector)

Returns a dense vector containing the counts received by the detector for the given duration t and LightSource
"""
function denseReadout(duration::Number,source::LightSource,bs::Beamsplitter,detect::Detector)
    field = eField(source)
    return denseReadout(duration,field,bs,detect)
end

"""
    denseReadout(t::Number,field::eField,bs::Beamsplitter,detect::Detector)

Returns a dense vector containing the counts received by the detector for the given duration t and EM field
"""
function denseReadout(t::Number,field::eField,bs::Beamsplitter,detect::Detector)
    nb = nbar(t,field.source)
    γint = γIntensity(nb,t,detect.resolution,field)
    return denseReadout(γint,bs,detect)
end

"""
    denseReadout(γint:γIntensity,bs::Beamsplitter,detect::Detector)

Returns a dense vector containing the counts received by the detector for the given photon intensity time-series
"""
function denseReadout(γint::γIntensity,bs::Beamsplitter,detect::Detector)
    beam1Dist = Poisson.(bs.t*γint.γvec)
    beam2Dist = Poisson.(bs.r*γint.γvec)
    beam1 = zeros(Int,length(γint))
    beam2 = zeros(Int,length(γint))
    ideadtime = convert(Int,ceil(detect.deadtime/detect.resolution))

    # counts for beam 1
    i = 1
    while i ≤ length(beam1)
        ct = rand(beam1Dist[i])
        if ct != 0
            beam1[i] = ct
            i+=ideadtime
        end
        i+=1
    end

    # counts for beam 2
    i = 1
    while i ≤ length(beam2)
        ct = rand(beam2Dist[i])
        if ct != 0
            beam2[i] = ct
            i+=ideadtime
        end
        i+=1
    end

    return DenseReadout(beam1,beam2)
end

export denseReadout
