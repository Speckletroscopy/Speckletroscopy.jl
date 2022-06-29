"""
	function σTemp(ω0::Real,temp::Number)

Returns the standard deviation due to Doppler broadening of frequency ω0 at temp
"""
function σTemp(ω0::Real,temp::Number)
	kbOverMhC2 = 9.178e-14;
	return sqrt(kbOverMhC2*temp)*ω0
end

export σTemp
#-------------------------------------------------------------------------------

struct LightSource
    n::Integer # number of atoms
    Em::Vector # line magnitudes
    νm::Vector # central frequencies of lines in GHz
    σ::Number # Doppler broadening in GHz
    νMin::Number # minimum bandpass frequency in GHz
    νMax::Number # maximum bandpass frequency in GHz
    γRate::Number # photon count rate in GHz
    function LightSource(n::Integer,Em::Vector,νm::Vector,σ::Number,νMin::Number,νMax::Number,γRate::Number)
        @assert length(Em) == length(νm) "Vectors for line magnitudes and frequencies must have the same length"
        @assert νMin < νMax "Bandpass minimum must be less than maximum"
        new(
            n,
            convert(Vector{Complex},Em),
            convert(Vector{Real},νm),
            convert(Real,σ),
            convert(Real,νMin),
            convert(Real,νMax),
            convert(Real,γRate)
        )
    end
end

function LightSource(n::Integer,Em::Vector,νm::Vector,σ::Number, γRate::Number)
    νMin = min(νm...) - 5*σ # automatically set the bandpass minimum to 5σ below the lowest line frequency
    νMax = max(νm...) + 5*σ # automatically set the bandpass maximum to 5σ above the highest line frequency
    return LightSource(n,Em,νm,σ,νMin,νMax,γRate)
end

function LightSource(params::Dict)
    return LightSource(
        params[:n],
        params[:Em],
        params[:νm],
        params[:σ],
        params[:fγ]
    )
end

export LightSource

#-------------------------------------------------------------------------------

"""
    ν0(source::LightSource)

Returns the average frequency of spectral lines in source
"""
function ν0(source::LightSource)
    return mean(source.νm)
end

export ν0

#-------------------------------------------------------------------------------

"""
    Δm(source::LightSource)

Returns vector of the sources line separation from the average frequency of all lines
"""
function Δm(source::LightSource)
    return source.νm .- ν0(source) 
end

export Δm

#-------------------------------------------------------------------------------

"""
    lineShifts(source::LightSource)

Returns a vector with the frequency differences between all spectral lines
"""
function lineShifts(source::LightSource)
    return lineShifts(source.νm)
end

"""
    lineShifts(source::LightSource)

Returns a vector with the frequency differences between all spectral lines
"""
function lineShifts(lines::Vector{T}) where {T<:Number}
    @assert length(lines)>1 "Must have at least two spectral lines"
    return map(x->abs(x[2]-x[1]),subsets(lines,2))
end

export lineShifts

#-------------------------------------------------------------------------------

"""
    nbar(t::Real,source::LightSource)

Returns the average photon counts received over a duration t from source
"""
function nbar(t::Real,source::LightSource)
    return t*source.γRate
end

export nbar

################################################################################
# Electric field 
################################################################################
"""
    eField(νn::Vector,ϕmn::Matrix,source::LightSource)

Container holding frequencies and phases for one realization of the electric field.
"""
struct eField
    νn::Vector
    ϕmn::Matrix
    source::LightSource

    function eField(νn::Vector,ϕmn::Matrix,source::LightSource)
        @assert (size(ϕmn)[1] == length(source.νm) && size(ϕmn)[2] == length(νn)) "Phase array must have shape (m,n)"
        new(νn,ϕmn,source)
    end
end

"""
    eField(source::LightSource, seed::Integer = -1)

Generate frequencies and phases for a single instance of the electric field.
"""
function eField(source::LightSource, seed::Integer = -1)
    if seed != -1
       Random.seed!(seed) 
    end
    νDist = Normal(ν0(source),source.σ)
    νn    = rand(νDist,source.n)
    ϕmn   = 2*π*rand(Float64,(length(source.νm),source.n))
    return eField(νn,ϕmn,source)
end

export eField

#-------------------------------------------------------------------------------

"""
	function eField(t::Number,field::eField)

Returns the electric field value at time t(ns)
"""
function eFieldT(t::Number,field::eField)
	# generate frequencies
	ωmn = transpose(field.νn) .+ Δm(field.source)
    ωmn .*= 2*π
	# add the phase
	exponentmn = -im*(t*ωmn+field.ϕmn)
	# put them in the exponent
	enm = exp.(exponentmn)
	# multiply by the field magnitude
	fieldnm = field.source.Em .* enm
	# add it all together
	return sum(ivec(fieldnm))
end

export eFieldT

################################################################################
# Calculate intensity
################################################################################
"""
    intensity(et::Number)

Returns the intensity given the value of the electric field at a particular point in time and space
"""
function intensity(et::Number)
    return real(et*conj(et))
end

"""
    intensity(t::Real,field::eField)

Returns the intensity given an instance of the EM field and a time t(ns).
"""
function intensity(t::Number,field::eField)
    efield = eFieldT(t,field)
    return intensity(efield)
end

export intensity

#-------------------------------------------------------------------------------

"""
    γIntensity(nbar::T,γvec::Vector{T}) where T<:Real

Returns γIntensity object which contains the average photon count rate and the renormalized intensity time series
"""
struct γIntensity{T<:Real}
    nbar::T
    γvec::Vector{T}
    function γIntensity(nbar::T,γvec::Vector{T}) where T<:Real
        total = sum(γvec)
        coeff = nbar/total
        return new{T}(nbar,coeff*γvec)
    end
end

function γIntensity(nbar::Real,t::Real,dt::Real,source::LightSource)
    field = eField(source)
    return γIntensity(nbar,t,dt,field)
end

function γIntensity(nbar::Real,t::Real,dt::Real,field::eField)
    intensitySeries = map(time->intensity(time,field),0:dt:t)
    return γIntensity(nbar,intensitySeries)
end

export γIntensity

#-------------------------------------------------------------------------------

# array indexing and iteration interface for γIntensity object
function Base.length(γint::γIntensity)
    return length(γint.γvec)
end

#-------------------------------------------------------------------------------

function Base.getindex(γint::γIntensity,i::Int)
    return γint.γvec[i]
end

#-------------------------------------------------------------------------------

function Base.firstindex(γint::γIntensity)
    return γint.γvec[1]
end

#-------------------------------------------------------------------------------

function Base.lastindex(γint::γIntensity)
    return γint.γvec[end]
end

#-------------------------------------------------------------------------------

function Base.iterate(γint::γIntensity)
    return (Base.firstindex(γint),1)
end

function Base.iterate(γint::γIntensity,state::Int)
    return state >= length(γint) ? nothing : (γint[state+1],state+1)
end

#-------------------------------------------------------------------------------

"""
    Beamsplitter(r::Number,t::Number)

Struct to hold beam splitter coefficients.
"""
struct Beamsplitter
    r::Number
    t::Number
    function Beamsplitter(r::Number,t::Number)
        beamNorm = sqrt(r^2+t^2)
        new(r/beamNorm,t/beamNorm)
    end
end

export Beamsplitter

