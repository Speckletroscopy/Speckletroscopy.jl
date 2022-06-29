################################################################################
# Correlation Data Structures
################################################################################

abstract type Correlation end
abstract type CorrelationVector <: Correlation end
abstract type CorrelationMatrix <: Correlation end

struct DenseCorrelationVector <: CorrelationVector
    n::Int # column of correlation matrix
    data::Vector
end

struct SparseCorrelationVector <: CorrelationVector
    n::Int # column of correlation matrix
    data::SparseVector
end

struct DenseCorrelationMatrix <: CorrelationMatrix
    data::Matrix
end

struct SparseCorrelationMatrix <: CorrelationMatrix
    data::SparseMatrixCSC
end

function Base.:+(x::T...) where {T<:Correlation}
    if length(x) == 1
        return x[1]
    end
    newdata = x[1].data
    for corr in x[2:end]
        newdata += corr.data
    end
    if T<:CorrelationVector
        icheck = true
        for corr in x[2:end]
            icheck *= x[1].n == corr.n
        end
        newn = icheck ? x[1].n : -1
        return typeof(x[1])(newn,newdata)
    else
        return typeof(x[1])(newdata)
    end
end
################################################################################
# Correlation Functions
################################################################################

"""
    function correlate2d(readout::DenseReadout, window::Int)

Calculates the entire correlation matrix between beam 1 and 2
"""
function correlate2d(readout::DenseReadout, window::Int)
    @assert window < length(readout.beam1) "Window must be smaller than length of beam vector"

    mnMax = length(readout.beam1) - window
    out = map(mn->correlate(readout,mn[1],mn[2],window), Iterators.product(1:mnMax,1:mnMax))
    
    return DenseCorrelationMatrix(out)
end
export correlate2d

#-------------------------------------------------------------------------------

"""
    function correlate1d(readout::DenseReadout, n::Int, window::Int)

Calculates the nth column of the correlation matrix between beam 1 and 2.
"""
function correlate1d(readout::DenseReadout, n::Int, window::Int)
    @assert window < length(readout.beam1) "Window must be smaller than length of beam vector"

    mMax = length(readout.beam1) - window
    out = map(m->correlate(readout, m, n, window), 1:mMax)

    return DenseCorrelationVector(n,out)
end
export correlate1d

#-------------------------------------------------------------------------------

"""
    function correlate(readout::DenseReadout, m::Int, n::Int, window::Int)

Calculates the mn-th element of the correlation matrix between beam 1 and 2.
"""
function correlate(readout::DenseReadout, m::Int, n::Int, window::Int)
    @assert m <= length(readout.beam1) - window "m index out of bounds"
    @assert n <= length(readout.beam1) - window "n index out of bounds"

    v1 = view(readout.beam1,m:(m+window))
    v2 = view(readout.beam2,n:(n+window))
    
    num = dot(v1,v2)/window
    den = mean(v1)*mean(v2)
    return num/den
end

"""
	function correlate(u::Vector{T},v::Vector{T},offset::Integer,window::Integer = -1) where {T<:Number}

Calculates correlation between vectors u and v with given offset. Specify averaging window to limit range of correlation. If the window extends beyond the end of one vector, it treats out-of-bounds indices as zero.
"""
function correlate(u::Vector{T},v::Vector{T},offset::Integer,window::Integer = -1) where {T<:Number}

	@assert offset <= length(u) "Offset out of bounds"
	@assert window <= length(u) && window <= length(v) "Window must be smaller than input vector lengths"

	if window == -1
		window = length(u)
	end

	v1 = view(u,1:window)
	v2 = view(v,1+offset:min(window+offset,length(v)))
	if window+offset > length(v)
		v2 = vcat(v2,zeros(window+offset-length(v)))
	end

	return dot(v1,v2)/window
end

export correlate

#-------------------------------------------------------------------------------

"""
	function autocorrelate(u::Vector{T},offset::Integer, window::Integer = -1) where {T<:Number}

Calculates correlation of vector u with itself.
"""
function autocorrelate(u::Vector{T},offset::Integer, window::Integer = -1) where {T<:Number}
	correlate(u,u,offset,window)
end

export autocorrelate

#-------------------------------------------------------------------------------

"""
	ncorrelate(u::Vector{T},v::Vector{T},offset::Integer,window::Integer = -1) where {T<:Number}

Calculates normalized correlation between vectors u and v with given offset. Specify averaging window to limit range of correlation. If the window extends beyond the end of one vector, it treats out-of-bounds indices as zero.
"""
function ncorrelate(u::Vector{T},v::Vector{T},offset::Integer,window::Integer = -1) where {T<:Number}

	@assert offset <= length(u) "Offset out of bounds"
	@assert window <= length(u) && window <= length(v) "Window must be smaller than input vector lengths"

	if window == -1
		window = length(u)
	end

	v1 = view(u,1:window)
	v2 = view(v,1+offset:min(window+offset,length(v)))
	if window+offset > length(v)
		v2 = vcat(v2,zeros(window+offset-length(v)))
	end

	return dot(v1,v2)/(window*mean(v1)*mean(v2))
end

export ncorrelate

#-------------------------------------------------------------------------------

"""
    function corrTimes(τ::Vector,γCounts::Vector)

Returns an array of τ values for which the γCounts autocorrelation is non-zero.
"""
function corrTimes(τ::Vector,γCounts1::Vector,γCounts2::Vector)
    @assert length(γCounts1) == length(γCounts2) "Count vectors must have the same length"
    return τ[map(i->correlate(γCounts1,γCounts2,i,length(τ)) > 0 ? true : false,collect(0:length(τ)-1))]
end

export corrTimes

#-------------------------------------------------------------------------------

"""
    function autocorrTimes(τ::Vector,γCounts::Vector)

Returns an array of τ values for which the γCounts autocorrelation is non-zero.
"""
function autocorrTimes(τ::Vector,γCounts::Vector)
    return τ[map(i->autocorrelate(γCounts,i,length(τ)) > 0 ? true : false,collect(0:length(τ)-1))]
end

export autocorrTimes

#-------------------------------------------------------------------------------

"""
    γCoincidentτ(beam1::Vector{T}, beam2::Vector{T}) where {T<:Number}

Returns the separation index between the first coincident counts of beam 1 and 2.
Returns -1 if no coincident counts are calculated.
"""
function γCoincidentIndex(beam1::Vector{T}, beam2::Vector{T}) where {T<:Number}
    @assert length(beam1) == length(beam2) "Beam vectors must have equal length"
    # iterate through beam 1
    for (i,nbar) in enumerate(beam1)
        if poissonCount(nbar) > 0
            for (j,mbar) in enumerate(beam2[i:end])
                if poissonCount(mbar) > 0
                    # return the separation index between coincident counts
                    return j-1 # subtract 1 since that's where Julia indexing starts
                end
            end
            return -1
        end
    end
    return -1
end

export γCoincidentIndex

#-------------------------------------------------------------------------------

"""
    firstNonzero(a::Vector)

Returns the index of the first non-zero element in a.
Returns -1 if all values are zero.
"""
function firstNonzero(a::Vector)
    for (i,val) in enumerate(a)
        if val != 0
            return i
        end
    end
    return -1
end

export firstNonzero
################################################################################
# Counting Related Functions
################################################################################
"""
    function countDeltaTimes(τ::Vector,γCounts::Vector)

Returns an array of τ values for which the γCounts autocorrelation is non-zero.
"""
function countTimes(times::Vector,γCounts::Vector)
    out = Vector{Real}(undef,0)
    for (i,counts) in enumerate(γCounts)
        if counts != 0
            countTimes = times[i]*ones(counts)
            out = vcat(out,countTimes)
        end
    end
    return out
end

export countTimes

################################################################################
# Functions from MGST2021
################################################################################

"""
    stauAvg(τ::Number,params::eFieldParams,n::Integer)

Returns the average of the Doppler noise term
"""
function stauAvg(τ::Number,source::LightSource)
    return stauAvg(τ,source.σ,source.n)
end

function stauAvg(τ::Number,params::SpeckleParams)
    return stauAvg(τ,params.σ,params.n)
end
"""
    stauAvg(τ::Number,σ::Number,n::Integer)

Returns the average of the Doppler noise term
"""
function stauAvg(τ::Number,σ::Number,n::Integer)
    term1 = n
    term2 = n*(n-1)
    term2 *= exp(-σ^2*τ^2)
    return term1 + term2
end

export stauAvg

"""
    stauVar(τ::Number,source::LightSource)

Returns the variance of the Doppler noise term
"""
function stauVar(τ::Number,source::LightSource)
    return stauVar(τ,source.σ,source.n)
end

"""
    stauVar(τ::Number,σ::Number,n::Integer)

Returns the variance of the Doppler noise term
"""
function stauVar(τ::Number,σ::Number,n::Integer)
    στ2 = σ^2*τ^2

    prod1 = 8*n*(n-1)
    prod1 *= exp(-2*στ2)

    prod2 = n-1+cosh(στ2)

    prod3 = sinh(στ2/2)^2

    return prod1*prod2*prod3
end

export stauVar

"""
    stau(τ::Number,field::eField)

Calculates the value of the Doppler noise term for the given eField
"""
function stau(τ::Number,field::eField)
    return stau(τ,field.νn)
end

"""
    stau(τ::Number,νn::Vector)

Calculates the value of the Doppler noise term for the given τ and frequencies
"""
function stau(τ::Number,νn::Vector)
    terms = exp.(-2.0*π*im*τ*νn)
    sumterms = sum(terms)
    return real(sumterms*conj(sumterms))
end

export stau

function Δm(params::SpeckleParams)
    return params.νm .- mean(params.νm)
end

"""
    g2Calc(τ::Number,params::SpeckleParams)

Returns the average calculated value of g2(τ) from MGST2021.
"""
function g2Calc(τ::Number,params::SpeckleParams)
    n = params.n
    em2 = real.(params.Em .* conj.(params.Em))
    em4 = em2 .* em2
    sumEm2 = sum(em2)
    sumEm4 = sum(em4)

    g2τ = 1.0

    term2 = -sumEm4/(n*sumEm2^2)

    g2τ += term2

    term3 = sum(em2 .* exp.(-im*2*π*τ*Δm(params)))/sumEm2
    term3 *= conj(term3)
    term3 = real(term3)
    term3 *= stauAvg(τ,params)/n^2

    return g2τ+term3
end
"""
    g2Calc(τ::Number,n::Integer,source::LightSource)

Returns the average calculated value of g2(τ) from MGST2021.
"""
function g2Calc(τ::Number,source::LightSource)
    n = source.n
    em2 = real.(source.Em .* conj.(source.Em))
    em4 = em2 .* em2
    sumEm2 = sum(em2)
    sumEm4 = sum(em4)

    g2τ = 1.0

    term2 = -sumEm4/(n*sumEm2^2)

    g2τ += term2

    term3 = sum(em2 .* exp.(-im*2*π*τ*Δm(source)))/sumEm2
    term3 *= conj(term3)
    term3 = real(term3)
    term3 *= stauAvg(τ,source)/n^2

    return g2τ+term3
end

export g2Calc
