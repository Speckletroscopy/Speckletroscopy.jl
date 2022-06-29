################################################################################
# Fold Iterator
################################################################################

"""
    Fold(v::AbstractVector,window::Int,start::Int)

Returns an iterator to step through slices of vector v with size window.
Optionally include a starting index to begin folding at a value other than index = 1.
"""
struct Fold{T<:AbstractVector}
    v::T
    window::Int
    start::Int
end

function Fold(v::T,window::Int) where {T<:AbstractVector}
    return Fold(v,window,1)
end

export Fold

Base.iterate(F::Fold) = (F.v[F.start:(F.start+F.window-1)], F.start+F.window)
Base.iterate(F::Fold, state) = state+F.window-1 > length(F.v) ? nothing : (F.v[state:(state+F.window-1)],state+F.window)

################################################################################
# Fourier Transforms
################################################################################
"""
    fftPositiveFreq(f,ft)

Returns (frequencies, coefficients) for positive frequencies of the fourier transform
"""
function fftPositiveFreq(f::AbstractVector{T},ft::AbstractVector{T}) where {T<:Real}
	@assert length(ft) == length(f) "Fourier transform and frequency vectors much have matching length"
	ftRange = f .>= 0
	return (f[ftRange],ft[ftRange])
end

export fftPositiveFreq

"""
    meanFFT(fft::Vector{S},cuts::Tuple{T,T}) where {S<:Real, T<:Integer}

Takes mean before taking fourier transform
"""
function meanFFT(timeSeries::Vector{S},cuts::Tuple{T,T}) where {S<:Real, T<:Integer}
	tsPrep = timeSeries[cuts[1]:cuts[2]]
	tsPrepMean = mean(tsPrep)
	tsPrep = tsPrep .- tsPrepMean
	return FFTW.fft(tsPrep)
end

export meanFFT

"""
    fftFreq(tres::Real,τ::Vector,cuts::Tuple{T,T}) where {T<:Integer}

Returns frequencies associated with Fourier transform 
"""
function fftFreq(tres::Real,cuts::Tuple{T,T}) where {T<:Integer}
	τlen = cuts[2]-cuts[1]+1
	freqFFT = FFTW.fftfreq(τlen,1/tres)
end

export fftFreq

"""
    niceFFT(timeSeries::Vector{S},cuts::Tuple{T,T},params::SpeckleParams) where {S<:Real, T<:Integer}

Performs some light data processing for a nice looking correlation fft.
    - Does a low cut on the time series to get past the non-periodic part of g2τ
    - Subtracts the mean of the time series to minimize the constant term coefficient

Returns (frequencies, coefficients)
"""
function niceFFT(timeSeries::Vector{S},params::SpeckleParams) where {S<:Real}
    # low cut 5 efolds beyond the non-periodic part of g2τ
    τlow = sqrt(5)/params.σ
    ilow = convert(Int,ceil(τlow/params.resolution))
    cuts = (ilow,length(timeSeries))
    mfft = meanFFT(timeSeries,cuts)
    mfft = abs.(mfft)
    mfft = convert(Vector{Float64},mfft)
    freqs = fftFreq(params.resolution,cuts)
    return fftPositiveFreq(freqs,mfft)
end
export niceFFT

struct SpeckleFFT{V<:Vector}
    id::UUID
    freqs::Vector
    singles::Vector{V}
    sumFFT::Vector
    FFTsum::Vector
end
function SpeckleFFT(sim::SpeckleSim{U,V}) where {U<:SpeckleReadout,V<:Correlation}
    corrFFTvec = map(corr->niceFFT(corr.data,sim.params)[2],sim.corr)
    sumCorrFFT = +(corrFFTvec...)
    corrSum = +(getproperty.(sim.corr,:data)...)
    freqs,corrSumFFT = niceFFT(corrSum,sim.params)

    return SpeckleFFT(sim.id,freqs,corrFFTvec,corrSumFFT,sumCorrFFT)
end
export SpeckleFFT
################################################################################
# Linear Prediction
################################################################################
function LPCoeff(corr::CorrelationVector)
    p = 20
    return DSP.LPC.lpc(corr.data,p)
end
################################################################################
# Analysis functions
################################################################################

function snr(ft::Vector{T},freqs::Vector{T},ishifts::Vector{U},sym::Vector{V}) where {T<:Number,U<:Tuple,V<:Symbol}
    @assert length(ishifts) == length(sym)+1 "Slicing went wrong!"
    # slice out the peaks and put the noise in a single view
    freqviews = map(slice->view(freqs,UnitRange(slice[1],slice[2])),ishifts)
    freqviews = CatView(freqviews...)
    noiseviews = map(slice->view(ft,UnitRange(slice[1],slice[2])),ishifts)
    noiseviews = CatView(noiseviews...)

    # linear least squares fit for the noise
    linearModel(ν,p) = p[1]*ν .+ p[2] # specify a linear model
    p0 = [0.0,mean(noiseviews)] # initial guess
    fit = curve_fit(linearModel,freqviews,noiseviews,p0) 
    # shift noise to horizontal axis
    noiseShift =linearModel(freqviews,fit.param)
    noiseviews = noiseviews - noiseShift

    # make a histogram of the noise
    noisehist = StatsBase.fit(Histogram,noiseviews,nbins = convert(Int,floor(length(noiseviews)/10)))

    # get max of histogram
    maxnoise = findmax(noisehist.weights)

    # find width of peak at half max value
    halfmax = maxnoise[1]/2
    ilhalf = maxnoise[2]
    irhalf = maxnoise[2]
    while ilhalf > 1 && noisehist.weights[ilhalf] > halfmax
        ilhalf -= 1
    end
    while irhalf < length(noisehist.weights) && noisehist.weights[irhalf] > halfmax
        irhalf += 1
    end
    edges = collect(noisehist.edges[1])
    fwhm = edges[irhalf] - edges[ilhalf]

    σ = fwhm/2.355
    # select the max value in each peak window as the peak height
    peaks = Float64[]
    ipeaks = Int[]
    for k=1:length(sym)
        peakview = view(ft,UnitRange(ishifts[k][2],ishifts[k+1][1]))
        peakval,ipeak = findmax(peakview)
        ipeak = ishifts[k][2]+ipeak-1
        push!(peaks,peakval)
        push!(ipeaks,ipeak)
    end
    # subtract the mean value of the noise from the peak height
    peaks = peaks - map(ind->linearModel(freqs[ind],fit.param),ipeaks)
    # calculate snr
    peaks = peaks/σ

    # we don't have a theoretical value for the peak height, so take snr=1 if fluctuations are larger than the signal
    peaks = map(peak-> peak < 1 ? 1 : peak,peaks)

    # zip together shifts and snrs for reference
    return collect(zip(sym,peaks))
end

"""
    snr(sfft::SpeckleFFT,params::SpeckleParams)

Calculates the signal-to-noise ratio for each peak in sfft
"""
function snr(sfft::SpeckleFFT,params::SpeckleParams)

    if length(params.νm) > 1
        halfwindow = 2
        # find frequency shifts
        shifts = map(x->abs(x[2]-x[1]),subsets(params.νm,2))
        shifts = sort(shifts)
       
        # creat a dictionary to hold everything
        shiftSym = Symbol.(shifts)
        snrDict = Dict{Symbol,Union{Vector{String},Vector{Float64}}}()
        for symbol in shiftSym
            snrDict[symbol] = Float64[]
        end
        snrDict[:type] = String[]

        # find fft indices where frequency shifts are supposed to be
        ishifts = Tuple{Int64,Int64}[]
        i = 1
        istart = i
        iend = i
        for shift in shifts
            # iterate until we hit a shift or run out of frequencies
            while i > length(sfft.freqs) || sfft.freqs[i] < shift 
                i+=1
            end

            # back up one index plus half the size of the cut around the peak
            itest = i-1-halfwindow

            # check if itest is valid
            if i > length(sfft.freqs)
                iend = length(sfft.freqs)
            elseif itest < istart
                iend = istart
            else
                iend = itest
            end
            # store range for later use
            push!(ishifts,(istart,iend))
            # new starting point is on the other side of the peak
            istart = iend + 2*halfwindow
        end
        @assert i < length(sfft.freqs) "Counting should have stopped before the final shift"
        # store final slice 
        push!(ishifts,(istart,length(sfft.freqs)))
        @assert length(ishifts) == length(shifts)+1 "Slicing went wrong!"
        
        # calculate the snr for each fourier transform
        # collect results in a single dictionary
        singleSnr = map(single->snr(single,sfft.freqs,ishifts,shiftSym),sfft.singles)
        for snrVec in singleSnr
            push!(snrDict[:type],"single")
            for (name,val) in snrVec
                push!(snrDict[name],val)
            end
        end
        sumFFTsnr = snr(sfft.sumFFT,sfft.freqs,ishifts,shiftSym)
        push!(snrDict[:type],"sumFFT")
        for (name,val) in sumFFTsnr
            push!(snrDict[name],val)
        end
        FFTsumSnr = snr(sfft.FFTsum,sfft.freqs,ishifts,shiftSym)
        push!(snrDict[:type],"FFTsum")
        for (name,val) in FFTsumSnr 
            push!(snrDict[name],val)
        end

        # convert dictionary to table and return
        return DataFrame(snrDict)

    else
        return nothing
    end
end
