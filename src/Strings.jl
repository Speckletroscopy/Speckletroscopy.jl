################################################################################
# String management for Speckletroscopy.jl
################################################################################


"""
    makeName(params::Dict; excpt::Dict = Dict())

Returns a string based on parameters in the simulation dictionary.
"""
function makeName(params::Dict; excpt::Dict = Dict())
    date = today()
    mm = length(month(date)) == 1 ? string(0,month(date)) : month(date)
    dd = length(day(date)) == 1 ? string(0,day(date)) : day(date)
    out = string(year(date),mm,dd,"_")
    for key in keys(params)
        keyname = key
        val = params[key]
        if isa(val,Array)
            val = length(val)
            keyname = string("len-",keyname)
        end
        out = string(out,keyname,"=",val,"_")
    end
    return out
end

export makeName

"""
    function paramTable(params::Dict)

Returns a markdown table with given parameters
"""
function paramTable(params::Dict,names::Dict = Dict())
    out = "| Description | Value(s) |\n|:---:|:---:|\n"
    for (key,value) in names 
        @assert key in keys(params) "Key $key not found in parameters"
        pval = params[key]
        if length(pval) == 1
            pval = pval[1]
        end
        out = string(out,"|$value|$pval|\n")
    end
    return out
end

export paramTable
