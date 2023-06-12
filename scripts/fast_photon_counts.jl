using Distributed
@everywhere using DrWatson

@everywhere begin
    @quickactivate :Speckletroscopy
    using SharedArrays
end


light_source_dict = Dict(
    :n  => 100,
    :Em => [1.0,1.0,1.0],
    :νm => [456810,456813,456815], #GHz
    :σ  => 55.0,
    :fγ => 2.0e6
)

light_source = LightSource(light_source_dict)
