module Micropalaz

using CSV
using DataFrames
using Query
using OnlineStats
using StatsBase
using Plots 
theme(:mute)

function diversityData(hole_id::String, fossil_group::String, div_file::String)
    """
    Function taking three parameters.
    hole_id : Hole ID in the form of e.g. "177_1090B" 
    fossil_group : As it is in Neptune Database, "D" for diatoms, "R" for Rads
    div_file : Exported diversity file from Neptune
    """
    div_data = CSV.read(div_file, DataFrame)
    data = @from i in div_data begin
        @where i.hole_id == hole_id 
        @where i.fossil_group == fossil_group
        #@se
        @select {i.taxon_name, i.sample_age_ma, i.value}
        @collect DataFrame
    end
    return data
    #return("Total number of taxa in this site is $(length(data[!, :taxon_name]))")
end

function dataPlot(data)
    """
    Plotting diveristy data age - total num. of species 
    parameter is the output of diversityData function
    """
    age = data[!, :sample_age_ma]
    a = sort(countmap(age::Array{Float64,1}))
    scatter(collect(keys(a)), collect(values(a)), alpha = 0.7)
    plot!(collect(keys(a)), collect(values(a)), alpha = 0.5)
    return collect(keys(a)), collect(values(a))
    #plt.plot.vlines(33.9)
end

function diversity(age, data::DataFrame)
    """
    Return single sample (age) data
    """
    siteData = @from i in data begin
        @where i.sample_age_ma == age
        @select {i.taxon_name, i.sample_age_ma, i.value}
        @collect DataFrame
    end
    return siteData
end

function SiteChao1(hole_id::String, fossil_group::String, div_file::String)
    """
    Return Age-Chao1 estimation for Site -of interest.
    """
    divdata = diversityData(hole_id, fossil_group, div_file)
    my_vec = []
    for i in sort(unique(divdata.sample_age_ma)) 
        my = diversity(i, divdata)
        push!(my_vec, neptuneChao1(my))   
    end
    return sort(unique(divdata.sample_age_ma)), my_vec
end


function plot_Chao1(c1_site_data, raw_site_data, leg_pos = :bottomleft)
    """
    Plot Chao1 estimations and raw diversity.
    """
    plot(raw_site_data[1], raw_site_data[2], fillrange = c1_site_data[2], 
    fillalpha = 0.35, c = 1, label = "Chao1 richness estimation",
    legend = leg_pos, xlabel = "Age (Ma)", ylabel = "Species Richenss", dpi = 300)
    scatter!(raw_site_data[1], raw_site_data[2], label = "Raw richness", alpha = 0.1)
    plot!(size=(500, 300))
end


function neptuneObserved(data::DataFrame)
    """
    return observed diversity
    a: is the example data
    """
    return nrow(data)
end

function neptuneGoodsCov(data)
    """
    Good's coverage of a sample
    """
    return 1 - (sum(data.value .== 1) / sum(data.value))
end

function neptuneChao1(data)
    """
    Return Chao1 richness estimation
    """
    return neptuneObserved(data) + (sum(data[!, :value] .== 1))^2 / (2 * sum(data[!, :value] .== 2))  
end

function neptuneACE(data, treshold::Int64)
    """
    Returns ACE richness estimation.
    """
    data = data[data.value .> 0, :]
    S_abun = sum(data.value .> treshold) #rechness of abundant taxa
    S_rare = sum(data.value .<= 10) #richness of rare taxa
    singleton = sum(data.value .== 1) #num. of singleton taxa
    rare_data = filter(x -> x.value .<= treshold, data)
    N_rare = sum(rare_data.value) #abundance of rare individuals
    C_ace = 1 - (singleton / N_rare)
    z = [1:10;]
    a_1 = [sum(a.value .== i) for i in 1:10]
    f_1 = (z .* (z .- 1)) .* a_1
    G_ace = (S_rare / C_ace) * sum(f_1) / (N_rare * (N_rare - 1))
    S_ace = S_abun + (S_rare / C_ace) + (singleton / C_ace) * max(G_ace, 0)
end

function RAC(data)
    """
    Returns Rank Abundance for given sample 
    """
    vals = data.value
    vals_ab_rank = sort(data.value, rev = true)
    return vals_ab_rank
end


"""
Return Rank Abundance plot
"""
RACplot(x) = scatter(x, xlabel = "Rank in abundance", ylabel = "Abundance", color = x, alpha =0.7, legend = false, theme= :mute)



end # module
