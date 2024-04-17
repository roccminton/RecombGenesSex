using CairoMakie
using JLD
using Statistics
using DataFrames
using LaTeXStrings
using BasicInterpolators

include("tools.jl")

#Change to your local directory
abs_path = "/home/larocca/github/RecombGenesSex/"

#=
Remark: The data required to generate the figure is not included in this repository
due to its large size. However, you can generate the data on your own by executing
the following for loop, while including the 'runandsaveIBM.jl' file.

    for N in 100:100:1000
        for dni in 0.01:0.01:0.1
            dosimulation(N,dni,"N=$N,dni=$dni",abs_path)
        end
    end
    
However this is not recommended to execute this simulations in series, but rather in parallel.
Otherwise executing the above for loop takes a very very long time.
=#


#---

#Font size 
fontsize_theme = Theme(fontsize = 25)
set_theme!(fontsize_theme)

Ns_data = 100:100:1000
μs_data = 0.01:0.01:0.1


#loading the data will not work on other machines, see above remark
data_path = "/media/larocca/PortableSSD/Data/NoRecombination/BigN/"
#make this true to re-load data
if false
    #setup empty DataFrame to store Data
    df = DataFrame(N=Int[],μ=Float64[],i=Int[],r=Union{Nothing,Float64}[])

    #populate DataFrame with statistics
    R = Vector{Int}(undef,3)
    println("Preparing Data...")
    for N in Ns_data
        println("Currently at N=$N")
        @showprogress for μ in μs_data
            for i in 1:3
                data = load(data_path * "N=$N,dni=$(μ)_$i.jld")
                v = data["ML"] ./ data["PopSize"]
                R[i] = findinc(v)
            end
            n = 3-count(iszero,R)
            t = (iszero(n) ? nothing : sum(R)/n)
            push!(df,[N,μ,n,t])
        end
    end
end

#set data
Ns = 100:1050
μs = 0.01:0.0001:0.105

qt = 0.99

mC0 = load(abs_path * "data/meanC0extimes.jld")

#C = [C₀(N,μ) for N in Ns, μ in μs]
p = BicubicInterpolator(mC0["Ns"],mC0["dnis"],mC0["C0"],NoBoundaries())
C = [p(N,μ) for N in Ns, μ in μs]
x = sort!(reduce(vcat,[fill(C₀(row.N,row.μ),row.i) for row in eachrow(df)]))
q = quantile(x,qt,sorted=true)

#generate figure

f = Figure(size=(1200,1000))

ax, hm = heatmap(
    f[1,1],mC0["Ns"],mC0["dnis"],log10.(mC0["C0"]./(2*mC0["K"])),
    )
ax.title = L"$ $Equilibrium fraction of $C_0$ and its extinction events"
ax.xlabel = L"N"
ax.ylabel = L"2\mu"
ax.xtickformat = vs -> map(v->L"%$(round(Int,v))",vs)
ax.ytickformat = vs -> map(v->L"%$(round(v,digits=2))",vs)
ylims!(ax,(minimum(μs),maximum(μs)))

tickrange = -1:-1:-5

cb = Colorbar(
    f[:,end+1],hm;
    ticks = (tickrange,[L"10^{%$x}" for x in tickrange])
    )

lines!(ax,Ns,μ_from_C.(q,Ns),linestyle=:dashdot,color=:darkblue)

scatter!(ax,Point2f(600,0.05),color=:darkblue,markersize=3*17)

for row in eachrow(df)
    if !iszero(row.i)
        scatter!(
            ax,Point2f(row.N,row.μ),
            color = [round(Int,row.r/10^5 * 100)],
            colormap = :buda,
            colorrange = (1,100),
            markersize=row.i*15
            )
    end
end


# tN = 420
# tμ = 0.051
# text!(ax,tN,tμ;
#     text=L"C_0 = %$(round(q*10^3;digits=2)) \times 10^{-3}",
#     rotation = -0.55,
# #    offset = (-10,0),
#     color = :darkblue,
#     fontsize = 12
#     )

sub_f = f[end+1,1:end-1] = GridLayout()

Label(sub_f[1,1],L"$ $Mean Extinction")
Label(sub_f[2,1],L"Time in $10^3$")
Label(sub_f[3,1],L"$ $generations")

Colorbar(sub_f[:,2],
    limits = (0,10^5), colormap = :buda,
    flipaxis=false,ticks=(range(0,10^5;length=3),[L"0",L"50",L"100"]),
    vertical = false
    )

rowgap!(sub_f,1,0)
rowgap!(sub_f,2,0)
colsize!(sub_f,2,Relative(0.9))

rowgap!(f.layout,1,-10)

#uncomment to save the output
#save(abs_path * "figures/figure2.pdf", f, pt_per_unit=1)

f
