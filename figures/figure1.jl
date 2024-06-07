using CairoMakie
using Statistics
using Distributions
using JLD

include("tools.jl")

#Change to your local directory
abs_path = "/home/larocca/github/RecombGenesSex/"

#=
Remark: The data required to generate the figure is not included in this repository
due to its large size. However, you can generate the data on your own by executing
the following function, while including the 'runandsaveIBM.jl' file.

    dosimulation(600,0.05,"N=$N,dni=$dni",abs_path)

Be aware though that executing this function may take some time.
=#

fontsize_theme = Theme(fontsize = 30)
set_theme!(fontsize_theme)

#Plot generating functions
#---

function plot_A!(f,ts,te,ml,prev,extimes,times,t)
    #create two axis which share a x-axis
    axprev = Axis(
        f[1, 1],
        yticklabelcolor = :orange,
        yaxisposition = :right,
        xaxisposition = :top,
        xticksvisible = false,
        ytickformat = vs -> map(v->L"$%$(round(Int,v*100))$ %",vs),
        ylabel = L"$ $Prevalence",
    )
    axprev.xticks = (extimes, [L"\dagger" for i = 0:length(extimes)-1])

    axload = Axis(
        f[1, 1],
        yticklabelcolor = :red,
        ylabel = L"$ $Mutation Burden",
        xlabel = L"time in $10^3$ generations",
        ytickformat = vs -> map(v->L"%$(round(Int,v))",vs)
    )
    axload.xticks = (
        range(ts, te; length = 5),
        [
            L"%$(round.(Integer, x / 1000))" for
            x in range(ts, te; length = 5)
        ],
    )

    #set the xlims and ylims
    linkxaxes!(axload, axprev)
    xlims!(axload, (ts, te))
    xlims!(axprev, (ts, te))
    ylims!(axload, low = 0)
    ylims!(axprev, low = 0)

    #add vline at t for middle corr matrix and for frames T
    vlines!(axprev,[t],color=:black,linestyle=:dashdot,linewidth=1)
    for (i,T) in enumerate(times)
        vlines!(axprev,[T[1],T[end]],linestyle=:dashdot,color=:black,linewidth=1)
    end

    #plot data
    prevline = lines!(axprev,ts:te,prev[ts+1:te+1],color = :orange,label = L"$ $Prevalence")
    loadline = lines!(axload, ts:te, ml[ts+1:te+1], color = :red, label = L"$ $Mutation Burden")
end

function plot_B!(f,ts,te,d,N_max,extimes,times)
    #iterate over all time frames
    for (i,T) in enumerate(times)
        tbl = (x = repeat(T, inner=N_max),
            height = reduce(vcat, [view(d["LoadHist"],1:N_max,t) ./ sum(view(d["LoadHist"],1:N_max,t)) for t in T]),
            grp = repeat(1:N_max,outer=length(T))
        )
        #create axis
        ax = Axis(
            f[1, i],
            xticksvisible = false,
            ytickformat = vs -> map(v->L"$%$(round(Int,v*100))$ %",vs),
            yaxisposition = :right,
        )
        #filter for extimes in frame
        ex = filter(t->t in T, extimes)
        if isempty(ex)
            hidexdecorations!(ax)
        else
            fex = findfirst(x->x==ex[1],extimes)
            ax.xticks = (ex, [L"\dagger" for i in fex-1:fex+length(ex)-2])
        end
        #set ax lims
        xlims!(ax,(T[1],T[end]))
        ylims!(ax,(0,1))
        #square aspect ratio
        colsize!(f, i, Aspect(1, 1.0))
        #plot data
        barplot!(ax,tbl.x,tbl.height,
                stack = tbl.grp,
                color = tbl.grp,
                colormap = :nipy_spectral,
                colorrange = (1,N_max),
        )
        #add lines for extimes
        vlines!(ax,ex,color=:gray,linewidth=1,alpha=0.5)
        #add frame name
        poly!(ax,Rect(T[1]+70,0.92,375,0.07),color=:lightgray,strokecolor=:black,strokewidth=1)
        text!(ax,0,1,text=L"T_{%$i}",align=(:left,:top),space=:relative,offset = (5,-5),fontsize=25)
        #add decorations if applicable
        if isone(i)
            ax.ylabel = L"$ $Haploid Load Class Distribution"
            ax.yaxisposition = :left
            hideydecorations!(ax,label=false)
        elseif i==3
            hideydecorations!(ax,ticklabels=false,ticks=false)
        else
            hideydecorations!(ax)
        end
        hidespines!(ax)
    end
    #add color bar
    Colorbar(
        f[2, :], limits = (0, N_max),
        colormap=:nipy_spectral,label=L"$ $Mutation Burden per Genome",
        vertical = false, flipaxis = false,
        tickformat = vs -> map(v->L"%$(round(Int,v))",vs)
        )
end

function plot_C!(f,ts,te,cormatrixs,loadpos,times,t)

    titles = [
        L"Average over $t \in T_1",
        L"At time $t = %$(round(t * 10^(-4),digits=1)) \times 10^4 $",
        L"Average over $t \in T_3",
    ]

    max_freq = maximum(maximum.(loadpos))*1.1

    for (i,T) in enumerate(times)
        ax = f[1, i] = GridLayout()
        colsize!(f, i, Aspect(1, 1.0))
        ax1 = Axis(
            ax[1:2, 1],
            title = titles[i],
            ylabel = isone(i) ? L"$ $Correlation Matrix" : "",
        )
        heatmap!(ax1, cormatrixs[i], clim = (0, 1))
        hidexdecorations!(ax1)
        hideydecorations!(ax1, label = false)

        ax2 = Axis(
            ax[3, 1],
            ytickformat = vs -> map(v->L"$%$(round(Int,v*100))$ %",vs),
        )
        hidespines!(ax2)
        xlims!(ax2, (0.5, N + 0.5))
        ylims!(ax2, (0,max_freq))
        barplot!(ax2,loadpos[i],color = :darkblue)

        isone(i) ? hideydecorations!(ax2, ticklabels = false, ticks = false,grid=false) :
        hideydecorations!(ax2,grid=false)
        hidexdecorations!(ax2, label = false)

        rowgap!(ax, 1)
    end

    Colorbar(f[1, 4], limits = (0, 1),tickformat=vs->map(v->L"%$v",vs))

    #Add legend
    m = :rect
    Legend(
        f[2, 3:4],
        [MarkerElement(color = :darkblue, marker = m)],
        markersize = 20,
        [L"$ $Allele Frequency of Lethal Equivalent per Gene"],
        tellwidth = false,
        orientation = :horizontal,
        halign = :right,
        framevisible = false,
    )
end

function addLabel!(f)
    #Add A/B/C Labels
    axtext = Axis(f,bbox=BBox(0,1500,0,1600))
    xlims!(axtext,(0,1))
    ylims!(axtext,(0,1))
    hidedecorations!(axtext)
    hidespines!(axtext)

    for (h,l) in zip([0.995,0.65,0.305],[L"$ $A",L"$ $B",L"$ $C"])
        text!(axtext,0.02,h;text=l,font=:bold,fontsize=25,align=(:left,:top))
    end
end

function addBrackets!(f,ts,te,times)
    #add curly brackets and names
    axbox = Axis(f,bbox=BBox(105,1375,1020,1170))
    xlims!(axbox,(ts,te))
    ylims!(axbox,(0,1))
    hidespines!(axbox)
    hidedecorations!(axbox)
    for (i,T) in enumerate(times)
        bracket!(axbox,T[1],.5,T[end],.5;orientation=:down,text=L"T_{%$i}")
    end
end

#parameters
N=600
#start and end of plot window
ts = 0
te = 80_000
#reload?
reload = false

#Load Data
#---
if reload
    #loading the data will not work on other machines, see above remark
    d = load("/media/larocca/PortableSSD/Data/NoRecombination/LoadClassDistribution/N=600,dni=0.05_long_3.jld")

    #Collect Data
    #find extinction events
    extimes = findextimes(d)
    #set time frames for zoom and correlation matrixes
    t_C2 = 40_000
    times_AB = [2500:7500,12000:17000,extimes[end]-2500:extimes[end]+2500]
    times_C = copy(times_AB) ; times_C[2]=t_C2:t_C2

    #maximum number of mutation per haploid genome
    N_max = maximum(maximum(findlast(!iszero,d["LoadHist"][:,t]) for t in T) for T in times_AB)
    #collect corrmatrix and loadpositions
    println("Calculating covariance matrices ... ")
    cormatrixs = [meanr2(d,N,time_to_sampleinds(int,d["savesnap"])) for int in times_C]
    loadpos = vec.([mean(d["LoadPos"][:, int], dims = 2) ./ d["K"] for int in times_C])
    #reorder matrix and vector to better see cluster
    for i in 2:3
        cluster = findcluster(cormatrixs[i])
        cormatrixs[i] = sortcormatrix(cormatrixs[i],cluster)
        loadpos[i] = sortloadpos(loadpos[i],cluster)
    end
end
#Generate Figure
#---
println("Generating figure...")
f = Figure(size=(1500,1600))

fA = f[1,:] = GridLayout()
fB = f[2,:] = GridLayout()
fC = f[3,:] = GridLayout()

plot_A!(fA,ts,te,
    replace_NaN(d["ML"] ./ d["PopSize"]),
    replace_NaN(d["Ill"] ./ d["PopSize"]),
    extimes,times_AB,t_C2
    )
plot_B!(fB,ts,te,d,N_max,extimes,times_AB)
plot_C!(fC,ts,te,cormatrixs,loadpos,times_C,t_C2)
addLabel!(f)
addBrackets!(f,ts,te,times_AB)

#uncomment to save the output
#save(abs_path * "figures/figure1.pdf", f, pt_per_unit=1)

f
