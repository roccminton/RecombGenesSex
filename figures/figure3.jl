using CairoMakie
using JLD
using Statistics
using LaTeXStrings
using CurveFit

include("tools.jl")

#Change to your local directory
abs_path = "/home/larocca/github/RecombGenesSex/"

#Font size 
fontsize_theme = Theme(fontsize = 17)
set_theme!(fontsize_theme)

#---

data0 = load(abs_path * "data/data_r=0.jld")
data1 = load(abs_path * "data/data_r=1.jld")
dataMB0 = load(abs_path * "data/data_MB_r=0.jld")

Ns = data0["Ns"]
μ = data0["dni"] / 2


function addtext(g,i,j,text,diff)
	axtext = Axis(g[i,j])
	xlims!(axtext,(0,1))
	ylims!(axtext,(0,1))
	hidedecorations!(axtext)
	hidespines!(axtext)
	text!(axtext,diff,1-diff;text=text,font=:bold,fontsize=26,align=(:left,:top))
end

function splityaxis!(f,Ns,i,j,dataMB0,data1)

	diff = 2.5

	topmin = 40 - diff
	topmax = 80 + diff

	s_topmin = 38.5
	s_topmax = 50

	bottommax = 11

	scale(x) = x/2

	#shrink interval [topmin,topmax] to [s_topmin,s_topmax]
	trans(x) = x * ((s_topmax - s_topmin)/(topmax-topmin)) + s_topmin - ((s_topmax - s_topmin)/(topmax-topmin)) * topmin

	lims = Observable(((0.0, bottommax), (trans(topmin), trans(topmax))))

	g = f[i, j] = GridLayout()

	ax_top = Axis(
		f[i, j][1, 1],
		yticks = (
			range(trans(topmin+diff),trans(topmax-diff),length=3),
			[L"%$(round(Int,topmin+diff))",L"%$(round(Int,median([topmax,topmin])))",L"%$(round(Int,topmax-diff))"]
			)
		)
	ax_bottom = Axis(
		f[i, j][2, 1],
		yticks = (0:5:10,[L"0",L"5",L"10"]),
		xticks = (0:500:1000,[L"0",L"500",L"1000"]),
		xlabel = L"$ $Number of Genes",
		ylabel = L"$ $Haploid Mutation Burden"
		)

	on(lims) do (bottom, top)
	    ylims!(ax_bottom, bottom)
	    ylims!(ax_top, top)
	    rowsize!(g, 1, Auto(top[2] - top[1]))
	    rowsize!(g, 2, Auto(bottom[2] - bottom[1]))
	end

	hidexdecorations!(ax_top, grid = false)
	ax_top.bottomspinevisible = false
	ax_bottom.topspinevisible = false

	xlims!(ax_bottom,(0,Ns[end]+50))
	linkxaxes!(ax_top, ax_bottom)
	rowgap!(g, 10)

	angle = pi/8
	linelength = 30

	segments = lift(
	        @lift($(ax_top.yaxis.attributes.endpoints)[1]),
	        @lift($(ax_bottom.yaxis.attributes.endpoints)[2]),
	        @lift($(ax_top.elements[:yoppositeline][1])[1]),
	        @lift($(ax_bottom.elements[:yoppositeline][1])[2]),
	    ) do p1, p2, p3, p4
	    ps = Point2f[p1, p2, p3, p4]

	    map(ps) do p
	        a = p + Point2f(cos(angle), sin(angle)) * 0.5 * linelength
	        b = p - Point2f(cos(angle), sin(angle)) * 0.5 * linelength
	        (a, b)
	    end
	end

	linesegments!(f.scene, segments,color=:black)

	#plot burden line at 2N*h
	fine_Ns = 0:Ns[end]+100
	lines!(ax_bottom,fine_Ns,scale.(ML.(fine_Ns,μ)),color=:red)
	#plot high mutation burden if exists
	for (n,N) in enumerate(Ns)
		loads = collect(skipmissing(dataMB0["highMB"][n,:]))
		!isempty(loads) && scatter!(
			ax_top,fill(N,length(loads)),trans.(scale.(loads)),
			color = :darkblue,
			marker = :diamond
		)
	end
	#plot regular mutation burden
	tend = data1["extimes"][1]
    scatter!(
        ax_bottom,Ns,
		scale.([
			mean(
				dataMB0["meanMB"][n,i] for i in findmaxindex(view(dataMB0["extimes"],n,:),tend)
			) for n in 1:length(Ns)
		]),
        color = :darkblue,
    )
	scatter!(ax_bottom,Ns,scale.(data1["meanMB"]),color = :darkorange)
	lines!(
		ax_bottom,Ns,scale.(data1["meanMB"]),
		color = :darkorange,linestyle = :dash
	)

	notify(lims)

	addtext(f,i,j,L"$ $A",0.01)

	return f
end

function c0plot!(g,Ns,i,j,data0,data1)
	#set the y scale (remember to adjust yticks accordingly)
	scale(x) = log10(x)

	ax_C0 = Axis(
	    g[i,j],
	    ylabel=L"$ $Fraction of mutation free gamete",xlabel=L"$ $Number of genes",
		yaxisposition = :right,
		ytickformat = vs -> [L"10^{%$(round(v,digits=1))}" for v in vs],
		xtickformat = vs -> [L"%$(round(Int,v))" for v in vs]
	)

	hidexdecorations!(ax_C0,grid=false)

	tend = data1["extimes"][1]
	#scatter the r=0 points
	scatter!(
		ax_C0,Ns,
		scale.([
			mean(
				data0["meanC0"][n,i] for i in findmaxindex(view(data0["extimes"],n,:),tend)
			) for n in 1:length(Ns)
		]),
		color = :darkblue
		)
	#scatter the r=1 points
	scatter!(
			ax_C0,Ns,scale.(data1["meanC0"]),
			color = :darkorange
			)
	#dashed line for r=0
	lines!(
		ax_C0,Ns,
		scale.([
			mean(
				data0["meanC0"][n,i] for i in findmaxindex(view(data0["extimes"],n,:),tend)
			) for n in 1:length(Ns)
		]),
		color = :darkblue,linestyle=:dash
		)
	#dashed line for r=1
	lines!(
			ax_C0,Ns,scale.(data1["meanC0"]),
			color = :darkorange,linestyle=:dash
			)
	#set the axis lims
	xlims!(ax_C0,(0,Ns[end]+50))
	ylims!(ax_C0,high=log10(1.0))

	addtext(g,i,j,L"$ $B",0.01)
end

function varplot!(g,Ns,i,j,data0,data1)
	scale(x) = x

	ax_var = Axis(
	    g[i,j],
	    ylabel=L"Variance of Hapolid Mutation Burden$",xlabel=L"$ $Number of genes",
		yaxisposition= :right,
		xtickformat = vs -> [L"%$(round(Int,v))" for v in vs],
		ytickformat = vs -> [L"%$(round(Int,v))" for v in vs]
	)

	#hidexdecorations!(ax_var,grid=false)

	fine_Ns = 0:Ns[end]+100
	tend = data1["extimes"][1]
	var0 = scale.([mean(
			data0["varhMB"][n,i] for i in findmaxindex(view(data0["extimes"],n,:),tend)
		) for n in 1:length(Ns)])

	scatter!(ax_var,Ns,var0,color = :darkblue)


	linear(a,b,x) = a+b*x
	lines!(
		ax_var,fine_Ns,linear.(linear_fit(Ns,var0)...,fine_Ns),
		color = :darkblue, linestyle = :dash
		)

	scatter!(ax_var,Ns,scale.(data1["varhMB"]),color = :darkorange)
	lines!(ax_var,Ns,scale.(data1["varhMB"]),color = :darkorange,linestyle = :dash)

	ylims!(ax_var,(0,11))
	xlims!(ax_var,(0,Ns[end]+50))

	addtext(g,i,j,L"$ $C",0.01)
end

function addLegend(g,i,j)
	elem_1 = [LineElement(color=:darkblue,linestyle=:dash),
			MarkerElement(color=:darkblue,marker=:circle)
		]
	elem_2 = [LineElement(color=:darkorange,linestyle=:dash),
			MarkerElement(color=:darkorange,marker=:circle)
		]
	elem_3 = [LineElement(color=:red)]
	Legend(g[i,j],
		[elem_1,elem_2,elem_3],
		[L"$ $No Recombination", L"$ $Full Recombination", L"N\sqrt{1-e^{-\frac{\mu}{N}}}"],
		tellheigth = false, tellwidth = false,
		margin = (40,10,10,40), halign = :left, valign =:top,
		patchsize = (35,35)
	)
end


#---

g = Figure(size=(1200,600))

splityaxis!(g,Ns,1:2,1,dataMB0,data1)
c0plot!(g,Ns,1,2,data0,data1)
varplot!(g,Ns,2,2,data0,data1)
addLegend(g,1:2,1)

#uncomment to save output
#save(abs_path * "figures/figure3.pdf",g)

g
