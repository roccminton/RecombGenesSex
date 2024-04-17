#Helper Functions for Plotting

using ProgressMeter

replace_NaN(v) = map(x -> isnan(x) ? zero(x) : x, v)

extime(data, loadclass, birthtime, tend) = findfirst(
    x -> iszero(data["LoadHist"][loadclass+1,x]),
    birthtime:tend,
)

function findextimes(data, tend = 0)
    iszero(tend) && (tend = data["historylength"])
    extimes = [extime(data, 0, 1, tend)]
    lc = 1
    while !isnothing(extimes[end])
        ext = extime(data, lc, extimes[end], tend)
        isnothing(ext) && return extimes
        push!(extimes, ext + extimes[end])
        lc += 1
    end
end

function findcluster(corr,ε=0.3)
    N = size(corr)[1]
    cluster = []
    for c in 1:N-1
        for r in c+1:N
            if corr[c,r] ≥ 1-ε
                i = findfirst(x->c∈x,cluster)
                if isnothing(i)
                    push!(cluster,[c,r])
                elseif r ∉ cluster[i]
                    push!(cluster[i],r)
                end
            end
        end
    end
    push!(cluster,filter(x->x∉reduce(vcat,cluster),1:N))
    return reduce(vcat,cluster)
end

sortcormatrix(corr,cluster) = corr[cluster,cluster]
sortcormatrix(corr) = sortcormatrix(corr,findcluster(corr))
sortloadpos(loadpos,cluster) = loadpos[cluster]
time_to_sampleinds(T,times) = searchsortedlast(times,T[1]):searchsortedfirst(times,T[end])

function meanr2!(L,P,d,N,l,Ts)
    L .= 0.0
    @showprogress map(t->addr2!(L,P,d,t,N,l),Ts)
    L .= L ./ length(Ts)
end

function meanr2(d,N,Ts)
    L = zeros(N,N)
    P = zeros(N)
    meanr2!(L,P,d,N,d["samplesize"],Ts)
    return L
end

function addr2!(L,P,d,t,N,l)
    get_Ps!(P,d,t,N,l)
    for n in 1:N
        if iszero(P[n])
            rnn = 0.0
        else
            rnn = D(P[n],P[n],P[n],norm=false)^2 / (P[n]*(1-P[n]))^2
            for m in n+1:N
                if iszero(P[m])
                    rnm = 0.0
                else
                    pnm = sum(view(d["SnapShots"],n,:,t) .* view(d["SnapShots"],m,:,t)) ./ l
                    rnm = D(P[n],P[m],pnm,norm=false)^2 / (P[n]*(1-P[n]) * P[m]*(1-P[m]))
                end
                L[n,m] += rnm
                L[m,n] += rnm
            end
        end
        L[n,n] += rnn
    end
end

function get_Ps!(P,d,t,N,l)
    for n in 1:N
        P[n] = sum(view(d["SnapShots"],n,:,t)) ./ l
    end
end

function D(pn,pm,pnm;norm=true)
    D = pnm-pn*pm
    iszero(D) && return D
    return (norm ? D/D_max(D,pn,pm) : D)
end

#calculates the forward differences between two neighbouring values of v
#and safes the result in r
function ∇!(r,v)
    n=length(v)
    !isequal(length(r),n-1) && error("Error:Dimension mismatch!")
    for (i,x) in enumerate(view(v,1:n-1))
        r[i] = v[i+1]-x
    end
end
function ∇(v)
    r = Vector{Float64}(undef,length(v)-1)
    ∇!(r,v)
    return r
end

#calculates the second order forward differences of v
#and safes the result in r
function Δ!(r,v)
    n=length(v)
    !isequal(length(r),n-2) && error("Error:Dimension mismatch!")
    for (i,x) in enumerate(view(v,1:n-2))
        r[i] = v[i+2]-2*v[i+1]+x
    end
end
function Δ(v)
    r = Vector{Float64}(undef,length(v)-2)
    Δ!(r,v)
    return r
end

#smoothes the input vector `v` by a moving mean of length `2*d` and saves the result
#in `r`
function mm!(r,v,d)
    n = length(v)
    !isequal(n-2*d,length(r)) && error("Error:Dimension Mismatch")
    for i in d+1:n-d
        r[i-d] = mean(view(v,i-d:i+d))
    end
end
function mm(v,d)
    r = Vector{Float64}(undef,length(v)-2*d)
    mm!(r,v,d)
    return r
end

function findinc(v;δ=5*10^3,ε=10)
    #calculate the first order forward difference, smooth is with the moving
    #mean and look for the location of the maximum value
    r = argmax(mm(∇(v),δ))
    m₁ = mean(view(v,1:r))
    m₂ = mean(view(v,r+1:length(v)))
    #if the rise is at the beginning and the means befor and after the rise
    #are close, then it was only the initla rise
    r < 100 && m₂-m₁ < 10.0 && return 0
    return r
end

#Function to calculate the equilibrium size of C₀ from N and μ
h(N,μ) = sqrt(1-exp(-μ/N))
C₀(N,μ) = (1-h(N,μ))^N
μ_from_C(C,N) = -N*log(2*C^(1/N)-C^(2/N))
ML(N,μ) = 2N*h(N,μ)

#finds all indices of extimes that are equal to tend, if non of them is equal
#to tend it returns the index of the maximum value in extimes
function findmaxindex(extimes,tend)
	is = findall(x->x==tend,extimes)
	isempty(is) && return [argmax(extimes)]
	return is
end

