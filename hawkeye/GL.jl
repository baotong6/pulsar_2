using Trapz
using Pkg
const logf = Ref{Vector{Float64}}()

function compute_bin(Tlist, m, ω, ϕ) #OK
    n = zeros(Int, m)
    j = @. floor(Int, m*mod2pi(ω*Tlist+ϕ)/(2π)) + 1
    @inbounds for u in 1:m
        n[u] = count(j.==u)
    end
    n
end

function compute_factor(N, m, ν) #OK
    f1 = 1/(2π*ν)
    logf2 = sum(log.(1:N+m-1))-sum(log.(1:m-1))
    logf3 = N*log(m)
    log(f1)-logf2+logf3
end

function precompute_factorial(N::Int) #OK
    logf = zeros(N)
    @inbounds for i in 2:N
        logf[i] = logf[i-1]+log(i)
    end
    logf
end

function get_T_in_mbins(m, ω, ϕ, gtis)
    T = 2π/ω
    dT = T/m
    T_in_perbin = zeros(m)
    for gti in gtis
        t_start = gti[1]
        t_end = gti[2]
        ϕ_start = mod2pi(ω*t_start+ϕ)/(2π)
        pm = floor(Int, m*ϕ_start) + 1
        temp = pm*dT-ϕ_start*T
        T_in_perbin[pm] += temp
        temp += t_start
        pm = mod(pm+1, 1:m)
        while temp < t_end
            if t_end > temp+dT
                T_in_perbin[pm] += dT
            else
                T_in_perbin[pm] += t_end-temp
            end
            pm = mod(pm+1, 1:m)
            temp += dT
        end
    end
    T_in_perbin
end

function compute_S(Tlist, m, ω, ϕ, gtis) #OK
    nj = compute_bin(Tlist, m, ω, ϕ)
    τj = get_T_in_mbins(m, ω, ϕ, gtis)
    sj = τj/(sum(τj)/m)
    lnS = -sum(nj.*log.(sj))
    lnS
end

function compute_W_scaled(Tlist, m, ω, ϕ, factor, gtis) #OK
    n = compute_bin(Tlist, m, ω, ϕ)
    logfac = 0.
    @inbounds for i in 1:m
        if n[i] > 0
            logfac += logf[][n[i]]
        end
    end
    lnS = compute_S(Tlist, m, ω, ϕ, gtis)
    exp(logfac+factor+lnS)
end

function compute_Om1(Tlist, m, ω, factor, ni, gtis) #OK
    ϕ = (0:ni-1)/ni*2π/m
    y = zeros(length(ϕ))
    @inbounds for i in 1:length(y)
        y[i] = compute_W_scaled(Tlist, m, ω, ϕ[i], factor[m], gtis)
    end
    trapz(ϕ, y)*m
end

function compute_Om1w(Tlist, m_max, ω, factor, ni, gtis) #OK
    Om1w = zeros(length(ω), m_max)
    @inbounds Threads.@threads for m in 1:m_max
        Om1x(ω) = compute_Om1(Tlist, m, ω, factor, ni, gtis)
        @views Om1w[:, m] = Om1x.(ω)
    end
    Om1w
end

function compute_GL(Tlist, ω_range=nothing, m_max=12, ni=10; gtis=[]) # gtis only used for computing lnS
    N = length(Tlist)
    global logf[] = precompute_factorial(N)
    ν = m_max - 1
    if ω_range == nothing
        T = maximum(Tlist)-minimum(Tlist)
        ω_hi = min(π*N/T, 2π/10)
        ω_lo = 20π/T
        dω = π/T
        ω = ω_lo:dω:ω_hi
    else
        ω = ω_range
        ω_hi = maximum(ω_range)
        ω_lo = minimum(ω_range)
        dω = ω[2]-ω[1] # ?
    end
    factor = zeros(m_max)
    @inbounds for m in 1:m_max
        factor[m] = compute_factor(N, m, ν)
    end
    if gtis == []
        Om1w = compute_Om1w(Tlist, m_max, ω, factor, ni, [[minimum(Tlist), maximum(Tlist)]])
    else
        Om1w = compute_Om1w(Tlist, m_max, ω, factor, ni, gtis)
    end
    pw = @. 1/ω/log(ω_hi/ω_lo)
    O1m = zeros(m_max)
    @inbounds for i in 1:m_max
        @views O1m[i] = trapz(ω, pw.*Om1w[:,i])
    end
    O_period = @views sum(O1m[2:end])
    p_period = O_period/(1+O_period)
    m_opt = @views argmax(O1m[2:end]) + 1
    S = @views Om1w[:,m_opt]./ω
    S ./= trapz(ω, S)
    ω_peak = ω[argmax(S)]
    cdf = copy(S)
    for i in 1:length(cdf)
        cdf[i] = @views trapz(ω[1:i], S[1:i])
    end
    ω_peak_range = @view ω[0.025 .<cdf.<0.975]
    if length(ω_peak_range) <= 1
        ω_peak_range = [ω_peak-dω, ω_peak+dω]
    else
        ω_peak_range = [ω_peak_range[1], ω_peak_range[end]]
    end
    p_period, ω_peak, m_opt, ω_peak_range
end

