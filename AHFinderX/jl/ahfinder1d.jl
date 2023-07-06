module AHFinder1d

const M = 1.0

function metric(r)
    Ψ = (1 + M / 2r)^4
    dr_Ψ = 4 * (1 + M / 2r)^3 * M / 2 * (-1 / r^2)
    grr = Ψ
    dr_grr = dr_Ψ
    gθθ = Ψ * r^2
    dr_gθθ = dr_Ψ * r^2 + Ψ * 2r
    return grr, gθθ, dr_grr, dr_gθθ
end

function expansion(h)
    r = h

    grr, gθθ, dr_grr, dr_gθθ = metric(r)
    gurr = 1 / grr
    guθθ = 1 / gθθ
    dr_gurr = -1 / grr^2 * dr_grr
    dr_guθθ = -1 / gθθ^2 * dr_gθθ
    Gammarrr = 1 / 2 * gurr * (dr_grr + dr_grr - dr_grr)
    Gammarθθ = 1 / 2 * gurr * (-dr_gθθ)

    gradFr = 1
    dr_gradFr = 0
    gradFur = gurr * gradFr
    dr_gradFur = dr_gurr * gradFr + gurr * dr_gradFr
    len_gradF = sqrt(gradFur * gradFr)
    dr_len_gradF = 1 / (2 * sqrt(gradFur * gradFr)) * (dr_gradFur * gradFr + gradFur * dr_gradFr)
    sr = gradFr / len_gradF
    dr_sr = dr_gradFr / len_gradF - gradFr / len_gradF^2 * dr_len_gradF
    sur = gurr * sr
    grad_srr = dr_sr - Gammarrr * sr
    grad_sθθ = -Gammarθθ * sr

    ρ = 2 * r^2 * len_gradF / ((gurr - sur * sur) * (1 - 1 * 1) + 2 * (guθθ - 0 * 0) * (r^2 - 0 * 0))

    qurr = gurr - sur * sur
    quθθ = guθθ
    Θ = qurr * grad_srr + 2 * quθθ * grad_sθθ

    return ρ, Θ
end

function step(h)
    # α = 1.0
    # β = 0.5
    # lmax = 1
    # A = α / (lmax * (lmax + 1)) + β
    # B = β / α

    A = 1.0
    B = 0.0

    ρ, Θ = expansion(h)

    l = 0
    hnew = h - A / (1 + B * l * (l + 1)) * (ρ * Θ)

    hnew = clamp(hnew, 0.9 * h, 1.1 * h)

    return hnew
end

function solve(h)
    n = 0
    while true
        ρ, Θ = expansion(h)
        println("$n: $h: $Θ ($ρ)")
        abs(Θ) ≤ sqrt(eps()) && break
        n ≥ 100 && break
        n += 1
        h = step(h)
    end
    return h
end

end
