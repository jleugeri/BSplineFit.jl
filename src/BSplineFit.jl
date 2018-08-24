module BSplineFit
export Bspline1D, BSpline1DBasis, fit, decompose

_augment_knots(knots, order) = [fill(knots[1],order);knots;fill(knots[end],order)]

struct BSpline1D{O,KT,CT}
    knots::Vector{KT}
    coeffs::Vector{CT}
end

BSpline1D(order, knots, coeffs=zeros(length(knots)+order)) = 
    BSpline1D{order, eltype(knots), eltype(coeffs)}(_augment_knots(knots, order), coeffs)

function BSpline1DBasis(order, knots, id)
    b = BSpline1D(order, knots)
    b.coeffs[id] = 1
    return b
end

function (s::BSpline1D{O,KT,CT})(x) where {O,KT,CT}
    if s.knots[1] < x < s.knots[end]
        k = searchsortedlast(s.knots, x)-1
        d = s.coeffs[(1+k-O):(1+k)]
        for r ∈ 1:O
            for j ∈ reverse(r:O)
                α = (x-s.knots[1+j+k-O]) / (s.knots[1+j+k+1-r] - s.knots[1+j+k-O])
                d[j+1] = (one(α)-α) * d[j] + α * d[j+1]
            end
        end
        
        return d[O+1]
    else
        return zero(KT)
    end
end

function fit(x, y, order, knots)
    num_basis = length(knots)+order-1
    basis = [BSpline1DBasis(order, knots, id) for id ∈ 1:num_basis]
    M = [b(px) for px ∈ x, b ∈ basis]
    
    c = M\y
    b = BSpline1D(order, knots)
    b.coeffs[1:num_basis] .= c
    return b
end

function decompose(x, y, order, knots; include_component=1:(length(knots)-1+order))
    basis = [BSpline1DBasis(order, knots, id) for id ∈ include_component]
    
    comp = zeros(eltype(y), length(include_component), size(y, 2))
    
    for (i,b) ∈ enumerate(basis)
        bx = b.(x)
        comp[i,:] .= sum(bx .* y, dims=1) ./ sum(bx, dims=1)
    end
    return _augment_knots(knots, order)[order.+include_component], comp
end

end # module
