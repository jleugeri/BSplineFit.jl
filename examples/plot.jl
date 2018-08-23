using BSplineFit, Plots

knots = [0.0,1.0,2.0,4.0,7.0]
data_x = sort(rand(100)*7)
data_y = rand(100)*0.1 .+ 0.5 .+ sin.(data_x)./2

plts = []
for order ∈ 1:3
    p = scatter(data_x, data_y, label="", title="Fit with BSpline of order $(order).")
    for id ∈ 1:length(knots)+order
        b = BSplineFit.BSpline1DBasis(order, knots, id)
        plot!(t->b(t), 0, 7, lw=2, ls=:dash, label="")
    end
    
    bf = fit(data_x, data_y, order, knots)
    plot!(t->bf(t), 0, 7, lc=:black, lw=2, label="", ylims=(-0.5,1.5))
        
    push!(plts, p)
end

p = plot(plts..., layout=grid(3,1), size=(800,800))

savefig(p, "examples/fit.svg")
