points = unpack_solution.(sols.u)
# reshape into the same dimensions as the image
points = reshape(points, (100, 100))

times = map(points) do gp
    # check if went off the integration chart on the inner boundary
    if gp.status == StatusCodes.WithinInnerBoundary
        # get the time coordinate
        gp.x[1]
    else
        NaN
    end
end

heatmap(α, β, times, aspect_ratio = 1)