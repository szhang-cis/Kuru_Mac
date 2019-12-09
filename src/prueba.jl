function NodeArrangementHex(C)

    linear_bases_idx = [1,(C+2),(C+2)^2,(C+2)^2-(C+1)]
    quad_aranger = vcat(linear_bases_idx,deleteat!(collect(1:(C+2)^2),sort(linear_bases_idx)))
    element_numbering = copy(quad_aranger)
    faces_z = quad_aranger .+ maximum(element_numbering)
    element_numbering = vcat(element_numbering,faces_z)
    println(!(element_numbering.==faces_z[1:4]))
    #element_numbering = element_numbering[!(element_numbering.==faces_z[1:4])]
    return element_numbering
end

numbering = NodeArrangementHex(1)
println(numbering)
