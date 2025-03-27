M = matroid(Oscar.matrix(QQ, [1 1 0 0 0 1 1 0 0 0; -1 0 -1 -1 0 -1 0 -1 -1 0; 0 -1 1 0 -1 0 -1 1 0 -1; 0 0 0 1 1 0 0 0 1 1]))
targetSupports = Support[]
for i in [1,2,3,4,5]
        pi = point([n in [i,i+5] ? 1 : 0 for n in 1:10])
        p0 = point([0 for n in 1:10])
        push!(targetSupports, support([pi,p0],[0,0]))
end

targetSupport = mixed_support(targetSupports)

for S in supports(targetSupport)
    println(S)
end