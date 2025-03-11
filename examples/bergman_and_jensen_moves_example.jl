M = uniform_matroid(2,4)

p1 = point(0,0,0,1)
p2 = point(0,0,0,0)

p3 = point(0,0,0,0)
p4 = point(0,1,1,0)
p5 = point(0,2,0,0)
p6 = point(0,0,2,0)
p7 = point(0,2,2,0)

f1 = support([p1,p2],[0,0])
f2 = support([p3,p4,p5,p6,p7],[0,-6,-6,-4,-10])


f2Target = support([p3,p4,p5,p6,p7],[0,0,6,-4,2])

mixedSupport = mixed_support((f1,f2))
targetSupport = mixed_support((f1,f2Target))

f2Active = support([p3, p6], [0, -4])
chainOfFlats = chain_of_flats(M, [[3]])
candidateOne = mixed_cell(mixed_support((f1, f2Active)), chainOfFlats)

f2Active = support([p6, p7], [-4, -10])
chainOfFlats = chain_of_flats(M, [[4]])
candidateTwo = mixed_cell(mixed_support((f1, f2Active)), chainOfFlats)

T = tracker(mixedSupport, [candidateOne, candidateTwo], [targetSupport])