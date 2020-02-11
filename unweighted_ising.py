from dwave.system import DWaveSampler
from dwave.system.composites import FixedEmbeddingComposite
from minorminer import find_embedding
import dimod
import math
import dwave.inspector


def nCr(n, r):
    f = math.factorial
    return f(n) / f(r) / f(n-r)


N = 16
chainstrength = 40
numruns = 1000

# Ising model is 0.5 x_i x_j - 0.5
# There is a constant offset of -1/2 per pair in the clique
# The number of pairs is (N 2)
h = {}
J = {(i, j): 0.5 for i in range(N) for j in range(i + 1, N)}

offset = -nCr(N, 2)/2
bqm = dimod.BinaryQuadraticModel.from_ising(h, J, offset=offset)

# Do the embedding
dwave_sampler = DWaveSampler()
A = dwave_sampler.edgelist
embedding = find_embedding(J, A)

sampler = FixedEmbeddingComposite(DWaveSampler(), embedding)
sampleset = sampler.sample(bqm, chain_strength=chainstrength, num_reads=numruns)
dwave.inspector.show(bqm, sampleset, sampler)
