import dimod
from dimod import BinaryQuadraticModel
from dwave.system import DWaveSampler, EmbeddingComposite
import dwave.inspector
import numpy as np

with open("willxujun.data", "r") as myfile:
    data = myfile.read()
z2 = [int(k) for k in data.split()]
arr = np.array(z2)
numpy_matrix = arr.reshape(64, 64)

bqm = dimod.BinaryQuadraticModel.from_numpy_matrix(numpy_matrix)
sampler = EmbeddingComposite(DWaveSampler())
sampleset = sampler.sample(bqm, chain_strength=10000, num_reads=1000)
dwave.inspector.show(bqm, sampleset, sampler)
