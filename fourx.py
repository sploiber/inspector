from dijkstra_solver import dijkstra_solver
import networkx as nx
from dimod import BinaryQuadraticModel
import dimod
import dwave.inspector
from dwave.system import DWaveSampler, EmbeddingComposite


G = nx.Graph()
G.add_weighted_edges_from([('S', 'a', 5.0), ('a', 'G', 6.0), ('G', 'c', 3.0), ('c', 'S', 17.0), ('a', 'c', 8)])
lagrange = 4
chainstrength = 20
numruns = 1000

model = dijkstra_solver(G, 'S', 'G', lagrange)
Q = model.to_qubo()
bqm = BinaryQuadraticModel.from_qubo(Q[0], offset=Q[1])
sampler = EmbeddingComposite(DWaveSampler())
sampleset = sampler.sample(bqm, chain_strength=chainstrength, num_reads=numruns)
dwave.inspector.show(bqm, sampleset, sampler)
