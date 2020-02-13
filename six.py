from dijkstra_solver import dijkstra_solver
import networkx as nx
from dimod import BinaryQuadraticModel
import dimod
import dwave.inspector
from dwave.system import DWaveSampler, EmbeddingComposite


G = nx.Graph()
G.add_weighted_edges_from([('S', 'a', 4.0), ('a', 'b', 7.0), ('b', 'G', 2.0), ('G', 'd', 3.0), ('c', 'd', 8.0), ('c', 'S', 6.0), ('c', 'b', 9.0)])
lagrange = 20
chainstrength = 20
numruns = 1000

model = dijkstra_solver(G, 'S', 'G', lagrange)
Q = model.to_qubo()
bqm = BinaryQuadraticModel.from_qubo(Q[0], offset=Q[1])
sampler = EmbeddingComposite(DWaveSampler())
sampleset = sampler.sample(bqm, chain_strength=chainstrength, num_reads=numruns)
dwave.inspector.show(bqm, sampleset, sampler)
