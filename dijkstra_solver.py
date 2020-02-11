from pyqubo import Binary, Constraint
from itertools import combinations
import math
import numpy as np


def nCr(n, r):
    f = math.factorial
    return f(n) / f(r) / f(n-r)


def qubo_c(n):

    if n == 1:
        return 1
    if n == 2:
        return -2
    sum_all = 1 - nCr(n, 1) + 2 * nCr(n, 2) - sum([-qubo_c(k) * nCr(n, k) for k in range(3, n)])
    return sum_all


def dijkstra_solver(G, node1, node2, lagrange):
    H = 0
    pq_vars = {}

    # Constraints
    # Find the complete set of edges, and create PyQUBO variables
    for edge in G.edges():
        edge_label = sorted(edge)[0] + sorted(edge)[1]
        pq_vars[edge_label] = Binary(edge_label)

    for node in G:
        neighbors = G.edges(node)
        incident_edges = [pq_vars[sorted(edge)[0] + sorted(edge)[1]] for edge in neighbors]
        if node == node1 or node == node2:
            # Starting or ending node - 1 in N condition
            H += lagrange * Constraint((sum(incident_edges) - 1) ** 2, label=' special ' + node)
        else:
            # Other node - may be 0 or 2 in N
            # In QUBO, the 0 case does not contribute, assuming the sum is
            # 0 (valid solution)
            sum_z = sum(incident_edges)
            for n in range(2, len(incident_edges) + 1):
                for c in combinations(incident_edges, n):
                    m_prod = np.prod(c)
                    sum_z += qubo_c(n) * m_prod
            H += lagrange * Constraint(sum_z, ' usual ' + node)

    # Objective terms
    for edge in G.edges():
        edge_label = sorted(edge)[0] + sorted(edge)[1]
        H += pq_vars[edge_label] * G[sorted(edge)[0]][sorted(edge)[1]]['weight']

    model = H.compile()
    return model
