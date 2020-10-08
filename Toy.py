import math
import networkx as nx
import pandas as pd
import matplotlib.pyplot as plt
import random
import queue
import numpy as np
from operator import itemgetter

G = nx.Graph()
G.add_nodes_from(['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M'])
G.add_edge('A', 'B', weight=1)
G.add_edge('A', 'C', weight=9)
G.add_edge('A', 'G', weight=8)
G.add_edge('A', 'H', weight=18)
G.add_edge('G', 'F', weight=4)
G.add_edge('F', 'C', weight=2)
G.add_edge('B', 'I', weight=5)
G.add_edge('I', 'H', weight=3)
G.add_edge('B', 'E', weight=7)
G.add_edge('B', 'D', weight=5)
G.add_edge('B', 'J', weight=1)
G.add_edge('B', 'K', weight=1)
G.add_edge('E', 'D', weight=2)
G.add_edge('H', 'L', weight=15)
G.add_edge('H', 'M', weight=3)


def degree(G):
    degree = nx.degree_centrality(G)
    sw={}
    wd={}
    for node in G.nodes():
        sw[node]=0
        for neigh in G.neighbors(node):
            sw[node] += G[neigh][node]['weight']
        wd[node] = (degree[node]**(1/2))*sw[node]
    return wd


def betweeness(G):
    between = nx.betweenness_centrality(G, weight='weight')
    return between


def closeness_centrality(G):
    close = nx.closeness_centrality(G, distance='weight')
    return close


def WH(G):
    Wh = {}
    WH = {}
    # for i in range(df.shape[0] - 1):
    #   G[df['Node1'][i]][df['Node2'][i]]['weight'] = G.degree[df['Node1'][i]]*G.degree[df['Node2'][i]]

    for node in G.nodes():
        Wh[node] = 0
        n_deg = 0
        dic1 = {}
        dic2 = {}
        for j in G.neighbors(node):
            n_deg += G.degree[j]
            dic1[j] = G[node][j]['weight']
            dic2[j] = G.degree[j]
        x = 0
        for m in range(n_deg):
            x = n_deg - m
            flag = 0

            if x > min(dic1.values()):
                a = list(dic1.keys())[list(dic1.values()).index(min(dic1.values()))]
                if dic2[a] == 1:
                    del dic1[a]
                    del dic2[a]
                else:
                    dic2[a] -= 1
            else:
                Wh[node] = x
                break

        for node in G.nodes():
            x = 0
            for neigh in G.neighbors(node):
                if neigh in Wh.keys():
                    x += Wh[neigh] * G.degree[neigh]
            WH[node] = x
    return WH


def WvoteRank(G, n):
    capability = {}
    Nodes = list(G.nodes())
    totdeg = 0
    for node in Nodes:
        capability[node] = [0, 1]
        totdeg += G.degree(node)
    totdeg /= len(Nodes)
    val = 1 / totdeg
    selected = []
    for i in range(n):
        for node in Nodes:
            capability[node][0] = 0
        bes = 0
        bestnode = -1
        for node in Nodes:
            allnodes = list(G.neighbors(node))
            for neig in allnodes:
                if capability[neig][1] > 0:
                    capability[node][0] += capability[neig][1] * G[neig][node]['weight']
            capability[node][0] = int(math.sqrt(capability[node][0] * G.degree(node)))
            if capability[node][0] > bes and node not in selected:
                bes = capability[node][0]
                bestnode = node
        if bestnode == -1:
            break
        selected.append(bestnode)
        neighbors = list(G.neighbors(bestnode))
        for node in neighbors:
            capability[node][1] -= val
        caps = {}
        for i in G.nodes():
            caps[i] = capability[i][0]
    return selected, caps


def SIRModel(g, spreaders, beta, S):
    infected_scale = np.array(np.zeros(50))
    ftc = 0
    while S > 0:
        S = S - 1
        infected = spreaders
        status = {}
        for i in g.nodes():
            status[i] = 0
        for i in infected:
            status[i] = 1
        n = g.number_of_nodes()
        infected_nodes = len(infected)
        recovered_nodes = 0
        time_stamp = 0
        infected_scale[time_stamp] = infected_scale[time_stamp] + (infected_nodes + recovered_nodes) / n
        infected = spreaders
        while len(infected) > 0:
            susceptible_to_infected = []
            time_stamp = time_stamp + 1
            for i in infected:
                susceptible = []
                status[i] = 2
                for neighbor in g.neighbors(i):
                    if status[neighbor] == 0:
                        susceptible.append(neighbor)
                total_susceptible = len(susceptible)
                no_of_susceptible_to_infected = round(beta * total_susceptible)
                while no_of_susceptible_to_infected > 0:
                    random_index = random.randint(0, total_susceptible - 1)
                    if susceptible[random_index] not in susceptible_to_infected:
                        susceptible_to_infected.append(susceptible[random_index])
                        status[susceptible[random_index]] = 1
                        no_of_susceptible_to_infected = no_of_susceptible_to_infected - 1
            infected_nodes = len(susceptible_to_infected)
            recovered_nodes = len(infected)
            ftc = ftc + recovered_nodes / n
            infected_scale[time_stamp] = infected_scale[time_stamp] + (infected_nodes + recovered_nodes) / n
            infected = susceptible_to_infected
    return infected_scale, ftc


def Find_Spreaders(g, centrality, no_of_spreaders):
    spreaders_list = centrality
    # print(spreaders_list)
    # print(no_of_spreaders)
    spreaders = []
    while no_of_spreaders > 0:
        no_of_spreaders = no_of_spreaders - 1
        matches = sorted(spreaders_list.items(), key=lambda kv: kv[1], reverse=True)
        # print( matches)
        # print( spreaders)
        key = (matches[0][0])
        spreaders.append(matches[0][0])
        del spreaders_list[key]
    # print(spreaders)
    return spreaders


def E_W_Voterank(G, n):
    capability={}
    Nodes = list(G.nodes())
    totdeg = 0
    for node in Nodes:
        capability[node] = [0, 1]
        totdeg += G.degree(node)
    totdeg /= len(Nodes)
    val = 1 / totdeg
    selected = []
    for i in range(n):
        for node in Nodes:
            capability[node][0] = 0
        bes = 0
        bestnode = -1
        for node in Nodes:
            allnodes = list(G.neighbors(node))
            for neig in allnodes:
                if capability[neig][1] > 0:
                    capability[node][0] += capability[neig][1] * G[neig][node]['weight']
                    for n2p in list(G.neighbors(neig)):
                        if capability[n2p][1] > 0:
                            capability[node][0] += capability[n2p][1] * G[n2p][neig]['weight']
            capability[node][0] = int(math.sqrt(capability[node][0] * G.degree(node)))
            if capability[node][0] > bes and node not in selected:
                bes = capability[node][0]
                bestnode = node
        if bestnode == -1:
            break
        selected.append(bestnode)
        neighbors = list(G.neighbors(bestnode))
        for node in neighbors:
            capability[node][1] -= val
            for neig in list(G.neighbors(node)):
                capability[neig][1] -= val/10
        caps = {}
        for i in G.nodes():
            caps[i] = capability[i]
    return selected, caps



def dict2list(d, k):
    selected = []
    i = 0
    for node in sorted(d.items(), key=itemgetter(1), reverse=True):
        if i < k:
            selected.append(node[0])
            i += 1
    return selected


def SPL(G, spreaders):
    total = 0
    n = len(spreaders)
    c = 0
    for i in range(0, n):
        node_i = spreaders[i]
        for j in range(i + 1, n):
            node_j = spreaders[j]
            if nx.has_path(G, spreaders[i], spreaders[j]):
                c = c + 1

                total = total + nx.shortest_path_length(G, node_i, node_j, 'weight')
    # print(total)
    # print(n)
    total = total / (n * n - n)
    total = 2 * total
    return total





nodes = list(G.nodes())
edges = list(G.edges())

k = 1
#S = int(input("Input the number of times SIR should average out: "))
#T = int(input("Input the timestamps: "))
beta = 0.01
deg = degree(G)
close = closeness_centrality(G)
between = betweeness(G)
W_Hindex = WH(G)
W_Hindex_selected = dict2list(W_Hindex, k)
WVoterank_selected, caps1 = WvoteRank(G, k)
E_WVoterank_selected, caps = E_W_Voterank(G, k)

print("Degree")
print(deg)

print("close")
print(close)

print("between")
print(between)

print("WHI")
print(W_Hindex)

print("WVP")
print(caps)

print("WV")
print(caps1)