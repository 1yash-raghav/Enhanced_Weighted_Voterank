import math
import networkx as nx
import pandas as pd
import matplotlib.pyplot as plt
import random
import queue
import numpy as np
from operator import itemgetter
import csv

df = pd.read_csv('New.csv', squeeze=True)
G = nx.Graph()
nod = []
for i in range(df.shape[0] - 1):
    nod.append(df['Node1'][i])
    nod.append(df['Node2'][i])
tmp = set(nod)
nod = list(tmp)
G.add_nodes_from(nod)
print(len(G.nodes))
for i in range(df.shape[0] - 1):
    G.add_edge(df['Node1'][i], df['Node2'][i], weight=df['Weight'][i]+15)
G.remove_edges_from(nx.selfloop_edges(G))


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
        capability[bestnode][1] = 0
        capability[bestnode][0] = 0
        neighbors = list(G.neighbors(bestnode))
        for node in neighbors:
            capability[node][1] -= val
    return selected


def SIRModel(g, spreaders, beta, S):
    infected_scale = np.array(np.zeros(T))
    ftc = np.array(np.zeros(T))
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
            ftc[time_stamp] = ftc[time_stamp] + recovered_nodes / n
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


def E_W_Voterank(G, n, core):
    capability={}
    Nodes = list(G.nodes())
    totdeg = 0.0
    for node in Nodes:
        capability[node] = [0.0, 1.0]
        totdeg += G.degree(node)
    totdeg /= len(Nodes)
    val = 1.0 / totdeg
    selected = []
    for i in range(n):
        for node in Nodes:
            capability[node][0] = 0.0
        bes = 0.0
        bestnode = -1
        for node in Nodes:
            allnodes = list(G.neighbors(node))
            for neig in allnodes:
                if capability[neig][1] > 0.0:
                    capability[node][0] += core[neig] * capability[neig][1] * G[neig][node]['weight']
                    for n2p in list(G.neighbors(neig)):
                        if capability[n2p][1] > 0.0:
                            capability[node][0] +=core[n2p] * capability[n2p][1] * G[n2p][neig]['weight']
            capability[node][0] = int(math.sqrt(capability[node][0] * core[node] * G.degree(node)))
            if capability[node][0] > bes and node not in selected:
                bes = capability[node][0]
                bestnode = node
        if bestnode == -1:
            break
        selected.append(bestnode)
        capability[bestnode][1] = 0.0
        capability[bestnode][0] = 0.0
        neighbors = list(G.neighbors(bestnode))
        for node in neighbors:
            capability[node][1] -= val
            for neig in list(G.neighbors(node)):
                capability[neig][1] -= val
    return selected



def dict2list(d, k):
    selected = []
    i = 0
    for node in sorted(d.items(), key=itemgetter(1), reverse=True):
        if i < k:
            selected.append(node[0])
            i += 1
    return selected


def ftToftc(ft):
    ftc = [0 for i in range(0, T + 1)]
    for i in range(T):
        ftc[i + 1] = ft[i] + ftc[i]
    return ftc


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


def SIR(G, selected, beta, ft, ftc):
    Q = queue.Queue(maxsize=len(list(G.nodes())))
    recovered = []
    infected = []
    for node in selected:
        Q.put(node)
        infected.append(node)
    t = 0
    while t < T:
        tot = []
        while not Q.empty():
            curr = Q.get()
            if curr in recovered:
                continue
            recovered.append(curr)
            ne = list(G.neighbors(curr))
            random.shuffle(ne)
            neighbors = int(len(ne) * beta) + 1
            count = 0
            idx = 0
            while count < neighbors and idx < len(ne):
                if ne[idx] not in recovered and ne[idx] not in tot:
                    tot.append(ne[idx])
                    count += 1
                idx += 1
        for i in tot:
            Q.put(i)
        tot.clear()
        ftc[t] = ftc[t]+ len(recovered)
        ft[t] = ft[t]+len(recovered)+Q.qsize()
        t += 1
    return ft, ftc


def SIRcall(G, selected, beta, S):
    ft = [0.0 for i in range(T)]
    ftc = [0.0 for i in range(T)]
    for i in range(S):
        ft, ftc = SIR(G, selected, beta, ft, ftc)
    ft[:] = [x/S for x in ft]
    ftc[:] = [x/S for x in ftc]
    return ft, ftc

edges = list(G.edges())
nodes = list(G.nodes())

k = int(input("Input the number of infectors: "))
S = int(input("Input the number of times SIR should average out: "))
T = int(input("Input the timestamps: "))
beta = float(input("Input beta: "))

core = nx.core_number(G)
deg = degree(G)
print("Degree Done \n")
close = closeness_centrality(G)
print("Closeness Done \n")
between = betweeness(G)
print("Betweeness Done\n")
W_Hindex = WH(G)
print("W_Hindex Done\n")

W_Hindex_selected = dict2list(W_Hindex, k)
WVoterank_selected = WvoteRank(G, k)
print("WVoteRank Done \n")
E_WVoterank_selected = E_W_Voterank(G, k, core)
print("E_WVoteRank Done \n")
deg_sel = dict2list(deg, k)
close_sel = dict2list(close, k)
between_sel = dict2list(between, k)

ft_deg, ftc_deg = SIRcall(G, deg_sel, beta, S)
print("Deg SIR Done")
ft_bet, ftc_bet = SIRcall(G, between_sel, beta, S)
print("Bet SIR Done \n")
ft_close, ftc_close = SIRcall(G, close_sel, beta, S)
print("Clos SIR Done \n")
ft_whindex, ftc_whindex = SIRcall(G, W_Hindex_selected, beta, S)
print("WHI SIR Done \n")
ft_wvote, ftc_wvote = SIRcall(G, WVoterank_selected, beta, S)
print("WVo SIR Done \n")
ft_e_wvote, ftc_e_wvote = SIRcall(G, E_WVoterank_selected, beta, S)
print("EWV SIR Done \n")

fileName= 'New'
with open(fileName+'ft.csv', 'w', newline='') as f:
    thewriter=csv.writer(f)
    thewriter.writerow(['Degree', 'Between', 'Closeness', 'WHI', 'WV', 'EWV'])
    for i in range(T):
        thewriter.writerow([str(ft_deg[i]), str(ft_bet[i]), str(ft_close[i]), str(ft_whindex[i]), str(ft_wvote[i]), str(ft_e_wvote[i])])

with open(fileName+'ftc.csv', 'w', newline='') as f:
    thewriter=csv.writer(f)
    thewriter.writerow(['Degree', 'Between', 'Closeness', 'WHI', 'WV', 'EWV'])
    for i in range(T):
        thewriter.writerow([str(ftc_deg[i]), str(ftc_bet[i]), str(ftc_close[i]), str(ftc_whindex[i]), str(ftc_wvote[i]), str(ftc_e_wvote[i])])

#F(t) vs  t
plt.plot(range(T), ft_deg, '-g*', label="Degree ", markersize=5)
plt.plot(range(T), ft_close, '-r*', label="Closeness", markersize=5)
plt.plot(range(T), ft_bet, '-b*', label="Betweeness", markersize=5)
plt.plot(range(T), ft_whindex, '-yp', label="W_Hindex", markersize=5)
plt.plot(range(T), ft_wvote, '-mo', label="WVoteRank", markersize=5)
plt.plot(range(T), ft_e_wvote, label="Improved-WVoterank", marker='X', markerfacecolor='black')
plt.xlabel('Time')
plt.ylabel('F(t)')
plt.legend(loc='best')
n = len(list(G.nodes()))
plt.show()


#F(tc) vs t
plt.plot(range(T), ftc_deg, '-g*', label="Degree ", markersize=5)
plt.plot(range(T), ftc_close, '-r*', label="Closeness", markersize=5)
plt.plot(range(T), ftc_bet, '-b*', label="Betweeness", markersize=5)
plt.plot(range(T), ftc_whindex, '-yp', label="W_Hindex", markersize=5)
plt.plot(range(T), ftc_wvote, '-mo', label="WVoteRank", markersize=5)
plt.plot(range(T), ftc_e_wvote, label="Improved-WVoterank", marker='X', markerfacecolor='black')
plt.xlabel('Time')
plt.ylabel('F(tc)')
plt.legend(loc='best')
n = len(list(G.nodes()))
plt.show()


