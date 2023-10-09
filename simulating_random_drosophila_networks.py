# -*- coding: utf-8 -*-
"""
Created on Thu Sep 28 16:50:43 2023

@author: eeguia
"""

# Import necessary libraries
import time
start_time = time.time()
import pandas as pd
import networkx as nx
from random import sample
import seaborn as sns
import matplotlib.pyplot as plt

# User declaration of input genes, number of trials to run, and experimental comparison
n_input_genes = 44
n_trials = 10000
exp_comparison = 'SI'

# Import interaction data from BioGRID
print('Loading BioGRID Data')
interaction_data = pd.read_csv(
    r"BIOGRID-ORGANISM-Drosophila_melanogaster-4.4.225.tab3.txt", 
    usecols = ['Entrez Gene Interactor A', 'Entrez Gene Interactor B'], 
    sep = '\t', 
    low_memory = False) 

print('Creating BioGRID Network using NetworkX')
# Make pairs of columns into a list of tuples
interaction_tuple_list = [(row['Entrez Gene Interactor A'], row['Entrez Gene Interactor B']) for index, row in interaction_data.iterrows()]
# Create graph object in NetworkX
biogrid_int_graph = nx.Graph()
biogrid_int_graph.add_edges_from(interaction_tuple_list) 
biogrid_int_graph_nodes = list(biogrid_int_graph.nodes()) 

def first_adj_network(input_genes, remove_inputs_from_neighbor_analysis = True):
    '''
    
    Parameters
    ----------
    input_genes : The list of Entrez Gene IDs from which the function will find the first neighbors.
    remove_inputs_from_neighbor_analysis : if any of the input genes are adjacent to each other, they will be removed from the neighbor directory unless this is set to False.

    Raises
    ------
    ValueError : If the input_genes are not in a list.

    Returns
    -------
    connections : A dictionary with the first neighbors

    '''
    if type(input_genes) != list:
        raise ValueError('Input Genes are not in list format')

    connections = {}    
    #add first neighbors
    first_neighbors = []
    for gene in input_genes:
        for subgene in list(biogrid_int_graph[gene]):
            first_neighbors.append(subgene)
        
    #remove duplicates to only get unique first neighbors
    first_neighbors = list(set(first_neighbors))
    
    if remove_inputs_from_neighbor_analysis == True:
        #removing any input genes from the first neighbor list
        first_neighbors = [x for x in first_neighbors if x not in input_genes]
        
    # go through the first neighbors
    for neighbor in first_neighbors:
        # find the genes it connects to
        adj = list(biogrid_int_graph[neighbor])
        connections_to_input = 0
        for gene in adj:
            # if it connects to an input gene, count that
            if gene in input_genes:
                connections_to_input += 1
        connections[neighbor] = connections_to_input
        
    return connections

print('Running trial')

output_table_2plus = pd.DataFrame(columns = ['Trial', 
                                             'No. of input genes', 
                                             'No. of neighbors', 
                                             'No. of neighbors w/ 1 int', 
                                             'No. of neighbors w/ 2+ ints', 
                                             'Percentage of 2+ ints'])

output_table_3plus = pd.DataFrame(columns = ['Trial', 
                                             'No. of input genes', 
                                             'No. of neighbors', 
                                             'No. of neighbors w/ 1 int', 
                                             'No. of neighbors w/ 2 ints', 
                                             'No. of neighbors w/ 3+ ints', 
                                             'Percentage of 3+ ints'])

for i in range(n_trials):
    trial = i+1
    rand_input_genes = sample(biogrid_int_graph_nodes, n_input_genes)
    connections_dict = first_adj_network(rand_input_genes)
    
    n_neighbors = len(connections_dict)
    n_1_connections = list(connections_dict.values()).count(1)
    n_2_connections = list(connections_dict.values()).count(2)
    n_2plus_connections = n_neighbors - n_1_connections
    n_3plus_connections = n_neighbors - n_1_connections - n_2_connections
    percentage_2plus_connections = n_2plus_connections / n_neighbors
    percentage_3plus_connections = n_3plus_connections / n_neighbors
    
    new_2plus_entry = pd.Series({'Trial': trial, 
                                 'No. of input genes': n_input_genes, 
                                 'No. of neighbors': n_neighbors, 
                                 'No. of neighbors w/ 1 int': n_1_connections, 
                                 'No. of neighbors w/ 2+ ints': n_2plus_connections,
                                 'Percentage of 2+ ints': percentage_2plus_connections})

    new_3plus_entry = pd.Series({'Trial': trial, 
                                'No. of input genes': n_input_genes, 
                                'No. of neighbors': n_neighbors, 
                                'No. of neighbors w/ 1 int': n_1_connections, 
                                'No. of neighbors w/ 2 ints': n_2_connections, 
                                'No. of neighbors w/ 3+ ints': n_3plus_connections, 
                                'Percentage of 3+ ints': percentage_3plus_connections})
    
    new_idx = len(output_table_2plus)
    output_table_2plus.loc[new_idx] = new_2plus_entry
    output_table_3plus.loc[new_idx] = new_3plus_entry
    
    if trial %100 == 0:
        print(f'\t{trial}')

print('Plotting result')
sns.histplot(output_table_2plus['Percentage of 2+ ints'])
plt.title(f'Percentage of 2+ ints (w/ input genes removed)\n{n_trials} trials of a random network with {n_input_genes} input genes')
plt.savefig(f'{exp_comparison} Simulation Results 2+ with input genes removed.png', dpi = 300)
plt.show()

sns.histplot(output_table_3plus['Percentage of 3+ ints'])
plt.title(f'Percentage of 3+ ints (w/ input genes removed)\n{n_trials} trials of a random network with {n_input_genes} input genes')
plt.savefig(f'{exp_comparison} Simulation Results 3+ with input genes removed', dpi = 300)
plt.show()

output_table_2plus.to_excel(f'{exp_comparison} Output Table 2 ints+ with input genes removed.xlsx', index = False)
output_table_3plus.to_excel(f'{exp_comparison} Output Table 3 ints+ with input genes removed.xlsx', index = False)

print("Code completed in %s seconds" % (time.time() - start_time))