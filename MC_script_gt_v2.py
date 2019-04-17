#!/usr/bin/env python3
"""
the script takes in input a file path and the number of randomization to run:
> python3 MC_script.py input_file n
- input_file: str, the path to a .graphml file
- n: int, number of randomized graphs to run
"""
import os, shutil
import argparse
import gc
import datetime
import numpy.random

import networkx as nx
import numpy as np
import subprocess as sp
import graph_tool as gt

from graph_tool.generation import random_rewire


class RandomizedMotifs():

    def __init__(self, inp_path, out_dir, seed, reps, motif_notation):

        if os.path.isfile(inp_path):
            self.in_path = inp_path
        else:
            raise ValueError('Input file does not exist.')

        if os.path.isdir(out_dir):
            self.out_dir = out_dir
        else:
            raise ValueError('Output folder does not exist.')

        if None == seed:
            self.SEED = seed
        else:
            self.SEED = (int(seed))

        if int(reps)>0:
            self.repetitions = int(reps)
        else:
            raise ValueError('Number of repetitions must be an integer >= 1.')


        #load motif number to adj matrix dict
        self.d_motif_adj = {}
        with open(motif_notation, 'r') as fh:
            for line in fh:
                adj, motif = line.split('\t')
                self.d_motif_adj[int(motif)] = adj


        #load input graph (infer type of input)
        fmt = self.in_path.split('.')[-1]

        if 'graphml' == fmt:
            with open(self.in_path, 'r') as fh:
                self.G = nx.read_graphml(fh, node_type=int)
        elif 'csv' == fmt:
            with open(self.in_path, 'rb') as fh:
                self.G = nx.read_edgelist(fh, delimiter=',', create_using=nx.DiGraph())
        else:
            raise ValueError('Unknown input file format.')


        self.N = self.G.number_of_nodes()

        if 0 != self.N:

            # relabel nodes
            d_relabel = dict(zip(self.G.nodes(), range(self.N)))
            self.G = nx.relabel_nodes(self.G, d_relabel)

            #convert to gt graph
            edges_array = np.array(list(self.G.edges(data=False)))

            self.G = gt.Graph()
            self.G.add_edge_list(edges_array)

            del edges_array
            gc.collect()

            #create dir for tmp files
            self.tmp_dir = self.out_dir+'/tmp_'+str(self.SEED)
            os.mkdir(self.tmp_dir)


    def __parse_adj(self):
        
        next_id = 3
        adj = ''
        d = {}

        with open(self.tmp_dir+'/adjMatrix.txt', 'r') as fh:

            for i, line in enumerate(fh):

                if i == next_id:
                    d[adj] = line.split(':')[1].strip()
                    next_id += 5
                    adj = ''
                else:
                    adj = adj+line.strip()

        return d


    def __read_counts(self):

        d_motif_read = {}
        d_motif_count = {}

        with open(self.tmp_dir+'/Motif_count.txt', 'r') as fh:
            fh.readline()
            for line in fh:
                ID, count = line.strip().split('\t')
                d_motif_read[ID] = float(count)

        if len(d_motif_read.keys()) > 0:
            d_adj_ID = self.__parse_adj()

            for i in range(1,14):
                try:
                    d_motif_count[i] = int(d_motif_read[d_adj_ID[self.d_motif_adj[i]]])
                except KeyError:
                    d_motif_count[i] = 0

        else:
            for i in range(1,14):
                d_motif_count[i] = 0

        return d_motif_count


    def save_empty(self, it=''):

        lines = [str(i)+'\t'+'0\n' for i in range(1,14)]

        with open(self.out_dir+'/'+self.in_path.rsplit('/',1)[-1].split('.')[0]+'_MC_'+str(self.SEED)+'_'+str(it)+'.tsv', 'w') as fh:                
            fh.writelines(lines)


    def randomize_and_dump(self):
        
        # randomization
        random_rewire(self.G, n_iter=4)

        # DUMPING
        edges = self.G.get_edges()[:, :2]
        with open(self.tmp_dir+'/tmp_kavosh.tsv', 'w+') as fh:
            np.savetxt(fh, edges, fmt='%s', delimiter='\t', header=str(self.N), comments='')

        #free some memory
        del edges
        gc.collect()

        return

    
    def run_kavosh(self):

        os.chdir('Kavosh_src')

        Kavosh_args = ['-i', self.tmp_dir+'/tmp_kavosh.tsv', '-o', self.tmp_dir, '-s', '3']
        sp.check_call(['./Kavosh']+Kavosh_args, stdout=sp.DEVNULL)

        os.chdir('..')

        return


    def get_bowtie_components(self):
        '''Classifying the nodes of a network into a bow-tie structure.
        Here we follow the definition of a bow-tie as in: 
        "Bow-tie Decomposition in Directed Graphs" - Yang et al. IEEE (2011) 
        
        input:  NetworkX directed graph or numpy adjacency matrix
        output: sets of nodes in the specified partitions (following the 
            NetworkX input graph node labelling or labelled according to
            the order of the adjacency matrix [0, n-1])
        '''
        # input G from self.tmp_dir + '/tmp_kavosh.tsv'
        G = nx.read_edgelist(self.tmp_dir + '/tmp_kavosh.tsv', delimiter = '\t', create_using = nx.DiGraph())
        # nx.reverse: reverse the direction of links
        GT = nx.reverse(G, copy=True) 
    
        strongly_con_comp = list(nx.strongly_connected_components(G))    
        strongly_con_comp = max(strongly_con_comp, key=len)

        S = strongly_con_comp

        v_any = list(S)[0]
        DFS_G = set(nx.dfs_tree(G, v_any).nodes())
        DFS_GT = set(nx.dfs_tree(GT, v_any).nodes())
        OUT = DFS_G - S
        IN = DFS_GT - S
        V_rest = set(G.nodes()) - S - OUT - IN

        TUBES = set()
        INTENDRILS = set()
        OUTTENDRILS = set()
        OTHER = set()
        for v in V_rest:
            irv = len(IN & set(nx.dfs_tree(GT,v).nodes())) is not 0
            vro = len(OUT & set(nx.dfs_tree(G,v).nodes())) is not 0
            if irv and vro:
                TUBES.add(v)
            elif irv and not vro:
                INTENDRILS.add(v)
            elif not irv and vro:
                OUTTENDRILS.add(v)
            elif not irv and not vro:
                OTHER.add(v)
            
        d = {'S': len(S),
             'IN': len(IN),
             'OUT': len(OUT),
             'TUBES': len(TUBES),
             'INTENDRILS': len(INTENDRILS),
             'OUTTENDRILS': len(OUTTENDRILS),
             'OTHER': len(OTHER)
            }
        return d


    def parse_and_save(self, it=''):

        #PARSE KAVOSH OUT
        d_motif_count = self.__read_counts()
        
        # PARSE BOW-TIE
        d_bowtie = self.get_bowtie_components()

        ### SAVE OUTPUT
        lines = []
        for m in sorted(d_motif_count.keys()):
            lines.append(str(m)+'\t'+str(d_motif_count[m])+'\n')
        for m in sorted(d_bowtie.keys()):
            lines.append(str(m) + '\t' + str(d_bowtie[m]) + '\n')

        with open(self.out_dir+'/'+self.in_path.rsplit('/',1)[-1].split('.')[0]+'_MC_'+str(self.SEED)+'_'+str(it)+'.tsv', 'w') as fh:                
            fh.writelines(lines)

        return


    def clean(self, dir=False):

        ### DELETE TMP FILES
        if dir:
            shutil.rmtree(self.tmp_dir)
        else:
            rm = ['/tmp_kavosh.tsv', '/Motif_count.txt', '/adjMatrix.txt']

            for f in rm:
                os.remove(self.tmp_dir+f)

        return


    def randomize_and_count(self):

        numpy.random.seed(self.SEED)

        for i in range(self.repetitions):

            if 0 == self.N:
                self.save_empty(i+1)
            else:
                start = datetime.datetime.now()

                self.randomize_and_dump()
                print('Rand time: {}'.format(datetime.datetime.now()-start))

                self.run_kavosh()
                self.parse_and_save(i+1)
                self.clean()

        if 0 != self.N:
            #free some memory and clean all tmp
            self.clean(dir=True)

            del self.G
            gc.collect()

        return




if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('-i','--input_file', help='Path to the .graphml file containing the network to be randomized')
    parser.add_argument('-o', '--output_folder', default='/home/ubuntu/Output', help='Path to the folder where the file with motif counts will be stored')
    parser.add_argument('-s', '--SEED', default=None, help='Random seed (needed for network randomization)')
    parser.add_argument('-n', '--repetitions', default=1, help='Number of randomizations to perform')
    parser.add_argument('-m', '--motif_notation', default="3Motif_notation.tsv", help='Motif Notation .tsv file')

    args = parser.parse_args()

    assert os.path.isfile(args.input_file), "Input file {0} not found.".format(args.input_file)
    assert os.path.isfile(args.motif_notation), "Motif file {0} not found.".format(args.motif_notation)

    input_file = os.path.abspath(args.input_file)
    output_folder = os.path.abspath(args.output_folder)
    motif_notation = os.path.abspat(args.motif_notation)

    simul = RandomizedMotifs(input_file, output_folder, args.SEED, args.repetitions, motif_notation)

    simul.randomize_and_count()
