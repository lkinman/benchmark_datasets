import numpy as np
import pandas as pd
from datetime import datetime as dt
import argparse
import os
import re

def add_args(parser):
    parser.add_argument('--infile', type = str, required=True, help='input .star file for filtering')
    parser.add_argument('--two', type = int, required=True, help = 'number of 2bp particles to be randomly subsampled')
    parser.add_argument('--eight', type = int, required=True, help = 'number of 8 bp particles to be randomly subsampled')
    parser.add_argument('--fourteen', type = int, default = None, help = 'number of 14 bp particles to be randomly subsampled')
    parser.add_argument('--outfile', default = './', help = 'output .star file')
    return parser 

class Starfile():
    
    def __init__(self, headers, df, data_optics=None, relion31=False):
        assert headers == list(df.columns), f'{headers} != {df.columns}'
        self.headers = headers
        self.df = df
        self.data_optics = data_optics
        self.relion31 = relion31

    def __len__(self):
        return len(self.df)

    @classmethod
    def load(cls, starfile):
        # detect star file type
        f = open(starfile,'r')
        BLOCK = 'data_'
        while 1:
            for line in f:
                if line.startswith(BLOCK):
                    break
            break
        if line.startswith('data_optics'):
            return cls._parse_relion31(starfile)
        else:
            return cls._parse_block(starfile, block_header='data_')
  
    @classmethod
    def _parse_relion31(cls, starfile):
        data_optics = cls._parse_block(starfile, block_header='data_optics')
        s = cls._parse_block(starfile, block_header='data_particles')
        s.data_optics = data_optics
        s.relion31 = True
        return s

    @classmethod
    def _parse_block(self, starfile, block_header='data_'):
        f = open(starfile,'r')
        # get to data block
        while 1:
            for line in f:
                if line.startswith(block_header):
                    break
            break
        # get to header loop
        while 1:
            for line in f:
                if line.startswith('loop_'):
                    break
            break
        # get list of column headers
        while 1:
            headers = []
            for line in f:
                if line.startswith('_'):
                    headers.append(line)
                else:
                    break
            break 
        # assume all subsequent lines until empty line is the body
        headers = [h.strip().split()[0] for h in headers]
        body = [line]
        for line in f:
            if line.strip() == '':
                break
            body.append(line)
        # put data into an array and instantiate as dataframe
        words = [l.strip().split() for l in body]
        words = np.array(words)
        assert words.ndim == 2, f"Error in parsing. Uneven # columns detected in parsing {set([len(x) for x in words])}." 
        assert words.shape[1] == len(headers), f"Error in parsing. Number of columns {words.shape[1]} != number of headers {len(headers)}" 
        data = {h:words[:,i] for i,h in enumerate(headers)}
        df = pd.DataFrame(data=data)
        return self(headers, df)

    def _write_block(self, f, headers, df, block_header='data_'):
        f.write(f'{block_header}\n\n')
        f.write('loop_\n')
        f.write('\n'.join(headers))
        f.write('\n')
        for i in df.index:
            # TODO: Assumes header and df ordering is consistent
            f.write(' '.join([str(v) for v in df.loc[i]]))
            f.write('\n')

    def write(self, outstar):
        f = open(outstar,'w')
        f.write('# Created {}\n'.format(dt.now()))
        f.write('\n')
        
        if self.relion31:
            self._write_block(f, self.data_optics.headers, self.data_optics.df, block_header='data_optics')
            f.write('\n\n')
            self._write_block(f, self.headers, self.df, block_header='data_particles')
        else:
            self._write_block(f, self.headers, self.df, block_header='data_')

def return_distribution_star(in_star, distribution_dict, out_star):
    #read in refinement .star file with all particles
    s = Starfile.load(in_star)
    pattern = r'([1-9][0-4]*)(bp)'
    origins = s.df['_rlnImageName'].str.extract(pattern)[0].astype('int')

    #select random indices from each dataset
    selected_inds = {}
    for dataset in distribution_dict:
        subset_inds = origins[origins == dataset].index
        assert len(subset_inds) >= distribution_dict[dataset], 'number of particles to be selected must be less than or equal to real number of particles'
    
        selected_inds[dataset] = np.random.choice(subset_inds, distribution_dict[dataset], replace = False)

    #shuffle indices to randomize .star file order
    all_selected_inds = np.concatenate([selected_inds[dataset] for dataset in distribution_dict])
    all_selected_inds_shuffle = np.random.permutation(all_selected_inds)
    s.df = s.df.loc[all_selected_inds_shuffle]
    
    #reassign random subsets
    sampling_distr = np.ones(int(np.ceil((len(all_selected_inds))/2)*2))
    sampling_distr[(len(all_selected_inds))//2:] = 2
    sampling_distr = sampling_distr.astype('int')

    s.df['_rlnRandomSubset'] = np.random.choice(sampling_distr, len(s.df), replace = False)

    #write new .star file
    s.write(out_star)
    
    return 

def main(args):
    desired_num_particles = {2: args.two, 8: args.eight, 14: args.fourteen}
    return_distribution_star(args.infile, desired_num_particles, args.outfile)
    return

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    main(add_args(parser).parse_args())