import itertools
import tempfile, os
from pysimm import system, lmps, forcefield
from pysimm.apps.random_walk import random_walk, copolymer
from multiprocessing import Pool
import argparse

parser = argparse.ArgumentParser(description='Generation of amorphous cell')
parser.add_argument('monomer', type=str, help='Number of chains in unit cell')
parser.add_argument('--chains', type=int, default=4, help='Number of chains in unit cell')
parser.add_argument('--size', type=str, default='6 3 1 0', help='Number of chains in unit cell')
parser.add_argument('--density', type=float, default=1.1, help='Target density of polymer')


args = parser.parse_args()
print('Requesed arguments:\n')
print('Number of chains: ', args.size)
print('Density: ', args.density)
print('File with monomer structure: ', args.monomer)
print('Chain size pattern: ', args.chains)

assert len(args.size.split()) in (1, args.chains), 'number of chain and pattern of their size must match'
assert os.path.exists(args.monomer), f'file {args.monomer} not exists'

name = 'uniform_polymer'

def _monomer(molfile,ih,it,position):
    sys = system.read_mol(molfile)
    ff = forcefield.Dreiding()
    sys.apply_forcefield(ff)
    head = sys.particles[ih]
    tail = sys.particles[it]
    head.linker = 'head'
    tail.linker = 'tail'
    if position=='tail' or position=='middle':
        for bond in head.bonds: # tail
            if bond.a.elem == 'H' or bond.b.elem == 'H':
                hatom = bond.a if bond.b is tail else bond.b
                sys.particles.remove(hatom.tag, update=False)
                break
    if position=='head' or position=='middle':  
        for bond in tail.bonds: # head 
            if bond.a.elem == 'H' or bond.b.elem == 'H':
                hatom = bond.a if bond.b is tail else bond.b
                sys.particles.remove(hatom.tag, update=False)
                break
    sys.remove_spare_bonding()
    sys.pair_style = 'lj/cut'
    lmps.quick_min(sys, min_style='fire')
    sys.add_particle_bonding()
    return sys


class Monomer():
    
    def __init__(self, id, molfile, ihead, itail):
        self.id = id # integer representing monomer type
        self.molfile = molfile
        with open(molfile,'r') as fi:
            self._mol_cont = fi.read()
        self.ihead = ihead
        self.itail = itail
        
    def get_middle(self):
        return _monomer(self.molfile,self.ihead,self.itail,'middle')
    
    def get_head(self):
        return _monomer(self.molfile,self.ihead,self.itail,'head')
    
    def get_tail(self):
        return _monomer(self.molfile,self.ihead,self.itail,'tail')
    
    def get_single(self):
        return _monomer(self.molfile,self.ihead,self.itail,'')

class Chain():
    def __init__(self, monomers, size, **kwargs):
        self.monomers = monomers # [Monomer1, Monomer2, ...]
        assert size >=2, 'size of chain must be >1'
        self.size = size
        self._kwargs = kwargs
        self.pattern = kwargs.get('pattern', None) #TODO: check if input is valid
        self.ratio = kwargs.get('ratio', None) #TODO: check if input is valid
        assert self.pattern != None, 'pattern must be explicitly specified in current implementation. For example:\n\
            1. Monomers A and B: [1,2,4,2,3,4] means [tailA, A, 2B, 4A, 2B, 3A, 4B ,headA]\n\
            2. Monomers A B C: [1,2,4,2,3,4] means [tailA, A, 2B, 4C, 2A, 3B, 4C ,headA]\n\
            3. Monomer A: [11] means [tailA, 11A, headA]\n'
        assert self.size == sum(self.pattern) + 2

    def build(self):
        '''building monomer sequence for pysimm.apps.random_walk.copolymer() method'''
        _seq = []
        _cycle = itertools.cycle(self.monomers)
        A = next(_cycle)
        head = A.get_head()
        tail = A.get_tail()
        _seq.append(tail)
        for i in self.pattern:
            mono = next(_cycle)
            mid = mono.get_middle() #TODO optimize to process only unique monomers
            _seq.append(mid)
        _seq.append(head)
        # FF setup
        f = forcefield.Dreiding()
        for _ch in _seq:
            _ch.pair_style = 'lj'
            _ch.apply_charges(f, charges='gasteiger')
        # build chain
        chain = copolymer(_seq, self.size, pattern=[1, *self.pattern, 1], forcefield=f)
        chain.apply_forcefield(f)
        chain.apply_charges(f, charges='gasteiger')
        return chain

if __name__ == '__main__':
  mol_path, ih, it = args.monomer.split()
  ih, it = int(ih), int(it)
  monomer = Monomer(0, mol_path, ih, it)
  
  _sizes = [int(i) for i in args.size.split()]
  _chains = [Chain([monomer], _sizes[i], pattern=[_sizes[i]-2]) for i in range(args.chains)]
  chains = [_ch.build() for _ch in _chains]
  uniform_polymer = system.replicate(chains, [1 for i in chains], density=args.density, rand=True)
  
  uniform_polymer.write_xyz(f'{name}.xyz')
  uniform_polymer.write_yaml(f'{name}.yaml')
  uniform_polymer.write_lammps(f'{name}.lmps')
  uniform_polymer.write_chemdoodle_json(f'{name}.json')
