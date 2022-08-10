import numpy as np
import logging

# Class reworked from C. Hofer
class Site:    
    def __init__(self, x, y, z = 0, site_id = None):
        self.id = site_id
        self.coords = np.array([x, y])
        self.neighbors = []
        self.candidates = []
        self.heavy = False

    def relocate(self, x, y, z = 0):
        self.coords = np.array([x, y, z])                 
        
    def distance(self, site):
        x = self.coords - site.coords
        return np.sqrt(x[0]**2 + x[1]**2)

    def second_nearest_neighbors(self):
        out = []

        for n in self.neighbors:
            for nn in n.neighbors:
                if nn is not self:
                    out.append(nn)
        
        return out
    
    def add_candidate(self, site2):
        if self == site2:
            # An atom itself must not be a candidate for its neighbor
            return False
        if not (self in site2.candidates or site2 in self.candidates):
            self.candidates.append(site2)
            site2.candidates.append(self)
            return True
        else:
            return False
    
    def top_candidate(self):
        if self.candidates == []:
            return None
        diff = self.coords - np.array(list(map(lambda x: x.coords, self.candidates)))
        dist = np.sqrt(diff[:,0]**2 + diff[:,1]**2)
        for idx in np.argsort(dist):
            if self.eval_candidate(self.candidates[idx]):
                return self.candidates[idx]
        return None 

    def add_neighbor(self, bond_candidate):
        try:
            self.candidates.remove(bond_candidate)
            bond_candidate.candidates.remove(self)
        except: pass
        
        if self not in bond_candidate.neighbors and bond_candidate not in self.neighbors:
          # if self.eval_candidate(bond_candidate) and bond_candidate.eval_candidate(self): # done in top_candidate() else
            self.neighbors.append(bond_candidate)
            bond_candidate.neighbors.append(self)
            return True 
        else:
            return False

    def eval_candidate(self, bond_candidate): 
        # No more than 4 neighbors
        if len(self.neighbors)>=4:
            return False
        # No triangles
        return not any( [bond_candidate in x.neighbors for x in self.neighbors] )
    
    def output_info(self):
        a_id = self.id
        a_coords = self.coords
        return "site ID " + str(a_id) + " with the coordinates " + str(a_coords)


class Atom(object):
    
    def __init__(self, site, element, defined_by_user=False, main_path=None):
        self.debug_print = False # Some lines with print commands are inserted for debugging
        
        self.site = site
        self.origin = site
        self.element = element
        self.defined_by_user = defined_by_user
        self.main_path = main_path # ##new

    def move(self, new_site):
        old_site = self.site
        self.site = new_site
        if new_site != old_site and self.debug_print:
                txt_old = old_site.output_info()
                txt_new = self.site.output_info()
                print(" Info: Atom to be displaced from")
                print(txt_old)
                print("        to")
                print(txt_new)
    
    def output_info(self):
        a_id = self.site.id
        a_coords = self.site.coords
        return "Atom at site ID " + str(a_id) + " with the coordinates " + str(a_coords)

# Class reworked from C. Hofer  
class Bond:
    
    def __init__(self, a1, a2, id_bond = None):
        self.site1 = a1
        self.site2 = a2
        self.id_bond = None
    
    def length(self):
        return np.linalg.norm(self.site1.coords - self.site2.coords)     
        
    def coords(self):
        return (self.site1.coords, self.site2.coords)

# Class reworked from C. Hofer        
class Bonds:
    # sites ... list of class members of "Atom"
    # members ... list of class members of "Bond"
    
    def __init__(self, sites, max_bond_length, build_bonds_on_init_flag = True):
        self.members = []
        self.max_bond_length = max_bond_length # in pixels
        
        if build_bonds_on_init_flag:
            self.members = Bonds.build_bonds(sites, self.max_bond_length)
        
    def build_bonds(sites, max_bond_length):
        N_candidates = 0
        # Get all coords and delete old candidates and neighbors
        def func(x):
            x.candidates = []
            x.neighbors = []
            return x.coords
        all_coords = np.array(list(map(func, sites)))
        distances = np.linalg.norm(all_coords-all_coords[:, np.newaxis, :], axis=2)
        ind_0, ind_1 = np.where( (distances<=max_bond_length) )
        for k in range(len(ind_0)):
            if sites[ind_0[k]].add_candidate(sites[ind_1[k]]):
                N_candidates += 1
        print('Total number of candidate bonds: %d' % N_candidates)
        
        # Algorithm could be made faster
        bonds = []
        num_bonds = -1
        sites_copy = sites
        
        it = 0
        while len(bonds) > num_bonds and it < 50:
            it+=1
            num_bonds = len(bonds)
            #sites_check_omit = set()
            for i in range(len(sites_copy)):
                atc = sites_copy[i].top_candidate()
                if atc and sites_copy[i].add_neighbor(atc): #NEW functionality instead of just "if atc"
                    bonds.append(Bond(sites_copy[i], atc))   
                        
        logging.info("Number of bonds set: %d" % len(bonds))
        return bonds