import numpy as np
import logging
from scipy.optimize import linear_sum_assignment as lsa

class Path(object):
    
    def __init__(self, site1, site2, a_id=0, list_blockers = [], list_banned = [], is_subpath=False):
        self.debug_print = False # Some lines with print commands are inserted for debugging
        
        self.start = site1
        self.end = site2
        self.id = a_id
        # array of atoms
        self.list_blockers = np.array(list_blockers)
        # array of sites
        self.list_banned = np.array(list_banned)
        self.sitelist = np.array([])
        self.blocked_by = np.array([])
        # flags
        self.is_subpath = is_subpath
        self.is_valid = None
        
    def print_sitelist(self):
        txt = "Current path is (site ids): "
        for site in self.sitelist:
            txt += str(site.id)+" "
        print(txt)
            
    def blocked(self):
        out = np.array([])
        caused_by = np.array([])
        # Neighbors are also configured to be blocking
        for b in self.list_blockers:
            tmp = np.append(b.site, b.site.neighbors)
            out = np.append(out, tmp) 
            caused_by = np.append(caused_by, np.repeat(b, len(tmp)))
        return out, caused_by
    
    def add_to_path(self, site):
        if site in self.list_banned:
            print("  FAILURE: Banned path successor")
            return False
        else:
            self.sitelist.append(site)
            blocking, caused_by = self.blocked()
            if site in blocking:
                blocker = caused_by[np.where(site == blocking)][0]
                if not blocker in self.blocked_by:
                    self.blocked_by = np.append(self.blocked_by, blocker)
                    txt = blocker.output_info()
                    print("     Path blocked by %r" % txt)
            return True
        
    def ban_site(self, site):
        self.list_banned.append(site)
        
    def determine_path(self):
        self.sitelist = list([self.start]) # No legality check for the starting point
        if self.start == self.end:
            logging.info("  Info: Starting site and end site are equal.")
            self.is_valid = False
            return None

        if self.debug_print:
            print(self.start.neighbors) #DEBUG

        it = 0
        while (self.sitelist[-1].id != self.end.id) and (it <= 1000):
            #print("DEBUG1a: current last path member is "+str(self.sitelist[-1].id))
            #print("DEBUG1b: end of path should be "+str(self.end.id))
            #print("DEBUG")
            #self.print_sitelist()
            it += 1
            d = np.infty
            highest_rated_neighbor = None
            for candidate in self.sitelist[-1].neighbors:
                # Exclude any double occurence and forbid the banned sites
                # CAUTION:  Probably not working well at every instant
                #           Should probably be done in add_to_path()
                if (candidate in self.sitelist) or (candidate in self.list_banned):
                    continue
                
                # d_tilde ... Cost equivalent
                d_tilde = candidate.distance(self.end)
                # Downrate the path at an edge by 10 percent
                if len(candidate.neighbors) < 3:
                    d_tilde *= 1.1
                # Downrate an approach to a blocked position
                if candidate in self.list_blockers:
                    d_tilde *= 10 #*candidate.distance(blocker.site)
                # Refresh highest rated neighbor
                if d_tilde < d:
                    d = d_tilde
                    highest_rated_neighbor = candidate
            if not highest_rated_neighbor:
                print("  Failure: No allowed neighbor found")
                self.is_valid = False
                return None
            else:
                self.add_to_path(highest_rated_neighbor)
        self.is_valid = True

            
class Paths(object):
    
    def __init__(self, atoms, target_sites):
        self.debug_print = False # Some lines witch print commands are inserted for debugging
        
        self.members = np.array([])
        self.target_sites = np.array(target_sites) # numpy array of class member "Site"
        for atom in atoms:
            atom.site = atom.origin # reinit position of the atoms when calculating paths freshly
        self.atoms = np.array(atoms)
        
        # First determine atom-target assignment
        self.atoms_ordered_idx, self.target_sites_ordered_idx, cost = self.hungarian_lap()
        
        # Then sort the atom-target-pairs in ascending order w.r.t. distance
        self.build_succession(cost)
        self.atoms = self.atoms[self.atoms_ordered_idx]
        self.target_sites = self.target_sites[self.target_sites_ordered_idx]
        
        self.print_atoms_and_targets()
        
        
    def print_atoms_and_targets(self):
        print("=====")
        print("Atoms (sorted):")
        for atom in self.atoms:
            txt = atom.output_info()
            print(txt)
        print("=====")
        print("Target sites (sorted):")
        for site in self.target_sites:
            txt = site.output_info()
            print(txt)
        
    def hungarian_lap(self):
        #  The Hungarian algorithm for the (linear) assignment problem
        #  also known as the Munkres or Kuhn-Munkres algorithm
        #  Kuhn (1955), Munkres (1957)
        # TODO: It is probably not the least total path length
        # (sum of all path lengths) at every instant
        
        N = min(len(self.target_sites), len(self.atoms))
        if len(self.atoms) != len(self.target_sites):
            print(" Warning: Number of atoms is not equal to the number of target sites." +
                  " Reliable functionality for this case is not implemented." + 
                  " This could lead to unintended behavior and/or exceptions.")
        # Cost matrix
        C = np.zeros((N, N))
        for i in range(N):
            for j in range(N):
                C[i, j] = self.atoms[i].site.distance(self.target_sites[j])
        row_ind, col_ind = lsa(C)

        return np.array(row_ind, dtype=int), np.array(col_ind, dtype=int), C[row_ind, col_ind]
        
    def build_succession(self, cost):
        cost_ordered_idx = np.argsort(cost) # indices sorted by cost
        self.atoms_ordered_idx = self.atoms_ordered_idx[cost_ordered_idx]
        self.target_sites_ordered_idx = self.target_sites_ordered_idx[cost_ordered_idx]
        return True

    def determine_paths(self):
        for k in range(len(self.atoms)):       
            target = self.target_sites[k]
            source = self.atoms[k]
            self.members = np.append(self.members, Path(source.site, target))
            self.members[-1].determine_path()
        logging.info("%d paths determined." % len(self.members))
    
    def determine_paths_nooverlap(self):
        for k in range(len(self.atoms)):
        
            a_lb = np.delete(self.atoms, k) # list of (potential) blockers
            
            if self.debug_print:
                print("=== Potential blockers ===")
                for atom in a_lb:
                    print(atom.output_info())
                print("===")
                
            path_to_be_evaluated = Path(self.atoms[k].site, self.target_sites[k], list_blockers = a_lb)
            path_to_be_evaluated.determine_path()
            
            blocked_by = path_to_be_evaluated.blocked_by
            
            if len( blocked_by ) > 0:
                # write blocking sources into an array
                blockers = np.array([])
                for kk in range( len(blocked_by) ):
                    blocker_idx = np.where(blocked_by[kk] == self.atoms)[0][0]
                    blocker = self.atoms[blocker_idx]
                    blockers = np.append(blockers, blocker)
            
                # evaluate new paths
                new_target = self.target_sites[k]
                while len(blockers) > 0:
                    blocker = blockers[-1]
                    subpath = Path(blocker.site, new_target, is_subpath = True)
                    subpath.determine_path()
                
                    # Move atom, print new atom and targets list, add member to pathlist, print sitelist
                    new_target = blocker.site
                    blocker.move(subpath.sitelist[-1])
                    if self.debug_print:
                        self.print_atoms_and_targets()
                    self.members = np.append(self.members, subpath)
                    self.members[-1].print_sitelist()
                    
                    # delete blocker from list         
                    blockers = np.delete(blockers, -1)
                
                path_to_be_evaluated = Path(self.atoms[k].site, new_target, list_blockers = a_lb)
                path_to_be_evaluated.determine_path()
                if self.debug_print:
                    print("Length of blockers after rearrangement: "+str(len(blockers)))
                
            self.atoms[k].move(path_to_be_evaluated.sitelist[-1])
            if self.debug_print:
                self.print_atoms_and_targets()
            self.members = np.append(self.members, path_to_be_evaluated)
            self.members[-1].print_sitelist()
        
        N = 0
        for path in self.members:
            if not path.is_subpath and path.is_valid: N+=1
        logging.info("%d valid paths determined." % N)