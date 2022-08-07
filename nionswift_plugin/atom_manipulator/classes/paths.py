import numpy as np
import logging
import copy
from scipy.optimize import linear_sum_assignment as lsa

class Path(object):
    
    def __init__(self, site1, site2, a_id=0, list_blockers = [], list_banned = [], is_subpath=False,
                avoid_1nn=True, avoid_2nn=True):
        self.debug_print = False # Some lines with print commands are inserted for debugging.
        
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
        self.avoid_1nn = avoid_1nn
        self.avoid_2nn = avoid_2nn
        
    def print_sitelist(self):
        txt = "Current path is (site ids): "
        for site in self.sitelist:
            txt += str(site.id)+" "
        print(txt)

    def ban_site(self, site):
        self.list_banned.append(site)

    def blocking_sites(self):
        out1 = []
        out2 = []

        # Sites of atoms directly
        out0 = [x.site for x in self.list_blockers]

        # If nearest neighbors are also configured to block.
        if self.avoid_1nn:
            out1 = [x.neighbors for x in out0]

            # If second-nearest neighbors are also configured to block.
            if self.avoid_2nn:
                for i, nearest_neighbors in enumerate(out1):
                    second_neighbors = []
                    for nearest_neighbor in nearest_neighbors:
                        for n in nearest_neighbor.neighbors:
                            if n is not out0[i]:
                                second_neighbors.append(n)
                    out2.append(second_neighbors)

        return out0, out1, out2
            
    def blocked(self):
        out = np.array([])
        caused_by = np.array([])
       
        for b in self.list_blockers:
            tmp = np.array([])
            tmp = np.append(tmp, b.site)

            # If nearest neighbors are also configured to block.
            if self.avoid_1nn:
                tmp = np.append(tmp, b.site.neighbors)

                # If second-nearest neighbors are also configured to block.
                if self.avoid_2nn:
                    for n in b.site.neighbors:
                        tmp = np.append(tmp, n.neighbors)

            # Condition output.
            out = np.append(out, tmp)
            caused_by = np.append(caused_by, np.repeat(b, len(tmp)))

        return out, caused_by
           
    def direct_path_blocked_old(self):
        blocker0, blocker1, blocker2 = self.blocking_sites()

        block_codes_and_sites = []

        for site in self.sitelist_direct:

            if site in blocker0:
                block_codes_and_sites.append([0, site])


            for i, neighbors in enumerate(blocker1):
                if site in neighbors:
                    block_codes_and_sites.append([1, blocker0[i]])


            for i, second_neighbors in enumerate(blocker2):
                if site in second_neighbors:
                    block_codes_and_sites.append([2, blocker0[i]])

        return block_codes_and_sites

    def direct_path_blocked(self):
        blocker0, blocker1, blocker2 = self.blocking_sites()

        block_codes_and_sites = []

        for i, b0 in enumerate(blocker0):
            # Direct site block.
            if b0 in self.sitelist_direct:
                block_codes_and_sites.append([0, b0])
                continue
            
            # Block by nearest neighbor, if configured.
            neighbors = blocker1[i]
            is_in_path = [b1 in self.sitelist_direct for b1 in neighbors]
            if any(is_in_path):
                block_codes_and_sites.append([1, b0])
                continue

            # Block by second-nearest neighbor, if configured.
            second_neighbors = blocker2[i]
            is_in_path = [b2 in self.sitelist_direct for b2 in second_neighbors]
            if any(is_in_path):
                block_codes_and_sites.append([2, b0])
                continue

        return block_codes_and_sites
        
    def determine_direct_path(self):
        self.sitelist_direct = list([self.start]) # No legality check for the starting point.
        
        if self.start == self.end:
            logging.info("  Info: Starting site and end site are equal.")
            self.is_valid = False
            return None

        if self.debug_print:
            print(self.start.neighbors)

        it = 0
        while (self.sitelist_direct[-1].id != self.end.id) and (it <= 1000):
            it += 1
            d = np.infty
            highest_rated_neighbor = None
            for candidate in self.sitelist_direct[-1].neighbors:
                # Exclude any double occurence and forbid the banned sites.

                if (candidate in self.sitelist_direct) or (candidate in self.list_banned):
                    continue
                
                # d_tilde ... Cost equivalent
                d_tilde = candidate.distance(self.end)

                # Downrate the path at an edge or through a two-coordinated atom (lattice defect?) by 10 percent.
                if len(candidate.neighbors) < 3:
                    d_tilde *= 1.1

                # Refresh highest rated neighbor.
                if d_tilde < d:
                    d = d_tilde
                    highest_rated_neighbor = candidate

            if not highest_rated_neighbor:
                print("  Failure: No allowed neighbor found")
                self.is_valid = False
                return None

            else:
                self.sitelist_direct.append(highest_rated_neighbor)

        self.is_valid = True

    def determine_unblocked_path(self):

        self.sitelist = list([self.start]) # No legality check for the starting point.

        self.candidates_d_tilde = list() # List with cost equivalents (= d_tilde).
        self.calculated_Flags = list() # Flags if calculations were already done.

        if self.start == self.end:
            logging.info("  Info: Starting site and end site are equal.")
            self.is_valid = False
            return None

        if self.debug_print:
            print(self.start.neighbors)

        it = 0

        while (self.sitelist[-1].id != self.end.id) and (it <= 100):
            it += 1 # Endless loop protection.

            path_depth = len(self.sitelist)-1 # Current path depth.
            if path_depth == -1:
                print(" No allowed unblocked path")
                self.is_valid = False
                return None

            if len(self.calculated_Flags) <= path_depth:
                self.calculated_Flags.append(False)
                self.candidates_d_tilde.append(None)
            
            if self.calculated_Flags[path_depth] is not True:

                self.candidates_d_tilde[path_depth] = np.full(len(self.sitelist[-1].neighbors), np.infty)

                for i, candidate in enumerate(self.sitelist[-1].neighbors):

                    # Here, do not allow to go back.
                    if (candidate in self.sitelist) or (candidate in self.list_banned):
                        continue

                    blocking_sites_tmp, caused_by = self.blocked()

                    # Do not allow approaches to a blocked position, except if end is a blocking site.
                    if candidate in blocking_sites_tmp:
                        
                        idx = np.nonzero(self.end==blocking_sites_tmp)

                        if len(idx[0]) > 0: # Here, self.end is in blocking_sites_tmp.
                            remaining_distance = candidate.distance(self.end)
                            if remaining_distance < candidate.distance(caused_by[idx][0].site)*0.9:
                                pass
                            else:
                                continue
                        else:
                            continue
                
                # d_tilde ... Cost equivalent
                    self.candidates_d_tilde[path_depth][i] = candidate.distance(self.end)

                # Downrate the path at an edge or through a two-coordinated atom (lattice defect?) by 10 percent.
                    if len(candidate.neighbors) < 3:
                        self.candidates_d_tilde[path_depth][i] *= 1.1

                self.calculated_Flags[path_depth] = True

            # Refresh highest rated candidate.
            d_min = np.min(self.candidates_d_tilde[path_depth])
            if d_min < np.infty:
                highest_rated_candidate = self.sitelist[-1].neighbors[np.argmin(self.candidates_d_tilde[path_depth])]
                self.sitelist.append(highest_rated_candidate)
                try:
                    self.calculated_Flags[path_depth+1] = False
                except:
                    pass
            else:
                highest_rated_candidate = None
                print("  Going back: No allowed neighbor found")
                # Cost equivalents from previous candidates.
                d_tilde_parent = self.candidates_d_tilde[path_depth-1]
                # Do not allow the highest rated previous candidate any more.
                self.candidates_d_tilde[path_depth-1][np.argmin(d_tilde_parent)] = np.infty
                self.sitelist.remove(self.sitelist[-1])

        self.is_valid = True
    
class Paths(object):
    
    def __init__(self, atoms, target_sites):
        self.debug_print = False # Some lines with print commands are inserted for debugging.
        
        self.members = np.array([]) # numpy.ndarray of class member "Path".

        self.target_sites = np.array(target_sites) # numpy.ndarray of class member "Site".
        for atom in atoms:
            atom.site = atom.origin # Reinit position of the atoms before a fresh calculation of the paths.
       
        self.atoms = np.array(atoms)  # numpy.ndarray of class member "Atom".

        self.swapped_pairs = []
        
        # First determine atom-target assignment.
        self.atoms_ordered_idx, self.target_sites_ordered_idx, cost = self.hungarian_lap()
        
        # Then sort the atom-target-pairs in ascending order w.r.t. distance.
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
        #  ----------- scipy.optimize.linear_sum_assignment ---------
        #  https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.linear_sum_assignment.html
        #
        #  Scipy 0.18.1
        #  ------------
        #  The Hungarian algorithm for the (linear) assignment problem
        #  also known as the Munkres or Kuhn-Munkres algorithm
        #  Kuhn (1955), Munkres (1957)        
        #
        #  Scipy 1.8.1
        #  ------------
        #  Modified Jonker-Volgenant algorithm with no initialization. Described in:
        #  DF Crouse. On implementing 2D rectangular assignment algorithms. 
        #  IEEE Transactions on Aerospace and Electronic Systems, 52(4):1679-1696, August 2016,
        #  DOI: 10.1109/TAES.2016.140952


        ## Symmetrical problem
        # N = min(len(self.target_sites), len(self.atoms))
        # M = N
        #if len(self.atoms) != len(self.target_sites):
        #    print(" Warning: Number of atoms is not equal to the number of target sites." +
        #          " Reliable functionality for this case is not implemented." + 
        #          " This could lead to unintended behavior and/or exceptions.")


        ## Asymmetrical problem
        # Computational effort is raised by (max(M,N) choose min(M,N)) if M is uneqal to N
        M = len(self.atoms)
        N = len(self.target_sites)


        # Cost matrix
        C = np.zeros((M, N))
        for i in range(M):
            for j in range(N):
                C[i, j] = self.atoms[i].site.distance(self.target_sites[j])
        row_ind, col_ind = lsa(C)

        return np.array(row_ind, dtype=int), np.array(col_ind, dtype=int), C[row_ind, col_ind]
        
    def build_succession(self, cost):
        cost_ordered_idx = np.argsort(cost) # Indices sorted by cost.
        self.atoms_ordered_idx = self.atoms_ordered_idx[cost_ordered_idx]
        self.target_sites_ordered_idx = self.target_sites_ordered_idx[cost_ordered_idx]
        return True

    def determine_paths(self):
        for k in range(len(self.atoms)):       
            target = self.target_sites[k]
            source = self.atoms[k]
            self.members = np.append(self.members, Path(source.site, target))
            self.members[-1].determine_direct_path()
        logging.info("%d paths determined." % len(self.members))
    
    def determine_paths_no_collision(self, avoid_1nn=True, avoid_2nn=True):
        k = 0
        while k < len(self.atoms):

            ## EXPERIMENTAL
            swapped_path = False
        
            a_lb = np.delete(self.atoms, k) # List of (potential) blockers.
            
            if self.debug_print:
                print("=== Potential blockers ===")
                for atom in a_lb:
                    print(atom.output_info())
                print("===")
                
            path_to_be_evaluated = Path(self.atoms[k].site, self.target_sites[k], list_blockers = a_lb,
                                        avoid_1nn=avoid_1nn, avoid_2nn=avoid_2nn)
            self.atoms[k].main_path = path_to_be_evaluated # ## EXPERIMENTAL
            
            ## New no collision algorithm.
            path_to_be_evaluated.determine_direct_path()
            N_steps_direct = len(path_to_be_evaluated.sitelist_direct)-1
            block_codes_and_sites = path_to_be_evaluated.direct_path_blocked()

            if self.debug_print:
                print(f"Single direct path length: {N_steps_direct}")
                print(f"Block codes {block_codes_and_sites}")

            if block_codes_and_sites: # Direct path is blocked.

                path_to_be_evaluated.determine_unblocked_path()
                N_steps_unblocked_tot = len(path_to_be_evaluated.sitelist)-1

                N_steps_compound_path = 0
                is_first_iteration = True

                while block_codes_and_sites: # Iterate for as long as a collision is identified.

                    if not is_first_iteration:
                        path_to_be_evaluated.determine_unblocked_path()
                    else:
                        is_first_iteration = False

                    N_steps_unblocked = len(path_to_be_evaluated.sitelist)-1

                    if N_steps_unblocked > len(path_to_be_evaluated.sitelist_direct)-1:
                        #print("unblocked path is longer")

                        block_code, block_site = block_codes_and_sites[-1]
                        block_atom_idx = np.where([block_site == x.site for x in self.atoms])[0][0]
                        block_atom = self.atoms[block_atom_idx]
                        block_path = block_atom.main_path # ## EXPERIMENTAL
                        
                        if block_code == 0:
                            # If the blocker is directly a foreign atom, the sum of the paths
                            # with interchanged target sites will always be equally long.

                            #print("direct path block")

                            site_idx = np.where([block_site == x for x in path_to_be_evaluated.sitelist_direct])[0][0]

                            subpath = Path(block_atom.site, path_to_be_evaluated.sitelist_direct[-1], is_subpath=True,
                                            avoid_1nn=avoid_1nn, avoid_2nn=avoid_2nn)
                            subpath.sitelist = path_to_be_evaluated.sitelist_direct[site_idx:]

                            # Move atom in backend.
                            block_atom.move(subpath.sitelist[-1])
                            
                            # Add subpath to path members.
                            if self.debug_print:
                                self.print_atoms_and_targets()
                            self.members = np.append(self.members, subpath)
                            self.members[-1].print_sitelist()
                            N_steps_compound_path += len(subpath.sitelist)-1 

                            # Prepare for next loop iteration.
                            path_to_be_evaluated.sitelist_direct = path_to_be_evaluated.sitelist_direct[:site_idx+1]
                            path_to_be_evaluated.sitelist = path_to_be_evaluated.sitelist_direct
                            path_to_be_evaluated.end = path_to_be_evaluated.sitelist_direct[-1]
                            block_codes_and_sites.remove(block_codes_and_sites[-1])

                        else:
                            #print("indirect path block")
                            
                            ## EXPERIMENTAL, swapping paths ##
                            if block_path in self.members:
                                path_member_idx = np.where(block_path == np.array(self.members))[0][0]

                                # Avoid swapping back and forth.
                                swapped_pair = [self.atoms[k], self.atoms[block_atom_idx]]
                                if swapped_pair in self.swapped_pairs or swapped_pair[-1::-1] in self.swapped_pairs:
                                    pass
                                
                                else: # Here, the order of the paths is exchanged.
                                    self.swapped_pairs.append(swapped_pair)

                                    self.atoms[k], self.atoms[block_atom_idx] = self.atoms[block_atom_idx], self.atoms[k]
                                    self.target_sites[k], self.target_sites[block_atom_idx] = self.target_sites[block_atom_idx], self.target_sites[k]
                                
                                    # Reset the position of atoms beginning at block_atom_idx.
                                    for kk in range(block_atom_idx, len(self.atoms)):
                                        self.atoms[kk].site = self.atoms[kk].origin
                                    
                                    # Delete all path members beginning at the path_member_idx of block_path and continue loop there.
                                    self.members = self.members[:path_member_idx] # 
                                    k = block_atom_idx
                                    swapped_path = True
                                    break
                            ## EXPERIMENTAL END ##

                            # Try exchanging target sites between blocking atom and currently evaluated atom.
                            subpath = Path(block_atom.site, path_to_be_evaluated.sitelist_direct[-1], is_subpath=True,
                                            avoid_1nn=avoid_1nn, avoid_2nn=avoid_2nn)
                            subpath.determine_direct_path()
                            subpath.sitelist = subpath.sitelist_direct
                            N0 = len(subpath.sitelist)-1 

                            path_proposed = Path(self.atoms[k].site, block_site, list_blockers=a_lb,
                                                    avoid_1nn=avoid_1nn, avoid_2nn=avoid_2nn)
                            path_proposed.determine_direct_path()
                            path_proposed.sitelist = path_proposed.sitelist_direct
                            N1 = len(path_proposed.sitelist)-1

                            if N_steps_unblocked <= (N0+N1): # Here, unblocked path is smaller or equal than the compound path.
                                block_codes_and_sites = []

                            else: # Here, compound path will be picked.

                                # Move atom in backend. 
                                block_atom.move(subpath.sitelist[-1])
                                
                                # Add subpath to path members.
                                if self.debug_print:
                                    self.print_atoms_and_targets()
                                self.members = np.append(self.members, subpath)
                                self.members[-1].print_sitelist()
                                N_steps_compound_path += N0

                                # Prepare for next loop iteration.
                                path_to_be_evaluated = path_proposed
                                block_codes_and_sites = path_proposed.direct_path_blocked()

                    else: # Here, the unblocked path is equally long.
                        block_codes_and_sites = []
                        if is_first_iteration:
                            N_steps_compound_path = np.nan # No compound path was calculated.

                ## EXPERIMENTAL
                if swapped_path: # Start a new loop iteration with the loop counter "k" reset to the original block_atom_idx. 
                    continue

                N_steps_compound_path += len(path_to_be_evaluated.sitelist)-1

                if self.debug_print:
                    print(f"Single non-blocking path length: {N_steps_unblocked_tot}")
                    print(f"Sum of paths with rearranged target sites: {N_steps_compound_path}")

            else: # Direct path is not blocked.
                path_to_be_evaluated.sitelist = path_to_be_evaluated.sitelist_direct
            
            # Move atom in backend.
            self.atoms[k].move(path_to_be_evaluated.sitelist[-1])

            # Add path to path members.
            if self.debug_print:
                self.print_atoms_and_targets()
            self.members = np.append(self.members, path_to_be_evaluated)
            self.members[-1].print_sitelist()
            
            # Loop counter.
            k += 1
        
        N = 0
        for path in self.members:
            if not path.is_subpath and path.is_valid: N+=1
        logging.info("%d valid paths determined." % N)