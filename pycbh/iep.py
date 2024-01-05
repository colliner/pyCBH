import sys, os
from itertools import combinations


class fragweb:

    def __init__(self, node_ls):
        """
        Initializes a fragweb object with a list of nodes.

        Parameters:
        - node_ls (list): List of nodes.

        Returns:
        - None
        """
        self.n_nodes = len(node_ls)
        self.n_edges = 0
        self.node_ls = node_ls
        self.senders = list()
        self.receivers = list()
        self.neighbors = list()
        for _ in range(self.n_nodes):
            self.neighbors.append([])
        pass

    def query_edge(self, n1, n2):
        """
        Queries whether an edge between two nodes exists.

        Parameters:
        - n1 (int): Index of the first node.
        - n2 (int): Index of the second node.

        Returns:
        - bool: True if the edge exists, False otherwise.
        """
        for idx, idx1 in enumerate(self.senders):
            if idx1 == n1:
                if self.receivers[idx] == n2:
                    return True
            elif idx1 == n2:
                if self.receivers[idx] == n1:
                    return True
        return False

    def add_edge(self, n1, n2):
        """
        Adds an edge between two nodes.

        Parameters:
        - n1 (int): Index of the first node.
        - n2 (int): Index of the second node.

        Returns:
        - None
        """
        if not self.query_edge(n1, n2):
            n1, n2 = min([n1, n2]), max([n1, n2])
            self.senders.extend([n1])
            self.receivers.extend([n2])
            if n2 not in self.neighbors[n1]:
                self.neighbors[n1].append(n2)
        pass

    def groups_of_n(self):
        """
        Finds overlapping groups of nodes.

        Returns:
        - dict: Dictionary containing overlapping groups of nodes.
        """
        groups = dict()
        for idx, idx1 in enumerate(self.senders):
            idx2 = self.receivers[idx]
            for jdx in range(self.n_nodes):
                if jdx > idx1 and jdx > idx2:
                    if self.query_edge(idx1, jdx):
                        if self.query_edge(idx2, jdx):
                            #print('FOUND TRIPLE: {} {} {}'.format(idx1,idx2,jdx))
                            if 3 in groups:
                                groups[3].append(sorted([idx1, idx2, jdx]))
                            else:
                                groups[3] = [sorted([idx1, idx2, jdx])]
        #print(groups)
        done = False
        while not done:
            new_groups = dict()
            try:
                max_key = max(groups.keys())
            except:
                if len(groups.keys()) > 0:
                    max_key = groups.keys()[-1]
                else:
                    max_key = -1
            #for key in groups:
            if max_key != -1:
                for group in groups[max_key]:
                    #min_idx = max(group)
                    #for jdx in range(min_idx,self.n_nodes):
                    #print('trying candidates: {}'.format(self.neighbors[group[0]]))
                    for jdx in [
                            x for x in self.neighbors[group[0]]
                            if x > max(group)
                    ]:
                        how_many = 0
                        for x in group:
                            if x < jdx:
                                if self.query_edge(x, jdx):
                                    how_many += 1
                                else:
                                    break
                        if how_many == len(group):
                            size = len(group) + 1
                            if size in new_groups:
                                new_group = sorted(group + [jdx])
                                if new_group not in new_groups[size]:
                                    new_groups[size].append(new_group)
                            else:
                                new_groups[size] = [group + [jdx]]
            #print('\n\nnew groups: {}'.format(new_groups))
            how_many = 0
            for key in new_groups:
                if key not in groups:
                    groups[key] = list()
                if key in groups:
                    for group in new_groups[key]:
                        if sorted(group) not in groups[key]:
                            groups[key].append(sorted(group))
                            how_many += 1
                #print('{}-body overlap : {} possible groups'.format(key,len(new_groups[key])))
            if how_many == 0:
                done = True
        #print('\ngroups: {}'.format(groups))
        #print('\nsenders: {}\nreceivers: {}'.format(self.senders,self.receivers))
        return groups


def overlap(l1, l2):
    """
    Finds the overlap between two lists.
    O(m+n)

    Parameters:
    - l1 (list): First list.
    - l2 (list): Second list.

    Returns:
    - set: Set containing the overlapping elements.
    """
    return set(l1).intersection(l2)


def iep_graph(primary, graph):
    """
    Generates overlapping fragments in a graph based on the primary structure.

    Parameters:
    - primary (list): Primary structure.
    - graph: Graph information.

    Returns:
    - dict: Dictionary containing overlapping fragments.
    """
    print(primary)
    sys.exit()
    frag_list = {1: primary}
    web = fragweb(primary)
    #print(len(list(combinations(list(range(len(primary))),4))))
    #c=combinations(primary, 2)
    c = combinations(list(range(len(primary))), 2)
    #print(c)
    o_ls = list()
    for idx1, idx2 in c:
        o = overlap(primary[idx1], primary[idx2])
        if len(o) > 0:
            web.add_edge(idx1, idx2)
            #print('{} n {} = {}'.format(primary[idx1],primary[idx2],o))
            o_ls.append(list(o))
    frag_list[2] = o_ls

    #print(frag_list)
    groups = web.groups_of_n()
    for size in groups:
        if size not in frag_list:
            frag_list[size] = list()
        for group in groups[size]:
            setlist = list()
            for g in group:
                setlist.append(set(primary[g]))
            u = set.intersection(*setlist)
            if len(list(u)) > 0:
                frag_list[size].append(list(u))
        #print('{}-body overlap : found {} overlapping fragments'.format(size,len(frag_list[size])))
    #print(frag_list)
    return frag_list


def iep(primary):
    """
    Generates overlapping fragments based on the primary structure.

    Parameters:
    - primary (list): Primary structure.

    Returns:
    - dict: Dictionary containing overlapping fragments.
    """
    frag_list = {1: primary}
    web = fragweb(primary)
    #print(len(list(combinations(list(range(len(primary))),4))))
    #c=combinations(primary, 2)
    c = combinations(list(range(len(primary))), 2)
    #print(c)
    o_ls = list()
    for idx1, idx2 in c:
        o = overlap(primary[idx1], primary[idx2])
        if len(o) > 0:
            web.add_edge(idx1, idx2)
            #print('{} n {} = {}'.format(primary[idx1],primary[idx2],o))
            o_ls.append(list(o))
    frag_list[2] = o_ls

    #print(frag_list)
    groups = web.groups_of_n()
    for size in groups:
        if size not in frag_list:
            frag_list[size] = list()
        for group in groups[size]:
            setlist = list()
            for g in group:
                setlist.append(set(primary[g]))
            u = set.intersection(*setlist)
            if len(list(u)) > 0:
                frag_list[size].append(list(u))
        #print('{}-body overlap : found {} overlapping fragments'.format(size,len(frag_list[size])))
    #print(frag_list)
    return frag_list


if __name__ == '__main__':
    primary = [[0, 1, 2, 3], [1, 0, 2, 3, 4], [2, 1, 0, 3], [3, 1, 4, 0, 2, 5],
               [4, 3, 5, 1, 6], [5, 4, 6, 3], [6, 5, 4]]
    overlaps_ = [[0, 1, 2, 3], [1, 0, 2, 3], [1, 0, 2, 3, 4], [3, 1, 4, 5],
                 [4, 3, 5, 6], [5, 4, 6]]
    #primary= [[0, 1], [1, 0, 2, 3], [2, 1], [3, 1, 4], [4, 3, 5], [5, 4, 6], [6, 5]]
    #overlaps_ = [[0, 1], [1, 2], [1, 3], [3, 4], [4, 5], [5, 6]]
    primary = [[0, 1], [1, 0, 2], [2, 1, 3], [3, 2, 4], [4, 3]]
    overlaps_ = [[0, 1], [1, 2], [2, 3], [3, 4]]
    overlaps = iep(primary)
    print(overlaps)
    sys.exit()
