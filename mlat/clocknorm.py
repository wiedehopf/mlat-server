# -*- mode: python; indent-tabs-mode: nil -*-

# Part of mlat-server: a Mode S multilateration server
# Copyright (C) 2015  Oliver Jowett <oliver@mutability.co.uk>

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.

# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
Clock normalization routines.
"""

import pygraph.classes.graph
import pygraph.algorithms.minmax
import time
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from scipy.sparse.csgraph import minimum_spanning_tree
import numpy as np

from mlat import profile


class _Predictor(object):
    """Simple object for holding prediction state"""
    def __init__(self, predict, variance):
        self.predict = predict
        self.variance = variance


def _identity_predict(x):
    return x


def _make_predictors(clocktracker, station0, station1, now):
    """Return a tuple of predictors (p_01, p_10) where:

    p_01 will predict a station1 timestamp given a station0 timestamp
    p_10 will predict a station0 timestamp given a station1 timestamp

    Returns None if no suitable clock sync model is available for
    this pair of stations.
    """
    if station0 is station1:
        return None

    if station0.bad_syncs > 0 or station1.bad_syncs > 0:
        return None

    #if station0.epoch is not None and station0.epoch == station1.epoch:
    #    # Assume clocks are closely synchronized to the epoch (and therefore to each other)
    #    predictor = _Predictor(_identity_predict, station0.clock.jitter ** 2 + station1.clock.jitter ** 2)
    #    return (predictor, predictor)

    if station0 < station1:
        pairing = clocktracker.clock_pairs.get((station0, station1))
    else:
        pairing = clocktracker.clock_pairs.get((station1, station0))

    if pairing is None or not pairing.valid:
        return None

    variance = pairing.variance
    stale = now - pairing.updated

    # increase variance for stale pairings
    variance *= 1 + stale / 60

    # increase variance for pairing with fewer sync points
    if pairing.n < 10:
        variance *= 1 + (10 - pairing.n) / 20

    if station0 < station1:
        return (_Predictor(pairing.predict_peer, variance),
                _Predictor(pairing.predict_base, variance))
    else:
        return (_Predictor(pairing.predict_base, variance),
                _Predictor(pairing.predict_peer, variance))


def _label_heights(g, node, heights):
    """Label each node in the tree with a root of 'node'
    with its height, filling the map 'heights' which
    should be initially empty."""

    # we use heights as a visited-map too.
    heights[node] = 0
    for each in g.neighbors(node):
        if each not in heights:
            _label_heights(g, each, heights)
            mn = heights[each] + g.edge_weight((node, each))
            if mn > heights[node]:
                heights[node] = mn


def _tallest_branch(g, node, heights, ignore=None):
    """Find the edge in the tree rooted at 'node' that is part of
    the tallest branch. If ignore is not None, ignore that neighbour.
    Returns (pathlen,node)"""
    tallest = (0, None)

    for each in g.neighbors(node):
        if each is ignore:
            continue

        eh = heights[each] + g.edge_weight((node, each))
        if eh > tallest[0]:
            tallest = (eh, each)

    return tallest


def _convert_timestamps(g, timestamp_map, predictor_map, node, results, conversion_chain, variance):
    """Rewrite node and all unvisited nodes reachable from node using the
    chain of clocksync objects in conversion_chain, populating the results dict.

    node: the root node to convert
    timestamp_map: dict of node -> [(timestamp, utc), ...] to convert
    results: dict of node -> (variance, [(converted timestamp, utc), ...])
    conversion_chain: list of predictor tuples to apply to node, in order
    variance: the total error introduced by chain: sum([p.variance for p in chain])
    """

    # convert our own timestamp using the provided chain
    r = []
    results[node] = (variance, r)   # also used as a visited-map
    for ts, utc in timestamp_map[node]:
        for predictor in conversion_chain:
            ts = predictor.predict(ts)
        r.append((ts, utc))

    # convert all reachable unvisited nodes using a conversion to our timestamp
    # followed by the provided chain
    for neighbor in g.neighbors(node):
        if neighbor not in results:
            predictor = predictor_map[(neighbor, node)]
            _convert_timestamps(g, timestamp_map, predictor_map,
                                neighbor,
                                results,
                                [predictor] + conversion_chain, variance + predictor.variance)


@profile.trackcpu
def normalize(clocktracker, timestamp_map):
    """
    Given {receiver: [(timestamp, utc), ...]}

    return [{receiver: (variance, [(timestamp, utc), ...])}, ...]
    where timestamps are normalized to some arbitrary base timescale within each map;
    one map is returned per connected subgraph."""

    # Represent the stations as a weighted graph where there
    # is an edge between S0 and S1 with weight W if we have a
    # sufficiently recent clock correlation between S0 and S1 with
    # estimated variance W.
    #
    # This graph may have multiple disconnected components. Treat
    # each separately and do this:
    #
    # Find the minimal spanning tree of the component. This will
    # give us the edges to use to convert between timestamps with
    # the lowest total error.
    #
    # Pick a central node of the MST to use as the the timestamp
    # basis, where a central node is a node that minimizes the maximum
    # path cost from the central node to any other node in the spanning
    # tree.
    #
    # Finally, convert all timestamps in the tree to the basis of the
    # central node.

    # populate initial graph
    g = pygraph.classes.graph.graph()
    g.add_nodes(timestamp_map.keys())

    # build a weighted graph where edges represent usable clock
    # synchronization paths, and the weight of each edge represents
    # the estimated variance introducted by converting a timestamp
    # across that clock synchronization.

    # also build a map of predictor objects corresponding to the
    # edges for later use

    now = time.time()

    predictor_map = {}
    for si in timestamp_map.keys():
        for sj in timestamp_map.keys():
            if si < sj:
                predictors = _make_predictors(clocktracker, si, sj, now)
                if predictors:
                    predictor_map[(si, sj)] = predictors[0]
                    predictor_map[(sj, si)] = predictors[1]
                    g.add_edge((si, sj), wt=predictors[0].variance)

    # find a minimal spanning tree for each component of the graph
    mst_forest = pygraph.algorithms.minmax.minimal_spanning_tree(g)

    # rebuild the graph with only the spanning edges, retaining weights
    # also note the roots of each tree as we go
    g = pygraph.classes.graph.graph()
    g.add_nodes(mst_forest.keys())
    roots = []
    for edge in mst_forest.items():
        if edge[1] is None:
            roots.append(edge[0])
        else:
            g.add_edge(edge, wt=predictor_map[edge].variance)

    # for each spanning tree, find a central node and convert timestamps
    resultComponents = []
    for root in roots:
        # label heights of nodes, where the height of a node is
        # the length of the most expensive path to a child of the node
        heights = {}
        _label_heights(g, root, heights)

        # Find the longest path in the spanning tree; we want to
        # resolve starting at the center of this path, as this minimizes
        # the maximum path length to any node

        # find the two tallest branches leading from the root
        tall1 = _tallest_branch(g, root, heights)
        tall2 = _tallest_branch(g, root, heights, ignore=tall1[1])

        # Longest path is TALL1 - ROOT - TALL2
        # We want to move along the path into TALL1 until the distances to the two
        # tips of the path are equal length. This is the same as finding a node on
        # the path within TALL1 with a height of about half the longest path.
        target = (tall1[0] + tall2[0]) / 2
        central = root
        step = tall1[1]
        while step and abs(heights[central] - target) > abs(heights[step] - target):
            central = step
            _, step = _tallest_branch(g, central, heights, ignore=central)

        # Convert timestamps so they are using the clock units of "central"
        # by walking the spanning tree edges. Then finally convert to wallclock
        # times as the last step by dividing by the final clock's frequency
        results = {}
        conversion_chain = [_Predictor(lambda x: x/central.clock.freq, central.clock.jitter**2)]
        _convert_timestamps(g, timestamp_map, predictor_map, central, results,
                            conversion_chain, central.clock.jitter**2)

        resultComponents.append(results)

    return resultComponents

class _my_component(object):
    """Simple object for holding a graph component"""
    def __init__(self, label, receivers, size):
        self.label = label
        self.receivers = receivers
        self.size = size


@profile.trackcpu
def normalize2(clocktracker, timestamp_map):
    """
    Given {receiver: [(timestamp, utc), ...]}

    return [{receiver: (variance, [(timestamp, utc), ...])}, ...]
    where timestamps are normalized to some arbitrary base timescale within each map;
    one map is returned per connected subgraph."""

    # normalize2: use scipy.sparse.csgraph instead of pythongraph

    # Represent the stations as a weighted graph where there
    # is an edge between S0 and S1 with weight W if we have a
    # sufficiently recent clock correlation between S0 and S1 with
    # estimated variance W.
    #
    # This graph may have multiple disconnected components. Treat
    # each separately and do this:
    #
    # Find the minimal spanning tree of the component. This will
    # give us the edges to use to convert between timestamps with
    # the lowest total error.
    #
    # Pick a central node of the MST to use as the the timestamp
    # basis, where a central node is a node that minimizes the maximum
    # path cost from the central node to any other node in the spanning
    # tree.
    #
    # Finally, convert all timestamps in the tree to the basis of the
    # central node.

    receivers = list(timestamp_map.keys())
    receivers.sort() # to get a matrix with only entries above the diagonal when doing si < sj

    # populate initial graph
    g = pygraph.classes.graph.graph()
    g.add_nodes(receivers)

    # build a weighted graph where edges represent usable clock
    # synchronization paths, and the weight of each edge represents
    # the estimated variance introducted by converting a timestamp
    # across that clock synchronization.

    # also build a map of predictor objects corresponding to the
    # edges for later use

    now = time.time()

    reclen = len(receivers)
    predictor_count = 0
    predictor_map = {}

    row = []
    col = []
    data = []
    ri = 0
    for si in receivers:
        ci = 0
        for sj in receivers:
            if si < sj:
                predictors = _make_predictors(clocktracker, si, sj, now)
                if predictors:
                    predictor_map[(si, sj)] = predictors[0]
                    predictor_map[(sj, si)] = predictors[1]
                    g.add_edge((si, sj), wt=predictors[0].variance)
                    predictor_count += 1
                    data.append(predictors[0].variance)
                    row.append(ri)
                    col.append(ci)
            ci += 1
        ri += 1

    if predictor_count < 2:
        return []

    cm = csr_matrix((np.array(data), (np.array(row), np.array(col))), shape=(reclen, reclen))

    #for row in cm.toarray():
    #    print(row)

    mst = minimum_spanning_tree(csgraph=cm, overwrite=True)

    n_components, labels = connected_components(csgraph=mst, directed=False, return_labels=True)
    #print('labels: ' + str(labels))

    comps = {}
    for label in labels:
        if label not in comps:
            comps[label] = (_my_component(label=label, receivers=[], size=0))

    index = 0
    for label in labels:
        comps[label].size += 1 # increment component size
        comps[label].receivers.append(receivers[index])
        index += 1

    # make our dict a list so we can sort it
    comps = list(comps.values())
    comps.sort(key=lambda x: x.size, reverse=True)

    if len(comps) == 0:
        return [] # no results, return empty list

    bigComp = comps[0] # biggest component
    roots = [bigComp.receivers[0]] # let's just stay with a list for a moment ... doesn't hurt even if we only do one entry
    #print(bigComp.size)

    if bigComp.size < 3:
        return [] # too small, don't continue

    # rebuild the graph with only the spanning edges, retaining weights
    g = pygraph.classes.graph.graph()
    g.add_nodes(bigComp.receivers)

    coo = mst.tocoo(copy=False)
    for index in range(len(coo.data)):
        weight = coo.data[index]
        si = receivers[coo.row[index]]
        sj = receivers[coo.col[index]]
        if si in bigComp.receivers:
            if si < sj:
                edge = (si,sj)
            else:
                edge = (sj,si)
            #g.add_edge((si, sj), wt=predictor_map[edge].variance)
            g.add_edge((si, sj), wt=weight)

    # for each spanning tree, find a central node and convert timestamps
    # actually we're only searching the biggest spanning tree now
    resultComponents = []
    for root in roots:
        # label heights of nodes, where the height of a node is
        # the length of the most expensive path to a child of the node
        heights = {}
        _label_heights(g, root, heights)

        # Find the longest path in the spanning tree; we want to
        # resolve starting at the center of this path, as this minimizes
        # the maximum path length to any node

        # find the two tallest branches leading from the root
        tall1 = _tallest_branch(g, root, heights)
        tall2 = _tallest_branch(g, root, heights, ignore=tall1[1])

        # Longest path is TALL1 - ROOT - TALL2
        # We want to move along the path into TALL1 until the distances to the two
        # tips of the path are equal length. This is the same as finding a node on
        # the path within TALL1 with a height of about half the longest path.
        target = (tall1[0] + tall2[0]) / 2
        central = root
        step = tall1[1]
        while step and abs(heights[central] - target) > abs(heights[step] - target):
            central = step
            _, step = _tallest_branch(g, central, heights, ignore=central)

        # Convert timestamps so they are using the clock units of "central"
        # by walking the spanning tree edges. Then finally convert to wallclock
        # times as the last step by dividing by the final clock's frequency
        results = {}
        conversion_chain = [_Predictor(lambda x: x/central.clock.freq, central.clock.jitter**2)]
        _convert_timestamps(g, timestamp_map, predictor_map, central, results,
                            conversion_chain, central.clock.jitter**2)

        resultComponents.append(results)

    return resultComponents
