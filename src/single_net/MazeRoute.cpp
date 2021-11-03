#include "MazeRoute.h"

ostream &operator<<(ostream &os, const Solution &sol) {
    os << "cost=" << sol.cost << ", vertex=" << sol.vertex << ", prev=" << (sol.prev ? sol.prev->vertex : -1);
    return os;
}

void MazeRoute::constructGridGraph(const vector<gr::GrBoxOnLayer> &guides) {
    GuideGridGraphBuilder graphBuilder(grNet, graph, guides);
    mergedPinAccessBoxes = grNet.getMergedPinAccessBoxes(
        [](const gr::GrPoint &point) { return gr::PointOnLayer(point.layerIdx, point[X], point[Y]); });
    graphBuilder.run(mergedPinAccessBoxes);
}

void MazeRoute::constructGridGraph(const CongestionMap &congMap) {
    CoarseGridGraphBuilder graphBuilder(grNet, graph, congMap);
    mergedPinAccessBoxes =
        grNet.getMergedPinAccessBoxes([&](const gr::GrPoint &point) { return graphBuilder.grPointToPoint(point); });
    graphBuilder.run(mergedPinAccessBoxes);
}

db::RouteStatus MazeRoute::run() {
    vertexCosts.assign(graph.getNodeNum(), std::numeric_limits<db::CostT>::max());
    pinSols.assign(mergedPinAccessBoxes.size(), nullptr);
    const int startPin = 0;

    auto status = route(startPin);
    if (!db::isSucc(status)) {
        return status;
    }

    getResult();

    return status;
}

db::RouteStatus MazeRoute::runSweep() {
    vertexCosts.assign(graph.getNodeNum(), std::numeric_limits<db::CostT>::max());
    vertexPrev.assign(graph.getNodeNum(), -1);
    pinSolv.assign(mergedPinAccessBoxes.size(), -1);
    const int startPin = 0;

    auto status = sweepRoute(startPin);
    if (!db::isSucc(status)) {
        return status;
    }
    getSweepResult();
    return status;
}

void MazeRoute::getAllVertices(int u, std::unordered_set<int> *visit) {
    (*visit).insert(u);
    for (auto direction : directions)
        if (graph.hasEdge(u, direction)) {
            int v = graph.getEdgeEndPoint(u, direction);
            if ((*visit).find(v) != (*visit).end()) continue;
            getAllVertices(v, visit);
        }
}

db::RouteStatus MazeRoute::sweepRoute(int startPin) {
    // log() << "iter = " << iter << std::endl;
    // return route(startPin);

    auto sweep = [&](int start, EdgeDirection direction) {
        int u = start;
        db::CostT costLimit = std::numeric_limits<db::CostT>::max();
        while (graph.hasEdge(u, direction)) {
            int v = graph.getEdgeEndPoint(u, direction);
            db::CostT w = graph.getEdgeCost(u, direction);
            db::CostT newCost = w + vertexCosts[u];
            if (newCost < vertexCosts[v]) {
                vertexCosts[v] = newCost;
                vertexPrev[v] = u;
            }
            u = v;
        }
    };

    // preprocess each via line and wire line
    std::unordered_set<int> visit;
    for (auto vertex : graph.getVertices(startPin))
        if (visit.find(vertex) == visit.end()) getAllVertices(vertex, &visit);
    std::vector<int> viaStart, viaEnd, wireStart, wireEnd;
    for (auto vertex : visit) {
        if (!graph.hasEdge(vertex, DOWN)) viaStart.push_back(vertex);
        if (!graph.hasEdge(vertex, UP)) viaEnd.push_back(vertex);
        if (!graph.hasEdge(vertex, BACKWARD)) wireStart.push_back(vertex);
        if (!graph.hasEdge(vertex, FORWARD)) wireEnd.push_back(vertex);
    }

    // init from startPin
    for (auto vertex : graph.getVertices(startPin)) vertexCosts[vertex] = 0;

    std::unordered_set<int> visitedPin = {startPin};
    int nPinToConnect = mergedPinAccessBoxes.size() - 1;

    for (int iter = 0; nPinToConnect != 0; iter++) {
        int dstVertex = -1;
        int dstPinIdx = -1;
        db::CostT minCost = std::numeric_limits<db::CostT>::max();

        for (int round = 0, iterNum = std::max(5, 11 - iter * 2); round < iterNum; round++) {
            // via sweep
            for (auto vertex : viaStart) sweep(vertex, UP);
            for (auto vertex : viaEnd) sweep(vertex, DOWN);
            // wire sweep
            for (auto vertex : wireStart) sweep(vertex, FORWARD);
            for (auto vertex : wireEnd) sweep(vertex, BACKWARD);
        }

        for (int pinIdx = startPin; pinIdx < mergedPinAccessBoxes.size(); pinIdx++)
            if (visitedPin.find(pinIdx) == visitedPin.end()) {
                for (auto vertex : graph.getVertices(pinIdx))
                    if (vertexCosts[vertex] < minCost) {
                        minCost = vertexCosts[vertex];
                        dstVertex = vertex;
                        dstPinIdx = pinIdx;
                    }
            }

        if (dstVertex < 0) {
            printWarnMsg(db::RouteStatus::FAIL_DISCONNECTED_GRID_GRAPH, grNet.dbNet);
            printlog(visitedPin, nPinToConnect, mergedPinAccessBoxes.size(), mergedPinAccessBoxes);
            printlog(graph.checkConn());
            graph.writeDebugFile(grNet.getName() + ".graph");
            getchar();
            return db::RouteStatus::FAIL_DISCONNECTED_GRID_GRAPH;
        }

        // update pinSolv
        pinSolv[dstPinIdx] = dstVertex;

        // mark the path to be zero
        auto tmp = dstVertex;
        while (tmp >= 0 && vertexCosts[tmp] != 0) {
            vertexCosts[tmp] = 0;
            tmp = vertexPrev[tmp];
        }

        // mark all the accessbox of the pin to be almost zero
        for (auto vertex : graph.getVertices(dstPinIdx)) {
            if (vertex == dstVertex) continue;
            vertexCosts[vertex] = 0;
            vertexPrev[vertex] = -1;
        }

        visitedPin.insert(dstPinIdx);
        nPinToConnect--;
    }

    return db::RouteStatus::SUCC_NORMAL;
}

void MazeRoute::getSweepResult() {
    grNet.gridTopo.clear();

    std::unordered_map<int, std::shared_ptr<gr::GrSteiner>> visited;

    // back track from pin to source
    for (unsigned p = 0; p < mergedPinAccessBoxes.size(); p++) {
        std::unordered_map<int, std::shared_ptr<gr::GrSteiner>> curVisited;
        int cur = pinSolv[p];
        // for (auto vertex : graph.getVertices(p)) cur = (cur == -1 || (cur >= 0 && vertexCosts[vertex] < vertexCosts[cur])) ? vertex : cur;
        std::shared_ptr<gr::GrSteiner> prevS;
        while (cur >= 0) {
            auto it = visited.find(cur);
            if (it != visited.end()) {
                // graft to an existing node
                if (prevS) {
                    gr::GrSteiner::setParent(prevS, it->second);
                }
                break;
            } else {
                // get curS
                auto point = graph.getPoint(cur);
                auto curS = std::make_shared<gr::GrSteiner>(gr::GrPoint(point.layerIdx, point[X], point[Y]),
                                                            graph.getPinIdx(cur));
                if (prevS) {
                    gr::GrSteiner::setParent(prevS, curS);
                }
                if (curVisited.find(cur) != curVisited.end()) {
                    printlog("Warning: self loop found in a path for net", grNet.getName(), "for pin", p);
                }
                curVisited.emplace(cur, curS);
                // store tree root
                if (vertexPrev[cur] < 0) {
                    grNet.gridTopo.push_back(curS);
                    break;
                }
                // prep for the next loop
                prevS = curS;
                cur = vertexPrev[cur];
            }
        }
        for (const auto &v : curVisited) visited.insert(v);
    }

    // remove redundant Steiner nodes
    for (auto &tree : grNet.gridTopo) {
        gr::GrSteiner::mergeNodes(tree);
    }
}

db::RouteStatus MazeRoute::route(int startPin) {
    // define std::priority_queue
    auto solComp = [](const std::shared_ptr<Solution> &lhs, const std::shared_ptr<Solution> &rhs) {
        return rhs->cost < lhs->cost;
    };
    std::priority_queue<std::shared_ptr<Solution>, vector<std::shared_ptr<Solution>>, decltype(solComp)> solQueue(
        solComp);

    auto updateSol = [&](const std::shared_ptr<Solution> &sol) {
        solQueue.push(sol);
        if (sol->cost < vertexCosts[sol->vertex]) vertexCosts[sol->vertex] = sol->cost;
    };

    // init from startPin
    for (auto vertex : graph.getVertices(startPin)) updateSol(std::make_shared<Solution>(0, vertex, nullptr));

    std::unordered_set<int> visitedPin = {startPin};
    int nPinToConnect = mergedPinAccessBoxes.size() - 1;

    while (nPinToConnect != 0) {
        std::shared_ptr<Solution> dstVertex;
        int dstPinIdx = -1;

        // Dijkstra
        while (!solQueue.empty()) {
            auto newSol = solQueue.top();
            int u = newSol->vertex;
            solQueue.pop();

            // reach a pin?
            dstPinIdx = graph.getPinIdx(u);
            if (dstPinIdx != -1 && visitedPin.find(dstPinIdx) == visitedPin.end()) {
                dstVertex = newSol;
                break;
            }

            // pruning by upper bound
            if (vertexCosts[u] < newSol->cost) continue;

            const db::MetalLayer &uLayer = database.getLayer(graph.getPoint(u).layerIdx);

            for (auto direction : directions) {
                if (!graph.hasEdge(u, direction) ||
                    (newSol->prev && graph.getEdgeEndPoint(u, direction) == newSol->prev->vertex))
                    continue;

                // from u to v
                int v = graph.getEdgeEndPoint(u, direction);

                // edge cost
                db::CostT w = graph.getEdgeCost(u, direction);

                db::CostT newCost = w + newSol->cost;

                if (newCost < vertexCosts[v]) updateSol(std::make_shared<Solution>(newCost, v, newSol));
            }
        }

        if (!dstVertex) {
            printWarnMsg(db::RouteStatus::FAIL_DISCONNECTED_GRID_GRAPH, grNet.dbNet);
            printlog(visitedPin, nPinToConnect, mergedPinAccessBoxes.size(), mergedPinAccessBoxes);
            printlog(graph.checkConn());
            graph.writeDebugFile(grNet.getName() + ".graph");
            getchar();
            return db::RouteStatus::FAIL_DISCONNECTED_GRID_GRAPH;
        }

        // update pinSols
        pinSols[dstPinIdx] = dstVertex;

        // mark the path to be zero
        auto tmp = dstVertex;
        while (tmp && tmp->cost != 0) {
            updateSol(std::make_shared<Solution>(0, tmp->vertex, tmp->prev));
            tmp = tmp->prev;
        }

        // mark all the accessbox of the pin to be almost zero
        for (auto vertex : graph.getVertices(dstPinIdx)) {
            if (vertex == dstVertex->vertex) continue;
            updateSol(std::make_shared<Solution>(0, vertex, nullptr));
        }

        visitedPin.insert(dstPinIdx);
        nPinToConnect--;
    }

    return db::RouteStatus::SUCC_NORMAL;
}

// GrSteiner means a Steiner node, and each vertex corresponds to one of it.
// Trace back from each pin's solution vertex, and create a Steiner node for each of it.
// Link them to be a tree, finishing building the result tree.
void MazeRoute::getResult() {
    grNet.gridTopo.clear();

    std::unordered_map<int, std::shared_ptr<gr::GrSteiner>> visited;    

    // back track from pin to source
    for (unsigned p = 0; p < mergedPinAccessBoxes.size(); p++) {
        std::unordered_map<int, std::shared_ptr<gr::GrSteiner>> curVisited;
        auto cur = pinSols[p];
        std::shared_ptr<gr::GrSteiner> prevS;
        while (cur) {
            auto it = visited.find(cur->vertex);
            if (it != visited.end()) {
                // graft to an existing node
                if (prevS) {
                    gr::GrSteiner::setParent(prevS, it->second);
                }
                break;
            } else {
                // get curS
                auto point = graph.getPoint(cur->vertex);
                auto curS = std::make_shared<gr::GrSteiner>(gr::GrPoint(point.layerIdx, point[X], point[Y]),
                                                            graph.getPinIdx(cur->vertex));
                if (prevS) {
                    gr::GrSteiner::setParent(prevS, curS);
                }
                if (curVisited.find(cur->vertex) != curVisited.end()) {
                    printlog("Warning: self loop found in a path for net", grNet.getName(), "for pin", p);
                }
                curVisited.emplace(cur->vertex, curS);
                // store tree root
                if (!(cur->prev)) {
                    grNet.gridTopo.push_back(curS);
                    break;
                }
                // prep for the next loop
                prevS = curS;
                cur = cur->prev;
            }
        }
        for (const auto &v : curVisited) visited.insert(v);
    }

    // remove redundant Steiner nodes
    for (auto &tree : grNet.gridTopo) {
        gr::GrSteiner::mergeNodes(tree);
    }
}
