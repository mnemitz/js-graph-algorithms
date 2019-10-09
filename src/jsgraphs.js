const less = function(a1, a2, compare) {
    return compare(a1, a2) < 0;
};

const exchange = function(a, i, j) {
    var temp = a[i];
    a[i] = a[j];
    a[j] = temp;
};

const StackNode = function (value) {
    this.value = value;
    this.next = null;
};


const Stack = function() {
    this.N = 0;
    this.first = null;
};

Stack.prototype.push = function (a) {
    this.first = this._push(this.first, a);  
};

Stack.prototype._push = function(x, a) {
    if(x == null) {
        this.N++;
        return new StackNode(a);
    }  
    var oldX = x;
    this.N++;
    x = new StackNode(a);
    x.next = oldX;
    return x;
};

Stack.prototype.pop = function () {
    if (this.first == null) {
        return undefined;
    }  
    
    var oldFirst = this.first;
    var item = oldFirst.value;
    this.first = oldFirst.next;
    this.N--;
    
    return item;
};

Stack.prototype.size = function() {
    return this.N;  
};

Stack.prototype.isEmpty = function() {
    return this.N == 0;  
};

Stack.prototype.peep = function() {
    if (this.first == null) {
        return undefined;
    }  
    
    return this.first.value;
};

Stack.prototype.toArray = function() {
    var result = [];
    var x = this.first;
    while (x != null) {
        result.push(x.value);
        x = x.next;
    }
    return result;
};


var QueueNode = function(a) {
    this.value = a;  
    this.next = null;
};


var Queue = function() {
    this.first = null;
    this.last = null;
    this.N = 0;
};

Queue.prototype.enqueue = function(item) {
    var oldLast = this.last;
    this.last = new QueueNode(item);
    if(oldLast != null){
        oldLast.next = this.last;
    }
    if (this.first == null) {
        this.first = this.last;
    }
    this.N++;
};

Queue.prototype.dequeue = function() {
    if(this.first == null) {
        return undefined;
    }  
    
    var oldFirst = this.first;
    var item = oldFirst.value;
    this.first = oldFirst.next;
    
    if(this.first == null) {
        this.last = null;
    }
    
    this.N--;
    
    return item;
};

Queue.prototype.size = function() {
    return this.N;  
};

Queue.prototype.isEmpty = function() {
    return this.N == 0;
};

Queue.prototype.toArray = function() {
    var result = [];
    var x = this.first;
    while (x != null) {
        result.push(x.value);
        x = x.next;
    }
    return result;
};


var MinPQ = function(compare) {
    this.s = [];
    this.N = 0;
    if(!compare) {
        compare = function(a1, a2) {
            return a1 - a2;
        };
    }
    this.compare = compare;
};

MinPQ.prototype.enqueue = function(item) {
    while(this.s.lengh <= this.N+1) {
        this.s.push(0);
    }   
    this.s[++this.N] = item;
    this.swim(this.N);
};

MinPQ.prototype.swim = function(k) {
    while (k > 1){
        var parent = Math.floor(k / 2);
        if(less(this.s[k], this.s[parent], this.compare)){
            exchange(this.s, k, parent);
            k = parent;
        } else {
            break;
        }
    }
};

MinPQ.prototype.delMin = function() {
    if(this.N == 0) {
        return undefined;
    }  
    
    var item = this.s[1];
    exchange(this.s, 1, this.N--);
    this.sink(1);
    return item;
};

MinPQ.prototype.sink = function(k) {
    while(k * 2 <= this.N) {
        var child = 2 * k;
        if(child < this.N && less(this.s[child+1], this.s[child], this.compare)){
            child++;
        }
        if(less(this.s[child], this.s[k], this.compare)){
            exchange(this.s, child, k);
            k = child;
        } else {
            break;
        }
    }  
};

MinPQ.prototype.size = function(){
    return this.N;
};

MinPQ.prototype.isEmpty = function() {
    return this.N == 0;
};


var QuickUnion = function(V) {
    this.id = [];
    for (var v = 0; v < V; ++v) {
        this.id.push(v);
    }
};

QuickUnion.prototype.union = function(v, w) {
    var q = this.root(v);
    var p = this.root(w);
    
    if(p != q) {
        this.id[p] = q;
    }
};

QuickUnion.prototype.root = function(q) {
    while(this.id[q] != q) {
        q = this.id[q];
    }  
    return q;
};

QuickUnion.prototype.connected = function(v, w) {
    return this.root(v) == this.root(w);  
};


var IndexMinPQ = function(N, compare) {
    this.keys = [];
    this.pq = [];
    this.qp = []; // positions of key in pq
    for(var i = 0; i <= N; ++i) {
        this.keys.push(null);
        this.pq.push(0);
        this.qp.push(-1);
    }
    this.N = 0;
    
    if(!compare) {
        compare = function(a1, a2) {
            return a1 - a2;  
        };
    } 
    this.compare = compare;
};

IndexMinPQ.prototype.insert = function (index, key) {
    this.keys[index] = key;
    
    this.pq[++this.N] = index;
    this.qp[index] = this.N;
    this.swim(this.N);
};

IndexMinPQ.prototype.decreaseKey = function(index, key) {
    if(less(key, this.keys[index], this.compare)){
        this.keys[index] = key;
        this.swim(this.qp[index]);
    }
};

IndexMinPQ.prototype.minKey = function() {
    return this.keys[this.pq[1]];  
};

IndexMinPQ.prototype.min = function() {
    return this.pq[1];  
};

IndexMinPQ.prototype.delMin = function() {
    var key = this.pq[1];
    exchange(this.pq, 1, this.N);
    this.qp[this.pq[1]] = 1;
    
    this.qp[this.pq[this.N]] = -1;
    this.keys[this.pq[this.N]] = null;

    this.N--;
    
    this.sink(1);
    
    return key;
};

IndexMinPQ.prototype.swim = function (k) {
    while( k > 1) {
        var parent = Math.floor(k / 2);
        if(less(this.keys[this.pq[k]], this.keys[this.pq[parent]], this.compare)){
            exchange(this.pq, k, parent);
            this.qp[this.pq[k]] = k;
            this.qp[this.pq[parent]] = parent;
            k = parent;
        } else {
            break;
        }
    }  
};

IndexMinPQ.prototype.sink = function (k) {
    while(2 * k <= this.N) {
        var child = k * 2;
        if(child < this.N && less(this.keys[this.pq[child+1]], this.keys[this.pq[child]], this.compare)){
            child++;
        }
        if(less(this.keys[this.pq[child]], this.keys[this.pq[k]], this.compare)) {
            exchange(this.pq, k, child);
            this.qp[this.pq[k]] = k;
            this.qp[this.pq[child]] = child;
            k = child;
        } else {
            break;
        }
    }  
};

IndexMinPQ.prototype.containsIndex = function (index) {
    return this.qp[index] != -1;  
};

IndexMinPQ.prototype.isEmpty = function() {
    return this.N == 0;  
};

IndexMinPQ.prototype.size = function() {
    return this.N;
}


var Graph = function (V) {
    this.V = V;
    this.adjList = [];
    this.nodeInfo = [];
    this.edges = {};
    for (var i = 0; i < V; ++i) {
        this.adjList.push([]);
        this.nodeInfo.push({});
    }
};

Graph.prototype.addEdge = function(v, w){
    this.adjList[v].push(w);
    this.adjList[w].push(v);
    var edge_id = v + '_' + w;
    if(v > w) {
        edge_id = w + '_' + v;
    }
    this.edges[edge_id] = new Edge(v, w, 0);
};

Graph.prototype.adj = function(v) {
    return this.adjList[v];  
};

Graph.prototype.node = function(v) {
    return this.nodeInfo[v];  
};

Graph.prototype.edge = function(v, w) {
    var edge_id = v + '_' + w;
    if(v > w) {
        edge_id = w + '_' + v;
    }
    if (edge_id in this.edges) {
        return this.edges[edge_id];
    }
    return null;
};


var DiGraph = function(V) {
    this.V = V;
    this.adjList = [];
    this.nodeInfo = [];
    this.edges = {};
    for (var v = 0; v < V; ++v){
        this.adjList.push([]);
        this.nodeInfo.push({});
    }
};

DiGraph.prototype.addEdge = function(v, w){
    this.adjList[v].push(w);
    var edge_id = v + '_' + w;
    this.edges[edge_id] = new Edge(v, w, 0);
};

DiGraph.prototype.edge = function(v, w) {
    var edge_id = v + '_' + w;
    if(edge_id in this.edges) {
        return this.edges[edge_id];
    } else {
        return null;
    }
};

DiGraph.prototype.adj = function(v) {
    return this.adjList[v];  
};

DiGraph.prototype.node = function(v) {
    return this.nodeInfo[v];  
};

DiGraph.prototype.reverse = function(){
    var g = new DiGraph(this.V);
    for (var v = 0; v < this.V; ++v) {
        var adj_v = this.adjList[v];
        for (var i = 0; i < adj_v.length; ++i){
            var w = adj_v[i];
            g.addEdge(w, v);
        }
    }
    return g;
};


var Edge = function(v, w, weight) {
    this.v = v;
    this.w = w;
    this.weight = weight;
};

Edge.prototype.either = function() {
    return this.v;
};

Edge.prototype.other = function(x) {
    return x == this.v ? this.w : this.v;
};

Edge.prototype.from = function() {
    return this.v;
};

Edge.prototype.to = function() {
    return this.w;
};


var WeightedGraph = function(V) {
    this.V = V;
    this.adjList = [];
    this.nodeInfo = [];
    
    for ( var v = 0; v < V; ++v) {
        this.adjList.push([]);
        this.nodeInfo.push({});
    }
};

WeightedGraph.prototype.adj = function(v) {
    return this.adjList[v];  
};

WeightedGraph.prototype.edge = function(v, w) {
    var adj_v = this.adjList[v];
    for(var i=0; i < adj_v.length; ++i) {
        var x = adj_v[i].other(v);
        if(x == w) {
            return adj_v[i];
        }
    }
    return null;
};

WeightedGraph.prototype.node = function(v) {
    return this.nodeInfo[v];  
};

WeightedGraph.prototype.addEdge = function(e) {
    var v = e.either();
    var w = e.other(v);
    this.adjList[v].push(e);
    this.adjList[w].push(e);
};


var WeightedDiGraph = function(V) {
    WeightedGraph.call(this, V);
};

WeightedDiGraph.prototype = Object.create(WeightedGraph.prototype);

WeightedDiGraph.prototype.addEdge = function(e) {
    var v = e.from();
    this.adjList[v].push(e);
};

WeightedDiGraph.prototype.edge = function(v, w) {
    var adj_v = this.adjList[v];
    for(var i=0; i < adj_v.length; ++i) {
        var x = adj_v[i].other(v);
        if(x == w) {
            return adj_v[i];
        }
    }
    return null;
};

WeightedDiGraph.prototype.toDiGraph = function() {
    var g = new DiGraph(this.V);
    for(var v = 0; v < this.V; ++v) {
        var adj_v = this.adjList[v];
        for (var i =0; i < adj_v.length; ++i) {
            var e = adj_v[i];
            var w = e.other(v);
            g.addEdge(v, w);
        }
    }
    return g;
};


var FlowEdge = function(v, w, capacity) {
    this.v = v;
    this.w = w;
    this.capacity = capacity;
    this.flow = 0;
};

FlowEdge.prototype.residualCapacityTo = function (x) {
    if(x == this.v) {
        return this.flow;
    } else {
        return this.capacity - this.flow;
    }
};

FlowEdge.prototype.addResidualFlowTo = function (x, deltaFlow) {
    if(x == this.v) {
        this.flow -= deltaFlow;
    } else if(x == this.w) {
        this.flow += deltaFlow;
    }
};

FlowEdge.prototype.from = function() {
    return this.v;
};

FlowEdge.prototype.to = function() {
    return this.w;
};

FlowEdge.prototype.other = function(x) {
    return x == this.v ? this.w : this.v;
}


var FlowNetwork = function(V) {
    this.V = V;
    this.adjList = [];
    this.nodeInfo = [];
    for(var v = 0; v < V; ++v) {
        this.adjList.push([]);
        this.nodeInfo.push({});
    }
};

FlowNetwork.prototype.node = function(v) {
    return this.nodeInfo[v];
};

FlowNetwork.prototype.edge = function(v, w) {
    var adj_v = this.adjList[v];
    for(var i=0; i < adj_v.length; ++i) {
        var x = adj_v[i].other(v);
        if(x == w) {
            return adj_v[i];
        }
    }
    return null;
};

FlowNetwork.prototype.addEdge = function(e) {
    var v = e.from();
    this.adjList[v].push(e);
    var w = e.other(v);
    this.adjList[w].push(e);
};

FlowNetwork.prototype.adj = function(v) {
    return this.adjList[v];  
};


var DepthFirstSearch = function(G, s) {
    this.s = s;
    var V = G.V;
    this.marked = [];
    this.edgeTo = [];
    for (var v = 0; v < V; ++v) {
        this.marked.push(false);
        this.edgeTo.push(-1);
    }
    
    this.dfs(G, s);
};

DepthFirstSearch.prototype.dfs = function (G, v) {
    this.marked[v] = true;
    var adj_v = G.adj(v);
    for (var i = 0; i < adj_v.length; ++i){
        var w = adj_v[i];
        if (!this.marked[w]){
            this.edgeTo[w] = v;
            this.dfs(G, w);
        }
    }
};

DepthFirstSearch.prototype.hasPathTo = function(v) {
    return this.marked[v];
};

DepthFirstSearch.prototype.pathTo = function(v) {
    var path = new Stack();
    if(v == this.s) return [v];
    
    for(var x = v; x != this.s; x = this.edgeTo[x]) {
        path.push(x);
    }
    path.push(this.s);
    return path.toArray();
};


var BreadthFirstSearch = function(G, s) {
    var V = G.V;
    this.s = s;
    
    var queue = new Queue();
    queue.enqueue(s);
    this.marked = [];
    this.edgeTo = [];
    
    for(var v = 0; v < V; ++v) {
        this.marked.push(false);
        this.edgeTo.push(-1);
    }
    
    while (!queue.isEmpty()) {
        var v = queue.dequeue();
        this.marked[v] = true;
        var adj_v = G.adj(v);
        for (var i = 0; i < adj_v.length; ++i) {
            var w = adj_v[i];
            if(!this.marked[w]){
                this.edgeTo[w] = v;
                queue.enqueue(w);
            }
        }
    }
};

BreadthFirstSearch.prototype.hasPathTo = function(v) {
    return this.marked[v];
};

BreadthFirstSearch.prototype.pathTo = function(v) {
    var path = new Stack();
    if(v == this.s) return [v];
    
    for(var x = v; x != this.s; x = this.edgeTo[x]) {
        path.push(x);
    }
    path.push(this.s);
    return path.toArray();
};


var ConnectedComponents = function(G) {
    this.count = 0;
    var V = G.V;
    this.marked = [];
    this.id = [];
    for (var v = 0; v < V; ++v) {
        this.marked.push(false);
        this.id.push(-1);
    }
    
    for (var v = 0; v < V; ++v) {
        if(!this.marked[v]){
            this.dfs(G, v);
            this.count++;
        }
    }
};

ConnectedComponents.prototype.dfs = function(G, v) {
    this.marked[v] = true;
    this.id[v] = this.count;
    var adj_v = G.adj(v);
    
    for(var i = 0; i < adj_v.length; ++i){
        var w = adj_v[i];
        if(!this.marked[w]){
            this.dfs(G, w);
        }
    }
};

ConnectedComponents.prototype.componentId = function(v) {
    return this.id[v];
};

ConnectedComponents.prototype.componentCount = function(){
    return this.count;
};



var TopologicalSort = function(G) {
    this.postOrder = new Stack();
    this.marked = [];
    var V = G.V;
    for (var v = 0; v < V; ++v) {
        this.marked.push(false);
    }
    
    for (var v = 0; v < V; ++v) {
        if(!this.marked[v]) {
            this.dfs(G, v);
        }
    }
};

TopologicalSort.prototype.dfs = function(G, v) {
    this.marked[v] = true;
    var adj_v = G.adj(v);
    for (var i = 0; i < adj_v.length; ++i) {
        var w = adj_v[i];
        if(!this.marked[w]){
            this.dfs(G, w);
        }
    }
    this.postOrder.push(v);
};

TopologicalSort.prototype.order = function() {
    return this.postOrder.toArray();  
};


var StronglyConnectedComponents = function(G) {
    var V = G.V;
    this.count = 0;
    this.marked = [];
    this.id = [];
    
    for(var v = 0; v < V; ++v) {
        this.marked.push(false);
        this.id.push(-1);
    }
    
    var order = new TopologicalSort(G.reverse()).order();
    for( var i = 0; i < order.length; ++i) {
        var v = order[i];
        if(!this.marked[v]){
            this.dfs(G, v);
            this.count++;
        }
    }
};

StronglyConnectedComponents.prototype.dfs = function (G, v) {
    this.marked[v] = true;
    this.id[v] = this.count;
    var adj_v = G.adj(v);
    for (var i = 0; i < adj_v.length; ++i){
        var w = adj_v[i];
        if(!this.marked[w]){
            this.dfs(G, w);
        }
    }
};


StronglyConnectedComponents.prototype.componentId = function(v) {
    return this.id[v];
};

StronglyConnectedComponents.prototype.componentCount = function(){
    return this.count;
};


var KruskalMST = function(G) {
    var V = G.V;
    var pq = new MinPQ(function(e1, e2){
        return e1.weight - e2.weight;
    });
    for(var v = 0; v < G.V; ++v){
        var adj_v = G.adj(v);
        for (var i = 0; i < adj_v.length; ++i) {
            var e = adj_v[i];
            if(e.either() != v) {
                continue;
            }
            pq.enqueue(e);
        }
    }
    
    this.mst = [];
    
    var uf = new QuickUnion(V);
    while (!pq.isEmpty() && this.mst.length < V-1) {
        var e = pq.delMin();
        var v = e.either();
        var w = e.other(v);
        
        if(!uf.connected(v, w)){
            uf.union(v, w);
            this.mst.push(e);
        }
    }
};


var LazyPrimMST = function(G) {
    var V = G.V;
    this.marked = [];
    for( var v = 0; v < V; ++v) {
        this.marked.push(false);
    }
    
    this.pq = new MinPQ(function(e1, e2){
        return e1.weight - e2.weight;
    });
    
    this.mst = [];
    
    this.visit(G, 0);
    
    while(!this.pq.isEmpty() && this.mst.length < V-1) {
        var e = this.pq.delMin();
        var v = e.either();
        var w = e.other(v);
        if(this.marked[v] && this.marked[w]) continue;
        this.mst.push(e);
        if(!this.marked[v]) this.visit(G, v);
        if(!this.marked[w]) this.visit(G, w);
    }
};

LazyPrimMST.prototype.visit = function(G, v) {
    this.marked[v]  = true;
    var adj_v = G.adj(v);
    for (var i = 0; i < adj_v.length; ++i) {
        var e = adj_v[i];
        if(!this.marked[e.other(v)]){
            this.pq.enqueue(e);
        }
    }
};


var EagerPrimMST = function(G) {
    var V = G.V;
    this.pq = new IndexMinPQ(V, function(e1, e2) {
        return e1.weight - e2.weight;
    });
    this.marked = [];
    for(var v = 0; v < V; ++v) {
        this.marked.push(false);
    }
    this.mst = [];
    this.visit(G, 0);
    while(!this.pq.isEmpty()) {
        var e = this.pq.minKey();
        var w = this.pq.delMin();
        
        this.mst.push(e);
        
        if(!this.marked[w]){
            this.visit(G, w);
        }
        
    }
};

EagerPrimMST.prototype.visit = function(G, v) {
    this.marked[v]  = true;
    var adj_v = G.adj(v);
    for(var i = 0; i < adj_v.length; ++i) {
        var e = adj_v[i];
        var w = e.other(v);
        if(this.marked[w]) continue;
        if(this.pq.containsIndex(w)){
            this.pq.decreaseKey(w, e);
        } else {
            this.pq.insert(w, e);
        }
    }
};


var Dijkstra = function(G, s) {
    var V = G.V;
    this.s = s;
    this.marked = [];
    this.edgeTo = [];
    this.cost = [];
    this.pq = new IndexMinPQ(V, function(cost1, cost2){
        return cost1, cost2;
    });
    
    for(var v =0; v < V; ++v){
        this.marked.push(false);
        this.edgeTo.push(null);
        this.cost.push(Number.MAX_VALUE);
    }
    
    this.cost[s] = 0;
    
    this.pq.insert(s, this.cost[s]);
    
    while(!this.pq.isEmpty()) {
        var v = this.pq.delMin();
        this.marked[v] = true;
        var adj_v = G.adj(v);
        for(var i = 0; i < adj_v.length; ++i) {
            var e = adj_v[i];
            this.relax(e);
        }
    }
    
};

    


Dijkstra.prototype.relax = function(e) {
    
    var v = e.from();
    var w = e.to();
    
    if(this.cost[w] > this.cost[v] + e.weight) {
        this.cost[w] = this.cost[v] + e.weight;
        this.edgeTo[w] = e;
        if(this.pq.containsIndex(w)){
            this.pq.decreaseKey(w, this.cost[w]);
        } else {
            this.pq.insert(w, this.cost[w]);
        }
    }
};



Dijkstra.prototype.hasPathTo = function(v) {
    return this.marked[v];  
};


Dijkstra.prototype.pathTo = function(v) {
    var path = new Stack();
    for(var x = v; x != this.s; x = this.edgeTo[x].other(x)) {
        path.push(this.edgeTo[x]);
    }  
    return path.toArray();
};

Dijkstra.prototype.distanceTo = function(v) {
    return this.cost[v];  
};


var BellmanFord = function(G, s) {
    var V = G.V;
    this.s = s;
    this.marked = [];
    this.edgeTo = [];
    this.cost = [];
    
    
    for(var v =0; v < V; ++v){
        this.marked.push(false);
        this.edgeTo.push(null);
        this.cost.push(Number.MAX_VALUE);
    }
    
    this.cost[s] = 0;
    this.marked[s] = true;
    
    for(var j = 0; j < V; ++j) {
        for(var v = 0; v < V; ++v){
            var adj_v = G.adj(v);
            for(var i = 0; i < adj_v.length; ++i) {
                var e = adj_v[i];
                this.relax(e);
            }
        }
    }
    
};

BellmanFord.prototype.relax = function(e) {
    
    var v = e.from();
    var w = e.to();
    
    if(this.cost[w] > this.cost[v] + e.weight) {
        this.cost[w] = this.cost[v] + e.weight;
        this.marked[w] = true;
        this.edgeTo[w] = e;
    }
};

BellmanFord.prototype.hasPathTo = function(v) {
    return this.marked[v];  
};


BellmanFord.prototype.pathTo = function(v) {
    var path = new Stack();
    for(var x = v; x != this.s; x = this.edgeTo[x].other(x)) {
        path.push(this.edgeTo[x]);
    }  
    return path.toArray();
};

BellmanFord.prototype.distanceTo = function(v) {
    return this.cost[v];  
};


var TopologicalSortShortestPaths = function(G, s) {
    var V = G.V;
    this.s = s;
    this.marked = [];
    this.edgeTo = [];
    this.cost = [];
    
    
    for(var v =0; v < V; ++v){
        this.marked.push(false);
        this.edgeTo.push(null);
        this.cost.push(Number.MAX_VALUE);
    }
    
    this.cost[s] = 0;
    this.marked[s] = true;
    
    var order = new TopologicalSort(G.toDiGraph()).order();
    
    
    for(var j = 0; j < order.length; ++j){
        var v = order[j];
        var adj_v = G.adj(v);
        for(var i = 0; i < adj_v.length; ++i) {
            var e = adj_v[i];
            this.relax(e);
        }
    }
    
    
};

TopologicalSortShortestPaths.prototype.relax = function(e) {
    
    var v = e.from();
    var w = e.to();
    
    if(this.cost[w] > this.cost[v] + e.weight) {
        this.cost[w] = this.cost[v] + e.weight;
        this.marked[w] = true;
        this.edgeTo[w] = e;
    }
};

TopologicalSortShortestPaths.prototype.hasPathTo = function(v) {
    return this.marked[v];  
};


TopologicalSortShortestPaths.prototype.pathTo = function(v) {
    var path = new Stack();
    for(var x = v; x != this.s; x = this.edgeTo[x].other(x)) {
        path.push(this.edgeTo[x]);
    }  
    return path.toArray();
};

TopologicalSortShortestPaths.prototype.distanceTo = function(v) {
    return this.cost[v];  
};


var FordFulkerson = function(G, s, t) {
    this.value = 0;
    var V = G.V;
    var bottle = Number.MAX_VALUE;
    this.marked = null;
    this.edgeTo = null;
    this.s = s;
    this.t = t;
    while(this.hasAugmentedPath(G)){
        
        for(var x = this.t; x != this.s; x = this.edgeTo[x].other(x)) {
            bottle = Math.min(bottle, this.edgeTo[x].residualCapacityTo(x));
        }
        
        for(var x = this.t; x != this.s; x = this.edgeTo[x].other(x)) {
            this.edgeTo[x].addResidualFlowTo(x, bottle);
        }
        
        
        this.value += bottle;
    }
};

FordFulkerson.prototype.hasAugmentedPath = function(G) {
    var V = G.V;
    this.marked = [];
    this.edgeTo = [];
    for(var v = 0; v < V; ++v) {
        this.marked.push(false);
        this.edgeTo.push(null);
    }
    
    var queue = new Queue();
    queue.enqueue(this.s);
    
    this.marked[this.s] = true;
    while(!queue.isEmpty()){
        var v = queue.dequeue();
        var adj_v = G.adj(v);
        
        for (var i = 0; i < adj_v.length; ++i) {
            var e = adj_v[i];
            var w = e.other(v);
            if(!this.marked[w] && e.residualCapacityTo(w) > 0){
                this.edgeTo[w] = e;
                this.marked[w] = true;
                if(w == this.t){
                    return true;
                }
                
                queue.enqueue(w);
            }
        }
    }
    
    return false;
};

FordFulkerson.prototype.minCut = function(G) {
    var cuts = [];
    var V = G.V;
    for(var v = 0; v < V; ++v){
        var adj_v = G.adj(v);
        for(var i = 0; i < adj_v.length; ++i) {
            var e = adj_v[i];
            if(e.from() == v && e.residualCapacityTo(e.other(v)) == 0) {
                cuts.push(e);
            }
        }
    }
    
    return cuts;
};

export { StackNode };
export { Stack };
export { QueueNode };
export { Queue };
export { MinPQ };
export { QuickUnion };
export { IndexMinPQ };
export { Graph };
export { DiGraph };
export { Edge };
export { WeightedGraph };
export { WeightedDiGraph };
export { FlowEdge };
export { FlowNetwork };
export { DepthFirstSearch };
export { BreadthFirstSearch };
export { ConnectedComponents };
export { TopologicalSort };
export { StronglyConnectedComponents };
export { KruskalMST };
export { LazyPrimMST };
export { EagerPrimMST };
export { Dijkstra };
export { BellmanFord };
export { TopologicalSortShortestPaths };
export { FordFulkerson };
