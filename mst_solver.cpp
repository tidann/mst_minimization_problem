#include <iostream>
#include <vector>
#include <algorithm>
#include <map>
#include <queue>
#include <cmath>
#include <chrono>

using namespace std;

const long long INF = 1e18;

// Edge structure
struct Edge {
    int u, v, w, id;
};

// Flow Edge
struct FlowEdge {
    int to;
    long long capacity;
    long long flow;
    int rev; // index of reverse edge
};

// Dinic's Algorithm for Max Flow
class Dinic {
    int n;
    vector<vector<FlowEdge>> adj;
    vector<int> level;
    vector<int> ptr;

public:
    Dinic(int n) : n(n), adj(n), level(n), ptr(n) {}

    void add_edge(int from, int to, long long capacity) {
        FlowEdge forward = {to, capacity, 0, (int)adj[to].size()};
        FlowEdge backward = {from, 0, 0, (int)adj[from].size()};
        adj[from].push_back(forward);
        adj[to].push_back(backward);
    }

    bool bfs(int s, int t) {
        fill(level.begin(), level.end(), -1);
        level[s] = 0;
        queue<int> q;
        q.push(s);
        while (!q.empty()) {
            int v = q.front();
            q.pop();
            for (const auto& edge : adj[v]) {
                if (edge.capacity - edge.flow > 0 && level[edge.to] == -1) {
                    level[edge.to] = level[v] + 1;
                    q.push(edge.to);
                }
            }
        }
        return level[t] != -1;
    }

    long long dfs(int v, int t, long long pushed) {
        if (pushed == 0) return 0;
        if (v == t) return pushed;
        for (int& cid = ptr[v]; cid < adj[v].size(); ++cid) {
            auto& edge = adj[v][cid];
            int tr = edge.to;
            if (level[v] + 1 != level[tr] || edge.capacity - edge.flow == 0) continue;
            long long tr_pushed = dfs(tr, t, min(pushed, edge.capacity - edge.flow));
            if (tr_pushed == 0) continue;
            edge.flow += tr_pushed;
            adj[tr][edge.rev].flow -= tr_pushed;
            return tr_pushed;
        }
        return 0;
    }

    long long max_flow(int s, int t) {
        long long flow = 0;
        while (bfs(s, t)) {
            fill(ptr.begin(), ptr.end(), 0);
            while (long long pushed = dfs(s, t, INF)) {
                flow += pushed;
            }
        }
        return flow;
    }

    // Get nodes reachable from source in residual graph
    vector<bool> get_cut(int s) {
        vector<bool> visited(n, false);
        queue<int> q;
        q.push(s);
        visited[s] = true;
        while (!q.empty()) {
            int v = q.front();
            q.pop();
            for (const auto& edge : adj[v]) {
                if (edge.capacity - edge.flow > 0 && !visited[edge.to]) {
                    visited[edge.to] = true;
                    q.push(edge.to);
                }
            }
        }
        return visited;
    }
};

int N, M;
vector<Edge> all_edges;
vector<int> tree_edge_indices;
vector<int> non_tree_edge_indices;
vector<int> final_weights;
vector<int> distinct_weights;
vector<pair<int, int>> constraints; // u -> v means x_u <= x_v

// Tree adjacency list: node -> list of (neighbor, edge_index)
vector<vector<pair<int, int>>> tree_adj;

// Find path in tree from start to end
bool get_path(int u, int target, int p, vector<int>& path_edges) {
    if (u == target) return true;
    for (auto& edge : tree_adj[u]) {
        int v = edge.first;
        int idx = edge.second;
        if (v != p) {
            path_edges.push_back(idx);
            if (get_path(v, target, u, path_edges)) return true;
            path_edges.pop_back();
        }
    }
    return false;
}

// Solve function using Divide and Conquer
void solve(const vector<int>& nodes, int L_idx, int R_idx) {
    if (nodes.empty()) return;
    
    if (L_idx == R_idx) {
        for (int id : nodes) {
            final_weights[id] = distinct_weights[L_idx];
        }
        return;
    }

    int mid_idx = (L_idx + R_idx) / 2;
    int val = distinct_weights[mid_idx];

    // Build flow network
    // Nodes 0..M-1 correspond to edges in the current active set
    // Source S = num_active, Sink T = num_active + 1
    int num_active = nodes.size();
    map<int, int> node_map;
    for(int i=0; i<num_active; ++i) node_map[nodes[i]] = i;
    
    int S = num_active;
    int T = num_active + 1;
    Dinic dinic(T + 1);

    for (int i = 0; i < num_active; ++i) {
        int id = nodes[i];
        int w = all_edges[id].w;
        
        // Minimize |x - w|
        // If w > val: prefer value > val (S-side). Cost 1 if <= val (T-side).
        // If w <= val: prefer value <= val (T-side). Cost 1 if > val (S-side).
        if (w > val) {
            dinic.add_edge(S, i, 1);
        } else {
            dinic.add_edge(i, T, 1);
        }
    }

    // Add constraints
    // Only need to enforce constraints between two active nodes.
    // Dependencies with inactive nodes (already fixed to left or right) are implicitly satisfied
    // by the nature of the D&C process and DAG structure of constraints.
    for (const auto& p : constraints) {
        int u = p.first;
        int v = p.second;
        
        if (node_map.count(u) && node_map.count(v)) {
            // Constraint: x_u <= x_v.
            // Invalid state: x_u > val (S-side) AND x_v <= val (T-side).
            // Prevent cut separating u (S) and v (T) by adding edge u -> v with INF capacity.
            dinic.add_edge(node_map[u], node_map[v], INF);
        }
    }
    
    dinic.max_flow(S, T);
    vector<bool> cut = dinic.get_cut(S);
    
    vector<int> left_nodes, right_nodes;
    for (int i = 0; i < num_active; ++i) {
        int id = nodes[i];
        if (cut[i]) {
            // Reachable from S -> assigned > val
            right_nodes.push_back(id);
        } else {
            // Not reachable -> assigned <= val
            left_nodes.push_back(id);
        }
    }
    
    solve(left_nodes, L_idx, mid_idx);
    solve(right_nodes, mid_idx + 1, R_idx);
}

int main() {
    if (scanf("%d %d", &N, &M) != 2) return 0;

    for (int i = 0; i < M; ++i) {
        int u, v, w;
        scanf("%d %d %d", &u, &v, &w);
        all_edges.push_back({u, v, w, i});
        distinct_weights.push_back(w);
    }
    
    sort(distinct_weights.begin(), distinct_weights.end());
    distinct_weights.erase(unique(distinct_weights.begin(), distinct_weights.end()), distinct_weights.end());
    
    tree_adj.resize(N + 1);
    
    // Map edges to indices
    map<pair<int, int>, int> edge_map;
    for (int i = 0; i < M; ++i) {
        int u = all_edges[i].u;
        int v = all_edges[i].v;
        if (u > v) swap(u, v);
        edge_map[{u, v}] = i;
    }
    
    // Identify tree edges
    vector<bool> is_tree(M, false);
    for (int i = 0; i < N - 1; ++i) {
        int u, v;
        scanf("%d %d", &u, &v);
        if (u > v) swap(u, v);
        int idx = edge_map[{u, v}];
        is_tree[idx] = true;
        tree_edge_indices.push_back(idx);
        tree_adj[u].push_back({v, idx});
        tree_adj[v].push_back({u, idx});
    }
    
    for (int i = 0; i < M; ++i) {
        if (!is_tree[i]) {
            non_tree_edge_indices.push_back(i);
        }
    }
    
    // Build constraints from Cycle Property
    // For every non-tree edge f=(u,v), and every tree edge e in the cycle of f:
    // weight(e) <= weight(f) => x_e <= x_f
    for (int f_idx : non_tree_edge_indices) {
        int u = all_edges[f_idx].u;
        int v = all_edges[f_idx].v;
        
        vector<int> path_edges;
        get_path(u, v, -1, path_edges);
        
        for (int e_idx : path_edges) {
            constraints.push_back({e_idx, f_idx});
        }
    }
    
    final_weights.resize(M);
    vector<int> initial_nodes(M);
    for (int i = 0; i < M; ++i) initial_nodes[i] = i;
    
    auto start_time = std::chrono::high_resolution_clock::now();
    solve(initial_nodes, 0, (int)distinct_weights.size() - 1);
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> elapsed = end_time - start_time;
    
    long long total_cost = 0;
    for (int i = 0; i < M; ++i) {
        total_cost += abs(final_weights[i] - all_edges[i].w);
    }
    
    printf("%lld\n", total_cost);
    fprintf(stderr, "Time: %.3f ms\n", elapsed.count());
    
    return 0;
}
