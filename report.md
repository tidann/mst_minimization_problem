# Programming Project: Minimum Spanning Tree

## Introduction

We are given a connected graph $G$ with $N$ vertices and $M$ edges, and a spanning tree $T$ of $G$. The task is to modify edge weights (to new integer values) so that $T$ becomes a minimum spanning tree of the modified graph, while minimizing the total change: $\sum_{e \in E} |W'(e) - W(e)|$, where $W$ and $W'$ denote original and new weights respectively.

This objective is equivalent to **$L_1$ isotonic regression** on a directed acyclic graph of constraints derived from the MST cycle property. Stout [1] (*$L_p$ Isotonic Regression Algorithms Using an $L_0$ Approach*) shows that such $L_1$ regression with integer values can be solved by **partitioning**: at each step we solve a binary (threshold) subproblem, which reduces to a minimum cut in a flow network. Our method implements this idea using divide-and-conquer over the set of distinct weights and Dinic’s max-flow algorithm for each binary subproblem. Correctness is argued in the *Validity of the Algorithm* section below.

## Problem Analysis

### Cycle property and constraints

From the course (Section 5.3), a spanning tree $T$ is an MST if and only if it satisfies the **cycle property**: for every non-tree edge $f = (u,v)$, the weight of $f$ is at least the weight of every tree edge on the unique cycle in $T \cup \{f\}$. Denoting the new weight of edge $e$ by $W'(e)$, we thus require
$$ W'(e) \le W'(f) \qquad \text{for every non-tree edge } f \text{ and every tree edge } e \text{ on the cycle of } f. $$

Let $x_e = W'(e)$ and $w_e = W(e)$. We must choose integer $x_e$ so that the constraints $x_e \le x_f$ hold for the relevant pairs $(e,f)$, and we minimize $\sum_e |x_e - w_e|$.

### Formulation as $L_1$ isotonic regression

The constraints form a partial order on the set of edges (a DAG): an arc $e \to f$ means $x_e \le x_f$. Minimizing $\sum |x_e - w_e|$ under these constraints is **$L_1$ isotonic regression** on that DAG. It is well known (and used in [1]) that an optimal solution can be chosen so that every $x_e$ lies in the set of original weights. Let $v_1 < v_2 < \cdots < v_k$ be the distinct values among $\{w_e\}$.

We solve the problem by **divide and conquer** on the index range $[L,R]$ of these values, deciding at each step which variables are $\le v_{mid}$ and which are $> v_{mid}$ via a **minimum cut** in an auxiliary flow network.

### Algorithm

**Procedure** `Solve(nodes $U$, range $[L,R]$)`:

1. **Base case.** If $L = R$, assign $x_i = v_L$ for every $i \in U$; return.
2. **Pivot.** Set $mid = \lfloor (L+R)/2 \rfloor$ and $val = v_{mid}$.
3. **Flow network.** Build a graph with source $S$, sink $T$, and one node per element of $U$.
   - **Cost edges:** For each $i \in U$: if $w_i > val$, add $S \to i$ with capacity $1$ (putting $i$ on the $T$-side costs $1$); if $w_i \le val$, add $i \to T$ with capacity $1$ (putting $i$ on the $S$-side costs $1$).
   - **Constraint edges:** For each constraint $x_u \le x_v$ with $u,v \in U$, add $u \to v$ with capacity $\infty$ so the cut cannot separate $u$ (on $S$) from $v$ (on $T$).
4. **Min cut.** Run a max-flow algorithm (e.g. Dinic) from $S$ to $T$. In the residual graph, let $U_{>}$ be the set of $i \in U$ reachable from $S$, and $U_{\le} = U \setminus U_{>}$. Then variables in $U_{>}$ will be assigned values $> val$, and those in $U_{\le}$ values $\le val$.
5. **Recurse.** Call `Solve(U_{\le}, [L, mid])` and `Solve(U_{>}, [mid+1, R])`.

### Validity of the Algorithm

We minimize $\Phi(x) = \sum_{i \in V} |x_i - w_i|$ subject to $x_u \le x_v$ for all constraint edges, with $x_i \in S = \{v_1, \dots, v_k\}$. This is **$L_1$ isotonic regression** on a DAG with integer values (unit weights, i.e. weight 1 per term).

**Theoretical framework [1].** Stout’s paper shows that $L_1$ isotonic regression with $\ell$ distinct values can be solved by **partitioning**: repeatedly solving binary (two-label) $L_1$ subproblems, each equivalent to a minimum cut / maximum flow problem.

- **Theorem 1 (ii)** of [1]: On a DAG $G$ with integer weighted function $(f,w)$ and $\ell$ distinct values, an exact $L_1$ isotonic regression can be found in time $O(\mathcal{V}(G) + \mathcal{F}_2(\widehat{n},\widehat{m},U) \cdot \log \ell)$. **Notation:** $\mathcal{V}(G)$ = time to build the *violator dag* of $G$ (a directed graph encoding all violating pairs $u \prec v$ with $f(u) > f(v)$); $\mathcal{F}_2(\widehat{n},\widehat{m},U)$ = time of the max-flow / min-cost flow algorithm on the *binary* subproblem graph, which has $\widehat{n}$ vertices, $\widehat{m}$ edges, and maximum capacity/weight $U$; $\ell$ = number of distinct function values. So the total time is “build violator structure” plus $\log \ell$ times “solve one binary flow problem”.
- **Section 4.1** of [1] (partitioning for $L_1$): For two values $a < b$ with no data in $(a,b)$, define binary $g$ by $g(v)=a$ if $f(v)\le a$, $g(v)=b$ if $f(v)\ge b$. There is an optimal $L_1$ isotonic regression of $g$ with values in $\{a,b\}$ that partitions vertices into $V_a$ and $V_b$; no vertex in $V_b$ precedes one in $V_a$. Regressions on the subgraphs induced by $V_a$ and $V_b$ can be computed **independently**. By choosing $a,b$ as median (or quartiles) and iterating, at most $O(\log \ell)$ steps are needed for $\ell$ distinct values.

Our algorithm follows this scheme: at each step we choose a pivot $val = v_{mid}$, solve the binary subproblem “$\le val$” vs “$> val$”, partition the edges by the cut, then recurse on the two subsets.

**Binary subproblem and min-cut.** For a fixed pivot $val$, we decide for each variable $x_i$ whether $x_i \le val$ or $x_i > val$, minimizing the resulting $L_1$ cost under $x_u \le x_v$. This is the **binary** $L_1$ isotonic regression (equivalent to $L_0$ with two labels) in [1]. Define $y_i = 1$ iff $x_i > val$, $y_i = 0$ iff $x_i \le val$. The cost for variable $i$ is $|y_i - C_i|$ with $C_i = 1$ if $w_i > val$, $C_i = 0$ otherwise. The constraint $x_u \le x_v$ implies $y_u \le y_v$. The flow network we build (source $S$, sink $T$, capacity 1 edges $S \to i$ or $i \to T$ according to $w_i$ vs $val$, and $\infty$ capacity $u \to v$ for each constraint) is exactly the standard reduction of this binary problem to **minimum cut**. By the **max-flow min-cut theorem** (course), a minimum cut gives the optimal binary assignment; nodes reachable from $S$ in the residual graph are “$> val$”, the others “$\le val$”.

**Correctness.** (1) By Section 4.1 of [1], the binary subproblem at pivot $val$ yields a partition consistent with some global optimal $L_1$ isotonic regression (since there exists an optimal solution with values in $\{v_1,\dots,v_k\}$). (2) Constraints between the two sides are respected (we never have $x_u > val$ and $x_v \le val$), so the subproblems on $U_{\le}$ and $U_{>}$ are independent. (3) By Theorem 1 (ii) and Section 4.1 of [1], this partitioning procedure produces an **exact** $L_1$ isotonic regression. Our implementation solves each binary subproblem by max-flow min-cut; therefore the final assignment is optimal for the MST weight modification problem.

### Complexity Analysis

We use the course’s complexity bounds for Dinic’s max-flow algorithm (Section 7). The overall cost comes from (i) building constraints, (ii) recursion depth, and (iii) the cost of each max-flow.

- **Constraints.** For each of the $O(M)$ non-tree edges we find the fundamental cycle in $T$ in $O(N)$ time; the cycle has length $O(N)$. So we obtain $O(MN)$ constraints. Each flow instance has $O(M)$ nodes and $O(M) + O(MN) = O(MN)$ edges.

- **Recursion.** We do a binary search over the $k \le M$ distinct weights, so the depth is $O(\log M)$. At each level the active sets form a partition of the edge set, so the total “size” of all subproblems at one level is $O(M)$ nodes and $O(MN)$ edges. Because the cost of Dinic is super-linear in the graph size, the total work per level is at most the cost of one max-flow on a graph of $O(M)$ nodes and $O(MN)$ edges (convexity of $T(V,E)$).

- **Max-flow.** From the course: Dinic runs in $O(V^2 E)$ in general, and in $O(E\sqrt{V})$ for unit-capacity networks. Our source/sink edges have capacity 1; the constraint edges have infinite capacity but do not increase the value of the maximum flow (bounded by $M$). Using the unit-capacity bound, one max-flow is $O(MN\sqrt{M}) = O(M^{1.5} N)$.

**Total:** $O(\log M \cdot M^{1.5} N)$. For $M \le 800$, $N \le 50$ this is on the order of $10^7$ operations, which matches the observed run times (all tests under 0.03 s).

## Implementation and Results

The solver is implemented in C++17. We build the constraint DAG from the given tree and non-tree edges, collect distinct weights, then run the divide-and-conquer procedure with Dinic’s algorithm for each min-cut. Output is the optimal total cost $\sum_e |x_e - w_e|$.

**Test results** (10 provided cases):

| Test   | Cost | Time (s) |
|--------|------|----------|
| 0.in   | 12   | 0.004    |
| 1.in   | 0    | 0.003    |
| 2.in   | 712  | 0.007    |
| 3.in   | 852  | 0.011    |
| 4.in   | 930  | 0.023    |
| 5.in   | 1484 | 0.020    |
| 6.in   | 1244 | 0.021    |
| 7.in   | 367  | 0.024    |
| 8.in   | 1005 | 0.023    |
| 9.in   | 1214 | 0.023    |

All runs complete in under 0.03 seconds.

**Build and run:** from the project root, `mkdir build && cd build`, then `cmake ..`, `make`, and `./mst_solver < ../tests/0.in` (or the path to any test input).

## Bibliography

[1] Q. F. Stout, *$L_p$ Isotonic Regression Algorithms Using an $L_0$ Approach*, arXiv:2107.00251, 2021.  
https://arxiv.org/pdf/2107.00251
