#include <bits/stdc++.h>
using namespace std;

// COMPILE WITH: g++ -O2 sequential.cpp -o sequential
// RUN WITH: ./sequential N

// --- Type aliases ----------------------------------------------
using Perm = vector<int>;

// --- Helpers ------------------------------------------------

// Convert [1,2,10,3] -> "1,2,10,3"
string perm_to_str(const Perm &v) {
    ostringstream oss;
    oss << v[0];
    for (int i = 1; i < (int)v.size(); i++)
        oss << ',' << v[i];
    return oss.str();
}

// Check if v is the identity 1,2,...,n
bool is_identity(const Perm &v) {
    for (int i = 0; i < (int)v.size(); i++)
        if (v[i] != i+1) return false;
    return true;
}

// --- Core building-block: swap the symbol x with its right neighbor -------------------
//
// Find position i of symbol x, then swap v[i] <-> v[i+1].
Perm swap_symbol(const Perm &v, int x, const vector<int> &inv) {
    int i = inv[x];        // 0-based index of x
    Perm p = v;
    std::swap(p[i], p[i+1]);
    return p;
}

// The FindPosition(v) subroutine (rules 1.1-1.3)
Perm find_position(const Perm &v,
                   int t, int n,
                   const vector<int> &inv,
                   int rpos)         // rpos = first out-of-place index (0-based)
{
    // 1.1: if t=2 and swap(v,2)==identity
    if (t == 2) {
        Perm cand = swap_symbol(v, 2, inv);
        if (is_identity(cand))
            return swap_symbol(v, 1, inv);    // Swap(v, t-1)
    }
    // 1.2: else if v[n-2] in {t, n-1}
    if (v[n-2] == t || v[n-2] == n-1) {
        // swap at the rightmost out-of-place
        return swap_symbol(v, v[rpos], inv);
    }
    // 1.3: else swap(v, t)
    return swap_symbol(v, t, inv);
}

// The parent-finding routine (Algorithm 1)
Perm get_parent(const Perm &v, int t, int n) {
    // Precompute inverse and rpos once per vertex
    vector<int> inv(n+1);
    for (int i = 0; i < n; i++) inv[v[i]] = i;
    int rpos = -1;
    // rpos = first index from right where v[i] != i+1
    for (int i = n-1; i >= 0; i--) {
        if (v[i] != i+1) { rpos = i; break; }
    }

    //  Rule (1) & (2): if last symbol = n
    if (v[n-1] == n) {
        if (t != n-1) {
            return find_position(v, t, n, inv, rpos);
        } else {
            // Swap(v, v[n-2])
            return swap_symbol(v, v[n-2], inv);
        }
    }
    // Rule (3): if last two symbols are (n,n-1) in that order,
    //            and swapping n wouldn't go straight to the root
    if (v[n-1] == n-1 && v[n-2] == n) {
        Perm s = swap_symbol(v, n, inv);
        if (!is_identity(s)) {
            if (t == 1) return s;         // Swap(v,n)
            else        return swap_symbol(v, t-1, inv);
        }
    }
    // Rule (5)&(6):
    if (v[n-1] == t) {
        return swap_symbol(v, n, inv);    // (5)
    } else {
        return swap_symbol(v, t, inv);    // (6)
    }
}

//  Main sequential driver --------------------------------------------
int main(int argc, char* argv[]) {
    if (argc != 2) {
        cerr << "Usage: " << argv[0] << " N\n" << "  Enumerates all N! vertices of B_n and computes\n"
             << "  each vertex's parent in the n-1 ISTs, sequentially.\n";
        return 1;
    }
    int n = atoi(argv[1]);
    if (n < 2) {
        cerr << "Error: N must be ≥ 2.\n";
        return 1;
    }

    // Enumerate all permutations in lex order
    Perm p(n);
    iota(p.begin(), p.end(), 1);
    vector<Perm> perms;
    do {
        perms.push_back(p);
    } while (next_permutation(p.begin(), p.end()));
    int Nperm = perms.size();

    // Build a map from "perm-string" -> lex-rank ID
    unordered_map<string,int> id_map;
    id_map.reserve(Nperm*2);
    for (int i = 0; i < Nperm; i++) {
        id_map[ perm_to_str(perms[i]) ] = i;
    }

    // Allocate parent table: parent[v_id][t-1] = parent ID in T^n_t
    vector< vector<int> > parent(Nperm, vector<int>(n-1, -1));

    // Compute, vertex by vertex
    for (int vid = 0; vid < Nperm; vid++) {
        const Perm &v = perms[vid];
        // root (identity) has no parent
        if (is_identity(v)) continue;
        // for each spanning tree t=1..n-1
        for (int t = 1; t < n; t++) {
            Perm par = get_parent(v, t, n);
            int pid = id_map[ perm_to_str(par) ];
            parent[vid][t-1] = pid;
        }
    }

    // create CSV for the result
    ofstream ofs("Sequential_parents.csv");
    ofs << "vertex_id,perm";
    for (int t = 1; t < n; t++) ofs << ",T" << t;
    ofs << "\n";

    for (int i = 0; i < Nperm; i++) {
        ofs << i << "," << perm_to_str(perms[i]);
        for (int t = 1; t < n; t++) {
            ofs << "," << parent[i][t-1];
        }
        ofs << "\n";
    }
    ofs.close();

    cout << "Done. Written " << Nperm << " vertices × " << (n-1)
         << " trees -> parents.csv\n";
    return 0;
}


// EXPLANATION on OUTPUT CSV FILE:

/* 

FOR N = 3, this would be the output:

vertex_id,perm,T1,T2
0,1,2,3,-1,-1
1,1,3,2,4,0
2,2,1,3,3,3
3,2,3,1,2,5
4,3,1,2,5,1
5,3,2,1,3,4


 => First of all, T1 and T2 means there are 2 ISTs found for N=3.
 => Vertex_id is simply the ids assigned to all permutations.
 => perm column has the permutation of the sequence i.e 123
 => first row, 123 and T1 = -1 and T2 = -1 means:
        123 permutation has -1 parent in Tree1 and -1 parent in Tree2

    Similarly, second row, 132 permutation has T1=4 and T2=0
        meaning, its parent in Tree1 is vertex 4 (312) and in Tree2 is vertex0(123)

    And so on...
 => It is one representation of the ISTs as each vertex only needs to know its own parent.
 => For future uses of these ISTs, developers can process this CSV file to build the ISTs in graph form

*/