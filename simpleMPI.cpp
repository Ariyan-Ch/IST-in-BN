// COMPILE WITH: mpicxx -O2 parallel_mpi.cpp -o parallel_mpi
// RUN WITH: mpirun --hostfile machines -np <num_procs> ./parallel_mpi N

#include <mpi.h>
#include <bits/stdc++.h>
using namespace std;

// --- Type aliases ----------------------------------------------
using Perm = vector<int>;
std::vector<long long> fact;  // Global factorial cache
// --- Helpers ------------------------------------------------

// Given an index 0 ≤ idx < N! and size N, produce the idx-th
// permutation in lex order in O(N²) time and O(N) space.
Perm unrank(long long idx, int N) {
    std::vector<int> elems(N);
    std::iota(elems.begin(), elems.end(), 1);
    Perm perm;

    for (int i = 0; i < N; i++) {
        long long f = fact[N - 1 - i];
        int pos = idx / f;
        idx %= f;
        perm.push_back(elems[pos]);
        elems.erase(elems.begin() + pos);
    }
    return perm;
}


// And to get the lex‐rank of any perm:
int rank_of(const Perm &perm) {
    int N = perm.size();
    std::vector<int> elems(N);
    std::iota(elems.begin(), elems.end(), 1);
    long long idx = 0;

    for (int i = 0; i < N; i++) {
        auto it = std::find(elems.begin(), elems.end(), perm[i]);
        int slot = std::distance(elems.begin(), it);
        idx += slot * fact[N - 1 - i];
        elems.erase(it);
    }
    return idx;
}

// Convert [1,2,10,3] -> "1-2-10-3"
string perm_to_str(const Perm &v) {
    ostringstream oss;
    oss << v[0];
    for (int i = 1; i < (int)v.size(); i++)
        oss << '-' << v[i];
    return oss.str();
}

// Check if v is the identity 1,2,...,n
bool is_identity(const Perm &v) {
    for (int i = 0; i < (int)v.size(); i++)
        if (v[i] != i+1) return false;
    return true;
}

// Swap the symbol x with its right neighbor
Perm swap_symbol(const Perm &v, int x, const vector<int> &inv) {
    int i = inv[x];        // 0-based index of x
    Perm p = v;
    std::swap(p[i], p[i+1]);
    return p;
}

// FindPosition subroutine (rules 1.1-1.3)
Perm find_position(const Perm &v,
                   int t, int n,
                   const vector<int> &inv,
                   int rpos) {
    if (t == 2) {
        Perm cand = swap_symbol(v, 2, inv);
        if (is_identity(cand))
            return swap_symbol(v, 1, inv);
    }
    if (v[n-2] == t || v[n-2] == n-1) {
        int j = rpos + 1;
        return swap_symbol(v, j, inv);
    }
    return swap_symbol(v, t, inv);
}

// Parent-finding routine (Algorithm 1)
Perm get_parent(const Perm &v, int t, int n) {
    vector<int> inv(n+1);
    for (int i = 0; i < n; i++) inv[v[i]] = i;
    int rpos = -1;
    for (int i = n-1; i >= 0; i--) {
        if (v[i] != i+1) { rpos = i; break; }
    }

    if (v[n-1] == n) {
        if (t != n-1) {
            return find_position(v, t, n, inv, rpos);
        } else {
            return swap_symbol(v, v[n-2], inv);
        }
    }
    if (v[n-1] == n-1 && v[n-2] == n) {
        Perm s = swap_symbol(v, n, inv);
        if (!is_identity(s)) {
            if (t == 1) return s;
            else        return swap_symbol(v, t-1, inv);
        }
    }
    if (v[n-1] == t) {
        return swap_symbol(v, n, inv);
    } else {
        return swap_symbol(v, t, inv);
    }
}

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc != 2) {
        if (rank == 0)
            cerr << "Usage: " << argv[0] << " N\n";
        MPI_Finalize();
        return 1;
    }
    int n = atoi(argv[1]);
    if (n < 2) {
        if (rank == 0)
            cerr << "Error: N must be ≥ 2.\n";
        MPI_Finalize();
        return 1;
    }

    // 1) Precompute factorials
    fact.resize(n + 1);
    fact[0] = 1;
    for (int i = 1; i <= n; ++i)
        fact[i] = fact[i - 1] * i;

    const long long Nperm = fact[n];
    const int trees = n - 1;

    // 2) Determine this rank's range
    long long base = Nperm / size;
    long long rem = Nperm % size;
    long long start = rank * base + std::min<long long>(rank, rem);
    long long count = base + (rank < rem ? 1 : 0);

    // 3) Start compute timer
    double t0 = MPI_Wtime();

    // 4) Open per-rank output file
    std::ofstream ofs("/tmp/results/Parallel_parents_rank" + std::to_string(rank) + ".csv");
    ofs << "vertex_id,perm";
    for (int t = 1; t <= trees; ++t)
        ofs << ",T" << t;
    ofs << "\n";

    // 5) Compute and write
    for (long long vid = start; vid < start + count; ++vid) {
        Perm v = unrank(vid, n);
        ofs << vid << "," << perm_to_str(v);
        if (!is_identity(v)) {
            for (int t = 1; t <= trees; ++t) {
                Perm par = get_parent(v, t, n);
                int pid = rank_of(par);
                ofs << "," << pid;
            }
        } else {
            for (int t = 1; t <= trees; ++t)
                ofs << ",-1";
        }
        ofs << "\n";
    }
    ofs.close();

    // 6) Report timing
    double t1 = MPI_Wtime();
    double my_time = t1 - t0;
    if (rank == 0)
        std::cout << "Each rank wrote its own output file.\n";
    std::cout << "Rank " << rank << " finished in " << my_time << " seconds.\n";

    MPI_Finalize();
    
  
    
    return 0;
}

