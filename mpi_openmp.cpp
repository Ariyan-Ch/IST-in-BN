// COMPILE WITH: mpicxx -O2 -fopenmp parallel_mpi.cpp -o parallel_mpi
// RUN WITH:     mpirun --hostfile machines -np <num_procs> ./parallel_mpi N

#include <mpi.h>
#include <omp.h>
#include <bits/stdc++.h>
using namespace std;

// --- Type aliases & globals -----------------------------------
using Perm = vector<int>;
static vector<long long> fact;

// --- Helpers (unchanged) --------------------------------------
Perm unrank(long long idx, int N) {
    vector<int> elems(N);
    iota(elems.begin(), elems.end(), 1);
    Perm perm; perm.reserve(N);
    for (int i = 0; i < N; i++) {
        long long f = fact[N - 1 - i];
        int pos = idx / f;
        idx %= f;
        perm.push_back(elems[pos]);
        elems.erase(elems.begin() + pos);
    }
    return perm;
}
int rank_of(const Perm &perm) {
    int N = perm.size();
    vector<int> elems(N);
    iota(elems.begin(), elems.end(), 1);
    long long idx = 0;
    for (int i = 0; i < N; i++) {
        auto it = find(elems.begin(), elems.end(), perm[i]);
        int slot = distance(elems.begin(), it);
        idx += slot * fact[N - 1 - i];
        elems.erase(it);
    }
    return static_cast<int>(idx);
}
string perm_to_str(const Perm &v) {
    ostringstream oss;
    oss << v[0];
    for (int i = 1; i < (int)v.size(); i++)
        oss << '-' << v[i];
    return oss.str();
}
bool is_identity(const Perm &v) {
    for (int i = 0; i < (int)v.size(); i++)
        if (v[i] != i+1) return false;
    return true;
}
Perm swap_symbol(const Perm &v, int x, const vector<int> &inv) {
    int i = inv[x];
    Perm p = v;
    swap(p[i], p[i+1]);
    return p;
}
Perm find_position(const Perm &v, int t, int n,
                   const vector<int> &inv, int rpos) {
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
Perm get_parent(const Perm &v, int t, int n) {
    vector<int> inv(n+1);
    for (int i = 0; i < n; i++) inv[v[i]] = i;
    int rpos = -1;
    for (int i = n-1; i >= 0; i--)
        if (v[i] != i+1) { rpos = i; break; }

    if (v[n-1] == n) {
        if (t != n-1)
            return find_position(v, t, n, inv, rpos);
        else
            return swap_symbol(v, v[n-2], inv);
    }
    if (v[n-1] == n-1 && v[n-2] == n) {
        Perm s = swap_symbol(v, n, inv);
        if (!is_identity(s))
            return (t == 1 ? s : swap_symbol(v, t-1, inv));
    }
    if (v[n-1] == t)
        return swap_symbol(v, n, inv);
    else
        return swap_symbol(v, t, inv);
}

// --- main ------------------------------------------------------
int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc != 2) {
        if (rank==0) cerr << "Usage: " << argv[0] << " N\n";
        MPI_Finalize();
        return 1;
    }
    int n = atoi(argv[1]);
    if (n < 2) {
        if (rank==0) cerr << "Error: N must be ≥ 2.\n";
        MPI_Finalize();
        return 1;
    }

    // 1) Factorials
    fact.resize(n+1);
    fact[0] = 1;
    for (int i = 1; i <= n; ++i)
        fact[i] = fact[i-1] * i;
    const long long Nperm = fact[n];
    const int trees = n - 1;

    // 2) This rank’s chunk of the factorial space
    long long base = Nperm / size;
    long long rem  = Nperm % size;
    long long start = rank*base + min<long long>(rank, rem);
    long long count = base + (rank < rem ? 1 : 0);

    // 3) Timer start
    double t0 = MPI_Wtime();

    // 4) Open file once
    ofstream ofs("/tmp/results/Parallel_parents_rank" + to_string(rank) + ".csv");
    ofs << "vertex_id,perm";
    for (int t = 1; t <= trees; ++t) ofs << ",T" << t;
    ofs << "\n";

    // 5) Compute & buffer in one big OpenMP region
    vector<vector<string>> thread_lines;
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        int nthreads = omp_get_num_threads();

        #pragma omp single
        thread_lines.resize(nthreads);

        // static chunk for this thread
        long long b = count / nthreads;
        long long r = count % nthreads;
        long long st = start + tid*b + min<long long>(tid, r);
        long long cnt = b + (tid < r ? 1 : 0);

        auto &buf = thread_lines[tid];
        buf.reserve(cnt);

        // do all my vids
        for (long long vid = st; vid < st + cnt; ++vid) {
            Perm v = unrank(vid, n);

            // compute inv & rpos once
            vector<int> inv(n+1);
            for (int i = 0; i < n; ++i) inv[v[i]] = i;
            int rpos = -1;
            for (int i = n-1; i >= 0; --i)
                if (v[i] != i+1) { rpos = i; break; }

            // get all parents
            vector<int> pids(trees);
            for (int t = 1; t <= trees; ++t) {
                Perm par = get_parent(v, t, n);
                pids[t-1] = rank_of(par);
            }

            // build CSV line
            ostringstream oss;
            oss << vid << "," << perm_to_str(v);
            for (int x : pids) oss << "," << x;
            oss << "\n";
            buf.push_back(oss.str());
        }

        // barrier then one‐time file write
        #pragma omp barrier
        #pragma omp single
        {
            for (auto &thbuf : thread_lines)
                for (auto &line : thbuf)
                    ofs << line;
        }
    }

    ofs.close();

    // 6) Timing report
    double t1 = MPI_Wtime();
    if (rank == 0) cout << "Each rank wrote its own output file.\n";
    cout << "Rank " << rank << " finished in " << (t1 - t0) << " seconds.\n";

    MPI_Finalize();
    return 0;
}

