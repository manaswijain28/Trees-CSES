#include <bits/stdc++.h>
using namespace std;
typedef int64_t T;
#define ll long long
// const int MM = 1000000007;
const ll p = 1e13;
const int m = 1e6+1;

//***************************************************
const int mod = 998244353;
const int maxn = 1e6+1;
struct MOD {
    int64_t v; explicit operator int64_t() const { return v % mod; }
    MOD() { v = 0; }
    MOD(int64_t _v) {
        v = (-mod < _v && _v < mod) ? _v : _v % mod;
        if (v < 0) v += mod;
    }
    friend bool operator==(const MOD& a, const MOD& b) {
        return a.v == b.v;
    }
    friend bool operator!=(const MOD& a, const MOD& b) {
        return !(a == b);
    }
    friend bool operator<(const MOD& a, const MOD& b) {
        return a.v < b.v;
    }

    MOD& operator+=(const MOD& m) {
        if ((v += m.v) >= mod) v -= mod;
        return *this;
    }
    MOD& operator-=(const MOD& m) {
        if ((v -= m.v) < 0) v += mod;
        return *this;
    }
    MOD& operator*=(const MOD& m) {
        v = v * m.v % mod; return *this;
    }
    MOD& operator/=(const MOD& m) { return (*this) *= inv(m); }
    friend MOD pow(MOD a, int64_t p) {
        MOD ans = 1; assert(p >= 0);
        for (; p; p /= 2, a *= a) if (p & 1) ans *= a;
            return ans;
    }
    friend MOD inv(const MOD& a) {
        assert(a.v != 0);
        return pow(a, mod - 2);
    }

    MOD operator-() const { return MOD(-v); }
    MOD& operator++() { return *this += 1; }
    MOD& operator--() { return *this -= 1; }
    MOD operator++(int32_t) { MOD temp; temp.v = v++; return temp; }
    MOD operator--(int32_t) { MOD temp; temp.v = v--; return temp; }
    friend MOD operator+(MOD a, const MOD& b) { return a += b; }
    friend MOD operator-(MOD a, const MOD& b) { return a -= b; }
    friend MOD operator*(MOD a, const MOD& b) { return a *= b; }
    friend MOD operator/(MOD a, const MOD& b) { return a /= b; }
    friend ostream& operator<<(ostream& os, const MOD& m) {
        os << m.v; return os;
    }
    friend istream& operator>>(istream& is, MOD& m) {
        int64_t x; is >> x;
        m.v = x;
        return is;
    }
    friend void __print(const MOD &x) {
        cerr << x.v;
    }
};


vector<MOD> fct(maxn, 1), invf(maxn, 1);
void calc_fact() {
    for(int i = 1 ; i < maxn ; i++) {
        fct[i] = fct[i - 1] * i;
    }
    invf.back() = MOD(1) / fct.back();
    for(int i = maxn - 1 ; i ; i--)
        invf[i - 1] = i * invf[i];
}

MOD choose(int n, int r) { // choose r elements out of n elements
    if(r > n)   return MOD(0);
    assert(r <= n);
    return fct[n] * invf[r] * invf[n - r];
}

MOD place(int n, int r) { // x1 + x2 ---- xr = n and limit value of xi >= n
    assert(r > 0);
    return choose(n + r - 1, r - 1);
}


//***************************************************

ll expo2(ll a, ll b){
    ll ans = 1;
    while(b){
        if(b&1){
            ans = (ans * a);  
        }
        a = (a*a);
        b >>= 1;  //divide by 2
    }
    return ans;

}

// ll fact(ll n){
//     ll a = 1;
//     for(ll i = 1;i<=n;i++){
//         a *= i;
//     }

//     return a;
// }

// ll De(ll n){
//     if(n==0){
//         return 1;
//     }

//     ll a = fact(n);
//     ll ans = 0;
//     for(ll i=2;i<=n;i++){
//         if(i%2==0){
//             ans += a/fact(i);
//         }else{
//             ans -= a/fact(i);
//         }
//     }


//     return ans;


// }

// const int N = 2e5;

// vector<int> index(N+1);

// bool static cmp(int a,int b){
//     return index[a]<index[b];
// }


//bfs using a priority queue



// void dijkstra(ll n,vector<vector<pair<ll,ll>>> &adj,vector<ll> &dist,ll src){

//     //min heap of pairs
//     priority_queue<pair<ll,ll>,vector<pair<ll,ll>> ,greater<>> pq;

//     dist[src] = 0;
//     pq.push({0,src});
//     while(!pq.empty()){
//         ll currDist = pq.top().first;
//         ll u = pq.top().second;

//         pq.pop();

//         if(currDist > dist[u]){
//             continue;
//         }

//         for(auto &it : adj[u]){
//             ll v = it.first;
//             ll wt = it.second;

//             if(dist[v] > dist[u]+wt){
//                 dist[v] = dist[u]+wt;
//                 pq.push({dist[v],v});
//             }
//         }

//     }

//     //O((n+m)*logn)
// }


// class UnionFind{
// private: vector<int> p,rank,sz;
// public: 

//     //constructor to assign parent of every node as itself
//     UnionFind(int n){
//         rank.assign(n,0);
//         p.assign(n,0);
//         sz.assign(n,0);
//         for(int i=0;i<n;i++){
//             p[i] = i;
//             sz[i] = 1;
//         }


//     }
    
//     //finding the parent of each node
//     int findSet(int i){
//         return (p[i] == i) ? i : p[i] = findSet(p[i]);

//     }
    
//     //if the two nodes are in the same components then they have the same parent
//     bool isSameSet(int i,int j){
//         return findSet(i) == findSet(j);
//     }

//     void unionSet(int i,int j){

//         //if not in same component
//         if(!isSameSet(i,j)){

//             //finding parent
//             int x = findSet(i);
//             int y = findSet(j);

//             // if(rank[x] > rank[y]){
//             //     p[y] = x;
//             // }else{
//             //     p[x] = y;
//             //     if(rank[x] == rank[y]){
//             //         rank[y]++; 
//             //     }
//             // }

//             if(sz[i] > sz[j]){
//                 swap(i,j);
//             }

//             p[i] = j;
//             sz[j] += sz[i];
//         }
//     }


// };

void solve(){

    // cout<<"H"<<endl;

    int n;
    cin>>n;
 
    vector<vector<int>> adj(n+1);
    for(int i=1;i<n;i++){
        int a,b;
        cin>>a>>b;
 
        adj[a].push_back(b);
        adj[b].push_back(a);
    }

    vector<int> dist(n+1,0);

    queue<int> q;
    q.push(1);

    vector<int> visited(n+1,0);

    visited[1] = 1;

    while(!q.empty()){
        int node = q.front();
        q.pop();

        for(auto it : adj[node]){
            if(visited[it] == 0){
                visited[it] = 1;
                dist[it] = dist[node] +1;
                q.push(it);
            }
        }
    }

    ll sum = 0;
    for(int i=0;i<=n;i++){
        sum += dist[i];
    }


    vector<int> subtree(n+1,0);

    function<void(int,int)> func = [&](int node,int parent)->void{
        subtree[node] = 1;
        for(auto it : adj[node]){
            if(it != parent){
                func(it,node);
                subtree[node] += subtree[it];
            }
        }
    };

    func(1,0);


    vector<ll> answer(n+1,0);

    function<void(int,int,ll)> dfs = [&](int node,int parent,ll ans)->void{
        answer[node] = ans;
        for(auto it : adj[node]){
            if(it != parent){
                dfs(it,node,ans+n-2*subtree[it]);
            }
        }


    };

    dfs(1,0,sum);

    for(int i=1;i<=n;i++){
        cout<<answer[i]<<" ";
    }

    cout<<endl;


    

}

int main() {\



    calc_fact();
    int t=1;
    // cin>>t;
    while(t--){
       solve();


   }




}
