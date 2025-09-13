#include <bits/stdc++.h>
using namespace std;
#define ll long long
 
 
void solve() {
    int n;
    cin>>n;
 
 
 
    vector<int> adj[n+1];
    for(int i=1;i<n;i++){
        int a,b;
        cin>>a>>b;
 
        adj[a].push_back(b);
        adj[b].push_back(a);
 
    }
 
 
    vector<int> distance1(n+1,0);
    vector<int> distance2(n+1,0);
 
    vector<int> depth(n+1,-1);
 
    function<void(int,int,int)> dfs = [&](int node,int parent,int dist)->void{
        distance1[node] = dist;
        for(auto &it : adj[node]){
            if(it != parent){
                dfs(it,node,dist+1);
            }
        }
    };
 
 
 
    dfs(1,-1,0);
 
    int firstEnd = 0;
    for(int i=1;i<=n;i++){
        if(distance1[firstEnd]<distance1[i]){
            firstEnd = i;
        }
    }
 
    dfs(firstEnd,-1,0);
 
    for(int i=1;i<=n;i++){
        distance2[i] = distance1[i];
    }
 
    int secondEnd = 0;
    for(int i=1;i<=n;i++){
        if(distance1[secondEnd]<distance1[i]){
            secondEnd = i;
        }
    }
 
    dfs(secondEnd,-1,0);
 
 
    for(int i=1;i<=n;i++){
        int answer = max(distance1[i],distance2[i]);
        
        cout<<answer<<" ";
    }
    cout<<endl;
 
 
 
 
 
    
    
 
}
 
int main(){
    int t=1;
    // cin>>t;
    while(t--){
        solve();
    }
}
