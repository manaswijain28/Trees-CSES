#include <bits/stdc++.h>
using namespace std;
#define ll long long


void solve(){


    int n;
    cin>>n;


    set<int> distinctColor[n+1];
    vector<vector<int>> adj(n+1);
    vector<int> answer(n+1);



    int a[n+1];

    for(int i=1;i<=n;i++){
        cin>>a[i];

    }


    for(int i=1;i<n;i++){
        int x,y;
        cin>>x>>y;


        adj[x].push_back(y);
        adj[y].push_back(x);

    }


    function<void(int,int)> dfs = [&](int node,int parent)-> void{
        distinctColor[node].insert(a[node]);

        for(auto it : adj[node]){
            if(it != parent){
                dfs(it,node);

                if(distinctColor[it].size()>distinctColor[node].size()){
                    swap(distinctColor[it],distinctColor[node]);
                }


                for(auto i : distinctColor[it]){
                    distinctColor[node].insert(i);
                }
            }
        }


        answer[node] = distinctColor[node].size();
    };


    dfs(1,0);

    for(int i=1;i<=n;i++){
        cout<<answer[i]<<" ";
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
