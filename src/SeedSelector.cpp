#include "SeedSelector.h"
#include "utils.h"
#include <vector>
#include <queue>
#include <unordered_set>
#include <cmath>
#include <algorithm>

int selectSeed(
    const std::vector<double> &IP,
    const std::vector<std::vector<int>> &adj,
    const std::vector<int> &cac_id)
{
    int n=IP.size();
    int C = *std::max_element(cac_id.begin(),cac_id.end())+1;
    std::vector<std::vector<int>> comms(C);
    for(int i=0;i<n;++i) comms[cac_id[i]].push_back(i);

    int best_u=-1;
    double bestSize=-1, bestDistSum=1e18;

    for(int c=0;c<C;++c){
        for(int u:comms[c]){
            // BFS distances
            std::vector<int> dist(n,-1);
            std::queue<int>q;
            dist[u]=0; q.push(u);
            while(!q.empty()){
                int v=q.front(); q.pop();
                for(int w:adj[v]){
                    if(dist[w]<0){
                        dist[w]=dist[v]+1;
                        q.push(w);
                    }
                }
            }
            int maxD=0;
            for(int i=0;i<n;++i) if(dist[i]>maxD) maxD=dist[i];
            if(maxD<1) continue;

            // IL[L]=avg IP of nodes at exactly distance L
            std::vector<double> IL(maxD+1,0.0);
            std::vector<int> cnt(maxD+1,0);
            for(int i=0;i<n;++i){
                if(dist[i]>0){
                    IL[dist[i]]+=IP[i];
                    cnt[dist[i]]++;
                }
            }
            for(int L=1;L<=maxD;++L){
                if(cnt[L]>0) IL[L]/=cnt[L];
            }

            int L0=-1;
            for(int L=1;L<maxD;++L){
                if(IL[L]>IL[L+1]){L0=L; break;}
            }
            if(L0<1) continue;
            if(IP[u] <= IL[L0]) continue;

            // Influence-BFS: only traverse black nodes
            std::unordered_set<int> vis;
            std::queue<int> qb;
            vis.insert(u); qb.push(u);
            while(!qb.empty()){
                int v=qb.front(); qb.pop();
                for(int w:adj[v]){
                    if(cac_id[w]==c && IP[w]>IL[L0] && vis.insert(w).second){
                        qb.push(w);
                    }
                }
            }
            double size = vis.size();
            double distSum=0;
            for(int v:vis) distSum += dist[v];

            if(size>bestSize || (size==bestSize && distSum<bestDistSum)){
                bestSize=size; bestDistSum=distSum; best_u=u;
            }
        }
    }
    return best_u;
}
