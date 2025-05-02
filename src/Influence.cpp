#include "Influence.h"
#include "utils.h"
#include <cmath>
#include <mpi.h>
#include <omp.h>
#include <algorithm>

double action_weight[3]={0.2,0.5,0.3};
double jaccard(const std::vector<int>&A,const std::vector<int>&B){
    int inter=0,uni=0;
    for(size_t i=0;i<A.size();++i){
        if(A[i]||B[i]){uni++; if(A[i]&&B[i]) inter++;}
    }
    return uni? double(inter)/uni:0.0;
}

std::vector<double> computeInfluence(
    CSRGraph &G,
    const std::vector<std::vector<int>>&I,
    const std::vector<int>&levels,
    MPI_Comm comm){
    int n=G.n,rank;
    MPI_Comm_rank(comm,&rank);
    const double d=0.85,tol=1e-6,eps=1e-12;
    std::vector<int>Npy(n,0);
    for(auto&e:G.edges) for(int a=0;a<3;++a) Npy[e.src]+=e.counts[a];
    std::vector<double>psi_out(n,0.0);
    for(auto&e:G.edges){
        double sa=0;for(int a=0;a<3;++a) sa+=action_weight[a]*e.counts[a];
        double C=jaccard(I[e.src],I[e.dst]);
        e.psi= Npy[e.src]? sa*C/Npy[e.src]:0.0;
        psi_out[e.src]+=e.psi;
    }
    std::vector<double>IP(n,1.0/n),IP2(n);
    int maxL=*std::max_element(levels.begin(),levels.end());
    for(int L=0;L<=maxL;++L){
        std::vector<int>verts;
        for(int u=0;u<n;++u) if(levels[u]==L) verts.push_back(u);
        bool done=false;
        while(!done){
            #pragma omp parallel for
            for(int i=0;i<verts.size();++i){
                int u=verts[i];
                int inD=G.rev_xadj[u+1]-G.rev_xadj[u];
                double tele=inD? (1-d)/inD:0.0;
                double sum=0;
                for(int j=G.rev_xadj[u];j<G.rev_xadj[u+1];++j){
                    auto &e=G.rev_edges[j];
                    sum+= e.psi*IP[e.src]/(psi_out[e.src]+eps);
                }
                IP2[u]=tele + d*sum;
            }
            double diff=0;
            for(int u:verts) diff+=fabs(IP2[u]-IP[u]);
            if(diff<tol) done=true;
            else for(int u:verts) IP[u]=IP2[u];
        }
        MPI_Allreduce(MPI_IN_PLACE,IP.data(),n,MPI_DOUBLE,MPI_SUM,comm);
    }
    return IP;
}
