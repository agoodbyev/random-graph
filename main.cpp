#include <iostream>
#include <vector>
#include <algorithm>
#include <ctime>
#include <cmath>
#include <stdlib.h>
#include <iomanip>

using namespace std;

int check(vector<int> v,int a){
    for(unsigned int i=0;i<v.size();i++){
        if(v[i]==a) return 1;
    }
    return 0;
}

double product(vector<double> v, int k=0){
    double p=1;
    for(unsigned int i=k;i<v.size();i++){
        p*=v[i];
    }
    return p;
}

double UB=1, LB=0, upper, lower;
int counter=0,st_time;

class RandomGraph
{
    public:
    int dimension;
    int edgesNumber;
    vector<int> degree;
    double currentReliability;
    vector<double> multipliers;

    class Edge{
        public:
        int in;
        int out;
        double reliability;

        Edge();

        Edge(int a,int b,double r){
            in=a;
            out=b;
            reliability=r;
        }

        ~Edge(){
            in=0;
            out=0;
            reliability=0;
        }

        void write(){
            cout << in << " " << out << " " << reliability;
        }
    };

    vector<Edge>edges;

    RandomGraph(){
        dimension=0;
        edgesNumber=0;
    }

    RandomGraph(const RandomGraph& orig){
        dimension = orig.dimension;
        edgesNumber = orig.edgesNumber;
        currentReliability = orig.currentReliability;
        for(int i=0;i<dimension;i++) degree.push_back(orig.degree[i]);
        for(int i=0;i<edgesNumber;i++) edges.push_back(orig.edges[i]);
    }

    RandomGraph(int n, int m){
        int a,b;
        double r;
        dimension = n;
        edgesNumber = m;
        cout << "Creating graph(" << n << "," << m << ")" << endl;
        cout << "Insert edges by numbers of two connected vertices and reliability: " << endl;
        for(int i=0;i<n;i++) degree.push_back(0);
        for(int i=0;i<m;i++){
            cout << i+1 << ") ";
            cin >> a >> b >> r;
            edges.push_back(Edge(a,b,r));
            degree[a]++;
            degree[b]++;
        }
    }

    RandomGraph(int a,vector<vector<double> >m){
        dimension = a;
        for(int i=0;i<dimension;i++){degree.push_back(0);}
        int k=0;
        for(int i=0;i<dimension;i++){
            for(int j=i+1;j<dimension;j++) {
                if(m[i][j]>=0) {
                    edges.push_back(Edge(i,j,m[i][j]));
                    k++;
                    degree[i]++;
                    degree[j]++;
                }
            }
        }
        edgesNumber=k;
    }

    ~RandomGraph(){
        dimension = 0;
        edgesNumber = 0;
        currentReliability = 0;
        degree.clear();
        edges.clear();
        multipliers.clear();
    }

    void deleteGraph(){
        dimension = 0;
        edgesNumber = 0;
        degree.clear();
        edges.clear();
        multipliers.clear();
    }

    RandomGraph& operator=(const RandomGraph& orig){
        dimension = orig.dimension;
        edgesNumber = orig.edgesNumber;
        for(int i=0;i<dimension;i++) degree.push_back(orig.degree[i]);
        for(int i=0;i<edgesNumber;i++) edges.push_back(orig.edges[i]);
        return *this;
    }

    void printEdgesList(){
        cout << endl << "Edges list:" << endl;
        for(int i=0;i<edgesNumber;i++){
            cout << i+1 << ") ";
            edges[i].write();
            cout << endl;
        }
    }

    void printAdjacencyMatrix(){
        vector<vector<double> >m;
        vector<double>p;
        for(int i=0;i<dimension;i++) p.push_back(0);
        for(int i=0;i<dimension;i++){
            m.push_back(p);
        }
        for(int i=0;i<edgesNumber;i++){
            m[edges[i].in][edges[i].out]=edges[i].reliability;
            m[edges[i].out][edges[i].in]=edges[i].reliability;
        }
        for(int i=0;i<dimension;i++){
            for(int j=0;j<dimension;j++){
                cout << std::setprecision(3) << m[i][j] << " ";
            }
            cout << endl;
        }
    }

    int createConnectedRandomGraph(int v,int e,double p=-1){
        if(e>v*(v-1)/2) {
            cout << "Impossible graph";
            return 0;
        }
        if(v>e+1) {
            cout << "Disconnected graph";
            return 0;
        }
        srand ( time(NULL) );
        dimension = v;
        edgesNumber = e;
        for(int i=0;i<dimension;i++) degree.push_back(0);
        int r, s;
        int flag;
        vector<int> a(v);
        for(int i=0;i<v;i++) a[i]=i;
        swap(a[rand()%v],a[v-1]);
        flag = v-1;
        while(flag!=0){
            r = rand()%flag;
            s = flag+rand()%(v-flag);
            degree[a[r]]++;
            degree[a[s]]++;
            if(p==-1) edges.push_back(Edge(a[r],a[s],(double)(rand())/RAND_MAX/2+0.5));
            else edges.push_back(Edge(a[r],a[s],p));
            swap(a[r],a[flag-1]);
            flag--;
            e--;
        }
        while(e!=0){
            r=rand()%v;
            s=rand()%v;
            flag=0;
            for(unsigned int i=0;i<edges.size();i++){
                if((edges[i].in==r&&edges[i].out==s)||(edges[i].in==s&&edges[i].out==r)){
                    flag=1;
                    break;
                }
            }
            if((r!=s)&&(flag==0)){
                if(p==-1) edges.push_back(Edge(r,s,(double)(rand())/RAND_MAX/2+0.5));
                else edges.push_back(Edge(r,s,p));
                degree[r]++;
                degree[s]++;
                e--;
            }
            flag = 0;
        }
        return 1;
    }

    int createGridGraph(int a, int b,double p=-1){
        srand ( time(NULL) );
        dimension = a*b;
        for(int i=0;i<dimension;i++) degree.push_back(0);
        for(int j=0;j<b;j++){
            for(int i=0;i<a-1;i++){
                if(p==-1) edges.push_back(Edge(i+a*j,i+a*j+1,(double)(rand())/RAND_MAX/2+0.5));
                else edges.push_back(Edge(i+a*j,i+a*j+1,p));
                edgesNumber++;
                degree[i+a*j]++;
                degree[i+a*j+1]++;
            }
        }
        for(int j=0;j<b-1;j++){
            for(int i=0;i<a;i++){
                if(p==-1) edges.push_back(Edge(i+a*j,i+a*(j+1),(double)(rand())/RAND_MAX/2+0.5));
                else edges.push_back(Edge(i+a*j,i+a*(j+1),p));
                edgesNumber++;
                degree[i+a*j]++;
                degree[i+a*(j+1)]++;
            }
        }
        return 1;
    }

    void swapVertices(int a,int b){
        if(a!=b){
            for(int i=0;i<edgesNumber;i++){
                if(edges[i].in==a) edges[i].in=b;
                    else if(edges[i].in==b) edges[i].in=a;
                if(edges[i].out==a) edges[i].out=b;
                    else if(edges[i].out==b) edges[i].out=a;
            }
            swap(degree[a],degree[b]);
        }
    }

    void eliminateLeaves(int logs=-1){
        int s = 0, jsave = 0,flag = 1;

        while(flag==1){
            for(int i=dimension-1;i>=0;i--){
                if(degree[i]==1){
                    for(int j=0;j<edgesNumber;j++){
                        if(edges[j].in==i){
                            s=edges[j].out;
                            jsave=j;
                            break;
                        }
                        if(edges[j].out==i){
                            s=edges[j].in;
                            jsave=j;
                            break;
                        }
                    }
                    degree[s]--;
                    swapVertices(dimension-1,i);
                    if(logs!=-1) cout << "Vertex " << i << " was eliminated as a list" << endl;
                    //reliability[s]+=reliability[dimension-1];
                    //multipliers.push_back(eliminateEdge(s,dimension-1));
                    multipliers.push_back(edges[jsave].reliability);
                    edges.erase(edges.begin()+jsave);//можно функцией
                    degree.erase(degree.begin()+dimension-1);
                    //reliability.erase(reliability.begin()+dimension-1);
                    dimension--;
                    edgesNumber--;
                }
            }
            flag=0;
            for(int i=dimension-1;i>=0;i--) if(degree[i]==1) flag=1;
        }
    }

    void eliminateIsolated(int logs=-1){
        int flag = 1;
        while(flag==1){
            flag = 0;
            for(int i=0;i<dimension;i++){
                    if(degree[i]==0) {
                        swapVertices(i,dimension-1);
                        degree.erase(degree.end()-1);
                        //degree.erase(degree.begin()+dimension-1);
                        //reliability.erase(reliability.begin()+dimension-1);
                        dimension--;
                        if(logs!=-1) cout << "Vertex " << i << " was eliminated as an isolated one" << endl;
                    }
            }
            for(int i=0;i<dimension;i++) if(degree[i]==0) flag=1;
        }
    }

    /*void reliability_list(){
        cout << "List of vertices reliabilities: ";
        for(int i=0;i<dimension;i++) cout << reliability[i] << " ";
        cout << endl;
    }*/

    void printDegreesList(){
        cout << "List of vertices degrees: ";
        for(int i=0;i<dimension;i++) cout << degree[i] << " ";
        cout << endl;
    }

    void contractVertices(int a, int b){
        int flag = 0;
        eliminateEdge(a,b);
        for(int i=0;i<edgesNumber;i++){
            if(edges[i].in==b){
                for(int j=0;j<edgesNumber;j++){
                    if(((edges[i].out==edges[j].in)&&(edges[j].out==a))||((edges[i].out==edges[j].out)&&(edges[j].in==a))){
                        edges[j].reliability=edges[j].reliability+edges[i].reliability-edges[i].reliability*edges[j].reliability;
                        if(edges[j].in==a) degree[edges[j].out]--;
                        else degree[edges[j].in]--;
                        edges.erase(edges.begin()+i);
                        edgesNumber--;
                        flag=1;
                        break;
                    }
                }
                if(flag==0){
                    edges[i].in=a;
                    degree[a]++;
                }
                else i--;
                flag=0;
                degree[b]--;
            }
            else if(edges[i].out==b){
                for(int j=0;j<edgesNumber;j++){
                    if(((edges[i].in==edges[j].in)&&(edges[j].out==a))||((edges[i].in==edges[j].out)&&(edges[j].in==a))){
                        edges[j].reliability=edges[j].reliability+edges[i].reliability-edges[i].reliability*edges[j].reliability;
                        if(edges[j].in==a) degree[edges[j].out]--;
                        else degree[edges[j].in]--;
                        edges.erase(edges.begin()+i);
                        edgesNumber--;
                        flag=1;
                        break;
                    }
                }
                if(flag==0){
                    edges[i].out=a;
                    degree[a]++;
                }
                else i--;
                flag=0;
                degree[b]--;
            }
        }
    }

    double eliminateEdge(int a,int b){
        double buff = 0;
        for(int i=0;i<edgesNumber;i++){
            if(((edges[i].in==a)&&(edges[i].out==b))||((edges[i].in==b)&&(edges[i].out==a))){
                buff=edges[i].reliability;
                degree[a]--;
                degree[b]--;
                edges.erase(edges.begin()+i);
                edgesNumber--;
                break;
            }
        }
        return buff;
    }

    void eliminateChains(int logs = -1, double pr = 0){

        vector<int> chain1;
        vector<int> chain2;
        int s, psave = 0, isCycle = 0, flag = 1, exist = 0;
        double summ = 0, buff, mult = 1;

        while(flag==1){
            flag=0;
            for(int i=0;i<dimension;i++){
                if(degree[i]==2) {
                    chain2.push_back(i);
                    for(int j=0;j<edgesNumber;j++){
                        if((edges[j].in==i)||(edges[j].out==i)) {

                            if(edges[j].in==i) {
                                chain1.push_back(edges[j].out);
                                s=edges[j].out;
                            }
                            else {
                                chain1.push_back(edges[j].in);
                                s=edges[j].in;
                            }

                            while(degree[s]==2){
                                for(int p=0;p<edgesNumber;p++){
                                    if(((edges[p].in==s)&&(check(chain1,edges[p].out)==0)&&(check(chain2,edges[p].out)==0))||((edges[p].out==s)&&(check(chain1,edges[p].in)==0)&&(check(chain2,edges[p].in)==0))){
                                        if(edges[p].in==s) {
                                                chain1.push_back(edges[p].out);
                                                s=edges[p].out;
                                                psave=p;
                                        }
                                        else {
                                            chain1.push_back(edges[p].in);
                                            s=edges[p].in;
                                            psave=p;
                                        }
                                        break;
                                    }

                                }
                                for(int u=0;u<edgesNumber;u++){
                                    if((degree[s]==2)&&(psave!=u)&&(edges[u].in==s)&&((check(chain2,edges[u].out)==1)||(check(chain1,edges[u].out)==1))) {
                                        if(logs!=-1) cout << "Cycle ";
                                        isCycle=1;
                                        break;
                                    }
                                    if((degree[s]==2)&&(psave!=u)&&(edges[u].out==s)&&((check(chain2,edges[u].in)==1)||(check(chain1,edges[u].in)==1))) {
                                        if(logs!=-1) cout << "Cycle ";
                                        isCycle=1;
                                        break;
                                    }
                                }
                                if(isCycle==1) break;
                            }
                        break;
                        }
                    }
                    /////////////////////
                    for(int j=0;j<edgesNumber;j++){
                        if(((edges[j].in==i)&&(check(chain1,edges[j].out)==0)&&(check(chain2,edges[j].out)==0))||((edges[j].out==i)&&(check(chain1,edges[j].in)==0)&&(check(chain2,edges[j].in)==0))){

                            if(edges[j].in==i) {
                                chain2.push_back(edges[j].out);
                                s=edges[j].out;
                            }
                            else {
                                chain2.push_back(edges[j].in);
                                s=edges[j].in;
                            }
                            while(degree[s]==2){
                                for(int p=0;p<edgesNumber;p++){
                                    if(((edges[p].in==s)&&(check(chain1,edges[p].out)==0)&&(check(chain2,edges[p].out)==0))||((edges[p].out==s)&&(check(chain1,edges[p].in)==0)&&(check(chain2,edges[p].in)==0))){
                                        if(edges[p].in==s) {
                                                chain2.push_back(edges[p].out);
                                                s=edges[p].out;
                                                psave=p;
                                        }
                                        else {
                                            chain2.push_back(edges[p].in);
                                            s=edges[p].in;
                                            psave=p;
                                        }
                                        break;
                                    }
                                }
                                for(int u=0;u<edgesNumber;u++){
                                    if((degree[s]==2)&&(psave!=u)&&(edges[u].in==s)&&((check(chain2,edges[u].out)==1)||(check(chain1,edges[u].out)==1))) {
                                        if(logs!=-1) cout << "Cycle ";
                                        isCycle=1;
                                        break;
                                    }
                                    if((degree[s]==2)&&(psave!=u)&&(edges[u].out==s)&&((check(chain2,edges[u].in)==1)||(check(chain1,edges[u].in)==1))) {
                                        if(logs!=-1) cout << "Cycle ";
                                        isCycle=1;
                                        break;
                                    }
                                }
                                if(isCycle==1) break;
                            }
                        break;
                        }
                    }
                    if((isCycle==0)&&(logs!=-1)) cout << "Chain ";
                    for(unsigned int j=0;j<chain1.size();j++) chain2.insert(chain2.begin(),chain1[j]);
                    chain1.clear();
                    if(logs!=-1) for(unsigned int j=0;j<chain2.size();j++) cout << chain2[j] << " ";
                    if((isCycle==0)&&(logs!=-1)) cout << "was reduced to an edge with probability ";
                break;
                }

            }
        //подсчет веро€тности ребра, замен€ющего цепь
            if((isCycle==0)&&(chain2.size()>0)){
                summ = 0;
                mult = 1;
                for(unsigned int j=1;j<chain2.size();j++){
                    buff=eliminateEdge(chain2[j-1],chain2[j]);
                    summ+=1.0/buff;
                    mult*=buff;
                }
                summ = summ-chain2.size()+2;

                UB-=pr*(1-mult*summ);
                pr*=mult*summ;

                multipliers.push_back(mult*summ);
                summ=1.0/summ;
                exist=0;

                for(int p=0;p<edgesNumber;p++){
                    if(((edges[p].in==chain2[0])&&(edges[p].out==chain2.back()))||((edges[p].out==chain2[0])&&(edges[p].in==chain2.back()))) {
                        edges[p].reliability=edges[p].reliability+summ-edges[p].reliability*summ;
                        exist=1;
                        if(logs!=-1) cout << edges[p].reliability << endl;
                        break;
                    }
                }
                if(exist==0){
                    edges.push_back(Edge(chain2[0],chain2.back(),summ));
                    degree[chain2[0]]++;
                    degree[chain2.back()]++;
                    edgesNumber++;
                    if(logs!=-1) cout << summ << endl;
                }
            chain2.clear();
            }

            if((isCycle==1)&&(chain2.size()>0)){
                mult=1;
                summ=1;
                for(unsigned int j=1;j<chain2.size();j++){
                    buff=eliminateEdge(chain2[j-1],chain2[j]);
                    mult*=buff;
                    summ+=(1-buff)/buff;
                }
                buff=eliminateEdge(chain2[0],chain2.back());
                mult*=buff;
                summ+=(1-buff)/buff;
                UB-=pr*(1-mult*summ);
                pr*=mult*summ;
                multipliers.push_back(mult*summ);
                if(logs!=-1) cout << " was cutted out with reliability " << mult*summ << endl;
                if((int)chain2.size()==dimension) flag=-1;
            }
            chain2.clear();
            isCycle = 0;
            summ = 0;
            exist = 0;
            eliminateIsolated();
            eliminateLeaves();
            for(int i=0;i<dimension;i++)  if((degree[i]==2)&&(flag!=-1)) flag=1;
        }
    }

    double doTDF(){//Third Degree Factorization
        double c1=1,c2=1,c3=1,c4=1,c5=1,p1=-1,p2=-1,p3=-1;
        int v=-1,k=0,v1=-1,v2=-1,v3=-1,flag=0,minimalDegree=dimension;

        if(dimension==5) return product(multipliers)*fc5();
        if(dimension==4) return product(multipliers)*fc4();
        eliminateLeaves();
        if(dimension==5) return product(multipliers)*fc5();
        if(dimension==4) return product(multipliers)*fc4();
        if((dimension<2)||(edgesNumber==0)) return product(multipliers);
        eliminateChains();
        if((dimension<2)||(edgesNumber==0)) return product(multipliers);
        if(dimension==5) return product(multipliers)*fc5();
        if(dimension==4) return product(multipliers)*fc4();

        for(v=0;v<dimension;v++){
            if(degree[v]==3) break;
        }
        if((v!=-1)&&(v!=dimension)&&(edgesNumber>=3)&&(dimension>4)){
            for(int i=0;i<edgesNumber;i++){
                if(edges[i].in==v) {
                    v1=edges[i].out;
                    p1=edges[i].reliability;
                    k=i;
                    break;
                }
                if(edges[i].out==v) {
                    v1=edges[i].in;
                    p1=edges[i].reliability;
                    k=i;
                    break;
                }
            }
            for(int i=k+1;i<edgesNumber;i++){
                if(edges[i].in==v) {
                    v2=edges[i].out;
                    p2=edges[i].reliability;
                    k=i;
                    break;
                }
                if(edges[i].out==v) {
                    v2=edges[i].in;
                    p2=edges[i].reliability;
                    k=i;
                    break;
                }
            }
            for(int i=k+1;i<edgesNumber;i++){
                if(edges[i].in==v) {
                    v3=edges[i].out;
                    p3=edges[i].reliability;
                    break;
                }
                if(edges[i].out==v) {
                    v3=edges[i].in;
                    p3=edges[i].reliability;
                    break;
                }
            }

            RandomGraph G1(*this);
            G1.multipliers.clear();
            G1.multipliers.push_back(1);
            G1.contractVertices(v,v1);
            G1.contractVertices(v,v2);
            G1.contractVertices(v,v3);
            G1.eliminateIsolated();
            G1.eliminateLeaves();
            G1.eliminateChains();
            if(G1.dimension==5) G1.multipliers.push_back(G1.fc5());
            else if(G1.dimension==4) G1.multipliers.push_back(G1.fc4());
            else {
                c1=G1.doTDF();
                flag=1;
            }
            if(flag==0) c1*=product(G1.multipliers);
            flag=0;
            G1.deleteGraph();

            RandomGraph G2(*this);
            G2.multipliers.clear();
            G2.multipliers.push_back(1);
            G2.eliminateEdge(v,v3);
            G2.contractVertices(v,v1);
            G2.contractVertices(v,v2);
            G2.eliminateIsolated();
            if(G2.isConnected()==1){
                G2.eliminateLeaves();
                G2.eliminateChains();
                if(G2.dimension==5) G2.multipliers.push_back(G2.fc5());
                else if(G2.dimension==4) G2.multipliers.push_back(G2.fc4());
                else {
                    c2=G2.doTDF();
                    flag=1;
                }
                if(flag==0) c2*=product(G2.multipliers);
                flag=0;
            }
            else c2=0;
            G2.deleteGraph();

            RandomGraph G3(*this);
            G3.multipliers.clear();
            G3.multipliers.push_back(1);
            G3.eliminateEdge(v,v2);
            G3.contractVertices(v,v1);
            G3.contractVertices(v,v3);
            G3.eliminateIsolated();
            if(G3.isConnected()==1){
                G3.eliminateLeaves();
                G3.eliminateChains();
                if(G3.dimension==5) G3.multipliers.push_back(G3.fc5());
                else if(G3.dimension==4) G3.multipliers.push_back(G3.fc4());
                else {
                    c3=G3.doTDF();
                    flag=1;
                }
                if(flag==0) c3*=product(G3.multipliers);
                flag=0;
            }
            else c3=0;
            G3.deleteGraph();

            RandomGraph G4(*this);
            G4.multipliers.clear();
            G4.multipliers.push_back(1);
            G4.eliminateEdge(v,v1);
            G4.contractVertices(v,v2);
            G4.contractVertices(v,v3);
            G4.eliminateIsolated();
            if(G4.isConnected()==1){
                G4.eliminateLeaves();
                G4.eliminateChains();
                if(G4.dimension==5) G4.multipliers.push_back(G4.fc5());
                else if(G4.dimension==4) G4.multipliers.push_back(G4.fc4());
                else {
                    c4=G4.doTDF();
                    flag=1;
                }
                if(flag==0) c4*=product(G4.multipliers);
                flag=0;
            }
            else c4=0;
            G4.deleteGraph();

            RandomGraph G5(*this);
            G5.multipliers.clear();
            G5.multipliers.push_back(1);
            G5.eliminateEdge(v,v1);
            G5.eliminateEdge(v,v2);
            G5.eliminateEdge(v,v3);
            G5.eliminateIsolated();
            if(G5.isConnected()==1){
                G5.eliminateLeaves();
                G5.eliminateChains();
                if(G5.dimension==5) G5.multipliers.push_back(G5.fc5());
                else if(G5.dimension==4) G5.multipliers.push_back(G5.fc4());
                else {
                    c5=G5.doTDF();
                    flag=1;
                }
                if(flag==0) c5*=product(G5.multipliers);
                flag=0;
            }
            else c5=0;
            G5.deleteGraph();

            return product(multipliers)*(c1*p1*p2*p3+c2*p1*p2*(1-p3)+c3*p1*(1-p2)*p3+c4*(1-p1)*p2*p3+((1-p1)*(1-p2)*p3+p1*(1-p2)*(1-p3)+(1-p1)*p2*(1-p3))*c5);
        }
        else {
            minimalDegree=dimension;
            for(int i=0;i<dimension;i++){
                if(degree[i]<minimalDegree) minimalDegree=degree[i];
            }
            if((minimalDegree>3)&&(minimalDegree!=dimension)) return product(multipliers)*doAF();
            else return product(multipliers);
        }
    }

    void dfs(int a, int *used){
        used[a]=1;
        for(int i=0;i<edgesNumber;i++){
            if((edges[i].in==a)&&(used[edges[i].out]==0)) dfs(edges[i].out,used);
            if((edges[i].out==a)&&(used[edges[i].in]==0)) dfs(edges[i].in,used);
        }
    }

    double doMS(){//Moore-Shannon
        double p;
        int *used;
        if(dimension==5) return fc5();
        if(dimension==4) return fc4();
        if(dimension==3) return fc3();
        if(edges.size()>2){
            RandomGraph G1(*this),G2(*this);
            G1.contractVertices(G1.edges[0].in,G1.edges[0].out);
            G1.eliminateIsolated();
            p=G2.eliminateEdge(G2.edges[0].in,G2.edges[0].out);
            for(int i=0;i<G2.dimension;i++) if(G2.degree[i]==0) return p*G1.doMS();
            used=new int[G2.dimension];
            for(int i=0;i<G2.dimension;i++) used[i]=0;
            G2.dfs(G2.edges[0].in,used);
            for(int i=0;i<G2.dimension;i++) if(used[i]==0) {
                delete[] used;
                return p*G1.doMS();
            }
            delete[] used;
            return p*G1.doMS()+(1-p)*G2.doMS();
        }
        if(edges.size()==2) return edges[0].reliability*edges[1].reliability;
        if(edges.size()==1) return edges[0].reliability;
        else return 1;
    }

    double doAF(){//Adaptive Factorization
        double p;
        int *used,minimalDegree=dimension;

        if(dimension==5) return fc5();
        if(dimension==4) return fc4();
        if(dimension==3) return fc3();
        multipliers.clear();
        multipliers.push_back(1);
        for(int i=0;i<dimension;i++){
            if(degree[i]<minimalDegree) minimalDegree=degree[i];
        }
        if(minimalDegree>3){
            RandomGraph G1(*this),G2(*this);
            G1.contractVertices(G1.edges[0].in,G1.edges[0].out);
            G1.eliminateIsolated();
            p=G2.eliminateEdge(G2.edges[0].in,G2.edges[0].out);
            minimalDegree=G1.dimension;
            for(int i=0;i<G1.dimension;i++){
                if(G1.degree[i]<minimalDegree) minimalDegree=G1.degree[i];
            }
            for(int i=0;i<G2.dimension;i++) if(G2.degree[i]==0) {
                if(minimalDegree<=3) return p*G1.doTDF();
                else return p*G1.doAF();
            }
            used=new int[G2.dimension];
            for(int i=0;i<G2.dimension;i++) used[i]=0;
            G2.dfs(G2.edges[0].in,used);
            for(int i=0;i<G2.dimension;i++) if(used[i]==0) {
                delete[] used;
                if(minimalDegree<=3) return p*G1.doTDF();
                else return p*G1.doAF();
            }
            delete[] used;
            if(minimalDegree<=3){
                minimalDegree=G2.dimension;
                for(int i=0;i<G2.dimension;i++){
                    if(G2.degree[i]<minimalDegree) minimalDegree=G2.degree[i];
                }
                if(minimalDegree<=3) return p*G1.doTDF()+(1-p)*G2.doTDF();
                else return p*G1.doTDF()+(1-p)*G2.doAF();
            }
            return p*G1.doAF()+(1-p)*G2.doAF();
        }
        else return doTDF();
    }

    double fc5(){//Fast Calculation Formula for 5 vertices
        double p,*q,*k;
        q=new double[20];
        k=new double[15];
        for(int i=0;i<10;i++) q[i]=0;
        for(int i=0;i<edgesNumber;i++){
            if(((edges[i].in==0)&&(edges[i].out==1))||((edges[i].in==1)&&(edges[i].out==0))) q[0]=edges[i].reliability;
            else if(((edges[i].in==0)&&(edges[i].out==2))||((edges[i].in==2)&&(edges[i].out==0))) q[1]=edges[i].reliability;
            else if(((edges[i].in==0)&&(edges[i].out==3))||((edges[i].in==3)&&(edges[i].out==0))) q[2]=edges[i].reliability;
            else if(((edges[i].in==0)&&(edges[i].out==4))||((edges[i].in==4)&&(edges[i].out==0))) q[3]=edges[i].reliability;
            else if(((edges[i].in==1)&&(edges[i].out==2))||((edges[i].in==2)&&(edges[i].out==1))) q[4]=edges[i].reliability;
            else if(((edges[i].in==1)&&(edges[i].out==3))||((edges[i].in==3)&&(edges[i].out==1))) q[5]=edges[i].reliability;
            else if(((edges[i].in==1)&&(edges[i].out==4))||((edges[i].in==4)&&(edges[i].out==1))) q[6]=edges[i].reliability;
            else if(((edges[i].in==2)&&(edges[i].out==3))||((edges[i].in==3)&&(edges[i].out==2))) q[7]=edges[i].reliability;
            else if(((edges[i].in==2)&&(edges[i].out==4))||((edges[i].in==4)&&(edges[i].out==2))) q[8]=edges[i].reliability;
            else if(((edges[i].in==3)&&(edges[i].out==4))||((edges[i].in==4)&&(edges[i].out==3))) q[9]=edges[i].reliability;
        }

        for(int i=10;i<20;i++) q[i]=1-q[i-10];

        k[0]=1-q[14]*(q[15]*q[16]+q[17]*q[18]);
        k[1]=1-q[17]*(q[11]*q[18]+q[12]*q[19]);
        k[2]=1-q[19]*(q[12]*q[15]+q[13]*q[16]);
        k[3]=1-q[13]*(q[10]*q[11]+q[16]*q[18]);
        k[4]=1-q[10]*(q[11]*q[12]+q[14]*q[15]);
        k[5]=q[0]*q[7]*q[8]+q[0]*q[9]*(q[7]*q[18]+q[17]*q[8])+q[10]*q[17]*q[18]*(1-4*q[19]);
        k[6]=q[2]*q[3]*q[4]+q[4]*q[9]*(q[2]*q[13]+q[12]*q[3])+q[13]*q[14]*q[19];
        k[7]=q[0]*q[3]*q[7]+q[6]*q[7]*(q[0]*q[13]+q[10]*q[3])+q[10]*q[16]*q[17];
        k[8]=q[0]*q[1]*q[9]+q[4]*q[9]*(q[0]*q[11]+q[10]*q[1])+q[10]*q[14]*q[19];
        k[9]=q[3]*q[4]*q[5]+q[3]*q[7]*(q[4]*q[15]+q[14]*q[5])+q[13]*q[15]*q[17];
        k[10]=q[1]*q[5]*q[6]+q[1]*q[9]*(q[5]*q[16]+q[15]*q[6])+q[11]*q[16]*q[19];
        k[11]=q[2]*q[4]*q[6]+q[2]*q[8]*(q[4]*q[16]+q[14]*q[6])+q[12]*q[16]*q[18];
        k[12]=q[1]*q[3]*q[5]+q[5]*q[8]*(q[1]*q[13]+q[11]*q[3])+q[11]*q[13]*q[15];
        k[13]=q[1]*q[2]*q[6]+q[6]*q[7]*(q[1]*q[12]+q[11]*q[2])+q[11]*q[12]*q[16];
        k[14]=q[0]*q[2]*q[8]+q[5]*q[8]*(q[0]*q[12]+q[10]*q[2])+q[12]*q[15]*q[18];

        p=1-q[11]*q[12]*(q[10]*q[13]*k[0]+q[14]*q[15]*(q[13]*q[16]*k[5]+q[18]*q[19]*k[7]));
        p-=q[15]*q[16]*(q[10]*q[14]*k[1]+q[17]*q[18]*(q[10]*q[11]*k[6]+q[12]*q[13]*k[8]));
        p-=q[11]*q[17]*(q[14]*q[18]*k[2]+q[13]*q[19]*(q[10]*q[15]*k[11]+q[14]*q[16]*k[14]));
        p-=q[12]*q[19]*(q[15]*q[17]*k[3]+q[10]*q[16]*(q[11]*q[18]*k[9]+q[14]*q[17]*k[12]));
        p-=q[13]*q[18]*(q[16]*q[19]*k[4]+q[10]*q[14]*(q[12]*q[17]*k[10]+q[15]*q[19]*k[13]));
        delete[] q;
        delete[] k;
        edges.clear();
        degree.clear();
        dimension=0;
        edgesNumber=0;
        return p;
    }

    double fc4(){//Fast Calculation Formula for 4 vertices
        double p,*q;
        q=new double[6];
        for(int i=0;i<6;i++) q[i]=1;
        for(int i=0;i<edgesNumber;i++){
            if(((edges[i].in==0)&&(edges[i].out==1))||((edges[i].in==1)&&(edges[i].out==0))) q[0]=1-edges[i].reliability;
            else if(((edges[i].in==1)&&(edges[i].out==2))||((edges[i].in==2)&&(edges[i].out==1))) q[1]=1-edges[i].reliability;
            else if(((edges[i].in==2)&&(edges[i].out==3))||((edges[i].in==3)&&(edges[i].out==2))) q[2]=1-edges[i].reliability;
            else if(((edges[i].in==0)&&(edges[i].out==3))||((edges[i].in==3)&&(edges[i].out==0))) q[3]=1-edges[i].reliability;
            else if(((edges[i].in==1)&&(edges[i].out==3))||((edges[i].in==3)&&(edges[i].out==1))) q[4]=1-edges[i].reliability;
            else if(((edges[i].in==0)&&(edges[i].out==2))||((edges[i].in==2)&&(edges[i].out==0))) q[5]=1-edges[i].reliability;
        }
        p=1+2*(q[1]*q[3]*q[4]*q[5]*(q[0]+q[2]-0.5)+q[0]*q[2]*q[4]*q[5]*(q[1]+q[3]-0.5)+q[0]*q[1]*q[2]*q[3]*(q[4]+q[5]-0.5))-6*q[0]*q[1]*q[2]*q[3]*q[4]*q[5]-q[0]*q[1]*q[4]-q[0]*q[3]*q[5]-q[1]*q[2]*q[5]-q[2]*q[3]*q[4];
        delete[] q;
        edges.clear();
        degree.clear();
        dimension=0;
        edgesNumber=0;
        return p;
    }

    double fc3(){//Fast Calculation Formula for 3 vertices
        double a=edges[0].reliability,b=edges[1].reliability,c;
        if(edges.size()==3) c=edges[2].reliability;
        else c=0;
        edges.clear();
        degree.clear();
        dimension=0;
        edgesNumber=0;
        return a*b+b*c+a*c-2*a*b*c;
    }

    int isConnected(){
        int *used;
        used=new int[dimension];
        for(int i=0;i<dimension;i++) used[i]=0;
        dfs(edges[0].in,used);
        for(int i=0;i<dimension;i++) if(used[i]==0){
            delete[] used;
            return 0;
        }
        delete[] used;
        return 1;
    }

    /////////////////////////////////////////////////////

    double doTDFCB(int level=1){//Third Degree Factorization with Cumulative Bounds Method
        double p1=-1,p2=-1,p3=-1,buff;
        int v=-1,v1=-1,v2=-1,v3=-1,minimalDegree,k=-1;
        counter++;
        if(counter<=10000) {upper=UB; lower=LB; }
        if(counter==10000) cout << "under 10000 iterations: " << endl << "UB " << UB << endl << "LB " << LB << " time " << (long double)(clock()-st_time)/CLK_TCK << endl;
        if(level==0){
            multipliers.clear();
            eliminateLeaves();
            eliminateChains(-1,product(multipliers));
            currentReliability=product(multipliers);
            UB=currentReliability;
            if(dimension<2) return currentReliability;
            if(dimension==5) {
                buff=fc5();
                UB-=currentReliability*(1-buff);
                LB+=currentReliability*buff;
                return currentReliability*buff;
            }
            if(dimension==4) {
                buff=fc4();
                UB-=currentReliability*(1-buff);
                LB+=currentReliability*buff;
                return currentReliability*buff;
            }
        }
        else{
            multipliers.clear();
            eliminateLeaves();
            UB-=currentReliability*(1-product(multipliers));
            currentReliability*=product(multipliers);
            multipliers.clear();
            eliminateChains(-1,currentReliability);
            currentReliability*=product(multipliers);
            if(dimension==5) {
                buff=fc5();
                UB-=currentReliability*(1-buff);
                LB+=currentReliability*buff;
                return currentReliability*buff;
            }
            if(dimension==4) {
                buff=fc4();
                UB-=currentReliability*(1-buff);
                LB+=currentReliability*buff;
                return currentReliability*buff;
            }
        }

        for(v=0;v<dimension;v++){
            if(degree[v]==3) break;
        }
        if((v!=dimension)&&(edgesNumber>=3)&&(dimension>=4)){
            for(int i=0;i<edgesNumber;i++){
                if(edges[i].in==v) {
                    v1=edges[i].out;
                    p1=edges[i].reliability;
                    k=i;
                    break;
                }
                if(edges[i].out==v) {
                    v1=edges[i].in;
                    p1=edges[i].reliability;
                    k=i;
                    break;
                }
            }
            for(int i=k+1;i<edgesNumber;i++){
                if(edges[i].in==v) {
                    v2=edges[i].out;
                    p2=edges[i].reliability;
                    k=i;
                    break;
                }
                if(edges[i].out==v) {
                    v2=edges[i].in;
                    p2=edges[i].reliability;
                    k=i;
                    break;
                }
            }
            for(int i=k+1;i<edgesNumber;i++){
                if(edges[i].in==v) {
                    v3=edges[i].out;
                    p3=edges[i].reliability;
                    break;
                }
                if(edges[i].out==v) {
                    v3=edges[i].in;
                    p3=edges[i].reliability;
                    break;
                }
            }

            RandomGraph G1(*this);
            G1.multipliers.clear();
            G1.contractVertices(v,v1);
            G1.contractVertices(v,v2);
            G1.contractVertices(v,v3);
            G1.eliminateIsolated();
            G1.currentReliability*=p1*p2*p3;
            G1.eliminateLeaves();
            UB-=G1.currentReliability*(1-product(G1.multipliers));
            G1.currentReliability*=product(G1.multipliers);
            G1.multipliers.clear();
            G1.eliminateChains(-1,G1.currentReliability);
            G1.currentReliability*=product(G1.multipliers);
            if(G1.dimension==5) {
                buff=G1.fc5();
                UB-=G1.currentReliability*(1-buff);
                LB+=G1.currentReliability*buff;
                G1.currentReliability*=buff;
            }
            else if(G1.dimension==4) {
                buff=G1.fc4();
                UB-=G1.currentReliability*(1-buff);
                LB+=G1.currentReliability*buff;
                G1.currentReliability*=buff;
            }
            else if(G1.dimension==3) {
                buff=G1.fc3();
                UB-=G1.currentReliability*(1-buff);
                LB+=G1.currentReliability*buff;
                G1.currentReliability*=buff;
            }
            else if(G1.dimension==2) {
                buff=G1.edges[0].reliability;
                UB-=G1.currentReliability*(1-buff);
                LB+=G1.currentReliability*buff;
                G1.currentReliability*=buff;
            }
            else if(G1.dimension>5)G1.currentReliability*=G1.doTDFCB();
            else LB+=G1.currentReliability;

            RandomGraph G2(*this);
            G2.multipliers.clear();
            G2.eliminateEdge(v,v3);
            G2.contractVertices(v,v1);
            G2.contractVertices(v,v2);
            G2.eliminateIsolated();
            G2.currentReliability*=p1*p2*(1-p3);
            if(G2.isConnected()==1){
                G2.eliminateLeaves();
                UB-=G2.currentReliability*(1-product(G2.multipliers));
                G2.currentReliability*=product(G2.multipliers);
                G2.multipliers.clear();
                G2.eliminateChains(-1,G2.currentReliability);
                G2.currentReliability*=product(G2.multipliers);
                if(G2.dimension==5) {
                    buff=G2.fc5();
                    UB-=G2.currentReliability*(1-buff);
                    LB+=G2.currentReliability*buff;
                    G2.currentReliability*=buff;
                }
                else if(G2.dimension==4) {
                    buff=G2.fc4();
                    UB-=G2.currentReliability*(1-buff);
                    LB+=G2.currentReliability*buff;
                    G2.currentReliability*=buff;
                }
                else if(G2.dimension==3) {
                    buff=G2.fc3();
                    UB-=G2.currentReliability*(1-buff);
                    LB+=G2.currentReliability*buff;
                    G2.currentReliability*=buff;
                }
                else if(G2.dimension==2) {
                    buff=G2.edges[0].reliability;
                    UB-=G2.currentReliability*(1-buff);
                    LB+=G2.currentReliability*buff;
                    G2.currentReliability*=buff;
                }
                else if(G2.dimension>5) G2.currentReliability*=G2.doTDFCB();
                else LB+=G2.currentReliability;
            }
            else {UB-=G2.currentReliability; G2.currentReliability=0;}

            RandomGraph G3(*this);
            G3.multipliers.clear();
            G3.eliminateEdge(v,v2);
            G3.contractVertices(v,v1);
            G3.contractVertices(v,v3);
            G3.eliminateIsolated();
            G3.currentReliability*=p1*(1-p2)*p3;
            if(G3.isConnected()==1){
                G3.eliminateLeaves();
                UB-=G3.currentReliability*(1-product(G3.multipliers));
                G3.currentReliability*=product(G3.multipliers);
                G3.multipliers.clear();
                G3.eliminateChains(-1,G3.currentReliability);
                G3.currentReliability*=product(G3.multipliers);
                if(G3.dimension==5) {
                    buff=G3.fc5();
                    UB-=G3.currentReliability*(1-buff);
                    LB+=G3.currentReliability*buff;
                    G3.currentReliability*=buff;
                }
                else if(G3.dimension==4) {
                    buff=G3.fc4();
                    UB-=G3.currentReliability*(1-buff);
                    LB+=G3.currentReliability*buff;
                    G3.currentReliability*=buff;
                }
                else if(G3.dimension==3) {
                    buff=G3.fc3();
                    UB-=G3.currentReliability*(1-buff);
                    LB+=G3.currentReliability*buff;
                    G3.currentReliability*=buff;
                }
                else if(G3.dimension==2) {
                    buff=G3.edges[0].reliability;
                    UB-=G3.currentReliability*(1-buff);
                    LB+=G3.currentReliability*buff;
                    G3.currentReliability*=buff;
                }
                else if(G3.dimension>5) G3.currentReliability*=G3.doTDFCB();
                else LB+=G3.currentReliability;
            }
            else {UB-=G3.currentReliability; G3.currentReliability=0;}

            RandomGraph G4(*this);
            G4.multipliers.clear();
            G4.eliminateEdge(v,v1);
            G4.contractVertices(v,v2);
            G4.contractVertices(v,v3);
            G4.eliminateIsolated();
            G4.currentReliability*=(1-p1)*p2*p3;
            if(G4.isConnected()==1){
                G4.eliminateLeaves();
                UB-=G4.currentReliability*(1-product(G4.multipliers));
                G4.currentReliability*=product(G4.multipliers);
                G4.multipliers.clear();
                G4.eliminateChains(-1,G4.currentReliability);
                G4.currentReliability*=product(G4.multipliers);
                if(G4.dimension==5) {
                    buff=G4.fc5();
                    UB-=G4.currentReliability*(1-buff);
                    LB+=G4.currentReliability*buff;
                    G4.currentReliability*=buff;
                }
                else if(G4.dimension==4) {
                    buff=G4.fc4();
                    UB-=G4.currentReliability*(1-buff);
                    LB+=G4.currentReliability*buff;
                    G4.currentReliability*=buff;
                }
                else if(G4.dimension==3) {
                    buff=G4.fc3();
                    UB-=G4.currentReliability*(1-buff);
                    LB+=G4.currentReliability*buff;
                    G4.currentReliability*=buff;
                }
                else if(G4.dimension==2) {
                    buff=G4.edges[0].reliability;
                    UB-=G4.currentReliability*(1-buff);
                    LB+=G4.currentReliability*buff;
                    G4.currentReliability*=buff;
                }
                else if(G4.dimension>5) G4.currentReliability*=G4.doTDFCB();
                else LB+=G4.currentReliability;
            }
            else {UB-=G4.currentReliability; G4.currentReliability=0;}

            RandomGraph G5(*this);
            G5.multipliers.clear();
            G5.eliminateEdge(v,v1);
            G5.eliminateEdge(v,v2);
            G5.eliminateEdge(v,v3);
            G5.eliminateIsolated();
            UB-=G5.currentReliability*(1-p1)*(1-p2)*(1-p3);
            G5.currentReliability*=(1-p1)*(1-p2)*p3+(1-p1)*p2*(1-p3)+p1*(1-p2)*(1-p3);
            if(G5.isConnected()==1){
                G5.eliminateLeaves();
                UB-=G5.currentReliability*(1-product(G5.multipliers));
                G5.currentReliability*=product(G5.multipliers);
                G5.multipliers.clear();
                G5.eliminateChains(-1,G5.currentReliability);
                G5.currentReliability*=product(G5.multipliers);
                if(G5.dimension==5) {
                    buff=G5.fc5();
                    UB-=G5.currentReliability*(1-buff);
                    LB+=G5.currentReliability*buff;
                    G5.currentReliability*=buff;
                }
                else if(G5.dimension==4) {
                    buff=G5.fc4();
                    UB-=G5.currentReliability*(1-buff);
                    LB+=G5.currentReliability*buff;
                    G5.currentReliability*=buff;
                }
                else if(G5.dimension==3) {
                    buff=G5.fc3();
                    UB-=G5.currentReliability*(1-buff);
                    LB+=G5.currentReliability*buff;
                    G5.currentReliability*=buff;
                }
                else if(G5.dimension==2) {
                    buff=G5.edges[0].reliability;
                    UB-=G5.currentReliability*(1-buff);
                    LB+=G5.currentReliability*buff;
                    G5.currentReliability*=buff;
                }
                else if(G5.dimension>5) G5.currentReliability*=G5.doTDFCB();
                else LB+=G5.currentReliability;
            }
            else {UB-=G5.currentReliability; G5.currentReliability=0;}

            return G1.currentReliability+G2.currentReliability+G3.currentReliability+G4.currentReliability+G5.currentReliability;
        }
        else{
            minimalDegree=dimension;
            for(int i=0;i<dimension;i++){
                if(degree[i]<minimalDegree) minimalDegree=degree[i];
            }
            if((minimalDegree>3)&&(minimalDegree!=dimension)) {currentReliability*=doTDFCBAF(); /*LB+=currentReliability;*/ return currentReliability;}
            else {LB+=currentReliability; return 1;}
        }
    }

    double doTDFCBAF(){//Adaptive Factorization for Third Degree Factorization with Cumulative Bounds Method
        double p,buff;
        int *used,minimalDegree=dimension;

        counter++;
        if(counter<=10000) {upper=UB; lower=LB;}
        if(counter==10000) cout << "under 10000 iterations: " << endl << "UB " << UB << endl << "LB " << LB << " time " << (long double)(clock()-st_time)/CLK_TCK << endl;
        if(dimension==5) {
            buff=fc5();
            UB-=currentReliability*(1-buff);
            LB+=currentReliability*buff;
            return fc5();
        }
        if(dimension==4) {
            buff=fc4();
            UB-=currentReliability*(1-buff);
            LB+=currentReliability*buff;
            return fc4();
        }
        if(dimension==3) {
            buff=fc3();
            UB-=currentReliability*(1-buff);
            LB+=currentReliability*buff;
            return fc3();
        }
        for(int i=0;i<dimension;i++){
            if(degree[i]<minimalDegree) minimalDegree=degree[i];
        }
        if(minimalDegree>3){
            RandomGraph G1(*this),G2(*this);
            G1.contractVertices(G1.edges[0].in,G1.edges[0].out);
            G1.eliminateIsolated();
            p=G2.eliminateEdge(G2.edges[0].in,G2.edges[0].out);
            G1.currentReliability*=p;
            G2.currentReliability*=1-p;
            minimalDegree=G1.dimension;
            for(int i=0;i<G1.dimension;i++){
                if(G1.degree[i]<minimalDegree) minimalDegree=G1.degree[i];
            }
            for(int i=0;i<G2.dimension;i++) if(G2.degree[i]==0) {
                if(minimalDegree<=3) return G1.doTDFCB();
                else return G1.doTDFCBAF();
            }

            used=new int[G2.dimension];
            for(int i=0;i<G2.dimension;i++) used[i]=0;
            G2.dfs(G2.edges[0].in,used);
            for(int i=0;i<G2.dimension;i++) if(used[i]==0) {
                delete[] used;
                if(minimalDegree<=3) return G1.doTDFCB();
                else return G1.doTDFCBAF();
            }
            delete[] used;
            if(minimalDegree<=3){
                minimalDegree=G2.dimension;
                for(int i=0;i<G2.dimension;i++){
                    if(G2.degree[i]<minimalDegree) minimalDegree=G2.degree[i];
                }
                if(minimalDegree<=3) return G1.doTDFCB()+G2.doTDFCB();
                else return G1.doTDFCB()+G2.doTDFCBAF();
            }
            minimalDegree=G2.dimension;
            /*for(int i=0;i<G2.dimension;i++){
                if(G2.degree[i]<minimalDegree) minimalDegree=G2.degree[i];
            }*/
            if(minimalDegree<=3) return G1.doTDFCBAF()+G2.doTDFCB();
            return G1.doTDFCBAF()+G2.doTDFCBAF();
        }
        else return doTDFCB();
    }

    double doRMS(int level=1){//Moore-Shannon with Reduction Methods
        double p,prod;
        int *used;
        multipliers.clear();
        if(level==0){
            eliminateLeaves();
            eliminateChains();
        }
        prod=product(multipliers);
        if(dimension==5) return prod*fc5();
        if(dimension==4) return prod*fc4();
        if(dimension==3) return prod*fc3();
        if(edges.size()>2){
            RandomGraph G1(*this),G2(*this);
            G1.contractVertices(G1.edges[0].in,G1.edges[0].out);
            G1.eliminateIsolated();
            G1.eliminateLeaves();
            G1.eliminateChains();
            p=G2.eliminateEdge(G2.edges[0].in,G2.edges[0].out);
            G2.eliminateLeaves();
            G2.eliminateChains();
            for(int i=0;i<G2.dimension;i++) if(G2.degree[i]==0) return prod*product(G1.multipliers)*p*G1.doRMS();
            used=new int[G2.dimension];
            for(int i=0;i<G2.dimension;i++) used[i]=0;
            G2.dfs(G2.edges[0].in,used);
            for(int i=0;i<G2.dimension;i++) if(used[i]==0) {
                delete[] used;
                return prod*product(G1.multipliers)*p*G1.doRMS();
            }
            delete[] used;
            return prod*(product(G1.multipliers)*p*G1.doRMS()+product(G2.multipliers)*(1-p)*G2.doRMS());
        }
        if(edges.size()==2) return prod*edges[0].reliability*edges[1].reliability;
        if(edges.size()==1) return prod*edges[0].reliability;
        else return prod;
    }

    int createWheelGraph(int n, double p){
        dimension=n;
        edgesNumber=0;
        for(int i=0;i<n;i++) degree.push_back(0);
        for(int i=1;i<n-1;i++){
            edges.push_back(Edge(i,i+1,p));
            degree[i]++;
            degree[i+1]++;
            edgesNumber++;
        }
        edges.push_back(Edge(1,n-1,p));
        degree[1]++;
        degree[n-1]++;
        edgesNumber++;
        for(int i=1;i<n;i++){
            edges.push_back(Edge(0,i,p));
            degree[0]++;
            degree[i]++;
            edgesNumber++;
        }
        return 1;
    }

};

int main(){
    int n,e;
    unsigned long long int start_time, end_time;

    cin>>n>>e;

    RandomGraph G;

    G.createGridGraph(n,e,0.9);
    RandomGraph G1(G),G2(G);

    cout << "=====================================" << endl;
    start_time =  clock();
    cout << "RMS " << G.doRMS(0);
    end_time =  clock();
    cout << " time " << (long double)(end_time-start_time)/CLK_TCK << endl;

    start_time =  clock();
    cout << "TDF " << G1.doTDF();
    end_time =  clock();
    cout << " time " << (long double)(end_time-start_time)/CLK_TCK << endl;

    cout << "=====================================" << endl;
    UB=1;
    LB=0;
    counter=0;
    cout << "TDF with cumulative bounds:" << endl;
    start_time = clock();
    st_time = clock();
    G2.doTDFCB(0);
    end_time = clock();
    if(counter>10000) cout << endl << "After a whole process:" << endl;
    cout << "UB " << UB << endl << "LB " << LB;
    cout << " time " << (long double)(end_time-start_time)/CLK_TCK << endl;
    cout << "Taken iterations at all: " << counter << endl;
    cout << "=====================================" << endl;

    return 0;
}
