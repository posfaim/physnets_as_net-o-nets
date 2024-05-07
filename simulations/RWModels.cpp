#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include <iostream>
#include <chrono>
#include <random>
#include <array>
#include <algorithm>

#include <iomanip>

#define SMALL_NUMBER 1e-14
#define MAX_NODE_TRIAL 100000
#define M_PI           3.14159265358979323846  /* pi */


// construct a trivial random generator engine from a time-based seed:
unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
std::default_random_engine generator(seed);


int RandomStartPoint(int N_unoccupied, int *grid){
    int sp, idx; 
    
    //random start point
     std::uniform_int_distribution<int> rand_start(0, N_unoccupied-1);
     sp = rand_start(generator);
     idx=-1;
     do{
         idx++;
         if(grid[idx]==-1) sp--;
     }while(sp>=0);

    return idx;
}


int LERW(int n, int N, int fill, int l0, std::vector<std::vector<int>> &grid_al, int *grid, std::vector<int> &walk){
    //loop-erased random walk
    int next;
    int l;
    int rd;
    
    //start point
    walk[0]=RandomStartPoint(N-fill, grid);
    l=1;
    
    while(true){
        //get next location
        std::uniform_int_distribution<int> rand_dir(0, grid_al[walk[l-1]].size()-1); // uniform, unbiased
        rd = rand_dir(generator);
        next = grid_al[walk[l-1]][rd];
        
        //is this the end?
        if(grid[next]>-1){
            walk[l]=next;
            l++;
            break;
        }
        
        //is this a loop?
        for(int s=0;s<l;s++){
            if(walk[s]==next){
                l=s;
                break;
            }
        }
        walk[l]=next;
        l++;
        
        //initial condition
        if(n==0 && l>l0) break;   
    }
        
    return l;
}

int SRW(int n, int N, int fill, int l0, std::vector<std::vector<int>> &grid_al, int *grid, std::vector<int> &walk){
    // simple random walk
    int next;
    int l;
    int rd;
    
    //start point
    walk[0]=RandomStartPoint(N-fill, grid);
    l=1;
    
    while(true){
        //get next location
        std::uniform_int_distribution<int> rand_dir(0, grid_al[walk[l-1]].size()-1); // uniform, unbiased
        rd = rand_dir(generator);
        next = grid_al[walk[l-1]][rd];
        
        if(l>(int) walk.size()-1)
            walk.resize(walk.size()+10);
        walk[l]=next;
        l++;
        
        //is this the end?
        if(grid[next]>-1){
            break;
        }
        
        //initial condition
        if(n==0 && l>l0) break;
    }
    
    return l;
}

int RPC(int n, int N, int fill, int l0, std::vector<std::vector<int>> &grid_al, int *grid, std::vector<int> &walk){
    // simple random walk
    int next;
    int l;
    
    std::uniform_int_distribution<int> rand_dir(0, N-1); // uniform, unbiased
    
    //start point
    walk[0]=RandomStartPoint(N-fill, grid);
    l=1;
    
    while(true){
        //get next location
        
        next = rand_dir(generator);
        
        if(l>(int) walk.size()-1)
            walk.resize(walk.size()+10);
        walk[l]=next;
        l++;
        
        //is this the end?
        if(grid[next]>-1){
            break;
        }
        
        //initial condition
        if(n==0 && l>l0) break;
    }
    
    return l;
}


int MSAW(int n, int N, int fill, int l0, std::vector<std::vector<int>> &grid_al, int *grid, std::vector<int> &walk){
    // self-avoiding random walk
    int next=0, rd;
    int l;
    
    // preparation for collecting self-avoiding options
    int max_degree=0, num_options;
    std::vector<int> options;
    for(auto neighbors:grid_al)if(int(neighbors.size())>max_degree) max_degree=neighbors.size();
    options.resize(max_degree);
    
    //start point
    walk[0]=RandomStartPoint(N-fill, grid);
    l=1;
    
    while(true){
        //look for possible self-avoiding next steps
        num_options=0;
        for(auto site:grid_al[walk[l-1]])
            if (std::find(walk.begin(), walk.begin()+l, site) == walk.begin()+l) {
                options[num_options]=site;
                num_options++;
            }
        if(num_options==0)
            return -1;
        
        //get next location
        std::uniform_int_distribution<int> rand_dir(0, num_options-1); // uniform, unbiased
        rd = rand_dir(generator);
        next = options[rd];
        walk[l]=next;
        l++;
        
        //is this the end?
        if(grid[next]>-1){
            break;
        }
        
        //initial condition
        if(n==0 && l>l0) break;   
    }
        
    return l;
}


int RRAY(int n, int N, int fill, int L, int D, int *grid, std::vector<int> &walk){
    // start in a random direction, walk straight until you hit something
    // this can only work on a grid where a straight line is meaningful
    
    int next=0, rax,rd,c;
    int l,Lk;
    
    
    //get random direction
    std::uniform_int_distribution<int> rand_axis(0, D-1);
    std::uniform_real_distribution<> uniform(0.0, 1.0);
    rax = rand_axis(generator);
    if(uniform(generator)<.5)
        rd = 1;
    else
        rd = -1;
    Lk=std::pow(L,rax);
    
    //start point
    walk[0]=RandomStartPoint(N-fill, grid);
    l=1;
    
    while(true){
        c = (walk[l-1]/Lk)%L;
        next = walk[l-1] - c*Lk + ((c+rd+L)%L)*Lk;
        
        
        walk[l]=next;
        l++;
        
        //is this the end?
        if(grid[next]>-1 || l>L){
            break;
        }  
    }
        
    return l;
}


int Full_C(std::string type, std::vector<std::vector<int>> &grid_al, int l0, std::vector<std::array<int,2>> &vws, std::vector<std::vector<int>> &walks){
    int N = grid_al.size();
    
    //initialize the grid
    int *grid;
    grid = new int[N];
    for(int i=0;i<N;i++) grid[i] = -1;
    
    int fill, unique_l;
    std::vector<int> walk, unique_sites;
    std::vector<int>::iterator it;
    walk.resize(N+1);

    int l, n;
    n = 0;
    fill = 0;
   
    while(N>fill){
        //if (PyErr_CheckSignals() != 0)
        //    break;
        if(type.compare("LERW") == 0) {
            l = LERW(n, N, fill, l0, grid_al, grid, walk);
            unique_l = l-1;
        } else if(type.compare("SRW") == 0){
            l = SRW(n, N, fill, l0, grid_al, grid, walk);

            //get number of unique visited sites
            unique_sites.clear();
            std::copy(walk.begin(), walk.begin()+l-1, std::back_inserter(unique_sites));
            std::sort(unique_sites.begin(),unique_sites.end());
            it = std::unique(unique_sites.begin(),unique_sites.end());
            unique_sites.resize( std::distance(unique_sites.begin(),it) );
            
            unique_l = unique_sites.size();
            
        } else if(type.compare("RPC") == 0){
            l = RPC(n, N, fill, l0, grid_al, grid, walk);

            //get number of unique visited sites
            unique_sites.clear();
            std::copy(walk.begin(), walk.begin()+l-1, std::back_inserter(unique_sites));
            std::sort(unique_sites.begin(),unique_sites.end());
            it = std::unique(unique_sites.begin(),unique_sites.end());
            unique_sites.resize( std::distance(unique_sites.begin(),it) );
            
            unique_l = unique_sites.size();
              
        } else if(type.compare("MSAW") == 0){
            l = MSAW(n, N, fill, l0, grid_al, grid, walk);
            unique_l = l-1;
            
        } else if(type.compare("RRAY") == 0){
            int D = (int) grid_al[0].size()/2;
            int L = (int) std::round(std::pow(grid_al.size(),1./D));
            l = RRAY(n, N, fill, L, D, grid, walk);
            unique_l = l-1;
            
        } else
            break;
            
        if(l<0)
            continue;
        
        //add edge
        if(n>0)
            vws.push_back({n,grid[walk[l-1]]});
        
        //update grid
        for(int s=0;s<l-1;s++) grid[walk[s]]=n; //leave out the junction point
        fill += unique_l;
        
        //update walks
        if(n==0)
            walks.push_back({walk.begin(),walk.begin()+l-1});
        else
            walks.push_back({walk.begin(),walk.begin()+l});
        
        //std::cout<<"fill: "<<fill<<std::endl<<std::endl;
        n++;
        
    }

    
    //free memory    
    delete [] grid;

    return n;
}


// Exposing functions to python


static PyObject* Full(PyObject* self, PyObject* args){
    int l0;
    
    PyObject* seq;
    int seqlen, edgeslen;
    std::vector<std::vector<int>> grid_al;
    char *type;
    
    if (!PyArg_ParseTuple(args, "sOi", &type, &seq, &l0))
        return NULL;       
    
    seq = PySequence_Fast(seq, "argument must be iterable");
    if(!seq)
        return NULL;
    
    seqlen = PySequence_Fast_GET_SIZE(seq);

    
    for(int i=0; i < seqlen; i++) {
        PyObject *edges  = PySequence_Fast_GET_ITEM(seq, i);
        edges = PySequence_Fast(edges, "argument must be iterable");
        edgeslen = PySequence_Fast_GET_SIZE(edges);
        
        grid_al.push_back({});
        
        for(int j=0; j < edgeslen; j++) {
            PyObject *litem;
            PyObject *item = PySequence_Fast_GET_ITEM(edges, j);
            
            litem = PyNumber_Long(item);
            
            grid_al[i].push_back(PyLong_AS_LONG(litem));
            Py_DECREF(litem);
        
        }
        Py_DECREF(edges);
    }    

    Py_DECREF(seq);
    //end of parsing adjacency list    
    
    std::vector<std::array<int,2>> vws;
    std::vector<std::vector<int>> walks;

    Full_C(type, grid_al, l0, vws, walks);

    PyObject* python_vws = PyList_New(vws.size());
    for (int e = 0; e<int(vws.size()); e++)
        PyList_SetItem(python_vws, e, Py_BuildValue("(ii)", vws[e][0], vws[e][1]));
    
    PyObject* python_walks = PyList_New(walks.size());
    for (int i = 0; i<int(walks.size()); i++){
        PyObject* python_walk = PyList_New(walks[i].size());
        for (int j = 0; j<int(walks[i].size()); j++)
            PyList_SetItem(python_walk, j, PyLong_FromLong(walks[i][j]));
        PyList_SetItem(python_walks, i, python_walk);
        }
        
    PyObject* return_tuple = Py_BuildValue("(OO)",python_vws,python_walks);
    
    Py_DECREF(python_vws);
    Py_DECREF(python_walks);
    
    return return_tuple;
}


static PyMethodDef mainMethods[] = {
    {"Full",Full, METH_VARARGS,"Generate full rw net."},
    {NULL,NULL,0,NULL}
};

static PyModuleDef RWModelsCPP = {
    PyModuleDef_HEAD_INIT,
    "RWModelsCPP","RWModelsCPP",
    -1,
    mainMethods
};

PyMODINIT_FUNC PyInit_RWModelsCPP(void){
    return PyModule_Create(&RWModelsCPP);
}


int main(){
    
    return 0;
}


