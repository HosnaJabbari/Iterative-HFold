#include "sparse_tree.hh"
#include <iomanip>
#include <vector>

#define maxSize 14 // 2^14 for sparse table





sparse_tree::sparse_tree(std::string structure,int n){
    this->n = n;
    this->structure = structure;
    tree.resize(n+1,Node());
    tree[0] = Node(0);
    FAI.resize(2*(n+1),-1);
    level.resize((n+1),-1);
    up.resize(n+1);
    logn.resize(2*(n+1));
    create_tree(n,structure);
    preprocess();
    ptr = 0;

    int cur = 1;
    dfs(0,-1,0);
    // depthArr.resize(euler.size());
    makeArr();
    level.clear();
    sparse_table.resize(euler.size(),std::vector<int>(maxSize,-1));
    buildSparseTable(euler.size());
}

sparse_tree::~sparse_tree(){
}


// Create Level depthArray corresponding
// to the Euler walk Array
void sparse_tree::makeArr()
{
    for (auto x : euler) depthArr.push_back(level[x]);
}
    
/**
 *  Query finds the node which is at the maximum height within a section of the depth array. 
 * As the parent of the given l and r must be between the two, we can search between.
 * By using the sparse table which holds the the maximum height node for 2^dx ahead of the node,
 * we can cover the full area between l and r and find the max height node
 */
const int sparse_tree::query(int l,int r) const{
    int d = r-l;
    int dx = logn[d];
    if (l==r) return l;

    if (depthArr[sparse_table[l][dx]] > depthArr[sparse_table[r-p2[dx]][dx]])
        return sparse_table[r-p2[dx]][dx];
    else
        return sparse_table[l][dx];
}


/**
 * Return the least common ancestor for the two indices we are given
 * From query, we get the node with the maximum height between the two (the parent)
 * and we find the actual position in the dfs
*/
const int sparse_tree::LCA(int i,int j) const{
            // trivial case
    if (i==j) return i;

    if (FAI[i] > FAI[j]) std::swap(i,j);
    // doing RMQ in the required range
    return euler[query(FAI[i], FAI[j])];
}   

/**
 * Fill the array p2 with the powers of 2 so that we do not need to recalculate them
*/
void sparse_tree::preprocess(){
    // memorizing powers of 2
    p2[0] = 1;
    for (int i=1; i<maxSize; i++)
        p2[i] = p2[i-1]*2; 

    // memorizing all log(n) values
    int val = 1,ptr=0;
    for (int i=1; i<2*n; i++){
        logn[i] = ptr-1;
        if (val==i){
            val*=2;
            logn[i] = ptr;
            ptr++;
        }
    }          
}
/**
 * Create a tree from a structure
 * A stack is used to hold the opening base pairs indices
 * Every iteration we create a node with an index of i
 * We assign it a pointer to its parent if it has one
 * If the index is a closing pair, we set the pair class variable for the nodes and we pop from the stack
 * If the index is an opening pair, we push it back as one of the children of the parent and push the index into the stack
 * If the index is an x, the base cannot pair and we change the pair value to -1
*/
void sparse_tree::create_tree(int n, std::string structure){

    // std::vector<Node> stack;
    std::vector<int> stackI;
    int count = 0;
    stackI.push_back(0);
    for(int i = 1; i<=n;++i){
        
        if(structure[i-1] == 'x') tree[i].pair = -1;
        
        if(structure[i-1] == ')'){
            tree[i].pair = stackI.back();
            tree[stackI.back()].pair = i; 
            stackI.pop_back();
            count = 0;
        }

        tree[i].index = i;
        tree[i].parent = &tree[stackI.back()];
        tree[stackI.back()].bases.push_back(i);  
        
        
        if(structure[i-1]== '('){
            tree[stackI.back()].children.push_back(i); // only push back base pairs -> pushes back the start of the pair
            stackI.push_back(i);
            count = 0;
        }
        up[i] = count;
        ++count;
        
    }

}

/**
 * Euler Walk ( preorder traversal)
 * converting tree to linear depthArray
 * Time Complexity : O(n)
 * */
void sparse_tree::dfs(int cur, int prev, int dep){
    
    if (FAI[cur]==-1) FAI[cur] = ptr;
    level[cur] = dep;
    euler.push_back(cur);
    ptr++;
    for (auto x:tree[cur].bases){
        if (x != prev){
            dfs(x,cur,dep+1);
            euler.push_back(cur);
            ptr++;
        }
    }
}

/**
 * Build the sparse table which allows for constant time queries through precomputed saved values
*/
void sparse_tree::buildSparseTable(int n){

    // filling base case values
    for (int i=1; i<n; i++)
        sparse_table[i-1][0] = (depthArr[i]>depthArr[i-1])?i-1:i;

    //dp to fill sparse table
    for (int l=1; l<maxSize; l++){
        for (int i=0; i<n; i++){
            if (sparse_table[i][l-1]!=-1 && sparse_table[i+p2[l-1]][l-1]!=-1)
            sparse_table[i][l] =
                (depthArr[sparse_table[i][l-1]]>depthArr[sparse_table[i+p2[l-1]][l-1]])?
                sparse_table[i+p2[l-1]][l-1] : sparse_table[i][l-1];
            else
                break;
        }
    }
}
    
/**
 * Returns the left innermost pair in a band between i and l
*/
const int sparse_tree::bp(int i, int l ) const{
    if(tree[l].parent->index == 0 || tree[l].pair > -1) return -2;
    if ((tree[l].parent)->index < i) return -1;
    return (tree[l].parent)->index;
}
/**
 * Returns the right innermost pair in a band between l and j
*/
const int sparse_tree::Bp(int l, int j) const{
    if(tree[l].parent->index == 0 || tree[l].pair > -1) return -2;
    if ((tree[l].parent)->pair > j) return -1;
    return (tree[l].parent)->pair;
}
/**
 * Returns the right outermostpair in a band between l and j
*/
const int sparse_tree::B(int l, int j) const{
    if(tree[l].parent->index == 0 || tree[l].pair > -1) return -2;
    if ((tree[l].parent)->pair > j) return -1;
    else{
        int lca = LCA(l,j);
        if(j == lca) return j;
        for(int x:tree[lca].children){
            if(x<j && tree[x].pair > l) return tree[x].pair;
        }
        return -100;
    }
}
// Returns the left outermost pair in a band between i and l
const int sparse_tree::b(int i, int l) const{
    if(tree[l].parent->index == 0 || tree[l].pair > -1) return -2;
    if ((tree[l].parent)->index < i) return -1;
    else{
        int lca = LCA(i,l);
        if(i==lca) return i;

        for(int x:tree[lca].children){
            if(x<l && tree[x].pair > l) return x;
        }
        return -100;
    }
}
/**
 * Returns whether there the area between i and j is weakly closed, specifically if all pairs in the [i,j] stay within [i,j]
*/
const bool sparse_tree::weakly_closed(int i, int j) const{
    if(j<i) return 0;
    if((i > tree[i].pair && tree[i].pair > 0) || tree[j].pair> j) return 0;
    if(i==j) return !(tree[j].pair > 0);
    // if((tree[i].pair > i && tree[j].pair > j) ||(tree[i].pair < i && tree[j].pair < j && tree[i].pair != -1 && tree[j].pair != -1)) return 0;
    // if(depthArr[FAI[i]] == depthArr[FAI[j]] && tree[i].parent->index == tree[j].parent->index) return 1;

    if((void*)tree[i].parent == (void*) tree[j].parent) {
        if(depthArr[FAI[i]] == depthArr[FAI[j]]) 
            return 1;
    }
    return 0;
}