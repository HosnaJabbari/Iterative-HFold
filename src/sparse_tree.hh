#ifndef SPARSE_TREE
#define SPARSE_TREE

#include <iomanip>
#include <vector>
#include <string>

#define maxSize 14 // 2^14 for sparse table

class Node {

    public:

        // Pointer to parent
        Node* parent = nullptr;

        int index = -1;

        int pair = -2;

        // vector of the positions of the children (all)
        std::vector<int> bases;
        // vector of the positions of the children (base pairs only)
        std::vector<int> children;

        Node(int index){
            this->index = index;
        }
        Node(){

        }



};

class sparse_tree{

    public:

        sparse_tree(std::string structure,int n);
        ~sparse_tree();

        // const int maxSize = 14; 
        std::vector<Node> tree; // A vector of Nodes corresponding to each index in the structure
        std::vector<int> FAI; // The index of the First appearance of a node
        std::vector<int> level; // Stores depths for each node of the tree
        std::vector<int> euler; // euler walk
        std::vector<int> depthArr; // depths corresponding to euler
        std::vector<int> logn; // holds logn values
        std::vector<int> up; // vector holding unpaired bases
        uint16_t n;
        std::string structure;
        int ptr; // Pointer to euler walk
        // uint16_t** sparse_table;
        std::vector< std::vector<int> > sparse_table;
        int p2[maxSize];

        int bp(int i, int l );
        int Bp(int l, int j);
        int B(int l, int j);
        int b(int i, int l);
        bool weakly_closed(int i, int j);


    private:
        void makeArr();
        int query(int l,int r);
        int LCA(int i,int j);
        void preprocess();
        void create_tree(int n, std::string structure);
        void dfs(int cur, int prev, int dep);
        void buildSparseTable(int n);


};

#endif