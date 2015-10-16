//#include "phylolm.h"
#include <math.h>
#include <stdlib.h>
#include <R.h>

//TODO re-write this function in C++ using RCpp; using calloc/free/pointers is not a good programming style.
void effectiveSampleSize (int *Npo, int *npo, int *pNpo, int *rootpo, double *transa, double *transb, int *des, int *anc, int *edge, double *output){
    int N = *Npo;     // number of edges, excluding root edge
    int n = *npo;     // number of tips
    int pN=*pNpo;     // number of internal nodes
    int npN = n + pN; //      total # nodes, including root node. N+1 if binary tree.
    int r = *rootpo;  // root node ID (in R)
    r--;
    double rootEdge= *transa; // length of root edge
    // transb: vector of edge lengths
    // des:    vector of descendents (one for each edge)
    // anc:    vector of ancestor nodes
    // edge:   vector of indices (for R) of edges to cut
    // !WARNING! edge indices are assumed to be sorted (ascending):
    //           to visit them in order during the post-order traversal
    // output: vector of size #edges to cut: for each subtree
    //         n_e = max(V) * one' V^{-1} one, where V=BM covariance for subtree
    //         so that max(V) = subtree height

    double* Theight=(double*)calloc(npN, sizeof(double));
    // Theight[i] = tree height of subtree rooted node i (c-style: = node i+1 R-style)
    // fixit: later, extract those
    // of subtree at the root if i=0 (including root edge!)
    //            at child node of edge i otherwise
    double* vec11=(double*)calloc(npN, sizeof(double));
    int* zero =(int*)calloc(npN, sizeof(int));
    // zero[i-1] will be -1 if node i has no children edge of length 0
    //                    d if node i has exactly 1 child edge of length 0,
    //                      to child node d+1.
    // error returned if some node has 2 or more child edges of length 0
    for (int iedge=0; iedge<N+1; iedge++)
        zero[iedge] = -1;
    int nextE = 0; //index of next edge to cut, in cut edge vector

    // loop over all N+1 edges, with root edge last. Assumes postorder traversal.
    for (int iedge=0; iedge < N+1; iedge++){
        double len;         // edge length
        int di, anci=0;
        if (iedge<N){      // non-root edge
            len=transb[iedge];
            di= des[iedge]-1;// subtree at node di+1.
            anci= anc[iedge]-1;
            if (Theight[anci] <= 0.0) Theight[anci] = Theight[di]+len;
        } else {           // iedge=N -> root edge
            len = rootEdge;
            di = r; // but anci meaningless
        }

        if (iedge == (edge[nextE]-1) ){ // edge to cut: get n_e and do *not* contribute to parent node
            if (di<n || Theight[di] == 0.0) // external edge or polytomous tips
                output[nextE] = 1.0;
            else
                output[nextE] = Theight[di] * vec11[di];
            nextE++;
        }

        else { // edge *not* to cut
            if (di<n){ // external edge
                if (len>0)
                    vec11[di] = 1/len;
                else {
                    if (zero[anci] >= 0) // anci already found with 1 child edge of length 0.
                        error("two or more sister external edges have length 0, V is singular\n (node %d in pruning-wise-ordered tree)", anci+1);
                    else
                        zero[anci] = di; // which is >= 0
                }
            }
            else { // internal edge. contributions from descendants of di have already been collected.
                int goodchildren = 1;
                if (zero[di] >= 0) { // exactly 1 child edge of length 0, to descendant d0=zero[di]
                    if (len<=0)
                        error("One external edge and its parent both have length 0\n (node %d in pruning-wise-ordered tree). To avoid this situation,\n please make a polytomy or resolve it differently ",di+1);
                    goodchildren = 0;
                    // still assuming di has more than 1 child, with current vec11[di]>0.
                }
                // rescaling to get vec11 correct ("p" instead of "p_A"),
                if (goodchildren)
                    vec11[di] /= (1+ len * vec11[di]);
                else
                    vec11[di] = 1/len;
            }
            // next: collect contribution of iedge to its parent, i.e. to ancestor node anci.
            // *except* if root edge (anci meaningless) or if external edge of zero length
            if ((iedge < N) && ((di>=n) || (len>0))){
                vec11[anci] += vec11[di];
            }
        }
    }

    free(Theight);
    free(vec11);
    free(zero);
}
