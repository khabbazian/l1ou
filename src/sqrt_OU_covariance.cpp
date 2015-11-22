#include <algorithm>
#include <Rcpp.h>

using namespace Rcpp;
typedef Matrix<REALSXP> NumericMatrix;

#define RASSERT(condition){if(!(condition)){throw std::range_error(std::string("internal error!@")+__FILE__);}}

// [[Rcpp::plugins(cpp11)]]
void one_step(const int i1, const int i2, const int e1, const int e2,
        const int counter, const int nTips,
        Rcpp::NumericMatrix &edgeList, //the third column contains lengths
        Rcpp::NumericVector &tips, 
        Rcpp::NumericMatrix &F, Rcpp::NumericMatrix &G, 
        Rcpp::NumericMatrix &D, Rcpp::NumericMatrix &B,
        double &rootEdge){

    const int nEdges = edgeList.nrow();

    //int e1=-1, e2=-1;
    //for (int i=0; i<nEdges; ++i)
    //    if (edgeList(i,1) == i1)
    //        e1 = i;
    //    else if (edgeList(i,1) == i2)
    //        e2 = i;

    RASSERT( e1!=-1 && e2!=-1 );

    const int i3  = edgeList(e1,0);
    RASSERT( edgeList(e2,0) == i3 );
    
    const double t1 = edgeList(e1,2);
    const double t2 = edgeList(e2,2);

    int e3 = -1; // -1 means root index
    //for(int i=0; i<nEdges; ++i) 
    for(int i=e2+1; i<nEdges; ++i) 
        if( edgeList(i,1) == i3){
            e3 = i;
            break;
        }

    double t3;
    if (e3 == -1)
        t3 = rootEdge>0?rootEdge:0;
    else
        t3 = edgeList(e3,2);
    
    const double  u = t1+t2;
    const double us = std::sqrt(u);

    D(_,counter) = (F(_,i1) - F(_,i2))/us;
    B(_,counter) = (G(_,i1)*t1 - G(_,i2)*t2)/us;

    F(_,i3) = (F(_,i1)*t2 + F(_,i2)*t1)/u;
    G(_,i3) = G(_,i1) + G(_,i2);

    if( e3>=0 )
        edgeList(e3,2) = t3 + 1/(1/t1+1/t2);
    else
        rootEdge = t3 + 1/(1/t1+1/t2);

    // erase-remove idiom
    tips.erase( std::remove( tips.begin(), tips.end(), i1 ), tips.end() );
    tips.erase( std::remove( tips.begin(), tips.end(), i2 ), tips.end() );
    tips.push_back(i3);
}



// [[Rcpp::export]]
Rcpp::List cmp_sqrt_OU_covariance(Rcpp::NumericMatrix edgeList, int nTips, double rootEdge){

    //TODO assert( edgeList.ncol == 3);

    Rcpp::NumericMatrix F(nTips, 2*nTips-1);
    Rcpp::NumericMatrix G(nTips, 2*nTips-1);
    for(int i=0; i<nTips; ++i)
        F(i,i) = G(i,i) = 1;

    Rcpp::NumericMatrix D(nTips,nTips);
    Rcpp::NumericMatrix B(nTips,nTips);

    Rcpp::NumericVector tips(nTips);
    //std::iota(tips.begin(),tips.end(),1);
    for(int i=0; i<tips.size(); ++i)
        tips(i) = i+1;

    int counter = 0;
    for(int i=0; i<edgeList.nrow() && tips.size() > 1; i+=2)
        one_step( edgeList(i,1), edgeList(i+1,1), i, i+1, counter++, nTips, edgeList, tips, F, G, D, B, rootEdge);
    
    for(int i=0; i<F.nrow(); ++i){
        D(i,counter) = F(i,tips[0])/std::sqrt(rootEdge);
        B(i,counter) = G(i,tips[0])*std::sqrt(rootEdge);
    }

    //sqrtInvSigma # instead of D
    //sqrtSigma    # instead of B
    
    return( Rcpp::List::create( Rcpp::Named("sqrtInvSigma") = D, Rcpp::Named("sqrtSigma") = B) );
}


