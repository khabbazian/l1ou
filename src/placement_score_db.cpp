#include <algorithm>
#include <vector>
#include <map>

#include <Rcpp.h>

// [[Rcpp::plugins(cpp11)]]

typedef std::map<std::size_t, double> DataBase;

DataBase db;
std::hash<std::string> myHash;

// [[Rcpp::export]]
void add_configuration_score_to_db(std::string str_key, double value){
    auto key = myHash(str_key);
    db[key]  = value;
}


// [[Rcpp::export]]
void erase_configuration_score_db(){
    db.erase(db.begin(), db.end());
}
 
// [[Rcpp::export]]
Rcpp::List get_score_of_configuration(std::string str_key){
    bool   valid = 0;
    double value = 0;

    auto key = myHash(str_key);
    auto itr = db.find(key);
    if(  itr != db.end() ){
        valid = 1;
        value = itr->second;
    }
    return( Rcpp::List::create( 
                Rcpp::Named("value") = value,
                Rcpp::Named("valid") = valid ) );
}



