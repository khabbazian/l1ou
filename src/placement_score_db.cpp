#include <algorithm>
#include <vector>
#include <map>

#include <Rcpp.h>

// [[Rcpp::plugins(cpp11)]]

typedef std::map<std::size_t, double> DataBase;
typedef std::pair<double,std::string>  ScoreConfig;
typedef std::vector<ScoreConfig> ScoreConfigVec;

DataBase db;
std::hash<std::string> myHash;
ScoreConfigVec myConfigVec;


// [[Rcpp::export]]
void add_configuration_score_to_db(std::string str_key, double value){
    auto key = myHash(str_key);
    db[key]  = value;
    myConfigVec.push_back( ScoreConfig(value, str_key) );
}

// [[Rcpp::export]]
Rcpp::List get_stored_config_score(){
   
    std::sort(myConfigVec.begin(), myConfigVec.end(),
            [](const ScoreConfig& lhs, const ScoreConfig& rhs) {
            return lhs.first < rhs.first; } );

    std::vector<std::string> strVec;
    std::vector<double> doubleVec;
    for(auto &itr : myConfigVec){
        doubleVec.push_back(itr.first);
        strVec.push_back   (itr.second);
    }

    return Rcpp::List::create( 
                Rcpp::Named("scores") = doubleVec,
                Rcpp::Named("configurations") = strVec ) ;
}

// [[Rcpp::export]]
void erase_configuration_score_db(){
    db.erase(db.begin(), db.end());
    myConfigVec.erase(myConfigVec.begin(), myConfigVec.end() );
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



