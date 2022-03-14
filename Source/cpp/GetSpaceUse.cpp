// [[Rcpp::depends(RcppProgress)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <progress.hpp>
#include <progress_bar.hpp>
#include <unordered_map>

using namespace arma;
using namespace Rcpp;

// [[Rcpp::export]]
// IMPORT QUANTILE FUNCTION FROM R
// NumericVector quantileCpp(NumericVector x, NumericVector probs) {
//   Environment stats("package:stats");
//   Function quantile = stats["quantile"];
//   int npr = probs.size();
//   NumericVector ans(npr);
//   for(int i=0; i<npr; i++){
//     ans[i] = as<double>(quantile(x, probs[i]));
//   }
// return ans;
// }

NumericVector quantileCpp(NumericVector x, NumericVector q) {
  NumericVector y = clone(x);
  std::sort(y.begin(), y.end());
  return y[x.size()*(q - 0.000000001)];
}


// [[Rcpp::export]]
// IMPORT MODE FUNCTION stackoverflow.com/questions/55212746/rcpp-fast-statistical-mode-function-with-vector-input-of-any-type
int fastIntMode(NumericVector x, bool narm = false) {
  if (narm) x = x[!is_na(x)];
  int myMax = 1;
  int myMode = 0;
  std::unordered_map<int, int> modeMap;
  modeMap.reserve(x.size());
  
  for (std::size_t i = 0, len = x.size(); i < len; ++i) {
    auto it = modeMap.find(x[i]);
    
    if (it != modeMap.end()) {
      ++(it->second);
      if (it->second > myMax) {
        myMax = it->second;
        myMode = x[i];
      }
    } else {
      modeMap.insert({x[i], 1});
    }
  }
  
  return myMode;
}



// [[Rcpp::export]]
List  GetSpaceUse(NumericMatrix sx,                   // X COORDINATES 
                     NumericMatrix sy,           // Y COORDINATES 
                     NumericMatrix z,            // Z STATE 
                     NumericVector sigma,        // SIGMA VALUES 
                     NumericMatrix habitatxy,    // HABITAT COORDINATES
                     NumericVector aliveStates,  // ALIVE STATES 
                     NumericMatrix regionID,     // MATRIX WITH REGION ID, ONE ROW PER REGION WITH 1 AND 0 WHETHER IT BELONGS TO THE REGION OR NOT NEED TO PROVIDE ROWNAMES TO IDENTIFY REGIONS.
                     NumericVector probs = NumericVector::create(0.025,0.975), //CI TO BE RETURNED
                     bool display_progress = true, // DISPLAY PROGRESS BAR
                     bool returnPosteriorCells = true  // IF POSTERIORS SHOULD BE RETURNED
) {
  
  //INITIATE OBJECTS 
  int niter = sx.nrow(), nind = sx.ncol();
  int ncells = regionID.ncol(), nregions = regionID.nrow();
  double xMin = min(habitatxy(_,0)), xMax = max(habitatxy(_,0));
  double yMin = min(habitatxy(_,1)), yMax = max(habitatxy(_,1));
  
  NumericMatrix pTot(ncells, niter);
  
  //INITIATE PROGRESS BAR
  Progress prog(niter, display_progress);
  
  //LOOP FOR EACH MCMC ITERATION  
  for(int ite = 0; ite < niter; ite++){
    //UPDATE PROGRESS BAR
    if(Progress::check_abort())return -1.0;
    prog.increment(); 
    
    //LOOP FOR EACH INDIVIDUAL
    NumericVector p(ncells);
    NumericVector dist2(ncells);
    NumericVector pTemp(ncells);
    
    for (int i = 0; i < nind; i++) {
      //CHECK IF ALIVE OR NOT
      if(std::find(aliveStates.begin(), aliveStates.end(), z(ite,i))!=aliveStates.end()) {
        if(sx(ite,i) >= xMin && sx(ite,i) <= xMax && sy(ite,i) >= yMin && sy(ite,i) <= yMax) {
          double sigma1 = (2*sigma(ite,i)*sigma(ite,i));
          //LOOP FOR EACH HABITAT CELL
          for (int j = 0; j < ncells; j++) {
            //CALCULATE DISTANCE^2
            dist2(j) = (pow(sx(ite,i) - habitatxy(j, 0),2) + pow(sy(ite,i) - habitatxy(j, 1),2)) ;
          }
          //CALCULATE HALFNORMAL DENSITY
          pTemp = exp((-dist2)/sigma1);
          //NORMALIZE TO GET SPACE-USE
          p += pTemp/sum(pTemp);
        }
      }
    }
    pTot(_,ite) = p;
  }
  
  //CELL STATISTICS 
  NumericVector outMean(ncells);
  NumericVector outMedian(ncells);
  NumericVector outsd(ncells);
  NumericVector outCV(ncells);
  NumericVector outCIL(ncells);
  NumericVector outCIH(ncells);
  
  for (int c = 0; c < ncells; c++){
    NumericVector tmpMean = pTot(c,_);
    outMean(c) = mean(tmpMean);
    outMedian(c) = median(tmpMean);
    outsd(c) = sd(tmpMean);
    outCV(c) = (100*outsd(c))/outMean(c);
    NumericVector outCILr= quantileCpp(tmpMean ,probs);
    outCIL(c) = outCILr(0);
    outCIH(c) = outCILr(1);
  }
  
  //================
  //==== RETURN ====
  //================
  //==== IF REGION STATISTICS ARE ASKED ==== 
  //INITIATE OBJECTS
  NumericMatrix subsetSum(nregions, niter);
  NumericMatrix summary ((nregions+1),5);
  //ATTRIBUTE NAMES FOR THE SUMMARY TABLE 
  CharacterVector Names =  rownames(regionID);
  //ADD A TOTAL ROW
  CharacterVector Names1 =  CharacterVector::create("Total");
  CharacterVector Names2(Names.size() + Names1.size());
  std::copy(Names.begin(), Names.end(), Names2.begin());
  std::copy(Names1.begin(), Names1.end(), Names2.begin() + Names.size());
  rownames(summary) = Names2;
  colnames(summary) = CharacterVector::create("mean", "median", "mode", "95%CILow","95%CIHigh");
  
  //SUMMARY AND POSTERIOR FOR EACH REGIONS 
  for (int r = 0; r < nregions; r++){
    //SUBSET TO CELL WITHIN REGIONS 
    NumericVector T = regionID(r,_);
    mat Xmat(pTot.begin(), pTot.nrow(), pTot.ncol(), false);
    colvec tIdx(T.begin(), T.size(), false); 
    mat subMat = Xmat.rows(find(tIdx == 1));
    int nsub = sum(T);
    // SUM THE ABUNDANCE OF ALL CELLS OF THE REGION FOR ALL ITERATIONS 
    for (int ite = 0; ite < niter; ite++){
      double total = 0;
      for (int j = 0; j < nsub; j++) {
        total += subMat(j,ite);
      }
      subsetSum(r,ite) = total;
    }
    //FILL IN THE SUMMARY TABLE
    NumericVector tmp = subsetSum(r, _);
    summary(r,0) = mean(tmp);
    summary(r,1) = median(tmp);
    summary(r,2) = fastIntMode(tmp);
    NumericVector outCILr= quantileCpp(tmp, probs);
    summary(r,3) = outCILr(0);
    summary(r,4) = outCILr(1);
  }
  
  //SUMMARY AND POSTERIOR FOR ALL REGIONS
  //SUBSET TO CELLS WITHIN ALL REGIONS 
  NumericVector AllRegionsID= colSums(regionID);
  mat Xmat1(pTot.begin(), pTot.nrow(), pTot.ncol(), false);
  colvec tIdx1(AllRegionsID.begin(), AllRegionsID.size(), false); 
  mat subMat1 = Xmat1.rows(find(tIdx1 > 0));
  int nsub1 = sum(AllRegionsID);
  
  //SUM THE ABUNDANCE FOR ALL REGIONS 
  NumericVector subsetSum1(niter) ;
  for (int ite = 0; ite < niter; ite++){
    double total = 0;
    for (int j = 0; j < nsub1; j++) {
      total += subMat1(j, ite);
    }
    subsetSum1(ite) = total;
  }
  
  //FILL IN THE "Total" SUMMARY TABLE
  summary(nregions,0) = mean(subsetSum1);
  summary(nregions,1) = median(subsetSum1);
  summary(nregions,2) = fastIntMode(subsetSum1);
  NumericVector outCILr1 = quantileCpp(subsetSum1, probs);
  summary(nregions,3) = outCILr1(0);
  summary(nregions,4) = outCILr1(1);
  
  //OUTPUT
  if(returnPosteriorCells){
    return List::create(Named("MeanCell") = outMean,
                        Named("MedianCell") = outMedian,
                        Named("SDCell") = outsd,
                        Named("CVCell") = outCV,
                        Named("CILCell") = outCIL,
                        Named("CIHCell") = outCIH,
                        Named("summary") = summary,
                        Named("PosteriorRegions") = subsetSum,
                        Named("PosteriorCells") = pTot);
  }else{
    return List::create(Named("MeanCell") = outMean,
                        Named("MedianCell") = outMedian,
                        Named("SDCell") = outsd,
                        Named("CVCell") = outCV,
                        Named("CILCell") = outCIL,
                        Named("CIHCell") = outCIH,
                        Named("summary") = summary,
                        Named("PosteriorRegions") = subsetSum);
  }
}

