#include <Rcpp.h>

using namespace Rcpp;



class Offset3D {
private:
  int nrows, ncols, nmats ;

public:
  Offset3D( int nrows_, int ncols_, int nmats_) : nrows(nrows_), ncols(ncols_), nmats(nmats_){}
  Offset3D( IntegerVector ns) : nrows(ns[0]), ncols(ns[1]), nmats(ns[2]){}
  int operator()( int i, int j, int k){
    return i + j * nrows + k * ( nrows * ncols ) ;
  }

} ;

class Offset4D {
private:
  int nrows, ncols, nmats, ncubes;

public:
  Offset4D( int nrows_, int ncols_, int nmats_, int ncubes_) : nrows(nrows_), ncols(ncols_), nmats(nmats_), ncubes(ncubes_){}
  Offset4D( IntegerVector ns) : nrows(ns[0]), ncols(ns[1]), nmats(ns[2]), ncubes(ns[3]){}
  int operator()( int i, int j, int k, int l){
    return i + j * nrows + k * ( nrows * ncols ) + l * (nrows * ncols * nmats);
  }

} ;



// [[Rcpp::export]]
NumericVector calculate_mi_and_omes(int i, int j, NumericVector P2, NumericMatrix P1, int q, int N) {
  Offset4D offset(N, N, q, q);
  double M = 0;
  double O = 0;
  for(int alpha = 0; alpha < q; alpha++) {
    for(int beta = 0; beta < q; beta++) {
      if( P2(offset(i, j, alpha, beta)) > 0 ) {
        M = M + P2(offset(i, j, alpha, beta)) * log(P2(offset(i, j, alpha, beta)) / P1(i, alpha) / P1(j, beta));
        double p = (P2(offset(i, j, alpha, beta)) - (P1(i, alpha) * P1(j, beta)));
        O = O + p*p;
      }
    }
  }

  return NumericVector::create(M, O);
}

int mapkey(int i, int alpha, int q) {
  int A = (q-1) * i + alpha;
  return A;
}

// [[Rcpp::export]]
NumericMatrix Compute_C(NumericVector Pij, NumericMatrix Pi, int N, int q) {
  NumericMatrix C(N*(q-1), N*(q-1));
  Offset4D offset(N, N, q, q);


  for(int i = 0; i < N; i++) {
    for(int j = 0; j < N; j++) {
      for(int alpha = 0; alpha < (q - 1); alpha++) {
        for(int beta = 0; beta < (q - 1); beta++) {
          C(mapkey(i, alpha, q), mapkey(j, beta, q)) = Pij(offset(i, j, alpha, beta)) - Pi(i, alpha) * Pi(j, beta);
        }
      }
    }
  }
  return C;
}


// [[Rcpp::export]]
NumericMatrix hammingdist(NumericMatrix X) {
  int n = X.nrow();
  int c = X.ncol();
  NumericMatrix m(n, n);

  for(int i = 0; i < n-1; i++) {
    for(int j = i+1; j < n; j++) {
      for( int k = 0; k < c; k++) {
        if(X(i, k) == X(j, k)) {
          m(i, j)++;
        }
      }
      m(i, j) = 1 - m(i, j)/c;
      m(j, i) = m(i, j);
    }
  }
  return m;
}

// [[Rcpp::export]]
NumericMatrix AccPi_true(NumericMatrix align, NumericVector W, int N, int q, int M) {
  NumericMatrix Pi_true(N, q);

  for(int j = 0; j < M; j++) {
    for(int i = 0; i < N; i++) {
      // Minus 1 because of th 1-indexing scheme of R and zero-indexing scheme in CPP
      Pi_true(i, align(j,i)-1) += W[j];
    }
  }

  return Pi_true;
}

// [[Rcpp::export]]
NumericVector AccPij_true(NumericMatrix align, NumericVector W, int N, int q, int M) {
  NumericVector Pij_true(N*N*q*q);
  Offset4D offset(N, N, q, q) ;
  Pij_true.attr("dim") = IntegerVector::create(N, N, q, q);
  for(int l = 0; l < M; l++) {
    for(int i = 0; i < N-1; i++) {
      for(int j = i + 1; j < N; j++) {
        // Minus 1 because of th 1-indexing scheme of R and zero-indexing scheme in CPP
        Pij_true(offset(i, j, align(l,i) - 1, align(l,j) - 1)) =
          Pij_true(offset(i, j, align(l, i) - 1, align(l, j) - 1)) + W[l];
        Pij_true(offset(j, i, align(l,j) - 1, align(l,i) - 1)) =
          Pij_true(offset(i, j, align(l, i) - 1, align(l, j) - 1));
      }
    }
  }
  return Pij_true;
}

// [[Rcpp::export]]
NumericVector arrayC(NumericVector input, IntegerVector dim) {
  input.attr("dim") = dim;
  return input;
}

