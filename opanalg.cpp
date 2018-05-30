/* This file contains the main algorithm and its components */
#include <iostream>
#include <vector>
  using std::string;
  using std::vector;
#include "opanalg.hpp"
//#include "subset_sum.hpp"
#include <primecount.hpp>
#include <algorithm>


void OPAN(int d, bool verbose, bool steps, string fname){

  std::vector<Tree> forest;
  Stats s(verbose, fname);

  for (int factors = 3; factors <= d; ++factors){

    Reset(s);

    /* Create a new tree for current number of factors */
    forest.push_back(Tree());
    Tree& tree = forest.back();

    vector<Node*> leaves; /* For parallel processing*/
    vector<vector<ZZ> > exp_seqs; /* Exponent sequences for the prime sequence
                                     (The prime sequence is kept track of 
                                     by the tree) */ 

    bool growing = true; /* Flag that tells when to stop looking for prime sequences */
    std::vector<ZZ> primes; 
    
    while (growing){

      Update(s, tree);
      if (verbose) Show(s, tree);
      if (steps){
        char c;
        std::cout << "Hit enter for the next step\n";
        std::cin >> c;
      }
      

      grow(tree, RR(3)); /* Grow the tree by a single prime */
      primes = strip_primes(tree); /* Get the current prime sequence */

      if (primes.size() != factors){ /* Not enough prime factors right now */

        exp_seqs.clear(); /* Thow out old exponent sequence. Theorem that allows us to re-use exponents does not apply anymore! */
        if (b_1(primes) >= RR(2)){ /* p_1 p_2 p_3 ... p_d is already abundant and there are not enough divisors */
          primes.pop_back();
          ZZ r = min_deficient(primes); /* find minimum prime so that Y(primes + r, {1,1,1,1...}) is deficient */
          replace(tree, r);
        }

        continue;

      }

      if (b_inf(primes) <= RR(2)){ /*Deficient no matter what exponents are input */
        growing = fail(tree); /* Checks if there are more prime sequences to find */
        continue;
      }

      if (exp_find(primes, exp_seqs)){
        //leaves.push_back(tree.curr); /* For parallel proccessing */

        /* Applies the efficiency theorem */
        Write(s, primes, exp_seqs); 
        efficiency(primes, exp_seqs, s);
        Write(s, primes, exp_seqs);

        replace(tree, primes.back()); /* Result of the efficiency theorem */

        success(tree);
        Write(s, primes, exp_seqs);
        //fail(tree); /* For parallel proccessing, may seem out of place */
      } else {
        if (!cap_check(primes, forest, factors)){ /* Calculate cap_{factors} (P), determine if tree should stop growing */
          growing = fail(tree);
        }
      }
    }

    Write(s, primes, exp_seqs);
    //expand_sets(leaves, s); /* This line of code is used for paralell proccessing */
   
  }
  return;
}

bool exp_find(vector<ZZ>& primes, vector<vector<ZZ> >& exp_seqs){
  /* First, test all exponents from previous primes,
   * taking advantage of the Efficiency Lemma */

  if (!exp_seqs.empty()){
    int i;

    for (i = 0; i < exp_seqs.size(); ++i){
      vector<ZZ>& exps = exp_seqs[i];
      if (b(primes, exps) < RR(2)){ break; }            
      if (!primitive(primes, exps)) { break; }
    }
    if (i != exp_seqs.size()) {

      /* Efficiency Lemma had a false hypothesis, find exponents manually */
      exp_seqs = expAlg(primes);
      return (!exp_seqs.empty());

    } else {
      /* Efficiecy Lemma had a true hypothesis */
      return true;
    }
  }
  /* There were no given exponents to check. Find exponents manually*/
  exp_seqs = expAlg(primes);
  return (!exp_seqs.empty());
}

bool cap_check(vector<ZZ>& primes, vector<Tree>& forest, int& factors){

  /* Run backup algorithm on tree, retrieve last truncated prime */
  ZZ truncp = backup(forest.back());
  primes = strip_primes(forest.back());

  if (primes.size() == factors - 1){ /* Follows from the continuity of OPANS theorem */
    return false;
  }

  /* Find cap_d(P) */
  ZZ cap = findmax(primes, forest);

  /* Check theorem concering capd(P) */
  if (truncp <= cap) { 
    grow(forest.back(), conv<RR>(truncp + 1));
    return true;
  } else {
    /* Theorem tells us that we should modify the primes and cannot look here anymore */
    return false;
  }
}


void efficiency(vector<ZZ>& primes, vector<vector<ZZ> >& exp_seqs, Stats& s){
  /* Assumes that input primes work with exp_seqs. */
  
  ZZ increment(2);
  ZZ bad(0);
  bool comedown = false;

  while (increment >= ZZ(2)){//(last_exps == exp_seqs){        
    vector<ZZ> new_primes = primes;
    vector<vector<ZZ> > last_exps = exp_seqs;

    ZZ prev = primes[primes.size() - 1];
    ZZ next = NextPrime(prev + increment);
    new_primes.pop_back();
    new_primes.push_back(next);
 
    if (bad != next){ 
      exp_find(new_primes, exp_seqs);
    } else {
      exp_seqs.clear();
    }

    if (last_exps == exp_seqs){

      if (!comedown){ increment *= ZZ(4); }
      
      primes = new_primes; // Move on to new primes

    } else {
             
      comedown = true;
      increment /= ZZ(4);
      bad = next;
      exp_seqs = last_exps;
    }      
  }



}


void Reset(Stats& s){

  std::cout << "Reset\n";

  s.number_found = 0;
  s.product = 0;
  s.prev_core.clear();
  s.init_core.clear();
  s.prev_exps.clear();
  s.init_tail == 0;
  s.prev_tail == 0;
}

void Update(Stats& s, Tree& t){

  std::cout << "Update\n";

}

void Show(Stats& s, Tree& t){
  if (s.verbose){
    std::cout << "Show\n";
  }
}

void dump_primes(Stats& s){
  printvectos(s.prev_core, s.out);//s.out);
}

void Write(Stats& s, vector<ZZ>& primes, vector<vector<ZZ> >& exp_sets){
 
  
  std::cout << "Write\n";

  bool echange = false;
  vector<ZZ> core_primes = primes;
  core_primes.pop_back();

  if (exp_sets != s.prev_exps){ 
    echange = true;
  } else if (core_primes != s.prev_core){
    echange = true; /* Dump primes in form p1^e1 p2^e2 p3^e3 p4^e4 q^e5 for e's : and q's: */
  
  }

  
 
  /* If there was a change in exponents, 
   * we need to display all prime sequences
   * with old exponents */

  if (echange) {

    /* Just print one prime */
    if (s.prev_tail == s.init_tail){
      printwithprime(s.prev_core, s.prev_tail, s.out);
      printexponents(s.prev_exps, s.out);
      s.out << '\n';
    } else {
      printwithgeneric(s.prev_core, s.out);
      s.out << "q: " << s.init_tail << " -> " << s.prev_tail << '\n';
      printexponents(s.prev_exps, s.out);
      s.out << '\n';
    }

    s.init_tail = primes.back();
    s.init_core = core_primes;

  }
  s.prev_tail = primes.back();
  s.prev_core = core_primes;
  s.prev_exps = exp_sets;
  
}


