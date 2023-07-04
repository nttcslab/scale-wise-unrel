#include "mylib/common.hpp"
#include "mylib/graph.hpp"

#include "tdzdd/DdSpec.hpp"
#include "tdzdd/DdStructure.hpp"
#include "tdzdd/DdEval.hpp"
#include "tdzdd/spec/FrontierBasedSearch.hpp"
#include "tdzdd/util/Graph.hpp"

#include "SAPPOROBDD/include/BDD.h"

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <cassert>
#include <chrono>
#include <queue>
#include <unordered_map>
#include <unordered_set>

  
class ToBDD: public tdzdd::DdEval<ToBDD, BDD> {
  int const offset;
  
  public:
  ToBDD(int offset = 0)
  : offset(offset) {
  }
  
  bool isThreadSafe() const {
    return false;
  }
  
  void initialize(int topLevel) const {
    while (BDD_VarUsed() < topLevel + offset) {
      BDD_NewVar();
    }
  }
  
  void evalTerminal(BDD& f, int value) const {
    f = BDD(value);
  }
  
  void evalNode(BDD& f, int level, tdzdd::DdValues<BDD, 2> const& values) const {
    f = values.get(0);
    f &= ~(BDDvar(BDD_VarOfLev(level + offset)));
    if (level + offset > 0) {
      BDD f1 = values.get(1);
      f |= (BDDvar(BDD_VarOfLev(level + offset))) & f1;
    }
  }
};

class ProbEval : public tdzdd::DdEval<ProbEval, double> {
private:
  std::vector<double> prob_list_;
  
public:
  ProbEval(const std::vector<double>& prob_list) : prob_list_(prob_list) {}
  
  void evalTerminal(double& p, bool one) const { p = one ? 1.0 : 0.0; }
  
  void evalNode(double& p, int level,
                tdzdd::DdValues<double, 2> const& values) const {
    double pc = prob_list_[prob_list_.size() - level];
    p = values.get(0) * (1 - pc) + values.get(1) * pc;
  }
};

class BDDRel{
public:
  BDDRel(){};
          
  double compute(BDD f, const std::vector<double> &_wp){
    wp.reserve(_wp.size() + 1);
    std::copy(_wp.begin(), _wp.end(), wp.begin()+1);
    mc.reserve(f.Size() * 2);
    
    int v = f.Top();
    double res = computeInner(f);
    return res;
  }
  
private:
  std::vector<double> wp;
  std::unordered_map<bddword, double> mc;
  
  double computeInner(BDD f)
  {
    // base case: f is terminal
    if (f == BDD(0)) return 0.0;
    if (f == BDD(1)) return 1.0;
    assert(f.Top() != 0);
    
    // get ID, var and level
    BDD h = f;
    int v = h.Top();
    int lev = BDD_LevOfVar(v);
    bddword id = h.GetID();
    bool neg = false;
    // negative edge processing
    if (id & 1){
      h = ~h;
      id = h.GetID();
      neg = true;
    }
    // cache search
    auto it = mc.find(id);
    if (it != mc.end()){
      if (!neg) return it->second;
      else      return 1.0 - it->second;
    }
    
    double res = 0.0;
    // 0-edge processing
    {
      BDD h0 = h.At0(v);
      res = (1.0 - wp[v]) * computeInner(h0);
    }
    // 1-edge processing
    {
      BDD h1 = h.At1(v);
      res += wp[v] * computeInner(h1);
    }
    // cache entry
    mc.emplace(id, res);
    if (!neg) return res;
    else      return 1.0 - res;
  }
};

void print_usage(char *fil){
  fprintf(stderr, "Usage: %s [graph_file] [probability_file] [source_file] [order_file]\n", fil);
}

int main(int argc, char **argv){
  if(argc < 5){
    fprintf(stderr, "ERROR: too few arguments.\n");
    print_usage(argv[0]);
    exit(EXIT_FAILURE);
  }
  
  Graph G;
  int n, m;
  std::unordered_set<int> srcs;
  std::vector<double> pi;
  
  {
    Graph H;
    if(!H.readfromFile(argv[1])){
      fprintf(stderr, "ERROR: reading graph file %s failed.\n", argv[1]);
      print_usage(argv[0]);
      exit(EXIT_FAILURE);
    }
    
    n = H.numV();
    m = H.numE();
    
    std::vector<double> prob(m);
    pi.resize(m);
    {
      FILE *fp;
      if((fp = fopen(argv[2], "r")) == NULL){
        fprintf(stderr, "ERROR: reading probability file %s failed.\n", argv[2]);
        print_usage(argv[0]);
        exit(EXIT_FAILURE);
      }
      
      for(size_t i=0; i<m; ++i){
        fscanf(fp, "%lf", &prob[i]);
      }
      fclose(fp);
    }
    
    {
      FILE *fp;
      if((fp = fopen(argv[3], "r")) == NULL){
        fprintf(stderr, "ERROR: reading source vertices file %s failed.\n", argv[3]);
        print_usage(argv[0]);
        exit(EXIT_FAILURE);
      }
      
      int src;
      while(fscanf(fp, "%d", &src) != EOF){
        srcs.emplace(src);
      }
      fclose(fp);
    }
    
    if(!G.readfromFile(argv[4])){
      fprintf(stderr, "ERROR: reading order file %s failed.\n", argv[4]);
      print_usage(argv[0]);
      exit(EXIT_FAILURE);
    }
    
    for(size_t i=0; i<m; ++i){
      pi[i] = prob[H.etovar(G.e[i].first, G.e[i].second)];
    }
  }
  
  auto cstart = std::chrono::system_clock::now();
  
  std::vector<double> pirev;
  pirev.reserve(pi.size());
  for(auto itr=pi.rbegin(); itr!=pi.rend(); ++itr) pirev.emplace_back(*itr);
  
  BDD_Init(1ULL << 24, 1ULL << 32);
  
  tdzdd::Graph tG;
  for(const auto& edg : G.e){
    tG.addEdge(std::to_string(edg.first), std::to_string(edg.second));
  }
  
  std::vector<BDD> Bu(n);
  
  int cnt=0;
  for(int u=1; u<=n; ++u){
    if(srcs.count(u)){
      continue;
    }else{
      Bu[cnt] = BDD(0);
      for(const auto& src : srcs){
        tG.clearColors();
        tG.setColor(std::to_string(u), 1);
        tG.setColor(std::to_string(src), 1);
        tG.update();
        tdzdd::FrontierBasedSearch fbs(tG, -1, false, false);
        tdzdd::DdStructure<2> dd(fbs);
        dd.useMultiProcessors(false);
        
        BDD Bus = dd.evaluate(ToBDD());
        Bu[cnt] |= Bus;
      }
      ++cnt;
    }
  }
  
  BDDRel BR;
  std::vector<double> res(n+1);
  
  std::vector<BDD> Bret(cnt+1);
  int pcnt = 0;
  Bret[cnt] = BDD(1);
  for(int i=cnt-1; i>=0; --i){
    Bret[i] = Bret[i+1] & (~Bu[i]);
  }
  res[pcnt+(n-cnt)] += BR.compute(Bret[0], pirev);
  for(uint64_t vl = 1; vl < (1ULL << cnt); ++vl){
    int lsbpos = log2ton(static_cast<uint64_t>(((vl ^ (vl-1)) + 1ULL) >> 1));
    pcnt += 1 - lsbpos;
    Bret[lsbpos] = Bret[lsbpos+1] & Bu[lsbpos];
    for(int i=lsbpos-1; i>=0; --i){
      Bret[i] = Bret[i+1] & (~Bu[i]);
    }
    res[pcnt+(n-cnt)] += BR.compute(Bret[0], pirev);
  }
  
  for(size_t ind=0; ind<=n; ++ind){
    printf("%zu: %.15lf\n", ind, res[ind]);
  }
  
  auto cend = std::chrono::system_clock::now();
  double ctime = std::chrono::duration_cast<std::chrono::milliseconds>(cend-cstart).count();
  
  fprintf(stderr, "calc time: %.6lf ms\n", ctime);
  
  return 0;
}