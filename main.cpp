#include "mylib/common.hpp"
#include "mylib/graph.hpp"

#include <cstdio>
#include <cstdlib>
#include <vector>
#include <array>
#include <cstdint>
#include <cassert>
#include <cstring>
#include <chrono>
#include <unordered_map>
#include <unordered_set>

void print_usage(char *fil){
  fprintf(stderr, "Usage: %s [graph_file] [probability_file] [source_file] [order_file]\n", fil);
}

class State{
public:
  std::array<int8_t, 16> comp;
  std::array<uint8_t, 16> numv;
  int8_t cnum;
  
  State() {};
          
  bool operator==(const State& rhs) const {
    return !memcmp(comp.data(), rhs.comp.data(), 16) && !memcmp(numv.data(), rhs.numv.data(), 16);
  }
};

namespace std{
  template <>
  struct hash<State>{
    public:
    uint64_t operator()(const State& s) const{
      const uint64_t *pc = reinterpret_cast<const uint64_t*>(s.comp.data());
      const uint64_t *pn = reinterpret_cast<const uint64_t*>(s.numv.data());
      return *pc * 9088778391368939145ULL + *(pc+1) * 98696038758792049ULL + *pn * 314159257ULL + *(pn+1);
    }
  };
}

class DPBlock{
public:
  double p = 0.0;
  int64_t lo;
  int64_t hi;
  
  DPBlock() {};
};

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
  
  G.buildFrontiers();
  
  size_t src_final = 0;
  for(size_t i=0; i<m; ++i){
    for(const auto& pos : G.fent[i]){
      if(srcs.count(G.mfros[i][pos])) src_final = i;
    }
  }
  
  std::vector<double> res(n+1);
  std::vector<std::vector<DPBlock>> dp(m+1);
  std::vector<std::unordered_map<State, int64_t>> maps(m+1);
  std::vector<int64_t> ssizes(m+1);
  
  // root node
  {
    State root;
    root.cnum = 1;
    root.comp.fill(-1);
    root.numv.fill(0);
    maps[0].emplace(root, 0);
    dp[0].emplace_back();
    dp[0][0].p = 1.0;
  }
  ssizes[0] = 1;
  
  for(size_t i=0; i<m; ++i){
    const auto& now_fro = G.fros[i];
    const auto& med_fro = G.mfros[i];
    const auto& next_fro = G.fros[i+1];
    const auto& now_vpos = G.vpos[i];
    const auto& now_ent = G.fent[i];
    const auto& now_lve = G.flve[i];
    size_t kk = now_fro.size();
    size_t tt = med_fro.size();
    size_t ll = next_fro.size();
    maps[i+1].reserve(ssizes[i] * 2);
    ssizes[i+1] = 0;
    
    // build ast_ent and ast_lve
    std::vector<size_t> ast_ent, not_ent, ast_lve, not_lve;
    for(const auto& pos : now_ent){
      if(srcs.count(med_fro[pos])) ast_ent.emplace_back(pos);
      else                         not_ent.emplace_back(pos);
    }
    for(const auto& pos : now_lve){
      if(srcs.count(med_fro[pos])) ast_lve.emplace_back(pos);
      else                         not_lve.emplace_back(pos);
    }
    
    for(const auto& ent : maps[i]){
      int64_t now_id = ent.second;
      const State& now_state = ent.first;
      const auto& now_comp = now_state.comp;
      const auto& now_numv = now_state.numv;
      int8_t cc = now_state.cnum;
      
      // generate intermediate state
      State med_state;
      auto& med_comp = med_state.comp;
      auto& med_numv = med_state.numv;
      med_comp.fill(-1);
      med_numv.fill(0);
      memcpy(med_comp.data(), now_comp.data(), tt);
      memcpy(med_numv.data(), now_numv.data(), cc);
      for(const auto& pos : ast_ent){
        med_comp[pos] = 0;
      }
      for(const auto& pos : not_ent){
        med_comp[pos] = cc++;
      }
      
      // lo_state processing
      {
        // generate lo_state
        State lo_state;
        auto& lo_comp = lo_state.comp;
        auto& lo_numv = lo_state.numv;
        lo_comp.fill(-1);
        lo_numv.fill(0);
        std::array<uint8_t, 16> tmp_numv;
        tmp_numv.fill(0);
        memcpy(tmp_numv.data(), med_numv.data(), cc);
        memcpy(lo_comp.data(), med_comp.data(), ll);
        for(const auto& pos : now_lve){
          if(pos < ll) lo_comp[pos] = -1;
          ++tmp_numv[med_comp[pos]];
        }
        std::array<int8_t, 16> renum;
        renum.fill(-1);
        bool prune = (i >= src_final);
        renum[0] = 0;
        int8_t cc_new = 1;
        for(auto&& val : lo_comp){
          if(val < 0) continue;
          if(renum[val] < 0) renum[val] = cc_new++;
          val = renum[val];
          prune = prune && val;
        }
        for(size_t c=0; c<cc; ++c){
          if(renum[c] >= 0){
            lo_numv[renum[c]] = tmp_numv[c];
          }
        }
        lo_state.cnum = cc_new;
        
        // find or generate id
        int64_t lo_id;
        if(prune){
          lo_id = -(static_cast<int64_t>(lo_numv[0]+1));
        }else{
          auto it = maps[i+1].find(lo_state);
          if(it != maps[i+1].end()){ // maps[i+1] has already had entry
            lo_id = it->second;
          }else{                     // there is no entry
            maps[i+1].emplace(lo_state, ssizes[i+1]);
            lo_id = ssizes[i+1]++;
          }
        }
        dp[i][now_id].lo = lo_id;
      }
      
      // hi_state processing
      int8_t cat_to   = med_comp[now_vpos.first];
      int8_t cat_from = med_comp[now_vpos.second];
      {
        // generate hi_state
        State hi_state;
        auto& hi_comp = hi_state.comp;
        auto& hi_numv = hi_state.numv;
        hi_comp.fill(-1);
        hi_numv.fill(0);
        std::array<uint8_t, 16> tmp_numv;
        tmp_numv.fill(0);
        memcpy(tmp_numv.data(), med_numv.data(), cc);
        memcpy(hi_comp.data(), med_comp.data(), ll);
        for(const auto& pos : now_lve){
          if(pos < ll) hi_comp[pos] = -1;
          ++tmp_numv[med_comp[pos]];
        }
        std::array<int8_t, 16> renum;
        renum.fill(-1);
        bool prune = (i >= src_final);
        renum[0] = 0;
        if(cat_to == 0)   renum[cat_from] = 0;
        if(cat_from == 0) renum[cat_to]   = 0;
        int8_t cc_new = 1;
        for(auto&& val : hi_comp){
          if(val < 0) continue;
          if(renum[val] < 0){
            renum[val] = cc_new++;
            if(val == cat_to)        renum[cat_from] = renum[val];
            else if(val == cat_from) renum[cat_to]   = renum[val];
          }
          val = renum[val];
          prune = prune && val;
        }
        for(size_t c=0; c<cc; ++c){
          if(renum[c] >= 0){
            hi_numv[renum[c]] += tmp_numv[c];
          }
        }
        hi_state.cnum = cc_new;
        
        // find or generate id
        int64_t hi_id;
        if(prune){
          hi_id = -(static_cast<int64_t>(hi_numv[0]+1));
        }else{
          auto it = maps[i+1].find(hi_state);
          if(it != maps[i+1].end()){ // maps[i+1] has already had entry
            hi_id = it->second;
          }else{                     // there is no entry
            maps[i+1].emplace(hi_state, ssizes[i+1]);
            hi_id = ssizes[i+1]++;
          }
        }
        dp[i][now_id].hi = hi_id;
      }
    }
    maps[i].clear();
    dp[i+1].resize(ssizes[i+1]);
  }
  
  // dp calculation
  for(size_t i=0; i<m; ++i){
    for(const auto& ent : dp[i]){
      if(ent.lo >= 0) dp[i+1][ent.lo].p += (1.0 - pi[i]) * ent.p;
      else            res[-(ent.lo+1)]  += (1.0 - pi[i]) * ent.p;
      if(ent.hi >= 0) dp[i+1][ent.hi].p += pi[i]         * ent.p;
      else            res[-(ent.hi+1)]  += pi[i]         * ent.p;
    }
  }

  for(size_t ind=0; ind<=n; ++ind){
    printf("%zu: %.15lf\n", ind, res[ind]);
  }
  
  auto cend = std::chrono::system_clock::now();
  double ctime = std::chrono::duration_cast<std::chrono::milliseconds>(cend-cstart).count();
  
  fprintf(stderr, "calc time: %.6lf ms\n", ctime);
  
  int64_t snum;
  for(int i=0; i<m; ++i){
    snum += ssizes[i];
  }
  fprintf(stderr, "#(states): %lld\n", snum);
  
  return 0;
}