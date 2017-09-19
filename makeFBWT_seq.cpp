#include <memory>

#include <iostream>

using namespace std;

#define CHAR_N 5


// SAはすでにsparseされたもの
// nはひとつの断片の大きさ
void makeFBWT(unique_ptr<uint32_t[]>& SA,unique_ptr<unsigned char[]>& fbwt,unsigned char* ref,uint32_t fbwtnum,uint32_t n){

  unique_ptr<uint32_t[]> sa_workspace1(new uint32_t[n]{});
  unique_ptr<uint32_t[]> sa_workspace2(new uint32_t[n]{});

  for(uint32_t i = 0;i < n;i++){
    sa_workspace1[i] = SA[i];
  }


  for(uint32_t f = 0;f < fbwtnum;f += 2){
    uint32_t count[CHAR_N] = {0,0,0,0,0};
    for(uint64_t i = 0;i < n;i++){
      fbwt[f*n + i] = ref[sa_workspace1[i] == 0 ? n*fbwtnum - 1 : sa_workspace1[i] - 1];
      count[fbwt[f*n + i]]++;
    }
    
    for(int i = CHAR_N - 1;i > 0;i--){count[i] = count[i - 1];} count[0] = 0;

    for(int i = 1;i < CHAR_N;i++){count[i] += count[i - 1];}

    for(uint64_t i = 0;i < n;i++){
      sa_workspace2[count[fbwt[f*n + i]]]
        = (sa_workspace1[i] == 0 ? n*fbwtnum - 1 : sa_workspace1[i] - 1);
      count[fbwt[f*n + i]]++;
    }

    if(f + 1 >= fbwtnum) break;


    for(int i = 0;i < CHAR_N;i++){count[i] = 0;}
    for(uint64_t i = 0;i < n;i++){
      fbwt[(f + 1)*n + i] = ref[sa_workspace2[i] == 0 ? n*fbwtnum - 1 : sa_workspace2[i] - 1];
      count[fbwt[(f + 1)*n + i]]++;
    }

    for(int i = CHAR_N - 1;i > 0;i--){count[i] = count[i - 1];} count[0] = 0;
    for(int i = 1;i < CHAR_N;i++){count[i] += count[i - 1];}

    for(uint64_t i = 0;i < n;i++){
      sa_workspace1[count[fbwt[(f + 1)*n + i]]] 
        = (sa_workspace2[i] == 0 ? n*fbwtnum - 1 : sa_workspace2[i] - 1);
      count[fbwt[(f + 1)*n + i]]++;
    }

  }

}
