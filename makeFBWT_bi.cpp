#include <memory>

#include <iostream>

using namespace std;

#define CHAR_N 5


// SAはすでにsparseされたもの
// nはひとつの断片の大きさ
void makeFBWT(unique_ptr<int[]>& SA,unique_ptr<unsigned char[]>& fbwt,unsigned char* ref,int fbwtnum,int n){

  unique_ptr<int[]> sa_workspace1(new int[n]{});
  unique_ptr<int[]> sa_workspace2(new int[n]{});

  for(int i = 0;i < n;i++){
    sa_workspace1[i] = SA[i];
  }


  cout << sa_workspace1[0] - 1 << ", " << ref[sa_workspace1[0] - 1] << "," << (int)ref[sa_workspace1[0] - 1] << " a " << n*fbwtnum<< "\n";

  for(int f = 0;f < fbwtnum;f += 2){
    //    cout << "hogehoge : " << f << "\n";
    int count[CHAR_N] = {0,0,0,0,0};
    for(int i = 0;i < n;i++){
      fbwt[f*n + i] = ref[sa_workspace1[i] == 0 ? n*fbwtnum - 1 : sa_workspace1[i] - 1];
      count[(fbwt[f*n + i] + 1)%CHAR_N]++;
      //      count[fbwt[f*n + i]]++;
    }
    
    for(int i = CHAR_N - 1;i > 0;i--){count[i] = count[i - 1];} count[0] = 0;

    for(int i = 1;i < CHAR_N;i++){count[i] += count[i - 1];}

    // for(int i = 0;i < CHAR_N;i++) cout << count[i] << ","; cout << "\n";

    for(int i = 0;i < n;i++){
      sa_workspace2[count[(fbwt[f*n + i] + 1)%CHAR_N]]
        //      sa_workspace2[count[fbwt[f*n + i]]]
        = (sa_workspace1[i] == 0 ? n*fbwtnum - 1 : sa_workspace1[i] - 1);
      count[(fbwt[f*n + i] + 1)%CHAR_N]++;
      //      count[fbwt[f*n + i]]++;
    }



    // for(int i = 1;i < n;i++){
    //   int k = 0;
    //   while(sa_workspace2[i] + k < n && sa_workspace2[i - 1] + k < n){
    //     if(ref[sa_workspace2[i] + k] == ref[sa_workspace2[i - 1] + k]){
    //       k++;
    //     }else if(ref[sa_workspace2[i] + k] < ref[sa_workspace2[i - 1] + k]){
    //       printf("err\n");
    //       exit(1);
    //     }else{
    //       break;
    //     }
    //   }
    // }

    // for(int i = 0;i < n;i++){
    //   printf("%d : ",sa_workspace2[i]);
    //   for(int j = sa_workspace2[i],c = 0; c < 10;c++,j++){
    //     printf("%d",ref[j]);
    //   }
    //   printf("\n");
    // }


    if(f + 1 >= fbwtnum) break;

    //    cout << "hogehoge : " << f + 1<< "\n";

    for(int i = 0;i < CHAR_N;i++){count[i] = 0;}
    for(int i = 0;i < n;i++){
      fbwt[(f + 1)*n + i] = ref[sa_workspace2[i] == 0 ? n*fbwtnum - 1 : sa_workspace2[i] - 1];
      count[(fbwt[(f + 1)*n + i] + 1)%CHAR_N]++;
    }

    for(int i = CHAR_N - 1;i > 0;i--){count[i] = count[i - 1];} count[0] = 0;
    for(int i = 1;i < CHAR_N;i++){count[i] += count[i - 1];}

    for(int i = 0;i < n;i++){
      sa_workspace1[count[(fbwt[(f + 1)*n + i] + 1)%CHAR_N]] 
        = (sa_workspace2[i] == 0 ? n*fbwtnum - 1 : sa_workspace2[i] - 1);
      count[(fbwt[(f + 1)*n + i] + 1)%CHAR_N]++;
      //      count[fbwt[(f + 1)*n + i]]++;
    }

  }

}
