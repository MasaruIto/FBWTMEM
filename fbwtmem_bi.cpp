#include <iostream>
#include <fstream>
#include <memory>
#include <utility>
#include <vector>
#include <algorithm>
#include <stdint.h>
//#include <random>
#include <stdlib.h>

#include "sais/sais.h"
#include <string.h>
#include "makeFBWT_bi.h"
#include <assert.h>
#include <sys/mman.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/time.h>

#include <thread>

#include <stdio.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <stdlib.h>

long show_getrusage(){
  struct rusage r;
  if (getrusage(RUSAGE_SELF, &r) != 0) {
    fprintf(stderr,"err\n");
  }
  return r.ru_maxrss;
}


using namespace std;

#define CHAR_N 5

int countMEMS = 0;
int firstCandidateNum = 0;
int firstCountDistribution[10240];


#define INTERVAL_DIST_MAX_LEN 102400
uint64_t initMEMIntervalfLenDist[INTERVAL_DIST_MAX_LEN*32];
uint64_t updateMEMIntervalLenDist[INTERVAL_DIST_MAX_LEN*32];

long maxMemorySize = 0;


void convertATGC(unsigned char* ref,int n){
  for(int i = 0;i < n;i++){
    switch(ref[i]){
    case 'A':
      ref[i] = 0;
      break;
    case 'C':
      ref[i] = 1;
      break;
    case 'G':
      ref[i] = 2;
      break;
    case 'T':
      ref[i] = 3;
      break;
    default:
      ref[i] = 4;
    }
  }
}


class C_OCC{

private:
  unique_ptr< unsigned char[]> FBWT;
  unique_ptr< unsigned char[]> RFBWT;

  class C{
  private:
    unique_ptr<int[]> ctable;
    int fnum;

    //    const string fname = "c";

  public:
    C(const char* dirname,string fname){
      FILE* fp;
      string filename = string(dirname) + "/" + fname;


      if((fp = fopen(filename.c_str(), "rb")) == NULL ) {
        printf("cannot open output file %s\n",filename.c_str());
        exit(EXIT_FAILURE);
      }

      int n;
      if(fseek(fp, 0, SEEK_END) == 0) {
        n = ftell(fp)/sizeof(int);
        rewind(fp);
        if(n < 0) {
          fprintf(stderr, "%s: Cannot ftell `%s': ", "C", filename.c_str());
          perror(NULL);
          exit(EXIT_FAILURE);
        }
      } else {
        fprintf(stderr, "%s: Cannot fseek `%s': ", "C", filename.c_str());
        perror(NULL);
        exit(EXIT_FAILURE);
      }

      ctable = make_unique<int[]>(n);

      if(fread(ctable.get(), sizeof(int), (size_t)n, fp) != (size_t)n) {
        fprintf(stderr, "%s: %s `%s': ",
                "C",
                (ferror(fp) || !feof(fp)) ? "Cannot read from" : "Unexpected EOF in",
                filename.c_str());
        perror(NULL);
        exit(EXIT_FAILURE);
      }

      fnum = n/(CHAR_N + 1);


      fclose(fp);
    }


    C(unique_ptr< unsigned char[]>& FBWT,int n,int fn){
      ctable = make_unique<int[]>((CHAR_N + 1)*fn);
      this->fnum = fn;

      for(int i = 0;i < (CHAR_N + 1)*fn;i++){
        ctable[i] = 0;
      }

      for(int f = 0;f < fn;f++){
        for(int i = 0;i < n;i++){
          ctable[(FBWT[i + f*n] + 1)%CHAR_N + (CHAR_N + 1)*f]++;
        }
      }

      // for(int i = 0;i < CHAR_N;i++){
      //   cout << ctable[i] << ",";
      // }
      // cout << "\n";



      for(int f = 0;f < fn;f++){
        for(int i = 1;i < CHAR_N;i++){
          ctable[i + f*(CHAR_N + 1)] += ctable[i - 1 + f*(CHAR_N + 1)];
        }
        for(int i = CHAR_N;i > 0;i--){
          ctable[i + f*(CHAR_N + 1)] = ctable[i - 1 + f*(CHAR_N + 1)];
        }
        for(int i = 0;i < CHAR_N;i++){
          ctable[f*(CHAR_N + 1) + i] = ctable[f*(CHAR_N + 1) + i + 1];
        }
        ctable[f*(CHAR_N + 1) + 4 + 1] = 0;
      }
    }


    int getFnum(){
      return fnum;
    }

    void printC(){
      for(int f = 0;f < fnum;f++){
        for(int i = 0;i <= CHAR_N;i++){
          cout << ctable[i + f*(CHAR_N + 1)] << " ";
        }
        cout << "\n";
      }
    }

    int getC(int c,int fn){
      return ctable[fn*(CHAR_N + 1) + c];
    }

    void outputFile(const char* dirname,string fname){
      FILE* fp;
      string filename = string(dirname) + "/" + fname;
      if((fp = fopen(filename.c_str(), "wb")) == NULL ) {
        printf("cannot open output file %s\n",filename.c_str());
        exit(EXIT_FAILURE);
      }
      
      int n = fwrite(ctable.get() , sizeof(int), (CHAR_N + 1)*fnum, fp);
      cout << n << " output ctable\n";
      fclose(fp);
    }
  };


  class OCC{
    //  private:
  public:
    // i番目以前のものの個数
    //    unique_ptr<int[]> occtable;
    C_OCC& c_occ;
    int fnum;
    int n;

    int compressEntrySize;


    static const int oneEntryNum = 128;

    //512bitにつめられるだけつめる
    //96ko ATGC32bit*4 96ko //24byte 4byte*4 24byte
    //1エントリー192個// 2bitを一文字としたときであるが、3bit一文字にすることにした
    //1エントリー128個
    struct oneEntry{
      uint64_t pre[3];
      uint32_t atgc[4];
      uint64_t post[3];
    };

    unique_ptr<oneEntry[]> compressedOcctable;

    void constructOccTableStructure(unique_ptr<unsigned char[]>& FBWT,int fn,int n,int compressEntrySize){
      compressedOcctable = make_unique<oneEntry[]>(compressEntrySize*fnum);

      // cout << compressEntrySize << " #hoge " << fn << " \n";

      for(int f = 0;f < fn;f++){
        int count[5] = {0,0,0,0,0};

        for(int i = 0;i < compressEntrySize;i++){
          uint64_t tmp = 0;
          for(int j = 0;j < 21;j++){
            tmp |= ((uint64_t)(FBWT[i*oneEntryNum + j + f*n] & 0x7) << (j*3));
            count[FBWT[i*oneEntryNum + j + f*n]]++;
          }
          tmp |= ((uint64_t)(FBWT[i*oneEntryNum + 21 + f*n] & 0x7) << (21*3));
          compressedOcctable[i + compressEntrySize*f].pre[0] = tmp;
          count[FBWT[i*oneEntryNum + 21 + f*n] & 0x7]++;

          tmp = 0;
          tmp |= ((uint64_t)(FBWT[i*oneEntryNum + 21 + f*n] & 0x7) >> 1);
          for(int j = 1;j < 21;j++){
            tmp |= ((uint64_t)(FBWT[i*oneEntryNum + j + 1*21 + f*n] & 0x7) << ((j - 1)*3 + 2));
            count[FBWT[i*oneEntryNum + j + 21 + f*n]]++;
          }
          tmp |= ((uint64_t)(FBWT[i*oneEntryNum + 42 + f*n] & 0x7) << (21*3 - 1));
          compressedOcctable[i + compressEntrySize*f].pre[1] = tmp;
          count[FBWT[i*oneEntryNum + 42 + f*n] & 0x7]++;

          tmp = 0;
          tmp |= ((uint64_t)(FBWT[i*oneEntryNum + 42 + f*n] & 0x7) >> 2);
          for(int j = 1;j < 22;j++){
            tmp |= ((uint64_t)(FBWT[i*oneEntryNum + j + 42 + f*n] & 0x7) << ((j - 1)*3 + 1));
            count[FBWT[i*oneEntryNum + j + 42 + f*n]]++;
          }
          compressedOcctable[i + compressEntrySize*f].pre[2] = tmp;

	  if(f == 1 && ((i == 6174332/oneEntryNum) ||  (i == 6174332/oneEntryNum + 1))){
	    cout << count[0] << ", " << count[1] << ", " << count[2] << ", " << count[3] << "\n";
	    int _count[5] = {0,0,0,0,0};
	    for(int j = 0;j < oneEntryNum;j++){
	      cout << (int)FBWT[i*oneEntryNum + oneEntryNum/2 + j + f*n];
	      _count[FBWT[i*oneEntryNum + oneEntryNum/2 + j + f*n]]++;
	    }
	    cout << "\n";
	    cout << _count[0] << ", " << _count[1] << ", " << _count[2] << ", " << _count[3] << "\n";
	  }
	  
          compressedOcctable[i + compressEntrySize*f].atgc[0] = count[0];
          compressedOcctable[i + compressEntrySize*f].atgc[1] = count[1];
          compressedOcctable[i + compressEntrySize*f].atgc[2] = count[2];
          compressedOcctable[i + compressEntrySize*f].atgc[3] = count[3];

          tmp = 0;
          for(int j = 0;j < 21;j++){
            tmp |= ((uint64_t)(FBWT[i*oneEntryNum + oneEntryNum/2 + j + f*n] & 0x7) << (j*3));
            count[FBWT[i*oneEntryNum + oneEntryNum/2 + j + f*n]]++;
          }
          tmp |= ((uint64_t)(FBWT[i*oneEntryNum + oneEntryNum/2 + 21 + f*n] & 0x7) << (21*3));
          compressedOcctable[i + compressEntrySize*f].post[0] = tmp;
          count[FBWT[i*oneEntryNum + oneEntryNum/2 + 21 + f*n] & 0x7]++;

          tmp = 0;
          tmp |= ((uint64_t)(FBWT[i*oneEntryNum + oneEntryNum/2 + 21 + f*n] & 0x7) >> 1);
          for(int j = 1;j < 21;j++){
            tmp |= ((uint64_t)(FBWT[i*oneEntryNum + oneEntryNum/2 + j + 1*21 + f*n] & 0x7) << ((j - 1)*3 + 2));
            count[FBWT[i*oneEntryNum + oneEntryNum/2 + j + 1*21 + f*n]]++;
          }
          tmp |= ((uint64_t)(FBWT[i*oneEntryNum + oneEntryNum/2 + 42 + f*n] & 0x7) << (21*3 - 1));
          compressedOcctable[i + compressEntrySize*f].post[1] = tmp;
          count[FBWT[i*oneEntryNum + oneEntryNum/2 + 42 + f*n] & 0x7]++;

          tmp = 0;
          tmp |= ((uint64_t)(FBWT[i*oneEntryNum + oneEntryNum/2 + 42 + f*n] & 0x7) >> 2);
          for(int j = 1;j < 22;j++){
            tmp |= ((uint64_t)(FBWT[i*oneEntryNum + oneEntryNum/2 + j + 2*21 + f*n] & 0x7) << ((j - 1)*3 + 1));
            count[FBWT[i*oneEntryNum + oneEntryNum/2 + j + 2*21 + f*n]]++;
          }
          compressedOcctable[i + compressEntrySize*f].post[2] = tmp;

          // for(int k = 0;k < 3;k++){
          //   uint64_t tmp = 0;
          //   for(int j = 0;j < 32;j++){
          //     tmp |= ((uint64_t)(FBWT[i*oneEntryNum + j + k*32 + oneEntryNum/2 + f*n] & 0x7) << (j*3));
          //     count[FBWT[i*oneEntryNum + j + k*32 + oneEntryNum/2 + f*n]]++;
          //   }
          //   compressedOcctable[i + compressEntrySize*f].post[k] = tmp;
          // }
        }
      }
    }

    char vv[5] = {'N','A','C','G','T'};

    void print64bit(uint64_t h){
      for(int k = 0;k < 32;k++)
        cout << vv[(h >> (k*2)) & 0x3];
    }


    void printOneCompressedEntry(int i){

      for(int k = 0;k < 21;k++)
        cout << vv[(compressedOcctable[i].pre[0] >> (k*3)) & 0x7];
      int tmp = ((compressedOcctable[i].pre[0] >> (63)) & 0x1) | (((compressedOcctable[i].pre[1]) & 0x3) << 1);
      cout << vv[tmp];
      cout << " ";

      for(int k = 0;k < 20;k++)
        cout << vv[(compressedOcctable[i].pre[1] >> (k*3 + 2)) & 0x7];

      tmp = ((compressedOcctable[i].pre[1] >> (62)) & 0x3) | (((compressedOcctable[i].pre[2]) & 0x1) << 2);
      cout << vv[tmp];
      cout << " ";

      for(int k = 0;k < 21;k++)
        cout << vv[(compressedOcctable[i].pre[2] >> (k*3 + 1)) & 0x7];
      cout << " ";

      // for(int j = 0;j < 3;j++){
      //   //        printf(" %016lx ", compressedOcctable[i].pre[j]);
      //   cout << " ";
      //   for(int k = 0;k < 32;k++)
      //     cout << vv[(compressedOcctable[i].pre[j] >> (k*2)) & 0x3];
      // }

      for(int j = 0;j < 4;j++){
        cout << " " << compressedOcctable[i].atgc[j];
      }      

      cout << " ";

      for(int k = 0;k < 21;k++)
        cout << vv[(compressedOcctable[i].post[0] >> (k*3)) & 0x7];
      tmp = ((compressedOcctable[i].post[0] >> (63)) & 0x1) | (((compressedOcctable[i].post[1]) & 0x3) << 1);
      cout << vv[tmp];
      cout << " ";

      for(int k = 0;k < 20;k++)
        cout << vv[(compressedOcctable[i].post[1] >> (k*3 + 2)) & 0x7];

      tmp = ((compressedOcctable[i].post[1] >> (62)) & 0x3) | (((compressedOcctable[i].post[2]) & 0x1) << 2);
      cout << vv[tmp];
      cout << " ";

      for(int k = 0;k < 21;k++)
        cout << vv[(compressedOcctable[i].post[2] >> (k*3 + 1)) & 0x7];

      // for(int j = 0;j < 3;j++){
      //   //        printf(" %016lx ", compressedOcctable[i].post[j]);
      //   cout << " ";
      //   for(int k = 0;k < 32;k++)
      //     cout << vv[(compressedOcctable[i].post[j] >> (k*2)) & 0x3];
      // }

      cout << "\n";
    }


    uint64_t atgcMask[4] = {0x0000000000000000,0x9249249249249249,0x2492492492492492,0xb6db6db6db6db6db};
    // uint64_t atgcMask[4] = {0x9249249249249249,0x2492492492492492,0xb6db6db6db6db6db,0x4924924924924924};
    // uint64_t atgcMask[4] = {0x0000000000000000,0x5555555555555555,0xcccccccccccccccc,0xffffffffffffffff};

    uint64_t preMask[22] = {
      0xffffffffffffffff, 0xfffffffffffffff8, 0xffffffffffffffc0, 0xfffffffffffffe00,
      0xfffffffffffff000, 0xffffffffffff8000, 0xfffffffffffc0000, 0xffffffffffe00000,
      0xffffffffff000000, 0xfffffffff8000000, 0xffffffffc0000000, 0xfffffffe00000000,
      0xfffffff000000000, 0xffffff8000000000, 0xfffffc0000000000, 0xffffe00000000000,
      0xffff000000000000, 0xfff8000000000000, 0xffc0000000000000, 0xfe00000000000000,
      0xf000000000000000, 0x8000000000000000
      // 0xffffffffffffffff, 0xfffffffffffffffc, 0xfffffffffffffff0, 0xffffffffffffffc0,
      // 0xffffffffffffff00, 0xfffffffffffffc00, 0xfffffffffffff000, 0xffffffffffffc000,
      // 0xffffffffffff0000, 0xfffffffffffc0000, 0xfffffffffff00000, 0xffffffffffc00000,
      // 0xffffffffff000000, 0xfffffffffc000000, 0xfffffffff0000000, 0xffffffffc0000000,
      // 0xffffffff00000000, 0xfffffffc00000000, 0xfffffff000000000, 0xffffffc000000000,
      // 0xffffff0000000000, 0xfffffc0000000000, 0xfffff00000000000, 0xffffc00000000000,
      // 0xffff000000000000, 0xfffc000000000000, 0xfff0000000000000, 0xffc0000000000000,
      // 0xff00000000000000, 0xfc00000000000000, 0xf000000000000000, 0xc000000000000000
    };

    uint64_t postMask[22] = {
      0x0000000000000000, 0x0000000000000007, 0x000000000000003f, 0x00000000000001ff,
      0x0000000000000fff, 0x0000000000007fff, 0x000000000003ffff, 0x00000000001fffff,
      0x0000000000ffffff, 0x0000000007ffffff, 0x000000003fffffff, 0x00000001ffffffff,
      0x0000000fffffffff, 0x0000007fffffffff, 0x000003ffffffffff, 0x00001fffffffffff,
      0x0000ffffffffffff, 0x0007ffffffffffff, 0x003fffffffffffff, 0x01ffffffffffffff,
      0x0fffffffffffffff, 0x7fffffffffffffff
      // 0x0000000000000000, 0x0000000000000003, 0x000000000000000f, 0x000000000000003f,
      // 0x00000000000000ff, 0x00000000000003ff, 0x0000000000000fff, 0x0000000000003fff,
      // 0x000000000000ffff, 0x000000000003ffff, 0x00000000000fffff, 0x00000000003fffff,
      // 0x0000000000ffffff, 0x0000000003ffffff, 0x000000000fffffff, 0x000000003fffffff,
      // 0x00000000ffffffff, 0x00000003ffffffff, 0x0000000fffffffff, 0x0000003fffffffff,
      // 0x000000ffffffffff, 0x000003ffffffffff, 0x00000fffffffffff, 0x00003fffffffffff,
      // 0x0000ffffffffffff, 0x0003ffffffffffff, 0x000fffffffffffff, 0x003fffffffffffff,
      // 0x00ffffffffffffff, 0x03ffffffffffffff, 0x0fffffffffffffff, 0x3fffffffffffffff
    };





    void printCompressedOcctable(){

      for(int f = 0;f < fnum;f++){
        for(int i = 0;i < compressEntrySize;i++){
          printOneCompressedEntry(i + f*compressEntrySize);
        }
        cout << "-------------------------------------\n";
      }
    }

    void testGetOccFromCompressed(){
      int i;
      char c;
      int fn;
      int occ;

      i = 63; c = 4; fn = 0; occ = getOccFromCompressed(i,c,fn);
      cout << "i : " << i << " c : " << (int)c << " fn : " << fn << " occ : " << occ << "\n";

      i = 62; c = 4; fn = 0; occ = getOccFromCompressed(i,c,fn);
      cout << "i : " << i << " c : " << (int)c << " fn : " << fn << " occ : " << occ << "\n";

      i = 43; c = 1; fn = 0; occ = getOccFromCompressed(i,c,fn);
      cout << "i : " << i << " c : " << (int)c << " fn : " << fn << " occ : " << occ << "\n";

      i = 42; c = 1; fn = 0; occ = getOccFromCompressed(i,c,fn);
      cout << "i : " << i << " c : " << (int)c << " fn : " << fn << " occ : " << occ << "\n";

      i = 41; c = 1; fn = 0; occ = getOccFromCompressed(i,c,fn);
      cout << "i : " << i << " c : " << (int)c << " fn : " << fn << " occ : " << occ << "\n";

      i = 22; c = 3; fn = 0; occ = getOccFromCompressed(i,c,fn);
      cout << "i : " << i << " c : " << (int)c << " fn : " << fn << " occ : " << occ << "\n";

      i = 21; c = 3; fn = 0; occ = getOccFromCompressed(i,c,fn);
      cout << "i : " << i << " c : " << (int)c << " fn : " << fn << " occ : " << occ << "\n";

      i = 20; c = 3; fn = 0; occ = getOccFromCompressed(i,c,fn);
      cout << "i : " << i << " c : " << (int)c << " fn : " << fn << " occ : " << occ << "\n";

      i = 0; c = 3; fn = 0; occ = getOccFromCompressed(i,c,fn);
      cout << "i : " << i << " c : " << (int)c << " fn : " << fn << " occ : " << occ << "\n";

      i = 1; c = 3; fn = 0; occ = getOccFromCompressed(i,c,fn);
      cout << "i : " << i << " c : " << (int)c << " fn : " << fn << " occ : " << occ << "\n";

      i = 64; c = 1; fn = 0; occ = getOccFromCompressed(i,c,fn);
      cout << "i : " << i << " c : " << (int)c << " fn : " << fn << " occ : " << occ << "\n";

      i = 65; c = 1; fn = 0; occ = getOccFromCompressed(i,c,fn);
      cout << "i : " << i << " c : " << (int)c << " fn : " << fn << " occ : " << occ << "\n";

      i = 84; c = 3; fn = 0; occ = getOccFromCompressed(i,c,fn);
      cout << "i : " << i << " c : " << (int)c << " fn : " << fn << " occ : " << occ << "\n";

      i = 85; c = 3; fn = 0; occ = getOccFromCompressed(i,c,fn);
      cout << "i : " << i << " c : " << (int)c << " fn : " << fn << " occ : " << occ << "\n";

      i = 86; c = 3; fn = 0; occ = getOccFromCompressed(i,c,fn);
      cout << "i : " << i << " c : " << (int)c << " fn : " << fn << " occ : " << occ << "\n";

      i = 87; c = 3; fn = 0; occ = getOccFromCompressed(i,c,fn);
      cout << "i : " << i << " c : " << (int)c << " fn : " << fn << " occ : " << occ << "\n";

      i = 105; c = 2; fn = 0; occ = getOccFromCompressed(i,c,fn);
      cout << "i : " << i << " c : " << (int)c << " fn : " << fn << " occ : " << occ << "\n";

      i = 106; c = 2; fn = 0; occ = getOccFromCompressed(i,c,fn);
      cout << "i : " << i << " c : " << (int)c << " fn : " << fn << " occ : " << occ << "\n";

      i = 107; c = 2; fn = 0; occ = getOccFromCompressed(i,c,fn);
      cout << "i : " << i << " c : " << (int)c << " fn : " << fn << " occ : " << occ << "\n";

      i = 105; c = 1; fn = 0; occ = getOccFromCompressed(i,c,fn);
      cout << "i : " << i << " c : " << (int)c << " fn : " << fn << " occ : " << occ << "\n";

      i = 106; c = 1; fn = 0; occ = getOccFromCompressed(i,c,fn);
      cout << "i : " << i << " c : " << (int)c << " fn : " << fn << " occ : " << occ << "\n";

      i = 107; c = 1; fn = 0; occ = getOccFromCompressed(i,c,fn);
      cout << "i : " << i << " c : " << (int)c << " fn : " << fn << " occ : " << occ << "\n";

      i = 128; c = 4; fn = 0; occ = getOccFromCompressed(i,c,fn);
      cout << "i : " << i << " c : " << (int)c << " fn : " << fn << " occ : " << occ << "\n";

      i = 129; c = 4; fn = 0; occ = getOccFromCompressed(i,c,fn);
      cout << "i : " << i << " c : " << (int)c << " fn : " << fn << " occ : " << occ << "\n";

    }



    void printBit(uint64_t a){
      for(int i = 0;i < 64;i++){
        cout << ((a >> (63 - i)) & 0x1);
        if(i % 3 == 0) cout << " ";
      }
      cout << '\n';
    }



    int getOccFromCompressed(int i,unsigned char c,int fn){

      uint64_t _atgcMask = atgcMask[c];

      if(i == n){
        int p = (n - 1) / oneEntryNum;
        int ret = compressedOcctable[p + fn*compressEntrySize].atgc[c];
        uint64_t tmp = (~((compressedOcctable[p + fn*compressEntrySize].post[2] >> 0x1) ^ _atgcMask));

        tmp = tmp & (tmp >> 1) & (tmp >> 2) & 0x1249249249249249;
        tmp = (tmp + (tmp >> 3)) & 0x71c71c71c71c71c7;
        tmp = (tmp + (tmp >> 6)) & 0x703f03f03f03f03f;
        tmp = (tmp + (tmp >> 12)) & 0x0fff000fff000fff;
        tmp = (tmp + (tmp >> 24)) & 0x0fff000000ffffff;
        tmp = (tmp + (tmp >> 48)) & 0x0000ffffffffffff;

        ret += tmp;

        uint64_t tmp2 = (~(((compressedOcctable[p + fn*compressEntrySize].post[1] >> 2) 
                            | ((compressedOcctable[p + fn*compressEntrySize].post[2] & 0x1) << 62)) ^ _atgcMask));

        tmp2 = tmp2 & (tmp2 >> 1) & (tmp2 >> 2) & 0x1249249249249249;
        tmp2 = (tmp2 + (tmp2 >> 3)) & 0x71c71c71c71c71c7;
        tmp2 = (tmp2 + (tmp2 >> 6)) & 0x703f03f03f03f03f;
        tmp2 = (tmp2 + (tmp2 >> 12)) & 0x0fff000fff000fff;
        tmp2 = (tmp2 + (tmp2 >> 24)) & 0x0fff000000ffffff;
        tmp2 = (tmp2 + (tmp2 >> 48)) & 0x0000ffffffffffff;

        ret += tmp2;

        uint64_t tmp3 = (~((compressedOcctable[p + fn*compressEntrySize].post[0]) ^ _atgcMask));
        tmp3 = tmp3 & (tmp3 >> 1) & (tmp3 >> 2) & 0x1249249249249249;
        tmp3 = (tmp3 + (tmp3 >> 3)) & 0x71c71c71c71c71c7;
        tmp3 = (tmp3 + (tmp3 >> 6)) & 0x703f03f03f03f03f;
        tmp3 = (tmp3 + (tmp3 >> 12)) & 0x0fff000fff000fff;
        tmp3 = (tmp3 + (tmp3 >> 24)) & 0x0fff000000ffffff;
        tmp3 = (tmp3 + (tmp3 >> 48)) & 0x0000ffffffffffff;

        ret += (((compressedOcctable[p + fn*compressEntrySize].post[0] >> 63) 
                | ((compressedOcctable[p + fn*compressEntrySize].post[1] & 0x3) << 1)) == c);

        ret += tmp3;
        return ret;
      }

      int p = i / oneEntryNum;
      int ret = compressedOcctable[p + fn*compressEntrySize].atgc[c];
      int num = (i % oneEntryNum);


      if(num < 21){
        uint64_t tmp = (~(compressedOcctable[p + fn*compressEntrySize].pre[0] ^ _atgcMask)) & preMask[num];
        uint64_t tmp2 = (~(((compressedOcctable[p + fn*compressEntrySize].pre[1] << 1) 
                            | (compressedOcctable[p + fn*compressEntrySize].pre[0] >> 63)) ^ _atgcMask));
        uint64_t tmp3 = (~((compressedOcctable[p + fn*compressEntrySize].pre[2] >> 1) ^ _atgcMask));


        tmp = tmp & (tmp >> 1) & (tmp >> 2) & 0x1249249249249249;
        tmp = (tmp + (tmp >> 3)) & 0x71c71c71c71c71c7;
        tmp = (tmp + (tmp >> 6)) & 0x703f03f03f03f03f;
        tmp = (tmp + (tmp >> 12)) & 0x0fff000fff000fff;
        tmp = (tmp + (tmp >> 24)) & 0x0fff000000ffffff;
        tmp = (tmp + (tmp >> 48)) & 0x0000ffffffffffff;

        ret -= tmp;

        tmp2 = tmp2 & (tmp2 >> 1) & (tmp2 >> 2) & 0x1249249249249249;
        tmp2 = (tmp2 + (tmp2 >> 3)) & 0x71c71c71c71c71c7;
        tmp2 = (tmp2 + (tmp2 >> 6)) & 0x703f03f03f03f03f;
        tmp2 = (tmp2 + (tmp2 >> 12)) & 0x0fff000fff000fff;
        tmp2 = (tmp2 + (tmp2 >> 24)) & 0x0fff000000ffffff;
        tmp2 = (tmp2 + (tmp2 >> 48)) & 0x0000ffffffffffff;

        ret -= tmp2;

        tmp3 = tmp3 & (tmp3 >> 1) & (tmp3 >> 2) & 0x1249249249249249;
        tmp3 = (tmp3 + (tmp3 >> 3)) & 0x71c71c71c71c71c7;
        tmp3 = (tmp3 + (tmp3 >> 6)) & 0x703f03f03f03f03f;
        tmp3 = (tmp3 + (tmp3 >> 12)) & 0x0fff000fff000fff;
        tmp3 = (tmp3 + (tmp3 >> 24)) & 0x0fff000000ffffff;
        tmp3 = (tmp3 + (tmp3 >> 48)) & 0x0000ffffffffffff;

        ret -= tmp3;

        ret -= (((compressedOcctable[p + fn*compressEntrySize].pre[1] >> 62) | ((compressedOcctable[p + fn*compressEntrySize].pre[2] & 0x1) << 0x2)) == c);
      }else if(num == 21){

        uint64_t tmp2 = (~(((compressedOcctable[p + fn*compressEntrySize].pre[1] << 1) 
                            | (compressedOcctable[p + fn*compressEntrySize].pre[0] >> 63)) ^ _atgcMask));

        tmp2 = tmp2 & (tmp2 >> 1) & (tmp2 >> 2) & 0x1249249249249249;
        tmp2 = (tmp2 + (tmp2 >> 3)) & 0x71c71c71c71c71c7;
        tmp2 = (tmp2 + (tmp2 >> 6)) & 0x703f03f03f03f03f;
        tmp2 = (tmp2 + (tmp2 >> 12)) & 0x0fff000fff000fff;
        tmp2 = (tmp2 + (tmp2 >> 24)) & 0x0fff000000ffffff;
        tmp2 = (tmp2 + (tmp2 >> 48)) & 0x0000ffffffffffff;

        ret -= tmp2;

        uint64_t tmp3 = (~((compressedOcctable[p + fn*compressEntrySize].pre[2] >> 1) ^ _atgcMask));
        tmp3 = tmp3 & (tmp3 >> 1) & (tmp3 >> 2) & 0x1249249249249249;
        tmp3 = (tmp3 + (tmp3 >> 3)) & 0x71c71c71c71c71c7;
        tmp3 = (tmp3 + (tmp3 >> 6)) & 0x703f03f03f03f03f;
        tmp3 = (tmp3 + (tmp3 >> 12)) & 0x0fff000fff000fff;
        tmp3 = (tmp3 + (tmp3 >> 24)) & 0x0fff000000ffffff;
        tmp3 = (tmp3 + (tmp3 >> 48)) & 0x0000ffffffffffff;

        ret -= tmp3;

        ret -= (((compressedOcctable[p + fn*compressEntrySize].pre[1] >> 62) | ((compressedOcctable[p + fn*compressEntrySize].pre[2] & 0x1) << 2)) == c);

      }else if(num <= 42){
        uint64_t tmp2 = (~(((compressedOcctable[p + fn*compressEntrySize].pre[1] >> 2) 
                            | ((compressedOcctable[p + fn*compressEntrySize].pre[2] & 0x1) << 62)) ^ _atgcMask)) & preMask[num - 22];

        // printBits((compressedOcctable[p + fn*compressEntrySize].pre[1] >> 2) 
        //                     | ((compressedOcctable[p + fn*compressEntrySize].pre[2] & 0x1) << 62));
        // printBits(_atgcMask);
        // printBits(tmp2);


        tmp2 = tmp2 & (tmp2 >> 1) & (tmp2 >> 2) & 0x1249249249249249;
        tmp2 = (tmp2 + (tmp2 >> 3)) & 0x71c71c71c71c71c7;
        tmp2 = (tmp2 + (tmp2 >> 6)) & 0x703f03f03f03f03f;
        tmp2 = (tmp2 + (tmp2 >> 12)) & 0x0fff000fff000fff;
        tmp2 = (tmp2 + (tmp2 >> 24)) & 0x0fff000000ffffff;
        tmp2 = (tmp2 + (tmp2 >> 48)) & 0x0000ffffffffffff;

        ret -= tmp2;

        // cout << tmp2 << "tmp3\n";

        uint64_t tmp3 = (~((compressedOcctable[p + fn*compressEntrySize].pre[2] >> 1) ^ _atgcMask));
        tmp3 = tmp3 & (tmp3 >> 1) & (tmp3 >> 2) & 0x1249249249249249;
        tmp3 = (tmp3 + (tmp3 >> 3)) & 0x71c71c71c71c71c7;
        tmp3 = (tmp3 + (tmp3 >> 6)) & 0x703f03f03f03f03f;
        tmp3 = (tmp3 + (tmp3 >> 12)) & 0x0fff000fff000fff;
        tmp3 = (tmp3 + (tmp3 >> 24)) & 0x0fff000000ffffff;
        tmp3 = (tmp3 + (tmp3 >> 48)) & 0x0000ffffffffffff;

        // cout << tmp3 << "tmp3\n";

        ret -= tmp3;
      }else if(num < 64){

        uint64_t tmp3 = (~((compressedOcctable[p + fn*compressEntrySize].pre[2] >> 1) ^ _atgcMask)) & preMask[num - 43];
        // cout << "\n";
        // printBits(compressedOcctable[p + fn*compressEntrySize].pre[2] >> 1);
        // printBits(_atgcMask);
        // printBits(tmp3);

        tmp3 = tmp3 & (tmp3 >> 1) & (tmp3 >> 2) & 0x1249249249249249;
        tmp3 = (tmp3 + (tmp3 >> 3)) & 0x71c71c71c71c71c7;
        tmp3 = (tmp3 + (tmp3 >> 6)) & 0x703f03f03f03f03f;
        tmp3 = (tmp3 + (tmp3 >> 12)) & 0x0fff000fff000fff;
        tmp3 = (tmp3 + (tmp3 >> 24)) & 0x0fff000000ffffff;
        tmp3 = (tmp3 + (tmp3 >> 48)) & 0x0000ffffffffffff;

        ret -= tmp3;
      }else if(num <= 85){
        uint64_t tmp = (~(compressedOcctable[p + fn*compressEntrySize].post[0] ^ _atgcMask)) & postMask[num - 64];

        tmp = tmp & (tmp >> 1) & (tmp >> 2) & 0x1249249249249249;
        tmp = (tmp + (tmp >> 3)) & 0x71c71c71c71c71c7;
        tmp = (tmp + (tmp >> 6)) & 0x703f03f03f03f03f;
        tmp = (tmp + (tmp >> 12)) & 0x0fff000fff000fff;
        tmp = (tmp + (tmp >> 24)) & 0x0fff000000ffffff;
        tmp = (tmp + (tmp >> 48)) & 0x0000ffffffffffff;

        ret += tmp;

      }else if(num == 86){
        uint64_t tmp = (~(compressedOcctable[p + fn*compressEntrySize].post[0] ^ _atgcMask));

        tmp = tmp & (tmp >> 1) & (tmp >> 2) & 0x1249249249249249;
        tmp = (tmp + (tmp >> 3)) & 0x71c71c71c71c71c7;
        tmp = (tmp + (tmp >> 6)) & 0x703f03f03f03f03f;
        tmp = (tmp + (tmp >> 12)) & 0x0fff000fff000fff;
        tmp = (tmp + (tmp >> 24)) & 0x0fff000000ffffff;
        tmp = (tmp + (tmp >> 48)) & 0x0000ffffffffffff;

        ret += tmp;

        ret += (((compressedOcctable[p + fn*compressEntrySize].post[0] >> 63) | ((compressedOcctable[p + fn*compressEntrySize].post[1] & 0x3) << 1)) == c);

      }else if(num <= 106){
        uint64_t tmp2 = (~(((compressedOcctable[p + fn*compressEntrySize].post[1] << 1) 
                            | ((compressedOcctable[p + fn*compressEntrySize].post[0]) >> 63)) ^ _atgcMask)) & postMask[num - 85];

        tmp2 = tmp2 & (tmp2 >> 1) & (tmp2 >> 2) & 0x1249249249249249;
        tmp2 = (tmp2 + (tmp2 >> 3)) & 0x71c71c71c71c71c7;
        tmp2 = (tmp2 + (tmp2 >> 6)) & 0x703f03f03f03f03f;
        tmp2 = (tmp2 + (tmp2 >> 12)) & 0x0fff000fff000fff;
        tmp2 = (tmp2 + (tmp2 >> 24)) & 0x0fff000000ffffff;
        tmp2 = (tmp2 + (tmp2 >> 48)) & 0x0000ffffffffffff;

        ret += tmp2;

        uint64_t tmp3 = (~((compressedOcctable[p + fn*compressEntrySize].post[0]) ^ _atgcMask));
        tmp3 = tmp3 & (tmp3 >> 1) & (tmp3 >> 2) & 0x1249249249249249;
        tmp3 = (tmp3 + (tmp3 >> 3)) & 0x71c71c71c71c71c7;
        tmp3 = (tmp3 + (tmp3 >> 6)) & 0x703f03f03f03f03f;
        tmp3 = (tmp3 + (tmp3 >> 12)) & 0x0fff000fff000fff;
        tmp3 = (tmp3 + (tmp3 >> 24)) & 0x0fff000000ffffff;
        tmp3 = (tmp3 + (tmp3 >> 48)) & 0x0000ffffffffffff;

        ret += tmp3;
      }else if(num < 128){

        uint64_t tmp = (~((compressedOcctable[p + fn*compressEntrySize].post[2] >> 0x1) ^ _atgcMask)) & postMask[num - 107];

        tmp = tmp & (tmp >> 1) & (tmp >> 2) & 0x1249249249249249;
        tmp = (tmp + (tmp >> 3)) & 0x71c71c71c71c71c7;
        tmp = (tmp + (tmp >> 6)) & 0x703f03f03f03f03f;
        tmp = (tmp + (tmp >> 12)) & 0x0fff000fff000fff;
        tmp = (tmp + (tmp >> 24)) & 0x0fff000000ffffff;
        tmp = (tmp + (tmp >> 48)) & 0x0000ffffffffffff;

        ret += tmp;

        ret += (((compressedOcctable[p + fn*compressEntrySize].post[1] >> 62)
                 | ((compressedOcctable[p + fn*compressEntrySize].post[2] & 0x1) << 2)) == c);

        uint64_t tmp2 = (~(((compressedOcctable[p + fn*compressEntrySize].post[1] << 1) 
                            | ((compressedOcctable[p + fn*compressEntrySize].post[0]) >> 63)) ^ _atgcMask));

        tmp2 = tmp2 & (tmp2 >> 1) & (tmp2 >> 2) & 0x1249249249249249;
        tmp2 = (tmp2 + (tmp2 >> 3)) & 0x71c71c71c71c71c7;
        tmp2 = (tmp2 + (tmp2 >> 6)) & 0x703f03f03f03f03f;
        tmp2 = (tmp2 + (tmp2 >> 12)) & 0x0fff000fff000fff;
        tmp2 = (tmp2 + (tmp2 >> 24)) & 0x0fff000000ffffff;
        tmp2 = (tmp2 + (tmp2 >> 48)) & 0x0000ffffffffffff;

        ret += tmp2;

        uint64_t tmp3 = (~((compressedOcctable[p + fn*compressEntrySize].post[0]) ^ _atgcMask));
        tmp3 = tmp3 & (tmp3 >> 1) & (tmp3 >> 2) & 0x1249249249249249;
        tmp3 = (tmp3 + (tmp3 >> 3)) & 0x71c71c71c71c71c7;
        tmp3 = (tmp3 + (tmp3 >> 6)) & 0x703f03f03f03f03f;
        tmp3 = (tmp3 + (tmp3 >> 12)) & 0x0fff000fff000fff;
        tmp3 = (tmp3 + (tmp3 >> 24)) & 0x0fff000000ffffff;
        tmp3 = (tmp3 + (tmp3 >> 48)) & 0x0000ffffffffffff;

        ret += tmp3;
      }


      

      return ret;




      // int p = i / oneEntryNum;
      // int j = (i % oneEntryNum)/32;

      // int ret = compressedOcctable[p + fn*compressEntrySize].atgc[c];

      // cout << "p : " << p << " j : " << j << " i%oneEntrynum : " << i % oneEntryNum <<  "\n";

      // if(j < 3){
      //   uint64_t tmp = (~(compressedOcctable[p + fn*compressEntrySize].pre[j] ^ atgcMask[c])) & preMask[(i % oneEntryNum) % 32];
      //   printf("%d\n",(i % oneEntryNum) % 32);

      //   print64bit(compressedOcctable[p + fn*compressEntrySize].pre[j]); printf(" %016lx tmp\n",tmp);
      //   print64bit(compressedOcctable[p + fn*compressEntrySize].pre[j]); printf(" %016lx tmp\n",tmp >> 1);
      //   tmp = (tmp & (tmp >> 1)) & 0x5555555555555555;// 0xaaaaaaaaaaaaaaaa;
      //   print64bit(compressedOcctable[p + fn*compressEntrySize].pre[j]); printf(" %016lx tmp\n",tmp);
        
      //   tmp = (tmp & 0x3333333333333333) + ((tmp >> 2) & 0x3333333333333333);
      //   tmp = (tmp & 0x0f0f0f0f0f0f0f0f) + ((tmp >> 4) & 0x0f0f0f0f0f0f0f0f);
      //   tmp = (tmp & 0x00ff00ff00ff00ff) + ((tmp >> 8) & 0x00ff00ff00ff00ff);
      //   tmp = (tmp & 0x0000ffff0000ffff) + ((tmp >> 16) & 0x0000ffff0000ffff);
      //   tmp = (tmp & 0x00000000ffffffff) + ((tmp >> 32) & 0x00000000ffffffff);
      //   cout << tmp << " first\n";
      //   ret -= tmp;

      //   for(int ii = j + 1;ii < 3;ii++){
      //     uint64_t tmp = (~(compressedOcctable[p + fn*compressEntrySize].pre[ii] ^ atgcMask[c]));
      //     print64bit(compressedOcctable[p + fn*compressEntrySize].pre[ii]); printf(" %016lx tmp\n",tmp);
      //     tmp = (tmp & (tmp >> 1)) & 0x5555555555555555;// 0xaaaaaaaaaaaaaaaa;
      //     print64bit(compressedOcctable[p + fn*compressEntrySize].pre[ii]); printf(" %016lx tmp\n",tmp);


      //     tmp = (tmp & 0x3333333333333333) + ((tmp >> 2) & 0x3333333333333333);
      //     tmp = (tmp & 0x0f0f0f0f0f0f0f0f) + ((tmp >> 4) & 0x0f0f0f0f0f0f0f0f);
      //     tmp = (tmp & 0x00ff00ff00ff00ff) + ((tmp >> 8) & 0x00ff00ff00ff00ff);
      //     tmp = (tmp & 0x0000ffff0000ffff) + ((tmp >> 16) & 0x0000ffff0000ffff);
      //     tmp = (tmp & 0x00000000ffffffff) + ((tmp >> 32) & 0x00000000ffffffff);
      //     cout << tmp << " first\n";
      //     ret -= tmp;
      //   }
      // }else{
      //   uint64_t tmp = (~(compressedOcctable[p + fn*compressEntrySize].post[j - 3] ^ atgcMask[c])) & postMask[(i % oneEntryNum) % 32];
      //   tmp = (tmp & (tmp >> 1)) & 0x5555555555555555;
      //   tmp = (tmp & 0x3333333333333333) + ((tmp >> 2) & 0x3333333333333333);
      //   tmp = (tmp & 0x0f0f0f0f0f0f0f0f) + ((tmp >> 4) & 0x0f0f0f0f0f0f0f0f);
      //   tmp = (tmp & 0x00ff00ff00ff00ff) + ((tmp >> 8) & 0x00ff00ff00ff00ff);
      //   tmp = (tmp & 0x0000ffff0000ffff) + ((tmp >> 16) & 0x0000ffff0000ffff);
      //   tmp = (tmp & 0x00000000ffffffff) + ((tmp >> 32) & 0x00000000ffffffff);
      //   ret += tmp;

      //   for(int ii = 0;ii < j - 3;ii++){
      //     uint64_t tmp = (~(compressedOcctable[p + fn*compressEntrySize].pre[ii] ^ atgcMask[c]));          
      //     tmp = (tmp & (tmp >> 1)) & 0x5555555555555555;
      //     tmp = (tmp & 0x3333333333333333) + ((tmp >> 2) & 0x3333333333333333);
      //     tmp = (tmp & 0x0f0f0f0f0f0f0f0f) + ((tmp >> 4) & 0x0f0f0f0f0f0f0f0f);
      //     tmp = (tmp & 0x00ff00ff00ff00ff) + ((tmp >> 8) & 0x00ff00ff00ff00ff);
      //     tmp = (tmp & 0x0000ffff0000ffff) + ((tmp >> 16) & 0x0000ffff0000ffff);
      //     tmp = (tmp & 0x00000000ffffffff) + ((tmp >> 32) & 0x00000000ffffffff);
      //     ret += tmp;
      //   }
      // }


      // return ret;
    }




    // c以上の文字数を取る
    int getOccMorethanCFromCompressed(int i,unsigned char c,int fn){
      //      uint64_t _atgcMask = atgcMask[c];

      if(i == n){
        int p = (n - 1) / oneEntryNum;
        int ret = 0;

        for(int j = c;j < 4;j++){
          ret += compressedOcctable[p + fn*compressEntrySize].atgc[j];
        }

        
        uint64_t tmp = 0;//(~((compressedOcctable[p + fn*compressEntrySize].post[2] >> 0x1) ^ _atgcMask));
        uint64_t _line = (compressedOcctable[p + fn*compressEntrySize].post[2] >> 0x1);
        uint64_t tmp2 = 0;
        uint64_t _line2 = ((compressedOcctable[p + fn*compressEntrySize].post[1] >> 2) 
                           | ((compressedOcctable[p + fn*compressEntrySize].post[2] & 0x1) << 62));
        uint64_t tmp3 = 0;
        uint64_t  _line3 = (compressedOcctable[p + fn*compressEntrySize].post[0]);
        for(int j = c;j < 4;j++){
          uint64_t _tmp = (~(_line ^ atgcMask[j]));
          tmp |= _tmp & (_tmp >> 1) & (_tmp >> 2) & 0x1249249249249249;

          uint64_t _tmp2 = (~(_line2 ^ atgcMask[j]));
          tmp2 |= _tmp2 & (_tmp2 >> 1) & (_tmp2 >> 2) & 0x1249249249249249;

          uint64_t _tmp3 = (~(_line3 ^ atgcMask[j]));
          tmp3 |= _tmp3 & (_tmp3 >> 1) & (_tmp3 >> 2) & 0x1249249249249249;
        }

        //        tmp = tmp & (tmp >> 1) & (tmp >> 2) & 0x1249249249249249;
        // printBit(tmp);

        tmp = (tmp + (tmp >> 3)) & 0x71c71c71c71c71c7;
        tmp = (tmp + (tmp >> 6)) & 0x703f03f03f03f03f;
        tmp = (tmp + (tmp >> 12)) & 0x0fff000fff000fff;
        tmp = (tmp + (tmp >> 24)) & 0x0fff000000ffffff;
        tmp = (tmp + (tmp >> 48)) & 0x0000ffffffffffff;

        ret += tmp;

        //        tmp2 = tmp2 & (tmp2 >> 1) & (tmp2 >> 2) & 0x1249249249249249;

        tmp2 = (tmp2 + (tmp2 >> 3)) & 0x71c71c71c71c71c7;
        tmp2 = (tmp2 + (tmp2 >> 6)) & 0x703f03f03f03f03f;
        tmp2 = (tmp2 + (tmp2 >> 12)) & 0x0fff000fff000fff;
        tmp2 = (tmp2 + (tmp2 >> 24)) & 0x0fff000000ffffff;
        tmp2 = (tmp2 + (tmp2 >> 48)) & 0x0000ffffffffffff;

        ret += tmp2;

        //        tmp3 = tmp3 & (tmp3 >> 1) & (tmp3 >> 2) & 0x1249249249249249;

        tmp3 = (tmp3 + (tmp3 >> 3)) & 0x71c71c71c71c71c7;
        tmp3 = (tmp3 + (tmp3 >> 6)) & 0x703f03f03f03f03f;
        tmp3 = (tmp3 + (tmp3 >> 12)) & 0x0fff000fff000fff;
        tmp3 = (tmp3 + (tmp3 >> 24)) & 0x0fff000000ffffff;
        tmp3 = (tmp3 + (tmp3 >> 48)) & 0x0000ffffffffffff;


        uint64_t ccc = ((compressedOcctable[p + fn*compressEntrySize].post[0] >> 63) 
                        | ((compressedOcctable[p + fn*compressEntrySize].post[1] & 0x3) << 1));
        ret += ((ccc >= c) && (ccc < 4));

        ret += tmp3;
        return ret;
      }

      int p = i / oneEntryNum;
      int ret = 0;//compressedOcctable[p + fn*compressEntrySize].atgc[c];
      for(int j = c;j < 4;j++){
        ret += compressedOcctable[p + fn*compressEntrySize].atgc[j];
      }

      int num = (i % oneEntryNum);


      // cout << "num : " << num << " i : " << i << " c : " << (int)c << "\n";
      // cout << ret << " # ret\n";

      if(num < 21){
        uint64_t tmp = 0;//(~((compressedOcctable[p + fn*compressEntrySize].post[2] >> 0x1) ^ _atgcMask));
        uint64_t _line = compressedOcctable[p + fn*compressEntrySize].pre[0];
        uint64_t tmp2 = 0;
        uint64_t _line2 = ((compressedOcctable[p + fn*compressEntrySize].pre[1] << 1) 
                           | (compressedOcctable[p + fn*compressEntrySize].pre[0] >> 63));
        uint64_t tmp3 = 0;
        uint64_t  _line3 = (compressedOcctable[p + fn*compressEntrySize].pre[2] >> 1);

        for(int j = c;j < 4;j++){
          uint64_t _tmp = ((~(_line ^ atgcMask[j])) & preMask[num]);
          tmp |= _tmp & (_tmp >> 1) & (_tmp >> 2) & 0x1249249249249249;

          uint64_t _tmp2 = (~(_line2 ^ atgcMask[j]));
          tmp2 |= _tmp2 & (_tmp2 >> 1) & (_tmp2 >> 2) & 0x1249249249249249;

          uint64_t _tmp3 = (~(_line3 ^ atgcMask[j]));
          tmp3 |= _tmp3 & (_tmp3 >> 1) & (_tmp3 >> 2) & 0x1249249249249249;
        }

        // printBit(tmp);
        // printBit(tmp2);
        // printBit(tmp3);

        // uint64_t tmp = (~(compressedOcctable[p + fn*compressEntrySize].pre[0] ^ _atgcMask)) & preMask[num];
        // uint64_t tmp2 = (~(((compressedOcctable[p + fn*compressEntrySize].pre[1] << 1) 
        //                     | (compressedOcctable[p + fn*compressEntrySize].pre[0] >> 63)) ^ _atgcMask));
        // uint64_t tmp3 = (~((compressedOcctable[p + fn*compressEntrySize].pre[2] >> 1) ^ _atgcMask));


        //        tmp = tmp & (tmp >> 1) & (tmp >> 2) & 0x1249249249249249;
        tmp = (tmp + (tmp >> 3)) & 0x71c71c71c71c71c7;
        tmp = (tmp + (tmp >> 6)) & 0x703f03f03f03f03f;
        tmp = (tmp + (tmp >> 12)) & 0x0fff000fff000fff;
        tmp = (tmp + (tmp >> 24)) & 0x0fff000000ffffff;
        tmp = (tmp + (tmp >> 48)) & 0x0000ffffffffffff;

        // cout << tmp << " tmp\n";


        ret -= tmp;

        //        tmp2 = tmp2 & (tmp2 >> 1) & (tmp2 >> 2) & 0x1249249249249249;
        tmp2 = (tmp2 + (tmp2 >> 3)) & 0x71c71c71c71c71c7;
        tmp2 = (tmp2 + (tmp2 >> 6)) & 0x703f03f03f03f03f;
        tmp2 = (tmp2 + (tmp2 >> 12)) & 0x0fff000fff000fff;
        tmp2 = (tmp2 + (tmp2 >> 24)) & 0x0fff000000ffffff;
        tmp2 = (tmp2 + (tmp2 >> 48)) & 0x0000ffffffffffff;

        // cout << tmp2 << " tmp2\n";


        ret -= tmp2;

        // tmp3 = tmp3 & (tmp3 >> 1) & (tmp3 >> 2) & 0x1249249249249249;
        tmp3 = (tmp3 + (tmp3 >> 3)) & 0x71c71c71c71c71c7;
        tmp3 = (tmp3 + (tmp3 >> 6)) & 0x703f03f03f03f03f;
        tmp3 = (tmp3 + (tmp3 >> 12)) & 0x0fff000fff000fff;
        tmp3 = (tmp3 + (tmp3 >> 24)) & 0x0fff000000ffffff;
        tmp3 = (tmp3 + (tmp3 >> 48)) & 0x0000ffffffffffff;

        // cout << tmp3 << " tmp3\n";

        ret -= tmp3;


        // cout << ret << " returnval\n";

        uint64_t ccc = ((compressedOcctable[p + fn*compressEntrySize].pre[1] >> 62) | ((compressedOcctable[p + fn*compressEntrySize].pre[2] & 0x1) << 0x2));
        // cout << ccc << " ccc\n";

        ret -= (ccc >= c && ccc < 4);

        // cout << ret << " returnval\n";

      }else if(num == 21){
        uint64_t tmp2 = 0;
        uint64_t _line2 = ((compressedOcctable[p + fn*compressEntrySize].pre[1] << 1) 
                           | (compressedOcctable[p + fn*compressEntrySize].pre[0] >> 63));
        uint64_t tmp3 = 0;
        uint64_t  _line3 = (compressedOcctable[p + fn*compressEntrySize].pre[2] >> 1);

        for(int j = c;j < 4;j++){
          uint64_t _tmp2 = (~(_line2 ^ atgcMask[j]));
          tmp2 |= _tmp2 & (_tmp2 >> 1) & (_tmp2 >> 2) & 0x1249249249249249;

          uint64_t _tmp3 = (~(_line3 ^ atgcMask[j]));
          tmp3 |= _tmp3 & (_tmp3 >> 1) & (_tmp3 >> 2) & 0x1249249249249249;
        }

        // uint64_t tmp2 = (~(((compressedOcctable[p + fn*compressEntrySize].pre[1] << 1) 
        //                     | (compressedOcctable[p + fn*compressEntrySize].pre[0] >> 63)) ^ _atgcMask));

        //        tmp2 = tmp2 & (tmp2 >> 1) & (tmp2 >> 2) & 0x1249249249249249;
        tmp2 = (tmp2 + (tmp2 >> 3)) & 0x71c71c71c71c71c7;
        tmp2 = (tmp2 + (tmp2 >> 6)) & 0x703f03f03f03f03f;
        tmp2 = (tmp2 + (tmp2 >> 12)) & 0x0fff000fff000fff;
        tmp2 = (tmp2 + (tmp2 >> 24)) & 0x0fff000000ffffff;
        tmp2 = (tmp2 + (tmp2 >> 48)) & 0x0000ffffffffffff;

        ret -= tmp2;

        //        uint64_t tmp3 = (~((compressedOcctable[p + fn*compressEntrySize].pre[2] >> 1) ^ _atgcMask));
        //        tmp3 = tmp3 & (tmp3 >> 1) & (tmp3 >> 2) & 0x1249249249249249;
        tmp3 = (tmp3 + (tmp3 >> 3)) & 0x71c71c71c71c71c7;
        tmp3 = (tmp3 + (tmp3 >> 6)) & 0x703f03f03f03f03f;
        tmp3 = (tmp3 + (tmp3 >> 12)) & 0x0fff000fff000fff;
        tmp3 = (tmp3 + (tmp3 >> 24)) & 0x0fff000000ffffff;
        tmp3 = (tmp3 + (tmp3 >> 48)) & 0x0000ffffffffffff;

        ret -= tmp3;

        uint64_t ccc = ((compressedOcctable[p + fn*compressEntrySize].pre[1] >> 62) | ((compressedOcctable[p + fn*compressEntrySize].pre[2] & 0x1) << 2));
        ret -= (ccc >= c && ccc < 4);

      }else if(num <= 42){
        uint64_t tmp2 = 0;
        uint64_t _line2 = ((compressedOcctable[p + fn*compressEntrySize].pre[1] >> 2) 
                           | ((compressedOcctable[p + fn*compressEntrySize].pre[2] & 0x1) << 62));
        uint64_t tmp3 = 0;
        uint64_t  _line3 = (compressedOcctable[p + fn*compressEntrySize].pre[2] >> 1);

        for(int j = c;j < 4;j++){
          uint64_t _tmp2 = (~(_line2 ^ atgcMask[j])) & preMask[num - 22];
          tmp2 |= _tmp2 & (_tmp2 >> 1) & (_tmp2 >> 2) & 0x1249249249249249;

          uint64_t _tmp3 = (~(_line3 ^ atgcMask[j]));
          tmp3 |= _tmp3 & (_tmp3 >> 1) & (_tmp3 >> 2) & 0x1249249249249249;
        }

        // uint64_t tmp2 = (~(((compressedOcctable[p + fn*compressEntrySize].pre[1] >> 2) 
        //                     | ((compressedOcctable[p + fn*compressEntrySize].pre[2] & 0x1) << 62)) ^ _atgcMask)) & preMask[num - 22];

        // printBits((compressedOcctable[p + fn*compressEntrySize].pre[1] >> 2) 
        //                     | ((compressedOcctable[p + fn*compressEntrySize].pre[2] & 0x1) << 62));
        // printBits(_atgcMask);
        // printBits(tmp2);


        // tmp2 = tmp2 & (tmp2 >> 1) & (tmp2 >> 2) & 0x1249249249249249;
        tmp2 = (tmp2 + (tmp2 >> 3)) & 0x71c71c71c71c71c7;
        tmp2 = (tmp2 + (tmp2 >> 6)) & 0x703f03f03f03f03f;
        tmp2 = (tmp2 + (tmp2 >> 12)) & 0x0fff000fff000fff;
        tmp2 = (tmp2 + (tmp2 >> 24)) & 0x0fff000000ffffff;
        tmp2 = (tmp2 + (tmp2 >> 48)) & 0x0000ffffffffffff;

        ret -= tmp2;

        // cout << tmp2 << "tmp3\n";

        //        uint64_t tmp3 = (~((compressedOcctable[p + fn*compressEntrySize].pre[2] >> 1) ^ _atgcMask));
        // tmp3 = tmp3 & (tmp3 >> 1) & (tmp3 >> 2) & 0x1249249249249249;
        tmp3 = (tmp3 + (tmp3 >> 3)) & 0x71c71c71c71c71c7;
        tmp3 = (tmp3 + (tmp3 >> 6)) & 0x703f03f03f03f03f;
        tmp3 = (tmp3 + (tmp3 >> 12)) & 0x0fff000fff000fff;
        tmp3 = (tmp3 + (tmp3 >> 24)) & 0x0fff000000ffffff;
        tmp3 = (tmp3 + (tmp3 >> 48)) & 0x0000ffffffffffff;

        // cout << tmp3 << "tmp3\n";

        ret -= tmp3;
      }else if(num < 64){

        uint64_t tmp3 = 0;
        uint64_t  _line3 = (compressedOcctable[p + fn*compressEntrySize].pre[2] >> 1);

        for(int j = c;j < 4;j++){
          uint64_t _tmp3 = (~(_line3 ^ atgcMask[j])) & preMask[num - 43];
          tmp3 |= _tmp3 & (_tmp3 >> 1) & (_tmp3 >> 2) & 0x1249249249249249;
        }

        // uint64_t tmp3 = (~((compressedOcctable[p + fn*compressEntrySize].pre[2] >> 1) ^ _atgcMask)) & preMask[num - 43];
        // cout << "\n";
        // printBits(compressedOcctable[p + fn*compressEntrySize].pre[2] >> 1);
        // printBits(_atgcMask);
        // printBits(tmp3);

        // tmp3 = tmp3 & (tmp3 >> 1) & (tmp3 >> 2) & 0x1249249249249249;
        tmp3 = (tmp3 + (tmp3 >> 3)) & 0x71c71c71c71c71c7;
        tmp3 = (tmp3 + (tmp3 >> 6)) & 0x703f03f03f03f03f;
        tmp3 = (tmp3 + (tmp3 >> 12)) & 0x0fff000fff000fff;
        tmp3 = (tmp3 + (tmp3 >> 24)) & 0x0fff000000ffffff;
        tmp3 = (tmp3 + (tmp3 >> 48)) & 0x0000ffffffffffff;

        ret -= tmp3;
      }else if(num <= 85){
        uint64_t tmp = 0;
        uint64_t _line = compressedOcctable[p + fn*compressEntrySize].post[0];

        for(int j = c;j < 4;j++){
          uint64_t _tmp0 = (~(_line ^ atgcMask[j])) & postMask[num - 64];
          tmp |= _tmp0 & (_tmp0 >> 1) & (_tmp0 >> 2) & 0x1249249249249249;
        }

        // uint64_t tmp = (~(compressedOcctable[p + fn*compressEntrySize].post[0] ^ _atgcMask)) & postMask[num - 64];

        // tmp = tmp & (tmp >> 1) & (tmp >> 2) & 0x1249249249249249;

        tmp = (tmp + (tmp >> 3)) & 0x71c71c71c71c71c7;
        tmp = (tmp + (tmp >> 6)) & 0x703f03f03f03f03f;
        tmp = (tmp + (tmp >> 12)) & 0x0fff000fff000fff;
        tmp = (tmp + (tmp >> 24)) & 0x0fff000000ffffff;
        tmp = (tmp + (tmp >> 48)) & 0x0000ffffffffffff;

        ret += tmp;

      }else if(num == 86){
        uint64_t tmp = 0;
        uint64_t _line = compressedOcctable[p + fn*compressEntrySize].post[0];

        for(int j = c;j < 4;j++){
          uint64_t _tmp0 = (~(_line ^ atgcMask[j]));
          tmp |= _tmp0 & (_tmp0 >> 1) & (_tmp0 >> 2) & 0x1249249249249249;
        }

        // uint64_t tmp = (~(compressedOcctable[p + fn*compressEntrySize].post[0] ^ _atgcMask));
        // tmp = tmp & (tmp >> 1) & (tmp >> 2) & 0x1249249249249249;

        tmp = (tmp + (tmp >> 3)) & 0x71c71c71c71c71c7;
        tmp = (tmp + (tmp >> 6)) & 0x703f03f03f03f03f;
        tmp = (tmp + (tmp >> 12)) & 0x0fff000fff000fff;
        tmp = (tmp + (tmp >> 24)) & 0x0fff000000ffffff;
        tmp = (tmp + (tmp >> 48)) & 0x0000ffffffffffff;

        ret += tmp;

        uint64_t ccc = ((compressedOcctable[p + fn*compressEntrySize].post[0] >> 63) | ((compressedOcctable[p + fn*compressEntrySize].post[1] & 0x3) << 1));
        ret += (ccc >= c && ccc < 4);

      }else if(num <= 106){
        uint64_t tmp2 = 0;
        uint64_t _line2 = ((compressedOcctable[p + fn*compressEntrySize].post[1] << 1) 
                          | ((compressedOcctable[p + fn*compressEntrySize].post[0]) >> 63));
        uint64_t tmp3 = 0;
        uint64_t _line3 = (compressedOcctable[p + fn*compressEntrySize].post[0]);

        for(int j = c;j < 4;j++){
          uint64_t _tmp2 = (~(_line2 ^ atgcMask[j])) & postMask[num - 85];
          tmp2 |= _tmp2 & (_tmp2 >> 1) & (_tmp2 >> 2) & 0x1249249249249249;
          uint64_t _tmp3 = (~(_line3 ^ atgcMask[j]));
          tmp3 |= _tmp3 & (_tmp3 >> 1) & (_tmp3 >> 2) & 0x1249249249249249;
        }


        // uint64_t tmp2 = (~(((compressedOcctable[p + fn*compressEntrySize].post[1] << 1) 
        //                     | ((compressedOcctable[p + fn*compressEntrySize].post[0]) >> 63)) ^ _atgcMask)) & postMask[num - 85];

        // tmp2 = tmp2 & (tmp2 >> 1) & (tmp2 >> 2) & 0x1249249249249249;
        tmp2 = (tmp2 + (tmp2 >> 3)) & 0x71c71c71c71c71c7;
        tmp2 = (tmp2 + (tmp2 >> 6)) & 0x703f03f03f03f03f;
        tmp2 = (tmp2 + (tmp2 >> 12)) & 0x0fff000fff000fff;
        tmp2 = (tmp2 + (tmp2 >> 24)) & 0x0fff000000ffffff;
        tmp2 = (tmp2 + (tmp2 >> 48)) & 0x0000ffffffffffff;

        ret += tmp2;

        // uint64_t tmp3 = (~((compressedOcctable[p + fn*compressEntrySize].post[0]) ^ _atgcMask));
        // tmp3 = tmp3 & (tmp3 >> 1) & (tmp3 >> 2) & 0x1249249249249249;
        tmp3 = (tmp3 + (tmp3 >> 3)) & 0x71c71c71c71c71c7;
        tmp3 = (tmp3 + (tmp3 >> 6)) & 0x703f03f03f03f03f;
        tmp3 = (tmp3 + (tmp3 >> 12)) & 0x0fff000fff000fff;
        tmp3 = (tmp3 + (tmp3 >> 24)) & 0x0fff000000ffffff;
        tmp3 = (tmp3 + (tmp3 >> 48)) & 0x0000ffffffffffff;

        ret += tmp3;
      }else if(num < 128){

        uint64_t tmp = 0;
        uint64_t _line = (compressedOcctable[p + fn*compressEntrySize].post[2] >> 0x1);
        uint64_t tmp2 = 0;
        uint64_t _line2 = ((compressedOcctable[p + fn*compressEntrySize].post[1] << 1) 
                           | ((compressedOcctable[p + fn*compressEntrySize].post[0]) >> 63));
        uint64_t tmp3 = 0;
        uint64_t _line3 = (compressedOcctable[p + fn*compressEntrySize].post[0]);

        for(int j = c;j < 4;j++){
          uint64_t _tmp = (~(_line ^ atgcMask[j])) & postMask[num - 107];
          tmp |= _tmp & (_tmp >> 1) & (_tmp >> 2) & 0x1249249249249249;
          uint64_t _tmp2 = (~(_line2 ^ atgcMask[j]));
          tmp2 |= _tmp2 & (_tmp2 >> 1) & (_tmp2 >> 2) & 0x1249249249249249;
          uint64_t _tmp3 = (~(_line3 ^ atgcMask[j]));
          tmp3 |= _tmp3 & (_tmp3 >> 1) & (_tmp3 >> 2) & 0x1249249249249249;
        }


        // printBit(tmp);
        // printBit(tmp2);
        // printBit(tmp3);


        // uint64_t tmp = (~((compressedOcctable[p + fn*compressEntrySize].post[2] >> 0x1) ^ _atgcMask)) & postMask[num - 107];

        // tmp = tmp & (tmp >> 1) & (tmp >> 2) & 0x1249249249249249;
        tmp = (tmp + (tmp >> 3)) & 0x71c71c71c71c71c7;
        tmp = (tmp + (tmp >> 6)) & 0x703f03f03f03f03f;
        tmp = (tmp + (tmp >> 12)) & 0x0fff000fff000fff;
        tmp = (tmp + (tmp >> 24)) & 0x0fff000000ffffff;
        tmp = (tmp + (tmp >> 48)) & 0x0000ffffffffffff;

        ret += tmp;


        uint64_t ccc = ((compressedOcctable[p + fn*compressEntrySize].post[1] >> 62)
                        | ((compressedOcctable[p + fn*compressEntrySize].post[2] & 0x1) << 2));
        ret += (ccc >= c && ccc < 4);

        // uint64_t tmp2 = (~(((compressedOcctable[p + fn*compressEntrySize].post[1] << 1) 
        //                     | ((compressedOcctable[p + fn*compressEntrySize].post[0]) >> 63)) ^ _atgcMask));

        // tmp2 = tmp2 & (tmp2 >> 1) & (tmp2 >> 2) & 0x1249249249249249;
        tmp2 = (tmp2 + (tmp2 >> 3)) & 0x71c71c71c71c71c7;
        tmp2 = (tmp2 + (tmp2 >> 6)) & 0x703f03f03f03f03f;
        tmp2 = (tmp2 + (tmp2 >> 12)) & 0x0fff000fff000fff;
        tmp2 = (tmp2 + (tmp2 >> 24)) & 0x0fff000000ffffff;
        tmp2 = (tmp2 + (tmp2 >> 48)) & 0x0000ffffffffffff;

        ret += tmp2;

        // uint64_t tmp3 = (~((compressedOcctable[p + fn*compressEntrySize].post[0]) ^ _atgcMask));
        // tmp3 = tmp3 & (tmp3 >> 1) & (tmp3 >> 2) & 0x1249249249249249;
        tmp3 = (tmp3 + (tmp3 >> 3)) & 0x71c71c71c71c71c7;
        tmp3 = (tmp3 + (tmp3 >> 6)) & 0x703f03f03f03f03f;
        tmp3 = (tmp3 + (tmp3 >> 12)) & 0x0fff000fff000fff;
        tmp3 = (tmp3 + (tmp3 >> 24)) & 0x0fff000000ffffff;
        tmp3 = (tmp3 + (tmp3 >> 48)) & 0x0000ffffffffffff;

        ret += tmp3;
      }

      return ret;
    }


    OCC(const char* dirname,int fn,C_OCC& cocc,string fname):c_occ(cocc){
      string filename = string(dirname) + "/" + fname;
      FILE* fp;
      this->fnum = fn;

      if((fp = fopen(filename.c_str(), "rb")) == NULL ) {
        printf("cannot open output file %s\n",filename.c_str());
        exit(EXIT_FAILURE);
      }

      int n;
      if(fseek(fp, 0, SEEK_END) == 0) {
        n = ftell(fp)/sizeof(oneEntry);
        rewind(fp);
        if(n < 0) {
          fprintf(stderr, "%s: Cannot ftell `%s': ", "C", filename.c_str());
          perror(NULL);
          exit(EXIT_FAILURE);
        }
      } else {
        fprintf(stderr, "%s: Cannot fseek `%s': ", "C", filename.c_str());
        perror(NULL);
        exit(EXIT_FAILURE);
      }

      compressedOcctable = make_unique<oneEntry[]>(n);      

      if(fread(compressedOcctable.get(), sizeof(oneEntry), (size_t)n, fp) != (size_t)n) {
        fprintf(stderr, "%s: %s `%s': ",
                "C",
                (ferror(fp) || !feof(fp)) ? "Cannot read from" : "Unexpected EOF in",
                filename.c_str());
        perror(NULL);
        exit(EXIT_FAILURE);
      }

      this->fnum = fn;
      this->compressEntrySize = n/fn;
      cout << compressEntrySize << " compressEntrysize " << compressEntrySize*oneEntryNum << " n\n";

      this->n = this->compressEntrySize*oneEntryNum;

      fclose(fp);
    }

    // nはひとつの断片の大きさ
    OCC(unique_ptr< unsigned char[]>& FBWT,int n,int fn,C_OCC& cocc):c_occ(cocc){
    
      compressEntrySize = (n % oneEntryNum == 0) ? (n / oneEntryNum) : (n / oneEntryNum + 1);
      //      occtable = std::make_unique<int[]>((n + 1)*CHAR_N*fn);

      this->n = n;
      this->fnum = fn;

      constructOccTableStructure(FBWT,fn,n,compressEntrySize);

    }


    int getN(){
      return n;
    }

    void printOCCCompressed(int c){
      for(int i = 0;i <= n;i++){
        cout << i << ",";
        for(int f = 0;f < fnum;f++){
          cout << getOccFromCompressed(i,c,f) << ",";
        }
        cout << "\n";
      }
    }

    void printOCCCompressed(){
      for(int i = 0;i <= n;i++){
        cout << i << ",";
        for(int f = 0;f < fnum;f++){
          for(int c = 0;c < CHAR_N - 1;c++)
            cout << getOccFromCompressed(i,c,f) << ",";
        }
        cout << "\n";
      }
    }

    //getOccMorethanCFromCompressed(int i,unsigned char c,int fn){
    void printMorethanC(){
      for(int i = 0;i <= n;i++){
        cout << i << ",";
        for(int f = 0;f < fnum;f++){
          for(int c = 0;c < CHAR_N - 1;c++)
            cout << getOccMorethanCFromCompressed(i,c,f) << ",";
        }
        cout << "\n";
      }
    }

    void printMorethanC(int c){
      for(int i = 0;i <= n;i++){
        cout << i << ",";
        for(int f = 0;f < fnum;f++){
          cout << getOccMorethanCFromCompressed(i,c,f) << ",";
        }
        cout << "\n";
      }
    }


    // void printOCC(){
    //   for(int f = 0;f < fnum;f++){
    //     for(int i = 0;i <= n;i++){
    //       for(int c = 0;c < CHAR_N;c++)
    //         cout << occtable[f*(n + 1)*CHAR_N + i*CHAR_N + c] << " ";
    //       cout << "\n";
    //     }
    //     cout << "---------------------\n";
    //   }
    // }

    void outputFile(const char* dirname,string fname){
      FILE* fp;
      string filename = string(dirname) + "/" + fname;
      if((fp = fopen(filename.c_str(), "wb")) == NULL ) {
        printf("cannot open output file %s\n",filename.c_str());
        exit(EXIT_FAILURE);
      }

      fwrite(compressedOcctable.get() , sizeof(oneEntry), compressEntrySize*fnum, fp);
      fclose(fp);
    }
    
  };

  class tuple{
  public:
    int left;
    int right;
  };

  int power(int d,int i){
    if(i == 0) return 1;
    else return d*power(d,i - 1);
  }


  void _constructKMR(int key,int depth,int maxdepth,int fn){
    if(maxdepth == depth){
      KMR[key + kmerOneFnSize*fn].left = 1;
      KMR[key + kmerOneFnSize*fn].right = 0;
      return;
    }
    for(int i = 0;i < 4;i++){
      _constructKMR(key | (i << 2*depth), depth + 1, maxdepth,(fn + 1)%fnum);
    }
  }

  void _constructKMR(C_OCC& c_occ,const pair<int,int>& interval,int key,int depth,const int maxdepth,int fn){
    if(interval.first > interval.second){
      _constructKMR(key, depth, maxdepth,fn);
      return;
    }
    if(maxdepth == depth){

      // printf("fn : %d, key : %3x, interval : [%d:%d]\n",fn,key,interval.first,interval.second);

      KMR[key + kmerOneFnSize*fn].left = interval.first;
      KMR[key + kmerOneFnSize*fn].right = interval.second;

      return;
    }
    for(int i = 0;i < 4;i++){
      pair<int,int> in(interval);
      c_occ.updateIntervalCompressed(in, i, fn);
      _constructKMR(c_occ, in, key | (i << 2*depth), depth + 1, maxdepth,(fn + 1)%fnum);
    }
  }

  void constructKMR(C_OCC& c_occ,int num){
    kmerOneFnSize = power(4,num);
    KMR = make_unique<tuple[]>(kmerOneFnSize*fnum);
    
    for(int f = 0;f < fnum;f++){
      for(int i = 0;i < 4;i++){
        pair<int,int> interval(c->getC(i,f),c->getC(i + 1,f) - 1);

        // cout << "char : " << i << ", interval[" << interval.first << ":" << interval.second << "]\n";

        _constructKMR(c_occ, interval, i, 1, num,(f + 1)%fnum);
      }
    }
    // for(int f = 0;f < fnum;f++){
    //   for(int i = 1;i < 5;i++){
    //     pair<int,int> interval(c->getC(i,f),c->getC(i + 1,f) - 1);
    //     _constructKMR(c_occ, interval, i - 1, 1, num,(f + 1)%fnum);
    //   }
    // }
    // printKMR();
  }

  void _constructRKMR(int key,int depth,int maxdepth,int fn){
    if(maxdepth == depth){
      RKMR[key + kmerOneFnSize*fn].left = 1;
      RKMR[key + kmerOneFnSize*fn].right = 0;
      return;
    }
    for(int i = 0;i < 4;i++){
      _constructRKMR((key << 2) | i, depth + 1, maxdepth,(fn + 1)%fnum);
    }
  }

  void _constructRKMR(C_OCC& c_occ,const pair<int,int>& interval,int key,int depth,const int maxdepth,int fn){
    if(interval.first > interval.second){
      _constructRKMR(key, depth, maxdepth,fn);
      return;
    }
    if(maxdepth == depth){
      RKMR[key + kmerOneFnSize*fn].left = interval.first;
      RKMR[key + kmerOneFnSize*fn].right = interval.second;

      return;
    }
    for(int i = 0;i < 4;i++){
      pair<int,int> in(interval);
      c_occ.updateRIntervalCompressed(in, i, fn);
      _constructRKMR(c_occ, in, (key << 2) | i, depth + 1, maxdepth,(fn + 1)%fnum);
    }
  }

  void constructRKMR(C_OCC& c_occ,int num){
    kmerOneFnSize = power(4,num);
    RKMR = make_unique<tuple[]>(kmerOneFnSize*fnum);
    
    for(int f = 0;f < fnum;f++){
      for(int i = 0;i < 4;i++){
        pair<int,int> interval(rc->getC(i,f),rc->getC(i + 1,f) - 1);
        _constructRKMR(c_occ, interval, i, 1, num,(f + 1)%fnum);
      }
    }
    // printKMR();
  }

  void printKMR(){
    char vv[4] = {'A','C','G','T'};
    cout << "=========================== KMR ===========================\n";
    for(int i = 0;i < kmerOneFnSize;i++){
      for(int j = kmerSize - 1;j >= 0;j--){
        cout << vv[((i >> j*2) & 0x3)];
      }
      for(int f = 0;f < fnum;f++){
        cout << " : [" << KMR[i + kmerOneFnSize*f].left << ":" << KMR[i + kmerOneFnSize*f].right << "]\t";
      }
      cout << "\n";
    }
  }

  void outputKMR(const char* dirname){
    FILE* fp;
    string filename = string(dirname) + "/" + kmrfname;
    if((fp = fopen(filename.c_str(), "wb")) == NULL ) {
      printf("cannot open output file %s\n",filename.c_str());
      exit(EXIT_FAILURE);
    }

    fwrite(KMR.get() , sizeof(tuple), kmerOneFnSize*fnum, fp);
    fclose(fp);

    int o[2] = {kmerSize,kmerOneFnSize};

    filename = string(dirname) + "/" + kmrsizefname;
    if((fp = fopen(filename.c_str(), "wb")) == NULL ) {
      printf("cannot open output file %s\n",filename.c_str());
      exit(EXIT_FAILURE);
    }
    fwrite(o , sizeof(int), 2, fp);
    fclose(fp);
  }

  void outputRKMR(const char* dirname){
    FILE* fp;
    string filename = string(dirname) + "/" + rkmrfname;
    if((fp = fopen(filename.c_str(), "wb")) == NULL ) {
      printf("cannot open output file %s\n",filename.c_str());
      exit(EXIT_FAILURE);
    }
    fwrite(RKMR.get() , sizeof(tuple), kmerOneFnSize*fnum, fp);
    fclose(fp);
  }

  void inputKMR(const char* dirname){
    string filename = string(dirname) + "/" + kmrfname;
    FILE* fp;

    if((fp = fopen(filename.c_str(), "rb")) == NULL ) {
      printf("cannot open output file %s\n",filename.c_str());
      exit(EXIT_FAILURE);
    }

    int n;
    if(fseek(fp, 0, SEEK_END) == 0) {
      n = ftell(fp)/sizeof(tuple);
      rewind(fp);
      if(n < 0) {
        fprintf(stderr, "%s: Cannot ftell `%s': ", "inputSA", filename.c_str());
        perror(NULL);
        exit(EXIT_FAILURE);
      }
    } else {
      fprintf(stderr, "%s: Cannot fseek `%s': ", "inputSA", filename.c_str());
      perror(NULL);
      exit(EXIT_FAILURE);
    }

    KMR = make_unique<tuple[]>(n);

    if(fread(KMR.get(), sizeof(tuple), (size_t)n, fp) != (size_t)n) {
      fprintf(stderr, "%s: %s `%s': ",
              "inputSA",
              (ferror(fp) || !feof(fp)) ? "Cannot read from" : "Unexpected EOF in",
              filename.c_str());
      perror(NULL);
      exit(EXIT_FAILURE);
    }

    fclose(fp);

    int o[2];
    filename = string(dirname) + "/" + kmrsizefname;
    if((fp = fopen(filename.c_str(), "rb")) == NULL ) {
      printf("cannot open output file %s\n",filename.c_str());
      exit(EXIT_FAILURE);
    }
    if(fread(o, sizeof(int), 2, fp) != 2) {
      fprintf(stderr, "%s: %s `%s': ",
              "inputSA",
              (ferror(fp) || !feof(fp)) ? "Cannot read from" : "Unexpected EOF in",
              filename.c_str());
      perror(NULL);
      exit(EXIT_FAILURE);
    }

    kmerSize = o[0];
    kmerOneFnSize = o[1];

  }

  void inputRKMR(const char* dirname){
    string filename = string(dirname) + "/" + rkmrfname;
    FILE* fp;

    if((fp = fopen(filename.c_str(), "rb")) == NULL ) {
      printf("cannot open output file %s\n",filename.c_str());
      exit(EXIT_FAILURE);
    }

    int n;
    if(fseek(fp, 0, SEEK_END) == 0) {
      n = ftell(fp)/sizeof(tuple);
      rewind(fp);
      if(n < 0) {
        fprintf(stderr, "%s: Cannot ftell `%s': ", "inputSA", filename.c_str());
        perror(NULL);
        exit(EXIT_FAILURE);
      }
    } else {
      fprintf(stderr, "%s: Cannot fseek `%s': ", "inputSA", filename.c_str());
      perror(NULL);
      exit(EXIT_FAILURE);
    }

    RKMR = make_unique<tuple[]>(n);

    if(fread(RKMR.get(), sizeof(tuple), (size_t)n, fp) != (size_t)n) {
      fprintf(stderr, "%s: %s `%s': ",
              "inputSA",
              (ferror(fp) || !feof(fp)) ? "Cannot read from" : "Unexpected EOF in",
              filename.c_str());
      perror(NULL);
      exit(EXIT_FAILURE);
    }

    fclose(fp);

  }

  // int FBWTBackwardSearchCompressed(pair<int,int>* interval,const unsigned char* query,int q,int fn){
  //   interval->first = c->getC(query[q - 1],fn);
  //   interval->second = c->getC(query[q - 1] + 1,fn) - 1;

  //   int _q = q - 1;
  //   int _fn = (fn + 1) % fnum;

  //   if(verbose & LOG_INTERVAL){cout << "query [" << _q  << "] = " << (int)query[_q] << " fn : " << _fn << " [" << interval->first << ":" << interval->second << "]\n";}

  //   while(interval->first <= interval->second && _q > 0){
  //     _q -= 1;
  //     updateIntervalCompressed(interval, query[_q], _fn);
  //     if(verbose & LOG_INTERVAL){cout << "query [" << _q  << "] = " << (int)query[_q] << " fn : " << _fn << " [" << interval->first << ":" << interval->second << "]\n";}
  //     _fn = (_fn + 1) % fnum;
  //   }
  //   return q - _q + (interval->first <= interval->second) - 1;
  // }
public:

  unique_ptr<tuple[]> KMR;//到着点に入っている
  unique_ptr<tuple[]> RKMR;//到着点に入っている
  int kmerSize;
  int kmerOneFnSize;
  unique_ptr<OCC> occ;
  unique_ptr<C> c;
  unique_ptr<OCC> rocc;
  unique_ptr<C> rc;
  unique_ptr<int[]> sa;
  
  unsigned char* ref;

  int fnum;
  int verbose = 0;
  int n;

  int sparseMult;

  class MEMCandidate{
    vector<int> memcandidate;
    vector<pair<int,int>> memcandidate_r; // sa, len
    //    static vector<pair<int,int> > output;
    vector<pair<int,int> > output;
    int memCandidateNum = 0;

  public:

    void printMEMCandidate(){
      cout << "memcandidate : ";
      for(int i = 0;i < memCandidateNum;i++){
        cout << memcandidate[i] << ", \t";
        if(i % 9 == 0) cout << "\n";
      }
      cout << "\n";
    }

    void printMEMCandidate_r(){
      for(int i = 0;i < memCandidateNum;i++){
        cout << "sa : " << memcandidate_r[i].first << ", len : " << memcandidate_r[i].second << "\n";
      }
      cout << "\n";
    }

    void printOutput(){
      //      cout << "output [sa,len] : ";
      int i = 0;
      for(auto& m : output){
        // cout << "[" << m.first << ", " << m.second << "] ,";
        // i++;
        // if(i % 12 == 0) cout << "\n";
        cout << m.first << "  " << m.second << " #mem pos\n";
      }    
      //      cout << "\n";
    }

    void printSize(){
      cout << "output size is " << output.size() << "\n";
    }

    int getSize(){
      return output.size();
    }


    bool countForward(unsigned char* query,int queryLen,int pos,int sa,int len,int memLen,int s,C_OCC& c_occ){
      bool ret = false;
      for(int i = 0;i < c_occ.fnum*s;i++){
        // cout << (int) query[pos + i];

        if(pos + i >= queryLen || sa + len + i >= c_occ.n*c_occ.fnum){
          if(len + i >= memLen) output.push_back(make_pair(sa,len + i));
          return ret;
        }else if(query[pos + i] != c_occ.ref[sa + len + i]){
            // cout << __FUNCTION__ << len + i << "\n";
          // cout << __FUNCTION__ << len + i << " countForward\n";
          if(len + i >= memLen){
            // cout << __FUNCTION__ << len + i << "\n";
            output.push_back(make_pair(sa,len + i));
          }
          return i;
        }
        ret = true;
      }
      // cout << __FUNCTION__ << " " << c_occ.fnum*s << " occ.fnum*s\n";

      return ret;
    }

    bool countForwardAndBackward(unsigned char* query,int queryLen,int pos,int sa,int len,int memLen,int s,C_OCC& c_occ){
      bool ret = false;

      // static int cc = 0;
      // if(cc > 100) exit(1);
      // cc++;

      // backward
      int backwardLen = 0;
      while(true){
        if((pos - len - backwardLen - 1 < 0) || (sa - backwardLen - 1< 0)){
          break;
        }
        if(query[pos - len - backwardLen - 1] == 0){
          break;
        }
        if(query[pos - len - backwardLen - 1] != c_occ.ref[sa - backwardLen - 1]){
          break;
        }else{
          backwardLen++;
        }
      }

      // forward
      for(int i = 0;i < c_occ.fnum*s;i++){
        if(pos + i >= queryLen || sa + len + i >= c_occ.n*c_occ.fnum){
          if(len + i + backwardLen >= memLen) output.push_back(make_pair(sa - backwardLen,len + i + backwardLen));
          return ret;
        }else if(query[pos + i] != c_occ.ref[sa + len + i]){
          if(len + i + backwardLen >= memLen) output.push_back(make_pair(sa - backwardLen,len + i + backwardLen));
          return i;
        }
        ret = true;
      }
      return ret;
    }

    int linerlycountthreashold = 5;

    // 返り値は残った候補数
    int initMEMCandidate(pair<int,int>& interval,int _q,int len,unsigned char* query,int queryLen,int startPos,int memLen,int s,int linerlycountthreashold,C_OCC& c_occ){
      this->linerlycountthreashold = linerlycountthreashold;
      if(interval.second - interval.first + 1 <= linerlycountthreashold){
        for(int i = interval.first;i <= interval.second;i++){
          countForwardAndBackward(query, queryLen, startPos, c_occ.sa[i],len,memLen,s, c_occ);
        }
        memCandidateNum = 0;
        return 0;
      }

      int count = 0;

      if(_q < 0){
        for(int i = interval.first;i <= interval.second;i++){
          countForward(query, queryLen, startPos, c_occ.sa[i],len,memLen,s, c_occ);
        }
        return 0;
      }
      for(int i = interval.first;i <= interval.second;i++){
        if(c_occ.FBWT[i] != query[_q]){
          countForward(query, queryLen, startPos, c_occ.sa[i],len,memLen,s, c_occ);
        }else{
          memcandidate.push_back(c_occ.sa[i] == 0 ? c_occ.n - 1 : c_occ.sa[i] - 1);
          count++;
        }
      }
      memCandidateNum = count;
      return count;
    }

    int initMEMCandidate(pair<int,int>& interval,pair<int,int>& rinterval,int rfn,int _q,int len,unsigned char* query,int queryLen,int startPos,int memLen,int s,int th,C_OCC& c_occ){
      int count = 0;
      pair<int,int> _interval(interval);

      if(_q < 0){
        int i,_rfn;
        for(i = startPos,_rfn = rfn;
	    _interval.first <= _interval.second && i < queryLen && i < startPos + memLen - len && i < startPos + c_occ.fnum*s;
	    i++,_rfn = ((_rfn + 1) % c_occ.fnum)){
	  // cout << _interval.first << ":" << _interval.second << "," << len + i - startPos << "\n";
	  if(query[i] == 4) return 0;
          c_occ.updateForward(_interval, rinterval, query[i], _rfn);
        }

        for(;_interval.first <= _interval.second && i < queryLen && i < startPos + c_occ.fnum*s;
	    i++,_rfn = ((_rfn + 1) % c_occ.fnum)){
          pair<int,int> pre(_interval);

	  if(query[i] == 4){
	    for(int j = pre.first;j <= pre.second;j++){
	      output.push_back(make_pair(c_occ.sa[j],len + i - startPos + 1));
	    }
	    return 0;
	  }

	  c_occ.updateForward(_interval, rinterval, query[i], _rfn);
	  // cout << pre.first << ":" << pre.second << " hogemii pre " << output.size() << "\n";
	  // cout << _interval.first << ":" << _interval.second << " hogemii " << output.size() << "\n";
	  
          if(_interval.first > _interval.second){
	    for(int j = pre.first;j <= pre.second;j++){
	      output.push_back(make_pair(c_occ.sa[j],len + i - startPos + 1));
	    }
            break;
          }
          for(int j = pre.first;j < _interval.first;j++){
            output.push_back(make_pair(c_occ.sa[j],len + i - startPos + 1));
          }
          for(int j = _interval.second + 1;j <= pre.second;j++){
            output.push_back(make_pair(c_occ.sa[j],len + i - startPos + 1));
          }
        }
	
        return 0;
      }

      // int ccc = 0;
      // cout << "hoge : " << ccc++ << "\n";

      memcandidate_r = vector<pair<int,int>>(interval.second - interval.first + 1);

      // cout << "K*s " << c_occ.fnum * s << "\n";

      int iii,_rfn;
      for(iii = startPos,_rfn = rfn;iii < startPos + c_occ.fnum*s ;iii++,_rfn = ((_rfn + 1) % c_occ.fnum)){

        if(query[iii] == 4) {
          for(int j = _interval.first;j <= _interval.second;j++){
            memcandidate_r[j - interval.first] = make_pair(c_occ.sa[j],iii - startPos);
          }
          break;
        }

        if(iii == queryLen){
          for(int j = _interval.first;j <= _interval.second;j++){
            memcandidate_r[j - interval.first] = make_pair(c_occ.sa[j],iii - startPos);
          }
          break;
        }

        pair<int,int> pre(_interval);

        // cout << "hoge\n";

        // cout << __FUNCTION__ << " rfn : " << _rfn << " i : " << iii << " query[iii] : " << (int)query[iii] << " interval [" << interval.first << " : " << interval.second  << "], _interval [" << _interval.first << " : " << _interval.second << "], rinterval [" << rinterval.first << " : " << rinterval.second << "] bef\n";

        c_occ.updateForward(_interval, rinterval, query[iii], _rfn);

        // cout << __FUNCTION__ << " rfn : " << _rfn << " i : " << iii << " query[iii] : " << (int)query[iii] << " interval [" << interval.first << " : " << interval.second  << "], _interval [" << _interval.first << " : " << _interval.second << "], rinterval [" << rinterval.first << " : " << rinterval.second << "] af\n";

        if(_interval.first > _interval.second){
          for(int j = pre.first;j <= pre.second;j++){
	    if(j - interval.first >= interval.second - interval.first + 1){
	      cout << "out of range1! startPos " << startPos << "\n";
	      for(int k = 0;k < 100;k++){
		cout << (int)query[startPos + k];
	      }
	      cout << "\n";
	    }else if(j - interval.first < 0){
	      cout << "out of range2! startPos " << startPos << "\n";
	      for(int k = 0;k < 100;k++){
		cout << (int)query[startPos + k];
	      }
	      cout << "\n";
	    }
            memcandidate_r[j - interval.first] = make_pair(c_occ.sa[j],iii - startPos);
          }
          break;
        }
        for(int j = pre.first;j < _interval.first;j++){
          memcandidate_r[j - interval.first] = make_pair(c_occ.sa[j],iii - startPos);
        }
        if(iii == startPos + c_occ.fnum*s - 1){
          for(int j = _interval.first;j <= _interval.second;j++){
            memcandidate_r[j - interval.first] = make_pair(-1,-1);
          }
        }
        for(int j = _interval.second + 1;j <= pre.second;j++){
          memcandidate_r[j - interval.first] = make_pair(c_occ.sa[j],iii - startPos);
        }
      }

      // if((query[iii] != 4 && iii != queryLen) || (iii == startPos + c_occ.fnum*s)) {
      //   for(int j = _interval.first;j <= _interval.second;j++){
      //     memcandidate_r[j - interval.first] = make_pair(-1,-1);
      //   }
      // }

      // for(int j = 0;j <= interval.second - interval.first;j++){
      //   cout << memcandidate_r[j].first << " : " << memcandidate_r[j].second << "\n";
      // }


      // memCandidateNum = memcandidate_r.size();
      // printMEMCandidate_r();
      //cout << (a) << "," << (int)c_occ.FBWT[(interval.first/OCC::oneEntryNum)*OCC::oneEntryNum + pad] << "\n"; 

#define MEMCAN(a)                                                       \
      if((a) != query[_q]){                                             \
        if(memcandidate_r[i].second + len >= memLen){                   \
          memcandidate_r[i].second = memcandidate_r[i].second + len;    \
          output.push_back(memcandidate_r[i]);                          \
        }                                                               \
      }else{                                                            \
        memcandidate_r[count] = memcandidate_r[i];                      \
        count++;                                                        \
      }


      // initMEMIntervalfLenDist[th*INTERVAL_DIST_MAX_LEN + ((interval.second - interval.first + 1) < INTERVAL_DIST_MAX_LEN - 1 ? (interval.second - interval.first + 1) : INTERVAL_DIST_MAX_LEN - 1)]++;

      int pad = interval.first % OCC::oneEntryNum;
      OCC::oneEntry* e = &c_occ.occ->compressedOcctable[interval.first/OCC::oneEntryNum];
      for(int i = 0;i <= interval.second - interval.first;i++){
        if(pad < 21){
          uint64_t tmp = e->pre[0];
          MEMCAN((tmp >> (pad*3)) & 0x7);
        }else if(pad == 21){
          MEMCAN((e->pre[0] >> 63) | ((e->pre[1] & 0x3) << 1));
        }else if(pad < 42){
          uint64_t tmp = e->pre[1];
          MEMCAN((tmp >> ((pad - 22)*3 + 2)) & 0x7);
        }else if(pad == 42){
          MEMCAN((e->pre[1] >> 62) | ((e->pre[2] & 0x1) << 2));
        }else if(pad < 64){
          uint64_t tmp = e->pre[2];
          MEMCAN((tmp >> ((pad - 43)*3 + 1)) & 0x7);
        }else if(pad < 64 + 21){
          uint64_t tmp = e->post[0];
          MEMCAN((tmp >> ((pad - 64)*3)) & 0x7);
        }else if(pad == 21 + 64){
          MEMCAN((e->post[0] >> 63) | ((e->post[1] & 0x3) << 1));
        }else if(pad < 42 + 64){
          uint64_t tmp = e->post[1];
          MEMCAN((tmp >> ((pad - 22 - 64)*3 + 2)) & 0x7);
        }else if(pad == 42 + 64){
          MEMCAN((e->post[1] >> 62) | ((e->post[2] & 0x1) << 2));
        }else if(pad < 64 + 64){
          uint64_t tmp = e->post[2];
          MEMCAN((tmp >> ((pad - 43 - 64)*3 + 1)) & 0x7);
        }
        pad++;
        if(pad == OCC::oneEntryNum){
          pad = 0;
          e++;
        }
      }

#undef MEMCAN

      // int j = 0;
      // for(int i = interval.first;i <= interval.second;i++){
      //   cout << (int)c_occ.FBWT[i] << "," << (int)query[_q] << "\n";
      //   if(c_occ.FBWT[i] != query[_q]){
      //     if(memcandidate_r[j].second + len >= memLen){
      //       memcandidate_r[j].second = memcandidate_r[j].second + len;
      //       output.push_back(memcandidate_r[j]);
      //     }
      //   }else{
      //     memcandidate_r[count] = memcandidate_r[j];
      //     count++;
      //   }
      //   j++;
      // }
      // cout << count << "\n";
      memCandidateNum = count;
      return count;
    }

    int updateMEMCandidate_r(pair<int,int>& interval,unsigned char q,int len,int fn,int memLen,int th,C_OCC& c_occ){

      int count = 0;

      /*
      for(int i = interval.first,j = 0;i <= interval.second;i++,j++){
        if(c_occ.FBWT[fn*c_occ.n + i] != q){
          if(memcandidate_r[j].second >= 0){
            if(memcandidate_r[j].second + len >= memLen){
              // cout << memcandidate_r[j].first << ", " << memcandidate_r[j].second << " len " << len << " update\n";
              memcandidate_r[j].second += len;
              output.push_back(memcandidate_r[j]);
            }
          }
        }else{
          // if(memcandidate_r[j].first >= 0) isAllMinus = false;
          memcandidate_r[count].first = (memcandidate_r[j].first == 0 ? c_occ.n - 1: (memcandidate_r[j].first < 0 ? -1 : memcandidate_r[j].first - 1));
          memcandidate_r[count].second = memcandidate_r[j].second;
          count++;
        }
      }
      */


#define MEMCAN(a,_iii)                                                      \
      if((a) != q){                                                     \
        if(memcandidate_r[iii + _iii].second >= 0){                            \
          if(memcandidate_r[iii + _iii].second + len >= memLen){               \
            memcandidate_r[iii + _iii].second = memcandidate_r[_iii].second + len; \
            output.push_back(memcandidate_r[_iii]);                      \
          }                                                             \
        }                                                               \
      }else{                                                            \
        memcandidate_r[count].first = (memcandidate_r[iii].first == 0 ? c_occ.n - 1: (memcandidate_r[iii].first < 0 ? -1 : memcandidate_r[iii].first - 1)); \
        memcandidate_r[count].second = memcandidate_r[iii].second;      \
        count++;                                                        \
      }

#define MEMCAN21(_case)                                     \
        MEMCAN(((tmp >> ((pad + _case)*3)) & 0x7),_case);



      // updateMEMIntervalLenDist[th*INTERVAL_DIST_MAX_LEN + ((interval.second - interval.first + 1) < INTERVAL_DIST_MAX_LEN - 1 ? (interval.second - interval.first + 1) : INTERVAL_DIST_MAX_LEN - 1)]++;

      int pad = interval.first % OCC::oneEntryNum;
      OCC::oneEntry* e = &c_occ.occ->compressedOcctable[interval.first/OCC::oneEntryNum + fn*c_occ.occ->compressEntrySize];
      int iii = 0;
      while(true){
        
        if(pad < 21){
          uint64_t tmp = e->pre[0];


          for(;iii <= interval.second - interval.first - 7 && pad < 21 - 7;iii+=8,pad+=8){
            if(((tmp >> (pad*3)) & 0x7) != q){                                                     
              if(memcandidate_r[iii].second >= 0){                              
                if(memcandidate_r[iii].second + len >= memLen){                 
                  memcandidate_r[iii].second = memcandidate_r[iii].second + len; 
                  output.push_back(memcandidate_r[iii]);                        
                }                                                             
              }                                                               
            }else{                                                            
              memcandidate_r[count].first = (memcandidate_r[iii].first == 0 ? c_occ.n - 1: (memcandidate_r[iii].first < 0 ? -1 : memcandidate_r[iii].first - 1)); 
              memcandidate_r[count].second = memcandidate_r[iii].second;        
              count++;                                                        
            }
            if(((tmp >> ((pad + 1)*3)) & 0x7) != q){                                                     
              if(memcandidate_r[iii + 1].second >= 0){                              
                if(memcandidate_r[iii + 1].second + len >= memLen){                 
                  memcandidate_r[iii + 1].second = memcandidate_r[iii + 1].second + len; 
                  output.push_back(memcandidate_r[iii + 1]);                        
                }                                                             
              }                                                               
            }else{                                                            
              memcandidate_r[count].first = (memcandidate_r[iii + 1].first == 0 ? c_occ.n - 1: (memcandidate_r[iii + 1].first < 0 ? -1 : memcandidate_r[iii + 1].first - 1)); 
              memcandidate_r[count].second = memcandidate_r[iii + 1].second;        
              count++;                                                        
            }
            if(((tmp >> ((pad + 2)*3)) & 0x7) != q){                                                     
              if(memcandidate_r[iii + 2].second >= 0){                              
                if(memcandidate_r[iii + 2].second + len >= memLen){                 
                  memcandidate_r[iii + 2].second = memcandidate_r[iii + 2].second + len; 
                  output.push_back(memcandidate_r[iii + 2]);                        
                }                                                             
              }                                                               
            }else{                                                            
              memcandidate_r[count].first = (memcandidate_r[iii + 2].first == 0 ? c_occ.n - 1: (memcandidate_r[iii + 2].first < 0 ? -1 : memcandidate_r[iii + 2].first - 1)); 
              memcandidate_r[count].second = memcandidate_r[iii + 2].second;        
              count++;                                                        
            }
            if(((tmp >> ((pad + 3)*3)) & 0x7) != q){                                                     
              if(memcandidate_r[iii + 3].second >= 0){                              
                if(memcandidate_r[iii + 3].second + len >= memLen){                 
                  memcandidate_r[iii + 3].second = memcandidate_r[iii + 3].second + len; 
                  output.push_back(memcandidate_r[iii]);                        
                }                                                             
              }                                                               
            }else{                                                            
              memcandidate_r[count].first = (memcandidate_r[iii + 3].first == 0 ? c_occ.n - 1: (memcandidate_r[iii + 3].first < 0 ? -1 : memcandidate_r[iii + 3].first - 1)); 
              memcandidate_r[count].second = memcandidate_r[iii + 3].second;        
              count++;                                                        
            }
            if(((tmp >> ((pad + 4)*3)) & 0x7) != q){                                                     
              if(memcandidate_r[iii + 4].second >= 0){                              
                if(memcandidate_r[iii + 4].second + len >= memLen){                 
                  memcandidate_r[iii + 4].second = memcandidate_r[iii + 4].second + len; 
                  output.push_back(memcandidate_r[iii + 4]);                        
                }                                                             
              }                                                               
            }else{                                                            
              memcandidate_r[count].first = (memcandidate_r[iii + 4].first == 0 ? c_occ.n - 1: (memcandidate_r[iii + 4].first < 0 ? -1 : memcandidate_r[iii + 4].first - 1)); 
              memcandidate_r[count].second = memcandidate_r[iii + 4].second;        
              count++;                                                        
            }
            if(((tmp >> ((pad + 5)*3)) & 0x7) != q){                                                     
              if(memcandidate_r[iii + 5].second >= 0){                              
                if(memcandidate_r[iii + 5].second + len >= memLen){                 
                  memcandidate_r[iii + 5].second = memcandidate_r[iii + 5].second + len; 
                  output.push_back(memcandidate_r[iii + 5]);                        
                }                                                             
              }                                                               
            }else{                                                            
              memcandidate_r[count].first = (memcandidate_r[iii + 5].first == 0 ? c_occ.n - 1: (memcandidate_r[iii + 5].first < 0 ? -1 : memcandidate_r[iii + 5].first - 1)); 
              memcandidate_r[count].second = memcandidate_r[iii + 5].second;        
              count++;                                                        
            }
            if(((tmp >> ((pad + 6)*3)) & 0x7) != q){                                                     
              if(memcandidate_r[iii + 6].second >= 0){                              
                if(memcandidate_r[iii + 6].second + len >= memLen){                 
                  memcandidate_r[iii + 6].second = memcandidate_r[iii + 6].second + len; 
                  output.push_back(memcandidate_r[iii + 6]);                        
                }                                                             
              }                                                               
            }else{                                                            
              memcandidate_r[count].first = (memcandidate_r[iii + 6].first == 0 ? c_occ.n - 1: (memcandidate_r[iii + 6].first < 0 ? -1 : memcandidate_r[iii + 6].first - 1)); 
              memcandidate_r[count].second = memcandidate_r[iii + 6].second;        
              count++;                                                        
            }
            if(((tmp >> ((pad + 7)*3)) & 0x7) != q){                                                     
              if(memcandidate_r[iii + 7].second >= 0){                              
                if(memcandidate_r[iii + 7].second + len >= memLen){                 
                  memcandidate_r[iii + 7].second = memcandidate_r[iii + 7].second + len; 
                  output.push_back(memcandidate_r[iii]);                        
                }                                                             
              }                                                               
            }else{                                                            
              memcandidate_r[count].first = (memcandidate_r[iii + 7].first == 0 ? c_occ.n - 1: (memcandidate_r[iii + 7].first < 0 ? -1 : memcandidate_r[iii + 7].first - 1)); 
              memcandidate_r[count].second = memcandidate_r[iii + 7].second;        
              count++;                                                        
            }
          }

          for(;iii <= interval.second - interval.first - 3 && pad < 21 - 3;iii+=4,pad+=4){
            if(((tmp >> (pad*3)) & 0x7) != q){                                                     
              if(memcandidate_r[iii].second >= 0){                              
                if(memcandidate_r[iii].second + len >= memLen){                 
                  memcandidate_r[iii].second = memcandidate_r[iii].second + len; 
                  output.push_back(memcandidate_r[iii]);                        
                }                                                             
              }                                                               
            }else{                                                            
              memcandidate_r[count].first = (memcandidate_r[iii].first == 0 ? c_occ.n - 1: (memcandidate_r[iii].first < 0 ? -1 : memcandidate_r[iii].first - 1)); 
              memcandidate_r[count].second = memcandidate_r[iii].second;        
              count++;                                                        
            }
            if(((tmp >> ((pad + 1)*3)) & 0x7) != q){                                                     
              if(memcandidate_r[iii + 1].second >= 0){                              
                if(memcandidate_r[iii + 1].second + len >= memLen){                 
                  memcandidate_r[iii + 1].second = memcandidate_r[iii + 1].second + len; 
                  output.push_back(memcandidate_r[iii + 1]);                        
                }                                                             
              }                                                               
            }else{                                                            
              memcandidate_r[count].first = (memcandidate_r[iii + 1].first == 0 ? c_occ.n - 1: (memcandidate_r[iii + 1].first < 0 ? -1 : memcandidate_r[iii + 1].first - 1)); 
              memcandidate_r[count].second = memcandidate_r[iii + 1].second;        
              count++;                                                        
            }
            if(((tmp >> ((pad + 2)*3)) & 0x7) != q){                                                     
              if(memcandidate_r[iii + 2].second >= 0){                              
                if(memcandidate_r[iii + 2].second + len >= memLen){                 
                  memcandidate_r[iii + 2].second = memcandidate_r[iii + 2].second + len; 
                  output.push_back(memcandidate_r[iii + 2]);                        
                }                                                             
              }                                                               
            }else{                                                            
              memcandidate_r[count].first = (memcandidate_r[iii + 2].first == 0 ? c_occ.n - 1: (memcandidate_r[iii + 2].first < 0 ? -1 : memcandidate_r[iii + 2].first - 1)); 
              memcandidate_r[count].second = memcandidate_r[iii + 2].second;        
              count++;                                                        
            }
            if(((tmp >> ((pad + 3)*3)) & 0x7) != q){                                                     
              if(memcandidate_r[iii + 3].second >= 0){                              
                if(memcandidate_r[iii + 3].second + len >= memLen){                 
                  memcandidate_r[iii + 3].second = memcandidate_r[iii + 3].second + len; 
                  output.push_back(memcandidate_r[iii]);                        
                }                                                             
              }                                                               
            }else{                                                            
              memcandidate_r[count].first = (memcandidate_r[iii + 3].first == 0 ? c_occ.n - 1: (memcandidate_r[iii + 3].first < 0 ? -1 : memcandidate_r[iii + 3].first - 1)); 
              memcandidate_r[count].second = memcandidate_r[iii + 3].second;        
              count++;                                                        
            }
          }

          for(;iii <= interval.second - interval.first&& pad < 21;iii++,pad++){
            if(((tmp >> (pad*3)) & 0x7) != q){                                                     
              if(memcandidate_r[iii].second >= 0){                              
                if(memcandidate_r[iii].second + len >= memLen){                 
                  memcandidate_r[iii].second = memcandidate_r[iii].second + len; 
                  output.push_back(memcandidate_r[iii]);                        
                }                                                             
              }                                                               
            }else{                                                            
              memcandidate_r[count].first = (memcandidate_r[iii].first == 0 ? c_occ.n - 1: (memcandidate_r[iii].first < 0 ? -1 : memcandidate_r[iii].first - 1)); 
              memcandidate_r[count].second = memcandidate_r[iii].second;        
              count++;                                                        
            }
          }

            //            MEMCAN((tmp >> (pad*3)) & 0x7);

          if(iii == interval.second - interval.first + 1) break;
        }else if(pad == 21){

          if(((e->pre[0] >> 63) | ((e->pre[1] & 0x3) << 1)) != q){                                                     
            if(memcandidate_r[iii].second >= 0){                              
              if(memcandidate_r[iii].second + len >= memLen){                 
                memcandidate_r[iii].second = memcandidate_r[iii].second + len;  
                output.push_back(memcandidate_r[iii]);                        
              }                                                             
            }                                                             
          }else{                                                        
            memcandidate_r[count].first = (memcandidate_r[iii].first == 0 ? c_occ.n - 1: (memcandidate_r[iii].first < 0 ? -1 : memcandidate_r[iii].first - 1)); 
          memcandidate_r[count].second = memcandidate_r[iii].second;        
          count++;                                                        
        }

        //          MEMCAN((e->pre[0] >> 63) | ((e->pre[1] & 0x3) << 1));
        iii++;
        pad++;
        if(iii == interval.second - interval.first + 1) break;
      }else if(pad < 42){
        uint64_t tmp = e->pre[1];
        // - 3 
        for(;iii <= interval.second - interval.first - 7 && pad < 42 - 7;iii+=8,pad+=8){

          if(((tmp >> ((pad - 22)*3 + 2)) & 0x7) != q){                                                     
            if(memcandidate_r[iii].second >= 0){                              
              if(memcandidate_r[iii].second + len >= memLen){                 
                memcandidate_r[iii].second = memcandidate_r[iii].second + len;  
                output.push_back(memcandidate_r[iii]);                        
              }                                                             
            }                                                               
          }else{                                                            
            memcandidate_r[count].first = (memcandidate_r[iii].first == 0 ? c_occ.n - 1: (memcandidate_r[iii].first < 0 ? -1 : memcandidate_r[iii].first - 1)); 
            memcandidate_r[count].second = memcandidate_r[iii].second;        
            count++;                                                        
          }
          if(((tmp >> ((pad + 1 - 22)*3 + 2)) & 0x7) != q){                                                     
            if(memcandidate_r[iii + 1].second >= 0){                              
              if(memcandidate_r[iii + 1].second + len >= memLen){                 
                memcandidate_r[iii + 1].second = memcandidate_r[iii + 1].second + len;  
                output.push_back(memcandidate_r[iii + 1]);                        
              }                                                             
            }                                                               
          }else{                                                            
            memcandidate_r[count].first = (memcandidate_r[iii + 1].first == 0 ? c_occ.n - 1: (memcandidate_r[iii + 1].first < 0 ? -1 : memcandidate_r[iii + 1].first - 1)); 
            memcandidate_r[count].second = memcandidate_r[iii + 1].second;        
            count++;                                                        
          }
          if(((tmp >> ((pad + 2 - 22)*3 + 2)) & 0x7) != q){                                                     
            if(memcandidate_r[iii + 2].second >= 0){                              
              if(memcandidate_r[iii + 2].second + len >= memLen){                 
                memcandidate_r[iii + 2].second = memcandidate_r[iii + 2].second + len;  
                output.push_back(memcandidate_r[iii + 2]);                        
              }                                                             
            }                                                               
          }else{                                                            
            memcandidate_r[count].first = (memcandidate_r[iii + 2].first == 0 ? c_occ.n - 1: (memcandidate_r[iii + 2].first < 0 ? -1 : memcandidate_r[iii + 2].first - 1)); 
            memcandidate_r[count].second = memcandidate_r[iii + 2].second;        
            count++;                                                        
          }
          if(((tmp >> ((pad + 3 - 22)*3 + 2)) & 0x7) != q){                                                     
            if(memcandidate_r[iii + 3].second >= 0){                              
              if(memcandidate_r[iii + 3].second + len >= memLen){                 
                memcandidate_r[iii + 3].second = memcandidate_r[iii + 3].second + len;  
                output.push_back(memcandidate_r[iii + 3]);                        
              }                                                             
            }                                                               
          }else{                                                            
            memcandidate_r[count].first = (memcandidate_r[iii + 3].first == 0 ? c_occ.n - 1: (memcandidate_r[iii + 3].first < 0 ? -1 : memcandidate_r[iii + 3].first - 1)); 
            memcandidate_r[count].second = memcandidate_r[iii + 3].second;        
            count++;                                                        
          }
          if(((tmp >> ((pad + 4 - 22)*3 + 2)) & 0x7) != q){                                                     
            if(memcandidate_r[iii + 4].second >= 0){                              
              if(memcandidate_r[iii + 4].second + len >= memLen){                 
                memcandidate_r[iii + 4].second = memcandidate_r[iii + 4].second + len;  
                output.push_back(memcandidate_r[iii + 4]);                        
              }                                                             
            }                                                               
          }else{                                                            
            memcandidate_r[count].first = (memcandidate_r[iii + 4].first == 0 ? c_occ.n - 1: (memcandidate_r[iii + 4].first < 0 ? -1 : memcandidate_r[iii + 4].first - 1)); 
            memcandidate_r[count].second = memcandidate_r[iii + 4].second;        
            count++;                                                        
          }
          if(((tmp >> ((pad + 5 - 22)*3 + 2)) & 0x7) != q){                                                     
            if(memcandidate_r[iii + 5].second >= 0){                              
              if(memcandidate_r[iii + 5].second + len >= memLen){                 
                memcandidate_r[iii + 5].second = memcandidate_r[iii + 5].second + len;  
                output.push_back(memcandidate_r[iii + 5]);                        
              }                                                             
            }                                                               
          }else{                                                            
            memcandidate_r[count].first = (memcandidate_r[iii + 5].first == 0 ? c_occ.n - 1: (memcandidate_r[iii + 5].first < 0 ? -1 : memcandidate_r[iii + 5].first - 1)); 
            memcandidate_r[count].second = memcandidate_r[iii + 5].second;        
            count++;                                                        
          }
          if(((tmp >> ((pad + 6 - 22)*3 + 2)) & 0x7) != q){                                                     
            if(memcandidate_r[iii + 6].second >= 0){                              
              if(memcandidate_r[iii + 6].second + len >= memLen){                 
                memcandidate_r[iii + 6].second = memcandidate_r[iii + 6].second + len;  
                output.push_back(memcandidate_r[iii + 6]);                        
              }                                                             
            }                                                               
          }else{                                                            
            memcandidate_r[count].first = (memcandidate_r[iii + 6].first == 0 ? c_occ.n - 1: (memcandidate_r[iii + 6].first < 0 ? -1 : memcandidate_r[iii + 6].first - 1)); 
            memcandidate_r[count].second = memcandidate_r[iii + 6].second;        
            count++;                                                        
          }
          if(((tmp >> ((pad + 7 - 22)*3 + 2)) & 0x7) != q){                                                     
            if(memcandidate_r[iii + 7].second >= 0){                              
              if(memcandidate_r[iii + 7].second + len >= memLen){                 
                memcandidate_r[iii + 7].second = memcandidate_r[iii + 7].second + len;  
                output.push_back(memcandidate_r[iii + 7]);                        
              }                                                             
            }                                                               
          }else{                                                            
            memcandidate_r[count].first = (memcandidate_r[iii + 7].first == 0 ? c_occ.n - 1: (memcandidate_r[iii + 7].first < 0 ? -1 : memcandidate_r[iii + 7].first - 1)); 
            memcandidate_r[count].second = memcandidate_r[iii + 7].second;        
            count++;                                                        
          }

          // MEMCAN((tmp >> ((pad - 22)*3 + 2)) & 0x7);
        }
        for(;iii <= interval.second - interval.first - 3 && pad < 42 - 3;iii+=4,pad+=4){

          if(((tmp >> ((pad - 22)*3 + 2)) & 0x7) != q){                                                     
            if(memcandidate_r[iii].second >= 0){                              
              if(memcandidate_r[iii].second + len >= memLen){                 
                memcandidate_r[iii].second = memcandidate_r[iii].second + len;  
                output.push_back(memcandidate_r[iii]);                        
              }                                                             
            }                                                               
          }else{                                                            
            memcandidate_r[count].first = (memcandidate_r[iii].first == 0 ? c_occ.n - 1: (memcandidate_r[iii].first < 0 ? -1 : memcandidate_r[iii].first - 1)); 
            memcandidate_r[count].second = memcandidate_r[iii].second;        
            count++;                                                        
          }
          if(((tmp >> ((pad + 1 - 22)*3 + 2)) & 0x7) != q){                                                     
            if(memcandidate_r[iii + 1].second >= 0){                              
              if(memcandidate_r[iii + 1].second + len >= memLen){                 
                memcandidate_r[iii + 1].second = memcandidate_r[iii + 1].second + len;  
                output.push_back(memcandidate_r[iii + 1]);                        
              }                                                             
            }                                                               
          }else{                                                            
            memcandidate_r[count].first = (memcandidate_r[iii + 1].first == 0 ? c_occ.n - 1: (memcandidate_r[iii + 1].first < 0 ? -1 : memcandidate_r[iii + 1].first - 1)); 
            memcandidate_r[count].second = memcandidate_r[iii + 1].second;        
            count++;                                                        
          }
          if(((tmp >> ((pad + 2 - 22)*3 + 2)) & 0x7) != q){                                                     
            if(memcandidate_r[iii + 2].second >= 0){                              
              if(memcandidate_r[iii + 2].second + len >= memLen){                 
                memcandidate_r[iii + 2].second = memcandidate_r[iii + 2].second + len;  
                output.push_back(memcandidate_r[iii + 2]);                        
              }                                                             
            }                                                               
          }else{                                                            
            memcandidate_r[count].first = (memcandidate_r[iii + 2].first == 0 ? c_occ.n - 1: (memcandidate_r[iii + 2].first < 0 ? -1 : memcandidate_r[iii + 2].first - 1)); 
            memcandidate_r[count].second = memcandidate_r[iii + 2].second;        
            count++;                                                        
          }
          if(((tmp >> ((pad + 3 - 22)*3 + 2)) & 0x7) != q){                                                     
            if(memcandidate_r[iii + 3].second >= 0){                              
              if(memcandidate_r[iii + 3].second + len >= memLen){                 
                memcandidate_r[iii + 3].second = memcandidate_r[iii + 3].second + len;  
                output.push_back(memcandidate_r[iii + 3]);                        
              }                                                             
            }                                                               
          }else{                                                            
            memcandidate_r[count].first = (memcandidate_r[iii + 3].first == 0 ? c_occ.n - 1: (memcandidate_r[iii + 3].first < 0 ? -1 : memcandidate_r[iii + 3].first - 1)); 
            memcandidate_r[count].second = memcandidate_r[iii + 3].second;        
            count++;                                                        
          }

          // MEMCAN((tmp >> ((pad - 22)*3 + 2)) & 0x7);
        }
        for(;iii <= interval.second - interval.first && pad < 42;iii++,pad++){

          if(((tmp >> ((pad - 22)*3 + 2)) & 0x7) != q){                                                     
            if(memcandidate_r[iii].second >= 0){                              
              if(memcandidate_r[iii].second + len >= memLen){                 
                memcandidate_r[iii].second = memcandidate_r[iii].second + len;  
                output.push_back(memcandidate_r[iii]);                        
              }                                                             
            }                                                               
          }else{                                                            
            memcandidate_r[count].first = (memcandidate_r[iii].first == 0 ? c_occ.n - 1: (memcandidate_r[iii].first < 0 ? -1 : memcandidate_r[iii].first - 1)); 
            memcandidate_r[count].second = memcandidate_r[iii].second;        
            count++;                                                        
          }

          // MEMCAN((tmp >> ((pad - 22)*3 + 2)) & 0x7);
        }
        if(iii == interval.second - interval.first + 1) break;
      }else if(pad == 42){

        if(((e->pre[1] >> 62) | ((e->pre[2] & 0x1) << 2)) != q){                                                     
          if(memcandidate_r[iii].second >= 0){                              
            if(memcandidate_r[iii].second + len >= memLen){                 
              memcandidate_r[iii].second = memcandidate_r[iii].second + len;  
              output.push_back(memcandidate_r[iii]);                        
            }                                                             
          }                                                               
        }else{                                                            
          memcandidate_r[count].first = (memcandidate_r[iii].first == 0 ? c_occ.n - 1: (memcandidate_r[iii].first < 0 ? -1 : memcandidate_r[iii].first - 1)); 
          memcandidate_r[count].second = memcandidate_r[iii].second;        
          count++;                                                        
        }

      // MEMCAN((e->pre[1] >> 62) | ((e->pre[2] & 0x1) << 2));

      iii++;
      pad++;
      if(iii == interval.second - interval.first + 1) break;
      }else if(pad < 64){
        uint64_t tmp = e->pre[2];
        for(;iii <= interval.second - interval.first - 7 && pad < 64 - 7;iii+=8,pad+=8){

          if(((tmp >> ((pad - 43)*3 + 1)) & 0x7) != q){                                                     
            if(memcandidate_r[iii].second >= 0){                              
              if(memcandidate_r[iii].second + len >= memLen){                 
                memcandidate_r[iii].second = memcandidate_r[iii].second + len;
                output.push_back(memcandidate_r[iii]);                        
              }                                                             
            }                                                               
          }else{                                                            
            memcandidate_r[count].first = (memcandidate_r[iii].first == 0 ? c_occ.n - 1: (memcandidate_r[iii].first < 0 ? -1 : memcandidate_r[iii].first - 1)); 
            memcandidate_r[count].second = memcandidate_r[iii].second;        
            count++;                                                        
          }
          if(((tmp >> ((pad + 1 - 43)*3 + 1)) & 0x7) != q){                                                     
            if(memcandidate_r[iii + 1].second >= 0){                              
              if(memcandidate_r[iii + 1].second + len >= memLen){                 
                memcandidate_r[iii + 1].second = memcandidate_r[iii + 1].second + len;
                output.push_back(memcandidate_r[iii + 1]);                        
              }                                                             
            }                                                               
          }else{                                                            
            memcandidate_r[count].first = (memcandidate_r[iii + 1].first == 0 ? c_occ.n - 1: (memcandidate_r[iii + 1].first < 0 ? -1 : memcandidate_r[iii + 1].first - 1)); 
            memcandidate_r[count].second = memcandidate_r[iii + 1].second;        
            count++;                                                        
          }
          if(((tmp >> ((pad + 2 - 43)*3 + 1)) & 0x7) != q){                                                     
            if(memcandidate_r[iii + 2].second >= 0){                              
              if(memcandidate_r[iii + 2].second + len >= memLen){                 
                memcandidate_r[iii + 2].second = memcandidate_r[iii + 2].second + len;
                output.push_back(memcandidate_r[iii + 2]);                        
              }                                                             
            }                                                               
          }else{                                                            
            memcandidate_r[count].first = (memcandidate_r[iii + 2].first == 0 ? c_occ.n - 1: (memcandidate_r[iii + 2].first < 0 ? -1 : memcandidate_r[iii + 2].first - 1)); 
            memcandidate_r[count].second = memcandidate_r[iii + 2].second;        
            count++;                                                        
          }
          if(((tmp >> ((pad + 3 - 43)*3 + 1)) & 0x7) != q){                                                     
            if(memcandidate_r[iii + 3].second >= 0){                              
              if(memcandidate_r[iii + 3].second + len >= memLen){                 
                memcandidate_r[iii + 3].second = memcandidate_r[iii + 3].second + len;
                output.push_back(memcandidate_r[iii + 3]);                        
              }                                                             
            }                                                               
          }else{                                                            
            memcandidate_r[count].first = (memcandidate_r[iii + 3].first == 0 ? c_occ.n - 1: (memcandidate_r[iii + 3].first < 0 ? -1 : memcandidate_r[iii + 3].first - 1)); 
            memcandidate_r[count].second = memcandidate_r[iii + 3].second;        
            count++;                                                        
          }
          if(((tmp >> ((pad + 4 - 43)*3 + 1)) & 0x7) != q){                                                     
            if(memcandidate_r[iii + 4].second >= 0){                              
              if(memcandidate_r[iii + 4].second + len >= memLen){                 
                memcandidate_r[iii + 4].second = memcandidate_r[iii + 4].second + len;
                output.push_back(memcandidate_r[iii + 4]);                        
              }                                                             
            }                                                               
          }else{                                                            
            memcandidate_r[count].first = (memcandidate_r[iii + 4].first == 0 ? c_occ.n - 1: (memcandidate_r[iii + 4].first < 0 ? -1 : memcandidate_r[iii + 4].first - 1)); 
            memcandidate_r[count].second = memcandidate_r[iii + 4].second;        
            count++;                                                        
          }
          if(((tmp >> ((pad + 5 - 43)*3 + 1)) & 0x7) != q){                                                     
            if(memcandidate_r[iii + 5].second >= 0){                              
              if(memcandidate_r[iii + 5].second + len >= memLen){                 
                memcandidate_r[iii + 5].second = memcandidate_r[iii + 5].second + len;
                output.push_back(memcandidate_r[iii + 5]);                        
              }                                                             
            }                                                               
          }else{                                                            
            memcandidate_r[count].first = (memcandidate_r[iii + 5].first == 0 ? c_occ.n - 1: (memcandidate_r[iii + 5].first < 0 ? -1 : memcandidate_r[iii + 5].first - 1)); 
            memcandidate_r[count].second = memcandidate_r[iii + 5].second;        
            count++;                                                        
          }
          if(((tmp >> ((pad + 6 - 43)*3 + 1)) & 0x7) != q){                                                     
            if(memcandidate_r[iii + 6].second >= 0){                              
              if(memcandidate_r[iii + 6].second + len >= memLen){                 
                memcandidate_r[iii + 6].second = memcandidate_r[iii + 6].second + len;
                output.push_back(memcandidate_r[iii + 6]);                        
              }                                                             
            }                                                               
          }else{                                                            
            memcandidate_r[count].first = (memcandidate_r[iii + 6].first == 0 ? c_occ.n - 1: (memcandidate_r[iii + 6].first < 0 ? -1 : memcandidate_r[iii + 6].first - 1)); 
            memcandidate_r[count].second = memcandidate_r[iii + 6].second;        
            count++;                                                        
          }
          if(((tmp >> ((pad + 7 - 43)*3 + 1)) & 0x7) != q){                                                     
            if(memcandidate_r[iii + 7].second >= 0){                              
              if(memcandidate_r[iii + 7].second + len >= memLen){                 
                memcandidate_r[iii + 7].second = memcandidate_r[iii + 7].second + len;
                output.push_back(memcandidate_r[iii + 7]);                        
              }                                                             
            }                                                               
          }else{                                                            
            memcandidate_r[count].first = (memcandidate_r[iii + 7].first == 0 ? c_occ.n - 1: (memcandidate_r[iii + 7].first < 0 ? -1 : memcandidate_r[iii + 7].first - 1)); 
            memcandidate_r[count].second = memcandidate_r[iii + 7].second;        
            count++;                                                        
          }

          // MEMCAN((tmp >> ((pad - 43)*3 + 1)) & 0x7);
        }
        for(;iii <= interval.second - interval.first - 3 && pad < 64 - 3;iii+=4,pad+=4){

          if(((tmp >> ((pad - 43)*3 + 1)) & 0x7) != q){                                                     
            if(memcandidate_r[iii].second >= 0){                              
              if(memcandidate_r[iii].second + len >= memLen){                 
                memcandidate_r[iii].second = memcandidate_r[iii].second + len;
                output.push_back(memcandidate_r[iii]);                        
              }                                                             
            }                                                               
          }else{                                                            
            memcandidate_r[count].first = (memcandidate_r[iii].first == 0 ? c_occ.n - 1: (memcandidate_r[iii].first < 0 ? -1 : memcandidate_r[iii].first - 1)); 
            memcandidate_r[count].second = memcandidate_r[iii].second;        
            count++;                                                        
          }
          if(((tmp >> ((pad + 1 - 43)*3 + 1)) & 0x7) != q){                                                     
            if(memcandidate_r[iii + 1].second >= 0){                              
              if(memcandidate_r[iii + 1].second + len >= memLen){                 
                memcandidate_r[iii + 1].second = memcandidate_r[iii + 1].second + len;
                output.push_back(memcandidate_r[iii + 1]);                        
              }                                                             
            }                                                               
          }else{                                                            
            memcandidate_r[count].first = (memcandidate_r[iii + 1].first == 0 ? c_occ.n - 1: (memcandidate_r[iii + 1].first < 0 ? -1 : memcandidate_r[iii + 1].first - 1)); 
            memcandidate_r[count].second = memcandidate_r[iii + 1].second;        
            count++;                                                        
          }
          if(((tmp >> ((pad + 2 - 43)*3 + 1)) & 0x7) != q){                                                     
            if(memcandidate_r[iii + 2].second >= 0){                              
              if(memcandidate_r[iii + 2].second + len >= memLen){                 
                memcandidate_r[iii + 2].second = memcandidate_r[iii + 2].second + len;
                output.push_back(memcandidate_r[iii + 2]);                        
              }                                                             
            }                                                               
          }else{                                                            
            memcandidate_r[count].first = (memcandidate_r[iii + 2].first == 0 ? c_occ.n - 1: (memcandidate_r[iii + 2].first < 0 ? -1 : memcandidate_r[iii + 2].first - 1)); 
            memcandidate_r[count].second = memcandidate_r[iii + 2].second;        
            count++;                                                        
          }
          if(((tmp >> ((pad + 3 - 43)*3 + 1)) & 0x7) != q){                                                     
            if(memcandidate_r[iii + 3].second >= 0){                              
              if(memcandidate_r[iii + 3].second + len >= memLen){                 
                memcandidate_r[iii + 3].second = memcandidate_r[iii + 3].second + len;
                output.push_back(memcandidate_r[iii + 3]);                        
              }                                                             
            }                                                               
          }else{                                                            
            memcandidate_r[count].first = (memcandidate_r[iii + 3].first == 0 ? c_occ.n - 1: (memcandidate_r[iii + 3].first < 0 ? -1 : memcandidate_r[iii + 3].first - 1)); 
            memcandidate_r[count].second = memcandidate_r[iii + 3].second;        
            count++;                                                        
          }

          // MEMCAN((tmp >> ((pad - 43)*3 + 1)) & 0x7);
        }
        for(;iii <= interval.second - interval.first && pad < 64;iii++,pad++){

          if(((tmp >> ((pad - 43)*3 + 1)) & 0x7) != q){                                                     
            if(memcandidate_r[iii].second >= 0){                              
              if(memcandidate_r[iii].second + len >= memLen){                 
                memcandidate_r[iii].second = memcandidate_r[iii].second + len;
                output.push_back(memcandidate_r[iii]);                        
              }                                                             
            }                                                               
          }else{                                                            
            memcandidate_r[count].first = (memcandidate_r[iii].first == 0 ? c_occ.n - 1: (memcandidate_r[iii].first < 0 ? -1 : memcandidate_r[iii].first - 1)); 
            memcandidate_r[count].second = memcandidate_r[iii].second;        
            count++;                                                        
          }

          // MEMCAN((tmp >> ((pad - 43)*3 + 1)) & 0x7);
        }
        if(iii == interval.second - interval.first + 1) break;
      }else if(pad < 64 + 21){
        uint64_t tmp = e->post[0];
        for(;iii <= interval.second - interval.first - 7 && pad < 64 + 21 - 7;iii+=8,pad+=8){
          if(((tmp >> ((pad - 64)*3)) & 0x7) != q){                                                     
            if(memcandidate_r[iii].second >= 0){                              
              if(memcandidate_r[iii].second + len >= memLen){                 
                memcandidate_r[iii].second = memcandidate_r[iii].second + len;
                output.push_back(memcandidate_r[iii]);                        
              }                                                             
            }                                                               
          }else{                                                            
            memcandidate_r[count].first = (memcandidate_r[iii].first == 0 ? c_occ.n - 1: (memcandidate_r[iii].first < 0 ? -1 : memcandidate_r[iii].first - 1));
            memcandidate_r[count].second = memcandidate_r[iii].second;        
            count++;                                                        
          }
          if(((tmp >> ((pad + 1 - 64)*3)) & 0x7) != q){                                                     
            if(memcandidate_r[iii + 1].second >= 0){                              
              if(memcandidate_r[iii + 1].second + len >= memLen){                 
                memcandidate_r[iii + 1].second = memcandidate_r[iii + 1].second + len;
                output.push_back(memcandidate_r[iii + 1]);                        
              }                                                             
            }                                                               
          }else{                                                            
            memcandidate_r[count].first = (memcandidate_r[iii + 1].first == 0 ? c_occ.n - 1: (memcandidate_r[iii + 1].first < 0 ? -1 : memcandidate_r[iii + 1].first - 1));
            memcandidate_r[count].second = memcandidate_r[iii + 1].second;        
            count++;                                                        
          }
          if(((tmp >> ((pad + 2 - 64)*3)) & 0x7) != q){                                                     
            if(memcandidate_r[iii + 2].second >= 0){                              
              if(memcandidate_r[iii + 2].second + len >= memLen){                 
                memcandidate_r[iii + 2].second = memcandidate_r[iii + 2].second + len;
                output.push_back(memcandidate_r[iii + 2]);                        
              }                                                             
            }                                                               
          }else{                                                            
            memcandidate_r[count].first = (memcandidate_r[iii + 2].first == 0 ? c_occ.n - 1: (memcandidate_r[iii + 2].first < 0 ? -1 : memcandidate_r[iii + 2].first - 1));
            memcandidate_r[count].second = memcandidate_r[iii + 2].second;        
            count++;                                                        
          }
          if(((tmp >> ((pad + 3 - 64)*3)) & 0x7) != q){                                                     
            if(memcandidate_r[iii + 3].second >= 0){                              
              if(memcandidate_r[iii + 3].second + len >= memLen){                 
                memcandidate_r[iii + 3].second = memcandidate_r[iii + 3].second + len;
                output.push_back(memcandidate_r[iii + 3]);                        
              }                                                             
            }                                                               
          }else{                                                            
            memcandidate_r[count].first = (memcandidate_r[iii + 3].first == 0 ? c_occ.n - 1: (memcandidate_r[iii + 3].first < 0 ? -1 : memcandidate_r[iii + 3].first - 1));
            memcandidate_r[count].second = memcandidate_r[iii + 3].second;        
            count++;                                                        
          }
          if(((tmp >> ((pad + 4 - 64)*3)) & 0x7) != q){                                                     
            if(memcandidate_r[iii + 4].second >= 0){                              
              if(memcandidate_r[iii + 4].second + len >= memLen){                 
                memcandidate_r[iii + 4].second = memcandidate_r[iii + 4].second + len;
                output.push_back(memcandidate_r[iii + 4]);                        
              }                                                             
            }                                                               
          }else{                                                            
            memcandidate_r[count].first = (memcandidate_r[iii + 4].first == 0 ? c_occ.n - 1: (memcandidate_r[iii + 4].first < 0 ? -1 : memcandidate_r[iii + 4].first - 1));
            memcandidate_r[count].second = memcandidate_r[iii + 4].second;        
            count++;                                                        
          }
          if(((tmp >> ((pad + 5 - 64)*3)) & 0x7) != q){                                                     
            if(memcandidate_r[iii + 5].second >= 0){                              
              if(memcandidate_r[iii + 5].second + len >= memLen){                 
                memcandidate_r[iii + 5].second = memcandidate_r[iii + 5].second + len;
                output.push_back(memcandidate_r[iii + 5]);                        
              }                                                             
            }                                                               
          }else{                                                            
            memcandidate_r[count].first = (memcandidate_r[iii + 5].first == 0 ? c_occ.n - 1: (memcandidate_r[iii + 5].first < 0 ? -1 : memcandidate_r[iii + 5].first - 1));
            memcandidate_r[count].second = memcandidate_r[iii + 5].second;        
            count++;                                                        
          }
          if(((tmp >> ((pad + 6 - 64)*3)) & 0x7) != q){                                                     
            if(memcandidate_r[iii + 6].second >= 0){                              
              if(memcandidate_r[iii + 6].second + len >= memLen){                 
                memcandidate_r[iii + 6].second = memcandidate_r[iii + 6].second + len;
                output.push_back(memcandidate_r[iii + 6]);                        
              }                                                             
            }                                                               
          }else{                                                            
            memcandidate_r[count].first = (memcandidate_r[iii + 6].first == 0 ? c_occ.n - 1: (memcandidate_r[iii + 6].first < 0 ? -1 : memcandidate_r[iii + 6].first - 1));
            memcandidate_r[count].second = memcandidate_r[iii + 6].second;        
            count++;                                                        
          }
          if(((tmp >> ((pad + 7 - 64)*3)) & 0x7) != q){                                                     
            if(memcandidate_r[iii + 7].second >= 0){                              
              if(memcandidate_r[iii + 7].second + len >= memLen){                 
                memcandidate_r[iii + 7].second = memcandidate_r[iii + 7].second + len;
                output.push_back(memcandidate_r[iii + 7]);                        
              }                                                             
            }                                                               
          }else{                                                            
            memcandidate_r[count].first = (memcandidate_r[iii + 7].first == 0 ? c_occ.n - 1: (memcandidate_r[iii + 7].first < 0 ? -1 : memcandidate_r[iii + 7].first - 1));
            memcandidate_r[count].second = memcandidate_r[iii + 7].second;        
            count++;                                                        
          }


          // MEMCAN((tmp >> ((pad - 64)*3)) & 0x7);
        }
        for(;iii <= interval.second - interval.first - 3 && pad < 64 + 21 - 3;iii+=4,pad+=4){
          if(((tmp >> ((pad - 64)*3)) & 0x7) != q){                                                     
            if(memcandidate_r[iii].second >= 0){                              
              if(memcandidate_r[iii].second + len >= memLen){                 
                memcandidate_r[iii].second = memcandidate_r[iii].second + len;
                output.push_back(memcandidate_r[iii]);                        
              }                                                             
            }                                                               
          }else{                                                            
            memcandidate_r[count].first = (memcandidate_r[iii].first == 0 ? c_occ.n - 1: (memcandidate_r[iii].first < 0 ? -1 : memcandidate_r[iii].first - 1));
            memcandidate_r[count].second = memcandidate_r[iii].second;        
            count++;                                                        
          }
          if(((tmp >> ((pad + 1 - 64)*3)) & 0x7) != q){                                                     
            if(memcandidate_r[iii + 1].second >= 0){                              
              if(memcandidate_r[iii + 1].second + len >= memLen){                 
                memcandidate_r[iii + 1].second = memcandidate_r[iii + 1].second + len;
                output.push_back(memcandidate_r[iii + 1]);                        
              }                                                             
            }                                                               
          }else{                                                            
            memcandidate_r[count].first = (memcandidate_r[iii + 1].first == 0 ? c_occ.n - 1: (memcandidate_r[iii + 1].first < 0 ? -1 : memcandidate_r[iii + 1].first - 1));
            memcandidate_r[count].second = memcandidate_r[iii + 1].second;        
            count++;                                                        
          }
          if(((tmp >> ((pad + 2 - 64)*3)) & 0x7) != q){                                                     
            if(memcandidate_r[iii + 2].second >= 0){                              
              if(memcandidate_r[iii + 2].second + len >= memLen){                 
                memcandidate_r[iii + 2].second = memcandidate_r[iii + 2].second + len;
                output.push_back(memcandidate_r[iii + 2]);                        
              }                                                             
            }                                                               
          }else{                                                            
            memcandidate_r[count].first = (memcandidate_r[iii + 2].first == 0 ? c_occ.n - 1: (memcandidate_r[iii + 2].first < 0 ? -1 : memcandidate_r[iii + 2].first - 1));
            memcandidate_r[count].second = memcandidate_r[iii + 2].second;        
            count++;                                                        
          }
          if(((tmp >> ((pad + 3 - 64)*3)) & 0x7) != q){                                                     
            if(memcandidate_r[iii + 3].second >= 0){                              
              if(memcandidate_r[iii + 3].second + len >= memLen){                 
                memcandidate_r[iii + 3].second = memcandidate_r[iii + 3].second + len;
                output.push_back(memcandidate_r[iii + 3]);                        
              }                                                             
            }                                                               
          }else{                                                            
            memcandidate_r[count].first = (memcandidate_r[iii + 3].first == 0 ? c_occ.n - 1: (memcandidate_r[iii + 3].first < 0 ? -1 : memcandidate_r[iii + 3].first - 1));
            memcandidate_r[count].second = memcandidate_r[iii + 3].second;        
            count++;                                                        
          }


          // MEMCAN((tmp >> ((pad - 64)*3)) & 0x7);
        }
        for(;iii <= interval.second - interval.first && pad < 64 + 21;iii++,pad++){


          if(((tmp >> ((pad - 64)*3)) & 0x7) != q){                                                     
            if(memcandidate_r[iii].second >= 0){                              
              if(memcandidate_r[iii].second + len >= memLen){                 
                memcandidate_r[iii].second = memcandidate_r[iii].second + len;
                output.push_back(memcandidate_r[iii]);                        
              }                                                             
            }                                                               
          }else{                                                            
            memcandidate_r[count].first = (memcandidate_r[iii].first == 0 ? c_occ.n - 1: (memcandidate_r[iii].first < 0 ? -1 : memcandidate_r[iii].first - 1));
            memcandidate_r[count].second = memcandidate_r[iii].second;        
            count++;                                                        
          }


          // MEMCAN((tmp >> ((pad - 64)*3)) & 0x7);
        }
        if(iii == interval.second - interval.first + 1) break;
      }else if(pad == 21 + 64){

        if(((e->post[0] >> 63) | ((e->post[1] & 0x3) << 1)) != q){                                                     
          if(memcandidate_r[iii].second >= 0){                              
            if(memcandidate_r[iii].second + len >= memLen){                 
              memcandidate_r[iii].second = memcandidate_r[iii].second + len;
              output.push_back(memcandidate_r[iii]);                       
            }                                                             
          }                                                               
        }else{                                                            
          memcandidate_r[count].first = (memcandidate_r[iii].first == 0 ? c_occ.n - 1: (memcandidate_r[iii].first < 0 ? -1 : memcandidate_r[iii].first - 1)); 
          memcandidate_r[count].second = memcandidate_r[iii].second;        
          count++;                                                        
        }

        // MEMCAN((e->post[0] >> 63) | ((e->post[1] & 0x3) << 1));


        iii++;
        pad++;
        if(iii == interval.second - interval.first + 1) break;
      }else if(pad < 42 + 64){
        uint64_t tmp = e->post[1];
        for(;iii <= interval.second - interval.first - 7 && pad < 64 + 42 - 7;iii+=8,pad+=8){
          if(((tmp >> ((pad - 22 - 64)*3 + 2)) & 0x7) != q){              
            if(memcandidate_r[iii].second >= 0){                          
              if(memcandidate_r[iii].second + len >= memLen){             
                memcandidate_r[iii].second = memcandidate_r[iii].second + len;
                output.push_back(memcandidate_r[iii]);                        
              }                                                          
            }                                                            
          }else{                                                         
            memcandidate_r[count].first = (memcandidate_r[iii].first == 0 ? c_occ.n - 1: (memcandidate_r[iii].first < 0 ? -1 : memcandidate_r[iii].first - 1)); 
            memcandidate_r[count].second = memcandidate_r[iii].second;    
            count++;                                                      
          }
          if(((tmp >> ((pad + 1 - 22 - 64)*3 + 2)) & 0x7) != q){              
            if(memcandidate_r[iii + 1].second >= 0){                          
              if(memcandidate_r[iii + 1].second + len >= memLen){             
                memcandidate_r[iii + 1].second = memcandidate_r[iii + 1].second + len;
                output.push_back(memcandidate_r[iii + 1]);                        
              }                                                          
            }                                                            
          }else{                                                         
            memcandidate_r[count].first = (memcandidate_r[iii + 1].first == 0 ? c_occ.n - 1: (memcandidate_r[iii + 1].first < 0 ? -1 : memcandidate_r[iii + 1].first - 1)); 
            memcandidate_r[count].second = memcandidate_r[iii + 1].second;    
            count++;                                                      
          }
          if(((tmp >> ((pad + 2 - 22 - 64)*3 + 2)) & 0x7) != q){              
            if(memcandidate_r[iii + 2].second >= 0){                          
              if(memcandidate_r[iii + 2].second + len >= memLen){             
                memcandidate_r[iii + 2].second = memcandidate_r[iii + 2].second + len;
                output.push_back(memcandidate_r[iii + 2]);                        
              }                                                          
            }                                                            
          }else{                                                         
            memcandidate_r[count].first = (memcandidate_r[iii + 2].first == 0 ? c_occ.n - 1: (memcandidate_r[iii + 2].first < 0 ? -1 : memcandidate_r[iii + 2].first - 1)); 
            memcandidate_r[count].second = memcandidate_r[iii + 2].second;    
            count++;                                                      
          }
          if(((tmp >> ((pad + 3 - 22 - 64)*3 + 2)) & 0x7) != q){              
            if(memcandidate_r[iii + 3].second >= 0){                          
              if(memcandidate_r[iii + 3].second + len >= memLen){             
                memcandidate_r[iii + 3].second = memcandidate_r[iii + 3].second + len;
                output.push_back(memcandidate_r[iii + 3]);                        
              }                                                          
            }                                                            
          }else{                                                         
            memcandidate_r[count].first = (memcandidate_r[iii + 3].first == 0 ? c_occ.n - 1: (memcandidate_r[iii + 3].first < 0 ? -1 : memcandidate_r[iii + 3].first - 1)); 
            memcandidate_r[count].second = memcandidate_r[iii + 3].second;    
            count++;                                                      
          }
          if(((tmp >> ((pad + 4 - 22 - 64)*3 + 2)) & 0x7) != q){              
            if(memcandidate_r[iii + 4].second >= 0){                          
              if(memcandidate_r[iii + 4].second + len >= memLen){             
                memcandidate_r[iii + 4].second = memcandidate_r[iii + 4].second + len;
                output.push_back(memcandidate_r[iii + 4]);                        
              }                                                          
            }                                                            
          }else{                                                         
            memcandidate_r[count].first = (memcandidate_r[iii + 4].first == 0 ? c_occ.n - 1: (memcandidate_r[iii + 4].first < 0 ? -1 : memcandidate_r[iii + 4].first - 1)); 
            memcandidate_r[count].second = memcandidate_r[iii + 4].second;    
            count++;                                                      
          }
          if(((tmp >> ((pad + 5 - 22 - 64)*3 + 2)) & 0x7) != q){              
            if(memcandidate_r[iii + 5].second >= 0){                          
              if(memcandidate_r[iii + 5].second + len >= memLen){             
                memcandidate_r[iii + 5].second = memcandidate_r[iii + 5].second + len;
                output.push_back(memcandidate_r[iii + 5]);                        
              }                                                          
            }                                                            
          }else{                                                         
            memcandidate_r[count].first = (memcandidate_r[iii + 5].first == 0 ? c_occ.n - 1: (memcandidate_r[iii + 5].first < 0 ? -1 : memcandidate_r[iii + 5].first - 1)); 
            memcandidate_r[count].second = memcandidate_r[iii + 5].second;    
            count++;                                                      
          }
          if(((tmp >> ((pad + 6 - 22 - 64)*3 + 2)) & 0x7) != q){              
            if(memcandidate_r[iii + 6].second >= 0){                          
              if(memcandidate_r[iii + 6].second + len >= memLen){             
                memcandidate_r[iii + 6].second = memcandidate_r[iii + 6].second + len;
                output.push_back(memcandidate_r[iii + 6]);                        
              }                                                          
            }                                                            
          }else{                                                         
            memcandidate_r[count].first = (memcandidate_r[iii + 6].first == 0 ? c_occ.n - 1: (memcandidate_r[iii + 6].first < 0 ? -1 : memcandidate_r[iii + 6].first - 1)); 
            memcandidate_r[count].second = memcandidate_r[iii + 6].second;    
            count++;                                                      
          }
          if(((tmp >> ((pad + 7 - 22 - 64)*3 + 2)) & 0x7) != q){              
            if(memcandidate_r[iii + 7].second >= 0){                          
              if(memcandidate_r[iii + 7].second + len >= memLen){             
                memcandidate_r[iii + 7].second = memcandidate_r[iii + 7].second + len;
                output.push_back(memcandidate_r[iii + 7]);                        
              }                                                          
            }                                                            
          }else{                                                         
            memcandidate_r[count].first = (memcandidate_r[iii + 7].first == 0 ? c_occ.n - 1: (memcandidate_r[iii + 7].first < 0 ? -1 : memcandidate_r[iii + 7].first - 1)); 
            memcandidate_r[count].second = memcandidate_r[iii + 7].second;    
            count++;                                                      
          }
          // MEMCAN((tmp >> ((pad - 22 - 64)*3 + 2)) & 0x7);
        }
        for(;iii <= interval.second - interval.first - 3 && pad < 64 + 42 - 3;iii+=4,pad+=4){
          if(((tmp >> ((pad - 22 - 64)*3 + 2)) & 0x7) != q){              
            if(memcandidate_r[iii].second >= 0){                          
              if(memcandidate_r[iii].second + len >= memLen){             
                memcandidate_r[iii].second = memcandidate_r[iii].second + len;
                output.push_back(memcandidate_r[iii]);                        
              }                                                          
            }                                                            
          }else{                                                         
            memcandidate_r[count].first = (memcandidate_r[iii].first == 0 ? c_occ.n - 1: (memcandidate_r[iii].first < 0 ? -1 : memcandidate_r[iii].first - 1)); 
            memcandidate_r[count].second = memcandidate_r[iii].second;    
            count++;                                                      
          }
          if(((tmp >> ((pad + 1 - 22 - 64)*3 + 2)) & 0x7) != q){              
            if(memcandidate_r[iii + 1].second >= 0){                          
              if(memcandidate_r[iii + 1].second + len >= memLen){             
                memcandidate_r[iii + 1].second = memcandidate_r[iii + 1].second + len;
                output.push_back(memcandidate_r[iii + 1]);                        
              }                                                          
            }                                                            
          }else{                                                         
            memcandidate_r[count].first = (memcandidate_r[iii + 1].first == 0 ? c_occ.n - 1: (memcandidate_r[iii + 1].first < 0 ? -1 : memcandidate_r[iii + 1].first - 1)); 
            memcandidate_r[count].second = memcandidate_r[iii + 1].second;    
            count++;                                                      
          }
          if(((tmp >> ((pad + 2 - 22 - 64)*3 + 2)) & 0x7) != q){              
            if(memcandidate_r[iii + 2].second >= 0){                          
              if(memcandidate_r[iii + 2].second + len >= memLen){             
                memcandidate_r[iii + 2].second = memcandidate_r[iii + 2].second + len;
                output.push_back(memcandidate_r[iii + 2]);                        
              }                                                          
            }                                                            
          }else{                                                         
            memcandidate_r[count].first = (memcandidate_r[iii + 2].first == 0 ? c_occ.n - 1: (memcandidate_r[iii + 2].first < 0 ? -1 : memcandidate_r[iii + 2].first - 1)); 
            memcandidate_r[count].second = memcandidate_r[iii + 2].second;    
            count++;                                                      
          }
          if(((tmp >> ((pad + 3 - 22 - 64)*3 + 2)) & 0x7) != q){              
            if(memcandidate_r[iii + 3].second >= 0){                          
              if(memcandidate_r[iii + 3].second + len >= memLen){             
                memcandidate_r[iii + 3].second = memcandidate_r[iii + 3].second + len;
                output.push_back(memcandidate_r[iii + 3]);                        
              }                                                          
            }                                                            
          }else{                                                         
            memcandidate_r[count].first = (memcandidate_r[iii + 3].first == 0 ? c_occ.n - 1: (memcandidate_r[iii + 3].first < 0 ? -1 : memcandidate_r[iii + 3].first - 1)); 
            memcandidate_r[count].second = memcandidate_r[iii + 3].second;    
            count++;                                                      
          }
          // MEMCAN((tmp >> ((pad - 22 - 64)*3 + 2)) & 0x7);
        }
        for(;iii <= interval.second - interval.first && pad < 64 + 42;iii++,pad++){
          if(((tmp >> ((pad - 22 - 64)*3 + 2)) & 0x7) != q){              
            if(memcandidate_r[iii].second >= 0){                          
              if(memcandidate_r[iii].second + len >= memLen){             
                memcandidate_r[iii].second = memcandidate_r[iii].second + len;
                output.push_back(memcandidate_r[iii]);                        
              }                                                          
            }                                                            
          }else{                                                         
            memcandidate_r[count].first = (memcandidate_r[iii].first == 0 ? c_occ.n - 1: (memcandidate_r[iii].first < 0 ? -1 : memcandidate_r[iii].first - 1)); 
            memcandidate_r[count].second = memcandidate_r[iii].second;    
            count++;                                                      
          }
          // MEMCAN((tmp >> ((pad - 22 - 64)*3 + 2)) & 0x7);
        }
        if(iii == interval.second - interval.first + 1) break;
      }else if(pad == 42 + 64){
        if(((e->post[1] >> 62) | ((e->post[2] & 0x1) << 2)) != q){        
          if(memcandidate_r[iii].second >= 0){                            
            if(memcandidate_r[iii].second + len >= memLen){               
              memcandidate_r[iii].second = memcandidate_r[iii].second + len;
              output.push_back(memcandidate_r[iii]);                  
            }                                                         
          }                                                           
        }else{                                                        
          memcandidate_r[count].first = (memcandidate_r[iii].first == 0 ? c_occ.n - 1: (memcandidate_r[iii].first < 0 ? -1 : memcandidate_r[iii].first - 1)); 
          memcandidate_r[count].second = memcandidate_r[iii].second;   
          count++;                                                     
        }
        // MEMCAN((e->post[1] >> 62) | ((e->post[2] & 0x1) << 2));
        iii++;
        pad++;
        if(iii == interval.second - interval.first + 1) break;
        }else if(pad < 64 + 64){
          uint64_t tmp = e->post[2];
          for(;iii <= interval.second - interval.first - 7 && pad < 64 + 64 - 7;iii+=8,pad+=8){
            if(((tmp >> ((pad - 43 - 64)*3 + 1)) & 0x7) != q){           
              if(memcandidate_r[iii].second >= 0){                       
                if(memcandidate_r[iii].second + len >= memLen){          
                  memcandidate_r[iii].second = memcandidate_r[iii].second + len; 
                  output.push_back(memcandidate_r[iii]);                  
                }                                                         
              }                                                           
            }else{                                                        
              memcandidate_r[count].first = (memcandidate_r[iii].first == 0 ? c_occ.n - 1: (memcandidate_r[iii].first < 0 ? -1 : memcandidate_r[iii].first - 1)); 
              memcandidate_r[count].second = memcandidate_r[iii].second;   
              count++;                                                     
            }
            if(((tmp >> ((pad + 1 - 43 - 64)*3 + 1)) & 0x7) != q){           
              if(memcandidate_r[iii + 1].second >= 0){                       
                if(memcandidate_r[iii + 1].second + len >= memLen){          
                  memcandidate_r[iii + 1].second = memcandidate_r[iii + 1].second + len; 
                  output.push_back(memcandidate_r[iii + 1]);                  
                }                                                         
              }                                                           
            }else{                                                        
              memcandidate_r[count].first = (memcandidate_r[iii + 1].first == 0 ? c_occ.n - 1: (memcandidate_r[iii + 1].first < 0 ? -1 : memcandidate_r[iii + 1].first - 1)); 
              memcandidate_r[count].second = memcandidate_r[iii + 1].second;   
              count++;                                                     
            }
            if(((tmp >> ((pad + 2 - 43 - 64)*3 + 1)) & 0x7) != q){           
              if(memcandidate_r[iii + 2].second >= 0){                       
                if(memcandidate_r[iii + 2].second + len >= memLen){          
                  memcandidate_r[iii + 2].second = memcandidate_r[iii + 2].second + len; 
                  output.push_back(memcandidate_r[iii + 2]);                  
                }                                                         
              }                                                           
            }else{                                                        
              memcandidate_r[count].first = (memcandidate_r[iii + 2].first == 0 ? c_occ.n - 1: (memcandidate_r[iii + 2].first < 0 ? -1 : memcandidate_r[iii + 2].first - 1)); 
              memcandidate_r[count].second = memcandidate_r[iii + 2].second;   
              count++;                                                     
            }
            if(((tmp >> ((pad + 3 - 43 - 64)*3 + 1)) & 0x7) != q){           
              if(memcandidate_r[iii + 3].second >= 0){                       
                if(memcandidate_r[iii + 3].second + len >= memLen){          
                  memcandidate_r[iii + 3].second = memcandidate_r[iii + 3].second + len; 
                  output.push_back(memcandidate_r[iii + 3]);                  
                }                                                         
              }                                                           
            }else{                                                        
              memcandidate_r[count].first = (memcandidate_r[iii + 3].first == 0 ? c_occ.n - 1: (memcandidate_r[iii + 3].first < 0 ? -1 : memcandidate_r[iii + 3].first - 1)); 
              memcandidate_r[count].second = memcandidate_r[iii + 3].second;   
              count++;                                                     
            }
            if(((tmp >> ((pad + 4 - 43 - 64)*3 + 1)) & 0x7) != q){           
              if(memcandidate_r[iii + 4].second >= 0){                       
                if(memcandidate_r[iii + 4].second + len >= memLen){          
                  memcandidate_r[iii + 4].second = memcandidate_r[iii + 4].second + len; 
                  output.push_back(memcandidate_r[iii + 4]);                  
                }                                                         
              }                                                           
            }else{                                                        
              memcandidate_r[count].first = (memcandidate_r[iii + 4].first == 0 ? c_occ.n - 1: (memcandidate_r[iii + 4].first < 0 ? -1 : memcandidate_r[iii + 4].first - 1)); 
              memcandidate_r[count].second = memcandidate_r[iii + 4].second;   
              count++;                                                     
            }
            if(((tmp >> ((pad + 5 - 43 - 64)*3 + 1)) & 0x7) != q){           
              if(memcandidate_r[iii + 5].second >= 0){                       
                if(memcandidate_r[iii + 5].second + len >= memLen){          
                  memcandidate_r[iii + 5].second = memcandidate_r[iii + 5].second + len; 
                  output.push_back(memcandidate_r[iii + 5]);                  
                }                                                         
              }                                                           
            }else{                                                        
              memcandidate_r[count].first = (memcandidate_r[iii + 5].first == 0 ? c_occ.n - 1: (memcandidate_r[iii + 5].first < 0 ? -1 : memcandidate_r[iii + 5].first - 1)); 
              memcandidate_r[count].second = memcandidate_r[iii + 5].second;   
              count++;                                                     
            }
            if(((tmp >> ((pad + 6 - 43 - 64)*3 + 1)) & 0x7) != q){           
              if(memcandidate_r[iii + 6].second >= 0){                       
                if(memcandidate_r[iii + 6].second + len >= memLen){          
                  memcandidate_r[iii + 6].second = memcandidate_r[iii + 6].second + len; 
                  output.push_back(memcandidate_r[iii + 6]);                  
                }                                                         
              }                                                           
            }else{                                                        
              memcandidate_r[count].first = (memcandidate_r[iii + 6].first == 0 ? c_occ.n - 1: (memcandidate_r[iii + 6].first < 0 ? -1 : memcandidate_r[iii + 6].first - 1)); 
              memcandidate_r[count].second = memcandidate_r[iii + 6].second;   
              count++;                                                     
            }
            if(((tmp >> ((pad + 7 - 43 - 64)*3 + 1)) & 0x7) != q){           
              if(memcandidate_r[iii + 7].second >= 0){                       
                if(memcandidate_r[iii + 7].second + len >= memLen){          
                  memcandidate_r[iii + 7].second = memcandidate_r[iii + 7].second + len; 
                  output.push_back(memcandidate_r[iii + 7]);                  
                }                                                         
              }                                                           
            }else{                                                        
              memcandidate_r[count].first = (memcandidate_r[iii + 7].first == 0 ? c_occ.n - 1: (memcandidate_r[iii + 7].first < 0 ? -1 : memcandidate_r[iii + 7].first - 1)); 
              memcandidate_r[count].second = memcandidate_r[iii + 7].second;   
              count++;
            }
            // MEMCAN((tmp >> ((pad - 43 - 64)*3 + 1)) & 0x7);
          }
          for(;iii <= interval.second - interval.first - 3 && pad < 64 + 64 - 3;iii+=4,pad+=4){
            if(((tmp >> ((pad - 43 - 64)*3 + 1)) & 0x7) != q){           
              if(memcandidate_r[iii].second >= 0){                       
                if(memcandidate_r[iii].second + len >= memLen){          
                  memcandidate_r[iii].second = memcandidate_r[iii].second + len; 
                  output.push_back(memcandidate_r[iii]);                  
                }                                                         
              }                                                           
            }else{                                                        
              memcandidate_r[count].first = (memcandidate_r[iii].first == 0 ? c_occ.n - 1: (memcandidate_r[iii].first < 0 ? -1 : memcandidate_r[iii].first - 1)); 
              memcandidate_r[count].second = memcandidate_r[iii].second;   
              count++;                                                     
            }
            if(((tmp >> ((pad + 1 - 43 - 64)*3 + 1)) & 0x7) != q){           
              if(memcandidate_r[iii + 1].second >= 0){                       
                if(memcandidate_r[iii + 1].second + len >= memLen){          
                  memcandidate_r[iii + 1].second = memcandidate_r[iii + 1].second + len; 
                  output.push_back(memcandidate_r[iii + 1]);                  
                }                                                         
              }                                                           
            }else{                                                        
              memcandidate_r[count].first = (memcandidate_r[iii + 1].first == 0 ? c_occ.n - 1: (memcandidate_r[iii + 1].first < 0 ? -1 : memcandidate_r[iii + 1].first - 1)); 
              memcandidate_r[count].second = memcandidate_r[iii + 1].second;   
              count++;                                                     
            }
            if(((tmp >> ((pad + 2 - 43 - 64)*3 + 1)) & 0x7) != q){           
              if(memcandidate_r[iii + 2].second >= 0){                       
                if(memcandidate_r[iii + 2].second + len >= memLen){          
                  memcandidate_r[iii + 2].second = memcandidate_r[iii + 2].second + len; 
                  output.push_back(memcandidate_r[iii + 2]);                  
                }                                                         
              }                                                           
            }else{                                                        
              memcandidate_r[count].first = (memcandidate_r[iii + 2].first == 0 ? c_occ.n - 1: (memcandidate_r[iii + 2].first < 0 ? -1 : memcandidate_r[iii + 2].first - 1)); 
              memcandidate_r[count].second = memcandidate_r[iii + 2].second;   
              count++;                                                     
            }
            if(((tmp >> ((pad + 3 - 43 - 64)*3 + 1)) & 0x7) != q){           
              if(memcandidate_r[iii + 3].second >= 0){                       
                if(memcandidate_r[iii + 3].second + len >= memLen){          
                  memcandidate_r[iii + 3].second = memcandidate_r[iii + 3].second + len; 
                  output.push_back(memcandidate_r[iii + 3]);                  
                }                                                         
              }                                                           
            }else{                                                        
              memcandidate_r[count].first = (memcandidate_r[iii + 3].first == 0 ? c_occ.n - 1: (memcandidate_r[iii + 3].first < 0 ? -1 : memcandidate_r[iii + 3].first - 1)); 
              memcandidate_r[count].second = memcandidate_r[iii + 3].second;   
              count++;                                                     
            }
            // MEMCAN((tmp >> ((pad - 43 - 64)*3 + 1)) & 0x7);
          }
          for(;iii <= interval.second - interval.first && pad < 64 + 64;iii++,pad++){
            if(((tmp >> ((pad - 43 - 64)*3 + 1)) & 0x7) != q){           
              if(memcandidate_r[iii].second >= 0){                       
                if(memcandidate_r[iii].second + len >= memLen){          
                  memcandidate_r[iii].second = memcandidate_r[iii].second + len; 
                  output.push_back(memcandidate_r[iii]);                  
                }                                                         
              }                                                           
            }else{                                                        
              memcandidate_r[count].first = (memcandidate_r[iii].first == 0 ? c_occ.n - 1: (memcandidate_r[iii].first < 0 ? -1 : memcandidate_r[iii].first - 1)); 
              memcandidate_r[count].second = memcandidate_r[iii].second;   
              count++;                                                     
            }
            // MEMCAN((tmp >> ((pad - 43 - 64)*3 + 1)) & 0x7);
          }
          if(iii == interval.second - interval.first + 1) break;
          if(pad == 64 + 64){
            pad = 0;
            e++;
          }
        }
        // pad++;
        // if(pad == OCC::oneEntryNum){
        //   pad = 0;
        //   e++;
        // }
      }

#undef MEMCAN


      // bool isAllMinus = true;

      // for(int i = interval.first,j = 0;i <= interval.second;i++,j++){
      //   if(c_occ.FBWT[fn*c_occ.n + i] != q){
      //     if(memcandidate_r[j].second >= 0){
      //       if(memcandidate_r[j].second + len >= memLen){
      //         // cout << memcandidate_r[j].first << ", " << memcandidate_r[j].second << " len " << len << " update\n";
      //         memcandidate_r[j].second += len;
      //         output.push_back(memcandidate_r[j]);
      //       }
      //     }
      //   }else{
      //     // if(memcandidate_r[j].first >= 0) isAllMinus = false;
      //     memcandidate_r[count].first = (memcandidate_r[j].first == 0 ? c_occ.n - 1: (memcandidate_r[j].first < 0 ? -1 : memcandidate_r[j].first - 1));
      //     memcandidate_r[count].second = memcandidate_r[j].second;
      //     count++;
      //   }
      // }

      // if(isAllMinus){
      //   memCandidateNum = 0;
      //   return 0;
      // }else{
      //   memCandidateNum = count;
      //   return count;
      // }
      memCandidateNum = count;
      return count;
    }

    int updateMEMCandidate(pair<int,int>& interval,char q,int len,int fn,unsigned char* query,int queryLen,int startPos,int memLen,int s,C_OCC& c_occ){
      if(interval.second - interval.first + 1 <= linerlycountthreashold){
        for(int i = 0;i <= interval.second - interval.first;i++){
          countForwardAndBackward(query, queryLen, startPos, memcandidate[i],len,memLen,s, c_occ);
        }
        memCandidateNum = 0;
        return 0;
      }

      int count = 0;
      for(int i = interval.first,j = 0;i <= interval.second;i++,j++){
        if(c_occ.FBWT[fn*c_occ.n + i] != q){
          // cout << memcandidate[j] << " len " << len <<" update\n";
          countForward(query, queryLen, startPos, memcandidate[j],len,memLen,s, c_occ);
        }else{
          memcandidate[count] = (memcandidate[j] == 0 ? c_occ.n - 1: memcandidate[j] - 1);
          count++;
        }
      }
      memCandidateNum = count;
      return count;
    }

    void finalizeMEMCandidate(int len,unsigned char* query,int queryLen,int startPos,int memLen,int s,C_OCC& c_occ){
      for(int i = 0;i < memCandidateNum;i++){
        // cout << memcandidate[i] << " final\n";
        countForward(query, queryLen, startPos, memcandidate[i],len,memLen,s, c_occ);
      }
      memCandidateNum = 0;
    }

    void finalizeMEMCandidate_r(int len,int memLen){
      for(int i = 0;i < memCandidateNum;i++){
        if(memcandidate_r[i].second >= 0){
          if(memcandidate_r[i].second + len >= memLen){
            // cout << memcandidate_r[i].first << " : " << memcandidate_r[i].second << " final\n";
            memcandidate_r[i].second += len;            
            output.push_back(memcandidate_r[i]);
          }
        }
      }
      // for(int i = 0;i < memCandidateNum;i++){
      //   countForward(query, queryLen, startPos, memcandidate[i],len,memLen,s, c_occ);
      // }
      memCandidateNum = 0;
    }
  };

  void construct(unique_ptr< unsigned char[]>& FBWT,unique_ptr< unsigned char[]>& RFBWT,unique_ptr<int[]>& SA,int n,int fn,unsigned char* ref,int kmersize,int sparseMult){
    occ = make_unique<OCC>(OCC(FBWT,n,fn,*this));
    c = make_unique<C>(C(FBWT,n,fn));
    rc = make_unique<C>(C(RFBWT,n,fn)); // this is for outputIndex. not for search
    rocc = make_unique<OCC>(OCC(RFBWT,n,fn,*this));
    fnum = fn;
    this->kmerSize = kmersize;
    this->FBWT = move(FBWT);
    this->RFBWT = move(RFBWT);
    this->sa = move(SA);
    this->n = n;
    this->ref = ref;
    this->sparseMult = sparseMult;
    constructKMR(*this, kmersize);
    constructRKMR(*this, kmersize);
  }

  void construct(unique_ptr< unsigned char[]>& FBWT,unique_ptr< unsigned char[]>& RFBWT,unique_ptr<int[]>& SA,int n,int fn,unsigned char* ref,int kmersize,int sparseMult,int verbose){
    construct(FBWT,RFBWT,SA,n,fn,ref,kmersize,sparseMult);
    this->verbose = verbose;
  }

  static const int LOG_INTERVAL = 1;
  static const int LOG_MEM = 2;

  const string safname = "sa";
  const string reffname = "ref";
  const string kmrfname = "kmr";
  const string rkmrfname = "rkmr";
  const string kmrsizefname = "kmrsize";
  const string occfname = "occ";
  const string roccfname = "rocc";
  const string cfname = "c";
  const string rcfname = "rc";


  static const string fbwtfname;
  static const string rfbwtfname;

  void outputFBWT(const char* dirname){
    FILE* fp;
    string filename = string(dirname) + "/" + fbwtfname;
    if((fp = fopen(filename.c_str(), "wb")) == NULL ) {
      printf("cannot open output file %s\n",filename.c_str());
      exit(EXIT_FAILURE);
    }
    fwrite(FBWT.get() , sizeof(unsigned char), n*fnum, fp);
    fclose(fp);
  }

  void outputRFBWT(const char* dirname){
    FILE* fp;
    string filename = string(dirname) + "/" + rfbwtfname;
    if((fp = fopen(filename.c_str(), "wb")) == NULL ) {
      printf("cannot open output file %s\n",filename.c_str());
      exit(EXIT_FAILURE);
    }
    fwrite(RFBWT.get() , sizeof(unsigned char), n*fnum, fp);
    fclose(fp);
  }

  int inputFBWT(const char* dirname){
    string filename = string(dirname) + "/" + fbwtfname;
    FILE* fp;

    if((fp = fopen(filename.c_str(), "rb")) == NULL ) {
      printf("cannot open output file %s\n",filename.c_str());
      exit(EXIT_FAILURE);
    }

    int n;
    if(fseek(fp, 0, SEEK_END) == 0) {
      n = ftell(fp);
      rewind(fp);
      if(n < 0) {
        fprintf(stderr, "%s: Cannot ftell `%s': ", "inputFBWT", filename.c_str());
        perror(NULL);
        exit(EXIT_FAILURE);
      }
    } else {
      fprintf(stderr, "%s: Cannot fseek `%s': ", "inputFBWT", filename.c_str());
      perror(NULL);
      exit(EXIT_FAILURE);
    }

    FBWT = make_unique<unsigned char[]>(n);      

    if(fread(FBWT.get(), sizeof(unsigned char), (size_t)n, fp) != (size_t)n) {
      fprintf(stderr, "%s: %s `%s': ",
              "inputFBWT",
              (ferror(fp) || !feof(fp)) ? "Cannot read from" : "Unexpected EOF in",
              filename.c_str());
      perror(NULL);
      exit(EXIT_FAILURE);
    }

    fclose(fp);

    return n;
  }

  int inputRFBWT(const char* dirname){
    string filename = string(dirname) + "/" + rfbwtfname;
    FILE* fp;

    if((fp = fopen(filename.c_str(), "rb")) == NULL ) {
      printf("cannot open output file %s\n",filename.c_str());
      exit(EXIT_FAILURE);
    }

    int n;
    if(fseek(fp, 0, SEEK_END) == 0) {
      n = ftell(fp);
      rewind(fp);
      if(n < 0) {
        fprintf(stderr, "%s: Cannot ftell `%s': ", "inputFBWT", filename.c_str());
        perror(NULL);
        exit(EXIT_FAILURE);
      }
    } else {
      fprintf(stderr, "%s: Cannot fseek `%s': ", "inputFBWT", filename.c_str());
      perror(NULL);
      exit(EXIT_FAILURE);
    }

    RFBWT = make_unique<unsigned char[]>(n);      

    if(fread(RFBWT.get(), sizeof(unsigned char), (size_t)n, fp) != (size_t)n) {
      fprintf(stderr, "%s: %s `%s': ",
              "inputFBWT",
              (ferror(fp) || !feof(fp)) ? "Cannot read from" : "Unexpected EOF in",
              filename.c_str());
      perror(NULL);
      exit(EXIT_FAILURE);
    }

    fclose(fp);

    return n;
  }

  void outputSA(const char* dirname){
    FILE* fp;
    string filename = string(dirname) + "/" + safname;
    if((fp = fopen(filename.c_str(), "wb")) == NULL ) {
      printf("cannot open output file %s\n",filename.c_str());
      exit(EXIT_FAILURE);
    }
    fwrite(sa.get() , sizeof(int), n, fp);
    fclose(fp);
  }

  void outputRef(const char* dirname){
    FILE* fp;
    string filename = string(dirname) + "/" + reffname;
    if((fp = fopen(filename.c_str(), "wb")) == NULL ) {
      printf("cannot open output file %s\n",filename.c_str());
      exit(EXIT_FAILURE);
    }
    fwrite(ref , sizeof(char), n*fnum, fp);
    fclose(fp);
  }

  void inputSA(const char* dirname){
    string filename = string(dirname) + "/" + safname;
    FILE* fp;

    if((fp = fopen(filename.c_str(), "rb")) == NULL ) {
      printf("cannot open output file %s\n",filename.c_str());
      exit(EXIT_FAILURE);
    }

    int n;
    if(fseek(fp, 0, SEEK_END) == 0) {
      n = ftell(fp)/sizeof(int);
      rewind(fp);
      if(n < 0) {
        fprintf(stderr, "%s: Cannot ftell `%s': ", "inputSA", filename.c_str());
        perror(NULL);
        exit(EXIT_FAILURE);
      }
    } else {
      fprintf(stderr, "%s: Cannot fseek `%s': ", "inputSA", filename.c_str());
      perror(NULL);
      exit(EXIT_FAILURE);
    }

    sa = make_unique<int[]>(n);      

    if(fread(sa.get(), sizeof(int), (size_t)n, fp) != (size_t)n) {
      fprintf(stderr, "%s: %s `%s': ",
              "inputSA",
              (ferror(fp) || !feof(fp)) ? "Cannot read from" : "Unexpected EOF in",
              filename.c_str());
      perror(NULL);
      exit(EXIT_FAILURE);
    }

    fclose(fp);
  }

  void inputRef(const char* dirname){
    string filename = string(dirname) + "/" + reffname;
    FILE* fp;

    if((fp = fopen(filename.c_str(), "rb")) == NULL ) {
      printf("cannot open output file %s\n",filename.c_str());
      exit(EXIT_FAILURE);
    }

    int n;
    if(fseek(fp, 0, SEEK_END) == 0) {
      n = ftell(fp)/sizeof(char);
      rewind(fp);
      if(n < 0) {
        fprintf(stderr, "%s: Cannot ftell `%s': ", "inputSA", filename.c_str());
        perror(NULL);
        exit(EXIT_FAILURE);
      }
    } else {
      fprintf(stderr, "%s: Cannot fseek `%s': ", "inputSA", filename.c_str());
      perror(NULL);
      exit(EXIT_FAILURE);
    }

    ref = (unsigned char*)malloc(sizeof(char)*n);

    if(fread(ref, sizeof(unsigned char), (size_t)n, fp) != (size_t)n) {
      fprintf(stderr, "%s: %s `%s': ",
              "inputSA",
              (ferror(fp) || !feof(fp)) ? "Cannot read from" : "Unexpected EOF in",
              filename.c_str());
      perror(NULL);
      exit(EXIT_FAILURE);
    }

    fclose(fp);
  }

public:
  // n is one element of FBWT
  C_OCC(unique_ptr< unsigned char[]>& FBWT,unique_ptr< unsigned char[]>& RFBWT,unique_ptr<int[]>& SA,int n,int fn,int kmersize,unsigned char* ref,int sparseMult){
    construct(FBWT,RFBWT,SA,n,fn,ref,kmersize,sparseMult);
  }

  C_OCC(unique_ptr< unsigned char[]>& FBWT,unique_ptr< unsigned char[]>& RFBWT,unique_ptr<int[]>& SA,int n,int fn,unsigned char* ref,int kmersize,int sparseMult,int verbose){
    construct(FBWT,RFBWT,SA,n,fn,ref,kmersize,sparseMult,verbose);
  }

  C_OCC(const char* dirname,int sparseMult,int verbose){
    c = make_unique<C>(C(dirname,cfname));
    rc = make_unique<C>(C(dirname,rcfname));
    occ = make_unique<OCC>(OCC(dirname,c->getFnum(),*this,occfname));
    rocc = make_unique<OCC>(OCC(dirname,c->getFnum(),*this,roccfname));
    n = occ->getN();
    inputSA(dirname);
    // inputFBWT(dirname);
    // inputRFBWT(dirname);
    // inputRef(dirname);
    inputKMR(dirname);
    inputRKMR(dirname);
    this->verbose = verbose;
    this->fnum = c->getFnum();
    this->sparseMult = sparseMult;
  }

  void outputFile(const char* dirname){
    c->outputFile(dirname,cfname);
    rc->outputFile(dirname,rcfname);
    occ->outputFile(dirname,occfname);
    rocc->outputFile(dirname,roccfname);
    outputSA(dirname);
    outputFBWT(dirname);
    // outputRFBWT(dirname);
    outputRef(dirname);
    outputKMR(dirname);
    outputRKMR(dirname);
  }

  void setSparseMult(int s){
    sparseMult = s;
  }

  static int getFBWTSize(const char* dirname){
    string filename = string((char*)dirname) + "/" + fbwtfname;
    FILE* fp;
    if((fp = fopen(filename.c_str(), "rb")) == NULL ) {
      printf("cannot open output file %s\n",filename.c_str());
      exit(EXIT_FAILURE);
    }
    int n;
    if(fseek(fp, 0, SEEK_END) == 0) {
      n = ftell(fp);
      rewind(fp);
      if(n < 0) {
        fprintf(stderr, "%s: Cannot ftell `%s': ", "C", filename.c_str());
        perror(NULL);
        exit(EXIT_FAILURE);
      }
    } else {
      fprintf(stderr, "%s: Cannot fseek `%s': ", "C", filename.c_str());
      perror(NULL);
      exit(EXIT_FAILURE);
    }
    fclose(fp);

    return n;
  }


  // this is not good coding
  void updateRIntervalCompressed(pair<int,int>& interval,char q,int fn){
    interval.first = rc->getC(q,fn) + rocc->getOccFromCompressed(interval.first,q,fn);
    interval.second = rc->getC(q,fn) + rocc->getOccFromCompressed(interval.second + 1,q,fn) - 1;
  }
  // fnは逆順のもの
  void updateForward(pair<int,int>& interval,pair<int,int>& rinterval,unsigned char q,int fn){

    // this is pretty bad
    int lessthanq = rinterval.second + 1 - rinterval.first;
    // int lessthanq = 0;
    // for(int i = 1;i < q;i++){
    //   lessthanq += (rocc->getOccFromCompressed(rinterval.second + 1,i,fn) - rocc->getOccFromCompressed(rinterval.first,i,fn));
    // }

    int rinterval_first_occ = rocc->getOccFromCompressed(rinterval.first,q,fn);
    int rinterval_second_occ = rocc->getOccFromCompressed(rinterval.second + 1,q,fn);

    lessthanq -= (rinterval_second_occ - rinterval_first_occ);

    if(q + 1 < 4)
      lessthanq -= rocc->getOccMorethanCFromCompressed(rinterval.second + 1,q + 1,fn) - rocc->getOccMorethanCFromCompressed(rinterval.first,q + 1,fn);
    
    // cout << "fore lessthanq after\n";

    // for(int i = q + 1;i < 4;i++){
    //   lessthanq -= ((rocc->getOccFromCompressed(rinterval.second + 1,i,fn) - rocc->getOccFromCompressed(rinterval.first,i,fn)));
    // }

    //    cout << "fn : " << fn << ", q : " << (int)q << ", lessthanq : " << lessthanq << "\n";

    
    //    cout << "hoge\n";

    //    cout << "rinterval_first_occ : " << rinterval_first_occ << ", rinterval_second_occ : " << rinterval_second_occ << "\n";

    rinterval.first = rc->getC(q,fn) + rinterval_first_occ;
    rinterval.second = rc->getC(q,fn) + rinterval_second_occ - 1;

    interval.first = interval.first + lessthanq;
    interval.second = interval.first + rinterval_second_occ - rinterval_first_occ - 1;
  }

  // fnは正順のもの
  void updateBackward(pair<int,int>& interval,pair<int,int>& rinterval,unsigned char q,int fn){

    // this is not good
    int lessthanq = interval.second + 1 - interval.first;
    int interval_first_occ = occ->getOccFromCompressed(interval.first,q,fn);
    int interval_second_occ = occ->getOccFromCompressed(interval.second + 1,q,fn);

    // for(int i = 1;i < q;i++){
    //   lessthanq += (occ->getOccFromCompressed(interval.second + 1,i,fn) - occ->getOccFromCompressed(interval.first,i,fn));
    // }

    lessthanq -= (interval_second_occ - interval_first_occ);
    //getOccMorethanCFromCompressed(int i,unsigned char c,int fn)

    // for(int i = q + 1;i < 4;i++){
    //   lessthanq -= (occ->getOccFromCompressed(interval.second + 1,i,fn) - occ->getOccFromCompressed(interval.first,i,fn));
    // }
    // cout << "back lessthanq before\n";
    if(q + 1 < 4)
      lessthanq -= occ->getOccMorethanCFromCompressed(interval.second + 1,q + 1,fn) - occ->getOccMorethanCFromCompressed(interval.first,q + 1,fn);
    // }
    // cout << "back lessthanq after\n";

    interval.first = c->getC(q,fn) + interval_first_occ;
    interval.second = c->getC(q,fn) + interval_second_occ - 1;

    rinterval.first = rinterval.first + lessthanq;
    rinterval.second = rinterval.first + interval_second_occ - interval_first_occ - 1;
  }

  void testBackwardForward(){

    pair<int,int> interval(15,25);
    pair<int,int> rinterval(66,76);

    updateForward(interval, rinterval, 4, 2);

    cout << "interval [" << interval.first << " : " << interval.second << "], rinterval [" << rinterval.first << " : " << rinterval.second << "]\n";

    interval.first = 15; interval.second = 25;
    rinterval.first = 66; rinterval.second = 76;
    updateBackward(interval, rinterval, 2, 0);

    cout << "interval [" << interval.first << " : " << interval.second << "], rinterval [" << rinterval.first << " : " << rinterval.second << "]\n";
    
  }

  void testInitMem(){

    // printOCCCompressed(1);

    MEMCandidate mc;

    // int key = 0x0;
    // cout << "hash interval [" << KMR[key + 0*kmerOneFnSize].left << " : " << KMR[key + 0*kmerOneFnSize].right << "]  [" << RKMR[key + 0*kmerOneFnSize].left << " : " << RKMR[key + 0*kmerOneFnSize].right << "]\n";
    // cout << "hash interval [" << KMR[key + 1*kmerOneFnSize].left << " : " << KMR[key + 1*kmerOneFnSize].right << "]  [" << RKMR[key + 1*kmerOneFnSize].left << " : " << RKMR[key + 1*kmerOneFnSize].right << "]\n";
    // cout << "hash interval [" << KMR[key + 2*kmerOneFnSize].left << " : " << KMR[key + 2*kmerOneFnSize].right << "]  [" << RKMR[key + 2*kmerOneFnSize].left << " : " << RKMR[key + 2*kmerOneFnSize].right << "]\n";
    // cout << "hash interval [" << KMR[key + 3*kmerOneFnSize].left << " : " << KMR[key + 3*kmerOneFnSize].right << "]  [" << RKMR[key + 3*kmerOneFnSize].left << " : " << RKMR[key + 3*kmerOneFnSize].right << "]\n";

    
    pair<int,int> interval(15,25);
    pair<int,int> rinterval(66,76);

    //    string query = "ACAGGCAGA";
    unsigned char query[9] = {0,1,0,2,2,1,0,2,0};

    char vv[5] = {'A','C','G','T','N'};
    for(int i = 0;i < 9;i++){
      cout << vv[query[i]];
    }
    cout << "\n";

    //    int initMEMCandidate(pair<int,int>& interval,pair<int,int>& rinterval,int rfn,int _q,int len,const char* query,int startPos,int memLen,int s,C_OCC& c_occ){

    int num = mc.initMEMCandidate(interval, rinterval, 2, 1, 2, query,9, 4, 7, 1,0, *this);
    
    mc.printMEMCandidate_r();

    //  void updateIntervalCompressed(pair<int,int>& interval,char q,int fn,int cannum){
    updateIntervalCompressed(interval,query[1],0,num);
    cout << "updated interval [" << interval.first << " : " << interval.second << "]\n";

    mc.updateMEMCandidate_r(interval,query[0],3,1,7,0,*this);
    mc.printMEMCandidate_r();

    mc.finalizeMEMCandidate_r(3, 7);

    //updateMEMCandidate_r(pair<int,int>& interval,char q,int len,int fn,int memLen,C_OCC& c_occ){
    //    updateMEMCandidate_r(interval,);
    

  }

  void testFindMEM_r(){
    //  int FindMEMCompressed_r(unsigned char* query,int q,int startPos,int firstMatchnum,const int memLen,const int s,const int intervalMaxSize,int* memcount){

    //    string query = "ACAGGCAGA";
    unsigned char query[9] = {1,2,1,3,3,2,1,3,1};

    char vv[5] = {'A','C','G','T','N'};
    for(int i = 0;i < 9;i++){
      cout << vv[query[i]];
    }
    cout << "\n";

    uint64_t c = 0;
    FindMEMCompressed_r(query,9,5,3,5,1,100,&c,0);
    
  }




  void updateIntervalCompressed(pair<int,int>& interval,char q,int fn){
    interval.first = c->getC(q,fn) + occ->getOccFromCompressed(interval.first,q,fn);
    interval.second = c->getC(q,fn) + occ->getOccFromCompressed(interval.second + 1,q,fn) - 1;
  }

  void updateIntervalCompressed(pair<int,int>& interval,char q,int fn,int cannum){
    // cout << "interval.first " << interval.first << " q " << (int)q << " fn " << fn << " getC : " << c->getC(q,fn) << " getOcc : " << occ->getOccFromCompressed(interval.first,q,fn) << "\n";
    interval.first = c->getC(q,fn) + occ->getOccFromCompressed(interval.first,q,fn);
    interval.second = interval.first + cannum - 1;
  }


  void testOcc(){
    occ->printCompressedOcctable();
    occ->testGetOccFromCompressed();
  }

  void printOCCCompressed(int c){
    occ->printOCCCompressed(c);
  }

  void printOCCCompressed(){
    occ->printOCCCompressed();
  }

  void printMorethanC(){
    occ->printMorethanC();
  }
  void printMorethanC(int c){
    occ->printMorethanC(c);
  }

  void printC(){
    c->printC();
  }


  int FBWTBackwardSearchCompressed(pair<int,int>& interval,const unsigned char* query,int q,int fn){
    interval.first = c->getC(query[q - 1],fn);
    interval.second = c->getC(query[q - 1] + 1,fn) - 1;

    int _q = q - 1;
    int _fn = (fn + 1) % fnum;

    while(interval.first <= interval.second && _q > 0 && query[_q] != 0){
      _q -= 1;
      updateIntervalCompressed(interval, query[_q], _fn);
      _fn = (_fn + 1) % fnum;
    }
    return q - _q + (interval.first <= interval.second) - 1;
  }

  void FindAllMEMCompressed(unsigned char* query,int q,int startPos,int memLen,int linearcomp,const int intervalMaxSize,int* count){
    // int cou = 0;

    // cout << "fnum : " << fnum << ", sparseMult : " << sparseMult << "\n";
    // cout << __FUNCTION__ << " start\n";

    // FindMEMCompressed(query, q, 161809754,memLen - fnum*sparseMult + 1,memLen,sparseMult,linearcomp,intervalMaxSize,count);
    // exit(1);

    // cout << intervalMaxSize << "\n";
    for(int k = 0;k < fnum;k++){
      for(int i = startPos - k;i >= memLen - fnum*sparseMult + 1;i-=fnum*sparseMult){
        // cou++;
        //        if(cou % 1000 == 0 && cou >= 422000) 
        // if(cou < 423000 && cou >= 422000) 
	// cout << cou << "/" << i << "/" << startPos <<  " : " << *count << "\n";
        FindMEMCompressed(query, q, i,memLen - fnum*sparseMult + 1,memLen,sparseMult,linearcomp,intervalMaxSize,count);
        //        if(cou == 100) exit(1);
      }
    }
  }

  void FindAllMEMCompressed_r(unsigned char* query,int q,int startPos,int memLen,const int intervalMaxSize,uint64_t* count,int th){
    // cout << __FUNCTION__ << " start\n";

    if(sparseMult <= 0) {
      cerr << "sparseMult must be more than 0\n";
      exit(1);
    }
    // FindMEMCompressed_r(query, q, 524159,memLen - fnum*sparseMult + 1,memLen,sparseMult,intervalMaxSize,count,th);
    // exit(1);

    // int cou = 0;
    // cout << "fnum : " << fnum << ", sparseMult : " << sparseMult << "\n";
    // cout << intervalMaxSize << "\n";
    for(int k = 0;k < fnum;k++){
      //      cout << k << " times\n";
      for(int i = startPos - k;i >= memLen - fnum*sparseMult + 1;i-=fnum*sparseMult){
      // for(int i = 156997437;i >= memLen - fnum*sparseMult + 1;i-=fnum*sparseMult){
      //   cou++;
	// if(cou % 100 == 0)
	  // cout << cou << "/" << i << "/" << startPos <<  " : " << *count << "\n";
        FindMEMCompressed_r(query, q, i,memLen - fnum*sparseMult + 1,memLen,sparseMult,intervalMaxSize,count,th);
	long m = show_getrusage();
	if(maxMemorySize < m){
	  maxMemorySize = m;
	}
        //        if(cou == 100) exit(1);
      }
    }
  }

  int FindMEMCompressed(unsigned char* query,int queryLen,int startPos,int firstMatchnum,const int memLen,const int s,const int linearcomp,const int intervalMaxSize,int* memcount){

    pair<int,int> interval;
    if(startPos >= queryLen) return 0;

    int key = 0;
    for(int i = 0;i < kmerSize;i++){
      if(query[startPos - 1 - i] == 4) return 0;
      key |= ((query[startPos - 1 - i]) << 2*i);
    }

    int _q = startPos - 1 - kmerSize;
    int _fn = (fnum - (firstMatchnum % fnum) + kmerSize) % fnum;//(firstFn + 1) % fnum;

    interval.first = KMR[key + _fn*kmerOneFnSize].left;
    interval.second = KMR[key + _fn*kmerOneFnSize].right;
    //165691706 
    //    cout << startPos << "\n";
    // if(startPos == 165691738){
      // cout << "c_fn : " << _fn << " [" << interval.first << ":" << interval.second << "]\n";
    // }

    if(interval.first > interval.second){
      return 0;
    }

    for(int i = firstMatchnum - 1 - kmerSize;
        interval.first <= interval.second && _q >= 0 && i >= 0;
        _q--,i--,_fn = ((_fn + 1) % fnum)){
      if(query[_q] == 4) return 0;
      updateIntervalCompressed(interval, query[_q], _fn);
      // if(startPos == 165691738){
        // cout << "b_fn : " << _fn << " _q : " << _q << " [" << interval.first << ":" << interval.second << "]\n";
      // }
    }

    if(interval.first > interval.second || _q != startPos - 1 - firstMatchnum || _fn != 0) {
      return startPos - _q + (interval.first <= interval.second) - 1;
    }

    int susunda = 0;
    bool isBreak = false;
    while((interval.second - interval.first + 1 > intervalMaxSize) && (susunda < s - 1)){
      susunda++;
      for(int i = fnum - 1;
          interval.first <= interval.second && _q >= 0 && i >= 0;
          _q--,i--,_fn = ((_fn + 1) % fnum)){
        if(query[_q] == 4) {isBreak = true;break;}
        updateIntervalCompressed(interval, query[_q], _fn);
        // if(startPos == 165691738){
          // cout << "d_fn : " << _fn << " _q : " << _q << " [" << interval.first << ":" << interval.second << "]\n";
        // }
      }

      if(interval.first > interval.second || _fn != 0) {
        isBreak = true;
      }

      if(isBreak) break;
    }

    if(susunda > 0){
      FindMEMCompressed(query,queryLen,startPos + fnum*(s - susunda),firstMatchnum + (s - susunda)*fnum,memLen,susunda,linearcomp,intervalMaxSize,memcount);
    }

    if(interval.first > interval.second || _fn != 0 || isBreak) {
      return startPos - _q + (interval.first <= interval.second) - 1;
    }

    MEMCandidate mc;
    int cannum = mc.initMEMCandidate(interval, _q, firstMatchnum + susunda*fnum,query,queryLen,startPos,memLen,s - susunda,linearcomp,*this);
    // cout << "cannum : " << cannum << "\n";
    int count = susunda*fnum;
    while(1){
      if(cannum <= 0) break;
      if(interval.first > interval.second || _q < 0) break;
      updateIntervalCompressed(interval, query[_q], _fn,cannum);

      // if(startPos == 165691738){
        // cout << "a_fn : " << _fn << "_q : " << _q << " [" << interval.first << ":" << interval.second << "]\n";
      // }

      _q -= 1;
      _fn = (_fn + 1) % fnum;

      count++;
      if(interval.first > interval.second || _q < 0) break;
      if(query[_q] == 4) break;
      cannum = mc.updateMEMCandidate(interval, query[_q], firstMatchnum + count,_fn,query,queryLen,startPos,memLen,s - susunda,*this);
    }

    mc.finalizeMEMCandidate(firstMatchnum + count,query,queryLen,startPos,memLen,s - susunda,*this);
    //mc.printOutput();

    // mc.printSize();
    *memcount += mc.getSize();

    // static int cccc = 0;
    // if(cccc == 100){
    //   exit(0);
    // }
    // cccc++;
    // if(startPos == 165691738){
    //   exit(1);
    // }

    //    cout << startPos << " : " << *memcount << "\n";
    //    countMEMS = mc.getSize();

    return startPos - _q + (interval.first <= interval.second) - 1;
  }

  int FindMEMCompressed_r(unsigned char* query,int queryLen,int startPos,int firstMatchnum,const int memLen,const int s,const int intervalMaxSize,uint64_t* memcount,int th){

    pair<int,int> interval;
    pair<int,int> rinterval;

    if(startPos >= queryLen) return 0;
    int key = 0;
    for(int i = 0;i < kmerSize;i++){
      if(query[startPos - 1 - i] == 4) return 0;
      key |= ((query[startPos - 1 - i]) << 2*i);
    }

    int _q = startPos - 1 - kmerSize;
    int _fn = (fnum - (firstMatchnum % fnum) + kmerSize) % fnum;//(firstFn + 1) % fnum;
    int _rfn = ((fnum - _fn + kmerSize - 1) % fnum + 1)%fnum;//(fnum - _fn + 1 + (kmerSize - 1)%fnum) % fnum;

    interval.first = KMR[key + _fn*kmerOneFnSize].left;
    interval.second = KMR[key + _fn*kmerOneFnSize].right;

    rinterval.first = RKMR[key + _rfn*kmerOneFnSize].left;
    rinterval.second = RKMR[key + _rfn*kmerOneFnSize].right;

    // cout << startPos << "\n";

    // if(startPos == 165691738){
      // cout << "c_fn : " << _fn << " [" << interval.first << ":" << interval.second << "]\n";
    // }


    if(interval.first > interval.second){
      return 0;
    }

    for(int i = firstMatchnum - 1 - kmerSize;
        interval.first <= interval.second && _q >= 0 && i >= 0;
        _q--,i--,_fn = ((_fn + 1) % fnum)){
      if(query[_q] == 4) return 0;

      // cout << __FUNCTION__ << " a fn : " << _fn << ",interval : [" << interval.first << ":" << interval.second << "] rfn : " << _rfn  << ",rinterval : [" << rinterval.first << ":" << rinterval.second << "]\n" ;

      updateBackward(interval,rinterval, query[_q], _fn);
      // if(startPos == 165691738){
        // cout << "b_fn : " << _fn << " _q : " << _q << " [" << interval.first << ":" << interval.second << "]\n";
      // }
    }

    // cout << __FUNCTION__ << " b fn : " << _fn << ",interval : [" << interval.first << ":" << interval.second << "] rfn : " << _rfn  << ",rinterval : [" << rinterval.first << ":" << rinterval.second << "]\n" ;

    if(interval.first > interval.second || _q != startPos - 1 - firstMatchnum || _fn != 0) {
      return startPos - _q + (interval.first <= interval.second) - 1;
    }

    int susunda = 0;
    bool isBreak = false;
    while((interval.second - interval.first + 1 > intervalMaxSize) && (susunda < s - 1)){
      susunda++;
      for(int i = fnum - 1;
          interval.first <= interval.second && _q >= 0 && i >= 0;
          _q--,i--,_fn = ((_fn + 1) % fnum)){
        if(query[_q] == 4) {isBreak = true;break;}
        // cout << __FUNCTION__ << " c1 fn : " << _fn << ",interval : [" << interval.first << ":" << interval.second << "] rfn : " << _rfn  << ",rinterval : [" << rinterval.first << ":" << rinterval.second << "]\n" ;
        updateBackward(interval, rinterval, query[_q], _fn);
        // if(startPos == 165691738){
          // cout << "d_fn : " << _fn << " _q : " << _q << " [" << interval.first << ":" << interval.second << "]\n";
        // }
        // cout << __FUNCTION__ << " c2 fn : " << _fn << ",interval : [" << interval.first << ":" << interval.second << "] rfn : " << _rfn  << ",rinterval : [" << rinterval.first << ":" << rinterval.second << "]\n" ;
      }

      if(interval.first > interval.second || _fn != 0) {
        isBreak = true;
      }

      if(isBreak) break;
    }

    if(susunda > 0){
      FindMEMCompressed_r(query,queryLen,startPos + fnum*(s - susunda),firstMatchnum + (s - susunda)*fnum,memLen,susunda,intervalMaxSize,memcount,th);
    }

    if(interval.first > interval.second || _fn != 0 || isBreak) {
      return startPos - _q + (interval.first <= interval.second) - 1;
    }

    // cout << __FUNCTION__ << " d fn : " << _fn << ",interval : [" << interval.first << ":" << interval.second << "] rfn : " << _rfn  << ",rinterval : [" << rinterval.first << ":" << rinterval.second << "]\n" ;
    MEMCandidate mc;
    int cannum = mc.initMEMCandidate(interval, rinterval, _rfn,_q, firstMatchnum + susunda*fnum,query,queryLen,startPos,memLen,s - susunda,th,*this);
    // cout << "cannum : " << cannum << "\n";

    int count = susunda*fnum;
    while(1){
      if(cannum <= 0) break;
      if(interval.first > interval.second || _q < 0) break;
      if(query[_q] == 4) break;
      updateIntervalCompressed(interval, query[_q], _fn,cannum);
      // cout << __FUNCTION__ << " e fn : " << _fn << ",interval : [" << interval.first << ":" << interval.second << "] rfn : " << _rfn  << ",rinterval : [" << rinterval.first << ":" << rinterval.second << "]\n" ;

      // if(startPos == 165691738){
        // cout << "a_fn : " << _fn << " _q : " << _q << " [" << interval.first << ":" << interval.second << "]n";
      // }


      _q -= 1;
      _fn = (_fn + 1) % fnum;

      count++;
      if(interval.first > interval.second || _q < 0) break;
      if(query[_q] == 4) break;
      // cout << "[" << interval.first << ":" << interval.second << "] fn " << fnum << " n " << n << " query[" << _q << "] = " << query[_q] << " bef\n";
      cannum = mc.updateMEMCandidate_r(interval, query[_q], firstMatchnum + count,_fn,memLen,th,*this);
      // cout << "[" << interval.first << ":" << interval.second << "] fn " << fnum << " n " << n << " query[" << _q << "] = " << query[_q] << " cannum " << cannum << " af\n";
    }

    mc.finalizeMEMCandidate_r(firstMatchnum + count,memLen);
    *memcount += mc.getSize();

    // static int cccc = 0;
    // if(cccc == 100){
    //   exit(0);
    // }
    // cccc++;
    // if(startPos == 165691738){
    //   exit(1);
    // }

    //    cout << startPos << " : " << *memcount << "\n";
       // mc.printOutput();

    return startPos - _q + (interval.first <= interval.second) - 1;
  }
};



char vv[5] = {'N','A','C','G','T'};

void arrange(unsigned char* ref,int n){



  for(int i = 0;i < n;i++){
    switch((char)ref[i]){
    case 'a':
    case 'A':
      ref[i] = 'A';
      break;
    case 'c':
    case 'C':
      ref[i] = 'C';
      break;
    case 'g':
    case 'G':
      ref[i] = 'G';
      break;
    case 't':
    case 'T':
      ref[i] = 'T';
      break;
    default:
      ref[i] = 1;
      break;
    }
  }
}

void load_fasta(string filepath,string& ref){
  ref = "`";
  string line;
  uint64_t length = 1;

  string filename = string(filepath);

  ifstream ifs(filename.c_str());

  if(!ifs.is_open()){
    cerr << "unable to open " << filename << "\n"; exit(1);
  }

  cout << filepath << "\n";

  // int count = 0;
  
  while(!ifs.eof()){
    getline(ifs, line); // Load one line at a time.

    // count++;
    // if(count % 10000 == 0) cout << count << " : " << line << "\n";
    
    if(line.length() == 0) continue;
    if(line[line.length()-1]=='\r'){
        line.erase(--line.end());
        printf("shouldn't happen");
        exit(1);
    }

    long start = 0, end = line.length() - 1;

    // Meta tag line and start of a new sequence.
    if(line[0] == '>') {
      // Save previous sequence and meta data.
      if(length > 0) {
	ref += '`'; // ` character used to separate strings
        //	startpos.push_back(S.length());
	//lengths.push_back(length+1);
      }
      // Reset parser state.
    } else { // Collect sequence data.
      length += end - start + 1;
      for(long i = start; i <= end; i++) { 
	ref += line[i];
      }
    }
  }
  ref += "`";
}

void load_fastq(string filepath,vector<string>& ref){
  string line;
  uint64_t length = 0;

  string filename = string(filepath);

  ifstream ifs(filename.c_str());

  if(!ifs.is_open()){
    cerr << "unable to open " << filename << "\n"; exit(1);
  }

  while(!ifs.eof()){
    getline(ifs, line); // Load one line at a time.
    if(line.length() == 0) continue;
    if(line[0] != '@') continue;
    if(line[line.length()-1]=='\r'){
        line.erase(--line.end());
        printf("shouldn't happen");
        exit(1);
    }
    getline(ifs, line);
    ref.push_back(line);
  }
}



const int oneEntryNum = 128;

void outputIndex2(int argc, const char *argv[]) {
  if(argc != 5) { fprintf(stderr,"args : file fbwtnum outputDir kmer\n");exit(1); }

  string _ref;
  load_fasta(string(argv[1]), _ref);

  int fbwtnum = atoi(argv[2]);
  uint64_t nn = _ref.length();
  uint64_t n = ((nn + 1) % (fbwtnum*oneEntryNum) == 0 ? nn + 1 : nn + (fbwtnum*oneEntryNum) + 1 - ((nn + 1) % (fbwtnum*oneEntryNum)));
  const char* dirname = argv[3];
  int kmer = atoi(argv[4]);

  if(n % oneEntryNum != 0){
    fprintf(stderr,"you should input bytes multiple of %d\n",oneEntryNum);exit(1);
  }
  if(n %  fbwtnum != 0){
    fprintf(stderr,"you should input bytes multiple of %d\n",fbwtnum);exit(1);
  }

  for(uint64_t i = nn;i < n;i++){
    _ref += "`";
  }

  unsigned char* ref = (unsigned char*)strdup(_ref.c_str());

  arrange(ref, n);
  string _rref = _ref;
  reverse(_rref.begin(),_rref.end());
  unsigned char* rref = (unsigned char*)strdup(_rref.c_str());
  arrange(rref, n);

  unique_ptr<int[]> SA(new int[n]{});
  unique_ptr<int[]> RSA(new int[n]{});
  unique_ptr<unsigned char[]> fbwt(new unsigned char[n]{});
  unique_ptr<unsigned char[]> rfbwt(new unsigned char[n]{});
  
  if(sais(ref, SA.get(), (int)n) != 0) {
    fprintf(stderr, "%s: Cannot allocate memory.\n", argv[0]);
    exit(EXIT_FAILURE);
  }
  
  if(sais(rref, RSA.get(), (int)n) != 0) {
    fprintf(stderr, "%s: Cannot allocate memory.\n", argv[0]);
    exit(EXIT_FAILURE);
  }
  cout << n << " n\n";
  
  convertATGC(ref, n);
  convertATGC(rref, n);

  int count = 0;
  for(uint64_t i = 0;i < n;i++){
    if(SA[i] % fbwtnum == 0){
      SA[count] = SA[i];
      count++; 
   }
  }
  count = 0;
  for(uint64_t i = 0;i < n;i++){
    if(RSA[i] % fbwtnum == 0){
      RSA[count] = RSA[i];
      count++;
    }
  }


  cerr << "start makeFBWT " << n/fbwtnum << "\n";
  makeFBWT(SA,fbwt,ref,fbwtnum,n/fbwtnum);
  makeFBWT(RSA,rfbwt,rref,fbwtnum,n/fbwtnum);
  cerr << "end makeFBWT " << n/fbwtnum << "\n";

  C_OCC c_occ(fbwt,rfbwt,SA,n/fbwtnum,fbwtnum,ref,kmer,1,3);

  c_occ.outputFile(dirname);

  free(ref);
  free(rref);

}




void outputIndex(int argc, const char *argv[]) {
  if(argc != 6) { fprintf(stderr,"args : file bytes fbwtnum outputDir kmer\n");exit(1); }

  FILE *fp;
  const char* filename = argv[1];
  int nn = atoi(argv[2]);
  int fbwtnum = atoi(argv[3]);
  const char* dirname = argv[4];
  int kmer = atoi(argv[5]);

  cout << "file is : " << filename << ", datalen is " << nn <<  ", fbwtnum is : " << fbwtnum << "\n";


  int n = ((nn + 2) % (fbwtnum*oneEntryNum) == 0 ? nn + 2 : nn + (fbwtnum*oneEntryNum) + 2 - ((nn + 2) % (fbwtnum*oneEntryNum)));

  cout << "file is : " << filename << ", datalen is " << n << ", fbwtnum is : " << fbwtnum << "\n";

  printf("%d bytes\n",n);

  if(n % oneEntryNum != 0){
    fprintf(stderr,"you should input bytes multiple of %d\n",oneEntryNum);exit(1);
  }
  if(n %  fbwtnum != 0){
    fprintf(stderr,"you should input bytes multiple of %d\n",fbwtnum);exit(1);
  }

  assert(n % fbwtnum == 0);

  unsigned char* ref = (unsigned char*)malloc(sizeof(unsigned char)*n);
  unsigned char* rref = (unsigned char*)malloc(sizeof(unsigned char)*n);
  unique_ptr<int[]> SA(new int[n]{});
  unique_ptr<int[]> RSA(new int[n]{});
  unique_ptr<unsigned char[]> fbwt(new unsigned char[n]{});
  unique_ptr<unsigned char[]> rfbwt(new unsigned char[n]{});


  if((ref == NULL) || (SA == nullptr)) {
    fprintf(stderr, "%s: Cannot allocate memory.\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  if((fp = fopen(filename, "rb")) == NULL) {
    fprintf(stderr, "%s: Cannot open file `%s': ", argv[0], filename);
    perror(NULL);
    exit(EXIT_FAILURE);
  }

  fread(ref + 1, sizeof(unsigned char), (size_t)n, fp);
  fclose(fp);

  ref[0] = 1;
  nn++;// this is because of head padding
  for(;nn < n;nn++){
    ref[nn] = 1;
  }

  arrange(ref,n);
  cout << "arrange : " << n << "\n";

  for(int i = 0;i < n;i++){
    rref[n - 1 - i] = ref[i];
  }

  if(sais(ref, SA.get(), (int)n) != 0) {
    fprintf(stderr, "%s: Cannot allocate memory.\n", argv[0]);
    exit(EXIT_FAILURE);
  }
  
  if(sais(rref, RSA.get(), (int)n) != 0) {
    fprintf(stderr, "%s: Cannot allocate memory.\n", argv[0]);
    exit(EXIT_FAILURE);
  }
  
  convertATGC(ref, nn);
  convertATGC(rref, nn);



  // for(int i = 0;i < nn;i++){
  //   cout << (int)ref[i];
  //   if(i % 100 == 0) cout << "\n";
  // }



  int count = 0;
  for(int i = 0;i < n;i++){
    if(SA[i] % fbwtnum == 0){
      SA[count] = SA[i];
      count++; 
   }
  }
  count = 0;
  for(int i = 0;i < n;i++){
    if(RSA[i] % fbwtnum == 0){
      RSA[count] = RSA[i];
      count++;
    }
  }

  cout << "makeFBWT " << nn/fbwtnum << "\n";

  makeFBWT(SA,fbwt,ref,fbwtnum,nn/fbwtnum);
  makeFBWT(RSA,rfbwt,rref,fbwtnum,nn/fbwtnum);

  // for(int i = 0;i < nn/fbwtnum;i++){
  //   for(int j = 0;j < fbwtnum;j++){
  //    cout << vv[fbwt[j*(nn/fbwtnum) + i]] << ",";
  //   }
  //   cout << "\n";
  // }


  cout << "end makeFBWT " << nn/fbwtnum << "\n";

  C_OCC c_occ(fbwt,rfbwt,SA,nn/fbwtnum,fbwtnum,ref,kmer,1,3);
  cout << "after\n";

  c_occ.outputFile(dirname);
  free(ref);
}

double dtime(){
  struct timeval t;
  gettimeofday(&t,NULL);
  return (double)t.tv_sec + (double)t.tv_usec/1000000.0;
}


void test5(int argc, const char *argv[]){
  if(argc != 4) { fprintf(stderr,"args : inputDir query startPos\n");exit(1); }

  const char* dirname = argv[1];
  const unsigned char* query = (unsigned char*)argv[2];
  int startPos = atoi(argv[3]);
  //  const char* refname = argv[4];
  int queryLen = strlen((char*)query);
  unsigned char* _query = (unsigned char*)malloc(sizeof(unsigned char)*queryLen);
  memmove(_query, query, queryLen);
  convertATGC(_query, queryLen);

  cout << "dirname is : " << dirname << ", query is " << query <<  ", queryLen is : " << queryLen << "\n";

  C_OCC c_occ(dirname,1,0);

  int times = 100;
  for(int memLen = 20;memLen < 100;memLen++){
    countMEMS = 0;

    // memLen - K*sparse + 1 > 10
    int sparse = 10;
    while(memLen - 8*sparse + 1 < 10){
      sparse--;
    }
    if(sparse <= 0){
      cout << "illegal sparse num " << sparse << "\n";
      exit(1);
    }
    if(sparse >= 5){
      sparse = 4;
    }
    c_occ.setSparseMult(sparse);
    double tstart = dtime();
    for(int i = 0;i < times;i++)
      c_occ.FindAllMEMCompressed(_query, queryLen, startPos,memLen,5,100,&countMEMS);
    double tend = dtime();
    cerr << "time for mapping (gettimeofday): " << tend - tstart << " memLen " << memLen << " exec times " << times << " hits " << countMEMS << " sparse " << sparse << " #measure" << endl;
  }

  free(_query);
}

void test6(int argc, const char *argv[]){ // input from file
  if(argc != 7) { fprintf(stderr,"args : inputDir queryfile memlen linearcomp intervalthreashold sparseMult\n");exit(1); }

  const char* dirname = argv[1];
  const unsigned char* queryFile = (unsigned char*)argv[2];
  const int memlen = atoi(argv[3]);
  const int linearcomp = atoi(argv[4]);
  const int intervalthreashold = atoi(argv[5]);
  const int sparseMult = atoi(argv[6]);
  int queryLen;

  if(memlen + 1 - sparseMult*8 <= 10){
    return;
  }

  FILE* fp;
  if((fp = fopen((char*)queryFile, "rb")) == NULL ) {
    printf("cannot open output file %s\n",queryFile);
    exit(EXIT_FAILURE);
  }

  if(fseek(fp, 0, SEEK_END) == 0) {
    queryLen = ftell(fp)/sizeof(char);
    rewind(fp);
    if(queryLen < 0) {
      fprintf(stderr, "%s: Cannot ftell `%s': ", "inputSA", queryFile);
      perror(NULL);
      exit(EXIT_FAILURE);
    }
  } else {
    fprintf(stderr, "%s: Cannot fseek `%s': ", "inputSA", queryFile);
    perror(NULL);
    exit(EXIT_FAILURE);
  }

  unsigned char* query = (unsigned char*)malloc(sizeof(char)*queryLen);

  if(fread(query, sizeof(unsigned char), (size_t)queryLen, fp) != (size_t)queryLen) {
    fprintf(stderr, "%s: %s `%s': ",
            "inputSA",
            (ferror(fp) || !feof(fp)) ? "Cannot read from" : "Unexpected EOF in",
            queryFile);
    perror(NULL);
    exit(EXIT_FAILURE);
  }

  fclose(fp);

  convertATGC(query, queryLen);

  cerr << "dirname is : " << dirname << ", queryFile is " << queryFile <<  ", queryLen is : " << queryLen << "\n";

  C_OCC c_occ(dirname,sparseMult,0);

  cerr << "start findMEMs\n";
  double tstart = dtime();
  //  c_occ.FindMEMCompressed(interval, query, queryLen, 165145300,memLen - fnum + 1,memLen);
  c_occ.FindAllMEMCompressed(query, queryLen, queryLen,memlen,linearcomp,intervalthreashold,&countMEMS); //
  double tend = dtime();
  cout << "time for mapping (gettimeofday): " << tend - tstart << " memLen " << memlen << " linearcomp " << linearcomp << " hits " << countMEMS 
       << " intervalthreashold " << intervalthreashold << " indexFile " << dirname << " queryFile " << queryFile << " sparseMult " << sparseMult << " #measure" << endl;
  free(query);
}

void test7(int argc, const char *argv[]){ // input from file
  if(argc != 3) { fprintf(stderr,"args : inputDir queryfile\n");exit(1); }

  const char* dirname = argv[1];
  const unsigned char* queryFile = (unsigned char*)argv[2];
  int queryLen;

  FILE* fp;
  if((fp = fopen((char*)queryFile, "rb")) == NULL ) {
    printf("cannot open output file %s\n",queryFile);
    exit(EXIT_FAILURE);
  }

  if(fseek(fp, 0, SEEK_END) == 0) {
    queryLen = ftell(fp)/sizeof(char);
    rewind(fp);
    if(queryLen < 0) {
      fprintf(stderr, "%s: Cannot ftell `%s': ", "inputSA", queryFile);
      perror(NULL);
      exit(EXIT_FAILURE);
    }
  } else {
    fprintf(stderr, "%s: Cannot fseek `%s': ", "inputSA", queryFile);
    perror(NULL);
    exit(EXIT_FAILURE);
  }

  unsigned char* query = (unsigned char*)malloc(sizeof(char)*queryLen);

  if(fread(query, sizeof(unsigned char), (size_t)queryLen, fp) != (size_t)queryLen) {
    fprintf(stderr, "%s: %s `%s': ",
            "inputSA",
            (ferror(fp) || !feof(fp)) ? "Cannot read from" : "Unexpected EOF in",
            queryFile);
    perror(NULL);
    exit(EXIT_FAILURE);
  }

  fclose(fp);

  convertATGC(query, queryLen);

  cout << "dirname is : " << dirname << ", queryFile is " << queryFile <<  ", queryLen is : " << queryLen << "\n";

  C_OCC c_occ(dirname,4,0);

  int times = 10000;
  //  c_occ.setSparseMult(5);
  double tstart = dtime();
  for(int i = 0;i < times;i++){
    // unsigned char* _query = query;
    // c_occ.FindAllMEMCompressed(interval, _query, queryLen, queryLen,50);
    unsigned char* _query = query + rand()%(queryLen - 10000);
    c_occ.FindAllMEMCompressed(_query, 10000, 10000,50,5,100,&countMEMS);
  }
  double tend = dtime();
  cerr << "time for mapping (gettimeofday): " << tend - tstart << " exec times " << times 
       << " hits " << countMEMS << " #measure" << endl;


  // for(int i = 0;i < times;i++){
  //   unsigned char* _query = query + rand()%(queryLen - 100);
  //   double tstart = dtime();
  //   for(int memLen = 20;memLen < 100;memLen++){
  //     int sparse = 10;
  //     while(memLen - 8*sparse + 1 < 10){
  //       sparse--;
  //     }
  //     if(sparse <= 0){
  //       cout << "illegal sparse num " << sparse << "\n";
  //       exit(1);
  //     }
  //     c_occ.setSparseMult(sparse);
  //     c_occ.FindAllMEMCompressed(interval, _query, 100, 100,memLen);
  //   }
  //   double tend = dtime();
  //   cerr << "time for mapping (gettimeofday): " << tend - tstart << " exec times " << times 
  //        << " hits " << countMEMS << " #measure" << endl;
  // }

  free(query);
}


/*
  queryのほぼすべてのsuffixを検索するだけ
  traverseを速度を比べるために作った
*/
void test8(int argc, const char *argv[]){ // input from file
  if(argc != 3) { fprintf(stderr,"args : inputDir queryfile\n");exit(1); }

  const char* dirname = argv[1];
  const unsigned char* queryFile = (unsigned char*)argv[2];
  int queryLen;

  FILE* fp;
  if((fp = fopen((char*)queryFile, "rb")) == NULL ) {
    printf("cannot open output file %s\n",queryFile);
    exit(EXIT_FAILURE);
  }

  if(fseek(fp, 0, SEEK_END) == 0) {
    queryLen = ftell(fp)/sizeof(char);
    rewind(fp);
    if(queryLen < 0) {
      fprintf(stderr, "%s: Cannot ftell `%s': ", "inputSA", queryFile);
      perror(NULL);
      exit(EXIT_FAILURE);
    }
  } else {
    fprintf(stderr, "%s: Cannot fseek `%s': ", "inputSA", queryFile);
    perror(NULL);
    exit(EXIT_FAILURE);
  }

  unsigned char* query = (unsigned char*)malloc(sizeof(char)*queryLen);

  if(fread(query, sizeof(unsigned char), (size_t)queryLen, fp) != (size_t)queryLen) {
    fprintf(stderr, "%s: %s `%s': ",
            "inputSA",
            (ferror(fp) || !feof(fp)) ? "Cannot read from" : "Unexpected EOF in",
            queryFile);
    perror(NULL);
    exit(EXIT_FAILURE);
  }

  fclose(fp);

  convertATGC(query, queryLen);

  cout << "dirname is : " << dirname << ", queryFile is " << queryFile <<  ", queryLen is : " << queryLen << "\n";

  C_OCC c_occ(dirname,4,0);

  pair<int, int> interval;
  int times = queryLen - 100;
  //  c_occ.setSparseMult(5);
  double tstart = dtime();
  for(int i = 0;i < times;i++){
    c_occ.FBWTBackwardSearchCompressed(interval, query, 100 + i, 0);
  }
  double tend = dtime();
  cerr << "time for mapping (gettimeofday): " << tend - tstart << " exec times " << times 
       << " #measure_traverseOnly" << endl;

  free(query);
}


/*
  ヒトゲノムの解析用のもの
*/
void test9(int argc, const char *argv[]){ // input from file
  if(argc != 4) { fprintf(stderr,"args : inputDir queryfile threadNum\n");exit(1); }


  const char* dirname = argv[1];
  const unsigned char* queryFile = (unsigned char*)argv[2];
  const int THREADNUM = atoi(argv[3]);

  int queryLen;

  FILE* fp;
  if((fp = fopen((char*)queryFile, "rb")) == NULL ) {
    printf("cannot open output file %s\n",queryFile);
    exit(EXIT_FAILURE);
  }

  if(fseek(fp, 0, SEEK_END) == 0) {
    queryLen = ftell(fp)/sizeof(char);
    rewind(fp);
    if(queryLen < 0) {
      fprintf(stderr, "%s: Cannot ftell `%s': ", "inputSA", queryFile);
      perror(NULL);
      exit(EXIT_FAILURE);
    }
  } else {
    fprintf(stderr, "%s: Cannot fseek `%s': ", "inputSA", queryFile);
    perror(NULL);
    exit(EXIT_FAILURE);
  }

  unsigned char* query = (unsigned char*)malloc(sizeof(char)*queryLen);

  if(fread(query, sizeof(unsigned char), (size_t)queryLen, fp) != (size_t)queryLen) {
    fprintf(stderr, "%s: %s `%s': ",
            "inputSA",
            (ferror(fp) || !feof(fp)) ? "Cannot read from" : "Unexpected EOF in",
            queryFile);
    perror(NULL);
    exit(EXIT_FAILURE);
  }

  fclose(fp);

  convertATGC(query, queryLen);

  cout << "dirname is : " << dirname << ", queryFile is " << queryFile <<  ", queryLen is : " << queryLen << "\n";

  C_OCC c_occ(dirname,2,0);

  pair<int, int> interval;
  int times = queryLen - 100;
  //  c_occ.setSparseMult(5);
  double tstart = dtime();
  countMEMS = 0;

  vector<thread> threads(THREADNUM);

  int countmem[32];
  for(int i = 0;i < 32;i++){
    countmem[i] = 0;
  }


  for(int t = 0;t < THREADNUM;t++){
    threads[t] = thread([&queryLen,t,&query,tstart,&c_occ,THREADNUM,countmem]() mutable{
        int oneBlock = queryLen/77/THREADNUM;
        for(int i = oneBlock*t;i < oneBlock*(t+1);i++){
          unsigned char* _query = query + 77*i;
          c_occ.FindAllMEMCompressed(_query, 76, 76,20,5,5,countmem + t);
          if(i % 1000 == 0)
            cout << "thread : " << t << " i : " << i << "/" << queryLen/77 << " matches : " << countmem[t] << " time : " << dtime() - tstart<< "\n";
        }
      });
  }
  for(int t = 0;t < THREADNUM;t++){
    threads[t].join();
  }

  int totalmemcount = 0;
  for(int i = 0;i < 32;i++){
    totalmemcount += countmem[i];
  }


  double tend = dtime();
  cerr << "time for mapping (gettimeofday): " << tend - tstart << " exec times " << times << " totalmemcount " << totalmemcount
       << " #measure_humanMEM" << endl;

  free(query);
}

// forward searchのテスト
void test10(){
  C_OCC c_occ("index/test_forward",2,0);

  //  c_occ.testBackwardForward();
   // c_occ.testInitMem();
  // c_occ.testFindMEM_r();

  c_occ.printMorethanC(3);

}

/*
  ヒトゲノムの解析用のもの。RFBWTをやっている
*/
void test11(int argc, const char *argv[]){ // input from file
  if(argc != 4) { fprintf(stderr,"args : inputDir queryfile threadNum\n");exit(1); }


  const char* dirname = argv[1];
  const unsigned char* queryFile = (unsigned char*)argv[2];
  const int THREADNUM = atoi(argv[3]);

  int queryLen;

  FILE* fp;
  if((fp = fopen((char*)queryFile, "rb")) == NULL ) {
    printf("cannot open output file %s\n",queryFile);
    exit(EXIT_FAILURE);
  }

  if(fseek(fp, 0, SEEK_END) == 0) {
    queryLen = ftell(fp)/sizeof(char);
    rewind(fp);
    if(queryLen < 0) {
      fprintf(stderr, "%s: Cannot ftell `%s': ", "inputSA", queryFile);
      perror(NULL);
      exit(EXIT_FAILURE);
    }
  } else {
    fprintf(stderr, "%s: Cannot fseek `%s': ", "inputSA", queryFile);
    perror(NULL);
    exit(EXIT_FAILURE);
  }

  unsigned char* query = (unsigned char*)malloc(sizeof(char)*queryLen);

  if(fread(query, sizeof(unsigned char), (size_t)queryLen, fp) != (size_t)queryLen) {
    fprintf(stderr, "%s: %s `%s': ",
            "inputSA",
            (ferror(fp) || !feof(fp)) ? "Cannot read from" : "Unexpected EOF in",
            queryFile);
    perror(NULL);
    exit(EXIT_FAILURE);
  }

  fclose(fp);

  convertATGC(query, queryLen);

  cout << "dirname is : " << dirname << ", queryFile is " << queryFile <<  ", queryLen is : " << queryLen << "\n";

  C_OCC c_occ(dirname,2,0);

  pair<int, int> interval;
  int times = queryLen - 100;
  //  c_occ.setSparseMult(5);
  double tstart = dtime();
  countMEMS = 0;

  vector<thread> threads(THREADNUM);

  uint64_t countmem[32];
  for(int i = 0;i < 32;i++){
    countmem[i] = 0;
  }


  for(int t = 0;t < THREADNUM;t++){
    threads[t] = thread([&queryLen,t,&query,tstart,&c_occ,THREADNUM,countmem]() mutable{
        int oneBlock = queryLen/77/THREADNUM;
        // for(int i = 13000;i < 14000;i++){
        for(int i = oneBlock*t;i < oneBlock*(t+1);i++){
	// int i = 13861; 
          unsigned char* _query = query + 77*i;
          // c_occ.FindAllMEMCompressed_r(_query, 76, 76,20,10000000,countmem + t,t);
          c_occ.FindAllMEMCompressed_r(_query, 76, 76,20,10,countmem + t,t);
          if(i % 1000 == 0)
            cout << "thread : " << t << " i : " << i << " / " << queryLen/77 << " matches : " << countmem[t] << " time : " << dtime() - tstart << " expected end time : " << (double)oneBlock/(i - oneBlock*t)*(dtime() - tstart) << "\n";
          // if(i == 1000) exit(0);
        }
      });
  }
  for(int t = 0;t < THREADNUM;t++){
    threads[t].join();
  }

  int totalmemcount = 0;
  for(int i = 0;i < 32;i++){
    totalmemcount += countmem[i];
  }


  double tend = dtime();
  cerr << "time for mapping (gettimeofday): " << tend - tstart << " exec times " << times 
       << " #measure_humanMEM" << endl;

  free(query);
}

/*
  細菌用
 */
void test12(int argc, const char *argv[]){ // input from file
  if(argc != 6) { fprintf(stderr,"args : inputDir queryfile memlen intervalthreashold sparseMult\n");exit(1); }

  const char* dirname = argv[1];
  const unsigned char* queryFile = (unsigned char*)argv[2];
  const int memlen = atoi(argv[3]);
  const int intervalthreashold = atoi(argv[4]);
  const int sparseMult = atoi(argv[5]);
  int queryLen;

  if(memlen + 1 - sparseMult*8 <= 10){
    return;
  }

  FILE* fp;
  if((fp = fopen((char*)queryFile, "rb")) == NULL ) {
    printf("cannot open output file %s\n",queryFile);
    exit(EXIT_FAILURE);
  }

  if(fseek(fp, 0, SEEK_END) == 0) {
    queryLen = ftell(fp)/sizeof(char);
    rewind(fp);
    if(queryLen < 0) {
      fprintf(stderr, "%s: Cannot ftell `%s': ", "inputSA", queryFile);
      perror(NULL);
      exit(EXIT_FAILURE);
    }
  } else {
    fprintf(stderr, "%s: Cannot fseek `%s': ", "inputSA", queryFile);
    perror(NULL);
    exit(EXIT_FAILURE);
  }

  unsigned char* query = (unsigned char*)malloc(sizeof(char)*queryLen);

  if(fread(query, sizeof(unsigned char), (size_t)queryLen, fp) != (size_t)queryLen) {
    fprintf(stderr, "%s: %s `%s': ",
            "inputSA",
            (ferror(fp) || !feof(fp)) ? "Cannot read from" : "Unexpected EOF in",
            queryFile);
    perror(NULL);
    exit(EXIT_FAILURE);
  }

  fclose(fp);


  convertATGC(query, queryLen);

  cerr << "dirname is : " << dirname << ", queryFile is " << queryFile <<  ", queryLen is : " << queryLen << "\n";

  C_OCC c_occ(dirname,sparseMult,0);

  cerr << "start findMEMs\n";
  double tstart = dtime();
  //  c_occ.FindMEMCompressed(interval, query, queryLen, 165145300,memLen - fnum + 1,memLen);
  uint64_t countMEMS = 0;
  c_occ.FindAllMEMCompressed_r(query, queryLen, queryLen,memlen,intervalthreashold,&countMEMS,0); //
  double tend = dtime();
  cout << "time for mapping (gettimeofday): " << tend - tstart << " memLen " << memlen << " hits " << countMEMS 
       << " intervalthreashold " << intervalthreashold << " indexFile " << dirname << " queryFile " << queryFile << " sparseMult " << sparseMult << " #measure" << endl;
  free(query);
}

/*
  細菌用
 */
void test122(int argc, const char *argv[]){ // input from file
  if(argc != 6) { fprintf(stderr,"args : inputDir queryfile memlen intervalthreashold sparseMult\n");exit(1); }

  const char* dirname = argv[1];
  // const unsigned char* queryFile = (unsigned char*)argv[2];
  const int memlen = atoi(argv[3]);
  const int intervalthreashold = atoi(argv[4]);
  const int sparseMult = atoi(argv[5]);

  string _query;
  cout << "loading query...\n";
  load_fasta(string(argv[2]), _query);

  unsigned char* query = (unsigned char*)strdup(_query.c_str());
  int queryLen = _query.length();

  cout << "convert query...\n";
  convertATGC(query, queryLen);

  cerr << "dirname is : " << dirname << ", queryFile is " << argv[2] <<  ", queryLen is : " << queryLen << "\n";

  C_OCC c_occ(dirname,sparseMult,0);

  cerr << "start findMEMs\n";
  double tstart = dtime();
  //  c_occ.FindMEMCompressed(interval, query, queryLen, 165145300,memLen - fnum + 1,memLen);
  uint64_t countMEMS = 0;
  c_occ.FindAllMEMCompressed_r(query, queryLen, queryLen,memlen,intervalthreashold,&countMEMS,0); //
  double tend = dtime();
  cout << "time_for_mapping_(gettimeofday): " << tend - tstart << " memLen " << memlen << " hits " << countMEMS 
       << " intervalthreashold " << intervalthreashold << " indexFile " << dirname << " queryFile " << argv[2]
       << " sparseMult " << sparseMult << " fnum " << c_occ.fnum << " maxMemorySize " << maxMemorySize << " #measure" << endl;
  free(query);
}

struct Query{
  int queryLen;
  unsigned char* query;  
};

void testSSR(int argc, const char *argv[]){ // input from file
  if(argc != 6) { fprintf(stderr,"args : inputDir queryfileSSR memlen intervalthreashold sparseMult\n");exit(1); }

  const char* dirname = argv[1];
  // const unsigned char* queryFile = (unsigned char*)argv[2];
  const int memlen = atoi(argv[3]);
  const int intervalthreashold = atoi(argv[4]);
  const int sparseMult = atoi(argv[5]);
  if(sparseMult <= 0){
    cerr << "sparseMult must be over 0\n";
    exit(1);
  }

  vector<string> _query;
  cout << "loading query...\n";
  load_fastq(string(argv[2]), _query);

  cout << "convert query...\n";
  unique_ptr<Query[]> queries(new Query[_query.size()]);

  for(int i = 0;i < _query.size();i++){
    queries[i].queryLen = _query[i].length();
    queries[i].query = (unsigned char*)strdup(_query[i].c_str());
    convertATGC(queries[i].query,queries[i].queryLen);
  }

  cerr << "dirname is : " << dirname << ", queryFile is " << argv[2] <<  ", queryNum is : " << _query.size() << "\n";

  C_OCC c_occ(dirname,sparseMult,0);

  cerr << "start findMEMs\n";
  double tstart = dtime();
  //  c_occ.FindMEMCompressed(interval, query, queryLen, 165145300,memLen - fnum + 1,memLen);
  uint64_t countMEMS = 0;
  uint64_t co = 0;
  for(int i = 0;i < _query.size();i++){
    c_occ.FindAllMEMCompressed_r(queries[i].query, queries[i].queryLen, queries[i].queryLen,memlen,intervalthreashold,&countMEMS,0); //
    co++;
    if(co % 10000 == 0){
      cerr << co << " : expected time " << (double)_query.size()/i*(dtime() - tstart) << "\n";
    }
  }
  double tend = dtime();
  cout << "time_for_mapping_(gettimeofday): " << tend - tstart << " memLen " << memlen << " hits " << countMEMS 
       << " intervalthreashold " << intervalthreashold << " indexFile " << dirname << " queryFile " << argv[2]
       << " sparseMult " << sparseMult << " fnum " << c_occ.fnum << " maxMemorySize " << maxMemorySize << " hashsize " 
       << c_occ.kmerSize << " #measure" << endl;
  for(int i = 0;i < _query.size();i++){
    free(queries[i].query);
  }
}


const string C_OCC::fbwtfname = "fbwt";
const string C_OCC::rfbwtfname = "rfbwt";
//vector<pair<int,int> > C_OCC::MEMCandidate::output;

int main(int argc, const char *argv[]) {

  // for(int i = 0;i < INTERVAL_DIST_MAX_LEN*32;i++){
  //   initMEMIntervalfLenDist[i] = 0;
  //   updateMEMIntervalLenDist[i] = 0; 
  // }

  // outputIndex(argc,argv);
  // test6(argc,argv);
  // test7(argc,argv);
  // test8(argc,argv);
  // test9(argc,argv);
  // test10();
  // test11(argc,argv);
  // test12(argc,argv);
  // outputIndex2(argc,argv);
  //test122(argc,argv);
   testSSR(argc,argv);
  // for(int i = 0;i < INTERVAL_DIST_MAX_LEN;i++){
  //   for(int th = 1;th < 32;th++){
  //     initMEMIntervalfLenDist[i] += initMEMIntervalfLenDist[th*INTERVAL_DIST_MAX_LEN + i];
  //     updateMEMIntervalLenDist[i] += updateMEMIntervalLenDist[th*INTERVAL_DIST_MAX_LEN + i];
  //   }
  // }

  // cout << "index , initMEM , updateMEM\n";
  // for(int i = 0;i < INTERVAL_DIST_MAX_LEN;i++){
  //   cout << (i + 1) << " , " << initMEMIntervalfLenDist[i] << " , " << updateMEMIntervalLenDist[i] << " #index_initMEM_updateMEM\n";
  // }  

}

