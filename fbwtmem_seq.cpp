#include <iostream>
#include <memory>
#include <utility>
#include <vector>
#include <stdint.h>
#include <fstream>
//#include <random>
#include <stdlib.h>

#include "sais/sais.h"
#include <string.h>
#include "makeFBWT_seq.h"
#include <assert.h>
#include <sys/mman.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/time.h>

#include <thread>
#include <algorithm>
#include <array>

#include <stdio.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>

using namespace std;

#define ACCREF(ii) (c_occ.ref[(ii)/2] >> (((ii)&1) << 2) & 0x7)

double dtime(){
  struct timeval t;
  gettimeofday(&t,NULL);
  return (double)t.tv_sec + (double)t.tv_usec/1000000.0;
}


long show_getrusage(){
  struct rusage r;
  if (getrusage(RUSAGE_SELF, &r) != 0) {
    cerr << "err\n";
  }
  return r.ru_maxrss;
}


struct Match{
  Match(uint64_t r,uint64_t q,uint64_t l):refpos(r),querypos(q),len(l){}
  uint64_t refpos;
  uint64_t querypos;
  uint64_t len;
};

#define TEST_INTV_OCC_MAX 1024*1024

// #define TEST_INTV_OCC

#ifdef TEST_INTV_OCC
struct TestIntv{
  int num;
  double firstmatch;
  double initmemcandidate;
  double updatememcandidate;
  int nobi;
  static int backwardnum;
};

TestIntv testintv[TEST_INTV_OCC_MAX];
#endif

#define CHAR_N 5

uint64_t countMEMS = 0;
long maxMemorySize = 0;

class C_OCC{

public:
  unique_ptr< unsigned char[]> FBWT;

  class C{
  private:
    unique_ptr<uint32_t[]> ctable;
    int32_t fnum;

    const string fname = "c";

  public:
    C(const char* basename){
      FILE* fp;
      string filename = string(basename) + "." + fname;


      if((fp = fopen(filename.c_str(), "rb")) == NULL ) {
        cerr << "fopen fail " << filename << "\n";
        exit(EXIT_FAILURE);
      }
      int64_t n;
      if(fseek(fp, 0, SEEK_END) == 0) {
        n = ftell(fp)/sizeof(uint32_t);
        rewind(fp);
        if(n < 0) {
          cerr << "ftell fail in " << filename << "\n";
          exit(EXIT_FAILURE);
        }
      } else {
        cerr << "fseek fail in " << filename << "\n";
        exit(EXIT_FAILURE);
      }

      ctable = make_unique<uint32_t[]>(n);

      if(fread(ctable.get(), sizeof(uint32_t), (size_t)n, fp) != (size_t)n) {
        cerr << "read fail " << filename << "\n";
        exit(EXIT_FAILURE);
      }

      fnum = n/(CHAR_N + 1);
      fclose(fp);
    }


    C(unique_ptr< unsigned char[]>& FBWT,uint64_t n,int32_t fn){
      ctable = make_unique<uint32_t[]>((CHAR_N + 1)*fn);
      this->fnum = fn;

      for(int64_t i = 0;i < (CHAR_N + 1)*fn;i++){
        ctable[i] = 0;
      }

      for(int32_t f = 0;f < fn;f++){
        for(uint64_t i = 0;i < n;i++){
          ctable[FBWT[i + f*n] + (CHAR_N + 1)*f]++;
        }
      }

      for(int32_t f = 0;f < fn;f++){
        for(uint64_t i = 1;i < CHAR_N;i++){
          ctable[i + f*(CHAR_N + 1)] += ctable[i - 1 + f*(CHAR_N + 1)];
        }
        for(uint64_t i = CHAR_N;i > 0;i--){
          ctable[i + f*(CHAR_N + 1)] = ctable[i - 1 + f*(CHAR_N + 1)];
        }
        ctable[f*(CHAR_N + 1)] = 0;
      }
    }


    int getFnum(){
      return fnum;
    }

    void printC(){
      for(int32_t f = 0;f < fnum;f++){
        for(uint64_t i = 0;i <= CHAR_N;i++){
          cerr << ctable[i + f*(CHAR_N + 1)] << " ";
        }
        cerr << "\n";
      }
    }

    uint32_t getC(uint64_t c,int32_t fn){
      return ctable[fn*(CHAR_N + 1) + c];
    }

    void outputFile(const char* basename){
      FILE* fp;
      string filename = string(basename) + "." + fname;
      if((fp = fopen(filename.c_str(), "wb")) == NULL ) {
        cerr << "fopen fail " << filename << "\n";
        exit(EXIT_FAILURE);
      }
      fwrite(ctable.get() , sizeof(uint32_t), (CHAR_N + 1)*fnum, fp);
      fclose(fp);
    }
  };


  class OCC{
    //  private:
  public:
    C_OCC& c_occ;
    int32_t fnum;
    uint64_t n;

    uint64_t oneFBWTOcctablesize;

    const string fname = "occ";

    static const uint64_t oneEntryNumOfChar = 128;

    struct oneEntry{
      uint64_t pre[3];
      uint32_t atgc[4];
      uint64_t post[3];
    };

    unique_ptr<oneEntry[]> occtable;

    void constructOcctable(unique_ptr<unsigned char[]>& FBWT,int32_t fn,uint64_t n,uint64_t oneFBWTOcctablesize){
      occtable = make_unique<oneEntry[]>(oneFBWTOcctablesize*fnum);

      for(int64_t f = 0;f < fn;f++){
        uint32_t count[5] = {0,0,0,0,0};

        for(uint64_t i = 0;i < oneFBWTOcctablesize;i++){
          uint64_t tmp = 0;
          for(uint64_t j = 0;j < 21;j++){
            tmp |= ((uint64_t)(FBWT[i*oneEntryNumOfChar + j + f*n] & 0x7) << (j*3));
            count[FBWT[i*oneEntryNumOfChar + j + f*n]]++;
          }
          tmp |= ((uint64_t)(FBWT[i*oneEntryNumOfChar + 21 + f*n] & 0x7) << (21*3));
          occtable[i + oneFBWTOcctablesize*f].pre[0] = tmp;
          count[FBWT[i*oneEntryNumOfChar + 21 + f*n] & 0x7]++;

          tmp = 0;
          tmp |= ((uint64_t)(FBWT[i*oneEntryNumOfChar + 21 + f*n] & 0x7) >> 1);
          for(uint64_t j = 1;j < 21;j++){
            tmp |= ((uint64_t)(FBWT[i*oneEntryNumOfChar + j + 1*21 + f*n] & 0x7) << ((j - 1)*3 + 2));
            count[FBWT[i*oneEntryNumOfChar + j + 21 + f*n]]++;
          }
          tmp |= ((uint64_t)(FBWT[i*oneEntryNumOfChar + 42 + f*n] & 0x7) << (21*3 - 1));
          occtable[i + oneFBWTOcctablesize*f].pre[1] = tmp;
          count[FBWT[i*oneEntryNumOfChar + 42 + f*n] & 0x7]++;

          tmp = 0;
          tmp |= ((uint64_t)(FBWT[i*oneEntryNumOfChar + 42 + f*n] & 0x7) >> 2);

          for(uint64_t j = 1;j < 22;j++){
            tmp |= ((uint64_t)(FBWT[i*oneEntryNumOfChar + j + 42 + f*n] & 0x7) << ((j - 1)*3 + 1));
            count[FBWT[i*oneEntryNumOfChar + j + 42 + f*n]]++;

          }
          occtable[i + oneFBWTOcctablesize*f].pre[2] = tmp;


          occtable[i + oneFBWTOcctablesize*f].atgc[0] = count[1];
          occtable[i + oneFBWTOcctablesize*f].atgc[1] = count[2];
          occtable[i + oneFBWTOcctablesize*f].atgc[2] = count[3];
          occtable[i + oneFBWTOcctablesize*f].atgc[3] = count[4];

          tmp = 0;
          for(uint64_t j = 0;j < 21;j++){
            tmp |= ((uint64_t)(FBWT[i*oneEntryNumOfChar + oneEntryNumOfChar/2 + j + f*n] & 0x7) << (j*3));
            count[FBWT[i*oneEntryNumOfChar + oneEntryNumOfChar/2 + j + f*n]]++;
          }
          tmp |= ((uint64_t)(FBWT[i*oneEntryNumOfChar + oneEntryNumOfChar/2 + 21 + f*n] & 0x7) << (21*3));
          occtable[i + oneFBWTOcctablesize*f].post[0] = tmp;
          count[FBWT[i*oneEntryNumOfChar + oneEntryNumOfChar/2 + 21 + f*n] & 0x7]++;

          tmp = 0;
          tmp |= ((uint64_t)(FBWT[i*oneEntryNumOfChar + oneEntryNumOfChar/2 + 21 + f*n] & 0x7) >> 1);
          for(uint64_t j = 1;j < 21;j++){
            tmp |= ((uint64_t)(FBWT[i*oneEntryNumOfChar + oneEntryNumOfChar/2 + j + 1*21 + f*n] & 0x7) << ((j - 1)*3 + 2));
            count[FBWT[i*oneEntryNumOfChar + oneEntryNumOfChar/2 + j + 1*21 + f*n]]++;
          }
          tmp |= ((uint64_t)(FBWT[i*oneEntryNumOfChar + oneEntryNumOfChar/2 + 42 + f*n] & 0x7) << (21*3 - 1));
          occtable[i + oneFBWTOcctablesize*f].post[1] = tmp;
          count[FBWT[i*oneEntryNumOfChar + oneEntryNumOfChar/2 + 42 + f*n] & 0x7]++;

          tmp = 0;
          tmp |= ((uint64_t)(FBWT[i*oneEntryNumOfChar + oneEntryNumOfChar/2 + 42 + f*n] & 0x7) >> 2);

          for(uint64_t j = 1;j < 22;j++){
            tmp |= ((uint64_t)(FBWT[i*oneEntryNumOfChar + oneEntryNumOfChar/2 + j + 2*21 + f*n] & 0x7) << ((j - 1)*3 + 1));
            count[FBWT[i*oneEntryNumOfChar + oneEntryNumOfChar/2 + j + 2*21 + f*n]]++;
          }
          occtable[i + oneFBWTOcctablesize*f].post[2] = tmp;

        }

      }
    }

    char vv[5] = {'N','A','C','G','T'};

    void print64bit(uint64_t h){
      for(int k = 0;k < 32;k++)
        cerr << vv[(h >> (k*2)) & 0x3];
    }


    void printOneOcctableEntry(int i){

      for(int k = 0;k < 21;k++)
        cerr << vv[(occtable[i].pre[0] >> (k*3)) & 0x7];
      int tmp = ((occtable[i].pre[0] >> (63)) & 0x1) | (((occtable[i].pre[1]) & 0x3) << 1);
      cerr << vv[tmp];
      cerr << " ";

      for(int k = 0;k < 20;k++)
        cerr << vv[(occtable[i].pre[1] >> (k*3 + 2)) & 0x7];

      tmp = ((occtable[i].pre[1] >> (62)) & 0x3) | (((occtable[i].pre[2]) & 0x1) << 2);
      cerr << vv[tmp];
      cerr << " ";

      for(int k = 0;k < 21;k++)
        cerr << vv[(occtable[i].pre[2] >> (k*3 + 1)) & 0x7];
      cerr << " ";

      for(int j = 0;j < 4;j++){
        cerr << " " << occtable[i].atgc[j];
      }      

      cerr << " ";

      for(int k = 0;k < 21;k++)
        cerr << vv[(occtable[i].post[0] >> (k*3)) & 0x7];
      tmp = ((occtable[i].post[0] >> (63)) & 0x1) | (((occtable[i].post[1]) & 0x3) << 1);
      cerr << vv[tmp];
      cerr << " ";

      for(int k = 0;k < 20;k++)
        cerr << vv[(occtable[i].post[1] >> (k*3 + 2)) & 0x7];

      tmp = ((occtable[i].post[1] >> (62)) & 0x3) | (((occtable[i].post[2]) & 0x1) << 2);
      cerr << vv[tmp];
      cerr << " ";

      for(int k = 0;k < 21;k++)
        cerr << vv[(occtable[i].post[2] >> (k*3 + 1)) & 0x7];

      cerr << "\n";
    }


    uint64_t atgcMask[4] = {0x9249249249249249,0x2492492492492492,0xb6db6db6db6db6db,0x4924924924924924};

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


  public:

    void printOcctable(){

      for(int32_t f = 0;f < fnum;f++){
        for(uint64_t i = 0;i < oneFBWTOcctablesize;i++){
          printOneOcctableEntry(i + f*oneFBWTOcctablesize);
        }
        cerr << "-------------------------------------\n";
      }
    }

    void printBits(uint64_t a){
      for(int i = 0;i < 64;i++){
        cerr << ((a >> (63 - i)) & 0x1);
        if(i % 3 == 0) cerr << " ";
      }
      cerr << '\n';
    }


    uint32_t getOcc(uint64_t i,unsigned char c,int32_t fn){

      uint64_t _atgcMask = atgcMask[c - 1];

      if(i == n){
        uint32_t p = (n - 1) / oneEntryNumOfChar;
        uint32_t ret = occtable[p + fn*oneFBWTOcctablesize].atgc[c - 1];
        uint64_t tmp = (~((occtable[p + fn*oneFBWTOcctablesize].post[2] >> 0x1) ^ _atgcMask));

        tmp = tmp & (tmp >> 1) & (tmp >> 2) & 0x1249249249249249;
        tmp = (tmp + (tmp >> 3)) & 0x71c71c71c71c71c7;
        tmp = (tmp + (tmp >> 6)) & 0x703f03f03f03f03f;
        tmp = (tmp + (tmp >> 12)) & 0x0fff000fff000fff;
        tmp = (tmp + (tmp >> 24)) & 0x0fff000000ffffff;
        tmp = (tmp + (tmp >> 48)) & 0x0000ffffffffffff;

        ret += tmp;

        uint64_t tmp2 = (~(((occtable[p + fn*oneFBWTOcctablesize].post[1] >> 2) 
                            | ((occtable[p + fn*oneFBWTOcctablesize].post[2] & 0x1) << 62)) ^ _atgcMask));

        tmp2 = tmp2 & (tmp2 >> 1) & (tmp2 >> 2) & 0x1249249249249249;
        tmp2 = (tmp2 + (tmp2 >> 3)) & 0x71c71c71c71c71c7;
        tmp2 = (tmp2 + (tmp2 >> 6)) & 0x703f03f03f03f03f;
        tmp2 = (tmp2 + (tmp2 >> 12)) & 0x0fff000fff000fff;
        tmp2 = (tmp2 + (tmp2 >> 24)) & 0x0fff000000ffffff;
        tmp2 = (tmp2 + (tmp2 >> 48)) & 0x0000ffffffffffff;

        ret += tmp2;

        uint64_t tmp3 = (~((occtable[p + fn*oneFBWTOcctablesize].post[0]) ^ _atgcMask));
        tmp3 = tmp3 & (tmp3 >> 1) & (tmp3 >> 2) & 0x1249249249249249;
        tmp3 = (tmp3 + (tmp3 >> 3)) & 0x71c71c71c71c71c7;
        tmp3 = (tmp3 + (tmp3 >> 6)) & 0x703f03f03f03f03f;
        tmp3 = (tmp3 + (tmp3 >> 12)) & 0x0fff000fff000fff;
        tmp3 = (tmp3 + (tmp3 >> 24)) & 0x0fff000000ffffff;
        tmp3 = (tmp3 + (tmp3 >> 48)) & 0x0000ffffffffffff;

        ret += (((occtable[p + fn*oneFBWTOcctablesize].post[0] >> 63) 
                | ((occtable[p + fn*oneFBWTOcctablesize].post[1] & 0x3) << 1)) == c);

        ret += tmp3;
        return ret;
      }

      uint32_t p = i / oneEntryNumOfChar;
      uint32_t ret = occtable[p + fn*oneFBWTOcctablesize].atgc[c - 1];
      uint32_t num = (i % oneEntryNumOfChar);


      if(num < 21){
        uint64_t tmp = (~(occtable[p + fn*oneFBWTOcctablesize].pre[0] ^ _atgcMask)) & preMask[num];
        uint64_t tmp2 = (~(((occtable[p + fn*oneFBWTOcctablesize].pre[1] << 1) 
                            | (occtable[p + fn*oneFBWTOcctablesize].pre[0] >> 63)) ^ _atgcMask));
        uint64_t tmp3 = (~((occtable[p + fn*oneFBWTOcctablesize].pre[2] >> 1) ^ _atgcMask));


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

        ret -= (((occtable[p + fn*oneFBWTOcctablesize].pre[1] >> 62) | ((occtable[p + fn*oneFBWTOcctablesize].pre[2] & 0x1) << 0x2)) == c);
      }else if(num == 21){

        uint64_t tmp2 = (~(((occtable[p + fn*oneFBWTOcctablesize].pre[1] << 1) 
                            | (occtable[p + fn*oneFBWTOcctablesize].pre[0] >> 63)) ^ _atgcMask));

        tmp2 = tmp2 & (tmp2 >> 1) & (tmp2 >> 2) & 0x1249249249249249;
        tmp2 = (tmp2 + (tmp2 >> 3)) & 0x71c71c71c71c71c7;
        tmp2 = (tmp2 + (tmp2 >> 6)) & 0x703f03f03f03f03f;
        tmp2 = (tmp2 + (tmp2 >> 12)) & 0x0fff000fff000fff;
        tmp2 = (tmp2 + (tmp2 >> 24)) & 0x0fff000000ffffff;
        tmp2 = (tmp2 + (tmp2 >> 48)) & 0x0000ffffffffffff;

        ret -= tmp2;

        uint64_t tmp3 = (~((occtable[p + fn*oneFBWTOcctablesize].pre[2] >> 1) ^ _atgcMask));
        tmp3 = tmp3 & (tmp3 >> 1) & (tmp3 >> 2) & 0x1249249249249249;
        tmp3 = (tmp3 + (tmp3 >> 3)) & 0x71c71c71c71c71c7;
        tmp3 = (tmp3 + (tmp3 >> 6)) & 0x703f03f03f03f03f;
        tmp3 = (tmp3 + (tmp3 >> 12)) & 0x0fff000fff000fff;
        tmp3 = (tmp3 + (tmp3 >> 24)) & 0x0fff000000ffffff;
        tmp3 = (tmp3 + (tmp3 >> 48)) & 0x0000ffffffffffff;

        ret -= tmp3;

        ret -= (((occtable[p + fn*oneFBWTOcctablesize].pre[1] >> 62) | ((occtable[p + fn*oneFBWTOcctablesize].pre[2] & 0x1) << 2)) == c);

      }else if(num <= 42){
        uint64_t tmp2 = (~(((occtable[p + fn*oneFBWTOcctablesize].pre[1] >> 2) 
                            | ((occtable[p + fn*oneFBWTOcctablesize].pre[2] & 0x1) << 62)) ^ _atgcMask)) & preMask[num - 22];

        tmp2 = tmp2 & (tmp2 >> 1) & (tmp2 >> 2) & 0x1249249249249249;
        tmp2 = (tmp2 + (tmp2 >> 3)) & 0x71c71c71c71c71c7;
        tmp2 = (tmp2 + (tmp2 >> 6)) & 0x703f03f03f03f03f;
        tmp2 = (tmp2 + (tmp2 >> 12)) & 0x0fff000fff000fff;
        tmp2 = (tmp2 + (tmp2 >> 24)) & 0x0fff000000ffffff;
        tmp2 = (tmp2 + (tmp2 >> 48)) & 0x0000ffffffffffff;

        ret -= tmp2;

        uint64_t tmp3 = (~((occtable[p + fn*oneFBWTOcctablesize].pre[2] >> 1) ^ _atgcMask));
        tmp3 = tmp3 & (tmp3 >> 1) & (tmp3 >> 2) & 0x1249249249249249;
        tmp3 = (tmp3 + (tmp3 >> 3)) & 0x71c71c71c71c71c7;
        tmp3 = (tmp3 + (tmp3 >> 6)) & 0x703f03f03f03f03f;
        tmp3 = (tmp3 + (tmp3 >> 12)) & 0x0fff000fff000fff;
        tmp3 = (tmp3 + (tmp3 >> 24)) & 0x0fff000000ffffff;
        tmp3 = (tmp3 + (tmp3 >> 48)) & 0x0000ffffffffffff;

        ret -= tmp3;
      }else if(num < 64){

        uint64_t tmp3 = (~((occtable[p + fn*oneFBWTOcctablesize].pre[2] >> 1) ^ _atgcMask)) & preMask[num - 43];

        tmp3 = tmp3 & (tmp3 >> 1) & (tmp3 >> 2) & 0x1249249249249249;
        tmp3 = (tmp3 + (tmp3 >> 3)) & 0x71c71c71c71c71c7;
        tmp3 = (tmp3 + (tmp3 >> 6)) & 0x703f03f03f03f03f;
        tmp3 = (tmp3 + (tmp3 >> 12)) & 0x0fff000fff000fff;
        tmp3 = (tmp3 + (tmp3 >> 24)) & 0x0fff000000ffffff;
        tmp3 = (tmp3 + (tmp3 >> 48)) & 0x0000ffffffffffff;

        ret -= tmp3;
      }else if(num <= 85){
        uint64_t tmp = (~(occtable[p + fn*oneFBWTOcctablesize].post[0] ^ _atgcMask)) & postMask[num - 64];

        tmp = tmp & (tmp >> 1) & (tmp >> 2) & 0x1249249249249249;
        tmp = (tmp + (tmp >> 3)) & 0x71c71c71c71c71c7;
        tmp = (tmp + (tmp >> 6)) & 0x703f03f03f03f03f;
        tmp = (tmp + (tmp >> 12)) & 0x0fff000fff000fff;
        tmp = (tmp + (tmp >> 24)) & 0x0fff000000ffffff;
        tmp = (tmp + (tmp >> 48)) & 0x0000ffffffffffff;

        ret += tmp;

      }else if(num == 86){
        uint64_t tmp = (~(occtable[p + fn*oneFBWTOcctablesize].post[0] ^ _atgcMask));

        tmp = tmp & (tmp >> 1) & (tmp >> 2) & 0x1249249249249249;
        tmp = (tmp + (tmp >> 3)) & 0x71c71c71c71c71c7;
        tmp = (tmp + (tmp >> 6)) & 0x703f03f03f03f03f;
        tmp = (tmp + (tmp >> 12)) & 0x0fff000fff000fff;
        tmp = (tmp + (tmp >> 24)) & 0x0fff000000ffffff;
        tmp = (tmp + (tmp >> 48)) & 0x0000ffffffffffff;

        ret += tmp;

        ret += (((occtable[p + fn*oneFBWTOcctablesize].post[0] >> 63) | ((occtable[p + fn*oneFBWTOcctablesize].post[1] & 0x3) << 1)) == c);

      }else if(num <= 106){
        uint64_t tmp2 = (~(((occtable[p + fn*oneFBWTOcctablesize].post[1] << 1) 
                            | ((occtable[p + fn*oneFBWTOcctablesize].post[0]) >> 63)) ^ _atgcMask)) & postMask[num - 85];

        tmp2 = tmp2 & (tmp2 >> 1) & (tmp2 >> 2) & 0x1249249249249249;
        tmp2 = (tmp2 + (tmp2 >> 3)) & 0x71c71c71c71c71c7;
        tmp2 = (tmp2 + (tmp2 >> 6)) & 0x703f03f03f03f03f;
        tmp2 = (tmp2 + (tmp2 >> 12)) & 0x0fff000fff000fff;
        tmp2 = (tmp2 + (tmp2 >> 24)) & 0x0fff000000ffffff;
        tmp2 = (tmp2 + (tmp2 >> 48)) & 0x0000ffffffffffff;

        ret += tmp2;

        uint64_t tmp3 = (~((occtable[p + fn*oneFBWTOcctablesize].post[0]) ^ _atgcMask));
        tmp3 = tmp3 & (tmp3 >> 1) & (tmp3 >> 2) & 0x1249249249249249;
        tmp3 = (tmp3 + (tmp3 >> 3)) & 0x71c71c71c71c71c7;
        tmp3 = (tmp3 + (tmp3 >> 6)) & 0x703f03f03f03f03f;
        tmp3 = (tmp3 + (tmp3 >> 12)) & 0x0fff000fff000fff;
        tmp3 = (tmp3 + (tmp3 >> 24)) & 0x0fff000000ffffff;
        tmp3 = (tmp3 + (tmp3 >> 48)) & 0x0000ffffffffffff;

        ret += tmp3;
      }else if(num < 128){

        uint64_t tmp = (~((occtable[p + fn*oneFBWTOcctablesize].post[2] >> 0x1) ^ _atgcMask)) & postMask[num - 107];

        tmp = tmp & (tmp >> 1) & (tmp >> 2) & 0x1249249249249249;
        tmp = (tmp + (tmp >> 3)) & 0x71c71c71c71c71c7;
        tmp = (tmp + (tmp >> 6)) & 0x703f03f03f03f03f;
        tmp = (tmp + (tmp >> 12)) & 0x0fff000fff000fff;
        tmp = (tmp + (tmp >> 24)) & 0x0fff000000ffffff;
        tmp = (tmp + (tmp >> 48)) & 0x0000ffffffffffff;

        ret += tmp;

        ret += (((occtable[p + fn*oneFBWTOcctablesize].post[1] >> 62)
                 | ((occtable[p + fn*oneFBWTOcctablesize].post[2] & 0x1) << 2)) == c);

        uint64_t tmp2 = (~(((occtable[p + fn*oneFBWTOcctablesize].post[1] << 1) 
                            | ((occtable[p + fn*oneFBWTOcctablesize].post[0]) >> 63)) ^ _atgcMask));

        tmp2 = tmp2 & (tmp2 >> 1) & (tmp2 >> 2) & 0x1249249249249249;
        tmp2 = (tmp2 + (tmp2 >> 3)) & 0x71c71c71c71c71c7;
        tmp2 = (tmp2 + (tmp2 >> 6)) & 0x703f03f03f03f03f;
        tmp2 = (tmp2 + (tmp2 >> 12)) & 0x0fff000fff000fff;
        tmp2 = (tmp2 + (tmp2 >> 24)) & 0x0fff000000ffffff;
        tmp2 = (tmp2 + (tmp2 >> 48)) & 0x0000ffffffffffff;

        ret += tmp2;

        uint64_t tmp3 = (~((occtable[p + fn*oneFBWTOcctablesize].post[0]) ^ _atgcMask));
        tmp3 = tmp3 & (tmp3 >> 1) & (tmp3 >> 2) & 0x1249249249249249;
        tmp3 = (tmp3 + (tmp3 >> 3)) & 0x71c71c71c71c71c7;
        tmp3 = (tmp3 + (tmp3 >> 6)) & 0x703f03f03f03f03f;
        tmp3 = (tmp3 + (tmp3 >> 12)) & 0x0fff000fff000fff;
        tmp3 = (tmp3 + (tmp3 >> 24)) & 0x0fff000000ffffff;
        tmp3 = (tmp3 + (tmp3 >> 48)) & 0x0000ffffffffffff;

        ret += tmp3;
      }
      return ret;
    }

    OCC(const char* basename,int32_t fn,C_OCC& cocc):c_occ(cocc){
      string filename = string(basename) + "." + fname;
      FILE* fp;
      this->fnum = fn;

      if((fp = fopen(filename.c_str(), "rb")) == NULL ) {
        cerr << "fopen fail " << filename << "\n";
        exit(EXIT_FAILURE);
      }

      int64_t n;
      if(fseek(fp, 0, SEEK_END) == 0) {
        n = ftell(fp)/sizeof(oneEntry);
        rewind(fp);
        if(n < 0) {
          cerr << "ftell fail " << filename << "\n";
          exit(EXIT_FAILURE);
        }
      } else {
        cerr << "fseek fail " << filename << "\n";
        exit(EXIT_FAILURE);
      }
      occtable = make_unique<oneEntry[]>(n);      
      if(fread(occtable.get(), sizeof(oneEntry), (size_t)n, fp) != (size_t)n) {
        cerr << "fread fail " << filename << "\n";
        exit(EXIT_FAILURE);
      }

      this->fnum = fn;
      this->oneFBWTOcctablesize = n/fn;
      this->n = this->oneFBWTOcctablesize*oneEntryNumOfChar;
      fclose(fp);

      cerr << "Occurrence table size is " <<  oneFBWTOcctablesize*oneEntryNumOfChar << " \n";
    }

    // n is size of one FBWT
    OCC(unique_ptr< unsigned char[]>& FBWT,uint64_t n,int32_t fn,C_OCC& cocc):c_occ(cocc){
    
      oneFBWTOcctablesize = (n % oneEntryNumOfChar == 0) ? (n / oneEntryNumOfChar) : (n / oneEntryNumOfChar + 1);

      this->n = n;
      this->fnum = fn;

      constructOcctable(FBWT,fn,n,oneFBWTOcctablesize);

    }


    uint64_t getN(){
      return n;
    }

    void printOCCCompressed(int c){
      for(uint32_t i = 0;i <= n;i++){
        cerr << i << ",";
        for(int32_t f = 0;f < fnum;f++){
          cerr << getOcc(i,c,f) << ",";
        }
        cerr << "\n";
      }
    }

    void printOCCCompressed(){
      for(uint32_t i = 0;i <= n;i++){
        cerr << i << ",";
        for(int32_t f = 0;f < fnum;f++){
          for(uint32_t c = 1;c < CHAR_N;c++)
            cerr << getOcc(i,c,f) << ",";
        }
        cerr << "\n";
      }
    }

    void outputFile(const char* basename){
      FILE* fp;
      string filename = string(basename) + "." + fname;
      if((fp = fopen(filename.c_str(), "wb")) == NULL ) {
        cerr << "fopen fail " << filename << "\n";
        exit(EXIT_FAILURE);
      }
      fwrite(occtable.get() , sizeof(oneEntry), oneFBWTOcctablesize*fnum, fp);
      fclose(fp);
    }
    
  };

  class tuple{
  public:
    uint32_t left;
    uint32_t right;
  };

  uint32_t power(uint32_t d,uint32_t i){
    if(i == 0) return 1;
    else return d*power(d,i - 1);
  }


  void _constructKMR(uint32_t key,uint32_t depth,uint32_t maxdepth,int32_t fn){
    if(maxdepth == depth){
      KMR[key + kmerOneFnSize*fn].left = 1;
      KMR[key + kmerOneFnSize*fn].right = 0;
      return;
    }
    for(uint32_t i = 0;i < 4;i++){
      _constructKMR(key | (i << 2*depth), depth + 1, maxdepth,(fn + 1)%fnum);
    }
  }

  void _constructKMR(C_OCC& c_occ,const pair<uint32_t,uint32_t>& interval,uint32_t key,uint32_t depth,const uint32_t maxdepth,int32_t fn){
    if(interval.first > interval.second){
      _constructKMR(key, depth, maxdepth,fn);
      return;
    }
    if(maxdepth == depth){
      KMR[key + kmerOneFnSize*fn].left = interval.first;
      KMR[key + kmerOneFnSize*fn].right = interval.second;

      return;
    }
    for(uint32_t i = 0;i < 4;i++){
      pair<uint32_t,uint32_t> in(interval);
      c_occ.updateInterval(in, i + 1, fn);
      _constructKMR(c_occ, in, key | (i << 2*depth), depth + 1, maxdepth,(fn + 1)%fnum);
    }
  }

  void constructKMR(C_OCC& c_occ,uint32_t num){
    kmerOneFnSize = power(4,num);
    KMR = make_unique<tuple[]>(kmerOneFnSize*fnum);
    
    for(int32_t f = 0;f < fnum;f++){
      for(uint32_t i = 1;i < 5;i++){
        pair<uint32_t,uint32_t> interval(c->getC(i,f),c->getC(i + 1,f) - 1);
        _constructKMR(c_occ, interval, i - 1, 1, num,(f + 1)%fnum);
      }
    }

    // printKMR();
  }

  void printKMR(){
    char vv[4] = {'A','C','G','T'};
    cerr << "=========================== KMR ===========================\n";
    for(uint32_t i = 0;i < kmerOneFnSize;i++){
      for(int64_t j = kmerSize - 1;j >= 0;j--){
        cerr << vv[((i >> j*2) & 0x3)];
      }
      for(int32_t f = 0;f < fnum;f++){
        cerr << " : [" << KMR[i + kmerOneFnSize*f].left << ":" << KMR[i + kmerOneFnSize*f].right << "]\t";
      }
      cerr << "\n";
    }
  }

  void outputKMR(const char* basename){
    FILE* fp;
    string filename = string(basename) + "." + kmrfname;
    if((fp = fopen(filename.c_str(), "wb")) == NULL ) {
      cerr << "fopen fail " << filename << "\n";
      exit(EXIT_FAILURE);
    }
    fwrite(KMR.get() , sizeof(tuple), kmerOneFnSize*fnum, fp);
    fclose(fp);

    uint32_t o[2] = {kmerSize,kmerOneFnSize};

    filename = string(basename) + "." + kmrsizefname;
    if((fp = fopen(filename.c_str(), "wb")) == NULL ) {
      cerr << "fopen fail " << filename << "\n";
      exit(EXIT_FAILURE);
    }
    fwrite(o , sizeof(uint32_t), 2, fp);
    fclose(fp);
  }

  void inputKMR(const char* basename){
    string filename = string(basename) + "." + kmrfname;
    FILE* fp;

    if((fp = fopen(filename.c_str(), "rb")) == NULL ) {
      cerr << "fopen fail " << filename << "\n";
      exit(EXIT_FAILURE);
    }

    int32_t n;
    if(fseek(fp, 0, SEEK_END) == 0) {
      n = ftell(fp)/sizeof(tuple);
      rewind(fp);
      if(n < 0) {
        cerr << "ftell fail " << filename << "\n";
        exit(EXIT_FAILURE);
      }
    } else {
      cerr << "fseek fail " << filename << "\n";
      exit(EXIT_FAILURE);
    }

    KMR = make_unique<tuple[]>(n);

    if(fread(KMR.get(), sizeof(tuple), (size_t)n, fp) != (size_t)n) {
      cerr << "fread fail " << filename << "\n";
      exit(EXIT_FAILURE);
    }

    fclose(fp);

    uint32_t o[2];
    filename = string(basename) + "." + kmrsizefname;
    if((fp = fopen(filename.c_str(), "rb")) == NULL ) {
      cerr << "fopen fail " << filename << "\n";
      exit(EXIT_FAILURE);
    }
    if(fread(o, sizeof(int), 2, fp) != 2) {
      cerr << "read fail " << filename << "\n";
      exit(EXIT_FAILURE);
    }

    kmerSize = o[0];
    kmerOneFnSize = o[1];

  }

  // int FBWTBackwardSearch(pair<int,int>* interval,const unsigned char* query,int q,int fn){
  //   interval->first = c->getC(query[q - 1],fn);
  //   interval->second = c->getC(query[q - 1] + 1,fn) - 1;

  //   int _q = q - 1;
  //   int _fn = (fn + 1) % fnum;

  //   if(verbose & LOG_INTERVAL){cerr << "query [" << _q  << "] = " << (int)query[_q] << " fn : " << _fn << " [" << interval->first << ":" << interval->second << "]\n";}

  //   while(interval->first <= interval->second && _q > 0){
  //     _q -= 1;
  //     updateIntervalCompressed(interval, query[_q], _fn);
  //     if(verbose & LOG_INTERVAL){cerr << "query [" << _q  << "] = " << (int)query[_q] << " fn : " << _fn << " [" << interval->first << ":" << interval->second << "]\n";}
  //     _fn = (_fn + 1) % fnum;
  //   }
  //   return q - _q + (interval->first <= interval->second) - 1;
  // }

  unique_ptr<tuple[]> KMR;//到着点に入っている
  uint32_t kmerSize;
  uint32_t kmerOneFnSize;
  unique_ptr<OCC> occ;
  unique_ptr<C> c;
  unique_ptr<uint32_t[]> sa;
  vector<pair<uint32_t,string>> refdescription;
  unsigned char* ref;

  int32_t fnum;
  int verbose = 0;
  uint64_t n;

  int sparseMult;

  void printOutput(vector<Match>& match,bool fourcolumn){
    if(fourcolumn){
      for(auto& m : match){
        pair<uint32_t,string> tmp(m.refpos,"");
        int id = (int)(lower_bound(refdescription.begin(),refdescription.end(),tmp) - refdescription.begin() - 1);
        cout << refdescription[id].second << " " << m.refpos << "  " << m.querypos << " " << m.len << "\n";
      }
    }else{
      for(auto& m : match){
        cout << m.refpos << "  " << m.querypos << " " << m.len <<" # refpos querypos memlen\n";
      }
    }
  }


  class MEMCandidate{

    vector<int> memcandidate;
    //    static vector<pair<int,int> > output;
    // vector<pair<uint32_t,uint32_t> > output;
    vector<Match>& output;
    uint32_t memCandidateNum = 0;
    bool fourcolumn = false;

  public:
    MEMCandidate(vector<Match>&  match):output(match){}

    void setPrintFourcolumn(bool f){
      fourcolumn = f;
    }

    void printMEMCandidate(){
      cerr << "memcandidate : ";
      for(uint32_t i = 0;i < memCandidateNum;i++){
        cerr << memcandidate[i] << ", \t";
        if(i % 9 == 0) cerr << "\n";
      }
      cerr << "\n";
    }



    void printSize(){
      cerr << "output size is " << output.size() << "\n";
    }

    int getSize(){
      return output.size();
    }


    // void check(unsigned char* query,C_OCC& c_occ,pair<int,int> p,int queryPos){
    //   while(true){
        
    //   }
    // }


    bool countForward(unsigned char* query,uint32_t queryLen,uint32_t pos,uint32_t sa,uint32_t len,uint32_t memLen,int s,C_OCC& c_occ){
      bool ret = false;
      for(uint32_t i = 0;i < (uint32_t)c_occ.fnum*s;i++){
        if(pos + i >= queryLen || sa + len + i >= c_occ.n*c_occ.fnum){
          if(len + i >= memLen) {
            output.push_back(Match(sa,pos - len,len + i));
	    // cerr << "sa : " << sa << " len : " << len + i << "\n";
            // cerr << "countForward 1\n";
          }
          return ret;
        }else if(query[pos + i] == 0){
          if(len + i >= memLen) {
            output.push_back(Match(sa,pos - len,len + i));
	    // cerr << "sa : " << sa << " len : " << len + i << "\n";
            // cerr << "countForward 2, " << c_occ.fnum*s << "," << i << "," << len << " \n";
          }
          return i;
        }else if(query[pos + i] != ACCREF(sa + len + i)){ //(c_occ.ref[(sa + len + i)/2] >> (((sa + len + i)&1) << 2) & 0x15)
          if(len + i >= memLen) {
            output.push_back(Match(sa,pos - len,len + i));
	    // cerr << "sa : " << sa << " len : " << len + i << "\n";
            // cerr << "countForward 3, " << c_occ.fnum*s << "," << i << "," << len << " \n";
          }
          return i;
        }
        ret = true;
      }
      return ret;
    }

    bool countForwardAndBackward(unsigned char* query,uint64_t queryLen,uint32_t pos,uint32_t sa,uint32_t len,uint32_t memLen,int s,C_OCC& c_occ){

      // static int cc = 0;
      // if(cc > 100) exit(1);
      // cc++;

      // backward
      int64_t backwardLen = 0;
      while(true){
        if(((int64_t)pos - (int64_t)len - backwardLen - 1 < 0) || ((int64_t)sa - backwardLen - 1 < 0)){
          break;
        }
        if(query[pos - len - backwardLen - 1] == 0){
          break;
        }
        if(query[pos - len - backwardLen - 1] != ACCREF(sa - backwardLen - 1)){ //c_occ.ref[sa - backwardLen - 1]
          break;
        }else{
          backwardLen++;
        }
      }

#ifdef TEST_INTV_OCC
      TestIntv::backwardnum++;
#endif
      return countForward(query, queryLen, pos, sa - backwardLen, len + backwardLen, memLen, s, c_occ);
    }

    uint32_t linerlycountthreashold = 5;

    // 返り値は残った候補数
    uint32_t initMEMCandidate(pair<uint32_t,uint32_t>& interval,int _q,uint32_t len,unsigned char* query,uint32_t queryLen,uint32_t startPos,uint32_t memLen,int s,uint32_t linerlycountthreashold,C_OCC& c_occ){

      this->linerlycountthreashold = linerlycountthreashold;
      if(interval.second - interval.first + 1 <= linerlycountthreashold){
        for(uint32_t i = interval.first;i <= interval.second;i++){
          countForwardAndBackward(query, queryLen, startPos, c_occ.sa[i],len,memLen,s, c_occ);
        }
        memCandidateNum = 0;
        return 0;
      }

      uint32_t count = 0;

      if(_q < 0){
        // cerr << "A!\n";
        for(uint32_t i = interval.first;i <= interval.second;i++){
          countForward(query, queryLen, startPos, c_occ.sa[i],len,memLen,s, c_occ);
        }
        return 0;
      }
      // for(int i = interval.first;i <= interval.second;i++){
      //   if(c_occ.FBWT[i] != query[_q]){
      //     countForward(query, queryLen, startPos, c_occ.sa[i],len,memLen,s, c_occ);
      //   }else{
      //     memcandidate.push_back(c_occ.sa[i] == 0 ? c_occ.n - 1 : c_occ.sa[i] - 1);
      //     count++;
      //   }
      // }
#define MEMCAN(a)                                                       \
      if((a) != query[_q]){                                             \
        countForward(query,queryLen,startPos,c_occ.sa[i],len,memLen,s,c_occ); \
      }else{                                                            \
        memcandidate.push_back(c_occ.sa[i]);                            \
        count++;                                                        \
      }
      
      // initMEMIntervalfLenDist[th*INTERVAL_DIST_MAX_LEN + ((interval.second - interval.first + 1) < INTERVAL_DIST_MAX_LEN - 1 ? (interval.second - interval.first + 1) : INTERVAL_DIST_MAX_LEN - 1)]++;

      uint32_t pad = interval.first % OCC::oneEntryNumOfChar;
      OCC::oneEntry* e = &c_occ.occ->occtable[interval.first/OCC::oneEntryNumOfChar];
      for(uint32_t i = interval.first;i <= interval.second;i++){
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
        if(pad == OCC::oneEntryNumOfChar){
          pad = 0;
          e++;
        }
      }

#undef MEMCAN


      memCandidateNum = count;
      return count;
    }

    uint32_t updateMEMCandidate(pair<uint32_t,uint32_t>& interval,unsigned char q,uint32_t firstmatchnum,uint32_t backwardlen,int fn,unsigned char* query,uint64_t queryLen,uint32_t startPos,uint32_t memLen,int s,uint32_t firstintervallen,C_OCC& c_occ){
      uint32_t len = firstmatchnum + backwardlen;
      if(interval.second - interval.first + 1 <= linerlycountthreashold ||interval.second - interval.first + 1  <= firstintervallen/100){
        for(uint32_t i = 0;i <= interval.second - interval.first;i++){
          countForwardAndBackward(query, queryLen, startPos, memcandidate[i] - backwardlen,len,memLen,s, c_occ);
        }
        memCandidateNum = 0;
        return 0;
      }

      uint32_t count = 0;

      // for(int i = interval.first,j = 0;i <= interval.second;i++,j++){
      //   if(c_occ.FBWT[fn*c_occ.n + i] != q){
      //     countForward(query, queryLen, startPos, memcandidate[j],len,memLen,s, c_occ);
      //   }else{
      //     memcandidate[count] = (memcandidate[j] == 0 ? c_occ.n - 1: memcandidate[j] - 1);
      //     count++;
      //   }
      // }

      // int ret = (c_occ.occ->getOcc(interval.second + 1,q,fn) - c_occ.occ->getOcc(interval.first,q,fn));
      // int remain = (interval.second - interval.first + 1) - ret;

#define MEMCAN(a)                                                       \
      if((a) != q){                                             \
	countForward(query, queryLen, startPos, memcandidate[count] - backwardlen < 0 ? memcandidate[count] - backwardlen + c_occ.n : memcandidate[count] - backwardlen,len,memLen,s, c_occ); \
	memcandidate.erase(memcandidate.begin() + count); \
      }else{								\
	count++;							\
      }


      //	memcandidate[count] = (memcandidate[j] == 0 ? c_occ.n - 1: memcandidate[j] - 1); 
      //remain--;if(remain == 0) {memCandidateNum = ret; return ret;

      uint32_t pad = interval.first % OCC::oneEntryNumOfChar;
      OCC::oneEntry* e = &c_occ.occ->occtable[interval.first/OCC::oneEntryNumOfChar + fn*c_occ.occ->oneFBWTOcctablesize];
      for(uint32_t i = interval.first,j = 0;i <= interval.second;i++,j++){
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
        if(pad == OCC::oneEntryNumOfChar){
          pad = 0;
          e++;
        }
      }

#undef MEMCAN

      
      memCandidateNum = count;
      return count;
    }

    void finalizeMEMCandidate(uint32_t firstmatchnum,uint32_t backwardlen,unsigned char* query,uint64_t queryLen,uint32_t startPos,uint32_t memLen,int s,C_OCC& c_occ){
      uint32_t len = firstmatchnum + backwardlen;
      // cerr << "finaleize << " << memCandidateNum << "\n";
      for(uint32_t i = 0;i < memCandidateNum;i++){
        countForward(query, queryLen, startPos, memcandidate[i] - backwardlen,len,memLen,s, c_occ);
      }
      // cerr << "finaleize end\n";
      memCandidateNum = 0;
    }
  };

  void construct(unique_ptr< unsigned char[]>& FBWT,unique_ptr<uint32_t[]>& SA,uint64_t n,int32_t fn,unsigned char* ref,uint32_t kmersize,uint32_t sparseMult){
    occ = make_unique<OCC>(OCC(FBWT,n,fn,*this));
    c = make_unique<C>(C(FBWT,n,fn));
    fnum = fn;
    this->kmerSize = kmersize;
    this->FBWT = move(FBWT);
    this->sa = move(SA);
    this->n = n;
    this->ref = ref;
    this->sparseMult = sparseMult;
    constructKMR(*this, kmersize);
  }

  void construct(unique_ptr< unsigned char[]>& FBWT,unique_ptr<uint32_t[]>& SA,uint64_t n,int32_t fn,unsigned char* ref,uint32_t kmersize,uint32_t sparseMult,int verbose){
    construct(FBWT,SA,n,fn,ref,kmersize,sparseMult);
    this->verbose = verbose;
  }

  static const int LOG_INTERVAL = 1;
  static const int LOG_MEM = 2;

  const string safname = "sa";
  const string reffname = "ref";
  const string kmrfname = "kmr";
  const string kmrsizefname = "kmrsize";
  const string descriptionfname = "description";

  static const string fbwtfname;


  void outputDescription(const char* basename){
    FILE* fp;
    string filename = string(basename) + "." + descriptionfname;
    if((fp = fopen(filename.c_str(), "wb")) == NULL ) {
      cerr << "fopen fail " << filename << "\n";
      exit(EXIT_FAILURE);
    }
    uint32_t size = refdescription.size();
    fwrite(&size , sizeof(size), 1, fp);
    for(uint32_t i = 0;i < size;i++){
      fwrite(&refdescription[i].first,sizeof(refdescription[i].first),1,fp);
      uint32_t len = refdescription[i].second.length();
      fwrite(&len,sizeof(uint32_t),1,fp);
      fwrite(refdescription[i].second.c_str(),sizeof(char),refdescription[i].second.length(),fp);
    }
    fclose(fp);
  }
  void inputDescription(const char* basename){
    string filename = string(basename) + "." + descriptionfname;
    FILE* fp;

    if((fp = fopen(filename.c_str(), "rb")) == NULL ) {
      cerr << "fopen fail " << filename << "\n";
      exit(EXIT_FAILURE);
    }

    uint32_t size;
    if(fread(&size, sizeof(uint32_t), (size_t)1, fp) != (size_t)1){
      cerr << "fread fail " << filename << "\n";
      exit(1);
    }
    // refdescription.resize(size);
    for(uint32_t i = 0;i < size;i++){
      uint32_t first;
      if(fread(&first, sizeof(first), (size_t)1, fp) != (size_t)1){
        cerr << "fread fail " << filename << "\n";
        exit(1);
      }
      char buf[1024];
      uint32_t len;
      if(fread(&len, sizeof(len), (size_t)1, fp) != (size_t)1){
        cerr << "fread fail " << filename << "\n";
        exit(1);
      }
      if(fread(buf, sizeof(char), (size_t)len, fp) != (size_t)len){
        cerr << "fread fail " << filename << "\n";
        exit(1);
      }
      string second(buf,len);
      refdescription.push_back(make_pair(first, second));
    }
    fclose(fp);

  }

  void outputFBWT(const char* basename){
    FILE* fp;
    string filename = string(basename) + "." + fbwtfname;
    if((fp = fopen(filename.c_str(), "wb")) == NULL ) {
      cerr << "fopen fail " << filename << "\n";
      exit(EXIT_FAILURE);
    }
    fwrite(FBWT.get() , sizeof(unsigned char), n*fnum, fp);
    fclose(fp);
  }


  int inputFBWT(const char* basename){
    string filename = string(basename) + "." + fbwtfname;
    FILE* fp;

    if((fp = fopen(filename.c_str(), "rb")) == NULL ) {
      cerr << "fopen fail " << filename << "\n";
      exit(EXIT_FAILURE);
    }

    int64_t n;
    if(fseek(fp, 0, SEEK_END) == 0) {
      n = ftell(fp);
      rewind(fp);
      if(n < 0) {
        cerr << "ftell fail " << filename << "\n";
        exit(EXIT_FAILURE);
      }
    } else {
      cerr << "fseek fail " << filename << "\n";
      exit(EXIT_FAILURE);
    }

    FBWT = make_unique<unsigned char[]>(n);      

    if(fread(FBWT.get(), sizeof(unsigned char), (size_t)n, fp) != (size_t)n) {
      cerr << "fread fail " << filename << "\n";
      exit(EXIT_FAILURE);
    }

    fclose(fp);

    return n;
  }

  void outputSA(const char* basename){
    FILE* fp;
    string filename = string(basename) + "." + safname;
    if((fp = fopen(filename.c_str(), "wb")) == NULL ) {
      cerr << "fopen fail " << filename << "\n";
      exit(EXIT_FAILURE);
    }
    fwrite(sa.get() , sizeof(uint32_t), n, fp);
    fclose(fp);
  }

  void outputRef(const char* basename){
    FILE* fp;
    string filename = string(basename) + "." + reffname;
    if((fp = fopen(filename.c_str(), "wb")) == NULL ) {
      cerr << "fopen fail " << filename << "\n";
      exit(EXIT_FAILURE);
    }
    fwrite(ref , sizeof(char), n*fnum/2, fp);
    fclose(fp);
  }

  void inputSA(const char* basename){
    string filename = string(basename) + "." + safname;
    FILE* fp;

    if((fp = fopen(filename.c_str(), "rb")) == NULL ) {
      cerr << "fopen fail " << filename << "\n";
      exit(EXIT_FAILURE);
    }

    int64_t n;
    if(fseek(fp, 0, SEEK_END) == 0) {
      n = ftell(fp)/sizeof(int);
      rewind(fp);
      if(n < 0) {
        cerr << "ftell fail " << filename << "\n";
        exit(EXIT_FAILURE);
      }
    } else {
      cerr << "fseek fail " << filename << "\n";
      exit(EXIT_FAILURE);
    }

    sa = make_unique<uint32_t[]>(n);

    if(fread(sa.get(), sizeof(int), (size_t)n, fp) != (size_t)n) {
      cerr << "fread fail " << filename << "\n";
      exit(EXIT_FAILURE);
    }

    fclose(fp);
  }

  void inputRef(const char* basename){
    string filename = string(basename) + "." + reffname;
    FILE* fp;

    if((fp = fopen(filename.c_str(), "rb")) == NULL ) {
      cerr << "fopen fail " << filename << "\n";
      exit(EXIT_FAILURE);
    }

    int64_t n;
    if(fseek(fp, 0, SEEK_END) == 0) {
      n = ftell(fp)/sizeof(char);
      rewind(fp);
      if(n < 0) {
        cerr << "ftell fail " << filename << "\n";
        exit(EXIT_FAILURE);
      }
    } else {
      cerr << "fseek fail " << filename << "\n";
      exit(EXIT_FAILURE);
    }

    ref = (unsigned char*)malloc(sizeof(char)*n);

    if(fread(ref, sizeof(unsigned char), (size_t)n, fp) != (size_t)n) {
      cerr << "fread fail " << filename << "\n";
      exit(EXIT_FAILURE);
    }

    fclose(fp);
  }

public:
  // n is one element of FBWT
  C_OCC(unique_ptr< unsigned char[]>& FBWT,unique_ptr<uint32_t[]>& SA,uint64_t n,int32_t fn,uint32_t kmersize,unsigned char* ref,int sparseMult,vector<pair<uint32_t,string>>& refdescription):refdescription(refdescription){
    construct(FBWT,SA,n,fn,ref,kmersize,sparseMult);
  }

  C_OCC(unique_ptr< unsigned char[]>& FBWT,unique_ptr<uint32_t[]>& SA,uint64_t n,int32_t fn,uint32_t kmersize,unsigned char* ref,int sparseMult,vector<pair<uint32_t,string>>& refdescription,int verbose):refdescription(refdescription){
    construct(FBWT,SA,n,fn,ref,kmersize,sparseMult,verbose);
  }

  C_OCC(const char* basename,int sparseMult,int verbose){
    c = make_unique<C>(C(basename));
    occ = make_unique<OCC>(OCC(basename,c->getFnum(),*this));
    n = occ->getN();
    inputSA(basename);
    // inputFBWT(basename);
    inputRef(basename);
    inputKMR(basename);
    // inputDescription(basename);
    this->verbose = verbose;
    this->fnum = c->getFnum();
    this->sparseMult = sparseMult;
  }


  void outputFile(const char* basename){
    c->outputFile(basename);
    occ->outputFile(basename);
    outputSA(basename);
    // outputFBWT(basename);
    outputRef(basename);
    outputKMR(basename);
    outputDescription(basename);
  }

  void setSparseMult(int s){
    sparseMult = s;
  }

  static int getFBWTSize(const char* basename){
    string filename = string((char*)basename) + "." + fbwtfname;
    FILE* fp;
    if((fp = fopen(filename.c_str(), "rb")) == NULL ) {
      cerr << "fopen fail " << filename << "\n";
      exit(EXIT_FAILURE);
    }
    int n;
    if(fseek(fp, 0, SEEK_END) == 0) {
      n = ftell(fp);
      rewind(fp);
      if(n < 0) {
        cerr << "ftell fail " << filename << "\n";
        exit(EXIT_FAILURE);
      }
    } else {
      cerr << "fseek fail " << filename << "\n";
      exit(EXIT_FAILURE);
    }
    fclose(fp);

    return n;
  }

  void updateInterval(pair<uint32_t,uint32_t>& interval,char q,int32_t fn){
    interval.first = c->getC(q,fn) + occ->getOcc(interval.first,q,fn);
    interval.second = c->getC(q,fn) + occ->getOcc(interval.second + 1,q,fn) - 1;
  }

  void updateInterval(pair<uint32_t,uint32_t>& interval,char q,int32_t fn,uint64_t cannum){
    interval.first = c->getC(q,fn) + occ->getOcc(interval.first,q,fn);
    interval.second = interval.first + cannum - 1;
  }


  void testOcc(){
    occ->printOcctable();
  }

  void printOCCCompressed(int c){
    occ->printOCCCompressed(c);
  }

  void printOCCCompressed(){
    occ->printOCCCompressed();
  }

  void printC(){
    c->printC();
  }


  int FBWTBackwardSearch(pair<uint32_t,uint32_t>& interval,const unsigned char* query,int q,int fn){
    interval.first = c->getC(query[q - 1],fn);
    interval.second = c->getC(query[q - 1] + 1,fn) - 1;

    int _q = q - 1;
    int _fn = (fn + 1) % fnum;

    while(interval.first <= interval.second && _q > 0 && query[_q] != 0){
      _q -= 1;
      updateInterval(interval, query[_q], _fn);
      _fn = (_fn + 1) % fnum;
      cerr << "fn " << _fn << " [" << interval.first << ":" << interval.second << "]\n";
      // if(_fn == 0) cerr << "sa : " << sa[interval.first] << "\n";
    }
    return q - _q + (interval.first <= interval.second) - 1;
  }

  void FindAllMEM(unsigned char* query,uint32_t q,uint64_t startPos,int64_t memLen,int linearcomp,const uint32_t intervalMaxSize,unique_ptr<vector<Match>[]>& matches,int threadnum){
    int maxsparseMult = (memLen - kmerSize) / fnum;
    if(sparseMult > maxsparseMult){
      cerr << "skip is too large, so shrinked to " << maxsparseMult << "\n";
      sparseMult = maxsparseMult;
    }


    if(threadnum == 1){
      for(int32_t k = 0;k < fnum;k++){
	for(int64_t i = startPos - k;i >= (int)(memLen - fnum*sparseMult + 1);i-=(int)(fnum*sparseMult)){
	  FindMEM(query, q, i,memLen - fnum*sparseMult + 1,memLen,sparseMult,linearcomp,intervalMaxSize,matches[0]);
	  long m = show_getrusage();
	  if(maxMemorySize < m){
	    maxMemorySize = m;
	  }
	}
      }
    }else{
    
      vector<thread> threads(threadnum);

      for(int t = 0;t < threadnum;t++){
	// cerr << t << ", " << threadnum << " t\n";
	threads[t] = thread([&,t]{
	    int64_t oneblocksize = (fnum*sparseMult)*(((startPos) - (memLen - fnum*sparseMult + 1))/(fnum*sparseMult)/threadnum);
	    // cerr << "oneblocksize " << oneblocksize << "\n";
	    for(int32_t k = 0;k < fnum;k++){
	      // cerr << t << " : " << threadnum << " : " << (int)(memLen - fnum*sparseMult + 1 - 1) << " : " 
	      //      << (int64_t)(startPos - k - oneblocksize*(t + 1)) << " : "
	      //      << ((t != (threadnum - 1) ? (int64_t)(startPos - k - oneblocksize*(t + 1)) : 
	      //          (int)(memLen - fnum*sparseMult + 1 - 1))) << "\n";
	      // for(int64_t i = startPos - k;i >= (int)(memLen - fnum*sparseMult + 1);i-=(int)(fnum*sparseMult)){
	      for(int64_t i = startPos - k - oneblocksize*t;
		  i > (t != threadnum - 1 ? (int64_t)(startPos - k - oneblocksize*(t + 1)) : 
		       (int)(memLen - fnum*sparseMult + 1 - 1));
		  i-=(int)(fnum*sparseMult)){
		FindMEM(query, q, i,memLen - fnum*sparseMult + 1,memLen,sparseMult,linearcomp,intervalMaxSize,matches[t]);
		long m = show_getrusage();
		if(maxMemorySize < m){
		  maxMemorySize = m;
		}
	      }
	    }
	  });
      }
      for(int t = 0;t < threadnum;t++){
	threads[t].join();
      }
    }

  }


  int FindMEM(unsigned char* query,uint64_t queryLen,uint64_t startPos,uint32_t firstMatchnum,const uint32_t memLen,const int s,const int linearcomp,const uint32_t intervalMaxSize,vector<Match>& match){

    pair<uint32_t,uint32_t> interval;
    if(startPos > queryLen) return 0;

#ifdef TEST_INTV_OCC
    double ttstart = dtime();
#endif
    
    uint32_t key = 0;
    for(uint32_t i = 0;i < kmerSize;i++){
      if(query[startPos - 1 - i] == 0) {// cerr << "query is 0 in hash\n";
        return 0;}
      key |= ((query[startPos - 1 - i] - 1) << 2*i);
    }

    int64_t _q = startPos - 1 - kmerSize;
    uint32_t _fn = (fnum - (firstMatchnum % fnum) + kmerSize) % fnum;//(firstFn + 1) % fnum;

    interval.first = KMR[key + _fn*kmerOneFnSize].left;
    interval.second = KMR[key + _fn*kmerOneFnSize].right;

    if(interval.first > interval.second){
      // cerr << "query is 0 in hash2\n";
      return 0;
    }

    for(int i = firstMatchnum - 1 - kmerSize;
        interval.first <= interval.second && _q >= 0 && i >= 0;
        _q--,i--,_fn = ((_fn + 1) % fnum)){
      if(query[_q] == 0) {// cerr << "query is 0 in backwardsearch\n";
        return 0;}

      updateInterval(interval, query[_q], _fn);
      // cerr << "_q " << _q << " [" << interval.first << ":" << interval.second << "] #backwardsearch\n";
    }

    if(interval.first > interval.second || _q != (int64_t)(startPos - 1 - firstMatchnum) || _fn != 0) {
      // cerr << "there is no interval\n";
      return 0;
    }

    int susunda = 0;
    bool isBreak = false;
    while((interval.second - interval.first + 1 > intervalMaxSize) && (susunda < (s/2 + (s&1)) - 1) && (_q >= 0)){
    // while((interval.second - interval.first + 1 > intervalMaxSize) && (susunda < s - 1) && (_q >= 0)){
      susunda++;
      for(int i = fnum - 1;
          interval.first <= interval.second && _q >= 0 && i >= 0;
          _q--,i--,_fn = ((_fn + 1) % fnum)){
        if(query[_q] == 0) {isBreak = true;break;}
        updateInterval(interval, query[_q], _fn);
        // cerr << "_q " << _q << " [" << interval.first << ":" << interval.second << "] #susu " << susunda << "\n";
      }

      if(interval.first > interval.second || _fn != 0) {
        isBreak = true;
      }

      if(isBreak) break;
    }

#ifdef TEST_INTV_OCC
    double firstmatchtime = dtime() - ttstart;
#endif
    // cerr << "ababa : " << susunda << "\n";

    if(susunda > 0){
      FindMEM(query,queryLen,startPos + fnum*(s - susunda),firstMatchnum + (s - susunda)*fnum,memLen,susunda,linearcomp,intervalMaxSize,match);
    }

    if(interval.first > interval.second || _fn != 0 || isBreak) {
      return 0;
    }

#ifdef TEST_INTV_OCC
    double tstart = dtime();
    int diffnum = interval.second - interval.first + 1;
#endif


    int64_t __q = _q;
    uint32_t __fn = _fn;


    uint32_t firstintervallen = interval.second - interval.first + 1;
    
#define INTERVAL 1000
    for(uint32_t i = interval.first;i <= interval.second;i += INTERVAL) {
      _q = __q;
      _fn = __fn;
      pair<uint32_t,uint32_t> _interval(i,i + INTERVAL - 1 < interval.second ? i + INTERVAL - 1 : interval.second);
      MEMCandidate mc(match);
      // cerr << firstMatchnum << "," << susunda << "," << fnum << "\n";
      uint32_t cannum = mc.initMEMCandidate(_interval, _q, firstMatchnum + susunda*fnum,query,
                                       queryLen,startPos,memLen,s - susunda,linearcomp,*this);
      uint32_t count = susunda*fnum;
      int32_t _firstintervallen = _interval.second - _interval.first + 1;

#ifdef TEST_INTV_OCC
      double _tstart = dtime();
      if(diffnum < TEST_INTV_OCC_MAX - 1){testintv[diffnum].num++;}
      else{testintv[diffnum].num++;}
#endif


      while(1){
        // if(diffnum - cannum < 1024*1024 - 1){testintv[diffnum - cannum]++;}
        // else{testintv[1024*1024 - 1]++;}
      
        if(cannum <= 0) break;
        if(_interval.first > _interval.second || _q < 0) break;
        if(query[_q] == 0) break;
        updateInterval(_interval, query[_q], _fn,cannum);
	
        // cerr << "_q " << _q << " [" << _interval.first << ":" << _interval.second
        //      << "] " << mc.getSize() << " #updateintervalcompressed\n";

        _q -= 1;
        _fn = (_fn + 1) % fnum;

        count++;
        if(_interval.first > _interval.second || _q < 0) break;
        if(query[_q] == 0) break;

        // diffnum = interval.second - interval.first + 1;
        // if(diffnum < TEST_INTV_OCC_MAX - 1){testintv[diffnum]++;}
        // else{testintv[diffnum]++;}

        cannum = mc.updateMEMCandidate(_interval, query[_q], firstMatchnum + susunda*fnum ,count - susunda*fnum,_fn,query,queryLen,startPos,memLen,s - susunda,_firstintervallen,*this);
      }

      mc.finalizeMEMCandidate(firstMatchnum + susunda*fnum, count - susunda*fnum,query,queryLen,startPos,memLen,s - susunda,*this);


#ifdef TEST_INTV_OCC
      double tend = dtime();
      testintv[diffnum].firstmatch += firstmatchtime;
      testintv[diffnum].initmemcandidate += (_tstart - tstart);
      testintv[diffnum].updatememcandidate += (tend - _tstart);
      testintv[diffnum].nobi += count - susunda*fnum;
#endif
    }
    //    countMEMS = mc.getSize();
    return firstintervallen;
  }
};


unsigned char convertATGC(unsigned char a){
  switch(a){
  case 'A':
  case 'a':
    return 1;
  break;
  case 'C':
  case 'c':
    return 2;
  break;
  case 'G':
  case 'g':
    return 3;
  break;
  case 'T':
  case 't':
    return 4;
  break;
  default:
    return 0;
  }
}

void convertATGC(unsigned char* ref,uint64_t n){
  for(uint32_t i = 0;i < n;i++){
    ref[i] = convertATGC(ref[i]);
  }
}

char vv[5] = {'N','A','C','G','T'};

void arrange(unsigned char* ref,uint64_t n){
  for(uint32_t i = 0;i < n;i++){
    switch(ref[i]){
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

uint64_t convert_fasta(unsigned char* ref,uint64_t n,vector<pair<uint32_t,string>>& refdescription){
  uint64_t num = 0;
  bool isSkip = false;
  string description = "";
  for(uint64_t i = 0;i < n;i++){
    if(isSkip){
      if(ref[i] == '\n') {
        refdescription.push_back(make_pair(num,description));
        description = "";
        isSkip = false; continue;
      }
      else description += ref[i];
    }else{
      if(ref[i] == '>'){
	ref[num] = '`'; num++;
	isSkip = true; continue;
      }else if(ref[i] == '\n'){
	continue;
      }else{
	ref[num] = ref[i]; num++;
      }
    }
  }
  ref[num] = '`'; num++;
  return num;
}



const int oneEntryNumOfChar = 128;


void compressRef(unsigned char* a,uint64_t n){
  for(uint64_t i = 0;i < n;i+=2){
    a[i/2] = a[i] | (a[i + 1] << 4);
  }
}

void convertSortedIntArray(unsigned char* ref,unique_ptr<uint32_t[]>& a,unique_ptr<uint32_t[]>& a2,uint32_t n,uint32_t sparse,uint32_t& sigma){
  uint32_t num = 1;
  for(uint32_t i = 0;i < n - 1;i++){
    bool isSame = true;
    for(uint32_t j = 0;j < sparse;j++){
      if(ref[a[i] + j] != ref[a[i + 1] + j]) {isSame = false; break;}
    }
    a2[a[i]/sparse] = num;
    if(!isSame) num++;
  }
  a2[a[n - 1]/sparse] = num;
  sigma = num;
}

unique_ptr<uint32_t[]> convertRefbySparce(unsigned char* ref,uint64_t n,uint32_t sparse,uint32_t& sigma){
  cerr << "start convertrefbysparse\n";
  uint32_t newsize = n/sparse;
  unique_ptr<uint32_t[]> newref(new uint32_t[newsize]);
  unique_ptr<uint32_t[]> newref2(new uint32_t[newsize]);
  for(uint32_t i = 0;i < newsize;i++){
    newref[i] = i*sparse;
  }
  cerr << "start radix sort sparse " << sparse << "\n";
  
  // radix sort sort
  for(int radix = sparse - 1;radix >= 0 ;radix -= 2){
    uint64_t count[256];
    uint64_t count2[256];
    for(uint32_t i = 0;i < 256;i++) {count[i] = 0;count2[i] = 0;}
    for(uint64_t i = 0;i < newsize;i++) {
      count[ref[newref[i] + radix]]++;
    }
    for(uint32_t i = 1;i < 256;i++) count2[i] += count2[i - 1] + count[i - 1];
    for(uint64_t i = 0;i < newsize;i++) {
      newref2[count2[ref[newref[i] + radix]]] = newref[i];
      count2[ref[newref[i] + radix]]++;
    }

    if(radix == 0) {
      convertSortedIntArray(ref,newref2,newref,newsize,sparse,sigma);
      return move(newref);
    }

    for(uint32_t i = 0;i < 256;i++) {count[i] = 0;count2[i] = 0;}
    for(uint64_t i = 0;i < newsize;i++) count[ref[newref2[i] + radix - 1]]++;
    for(uint32_t i = 1;i < 256;i++) count2[i] += count2[i - 1] + count[i - 1];
    for(uint64_t i = 0;i < newsize;i++) {
      newref[count2[ref[newref2[i] + radix - 1]]] = newref2[i];
      count2[ref[newref2[i] + radix - 1]]++;
    }
  }

  cerr << "convert sorted int aarr\n";
  convertSortedIntArray(ref,newref,newref2,newsize,sparse,sigma);
  return move(newref2);
}


#define ADJUSTSIZE(N) (((N) + 1) % (fbwtnum*oneEntryNumOfChar) == 0 ? (N) + 1 : (N) + (fbwtnum*oneEntryNumOfChar) + 1 - (((N) + 1) % (fbwtnum*oneEntryNumOfChar)))
#define OPENFILEADJUSTSIZE(PATH,SIZE,POINTER,TYPE) {			\
    FILE* fp;								\
    if((fp = fopen((char*)PATH, "rb")) == NULL ) {			\
      cerr << "fopen fail " << PATH << "\n";                            \
      exit(EXIT_FAILURE);						\
    }									\
    if(fseek(fp, 0, SEEK_END) == 0) {					\
      SIZE = ftell(fp)/sizeof(TYPE);					\
      rewind(fp);							\
      if(SIZE < 0) {							\
        cerr << "ftell fail " << PATH << "\n";                          \
	exit(EXIT_FAILURE);						\
      }									\
    } else {								\
      cerr << "fseek fail " << PATH << "\n";                            \
      exit(EXIT_FAILURE);						\
    }									\
    POINTER = (TYPE*)malloc(sizeof(TYPE)*ADJUSTSIZE(SIZE + 1));		\
    if(fread(POINTER, sizeof(TYPE), (size_t)SIZE, fp) != (size_t)SIZE) { \
      cerr << "fread fail " << PATH << "\n";                            \
      exit(EXIT_FAILURE);						\
    }									\
    fclose(fp);								\
  }

#define OPENFILE(PATH,SIZE,POINTER,TYPE) {				\
    FILE* fp;								\
    if((fp = fopen((char*)PATH, "rb")) == NULL ) {			\
      cerr << "fopen fail " << PATH << "\n";                            \
      exit(EXIT_FAILURE);						\
    }									\
    if(fseek(fp, 0, SEEK_END) == 0) {					\
      SIZE = ftell(fp)/sizeof(TYPE);					\
      rewind(fp);							\
      if(SIZE < 0) {							\
        cerr << "ftell fail " << PATH << "\n";                          \
	exit(EXIT_FAILURE);						\
      }									\
    } else {								\
      cerr << "fseek fail " << PATH << "\n";                            \
      exit(EXIT_FAILURE);						\
    }									\
    POINTER = (TYPE*)malloc(sizeof(TYPE)*SIZE);				\
    if(fread(POINTER, sizeof(TYPE), (size_t)SIZE, fp) != (size_t)SIZE) { \
      cerr << "fread fail " << PATH << "\n";                            \
      exit(EXIT_FAILURE);						\
    }									\
    fclose(fp);								\
  }

void outputIndex(string inputref,string outputdir,uint32_t kmer,uint32_t fbwtnum) {
  uint64_t reffilesize;
  unsigned char* ref;

  OPENFILEADJUSTSIZE(inputref.c_str(),reffilesize,ref,unsigned char) ;

  vector<pair<uint32_t,string>> refdescription;
  uint64_t nn = convert_fasta(ref,reffilesize,refdescription);
  uint64_t n = ADJUSTSIZE(nn);

  if(n % oneEntryNumOfChar != 0){
    cerr << "you should input bytes multiple of " << oneEntryNumOfChar << "\n";exit(1);
  }
  if(n %  fbwtnum != 0){
    cerr << "you should input bytes multiple of " << fbwtnum << "\n";exit(1);
  }

  for(uint64_t i = nn;i < n;i++){
    ref[i] = '`';
  }

  uint64_t newsize = n/fbwtnum;

  unique_ptr<uint32_t[]> SA(new uint32_t[newsize]{});
  unique_ptr<unsigned char[]> fbwt(new unsigned char[n]{});

  arrange(ref,n);
  
  uint32_t sigma = 0;
  unique_ptr<uint32_t[]> newref = convertRefbySparce(ref,n,fbwtnum,sigma);

  cerr << "refsize : " << n << "charsize : " << sigma << "\n";
  cerr << "constructing suffix array\n";
  if(sais_int((int*)newref.get(), (int*)SA.get(), (int)newsize,sigma + 1) != 0) {
    cerr << "sais fail\n";
    exit(1);
  }

  for(uint32_t i = 0;i < newsize;i++){
    SA[i] *= fbwtnum;
  }

  cerr << "converting reference\n";
  convertATGC(ref, n);

  cerr << "makeFBWT " << newsize << "\n";
  makeFBWT(SA,fbwt,ref,fbwtnum,newsize);

  cerr << "conpressing reference " << newsize << "\n";
  compressRef(ref, n);

  C_OCC c_occ(fbwt,SA,newsize,fbwtnum,kmer,ref,1,refdescription);

  c_occ.outputFile(outputdir.c_str());
  free(ref);
}


void memFromFastaMain(const char* queryfilepath,int32_t memLen,int32_t linearcomp,uint32_t intervalthreashold,C_OCC& c_occ,double& dtimes,unique_ptr<vector<Match>[]>& matches,uint64_t& countMEM,bool print,bool fourcolumn,int threadnum){
  FILE* fp;
  uint64_t queryfilesize;

  cerr << "open queryfile : " << queryfilepath << "\n";
  if((fp = fopen(queryfilepath, "rb")) == NULL ) {			
    cerr << "fopen fail " << queryfilepath << "\n";
    exit(1);						
  }									
  if(fseek(fp, 0, SEEK_END) == 0) {					
    queryfilesize = ftell(fp)/sizeof(unsigned char);
    rewind(fp);							
    if(queryfilesize < 0) {						
      cerr << "ftell fail " << queryfilepath << "\n";
      exit(EXIT_FAILURE);						
    }									
  } else {								
    cerr << "fseek fail " << queryfilepath << "\n";
    exit(EXIT_FAILURE);						
  }


  const uint32_t bufsize = 1024*1024;
  unsigned char buf[bufsize];
  uint32_t max_memorysize = 1024*1024*1024;
  unsigned char* query = (unsigned char*)malloc(sizeof(unsigned char)*max_memorysize);
  uint64_t queryLen = 0;

  bool isSkip = false;

  cerr << "start " << queryfilesize << "\n";
  string desc = "";
  for(uint64_t i = 0;i < queryfilesize;i+=bufsize){
    uint32_t size = fread(buf,sizeof(unsigned char),bufsize,fp);
    for(uint32_t j = 0;j < size;j++){
      if(isSkip){
	if(buf[j] == '\n') {isSkip = false; continue;}
        else desc += buf[j];
      }else{
	if(buf[j] == '>'){
	  if(queryLen > 0){
	    cerr << desc << ", size is " << queryLen << "\n";
	    double tstart = dtime();
	    c_occ.FindAllMEM(query, queryLen, queryLen,memLen,linearcomp,intervalthreashold,matches,threadnum);
	    double tend = dtime();
	    dtimes += tend - tstart;
	    queryLen = 0;
            desc = "";

	    if(print){
              for(int i = 0;i < threadnum;i++){
                c_occ.printOutput(matches[i],fourcolumn);
                countMEM += matches[i].size();
                matches[i].clear();
              }
	    }
	  }
	  isSkip = true; continue;
	}else if(buf[j] == '\n'){
	  continue;
	}else{
	  if(queryLen == max_memorysize){
	    max_memorysize *= 2;
	    unsigned char* _query = (unsigned char*)malloc(sizeof(unsigned char)*max_memorysize);
	    memcpy(_query,query,queryLen);
	    free(query);
	    query = _query;
	  }
	  query[queryLen] = convertATGC(buf[j]);
	  queryLen++;
	}
      }
    }
  }
  if(queryLen > 0){
    cerr << "querySize : " << queryLen << "\n";
    double tstart = dtime();
    c_occ.FindAllMEM(query, queryLen, queryLen,memLen,linearcomp,intervalthreashold,matches,threadnum);
    double tend = dtime();
    dtimes += tend - tstart;
    if(print){
      for(int i = 0;i < threadnum;i++){
        c_occ.printOutput(matches[i],fourcolumn);
        countMEM += matches[i].size();
        matches[i].clear();
      }
    }
  }

  if(!print){
    for(int i = 0;i < threadnum;i++){
      countMEM += matches[i].size();
    }
  }

  fclose(fp);							       
}


void MEM(const char* basename, vector<char*> queryFiles,int memlen,int linearcomp,int intervalthreashold,int sparseMult,bool print,bool fourcolumn,const int threadnum){
  if(sparseMult <= 0){
    cerr << "sparseMult must be over 0\n";
    exit(1);
  }

  C_OCC c_occ(basename,sparseMult,0);

  cerr << "start findMEMs\n";
  double dtimes = 0;
  unique_ptr<vector<Match>[]> matches(new vector<Match>[threadnum]);
  uint64_t countMEM = 0;
  for(auto& queryFile : queryFiles){
    cerr << "index is : " << basename << ", queryFile is " << queryFile <<  "\n";
    memFromFastaMain(queryFile,memlen,linearcomp,intervalthreashold,c_occ,dtimes,matches,countMEM,print,fourcolumn,threadnum);
  }
  
  uint64_t matchessize = 0;
  for(int i = 0;i < threadnum;i++){
    matchessize += matches[i].size();
  }

  cerr << "time : " << dtimes << ", memLen : " << memlen << ", hits : " << matchessize
       << ", intervalthreashold : " << intervalthreashold << ", indexFile : " << basename << ", queryFiles : " << queryFiles.size()
       << ", skip : " << sparseMult << ", fbwtnum : " << c_occ.fnum << ", maxMemorySize : " << maxMemorySize << ", hashsize : " 
       << c_occ.kmerSize << " #measure" << endl;
}

const string C_OCC::fbwtfname = "fbwt";
//vector<pair<int,int> > C_OCC::MEMCandidate::output;

void usage(string programname){
  cerr << programname << " [options] reference_file_path query_file_path_1 ... [query_file_path_n]\n";
  cerr << "Options related to making index\n";
  cerr << "-save path          save the index\n";
  cerr << "-kmer num           set the hash key length. default is 10\n";
  cerr << "-k num              set the size of FBWT\n";
  cerr << "\n";
  cerr << "options related to finding MEMs\n";
  cerr << "-l num              set the minimal MEM length. default is 50\n";
  cerr << "-load path          load the index\n";
  cerr << "-directcompth num   set the threshold of interval size to switch computing left maximal length to directly compare with reference. default is 10\n";
  cerr << "-intervallenth num  set the threshold of interval size to decide if more exact mathing in first step of algorithm. default is 10\n";
  cerr << "-skip num           sparsify the MEM-finding algorithm. default is possible maximum value\n";
  cerr << "-print [01]         print out the result when -print 1, not to do when -print 0. default is 1\n";
  cerr << "-fourcolumn [01]    print out result by fourcolumn align\n";
  cerr << "-threads num        execute num threads\n";
  cerr << "\n";
  cerr << "Make index\n";
  cerr << "fbwtmem_seq -kmer 8 -save path ref.fa\n";
  cerr << "Finding MEMs from index\n";
  cerr << "fbwtmem_seq -l 20 -load index query.fa\n";
  exit(1);
}


#ifdef TEST_INTV_OCC
int  TestIntv::backwardnum = 0;
#endif
int main(int argc, char *argv[]) {

#ifdef TEST_INTV_OCC
  for(int i = 0;i < TEST_INTV_OCC_MAX;i++){
    testintv[i].num = 0;
    testintv[i].firstmatch= 0;
    testintv[i].initmemcandidate = 0;
    testintv[i].updatememcandidate = 0;
    testintv[i].nobi = 0;
  }
#endif

  int memLen = 50;
  int kmer = 10;
  int directcompth = 10;
  int intervallenth = 10;
  int skip = 40;
  int fnum = 1;
  string load = "";
  string save = "";
  bool print = true;
  bool fourcolumn = false;
  int threadnum = 1;

  struct option long_options[] = {
    {"l",required_argument,0,0},
    {"kmer",required_argument,0,0},
    {"directcompth",required_argument,0,0},
    {"intervallenth",required_argument,0,0},
    {"skip",required_argument,0,0},
    {"save",required_argument,0,0},
    {"load",required_argument,0,0},
    {"k",required_argument,0,0},
    {"print",required_argument,0,0},
    {"fourcolumn",required_argument,0,0},
    {"threads",required_argument,0,0},
    {0,0,0,0}
  };

  while(true){
    int index;
    int retval = getopt_long_only(argc,argv,"",long_options,&index);
    if(retval < 0) break; // end options
    if(retval == '?'){
      cerr << "Invalid argument.\n";
      usage(argv[0]);
      exit(1);
    }
    switch(index){
    case 0:
      memLen = atoi(optarg);
      cerr << "set minimal MEM length : " << memLen << "\n";
      break;
    case 1:
      kmer = atoi(optarg);
      cerr << "set kmer : " << kmer << "\n";
      break;
    case 2:
      directcompth = atoi(optarg);
      cerr << "set directcompth : " << directcompth << "\n";
      break;
    case 3:
      intervallenth = atoi(optarg);
      cerr << "set intervallenth : " << intervallenth << "\n";
      break;
    case 4:
      skip = atoi(optarg);
      cerr << "set skip parameter : " << skip << "\n";
      break;
    case 5:
      save = optarg;
      cerr << "save index to " << save << "\n";
      break;
    case 6:
      load = optarg;
      cerr << "load index from " << load << "\n";
      break;
    case 7:
      fnum = atoi(optarg);
      cerr << "set size of FBWT " << fnum << "\n";
      break;
    case 8:
      print = atoi(optarg);
      cerr << "set print flag " << print << "\n";
      break;
    case 9:
      fourcolumn = atoi(optarg);
      cerr << "set fourcolumn flag " << fourcolumn << "\n";
      break;
    case 10:
      threadnum = atoi(optarg);
      cerr << "set threads " << threadnum << "\n";
      break;
    default:
      break;
    }
  }


  if(argc - optind < 1) {usage(argv[0]);}

  //void outputIndex(string inputref,string outputdir,uint32_t kmer,uint32_t fbwtnum) {
  if(!save.empty()){
    string inputref = string(argv[optind]);
    outputIndex(inputref,save,kmer,fnum);
  }
  if(!load.empty()){
    vector<char*> queryfiles;
    for(int i = optind;i < argc;i++){
      queryfiles.push_back(argv[i]);
    }
    MEM(load.c_str(),queryfiles,memLen,directcompth,intervallenth,skip,print,fourcolumn,threadnum);
  }

#ifdef TEST_INTV_OCC
  for(int i = 0;i < TEST_INTV_OCC_MAX;i++){
    if(testintv[i].num != 0)
      cerr << i << "," << testintv[i].num << ", " << testintv[i].firstmatch
	   << ", " << testintv[i].initmemcandidate << ", " << testintv[i].updatememcandidate
	   << ", " << testintv[i].nobi << ", " << Testintv[::backwardnum << ", ave, "
           << testintv[i].firstmatch/testintv[i].num << ", "
           << testintv[i].initmemcandidate/testintv[i].num << ", "
           << testintv[i].updatememcandidate/testintv[i].num << ", "
           << (double)testintv[i].nobi/testintv[i].num << ", "
           << " #frequency\n";
  }
#endif  

}

