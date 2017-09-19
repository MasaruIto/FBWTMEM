#ifndef __MAKEFBWT_H__
#define __MAKEFBWT_H__

#include <memory>

void makeFBWT(std::unique_ptr<int[]>& SA,std::unique_ptr<unsigned char[]>& fbwt,unsigned char* ref,int fbwtnum,int n);


#endif
