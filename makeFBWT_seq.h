#ifndef __MAKEFBWT_H__
#define __MAKEFBWT_H__

#include <memory>
#include <stdint.h>

void makeFBWT(std::unique_ptr<uint32_t[]>& SA,std::unique_ptr<unsigned char[]>& fbwt,unsigned char* ref,uint32_t fbwtnum,uint32_t n);


#endif
