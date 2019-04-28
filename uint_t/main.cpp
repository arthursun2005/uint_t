//
//  main.cpp
//  uint_t
//
//  Created by Arthur Sun on 4/27/19.
//  Copyright © 2019 Arthur Sun. All rights reserved.
//

#include <iostream>
#include "uint_t.h"

int main(int argc, const char * argv[]) {
    uint_t j = 1;
    
    unsigned long a = clock();
    
    for(int i = 2; i <= 6000; ++i) {
        j *= i;
    }
    
    unsigned long b = clock();

    printf("%lu\n", b - a);

    //printf("%lu\n", j.uint2());
    //printf("%s\n", j.toString(10).c_str());
    //printf("%lu\n", 1L << 47);
    //printf("%llu\n", 11111111llu * 11111111);
    return 0;
}
