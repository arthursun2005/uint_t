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
    uint_t j = uint_t(111111L) * 111111;
    printf("%lu\n", uint_t::base * j[1] + j[0]);
    printf("%lu\n", 1L << 63);
    printf("%llu\n", 111111llu * 111111);
    return 0;
}
