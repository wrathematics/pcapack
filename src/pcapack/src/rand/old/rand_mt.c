/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

// Copyright 2013, Schmidt and Heckendorf

#include <time.h>
#include <math.h>
#include <stdlib.h>

#include "rand.h"


// ###################### //
// Mersenne Twister 19937 //
// ###################### //

unsigned int MT[MT_SIZE];
int mt_index = 0;
unsigned int mt_initialized = 0;


void mt_init_gen(const unsigned int seed)
{
    int i;
    
    MT[0]=seed;
    
    for(i=1;i<MT_SIZE;i++)
        MT[i]=0xFFFFFFFF&(0x6C078965 * (MT[i-1] ^ ((MT[i-1])>>30)) + i);
}

void mt_init(unsigned int seed)
{
    mt_init_gen(seed);
    mt_initialized = 1;
}


void mt_check()
{
    unsigned int seed;
    
    if (mt_initialized == 0)
    {
        seed = (unsigned int) time(NULL);
        mt_init_gen(seed);
        mt_initialized = 1;
    }
}

void mt_gen()
{
    int i;
    unsigned int y;
    
    for(i=0;i<MT_SIZE;i++)
    {
        y=(MT[i]&0x80000000) + (MT[(i+1)%MT_SIZE] & 0x7FFFFFFF);
        MT[i]=MT[(i+397)%MT_SIZE] ^ (y>>1);
        
        if(y%2!=0)
            MT[i]^=0x9908B0DF;
    }
}

unsigned int mt_extract()
{
    unsigned int y;
    
    if(mt_index==0)
        mt_gen();
    
    y=MT[mt_index];
    y^=y>>11;
    y^=(y<<7)&0x9D2C5680;
    y^=(y<<15)&0xEFC60000;
    y^=y>>18;
    
    mt_index=(mt_index+1)%MT_SIZE;
    
    return y;
}

