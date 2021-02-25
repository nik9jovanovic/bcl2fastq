//
//  BCLHeaderReader.h
//  bcl2fastq
//
//  Created by Nikola Jovanovic on 8/3/18.
//  Copyright © 2018 Nikola Jovanovic. All rights reserved.
//
#pragma once

#include "Semaphore.h"
#include <string>
#include <stdio.h>
#include <zlib.h>
using namespace std;

//Container structure holding the header information about each tile - only CBCLs have this
struct TileInfo{
    unsigned int tile_num;
    unsigned int num_clusters;
    unsigned int uncomp_size;
    unsigned int comp_size;
    bool non_pf;
    
    TileInfo(){
        tile_num = 0;
        num_clusters = 0;
        uncomp_size = 0;
        comp_size = 0;
        non_pf = 0;
    };
    
    TileInfo(unsigned int tn, unsigned int nc, unsigned int us, unsigned int cs, bool np){
        tile_num = tn;
        num_clusters = nc;
        uncomp_size = us;
        comp_size = cs;
        non_pf = np;
    };
};

//Class for reading from CBCL files
class CBCLReader{

    unsigned int header_size;
    unsigned int* bins;
    unsigned int num_tiles;
    TileInfo* tiles;
    
    unsigned char* buffer;
    
public:
    CBCLReader();
    CBCLReader(const unsigned int hs, unsigned int* bin, const unsigned int nt, TileInfo* t, unsigned char* buf);
    CBCLReader(const CBCLReader &c);
    ~CBCLReader();
    
    unsigned char* getBuffer() const;
    unsigned int getHeaderSize() const;
    TileInfo* getTiles() const;
    unsigned int getNumTiles() const;
    
    unsigned int* getBins() const;
    unsigned int getUInt(FILE* f) const;
    bool checkTile(unsigned int tile, unsigned int tile_num) const;
    void clearBuffer();
    
    unsigned int parseHeader(const string path, Semaphore* sem);
    unsigned int readCBCL(const string path, const unsigned char tile, const unsigned int startpos, const unsigned int number_of_reads, Semaphore* sem);
};
