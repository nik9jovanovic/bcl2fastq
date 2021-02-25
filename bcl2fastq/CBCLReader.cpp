//
//  BCLHeaderReader.cpp
//  bcl2fastq
//
//  Created by Nikola Jovanovic on 8/3/18.
//  Copyright © 2018 Nikola Jovanovic. All rights reserved.
//
#include "CBCLReader.h"

CBCLReader::CBCLReader(){
    header_size = 0;
    
    bins = NULL;
    tiles = NULL;
    buffer = NULL;
};

CBCLReader::CBCLReader(const unsigned int hs, unsigned int* bin, const unsigned int nt, TileInfo* t, unsigned char* buf){
    header_size = hs;
    bins = bin;
    num_tiles = nt;
    tiles = t;
    buffer = buf;
};

CBCLReader::CBCLReader(const CBCLReader &c){
    header_size = c.getHeaderSize();
    bins = c.getBins();
    num_tiles = c.getNumTiles();
    tiles = c.getTiles();
    buffer = c.getBuffer();
};

CBCLReader::~CBCLReader(){
    delete [] bins;
    delete [] tiles;
    delete [] buffer;
};

unsigned char* CBCLReader::getBuffer() const{
    return buffer;
};

unsigned int CBCLReader::getHeaderSize() const{
    return header_size;
};

TileInfo* CBCLReader::getTiles() const{
    return tiles;
};

unsigned int CBCLReader::getNumTiles() const{
    return num_tiles;
};

unsigned int* CBCLReader::getBins() const{
    return bins;
};

unsigned int CBCLReader::getUInt(FILE* f) const{
    return fgetc(f) + (fgetc(f) << 8) + (fgetc(f) << 16) + (fgetc(f) << 24);
};

bool CBCLReader::checkTile(unsigned int tile, unsigned int tn) const{
    if(tiles[tile].tile_num == tn) return 1;
    return 0;
};

void CBCLReader::clearBuffer(){
    if(buffer) delete [] buffer;
    buffer = NULL;
};

unsigned int CBCLReader::parseHeader(const string path, Semaphore* sem){
    FILE* cbcl = fopen(path.c_str(), "r");
    if(cbcl == NULL){
        sem->notify();
        return 1;
    };
    
    unsigned int version =  fgetc(cbcl) + (fgetc(cbcl) << 8);
    header_size = getUInt(cbcl);
    unsigned int bit_per_basecall = fgetc(cbcl);
    unsigned int bit_per_qscore = fgetc(cbcl);
    if((version != 1) || (bit_per_basecall != 2) || (bit_per_qscore != 2)){
        //the only supported combination right now
        sem->notify();
        return 1;
    };
    unsigned int num_bins = getUInt(cbcl);
    bins = new unsigned int [num_bins];
    for(int bin = 0; bin < num_bins; bin++){
        unsigned int num = getUInt(cbcl);
        bins[num] = getUInt(cbcl);
    };
    num_tiles = getUInt(cbcl);
    tiles = new TileInfo [num_tiles];
    for(int tile = 0; tile < num_tiles; tile++){
        tiles[tile].tile_num = getUInt(cbcl);
        tiles[tile].num_clusters = getUInt(cbcl);
        tiles[tile].uncomp_size = getUInt(cbcl);
        tiles[tile].comp_size = getUInt(cbcl);
        tiles[tile].non_pf = fgetc(cbcl);
    };
    
    fclose(cbcl);
    sem->notify();
    return 0;
};

unsigned int CBCLReader::readCBCL(const string path, const unsigned char tile, const unsigned int startpos, const unsigned int number_of_reads, Semaphore* sem){
    if(buffer != NULL){
        delete [] buffer;
        buffer = NULL;
    };
    if(header_size == 0) if(this->parseHeader(path, sem)) return 1;
        
    FILE* f = fopen(path.c_str(), "r");
    if(f == NULL){
        sem->notify();
        return 1;
    };
    
    char* buffer_in = new char [tiles[tile].comp_size];
    char* buffer_out = new char [tiles[tile].uncomp_size];
    
    fseek(f,header_size, SEEK_SET);
    unsigned int slide_pos = startpos;
    unsigned int tmp_tile = 0;
    for(;(tmp_tile < num_tiles) && (slide_pos > tiles[tmp_tile].comp_size);tmp_tile++){
        slide_pos -= tiles[tmp_tile].uncomp_size;
        fseek(f,tiles[tmp_tile].comp_size, SEEK_CUR);
    };
    size_t result = fread(buffer_in,1,tiles[tile].comp_size,f);
    if(result != tiles[tile].comp_size) return 1;
    fclose(f);
    
    //gzip decompress the block
    z_stream infstream;
    infstream.zalloc = Z_NULL;
    infstream.zfree = Z_NULL;
    infstream.opaque = Z_NULL;
    infstream.avail_in = tiles[tile].comp_size;
    infstream.next_in = (Bytef *)buffer_in;
    infstream.avail_out = tiles[tile].uncomp_size;
    infstream.next_out = (Bytef *)buffer_out;
    inflateInit2(&infstream, 16+MAX_WBITS);
    inflate(&infstream, Z_NO_FLUSH);
    inflateEnd(&infstream);
    
    //determine how many reads will be processed
    unsigned int num_reads = number_of_reads;
    if(num_reads == 0) num_reads = tiles[tile].uncomp_size;
    if((slide_pos + num_reads) > tiles[tile].uncomp_size) num_reads = tiles[tile].uncomp_size - slide_pos;
    
    //fetch only the amount of reads specified
    delete [] buffer_in;
    buffer = new unsigned char [num_reads];
    for(int cpy = startpos; cpy < (startpos + num_reads); cpy++) buffer[cpy - startpos] = buffer_out[cpy];
    delete [] buffer_out;
    
    sem->notify();
    return 0;
};
