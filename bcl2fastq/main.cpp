//
//  main.cpp
//  bcl2fastq
//
//  Created by Nikola Jovanovic on 7/11/18.
//  Copyright © 2018 Nikola Jovanovic. All rights reserved.
//
#include "BCLConverter.h"

int main(int argc, const char * argv[]) {
    if(!argv[1]){
        cerr << "No project path specified!\n";
        exit(1);
    };
    BCLConverter bcc = BCLConverter();
    bcc.setProjectPath(argv[1]);
    bcc.parseRunInfo();
    bcc.parseSampleSheet();
    bcc.checkBCLDirectory();
    bcc.convert();
    exit(0);
};

//TODO
//add check if tile fetched has the same tile_num, add exception

//TODO
//add outputdir

//TODO
//process only lanes that exist in dirs, ignore all N reads

//TODO
//delete already existing files on first execution

//TODO
//check if all of the files are closed - fastq and cbcls

//TODO
//allow for random readstructure format - vector instead of array

//TODO compress output

//TODO bgzf

//TODO change returnmessage with exception or exit and cerr

//TODO determine max_reads_in_mem

//add params

//add help

//non-pf from tileinfo for filters

//TODO try with no samplesheet

//TODO potentially change pointers for array

//TODO logging

//TODO far far in the future adapter trimming

//todo umis

