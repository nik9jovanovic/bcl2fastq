//
//  Semaphore.h
//  bcl2fastq
//
//  Created by Nikola Jovanovic on 8/10/18.
//  Copyright © 2018 Nikola Jovanovic. All rights reserved.
//
#pragma once

#include <mutex>
#include <condition_variable>

//Concurrency class for multithreading
class Semaphore{
    std::mutex mutex;
    std::condition_variable condition;
    unsigned int count;
    
public:
    Semaphore();
    Semaphore(const unsigned int c);
    Semaphore(const Semaphore &s);
    
    unsigned int getCount() const;
    void setCount(const unsigned int c);
    
    void wait();
    void notify();
};
