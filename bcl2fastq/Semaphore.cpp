//
//  Semaphore.cpp
//  bcl2fastq
//
//  Created by Nikola Jovanovic on 8/10/18.
//  Copyright © 2018 Nikola Jovanovic. All rights reserved.
//
#include "Semaphore.h"

Semaphore::Semaphore(){
    count = 0;
};

Semaphore::Semaphore(const unsigned int c){
    count = c;
};

Semaphore::Semaphore(const Semaphore &s){
    count = s.getCount();
};

void Semaphore::notify() {
    std::unique_lock<decltype(mutex)> lock(mutex);
    ++count;
    condition.notify_one();
};

void Semaphore::wait(){
    std::unique_lock<decltype(mutex)> lock(mutex);
    while(!count) 
        condition.wait(lock);
    --count;
};

void Semaphore::setCount(const unsigned int c){
    count = c;
};

unsigned int Semaphore::getCount() const{
    return count;
};
