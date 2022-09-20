/*
 * Authors: Gustavo V. Barroso 
 * Created: 19/09/2022
 * Last modified: 19/09/2022
 *
 */

 
#ifndef _HASHTABLE_
#define _HASHTABLE_

#include <random>
#include <iostream>
#include <utility>
#include <chrono>
#include <array>
#include <functional>

#include "HashNode.h"

template<typename ValueType>
class HashTable {
    
private:
  HashNode<ValueType>** table_;
    
  unsigned long long tableSize_;
  unsigned long long largePrime_; 
  
  unsigned long long a_;
  unsigned long long b_;
    
public:
  HashTable(unsigned long long tableSize, unsigned long long universeSize): 
  table_(),
  tableSize_(tableSize),
  largePrime_(0),
  a_(0),
  b_(0)
  {
    // constructs zero initialized hash table of size tableSize_
    table_ = new HashNode<ValueType>* [tableSize_]();
    largePrime_ = findLargePrime_(universeSize);
    
    std::mt19937 gen;  
    std::array<int, std::mt19937::state_size> seedData;
    unsigned sem = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine re(sem);
    std::generate_n(seedData.data(), seedData.size(), std::ref(re));
    std::seed_seq seq(std::begin(seedData), std::end(seedData));
    gen.seed(seq);  
    
    std::uniform_int_distribution<unsigned long long> dist(0, largePrime_ - 1); 
  
    a_ = dist(gen);
    b_ = dist(gen);    
  }
  
  HashTable(): 
  table_(nullptr),
  tableSize_(0),
  largePrime_(0),
  a_(0),
  b_(0)
  { } 
  
  HashTable(HashTable&& other):
  table_(nullptr),
  tableSize_(0),
  largePrime_(0),
  a_(0),
  b_(0)
  {
    table_ = other.table_;
    tableSize_ = other.tableSize_;
    largePrime_ = other.largePrime_;
    a_ = other.a_;
    b_ = other.b_;
    
    other.table_ = nullptr;
    other.tableSize_ = 0;
    other.largePrime_ = 0;
    other.a_ = 0;
    other.b_ = 0;
  }
  
  HashTable(HashTable& other):
  table_(nullptr),
  tableSize_(0),
  largePrime_(0),
  a_(0),
  b_(0)
  {
    // deep copy 
    table_ = new HashNode<ValueType>* [other.tableSize_]; 
    for(unsigned long long i = 0; i < other.tableSize_; ++i)
    { 
      if(other.table_[i] != NULL)
        table_[i] = new HashNode<ValueType>(*(other.table_[i]));
    }
    
    tableSize_ = other.tableSize_;
    largePrime_ = other.largePrime_;
    a_ = other.a_;
    b_ = other.b_;
  }
  
  HashTable& operator=(HashTable& other)
  {
    if(this == &other)
      return *this;
    else
    {
      delete [] table_;
      // deep copy 
      table_ = new HashNode<ValueType>* [other.tableSize_]; 
      for(unsigned long long i = 0; i < other.tableSize_; ++i)
      { 
        if(other.table_[i] != NULL)
          table_[i] = new HashNode<ValueType>(*(other.table_[i]));
      }
    
      tableSize_ = other.tableSize_;
      largePrime_ = other.largePrime_;
      a_ = other.a_;
      b_ = other.b_;
    }
    
    return *this;
  }
  
  ~HashTable() 
  {
    // destroys all buckets one by one
    for(unsigned long long i = 0; i < tableSize_; ++i)
    {
      HashNode<ValueType>* entry = table_[i];
     
      while(entry != NULL) 
      {
        HashNode<ValueType>* prev = entry;
        entry = entry->getNext();
      
        delete prev;
      }
      
      table_[i] = NULL;
    }

    delete [] table_;
  }

public:
  void reset(ValueType val); 
    
  void init(const std::vector<unsigned long long>& chrSnps, ValueType val); 
  
  void setValue(unsigned long long key, ValueType val); 
    
  void insert(unsigned long long key, ValueType val);

  void remove(unsigned long long key);
  
  ValueType getValue(unsigned long long key) const;
  
  unsigned long long getTableSize() const
  {
    return tableSize_;
  }
  
  HashNode<ValueType>**& getTable()
  {
    return table_;
  }
  
private:
  unsigned long long preHash_(unsigned long long key) const;
  
  unsigned long long findLargePrime_(unsigned long long tableSize);

};

#endif
