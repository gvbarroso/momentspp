/*
 * Authors: Gustavo V. Barroso 
 * Created: 19/09/2022
 * Last modified: 19/09/2022
 *
 * Modified after:
 * https://medium.com/@aozturk/simple-hash-map-hash-table-implementation-in-c-931965904250
 */


#include "HashTable.h"

template<typename ValueType>
void HashTable<ValueType>::reset(ValueType val)
{
  for(size_t i = 0; i < tableSize_; ++i)
  {
    HashNode<ValueType>* entry = table_[i];
    
    while(entry != NULL) 
    {
      entry->setValue(val);
      entry = entry->getNext();
    }
  }
}

template<typename ValueType>
void HashTable<ValueType>::init(const std::vector<unsigned long long>& keys, ValueType val)
{
  for(auto itKey = std::begin(keys); itKey != std::end(keys); ++itKey)
    insert(*itKey, val);
}

template<typename ValueType>
void HashTable<ValueType>::setValue(unsigned long long key, ValueType val) 
{
  unsigned long long hashValue = preHash_(key);
  
  HashNode<ValueType>* entry = table_[hashValue];

  while(entry->getKey() != key)
    entry = entry->getNext();

  entry->setValue(val);
}

template<typename ValueType>
void HashTable<ValueType>::insert(unsigned long long key, ValueType value) 
{
  unsigned long long hashValue = preHash_(key);
  
  HashNode<ValueType>* prev = NULL;
  HashNode<ValueType>* entry = table_[hashValue];

  while(entry != NULL && entry->getKey() != key)
  {
    prev = entry;
    entry = entry->getNext();
  }

  if(entry == NULL)
  {
    entry = new HashNode<ValueType>(key, value);
        
    if(prev == NULL) // insert as first bucket 
      table_[hashValue] = entry;
    
    else
      prev->setNext(entry);
  }
  
  else // just update the value
    entry->setValue(value);
}

template<typename ValueType>
void HashTable<ValueType>::remove(unsigned long long key) 
{
  unsigned long long hashValue = preHash_(key);
  
  HashNode<ValueType>* prev = NULL;
  HashNode<ValueType>* entry = table_[hashValue];

  while(entry != NULL && entry->getKey() != key) 
  {
    prev = entry;
    entry = entry->getNext();
  }

  if(entry == NULL) // key not found
    return;
  
  else
  {
    if(prev == NULL) // remove first bucket of the list
      table_[hashValue] = entry->getNext();
    
    else
      prev->setNext(entry->getNext());
      
    delete entry;
  }
}

template<typename ValueType>
ValueType HashTable<ValueType>::getValue(unsigned long long key) const
{
  unsigned long long hashValue = preHash_(key);
  HashNode<ValueType>* entry = table_[hashValue];
  
  while(entry->getKey() != key)
    entry = entry->getNext();

  return entry->getValue();
}

template<typename ValueType>
unsigned long long HashTable<ValueType>::findLargePrime_(unsigned long long number)
{
  auto isPrime = [=](unsigned long long n)  
  {  
    // corner cases  
    if(n <= 1)
      return false; 
    
    if(n <= 3)
      return true;  
    
    // this is checked so that we can skip middle five numbers in below loop  
    if(n % 2 == 0 || n % 3 == 0)
      return false;  
    
    for(unsigned long long i = 5; i * i <= n; i = i + 6)  
    {
      if(n % i == 0 || n % (i + 2) == 0)  
        return false;  
    }
    
    return true;  
  }; 

  unsigned long long prime = number; 
  bool found = false; 
  
  while(!found) 
  { 
    ++prime; 
    
    if(isPrime(prime)) 
      found = true; 
  }
  
  return prime; 
}

template<typename ValueType>
unsigned long long HashTable<ValueType>::preHash_(unsigned long long key) const
{
  unsigned long long v = (((a_ * key) + b_) % largePrime_) % tableSize_;
  //std::cout << key << "," << a_ << "," << b_ << "," << largePrime_ << "," << tableSize_ << "," << v << "\n";
  return v;
}

template class HashTable<size_t>;
template class HashTable<double>;
template class HashTable<std::pair<size_t, size_t>>;
