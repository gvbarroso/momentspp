/*
 * Authors: Gustavo V. Barroso 
 * Created: 19/09/2022
 * Last modified: 19/09/2022
 *
 */


#ifndef _HASHNODE_
#define _HASHNODE_

#include <vector>
#include <string>
#include <algorithm>
#include <utility>

// hashing with chaining
template<typename ValueType>
class HashNode {
    
private:
    std::string key_; // moment name
    
    ValueType value_; // moment value (double)

    std::vector<size_t> popIndices_; // that appear in moment name
    
    // next bucket with the same prehash
    HashNode<ValueType>* next_;
    
public:
  HashNode(unsigned long long key, ValueType value):
  key_(key),
  value_(value),
  popIndices_(0),
  next_(nullptr) 
  { }
  
  HashNode(const HashNode<ValueType>& other):
  key_(0),
  value_(),
  popIndices_(0),
  next_(nullptr) 
  {
    key_ = other.key_;
    value_ = other.value_;
    popIndices_ = other.popIndices_;

    if(other.next_ != nullptr)
      next_ = new HashNode<ValueType>(*other.next_); // deep copy
  }

  const std::string& getKey() const
  {
    return key_;
  }

  ValueType getValue() const
  {
    return value_;
  }

  const std::vector<size_t> getPopIndices_() const
  {
    return popIndices_;
  }
    
  void setValue(ValueType value)
  {
    HashNode::value_ = value;
  }

  HashNode<ValueType>* getNext() const
  {
    return next_;
  }

  void setNext(HashNode<ValueType>* next)
  {
    HashNode::next_ = next;
  }
  
};

#endif
