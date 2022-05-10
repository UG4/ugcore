// Felix Salfelder, 2016-2017, 2021
//
// This program is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the
// Free Software Foundation; either version 3, or (at your option) any
// later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, 51 Franklin Street - Suite 500, Boston, MA 02110-1335, USA.
//
//
// // contains code from boost/graph
// //=======================================================================
// // Copyright 1997, 1998, 1999, 2000 University of Notre Dame.
// // Authors: Andrew Lumsdaine, Lie-Quan Lee, Jeremy G. Siek
// //
// // Distributed under the Boost Software License, Version 1.0. (See
// // accompanying file LICENSE_1_0.txt or copy at
// // http://www.boost.org/LICENSE_1_0.txt)
// //=======================================================================
#ifndef TREEDEC_BUCKET_SORTER_HPP
#define TREEDEC_BUCKET_SORTER_HPP

// avoid including the the other one
#define BOOST_GRAPH_DETAIL_BUCKET_SORTER_HPP

#include <vector>
#include <cassert>
//#include <boost/limits.hpp>
#include "trace.h"

// the __FreeBSD__ condition is a wild guess.
#if defined(BOOST_CLANG) && (1 == BOOST_CLANG) && defined(__FreeBSD__)
#define ITERATOR_CONSTRUCTOR_WORKAROUND
#endif

namespace boost {

  template <class BucketType, class ValueType, class Bucket, 
            class ValueIndexMap>
  class bucket_sorter {
  public:
    typedef BucketType bucket_type;
    typedef ValueType value_type;
    typedef Bucket Bucket_type;
    typedef ValueIndexMap value_index_map;
    typedef typename std::vector<value_type>::size_type size_type;
    
    bucket_sorter(size_type _length, bucket_type _max_bucket, 
                  const Bucket& _bucket = Bucket(), 
                  const ValueIndexMap& _id = ValueIndexMap()) 
      : next(_length+_max_bucket, invalid_value()), 
        prev(_length+_max_bucket, invalid_value()),
        head(next.begin() + _length),
        tail(prev.begin() + _length),
        id_to_value(_length),
        bucket(_bucket), id(_id) {
          // trace2("created bs", _length, _max_bucket);
          assert(head==next.begin()+_length);
          assert(!next.size() || head[0]==invalid_value());


          auto l = _length;
          auto h = head;
          for(; l < _length + _max_bucket; ++l){
            *h = l;
            ++h;
          }
        }
    bucket_sorter(){untested();
    }

    void remove(const value_type& x) {
      size_type i = get(id, x);
      assert(i<size());
#if 0
      trace3("rm", x, i, bucket[x]);
      trace1("rm", head[bucket[x] ]);
      trace2("rm", prev[i], next[i]);
#endif
      auto next_node = next[i];
      auto prev_node = prev[i];
      assert(prev_node!=i);

      prev[next_node] = prev_node;
      next[prev_node] = next_node;

    } // remove

    void push_back(const value_type& x) {
      assert(x<next.size());
      id_to_value[get(id, x)] = x;
      (*this)[bucket[x]].push_back(x);
    }
    void push_front(const value_type& x) {
      // dont know how to do this right now, but i mean it.
      // assert( !is_in_bucket(x) );
      assert(x<next.size());
      id_to_value[get(id, x)] = x;
      (*this)[bucket[x]].push_front(x);
    }
    void push(const value_type& x) {
      return push_front(x);
    }

    // update_back?
    void update_front(const value_type& x) {
      // dont know how to do this right now, but i mean it.
      // assert( is_in_bucket(x) );
      remove(x);
      (*this)[bucket[x]].push(x);
    }
    void update(const value_type& x) {
      return update_front(x);
    }
    void update_back(const value_type& x) { itested();
      // dont know how to do this right now, but i mean it.
      // assert( is_in_bucket(x) );
      remove(x);
      (*this)[bucket[x]].push_back(x);
    }

    //  private: 
    //    with KCC, the nested stack class is having access problems
    //    despite the friend decl.
    static size_type invalid_value() {
      return (std::numeric_limits<size_type>::max)();
    }
    
    typedef typename std::vector<size_type>::iterator Iter;
    typedef typename std::vector<size_type>::const_iterator ConstIter;
    typedef typename std::vector<value_type>::iterator IndexValueMap;
    typedef typename std::vector<value_type>::const_iterator ConstIndexValueMap;
    
  public:

    template<class Iter_, class IndexValueMap_>
    class stack_ {
    public:
      typedef bucket_sorter base;
      typedef bucket_sorter::value_type value_type;
    public:
      class const_iterator{
        // bug? feature?
        const_iterator(){unreachable();}
      public:
        const_iterator(size_type t, stack_ const& s_)
           : s(s_), b(t) {}
        const_iterator(const const_iterator& p)
           : s(p.s), b(p.b) {}
        const_iterator(const const_iterator&& p)
           : s(p.s), b(p.b) {untested();}
        ~const_iterator(){}
      public:
        const_iterator& operator=(const const_iterator& o){
          assert(&s==&o.s); // how to compile time check?!
                            // or just fallback to pointer?
          b = o.b;
          return *this;
        }
        const_iterator& operator=(const const_iterator&& o){
          assert(&s==&o.s); // how to compile time check?!
                            // or just fallback to pointer?
          b = o.b;
          return *this;
        }
        value_type operator*() const{
          assert(b<s.size());
          return s.value[b];
        }
        const_iterator& operator++(){
          assert(b!=invalid_value());
          assert(b!=s.next[b]);
          b = s.next[b];
          return *this;
        }
        bool operator!=(const_iterator const& o){
          return o.b!=b;
        }
        bool operator==(const_iterator const& o)
        { return o.b==b; }
      private:
        stack_ const& s;
        size_type b;
      };
    public:
      stack_(const stack_& p)
        : bucket_id(p.bucket_id),
          head(p.head),
          next(p.next),
          prev(p.prev),
          tail(p.tail),
          value(p.value),
          id(p.id)
      {untested();
      }

      stack_& operator=(const stack_& p) {
        bucket_id = p.bucket_id;
        head = p.head;
        next = p.next;
        prev = p.prev;
        tail = p.tail;
        value = p.value;
        id = p.id;
        return *this;
      }
      // stack_(const stack_&&){untested();}
    public:
      stack_(bucket_type _bucket_id, Iter_ h, Iter_ n, Iter_ p, Iter_ t, IndexValueMap_ v,
            const ValueIndexMap& _id)
#ifdef ITERATOR_CONSTRUCTOR_WORKAROUND
        : bucket_id(_bucket_id), head(), next(), prev(), value(v), id(_id) { itested();
        head = h;
        next = n;
        prev = p;
#else
        : bucket_id(_bucket_id), head(h), next(n), prev(p), tail(t), value(v), id(_id) { // }
#endif
      }

      // Avoid using default arg for ValueIndexMap so that the default
      // constructor of the ValueIndexMap is not required if not used.
      stack_(bucket_type _bucket_id, Iter_ h, Iter_ n, Iter_ p, Iter_ t, IndexValueMap_ v)
#ifdef ITERATOR_CONSTRUCTOR_WORKAROUND
        : bucket_id(_bucket_id), head(), next(), prev(), value(v) { untested();
        head = h;
        next = n;
        prev = p;
#else
        : bucket_id(_bucket_id), head(h), next(n), prev(p), tail(t), value(v) { untested();
#endif
      }

      void push_front(const value_type& x) {
        const size_type new_head = get(id, x);
        assert(new_head < size());
        const size_type current = head[bucket_id];
        if(new_head == current){ untested();
//          assert(false);
        }
        if ( current != invalid_value() ){
          assert(current!=new_head);
          prev[current] = new_head;
        }else{ untested();
          tail[bucket_id] = new_head;
        }

        prev[new_head] = bucket_id + (head - next);
        next[new_head] = current;
        head[bucket_id] = new_head;
      }
      void push_back(const value_type& x) {
        const size_type new_tail = get(id, x);
        assert(new_tail < size());
        const size_type current = tail[bucket_id];

        if(new_tail == current){ untested();
//          assert(false);
        }
        if ( current != invalid_value() ){
          assert(current!=new_tail);
          next[current] = new_tail;
        }else{
          head[bucket_id] = new_tail;
        }

        prev[new_tail] = current;
        next[new_tail] = bucket_id + (head - next);
        tail[bucket_id] = new_tail;
      }
      void push(const value_type& x) {
        return push_front(x);
      }
      void pop() {
        assert(!empty());
        assert(bucket_id<size());
        size_type current = head[bucket_id];
        size_type next_node = next[current];
        head[bucket_id] = next_node;
        // prev[bucket_id] = bucket_id + (head - prev);
        if ( next_node != invalid_value() ){
//          assert(next_node != prev[current]);
          prev[next_node] = bucket_id + (head - next);
        }else{ untested();
          unreachable();
        }
        // assert(next_node != prev[current]);
      }
      value_type const& top() const { return value[ head[bucket_id] ]; }
      value_type& top() { return value[ head[bucket_id] ]; }
      value_type const& front() const { return value[ head[bucket_id] ]; }
      value_type& front() { return value[ head[bucket_id] ]; }
      bool empty() const {
        return begin() == end();
      }
    public: // iterator access
      const_iterator begin() const{
        return const_iterator(head[bucket_id], *this);
      }
      // BUG: template override in degree.hpp does not match (why?)
      const_iterator rbegin() const{ untested();
        return const_iterator(head[bucket_id], *this);
      }
      const_iterator end() const{
        return const_iterator(size()+bucket_id, *this);
      }
    public: // debug
      size_t size()const{return head-next;}
    private:
      bucket_type bucket_id;
      Iter_ head;
      Iter_ next;
      Iter_ prev;
      Iter_ tail;
      IndexValueMap_ value;
      ValueIndexMap id;
    }; // stack

    typedef stack_<Iter, IndexValueMap> stack;
    typedef stack_<ConstIter, ConstIndexValueMap> const_stack;
    
    const_stack operator[](const bucket_type& i) const{
      assert(i < next.size());
      return const_stack(i, head, next.begin(), prev.begin(), tail,
                   id_to_value.begin(), id);
    }
    stack operator[](const bucket_type& i) {
      assert(i < next.size());
      return stack(i, head, next.begin(), prev.begin(), tail,
                   id_to_value.begin(), id);
    }
    // number of buckets.
    size_t size() const{ return next.size() - (head-next.begin()); }
  protected:
    std::vector<size_type>   next;
    std::vector<size_type>   prev;
    typename std::vector<size_type>::iterator head;
    typename std::vector<size_type>::iterator tail;
    std::vector<value_type>  id_to_value;
    Bucket bucket;
    ValueIndexMap id;
  }; // stack_
  
}

#endif

// vim:ts=8:sw=2:et
