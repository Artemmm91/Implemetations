#include <iostream>
#include <memory>
#include <cassert>

#include "smart_pointers.h"

template<typename T>
struct MyAlloc {
  using value_type = T;

  MyAlloc() = default;

  template<typename U>
  MyAlloc(const MyAlloc<U>&) {}
    
  T* allocate(size_t n) {
    return (T*) ::operator new(n * sizeof(T));
  }

  void deallocate(T* p, size_t n) {
    ::operator delete((void*)p);
  }
};


struct A {
  int x;
  static int deleted;
  virtual ~A() {
    ++deleted;
  }
};
int A::deleted = 0;

int main() {
    {
        SharedPtr<A> a = allocateShared<A, MyAlloc<A>>(MyAlloc<A>());
    }
    /*auto alloc = MyAlloc<A>();
    using ControlBlockType = ControlBlock<A, MyAlloc<A>, EmptyDeleter<A, MyAlloc<A>>>;
    using AllocChar = typename std::allocator_traits<MyAlloc<A>>::template rebind_alloc<char>;
    AllocChar alloc_data = alloc;
    
    char* data_ptr = alloc_data.allocate(sizeof(ControlBlockType) + sizeof(A));
    auto counter = reinterpret_cast<BaseControlBlock*>(data_ptr);
    auto ptr = data_ptr + sizeof(BaseControlBlock*);
    
    A* aptr = new (ptr) A();
    new (counter) ControlBlockType(EmptyDeleter<A, MyAlloc<A>>(alloc), alloc);
    
    aptr->~A();*/
    assert(A::deleted == 1);
}
