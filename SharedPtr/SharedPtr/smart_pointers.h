#include <iostream>

template <typename T, typename Alloc>
struct EmptyDeleter {
    Alloc alloc;
    EmptyDeleter(Alloc alloc): alloc(alloc) {}
    void operator()(T* ptr) {
        std::allocator_traits<Alloc>::destroy(alloc, ptr);
    }
};

template <typename T, typename Alloc, typename Deleter>
struct ControlData;

struct BaseControlBlock {
    virtual ~BaseControlBlock() {}
    virtual void increment_shared_count() = 0;
    virtual void decrement_shared_count() = 0;
    virtual void increment_weak_count() = 0;
    virtual void decrement_weak_count() = 0;
    virtual size_t get_shared_count() = 0;
    virtual size_t get_weak_count() = 0;
    virtual void operator()(void*) = 0;
    virtual void destroy() = 0;
};

template <typename T, typename Alloc, typename Deleter>
struct ControlBlock: public BaseControlBlock {
    using ControlBlockType = ControlBlock<T, Alloc, Deleter>;
    
    size_t shared_count;
    size_t weak_count;
    Deleter deleter;
    Alloc alloc;
    
    ControlBlock(Deleter d = Deleter(), Alloc alloc = Alloc()): shared_count(1), weak_count(0),
    deleter(d), alloc(alloc) {}
    
    void increment_shared_count() override {
        ++shared_count;
    }
    
    void decrement_shared_count() override {
        --shared_count;
    }
    
    void increment_weak_count() override {
        ++weak_count;
    }
    
    void decrement_weak_count() override {
        --weak_count;
    }
    
    size_t get_shared_count() override {
        return shared_count;
    }
    
    size_t get_weak_count() override {
        return weak_count;
    }
    
    void operator()(void* ptr) override {
        deleter(reinterpret_cast<T*>(ptr));
    }
    
    void destroy_all() {
        using ControlDataType = ControlData<T, Alloc, EmptyDeleter<T, Alloc>>;
        using AllocData = typename std::allocator_traits<Alloc>::template rebind_alloc<ControlDataType>;
        AllocData alloc_data = alloc;
        this->~ControlBlock();
        std::allocator_traits<AllocData>::deallocate(alloc_data, reinterpret_cast<ControlDataType*>(this), 1);
    }
    
    void destroy() override {
        if (std::is_same<Deleter, EmptyDeleter<T, Alloc>>::value) {
            destroy_all();
            return;
        }
        using AllocControlBlock = typename std::allocator_traits<Alloc>::template rebind_alloc<ControlBlockType>;
        AllocControlBlock alloc_control_block = alloc;
        this->~ControlBlock();
        std::allocator_traits<AllocControlBlock>::deallocate(alloc_control_block, this, 1);
    }
    
    ~ControlBlock() = default;
};

template <typename T, typename Alloc, typename Deleter>
struct ControlData {
    ControlBlock<T, Alloc, Deleter> control_block;
    T object;
};

template <typename T>
class SharedPtr {
    template <typename U>
    friend class SharedPtr;
    template <typename U>
    friend class WeakPtr;
private:
    BaseControlBlock* counter;
    T* ptr;
    
    struct construct_tag {};
    struct construct_control_tag {};
    
    template <typename Alloc, typename... Args>
    SharedPtr(typename SharedPtr<T>::construct_tag, Alloc, Args&& ...args);
    
    SharedPtr(typename SharedPtr<T>::construct_control_tag, T* ptr, BaseControlBlock* counter);
    
public:
    template <typename U, typename Deleter, typename Alloc>
    SharedPtr(U* ptr, Deleter d, Alloc alloc);
    
    SharedPtr();

    template <typename U>
    SharedPtr(U* ptr);

    SharedPtr(const SharedPtr& other);

    template <typename U>
    SharedPtr(const SharedPtr<U>& other);

    template <typename U>
    SharedPtr(SharedPtr<U>&& other);

    SharedPtr(SharedPtr&& other);

    template <typename U>
    SharedPtr<T>& operator=(U&& other);
    
    template <typename U, typename Deleter>
    SharedPtr(U* ptr, Deleter d);

    ~SharedPtr<T>();

    size_t use_count() const;

    void reset();
    template <typename U>
    void reset(U*);

    void swap(SharedPtr<T>&);

    T* get() const;
    T& operator*() const;
    T* operator->() const;
    explicit operator bool() const;

    template <typename U, typename... Args>
    friend SharedPtr<U> makeShared(Args&&... args);

    template <typename U, typename Alloc, typename... Args>
    friend SharedPtr<U> allocateShared(const Alloc&, Args&&...);
    
};


template <typename T>
class WeakPtr {
    template <typename U>
    friend class WeakPtr;
    template <typename U>
    friend class SharedPtr;
private:
    BaseControlBlock* counter;
    T* ptr;
    
    
public:
    WeakPtr();
    
    WeakPtr(const WeakPtr& other);
    
    template <typename U>
    WeakPtr(const WeakPtr<U>& other);
    
    template <typename U>
    WeakPtr(WeakPtr<U>&& other);
    
    WeakPtr(WeakPtr<T>&& other);
    
    template <typename U>
    WeakPtr(const SharedPtr<U>& shared);
    
    ~WeakPtr<T>();
    
    template <typename U>
    WeakPtr<T>& operator=(const SharedPtr<U>&);
    template <typename Other>
    WeakPtr<T>& operator=(Other&& other);
    
    size_t use_count() const;
    
    void swap(WeakPtr<T>&);
    
    T* get() const;
    T& operator*() const;
    T* operator->() const;
    explicit operator bool() const;
    
    SharedPtr<T> lock() const;
    bool expired() const;
};

template <typename T>
SharedPtr<T>::SharedPtr(): counter(nullptr), ptr(nullptr) {}

template <typename T>
template <typename U>
SharedPtr<T>::SharedPtr(U* ptr): SharedPtr(ptr, std::default_delete<U>(), std::allocator<U>()) {}

template <typename T>
SharedPtr<T>::SharedPtr(const SharedPtr& other): counter(other.counter), ptr(other.ptr) {
    if (counter == nullptr) {
        return;
    }
    counter->increment_shared_count();
}

template <typename T>
template <typename U>
SharedPtr<T>::SharedPtr(const SharedPtr<U>& other): counter(other.counter), ptr(other.ptr) {
    if (counter == nullptr) {
        return;
    }
    counter->increment_shared_count();
}

template <typename T>
template <typename U>
SharedPtr<T>::SharedPtr(SharedPtr<U>&& other): counter(other.counter), ptr(other.ptr){
    other.ptr = nullptr;
    other.counter = nullptr;
}

template <typename T>
SharedPtr<T>::SharedPtr(SharedPtr<T>&& other): counter(other.counter), ptr(other.ptr) {
    other.ptr = nullptr;
    other.counter = nullptr;
}

template <typename T>
template <typename U>
SharedPtr<T>& SharedPtr<T>::operator=(U&& other) {
    SharedPtr<T>(std::forward<U>(other)).swap(*this);
    return *this;
}

template <typename T>
template <typename U, typename Deleter>
SharedPtr<T>::SharedPtr(U* ptr, Deleter d): SharedPtr(ptr, d, std::allocator<U>()) {}

template <typename T>
template <typename U, typename Deleter, typename Alloc>
SharedPtr<T>::SharedPtr(U* ptr, Deleter d, Alloc alloc): counter(nullptr), ptr(ptr) {
    using ControlBlockType = ControlBlock<U, Alloc, Deleter>;
    using AllocControlBlock = typename std::allocator_traits<Alloc>::template rebind_alloc<ControlBlockType>;
    AllocControlBlock alloc_control_block = alloc;
    counter = alloc_control_block.allocate(1);
    new (counter) ControlBlockType(d, alloc);
    //std::allocator_traits<AllocControlBlock>::construct(alloc_control_block, counter, d, alloc);
}



template <typename T>
SharedPtr<T>::~SharedPtr() {
    if (counter == nullptr) {
        return;
    }
    
    if (counter->get_shared_count() == 1) {
        counter->operator()(reinterpret_cast<void*>(ptr));
        if (counter->get_weak_count() == 0) {
            counter->destroy();
        } else {
            counter->decrement_shared_count();
        }
        return;
    }
    counter->decrement_shared_count();
}

template <typename T, typename... Args>
SharedPtr<T> makeShared(Args&&... args) {
    return SharedPtr<T>(typename SharedPtr<T>::construct_tag(), std::allocator<T>(), std::forward<Args>(args)...);
}

template <typename T, typename Alloc, typename... Args>
SharedPtr<T> allocateShared(const Alloc& alloc, Args&&... args) {
    return SharedPtr<T>(typename SharedPtr<T>::construct_tag(), alloc, std::forward<Args>(args)...);
}

template <typename T>
template <typename Alloc, typename... Args>
SharedPtr<T>::SharedPtr(typename SharedPtr<T>::construct_tag, Alloc alloc, Args&&... args) {
    using ControlDataType = ControlData<T, Alloc, EmptyDeleter<T, Alloc>>;
    using ControlBlockType = ControlBlock<T, Alloc, EmptyDeleter<T, Alloc>>;
    using AllocData = typename std::allocator_traits<Alloc>::template rebind_alloc<ControlDataType>;
    AllocData alloc_data = alloc;
    
    auto* data_ptr = alloc_data.allocate(1);
    counter = &(data_ptr->control_block);
    ptr = &(data_ptr->object);
    
    std::allocator_traits<Alloc>::construct(alloc, ptr, std::forward<Args>(args)...);
    new (counter) ControlBlockType(EmptyDeleter<T, Alloc>(alloc), alloc);
}

template <typename T>
SharedPtr<T>::SharedPtr(typename SharedPtr<T>::construct_control_tag, T* ptr, BaseControlBlock* counter): counter(counter), ptr(ptr) {
    counter->increment_shared_count();
}

template <typename T>
size_t SharedPtr<T>::use_count() const {
    return (counter ? counter->get_shared_count() : 0);
}

template <typename T>
void SharedPtr<T>::reset() {
    SharedPtr<T>().swap(*this);
}

template <typename T>
template <typename U>
void SharedPtr<T>::reset(U* ptr) {
    if (ptr == nullptr) {
        reset();
        return;
    }
    SharedPtr<T>(ptr).swap(*this);
}

template <typename T>
void SharedPtr<T>::swap(SharedPtr<T>& other) {
    std::swap(counter, other.counter);
    std::swap(ptr, other.ptr);
}

template <typename T>
T* SharedPtr<T>::get() const {
    return ptr;
}

template <typename T>
T& SharedPtr<T>::operator*() const {
    return *ptr;
}

template <typename T>
T* SharedPtr<T>::operator->() const {
    return ptr;
}

template <typename T>
WeakPtr<T>::WeakPtr() : counter(nullptr), ptr(nullptr) {}

template <typename T>
template <typename U>
WeakPtr<T>::WeakPtr(const WeakPtr<U>& other): counter(other.counter), ptr(other.ptr) {
    if (counter == nullptr) {
        return;
    }
    counter->increment_weak_count();
}

template <typename T>
WeakPtr<T>::WeakPtr(const WeakPtr<T>& other): counter(other.counter), ptr(other.ptr) {
    if (counter == nullptr) {
        return;
    }
    counter->increment_weak_count();
}

template <typename T>
template <typename U>
WeakPtr<T>::WeakPtr(const SharedPtr<U>& other): counter(other.counter), ptr(other.ptr){
    counter->increment_weak_count();
}

template <typename T>
template <typename U>
WeakPtr<T>::WeakPtr(WeakPtr<U>&& other): counter(other.counter), ptr(other.ptr) {
    other.ptr = nullptr;
    other.counter = nullptr;
}

template <typename T>
WeakPtr<T>::WeakPtr(WeakPtr<T>&& other): counter(other.counter), ptr(other.ptr) {
    other.ptr = nullptr;
    other.counter = nullptr;
}

template <typename T>
WeakPtr<T>::~WeakPtr() {
    if (counter == nullptr) {
        return;
    }
    counter->decrement_weak_count();
    if (counter->get_weak_count() == 0 && counter->get_shared_count() == 0) {
        counter->destroy();
    }
}

template <typename T>
template <typename U>
WeakPtr<T>& WeakPtr<T>::operator=(const SharedPtr<U>& other) {
    WeakPtr<T>(other).swap(*this);
    return *this;
}

template <typename T>
template <typename Other>
WeakPtr<T>& WeakPtr<T>::operator=(Other &&other) {
    WeakPtr<T>(std::forward<Other>(other)).swap(*this);
    return *this;
}

template <typename T>
size_t WeakPtr<T>::use_count() const {
    return (counter ? counter->get_shared_count() : 0);
}

template <typename T>
void WeakPtr<T>::swap(WeakPtr<T>& other) {
    std::swap(counter, other.counter);
    std::swap(ptr, other.ptr);
}

template <typename T>
T* WeakPtr<T>::get() const {
    return ptr;
}

template <typename T>
bool WeakPtr<T>::expired() const {
    return counter->get_shared_count() == 0;
}

template <typename T>
SharedPtr<T> WeakPtr<T>::lock() const {
    return expired() ? SharedPtr<T>() : SharedPtr<T>(typename SharedPtr<T>::construct_control_tag(), ptr, counter);
}

template <typename T>
T& WeakPtr<T>::operator*() const {
    return *ptr;
}

template <typename T>
T* WeakPtr<T>::operator->() const {
    return ptr;
}


