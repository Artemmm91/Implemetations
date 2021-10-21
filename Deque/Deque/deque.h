#include <vector>

template <typename ElementType>
class Deque {
private:
    std::vector<ElementType> elements = {};
    size_t _begin = 0;
    size_t _end = 0;
    int delta = 0;
    
    void changeBuffer(size_t newSize);
    void increaseBuffer();
    
    template <typename IteratorType>
    class baseIterator;
    
public:
    using iterator = baseIterator<ElementType>;
    using const_iterator = baseIterator<const ElementType>;
    
    Deque();
    Deque(int newSize);
    Deque(int newSize, ElementType element);
    Deque(const Deque<ElementType>& another) = default;
    Deque<ElementType>& operator=(const Deque<ElementType>& another) = default;
    
    size_t size() const;
    
    ElementType& operator[](size_t index);
    const ElementType& operator[](size_t index) const;
    ElementType& at(size_t index);
    const ElementType& at(size_t index) const;
    
    void push_back(const ElementType& newElement);
    void push_front(const ElementType& newElement);
    void pop_back();
    void pop_front();
                
    const_iterator begin() const;
    const_iterator cbegin() const;
    const_iterator end() const;
    const_iterator cend() const;
    iterator end();
    iterator begin();
    
    void insert(iterator elementIterator, const ElementType& newElement);
    void erase(iterator elementIterator);
};

template <typename ElementType>
template <typename IteratorType>
class Deque<ElementType>::baseIterator {
protected:
    Deque<ElementType>* thisDeque;
    int index;
    
public:
    baseIterator(Deque<ElementType>* ElementTypehisDeque, int index);
    baseIterator(const Deque<ElementType>* ElementTypehisDeque, int index);
    size_t currentIndex() const;
    
    baseIterator& operator++();
    baseIterator operator++(int);
    baseIterator& operator--();
    baseIterator operator--(int);
    
    int operator-(const baseIterator& another) const;
    
    bool operator<(const baseIterator& another) const;
    bool operator>(const baseIterator& another) const;
    bool operator<=(const baseIterator& another) const;
    bool operator>=(const baseIterator& another) const;
    bool operator==(const baseIterator& another) const;
    bool operator!=(const baseIterator& another) const;
    
    IteratorType& operator*() const;
    IteratorType* operator->() const;
    
    baseIterator operator+(int number) const {
        return baseIterator(thisDeque, number + index);
    }
    friend baseIterator operator+(int number, const baseIterator& thisConstIterator) {
        return baseIterator(thisConstIterator.thisDeque, thisConstIterator.index + number);
    }
    
    baseIterator operator-(int number) const {
        return baseIterator(thisDeque, -number + index);
    }
};



//-------------------------------- Definitions -------------------------------- //



template <typename ElementType>
void Deque<ElementType>::changeBuffer(size_t newSize) {
    std::vector<ElementType> newElements(newSize * 3);
    for (size_t i = 0; i < size(); ++i) {
        newElements[newSize + i] = elements[_begin + i];
    }
    _end = newSize + size();
    _begin = newSize;
    elements = newElements;
    delta = static_cast<int>(newSize);
}
template <typename ElementType>
void Deque<ElementType>::increaseBuffer() {
    size_t newSize = elements.size();
    if (size() * 3 <= elements.size()) {
        newSize = size();
    }
    changeBuffer(newSize);
    return;
}

template <typename ElementType>
Deque<ElementType>::Deque(): elements(std::vector<ElementType>(3)),
                             _begin(1), _end(1), delta(1) {}

template <typename ElementType>
Deque<ElementType>::Deque(int newSize): elements(std::vector<ElementType>(newSize * 3)),
                                        _begin(newSize), _end(newSize * 2), delta(newSize) {}

template <typename ElementType>
Deque<ElementType>::Deque(int newSize, ElementType element):
        elements(std::vector<ElementType>(newSize * 3, element)),
        _begin(newSize), _end(newSize * 2), delta(newSize) {}

template <typename ElementType>
size_t Deque<ElementType>::size() const {
    return _end - _begin;
}

template <typename ElementType>
ElementType& Deque<ElementType>::operator[](size_t index) {
    return elements[_begin + index];
}

template <typename ElementType>
const ElementType& Deque<ElementType>::operator[](size_t index) const {
    return elements[_begin + index];
}

template <typename ElementType>
ElementType& Deque<ElementType>::at(size_t index) {
    if (index + _begin >= _end || index < 0) {
        throw std::out_of_range("Deque::at out of range");
        return elements[_begin];
    }
    return elements[_begin + index];
}

template <typename ElementType>
const ElementType& Deque<ElementType>::at(size_t index) const {
    if (index + _begin >= _end || index < 0) {
        throw std::out_of_range("Deque::at out of range");
        return elements[_begin];
    }
    return elements[_begin + index];
}

template <typename ElementType>
void Deque<ElementType>::push_back(const ElementType& newElement) {
    if (_end == elements.size()) {
        increaseBuffer();
    }
    elements[_end] = newElement;
    ++_end;
}

template <typename ElementType>
void Deque<ElementType>::push_front(const ElementType& newElement) {
    if (_begin == 0) {
        increaseBuffer();
    }
    --_begin;
    elements[_begin] = newElement;
}

template <typename ElementType>
void Deque<ElementType>::pop_back() {
    --_end;
}

template <typename ElementType>
void Deque<ElementType>::pop_front() {
    ++_begin;
}

template <typename ElementType>
template <typename IteratorType>
Deque<ElementType>::baseIterator<IteratorType>::
                    baseIterator(const Deque<ElementType>* thisDeque, int index):
                            thisDeque(const_cast<Deque<ElementType>*>(thisDeque)),
                            index(index) {}

template <typename ElementType>
template <typename IteratorType>
Deque<ElementType>::baseIterator<IteratorType>::
                    baseIterator(Deque<ElementType>* thisDeque, int index):
                            thisDeque(thisDeque), index(index) {}


template <typename ElementType>
template <typename IteratorType>
size_t Deque<ElementType>::baseIterator<IteratorType>::currentIndex() const {
    return index + thisDeque->delta;
}

template <typename ElementType>
template <typename IteratorType>
typename Deque<ElementType>::template baseIterator<IteratorType>&
         Deque<ElementType>::baseIterator<IteratorType>::operator++() {
    ++index;
    return *this;
}

template <typename ElementType>
template <typename IteratorType>
typename Deque<ElementType>::template baseIterator<IteratorType>
         Deque<ElementType>::baseIterator<IteratorType>::operator++(int) {
    auto another = *this;
    ++(*this);
    return another;
}

template <typename ElementType>
template <typename IteratorType>
typename Deque<ElementType>::template baseIterator<IteratorType>&
         Deque<ElementType>::baseIterator<IteratorType>::operator--() {
    --index;
    return *this;
}
template <typename ElementType>
template <typename IteratorType>
typename Deque<ElementType>::template baseIterator<IteratorType>
         Deque<ElementType>::baseIterator<IteratorType>::operator--(int) {
    auto another = *this;
    --(*this);
    return another;
}

template <typename ElementType>
template <typename IteratorType>
IteratorType& Deque<ElementType>::baseIterator<IteratorType>::operator*() const {
    return thisDeque->elements[currentIndex()];
}

template <typename ElementType>
template <typename IteratorType>
IteratorType* Deque<ElementType>::baseIterator<IteratorType>::operator->() const {
    return &thisDeque->elements[currentIndex()];
}


template <typename ElementType>
template <typename IteratorType>
bool Deque<ElementType>::baseIterator<IteratorType>::
         operator<(const baseIterator& another) const {
    return index < another.index;
}

template <typename ElementType>
template <typename IteratorType>
bool Deque<ElementType>::baseIterator<IteratorType>::
        operator>(const baseIterator& another) const {
    return another < (*this);
}

template <typename ElementType>
template <typename IteratorType>
bool Deque<ElementType>::baseIterator<IteratorType>::
        operator<=(const baseIterator& another) const {
    return !(another < (*this));
}

template <typename ElementType>
template <typename IteratorType>
bool Deque<ElementType>::baseIterator<IteratorType>::
        operator>=(const baseIterator& another) const {
    return !((*this) < another);
}

template <typename ElementType>
template <typename IteratorType>
bool Deque<ElementType>::baseIterator<IteratorType>::
        operator==(const baseIterator& another) const {
    return !((*this) < another) && !(another < (*this));
}

template <typename ElementType>
template <typename IteratorType>
bool Deque<ElementType>::baseIterator<IteratorType>::
        operator!=(const baseIterator& another) const {
    return ((*this) < another) || (another < (*this));
}

template <typename ElementType>
template <typename IteratorType>
int Deque<ElementType>::baseIterator<IteratorType>::
        operator-(const baseIterator& another) const {
    return index - another.index;
}

template <typename ElementType>
typename Deque<ElementType>::iterator Deque<ElementType>::begin() {
    return iterator(this, static_cast<int>(_begin) - delta);
}

template <typename ElementType>
typename Deque<ElementType>::const_iterator Deque<ElementType>::begin() const {
    return cbegin();
}

template <typename ElementType>
typename Deque<ElementType>::const_iterator Deque<ElementType>::cbegin() const {
    return const_iterator(this, static_cast<int>(_begin) - delta);
}

template <typename ElementType>
typename Deque<ElementType>::iterator Deque<ElementType>::end() {
    return iterator(this, static_cast<int>(_end) - delta);
}

template <typename ElementType>
typename Deque<ElementType>::const_iterator Deque<ElementType>::end() const {
    return cend();
}

template <typename ElementType>
typename Deque<ElementType>::const_iterator Deque<ElementType>::cend() const {
    return const_iterator(this, static_cast<int>(_end) - delta);
}

template<typename ElementType>
void Deque<ElementType>::insert(iterator elementIterator, const ElementType& newElement) {
    push_back(newElement);
    for (auto i = end() - 1; i > elementIterator; --i) {
        std::swap(*i, *(i - 1));
    }
}
 
template<typename ElementType>
void Deque<ElementType>::erase(iterator elementIterator) {
    for (auto i = elementIterator; i < end(); ++i) {
        std::swap(*i, *(i + 1));
    }
    pop_back();
}
