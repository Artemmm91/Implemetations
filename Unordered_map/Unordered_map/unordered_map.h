#include <cmath>
#include <functional>
#include <cstddef>
#include <vector>
#include <iterator>
#include <optional>
#include <type_traits>


template <typename T, typename Allocator = std::allocator<T>>
class List {
private:
    struct Node {
        Node* next = nullptr;
        Node* prev = nullptr;
        T data;
        
        Node() = default;
        Node(const T& value): data(value) {}
        void unlink() {
            Node* nextNode = next;
            Node* prevNode = prev;
            if (next)
                next->prev = prevNode;
            if (prev)
                prev->next = nextNode;
            next = nullptr;
            prev = nullptr;
        }
    };
    
    template <typename IteratorType>
    class baseIterator;
    
    Node* beginNode = nullptr;
    Node* endNode = nullptr;
    size_t listSize = 0;
    Allocator alloc;
    
    void clear();
    void copy(const List& another);
    
    using nodeAllocatorType = typename std::allocator_traits<Allocator>::template rebind_alloc<Node>;
    nodeAllocatorType nodeAlloc;
    void push_back();
    
    template <typename U>
    Node* createElement(U&& value);
    
    void clearElement(Node* element);

public:
    using iterator = baseIterator<T>;
    using const_iterator = baseIterator<const T>;
    using reverse_iterator = std::reverse_iterator<baseIterator<T>>;
    using const_reverse_iterator = std::reverse_iterator<baseIterator<const T>>;
    
    explicit List(const Allocator& alloc = Allocator());
    List(size_t count, const T& value, const Allocator& alloc = Allocator());
    List(size_t count, const Allocator& alloc = Allocator());
    
    List(const List& another);
    List(List&& another);
    
    List& operator=(const List& right);
    List& operator=(List&& right);
    ~List();

    size_t size() const;
    void pop_front();
    void pop_back();
    
    
    template <typename U>
    void push_front(U&& value);
    
    template <typename U>
    void push_back(U&& value);
    
    Allocator get_allocator() const;
    
    const_iterator begin() const {return const_iterator(beginNode);}
    const_iterator cbegin() const {return const_iterator(beginNode);}
    const_iterator end() const {return const_iterator(endNode);}
    const_iterator cend() const {return const_iterator(endNode);}
    iterator end() {return iterator(endNode);}
    iterator begin() {return iterator(beginNode);}
    
    const_reverse_iterator rbegin() const {return const_reverse_iterator(endNode);}
    const_reverse_iterator crbegin() const {return const_reverse_iterator(endNode);}
    const_reverse_iterator rend() const {return const_reverse_iterator(beginNode);}
    const_reverse_iterator crend() const {return const_reverse_iterator(beginNode);}
    reverse_iterator rend() {return reverse_iterator(beginNode);}
    reverse_iterator rbegin() {return reverse_iterator(endNode);}
    
    template <typename U>
    void insert(const_iterator thisIterator, U&& value);
    void erase(const_iterator thisIterator);
};


template <typename T, typename Allocator>
void List<T, Allocator>::clear() {
    while (listSize > 0) {
        pop_back();
    }
}

template <typename T, typename Allocator>
void List<T, Allocator>::copy(const List& another) {
    Node* it = another.beginNode;
    while (it != another.endNode) {
        push_back(it->data);
        it = it->next;
    }
}

template <typename T, typename Allocator>
List<T, Allocator>::List(const Allocator& newAlloc) {
    endNode = std::allocator_traits<nodeAllocatorType>::allocate(nodeAlloc, 1);
    endNode->prev = nullptr;
    endNode->next = nullptr;
    beginNode = endNode;
    alloc = newAlloc;
}

template <typename T, typename Allocator>
List<T, Allocator>::List(size_t count, const T& value, const Allocator& alloc): List(alloc) {
    endNode = std::allocator_traits<nodeAllocatorType>::allocate(nodeAlloc, 1);
    endNode->prev = nullptr;
    endNode->next = nullptr;
    beginNode = endNode;
    for (size_t i = 0; i < count; ++i) {
        push_back(value);
    }
}

template <typename T, typename Allocator>
List<T, Allocator>::List(size_t count, const Allocator& alloc): List(alloc) {
    endNode = std::allocator_traits<nodeAllocatorType>::allocate(nodeAlloc, 1);
    endNode->prev = nullptr;
    endNode->next = nullptr;
    beginNode = endNode;
    for (size_t i = 0; i < count; ++i) {
        push_back();
    }
}


template <typename T, typename Allocator>
List<T, Allocator>::List(const List& another):
        List(std::allocator_traits<Allocator>::select_on_container_copy_construction(another.get_allocator())) {
    endNode = std::allocator_traits<nodeAllocatorType>::allocate(nodeAlloc, 1);
    endNode->prev = nullptr;
    endNode->next = nullptr;
    beginNode = endNode;
    copy(another);
}

template <typename T, typename Allocator>
List<T, Allocator>::List(List&& another):
        beginNode(another.beginNode), endNode(another.endNode),listSize(another.listSize),
        alloc(std::move(another.alloc)), nodeAlloc(std::move(another.nodeAlloc)) {
    another.beginNode = nullptr;
    another.endNode = nullptr;
    another.listSize = 0;
}

template <typename T, typename Allocator>
List<T, Allocator>& List<T, Allocator>::operator=(const List<T, Allocator>& another) {
    if (this != &another) {
        clear();
        if(std::allocator_traits<Allocator>::propagate_on_container_copy_assignment::value){
            alloc = another.get_allocator();
        }
        copy(another);
    }
    return *this;
}

template <typename T, typename Allocator>
List<T, Allocator>& List<T, Allocator>::operator=(List<T, Allocator>&& right) {
    if (this != &right) {
        beginNode = right.beginNode;
        endNode = right.endNode;
        listSize = right.listSize;
        right.beginNode = nullptr;
        right.endNode = nullptr;
        right.listSize = 0;
    }
    return *this;
}

template <typename T, typename Allocator>
List<T, Allocator>::~List() {
    clear();
    std::allocator_traits<nodeAllocatorType>::deallocate(nodeAlloc, endNode, 1);
}


template <typename T, typename Allocator>
size_t List<T, Allocator>::size() const {
    return listSize;
}

template <typename T, typename Allocator>
template <typename U>
typename List<T, Allocator>::Node* List<T, Allocator>::createElement(U&& value){
    Node* newNode = std::allocator_traits<nodeAllocatorType>::allocate(nodeAlloc, 1);
    std::allocator_traits<nodeAllocatorType>::construct(nodeAlloc, newNode, std::forward<U>(value));
    return newNode;
}

template <typename T, typename Allocator>
void List<T, Allocator>::clearElement(Node* element){
    std::allocator_traits<nodeAllocatorType>::destroy(nodeAlloc, element);
    std::allocator_traits<nodeAllocatorType>::deallocate(nodeAlloc, element, 1);
}

template <typename T, typename Allocator>
void List<T, Allocator>::pop_front() {
    if (listSize > 0) {
        Node* newBeginNode = beginNode->next;
        Node* oldBeginNode = beginNode;
        oldBeginNode->unlink();
        beginNode = newBeginNode;
        listSize--;
        
        clearElement(oldBeginNode);
    }
}

template <typename T, typename Allocator>
void List<T, Allocator>::pop_back() {
    if (listSize > 0) {
        Node* oldLastNode = endNode->prev;
        oldLastNode->unlink();
        listSize--;
        
        clearElement(oldLastNode);
    }
    if (listSize == 0) {
        beginNode = endNode;
    }
}


template <typename T, typename Allocator>
template <typename U>
void List<T, Allocator>::push_front(U&& value) {
    Node* newNode = createElement(std::forward<U>(value));

    beginNode->prev = newNode;
    newNode->next = beginNode;
    beginNode = newNode;
    ++listSize;
}

template <typename T, typename Allocator>
template <typename U>
void List<T, Allocator>::push_back(U&& value) {
    Node* newNode = createElement(std::forward<U>(value));

    newNode->next = endNode;
    newNode->prev = endNode->prev;
    endNode->prev = newNode;
    if(newNode->prev != nullptr) {
        newNode->prev->next = newNode;
    } else {
        beginNode = newNode;
    }
    
    ++listSize;
}

template <typename T, typename Allocator>
void List<T, Allocator>::push_back() {
    Node* newNode = std::allocator_traits<nodeAllocatorType>::allocate(nodeAlloc, 1);
    std::allocator_traits<nodeAllocatorType>::construct(nodeAlloc, newNode);

    newNode->next = endNode;
    newNode->prev = endNode->prev;
    endNode->prev = newNode;
    if(newNode->prev != nullptr) {
        newNode->prev->next = newNode;
    } else {
        beginNode = newNode;
    }
    
    ++listSize;
}

template <typename T, typename Allocator>
Allocator List<T, Allocator>::get_allocator() const {
    return alloc;
}


template <typename ElementType, typename Allocator>
template <typename IteratorType>
class List<ElementType, Allocator>::baseIterator {
public:
    Node* thisNode = nullptr;
    
public:
    using iterator_category = std::bidirectional_iterator_tag;
    using value_type = IteratorType;
    using reference = IteratorType&;
    using pointer = IteratorType*;
    using difference_type = std::ptrdiff_t;
    
    
    baseIterator(Node* node): thisNode(node) {}
    operator baseIterator<const IteratorType>() {
           return baseIterator<const IteratorType>(thisNode);
    }
    
    baseIterator& operator++();
    baseIterator operator++(int);
    baseIterator& operator--();
    baseIterator operator--(int);
    
    int operator-(const baseIterator& another) const;
    
    bool operator==(const baseIterator& another) const;
    bool operator!=(const baseIterator& another) const;
    
    IteratorType& operator*() const;
    IteratorType* operator->() const;
    
    baseIterator operator+(int number) const;
    baseIterator operator-(int number) const;
    
    friend baseIterator operator+(int number, const baseIterator& thisIterator) {
        return thisIterator + number;
    }
    
    bool is_end() {
        return thisNode == nullptr || thisNode->next == nullptr;
    }
};


template <typename ElementType, typename Allocator>
template <typename IteratorType>
typename List<ElementType, Allocator>::template baseIterator<IteratorType>&
        List<ElementType, Allocator>::baseIterator<IteratorType>::operator++() {
    thisNode = thisNode->next;
    return *this;
}

template <typename ElementType, typename Allocator>
template <typename IteratorType>
typename List<ElementType, Allocator>::template baseIterator<IteratorType>&
        List<ElementType, Allocator>::baseIterator<IteratorType>::operator--() {
    thisNode = thisNode->prev;
    return *this;
}

template <typename ElementType, typename Allocator>
template <typename IteratorType>
typename List<ElementType, Allocator>::template baseIterator<IteratorType>
        List<ElementType, Allocator>::baseIterator<IteratorType>::operator++(int) {
    auto another = *this;
    ++(*this);
    return another;
}

template <typename ElementType, typename Allocator>
template <typename IteratorType>
typename List<ElementType, Allocator>::template baseIterator<IteratorType>
        List<ElementType, Allocator>::baseIterator<IteratorType>::operator--(int) {
    auto another = *this;
    --(*this);
    return another;
}

template <typename ElementType, typename Allocator>
template <typename IteratorType>
IteratorType& List<ElementType, Allocator>::baseIterator<IteratorType>::operator*() const {
    return thisNode->data;
}

template <typename ElementType, typename Allocator>
template <typename IteratorType>
IteratorType* List<ElementType, Allocator>::baseIterator<IteratorType>::operator->() const {
    return &thisNode->data;
}

template <typename ElementType, typename Allocator>
template <typename IteratorType>
bool List<ElementType, Allocator>::baseIterator<IteratorType>::operator==(const baseIterator& another) const {
    return thisNode == another.thisNode;
}

template <typename ElementType, typename Allocator>
template <typename IteratorType>
bool List<ElementType, Allocator>::baseIterator<IteratorType>::operator!=(const baseIterator& another) const {
    return thisNode != another.thisNode;
}

template <typename ElementType, typename Allocator>
template <typename IteratorType>
typename List<ElementType, Allocator>::template baseIterator<IteratorType>
        List<ElementType, Allocator>::baseIterator<IteratorType>::operator+(int number) const {
    if (number < 0) {
        return *this - (-number);
    }
    Node* iterateNode = thisNode;
    while (number > 0) {
        iterateNode = thisNode->next;
        --number;
    }
    return iterateNode;
}

template <typename ElementType, typename Allocator>
template <typename IteratorType>
typename List<ElementType, Allocator>::template baseIterator<IteratorType>
        List<ElementType, Allocator>::baseIterator<IteratorType>::operator-(int number) const {
    if (number < 0) {
        return *this + (-number);
    }
    Node* iterateNode = thisNode;
    while (number > 0) {
        iterateNode = thisNode->prev;
        --number;
    }
    return iterateNode;
}


template <typename ElementType, typename Allocator>
template <typename U>
void List<ElementType, Allocator>::insert(const_iterator thisIterator, U&& value) {
    Node* newNode = createElement(std::forward<U>(value));
    
    Node* thisNode = thisIterator.thisNode;
    newNode->next = thisNode;
    newNode->prev = thisNode->prev;
    if (thisNode->prev == nullptr) {
        beginNode = newNode;
    } else {
        thisNode->prev->next = newNode;
    }
    thisNode->prev = newNode;
    ++listSize;
}

template <typename ElementType, typename Allocator>
void List<ElementType, Allocator>::erase(const_iterator thisIterator) {
    Node* thisNode = thisIterator.thisNode;
    thisNode->next->prev = thisNode->prev;
    if (thisNode->prev == nullptr) {
        beginNode = thisNode->next;
    } else {
        thisNode->prev->next = thisNode->next;
    }
    --listSize;
    clearElement(thisNode);
}


template <typename Key, typename Value,
    typename Hash = std::hash<Key>, typename Equal = std::equal_to<Key>,
    typename Alloc = std::allocator<std::pair<const Key, Value>>>
class UnorderedMap {
public:
    using NodeType = std::pair<const Key, Value>;
private:
    using list_alloc = typename std::allocator_traits<Alloc>::
            template rebind_alloc<NodeType*>;
    using list_iterator = typename List<NodeType*, list_alloc>::iterator;
    using const_list_iterator = typename List<NodeType*, list_alloc>::const_iterator;
    using iterator_alloc = typename std::allocator_traits<Alloc>::
            template rebind_alloc<typename List<NodeType*, list_alloc>::iterator>;
    Hash hash_function;
    Equal equal_function;
    Alloc node_alloc;
    
    List<NodeType*, list_alloc> nodes;
    std::vector<list_iterator, iterator_alloc> pointers_to_beginning;
    std::vector<list_iterator, iterator_alloc> pointers_to_ends;
    size_t chunks = 0;
    float _max_load_factor = 0.8;
    size_t last_hash_node = 0;
    
    void restore_map_to_default();
    size_t current_hash(const Key& key) const;
    
    template <typename IteratorType>
    class base_iterator;
    
public:
    using const_iterator = base_iterator<const NodeType>;
    using iterator = base_iterator<NodeType>;

    UnorderedMap();
    UnorderedMap(const UnorderedMap& other);
    UnorderedMap(UnorderedMap&& other);

    UnorderedMap& operator=(const UnorderedMap & other);
    UnorderedMap& operator=(UnorderedMap&& other);

    ~UnorderedMap() {
        for (auto it = nodes.begin(); it != nodes.end(); ++it) {
              std::allocator_traits<Alloc>::destroy(node_alloc, (*it));
              std::allocator_traits<Alloc>::deallocate(node_alloc, (*it), 1);
        }
    }
    
    Value& operator[](const Key& key);
    Value& operator[](Key&& key);
    Value& at(const Key& key);
    const Value& at(const Key& key) const;
    
    size_t size() const;
    
    iterator begin();
    const_iterator begin() const;
    const_iterator cbegin() const ;
    iterator end();
    const_iterator end() const;
    const_iterator cend() const;
    
    std::pair<iterator, bool> insert(const NodeType& node);
    template <typename U>
    std::pair<iterator, bool> insert(U&& value);
    
    template <typename InputIterator>
    void insert(InputIterator first, InputIterator last);
    
    template <typename... Args >
    std::pair<iterator, bool> emplace(Args&&... args);
    
    iterator erase(const_iterator position);
    iterator erase(const_iterator first, const_iterator last);
    
    template <typename U>
    iterator find(U&& key);
    
    template <typename U>
    const_iterator find(U&& key) const;
    
    void reserve(size_t count);
    size_t max_size() const;
    float load_factor() const;
    float max_load_factor() const;
    void max_load_factor(float new_max_load_factor);
    
private:
    iterator deconst_convert(const_iterator it);
    void rehash(size_t count);
    void check_rehash();
    std::pair<iterator, bool> insert_node(NodeType* new_node);
    bool is_iterators_equal(const_iterator map_iter, list_iterator list_iter);
    
    template <typename Iterator, typename U>
    std::optional<Iterator> base_find(U&& key, size_t hash_node) const;
    
    template <typename U>
    iterator try_find(U&& key, size_t& hash_node);
};


template <typename Key, typename Value, typename Hash, typename Equal, typename Alloc>
template <typename IteratorType>
class UnorderedMap<Key, Value, Hash, Equal, Alloc>::base_iterator {
private:
    friend class UnorderedMap;
    
    using list_iterator_type = typename std::conditional_t<std::is_const<IteratorType>::value,
                                            const_list_iterator, list_iterator>;
    list_iterator_type this_iterator;
    
    base_iterator(list_iterator it): this_iterator(it) {}
    base_iterator(const_list_iterator it): this_iterator(it) {}
    
    list_iterator_type get_list_iterator() {
        return this_iterator;
    }
    
    bool is_end() {
        return this_iterator.is_end();
    }
public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = IteratorType;
    using reference = IteratorType&;
    using pointer = IteratorType*;
    using difference_type = std::ptrdiff_t;
    
    operator base_iterator<const IteratorType>() {
           return base_iterator<const IteratorType>(this_iterator);
    }
    
    bool operator==(const base_iterator& another) const {
        return another.this_iterator == this_iterator;
    }
    bool operator!=(const base_iterator& another) const {
        return !(*this == another);
    }
    
    base_iterator& operator++() {
        ++this_iterator;
        return *this;
    }
    base_iterator operator++(int) {
        auto current_it = *this;
        ++this_iterator;
        return current_it;
    }
    
    base_iterator& operator--() {
        --this_iterator;
        return *this;
    }
    base_iterator operator--(int) {
        auto current_it = *this;
        --this_iterator;
        return current_it;
    }
    
    IteratorType& operator*() const {
        return *(*this_iterator);
    }
    IteratorType* operator->() const {
        return (*this_iterator);
    }
};

template <typename Key, typename Value, typename Hash, typename Equal, typename Alloc>
UnorderedMap<Key, Value, Hash, Equal, Alloc>::UnorderedMap():
        nodes(List<NodeType*, list_alloc>()),
        pointers_to_beginning(std::vector<list_iterator, iterator_alloc>(1, nodes.end())),
        pointers_to_ends(std::vector<list_iterator, iterator_alloc>(1, nodes.end())),
        chunks(1) {
}

template <typename Key, typename Value, typename Hash, typename Equal, typename Alloc>
UnorderedMap<Key, Value, Hash, Equal, Alloc>::UnorderedMap(const UnorderedMap& other):
        nodes(List<NodeType*, list_alloc>()),
        pointers_to_beginning(std::vector<list_iterator, iterator_alloc>(other.chunks, nodes.end())),
        pointers_to_ends(std::vector<list_iterator, iterator_alloc>(other.chunks, nodes.end())),
        chunks(other.chunks),
        _max_load_factor(other._max_load_factor),
        last_hash_node(other.last_hash_node) {
    if(std::allocator_traits<Alloc>::propagate_on_container_copy_assignment::value){
        node_alloc = other.node_alloc;
    }
    for (auto it = other.cbegin(); it != other.cend(); ++it) {
        insert(*it);
    }
}

template <typename Key, typename Value, typename Hash, typename Equal, typename Alloc>
UnorderedMap<Key, Value, Hash, Equal, Alloc>::UnorderedMap(UnorderedMap&& other):
        node_alloc(std::move(other.node_alloc)),
        nodes(std::move(other.nodes)),
        pointers_to_beginning(std::move(other.pointers_to_beginning)),
        pointers_to_ends(std::move(other.pointers_to_ends)),
        chunks(other.chunks),
        _max_load_factor(other._max_load_factor),
        last_hash_node(other.last_hash_node) {
    other.restore_map_to_default();
}

template <typename Key, typename Value, typename Hash, typename Equal, typename Alloc>
UnorderedMap<Key, Value, Hash, Equal, Alloc>&
        UnorderedMap<Key, Value, Hash, Equal, Alloc>::operator=(const UnorderedMap& other) {
    if(std::allocator_traits<Alloc>::propagate_on_container_copy_assignment::value){
        node_alloc = other.node_alloc;
    }
    nodes = other.nodes;
    chunks = other.chunks;
    pointers_to_beginning = std::vector<list_iterator, iterator_alloc>(chunks, nodes.end());
    pointers_to_ends = std::vector<list_iterator, iterator_alloc>(chunks, nodes.end());
    for (auto element: nodes) {
        insert(element);
    }
    _max_load_factor = other._max_load_factor;
    last_hash_node = other.last_hash_node;
    return *this;
}

template <typename Key, typename Value, typename Hash, typename Equal, typename Alloc>
UnorderedMap<Key, Value, Hash, Equal, Alloc>&
        UnorderedMap<Key, Value, Hash, Equal, Alloc>::operator=(UnorderedMap&& other) {
    node_alloc = std::move(node_alloc);
    nodes = std::move(other.nodes);
    pointers_to_beginning = std::move(other.pointers_to_beginning);
    pointers_to_ends = std::move(other.pointers_to_ends);
    chunks = other.chunks;
    _max_load_factor = other._max_load_factor;
    last_hash_node = other.last_hash_node;
    other.restore_map_to_default();
    return *this;
}


template <typename Key, typename Value, typename Hash, typename Equal, typename Alloc>
void UnorderedMap<Key, Value, Hash, Equal, Alloc>::restore_map_to_default() {
    nodes = List<NodeType*, list_alloc>();
    pointers_to_beginning = std::vector<list_iterator, iterator_alloc>(1, nodes.end());
    pointers_to_ends = std::vector<list_iterator, iterator_alloc>(1, nodes.end());
    chunks = 1;
    _max_load_factor = 0.8;
    last_hash_node = 0;
}

template <typename Key, typename Value, typename Hash, typename Equal, typename Alloc>
size_t UnorderedMap<Key, Value, Hash, Equal, Alloc>::current_hash(const Key& key) const {
    return hash_function(key) % chunks;
}

template <typename Key, typename Value, typename Hash, typename Equal, typename Alloc>
Value& UnorderedMap<Key, Value, Hash, Equal, Alloc>::at(const Key& key) {
    auto it = find(key);
    if (it != end()) {
        return it->second;
    }
    throw std::out_of_range("No such key");
}

template <typename Key, typename Value, typename Hash, typename Equal, typename Alloc>
const Value& UnorderedMap<Key, Value, Hash, Equal, Alloc>::at(const Key& key) const {
    auto it = find(key);
    if (it != cend()) {
        return it->second;
    }
    throw std::out_of_range("No such key");
}

template <typename Key, typename Value, typename Hash, typename Equal, typename Alloc>
Value& UnorderedMap<Key, Value, Hash, Equal, Alloc>::operator[](const Key& key) {
    return (insert(std::make_pair(key, Value())).first)->second;
}

template <typename Key, typename Value, typename Hash, typename Equal, typename Alloc>
Value& UnorderedMap<Key, Value, Hash, Equal, Alloc>::operator[](Key&& key) {
    std::pair<iterator, bool> result = insert(std::make_pair(std::move(key), Value()));
    return (result.first)->second;
}

template <typename Key, typename Value, typename Hash, typename Equal, typename Alloc>
size_t UnorderedMap<Key, Value, Hash, Equal, Alloc>::size() const {
    return nodes.size();
}

template <typename Key, typename Value, typename Hash, typename Equal, typename Alloc>
typename UnorderedMap<Key, Value, Hash, Equal, Alloc>::iterator
UnorderedMap<Key, Value, Hash, Equal, Alloc>::begin() {
    return nodes.begin();
}

template <typename Key, typename Value, typename Hash, typename Equal, typename Alloc>
typename UnorderedMap<Key, Value, Hash, Equal, Alloc>::const_iterator
UnorderedMap<Key, Value, Hash, Equal, Alloc>::begin() const {
    return nodes.cbegin();
}

template <typename Key, typename Value, typename Hash, typename Equal, typename Alloc>
typename UnorderedMap<Key, Value, Hash, Equal, Alloc>::const_iterator
UnorderedMap<Key, Value, Hash, Equal, Alloc>::cbegin() const {
    return const_iterator(nodes.cbegin());
}


template <typename Key, typename Value, typename Hash, typename Equal, typename Alloc>
typename UnorderedMap<Key, Value, Hash, Equal, Alloc>::iterator
UnorderedMap<Key, Value, Hash, Equal, Alloc>::end() {
    return nodes.end();
}

template <typename Key, typename Value, typename Hash, typename Equal, typename Alloc>
typename UnorderedMap<Key, Value, Hash, Equal, Alloc>::const_iterator
UnorderedMap<Key, Value, Hash, Equal, Alloc>::end() const {
    return nodes.cend();
}

template <typename Key, typename Value, typename Hash, typename Equal, typename Alloc>
typename UnorderedMap<Key, Value, Hash, Equal, Alloc>::const_iterator
UnorderedMap<Key, Value, Hash, Equal, Alloc>::cend() const {
    return const_iterator(nodes.cend());
}

template <typename Key, typename Value, typename Hash, typename Equal, typename Alloc>
template <typename U>
std::pair<typename UnorderedMap<Key, Value, Hash, Equal, Alloc>::iterator, bool>
UnorderedMap<Key, Value, Hash, Equal, Alloc>::insert(U&& node) {
    NodeType* new_node = std::allocator_traits<Alloc>::allocate(node_alloc, 1);
    std::allocator_traits<Alloc>::construct(node_alloc, new_node, std::forward<U>(node));
    
    return insert_node(new_node);
}

template <typename Key, typename Value, typename Hash, typename Equal, typename Alloc>
std::pair<typename UnorderedMap<Key, Value, Hash, Equal, Alloc>::iterator, bool>
UnorderedMap<Key, Value, Hash, Equal, Alloc>::insert(const NodeType& node) {
    NodeType* new_node = std::allocator_traits<Alloc>::allocate(node_alloc, 1);
    std::allocator_traits<Alloc>::construct(node_alloc, new_node, node);
    
    return insert_node(new_node);
}

template <typename Key, typename Value, typename Hash, typename Equal, typename Alloc>
std::pair<typename UnorderedMap<Key, Value, Hash, Equal, Alloc>::iterator, bool>
UnorderedMap<Key, Value, Hash, Equal, Alloc>::insert_node(NodeType* new_node) {
    check_rehash();
    size_t hash_node = 0;
    iterator node_iter = try_find(new_node->first, hash_node);
    if (!node_iter.is_end()) {
        return std::make_pair(node_iter, false);
    }
    auto iter = iterator(pointers_to_beginning[hash_node]);
    
    if (iter.is_end()) {
        nodes.push_back(new_node);
        pointers_to_beginning[hash_node] = (--nodes.end());
        if (size() > 1)
            pointers_to_ends[last_hash_node] = pointers_to_beginning[hash_node];
        last_hash_node = hash_node;
        return std::make_pair(iterator(pointers_to_beginning[hash_node]), true);
    }
    
    iter = iterator(pointers_to_ends[hash_node]);
    auto prev = iterator(pointers_to_ends[hash_node]); --prev;
    nodes.template insert(iter.get_list_iterator(), new_node); ++iter;
    return std::make_pair(prev, true);
}

template <typename Key, typename Value, typename Hash, typename Equal, typename Alloc>
template<class... Args>
std::pair<typename UnorderedMap<Key, Value, Hash, Equal, Alloc>::iterator, bool>
UnorderedMap<Key, Value, Hash, Equal, Alloc>::emplace(Args&&... args) {
    NodeType* new_node = std::allocator_traits<Alloc>::allocate(node_alloc, 1);
    std::allocator_traits<Alloc>::construct(node_alloc, new_node, std::forward<Args>(args)...);
    return insert_node(new_node);
}


template <typename Key, typename Value, typename Hash, typename Equal, typename Alloc>
template <typename InputIterator>
void UnorderedMap<Key, Value, Hash, Equal, Alloc>::insert(InputIterator first, InputIterator last) {
    for (auto it = first; it != last; ++it) {
        insert(*it);
    }
}

template <typename Key, typename Value, typename Hash, typename Equal, typename Alloc>
bool UnorderedMap<Key, Value, Hash, Equal, Alloc>::is_iterators_equal(const_iterator map_iter, list_iterator list_iter) {
    if (map_iter == end() && list_iter == nodes.end()) return true;
    if (map_iter == end() || list_iter == nodes.end()) return false;
    return equal_function((*list_iter)->first, map_iter->first);
}

template <typename Key, typename Value, typename Hash, typename Equal, typename Alloc>
typename UnorderedMap<Key, Value, Hash, Equal, Alloc>::iterator
UnorderedMap<Key, Value, Hash, Equal, Alloc>::erase(const_iterator position) {
    if (position == cend()) {
        return end();
    }
    size_t hash_node = current_hash(position->first);
    iterator next = deconst_convert(position); ++next;
    if (is_iterators_equal(position, pointers_to_beginning[hash_node])) {
        if (is_iterators_equal(next, pointers_to_ends[hash_node])) {
            pointers_to_beginning[hash_node] = nodes.end();
            pointers_to_ends[hash_node] = nodes.end();
        } else {
            pointers_to_beginning[hash_node] = next.get_list_iterator();
        }
    } else if (next != end() && is_iterators_equal(next, pointers_to_ends[hash_node])) {
        iterator prev = deconst_convert(position);
        pointers_to_ends[hash_node] = prev.get_list_iterator();
    }
    nodes.erase(position.get_list_iterator());
    return next;
}

template <typename Key, typename Value, typename Hash, typename Equal, typename Alloc>
typename UnorderedMap<Key, Value, Hash, Equal, Alloc>::iterator
UnorderedMap<Key, Value, Hash, Equal, Alloc>::erase(const_iterator first, const_iterator last) {
    for (auto it = first; it != last; it = erase(it));
    return deconst_convert(last);
}

template <typename Key, typename Value, typename Hash, typename Equal, typename Alloc>
typename UnorderedMap<Key, Value, Hash, Equal, Alloc>::iterator
UnorderedMap<Key, Value, Hash, Equal, Alloc>::deconst_convert(const_iterator it) {
    if (it == cend()) {
        return end();
    }
    return find(it->first);
}


template <typename Key, typename Value, typename Hash, typename Equal, typename Alloc>
template <typename Iterator, typename U>
std::optional<Iterator> UnorderedMap<Key, Value, Hash, Equal, Alloc>::base_find(U&& key, size_t hash_node) const {
    auto it = Iterator(pointers_to_beginning[hash_node]);
    auto it_end = Iterator(pointers_to_ends[hash_node]);
    while (it != it_end) {
        if (equal_function(it->first, key)) {
            return it;
        }
        ++it;
    }
    return std::nullopt;
}


template <typename Key, typename Value, typename Hash, typename Equal, typename Alloc>
template <typename U>
typename UnorderedMap<Key, Value, Hash, Equal, Alloc>::iterator
UnorderedMap<Key, Value, Hash, Equal, Alloc>::find(U&& key) {
    size_t hash_node = current_hash(key);
    auto opt_iter = base_find<iterator>(key, hash_node);
    if (opt_iter) {
        return opt_iter.value();
    }
    return end();
}


template <typename Key, typename Value, typename Hash, typename Equal, typename Alloc>
template <typename U>
typename UnorderedMap<Key, Value, Hash, Equal, Alloc>::const_iterator
UnorderedMap<Key, Value, Hash, Equal, Alloc>::find(U&& key) const {
    size_t hash_node = current_hash(key);
    auto opt_iter = base_find<const_iterator, U>(key, hash_node);
    if (opt_iter) {
        return opt_iter.value();
    }
    return cend();
}


template <typename Key, typename Value, typename Hash, typename Equal, typename Alloc>
template <typename U>
typename UnorderedMap<Key, Value, Hash, Equal, Alloc>::iterator
UnorderedMap<Key, Value, Hash, Equal, Alloc>::try_find(U&& key, size_t& hash_node) {
    hash_node = current_hash(key);
    auto opt_iter = base_find<iterator, U>(key, hash_node);
    if (opt_iter) {
        return opt_iter.value();
    }
    return end();
}


template <typename Key, typename Value, typename Hash, typename Equal, typename Alloc>
void UnorderedMap<Key, Value, Hash, Equal, Alloc>::reserve(size_t count) {
    if (count > size()) {
        rehash(std::ceil(count / _max_load_factor));
    }
}


template <typename Key, typename Value, typename Hash, typename Equal, typename Alloc>
size_t UnorderedMap<Key, Value, Hash, Equal, Alloc>::max_size() const {
    return pointers_to_beginning.max_size();
}

template <typename Key, typename Value, typename Hash, typename Equal, typename Alloc>
float UnorderedMap<Key, Value, Hash, Equal, Alloc>::load_factor() const {
    return static_cast<float>(size()) / chunks;
}

template <typename Key, typename Value, typename Hash, typename Equal, typename Alloc>
float UnorderedMap<Key, Value, Hash, Equal, Alloc>::max_load_factor() const {
    return _max_load_factor;
}

template <typename Key, typename Value, typename Hash, typename Equal, typename Alloc>
void UnorderedMap<Key, Value, Hash, Equal, Alloc>::max_load_factor(float new_max_load_factor) {
    _max_load_factor = new_max_load_factor;
    check_rehash();
}


template <typename Key, typename Value, typename Hash, typename Equal, typename Alloc>
void UnorderedMap<Key, Value, Hash, Equal, Alloc>::rehash(size_t count) {
    List<NodeType*, list_alloc> all_nodes = std::move(nodes);
    restore_map_to_default();
    chunks = count;
    pointers_to_beginning = std::vector<list_iterator, iterator_alloc>(chunks, nodes.end());
    pointers_to_ends = std::vector<list_iterator, iterator_alloc>(chunks, nodes.end());
    for (auto element: all_nodes) {
        insert_node(element);
    }
}


template <typename Key, typename Value, typename Hash, typename Equal, typename Alloc>
void UnorderedMap<Key, Value, Hash, Equal, Alloc>::check_rehash() {
    if (load_factor() > _max_load_factor) {
        reserve(chunks * 2 + 1); // to change the module
    }
}

