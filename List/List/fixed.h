#include <cstddef>
#include <vector>


template <size_t chunkSize>
class FixedAllocator {
private:
    struct Chunk {
        std::byte memory{chunkSize};
        Chunk* nextChunk = nullptr;
    };
    
    
    static FixedAllocator<chunkSize>* instance;
    std::vector<Chunk*> blockPool;
    Chunk* nextChunk = nullptr;
    const size_t maxBlockSize = 1024;
    size_t blockSize = 1;
    
    Chunk* makeBlock() {
        Chunk* newBlock = new Chunk[blockSize];
        for (size_t i = 0; i + 1 < blockSize; ++i) {
            (newBlock + i)->nextChunk = newBlock + (i + 1);
        }
        return newBlock;
    }
    
    void updateSize() {
        if (blockSize < maxBlockSize) {
            blockSize *= 2;
        }
    }
    
    void addBlock() {
        updateSize();
        Chunk* newBlock = makeBlock();
        blockPool.push_back(newBlock);
        nextChunk = newBlock;
    }
    
    FixedAllocator() = default;
public:
    static FixedAllocator<chunkSize>* getInstance() {
        if (instance) {
            return instance;
        }
        instance = new FixedAllocator<chunkSize>();
    }
    ~FixedAllocator() {
        for (auto block: blockPool) {
            delete block;
        }
    }
    void* allocate() {
        if (nextChunk == nullptr) {
            addBlock();
        }
        void* newChunk = static_cast<void*>(nextChunk);
        nextChunk = nextChunk->next;
        return newChunk;
    }
    void deallocate(void* chunk) {
        Chunk* chunkToDelte = static_cast<Chunk*>(chunk);
        chunkToDelte->next = nextChunk;
        nextChunk = chunkToDelte;
    }
};
