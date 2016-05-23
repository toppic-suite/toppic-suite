#include "threadpool.hpp"

namespace prot {

template <typename T>
ThreadPool<T>::ThreadPool(int threads, std::string file_name) : 
    terminate(false), stopped(false) {

      for(int i = 0; i < threads; i++) {
        ThreadPtr thread_ptr(new boost::thread(&ThreadPool::Invoke, this));
        threadPool.emplace_back(thread_ptr);
        std::string thread_file_name = file_name + "_" + std::to_string(i);
        std::shared_ptr<T> writer_ptr = std::make_shared<T>(thread_file_name);
        std::pair<boost::thread::id, PrsmXmlWriterPtr> id_writer(thread_ptr->get_id(), writer_ptr);
        writerPool.push_back(id_writer);
      }
    }

template <typename T>
void ThreadPool<T>::Enqueue(std::function<void()> f) {
  // Scope based locking.
  {
    // Put unique lock on task mutex.
    boost::unique_lock<boost::mutex> lock(tasksMutex);

    // Push task into queue.
    tasks.push(f);
  }

  // Wake up one thread.
  condition.notify_one();
}

template <typename T>
void ThreadPool<T>::Invoke() {
  /*
     mtx.lock();
     index ++;
     std::cout << "Thread started index " << index << std::endl; 
     mtx.unlock();
     */
  std::function<void()> task;
  while(true) {
    // Scope based locking.
    {
      // Put unique lock on task mutex.
      boost::unique_lock<boost::mutex> lock(tasksMutex);

      // Wait until queue is not empty or termination signal is sent.
      condition.wait(lock, [this]{ return !tasks.empty() || terminate; });

      // If termination signal received and queue is empty then exit else continue clearing the queue.
      if (terminate && tasks.empty()) {
        return;
      }

      // Get next task in the queue.
      task = tasks.front();

      // Remove it from the queue.
      tasks.pop();
    }

    // Execute the task.
    task();
  }
}

template <typename T>
void ThreadPool<T>::ShutDown() {
  // Scope based locking.
  {
    // Put unique lock on task mutex.
    boost::unique_lock<boost::mutex> lock(tasksMutex);

    // Set termination flag to true.
    terminate = true;
  }

  // Wake up all threads.
  condition.notify_all();

  // Join all threads.
  for(ThreadPtr thread_ptr : threadPool){
    if (thread_ptr->joinable()) thread_ptr->join();
  }

  // Empty workers vector.
  threadPool.empty();

  // Indicate that the pool has been shut down.
  stopped = true;

  for (size_t i = 0; i < writerPool.size(); i++) {
    std::shared_ptr<T> writer_ptr = writerPool[i].second; 
    writer_ptr->close();
  }
}

template <typename T>
std::shared_ptr<T> ThreadPool<T>::getWriter(boost::thread::id thread_id) {
  for (size_t i = 0; i < writerPool.size(); i++) {
    if (writerPool[i].first == thread_id) {
      return writerPool[i].second;
    }
  }
  std::cout << "Thread Error: no writer found" << std::endl;
  return nullptr;
}

// Destructor.
template <typename T>
ThreadPool<T>::~ThreadPool() {
  if (!stopped){
    ShutDown();
  }
}

}
