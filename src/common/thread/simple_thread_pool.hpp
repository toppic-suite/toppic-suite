//Copyright (c) 2014 - 2018, The Trustees of Indiana University.
//
//Licensed under the Apache License, Version 2.0 (the "License");
//you may not use this file except in compliance with the License.
//You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
//Unless required by applicable law or agreed to in writing, software
//distributed under the License is distributed on an "AS IS" BASIS,
//WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//See the License for the specific language governing permissions and
//limitations under the License.

#ifndef TOPPIC_COMMON_THREAD_SIMPLE_THREAD_POOL_HPP_
#define TOPPIC_COMMON_THREAD_SIMPLE_THREAD_POOL_HPP_

#include <vector>
#include <queue>

#include <boost/thread/thread.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/condition_variable.hpp>

namespace toppic {

typedef std::shared_ptr<boost::thread> ThreadPtr;

class ToppicThread {
 public:
  ToppicThread (int id, ThreadPtr thread_ptr):
      id_(id), 
      thread_ptr_(thread_ptr) {}
  int getId () {return id_;}
  ThreadPtr getThreadPtr() {return thread_ptr_;}

 private:
  int id_;
  ThreadPtr thread_ptr_;  
};

typedef std::shared_ptr<ToppicThread> ToppicThreadPtr;
typedef std::vector<ToppicThreadPtr> ToppicThreadPtrVec;

class SimpleThreadPool {
 public:
  // Constructor.
  SimpleThreadPool(int threads);

  // Destructor.
  ~SimpleThreadPool();

  // Adds task to a task queue.
  void Enqueue(std::function<void()> f);

  // Shut down the pool.
  void ShutDown();

  int getQueueSize() {return tasks_.size();}

  int getThreadNum() {return thread_ptr_vec_.size();}

  int getId(boost::thread::id thread_id);

 private:
  // Thread pool storage.
  ToppicThreadPtrVec thread_ptr_vec_;

  // Queue to keep track of incoming tasks.
  std::queue<std::function<void()> > tasks_;

  // Task queue mutex.
  boost::mutex tasks_mutex_;

  // Condition variable.
  boost::condition_variable condition_;

  // Indicates that pool needs to be shut down.
  bool terminate_;

  // Indicates that pool has been terminated.
  bool stopped_;

  // Function that will be invoked by our threads.
  void Invoke();
};

typedef std::shared_ptr<SimpleThreadPool> SimpleThreadPoolPtr;

}  // namespace toppic

#endif
