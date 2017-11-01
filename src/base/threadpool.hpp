//Copyright (c) 2014 - 2017, The Trustees of Indiana University.
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



#ifndef PROT_THREAD_POOL_HPP_
#define PROT_THREAD_POOL_HPP_

#include <vector>
#include <queue>
#include <iostream>
#include <utility>
#include <string>

#include <boost/thread/thread.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/condition_variable.hpp>

#include "prsm/prsm_xml_writer.hpp"

namespace prot {

typedef std::shared_ptr<boost::thread> ThreadPtr;

template <typename T>
class ThreadPool {
 public:
  // Constructor.
  ThreadPool(int threads, std::string output_file_name);

  // Destructor.
  ~ThreadPool();

  // Adds task to a task queue.
  void Enqueue(std::function<void()> f);

  // Shut down the pool.
  void ShutDown();

  int getQueueSize() {return tasks.size();}

  int getThreadNum() {return threadPool.size();}

  std::shared_ptr<T> getWriter(boost::thread::id thread_id);

 private:
  // Thread pool storage.
  std::vector<ThreadPtr> threadPool;

  // prsm writer pool
  std::vector<std::pair<boost::thread::id, std::shared_ptr<T>>> writerPool;

  // Queue to keep track of incoming tasks.
  std::queue<std::function<void()> > tasks;

  // Task queue mutex.
  boost::mutex tasksMutex;

  // Condition variable.
  boost::condition_variable condition;

  // Indicates that pool needs to be shut down.
  bool terminate;

  // Indicates that pool has been terminated.
  bool stopped;

  // Function that will be invoked by our threads.
  void Invoke();
};

}  // namespace prot

#include "threadpool_impl.hpp"

#endif
