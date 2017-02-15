// Copyright (c) 2014 - 2017, The Trustees of Indiana University.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
// list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation and/or
// other materials provided with the distribution.
//
// Indiana University provides no reassurances that the source code provided does
// not infringe the patent or any other intellectual property rights of any other
// entity. Indiana University disclaims any liability to any recipient for claims
// brought by any other entity based on infringement of intellectual property
// rights or otherwise.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
// ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.



#ifndef PROT_THREAD_POOL_HPP_
#define PROT_THREAD_POOL_HPP_

#include <vector>
#include <queue>
#include <iostream>
#include "boost/thread/thread.hpp"
#include "boost/thread/mutex.hpp"
#include "boost/thread/condition_variable.hpp"

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

  //boost::mutex mtx;

  //int index = 0;

};

//typedef std::shared_ptr<ThreadPool> ThreadPoolPtr;

}

#include "threadpool_impl.hpp"

#endif
