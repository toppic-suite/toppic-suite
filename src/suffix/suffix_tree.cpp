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

#include <string>
#include <algorithm>
#include <vector>

#include "base/logger.hpp"

#include "node.hpp"
#include "edge.hpp"
#include "suffix.hpp"
#include "suffix_tree.hpp"

namespace toppic {

namespace suffix {

void SuffixTree::init() {
  root = std::make_shared<Node>(shared_from_this(), nullptr);
  SuffixPtr active = std::make_shared<Suffix>(root, 0, -1);
  LOG_DEBUG("text length " << text_.length());
  for (size_t i = 0; i < text_.length(); i++) {
    addPrefix(active, i);
  }
}

void SuffixTree::addPrefix(SuffixPtr active, int endIndex) {
  NodePtr lastParentNode = nullptr;
  NodePtr parentNode;

  while (true) {
    EdgePtr edge;
    parentNode = active->getOriginNode();

    // Step 1 is to try and find a matching edge for the given node.
    // If a matching edge exists, we are done adding edges, so we break out of this big loop.
    if (active->isExplicit()) {
      edge = active->getOriginNode()->findEdge(text_[endIndex]);
      if (edge != nullptr) break;
    } else {
      // implicit node, a little more complicated
      edge = active->getOriginNode()->findEdge(text_[active->getBeginIndex()]);
      int span = active->getSpan();
      if (text_.at(edge->getBeginIndex() + span + 1) == text_.at(endIndex)) break;
      parentNode = edge->splitEdge(active);
    }

    // We didn't find a matching edge, so we create a new one, add it to the tree at the parent node position,
    // and insert it into the hash table.  When we create a new node, it also means we need to create
    // a suffix link to the new node from the last node we visited.
    // Edge newEdge = new Edge(endIndex, text.length() - 1, parentNode);
    EdgePtr newEdge = std::make_shared<Edge>(endIndex, text_.length() - 1, parentNode, nullptr);
    newEdge->setLeafEdge();
    newEdge->leaf_->setLeafIndex(leaf_created_this_step_);
    //  newEdge.createSuffixPosition(getSeqIndex()+1, getLeafCreatedThisStep()+1);
    newEdge->leaf_->setSeqNum(seq_idx_);
    leaf_created_this_step_++;
    determineSuffixPos();
    newEdge->insert();
    updateSuffixNode(lastParentNode, parentNode);
    lastParentNode = parentNode;

    // This final step is where we move to the next smaller suffix
    if (active->getOriginNode() == root) {
      active->incBeginIndex();
    } else {
      active->changeOriginNode();
    }
    active->canonize();
  }
  updateSuffixNode(lastParentNode, parentNode);
  active->incEndIndex();
  active->canonize();
}

void SuffixTree::updateSuffixNode(NodePtr node, NodePtr suffixNode) {
  if ((node != nullptr) && (node != root)) {
    node->setSuffixNode(suffixNode);
  }
}

void SuffixTree::determineSuffixPos() {
  if (leaf_created_this_step_ == static_cast<int>(database_->getSeqLength(seq_idx_) + 1)) {
    seq_idx_++;
    leaf_created_this_step_ = 0;
  }
}

NodePtr SuffixTree::getRoot() {
  return root;
}

bool comparesp(SuffixPosPtr s1, SuffixPosPtr s2) {
  if (s1->getSeqNum() == s2->getSeqNum()) {
    return s1->getPosInSeq() - s2->getPosInSeq() < 0;
  }
  return s1->getSeqNum() - s2->getSeqNum() < 0;
}

std::vector<SuffixPosPtr> SuffixTree::search(std::string target) {
  std::vector<SuffixPosPtr> startPosList;
  std::vector<NodePtr> stack;
  EdgePtr matchEdge = findMatchEdge(target);
  if (matchEdge == nullptr) {
    return startPosList;
  }
  NodePtr node = matchEdge->getEndNode();
  if (node == nullptr) {
    startPosList.push_back(matchEdge->leaf_->getSuffixPosition());
    return startPosList;
  }

  stack.push_back(node);
  while (stack.size() > 0) {
    std::vector<NodePtr> childNodes;
    for (size_t i = 0; i < stack.size(); i++) {
      NodePtr node2 = stack[i];
      for (int i = 0; i < 21; i++) {
        EdgePtr edge = node2->getEdge(i);
        if (edge == nullptr) {
          continue;
        }
        if (edge->hasEndNode()) {
          childNodes.push_back(edge->getEndNode());
        } else {
          startPosList.push_back(edge->leaf_->getSuffixPosition());
        }
      }
    }

    stack = childNodes;
  }

  std::sort(startPosList.begin(), startPosList.end(), comparesp);
  return startPosList;
}

// find the matched edge for a string to be searched
EdgePtr SuffixTree::findMatchEdge(std::string target) {
  if (target.length() == 0) return nullptr;

  size_t i = 0;
  NodePtr node = root;
  EdgePtr edge = nullptr;
  while (i < target.length()) {
    edge = node->findEdge(target[i]);
    if (edge == nullptr) return nullptr;
    i++;
    for (int j = edge->getBeginIndex() + 1; j <= edge->getEndIndex(); j++) {
      if (i == target.length()) return edge;

      std::string item = std::string(1, edge->getItemAt(j));
      char chr = target[i];
      std::string tmp = std::string(1, chr);
      if (tmp != item) return nullptr;

      i++;
    }

    node = edge->getEndNode();
  }
  return edge;
}

}  // namespace suffix

}  // namespace toppic
