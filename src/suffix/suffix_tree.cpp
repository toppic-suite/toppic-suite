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

#include <sstream>
#include <string>
#include <algorithm>
#include <vector>

#include "node.hpp"
#include "edge.hpp"
#include "leaf_edge.hpp"
#include "suffix.hpp"
#include "suffix_tree.hpp"

namespace prot {
namespace suffix {

SuffixTree::SuffixTree(std::string text) {
  seqIndex = 0;
  leafCreatedThisStep = 0;
  this->text = text;
  root = new Node(this, NULL);
  Suffix *active = new Suffix(root, 0, -1);
  for (size_t i = 0; i < text.length(); i++)
    addPrefix(active, i);
}

SuffixTree::SuffixTree(std::string text, ProteinDatabase *database) {
  seqIndex = 0;
  leafCreatedThisStep = 0;
  this->text = text;
  root = new Node(this, NULL);
  this->database = database;
  Suffix * active = new Suffix(root, 0, -1);
  std::cout << "text.length() " << text.length() << std::endl;
  for (size_t i = 0; i < text.length(); i++) {
    addPrefix(active, i);
  }
}

void SuffixTree::addPrefix(Suffix *active, int endIndex) {
  Node *lastParentNode = NULL;
  Node *parentNode;

  while (true) {
    Edge *edge;
    parentNode = active->getOriginNode();

    // Step 1 is to try and find a matching edge for the given node.
    // If a matching edge exists, we are done adding edges, so we break out of this big loop.
    if (active->isExplicit()) {
      edge = active->getOriginNode()->findEdge(text[endIndex]);
      if (edge != NULL)
        break;
    } else {
      // implicit node, a little more complicated
      edge = active->getOriginNode()->findEdge(text[active->getBeginIndex()]);
      int span = active->getSpan();
      if (text.at(edge->getBeginIndex() + span + 1) == text.at(endIndex))
        break;
      parentNode = edge->splitEdge(active);
    }

    // We didn't find a matching edge, so we create a new one, add it to the tree at the parent node position,
    // and insert it into the hash table.  When we create a new node, it also means we need to create
    // a suffix link to the new node from the last node we visited.
    // Edge newEdge = new Edge(endIndex, text.length() - 1, parentNode);
    Edge *newEdge = new Edge(endIndex, text.length()-1, parentNode, NULL);
    newEdge->setLeafEdge();
    newEdge->leaf->setLeafIndex(leafCreatedThisStep);
    //  newEdge.createSuffixPosition(getSeqIndex()+1, getLeafCreatedThisStep()+1);
    newEdge->leaf->setSeqNum(seqIndex);
    leafCreatedThisStep++;
    determineSuffixPos();
    newEdge->insert();
    updateSuffixNode(lastParentNode, parentNode);
    lastParentNode = parentNode;

    // This final step is where we move to the next smaller suffix
    if (active->getOriginNode() == root)
      active->incBeginIndex();
    else
      active->changeOriginNode();
    active->canonize();
  }
  updateSuffixNode(lastParentNode, parentNode);
  active->incEndIndex();
  active->canonize();
}

void SuffixTree::updateSuffixNode(Node *node, Node *suffixNode) {
  if ((node != NULL) && (node != root)) {
    node->setSuffixNode(suffixNode);
  }
}

void SuffixTree::determineSuffixPos() {
  if (leafCreatedThisStep == static_cast<int>(database->getSeqLength(seqIndex) + 1)) {
    seqIndex++;
    leafCreatedThisStep = 0;
  }
}

Node * SuffixTree::getRoot() {
  return root;
}

bool comparesp(SuffixPosition *s1, SuffixPosition *s2) {
  if (s1->getSeqNum() == s2->getSeqNum())
    return s1->getPosInSeq() - s2->getPosInSeq() < 0;
  return s1->getSeqNum() - s2->getSeqNum() < 0;
}

std::vector<SuffixPosition *> SuffixTree::search(std::string target) {
  std::vector<SuffixPosition *> startPosList;
  std::vector<Node *> stack;
  Edge * matchEdge = findMatchEdge(target);
  if (matchEdge == NULL) {
    return startPosList;
  }
  Node *node = matchEdge->getEndNode();
  if (node == NULL) {
    startPosList.push_back(matchEdge->leaf->getSuffixPosition());
    return startPosList;
  }

  stack.push_back(node);
  while (stack.size() > 0) {
    std::vector<Node*> childNodes;
    for (size_t i = 0; i < stack.size(); i++) {
      Node* node2 = stack[i];
      for (int i = 0; i < 21; i++) {
        Edge *edge = node2->getEdge(i);
        if (edge == NULL)
          continue;
        if (edge->hasEndNode())
          childNodes.push_back(edge->getEndNode());
        else
          startPosList.push_back(edge->leaf->getSuffixPosition());
      }
    }

    stack = childNodes;
  }

  std::sort(startPosList.begin(), startPosList.end(), comparesp);
  return startPosList;
}

// find the matched edge for a string to be searched
Edge *SuffixTree::findMatchEdge(std::string target) {
  if (target.length() == 0)
    return NULL;

  size_t i = 0;
  Node *node = root;
  Edge *edge = NULL;
  while (i < target.length()) {
    edge = node->findEdge(target[i]);
    if (edge == NULL)
      return NULL;
    i++;
    for (int j = edge->getBeginIndex() + 1; j <= edge->getEndIndex(); j++) {
      if (i == target.length())
        return edge;
      std::string item = getstring(edge->getItemAt(j));
      char chr = target[i];
      std::string tmp = getstring(chr);
      if (tmp != item)
        return NULL;
      i++;
    }

    node = edge->getEndNode();
  }
  return edge;
}

std::string SuffixTree::getstring(char c) {
  std::stringstream ss;
  ss<< c;
  return ss.str();
}
}  // namespace suffix
}  // namespace prot
