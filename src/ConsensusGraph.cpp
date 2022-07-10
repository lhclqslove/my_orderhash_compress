#include "ConsensusGraph.h"
#include "bsc_helper.h"
#include "DirectoryUtils.h"
#include <algorithm>
#include <cassert>
#include <deque>
#include <fstream>
#include <iostream>
#include <stack>
#include <zlib.h>
#include "minimap.h"

Edge::Edge(Node *source, Node *sink, read_t read) : source(source), sink(sink) {
    count = 1;
    // insert into sorted vector
    reads.push_back(read);
}
Edge::Edge(Node *source, Node *sink, read_t read,read_t w) : source(source), sink(sink) {
    count = w;
    // insert into sorted vector
    reads.push_back(read);
}
Edge::Edge(Node *source, Node *sink, std::vector<read_t> &reads)
    : source(source), sink(sink), reads(reads) {
    count = reads.size();
}

void Edge::addRead(read_t read) {
    count++;
    // insert into sorted vector
    reads.insert(std::lower_bound(reads.begin(), reads.end(), read), read);
}
void Edge::addRead(read_t read,size_t w) {
    count+=w;
    // insert into sorted vector
    reads.insert(std::lower_bound(reads.begin(), reads.end(), read), read);
}

Node::Node(const char base) : base(base) {}
size_t Node::getNodeInDegree() const {
    return edgesIn.size();
}
size_t Node::getNodeOutDergree() const {
    return  edgesOut.size();
}
//找回主链的边
Edge *Node::getEdgeTo(Node *n) {
    const auto &it = std::find_if(
        edgesOut.begin(), edgesOut.end(),
        [&](const Edge* p) { return (p->sink == n); });
    if (it == edgesOut.end()) {
        return nullptr;
    } else {
        return *it;
    }
}

Edge *Node::getEdgeToSide(char base) {
    const auto &end = edgesOut.end();
    for (auto it = edgesOut.begin(); it != end; it++) {
        Node *n = (*it)->sink;
        if (!n->onMainPath && n->base == base)
            return *it;
    }
    return nullptr;
}

Edge *Node::getBestEdgeOut() {
    Edge *bestEdge = nullptr;
    read_t bestCount = 0;
    const auto &end = edgesOut.end();
    for (auto it = edgesOut.begin(); it != end; ++it) {
        Edge *e = *it;
        if (e->count > bestCount) {
            bestCount = e->count;
            bestEdge = e;
        }
    }
    return bestEdge;
}

Edge *Node::getBestEdgeIn() {
    Edge *bestEdge = nullptr;
    read_t bestCount = 0;
    const auto &end = edgesIn.end();
    for (auto it = edgesIn.begin(); it != end; ++it) {
        Edge *e = *it;
        if (e->count > bestCount) {
            bestCount = e->count;
            bestEdge = e;
        }
    }
    return bestEdge;
}

Edge *Node::getEdgeInRead(read_t read) const {
    for (auto e : edgesOut) {
        if (std::binary_search(e->reads.begin(), e->reads.end(),
                               read)) {
            return e;
        }
    }
    return nullptr;
}

Node *Node::getNextNodeInRead(read_t read) const {
    Edge *e = getEdgeInRead(read);
    return e ? e->sink : nullptr;
}

double Path::getAverageWeight() {
    size_t totalWeight = 0;
    for (Edge *e : edges)
        totalWeight += e->count;
    return ((double)totalWeight) / edges.size();
}

void Path::clear() {
    size_t l = edges.size();
    if (l > 0) {
        Node *currentNode = edges[0]->source;
        currentNode->onMainPath = false;
    }
    for (size_t i = 0; i < l; i++) {
        Node *currentNode = edges[i]->sink;
        currentNode->onMainPath = false;
    }
    edges.clear();
    path.clear();
}

ConsensusGraphWriter::ConsensusGraphWriter(const std::string &filePrefix) {
    const std::string posFileName = filePrefix + ".pos";
    const std::string ref_id=filePrefix+".refid";
//    const std::string read_id=filePrefix+".readid";
    const std::string editTypeFileName = filePrefix + ".type";
    const std::string editBaseFileName = filePrefix + ".base";
    const std::string idFileName = filePrefix + ".id";
    const std::string complementFileName = filePrefix + ".complement";
    const std::string genomeFileName = filePrefix + ".genome";
    const std::string loneFileName = filePrefix + ".lone";
    const std::string loneidFileName=filePrefix+".loneid";
    posFile.open(posFileName, std::ios::binary);
    refidFile.open(ref_id,std::ios::binary);
//    readidFile.open(read_id,std::ios::binary);
    editTypeFile.open(editTypeFileName);
    editBaseFile.open(editBaseFileName);
    idFile.open(idFileName, std::ios::binary);
    complementFile.open(complementFileName);
//    genomeFile.open(genomeFileName);
    loneidFile.open(loneidFileName, std::ios::binary);
    loneFile.open(loneFileName);
}
void ConsensusGraphWriter::writeLoneReadtoFile() {
    read_t pasId = 0;
    for (auto it : readsLone) {
            read_t diffId = it.first - pasId;
            idFile.write((char*)&diffId, std::ios::binary);
            loneidFile ;
            DirectoryUtils::write_var_uint32(diffId, loneidFile);
            std::string readStr1;
            it.second->read->to_string(readStr1);
            loneFile<<readStr1<<std::endl;
            pasId = it.first;
        }
    readsLone.clear();
}
void ConsensusGraphWriter::closefilestream() {
    posFile.close();
    refidFile.close();
    editTypeFile.close();
    editBaseFile.close();
    idFile.close();
    complementFile.close();
    loneidFile.close();
    loneFile.close();
}
void ConsensusGraph::initialize(const std::string &seed, read_t readId,
                                long pos) {
    // pos is zero here
    size_t len = seed.length();
    Node *currentNode = createNode(seed[0]);
    // We create a read that points to this node
    readsInGraph.insert(std::make_pair(
        readId, ConsensusGraph::Read(pos, currentNode, seed.length(), false)));
    rightMostUnchangedNode = currentNode;
    rightMostUnchangedNodeOffset = 0;
    leftMostUnchangedNode = currentNode;
    leftMostUnchangedNodeOffset = 0;
    mainPath.path.push_back(currentNode->base);
    currentNode->onMainPath = true;
    currentNode->cumulativeWeight = 0;
    for (size_t i = 1; i < len; ++i) {
        Node *nextNode = createNode(seed[i]);
        createEdge(currentNode, nextNode, readId);
        currentNode = nextNode;
    }
    startPos = pos; // = 0
    endPos = pos + 1; 
    // endPos is set to pos+1=1 here because we only inserted one base to mainPath
    // Rest will be inserted later when calculateMainPathGreedy is called
}
std::shared_ptr<my_read::Read> ConsensusGraph::get_mainpathread(size_t cnt)
{
    //原子操作获取新id
    auto  ptr= std::make_shared<my_read::Read>(ReadData::getnewindex());
    ptr->read=std::make_unique<DnaBitset>(mainPath.path.c_str(),mainPath.path.size());
    //沿着走一遍获取base和weight
//#ifdef LOG
//   std::cout<<"len:"<<mainPath.path.size()<<" edge.size"<<mainPath.edges.size()<<std::endl;
//   std::cout<<mainPath.path<<std::endl;
//   bool flag=false;
//#endif

   for (auto &e : mainPath.edges) {
//#ifdef LOG
//        if(e== nullptr)
//        {
//            std::cout<<"遇到空节点"<<std::endl;
//        }
//#endif
        ptr->w.push_back(e->count);

//#ifdef LOG
//       if(!flag)
//       {
//           std::cout<<"成功插入一个"<<std::endl;
//           flag=true;
//       }
//#endif
    }
//#ifdef LOG
//
//    std::cout<<"重新获取权值，插入vector成功"<<std::endl;
//#endif
    ptr->cnt=cnt;
    return  ptr;
}
std::shared_ptr<my_read::Read> ConsensusGraph::get_newsubread(std::string &s){
    auto  ptr= std::make_shared<my_read::Read>(ReadData::getnewindex());
    ptr->read=std::make_unique<DnaBitset>(s.c_str(),s.size());
    ptr->w.resize(s.size()-1,1);
    ptr->cnt=ReadData::recover_cnt;
    return  ptr;
}
void ConsensusGraph::initialize(const std::shared_ptr<my_read::Read> &read, long pos) {
    // pos is zero here
    std::string  seed;
    read->read->to_string(seed);
    size_t len = seed.length();
    Node *currentNode = createNode(seed[0]);
    assert(currentNode);
    _firstReadId=read->id;
    // We create a read that points to this node
    readsInGraph.insert(std::make_pair(
            read->id, ConsensusGraph::Read(pos, currentNode, seed.length(), false)));
#ifdef LOG
    if(currentNode== nullptr)
    {
        std::cout<<"init insert nullptr"<<std::endl;
    }
#endif
    rightMostUnchangedNode = currentNode;
    rightMostUnchangedNodeOffset = 0;
    leftMostUnchangedNode = currentNode;
    leftMostUnchangedNodeOffset = 0;
    mainPath.path.push_back(currentNode->base);
    currentNode->onMainPath = true;
    currentNode->cumulativeWeight = 0;
    for (size_t i = 1; i < len; ++i) {
        Node *nextNode = createNode(seed[i]);
        createEdge(currentNode, nextNode, read->id,read->w[i-1]);
        currentNode = nextNode;
    }
    startPos = pos; // = 0
    endPos = pos + 1;
    // endPos is set to pos+1=1 here because we only inserted one base to mainPath
    // Rest will be inserted later when calculateMainPathGreedy is called
}
bool ConsensusGraph::alignRead(const std::string &s, std::vector<Edit> &editScript, ssize_t &relPos,
            ssize_t &beginOffset, ssize_t &endOffset,ssize_t &match_length, size_t &editdis,size_t &ref_st,size_t &ref_ed,size_t &que_st,size_t &que_ed,size_t m_k, size_t m_w, size_t max_chain_iter) {
    // General comments: 
    // 1. The whole idea behind startPos, endPos and Read.pos in ConsensusGraph
    //    is to provide a reference point for the main path when looking for the next 
    //    read in addRelatedReads in Consensus.cpp. This is needed because we extend the 
    //    mainPath to the left in some cases, and so curPos in Consensus::addRelatedReads
    //    doesn't make sense unless we offset it by startPos which essentially tells us
    //    where the current start position of mainPath is wrt to the start of the first read
    //    in the contig (so startPos is always <=0). endPos is the end of the mainPath wrt
    //    the start of the first read (e.g., it is equal to the length of the first read at 
    //    the very start of the contig when there is only one read). It generally increases 
    //    as we add reads to the right. endPos seems less important for us. Finally, for each 
    //    read in the graph, there is a pos variable (0 for first read) which tells us roughly 
    //    where on the mainPath the read lies (it is actually decided in Consensus::addRelatedReads
    //    itself after the sort-merge procedure). The only use for this pos variable seems 
    //    to be for computation for endPos and startPos in calculateMainPathGreedy after the 
    //    new read is added in. The pos variable is also defined wrt the start of the first read.


    auto &originalString = mainPath.path;
    /// We use either the head or tail of mainPath as a reference to obtain an
    /// offsetGuess

    // TODO: can we remove the Read.pos, startPos, endPos variables in the minimap case 
    //       to help simplify things significantly?
    match_length=0;
    size_t editDis;

    const char* Abegin = originalString.c_str();
    int hits;

    bool success = true;
    //initialize the local buffer
    mm_tbuf_t *b = mm_tbuf_init();
    //initialize the mapopt and iopt
    mm_idxopt_t iopt;
    mm_mapopt_t mopt;
    //0 correpons to map-ont
    mm_set_opt(0, &iopt, &mopt);
    mopt.flag |= MM_F_CIGAR;
    mopt.flag |= MM_F_FOR_ONLY;
    mopt.max_chain_iter = max_chain_iter; 
    // mopt.max_chain_iter - we set this lower because it helps avoid super slow times
    // when we encounter highly repetitive sequences (for whole genome human data).
    // Note: setting it too small leads to worse compression for human datasets

    // only forward alignment, no reverse complement (which is handled elsewhere)
    //call the mm_idx_str to return the index for the reference read    
    // std::cout<<"k:"<<iopt.k<<"w:"<<iopt.w<<std::endl;
    // std::cout<<"flag:"<<iopt.flag<<"bits:"<<iopt.bucket_bits<<hits<<std::endl;      
    //the defalut parameters are: 15 10 false 14
    mm_idx_t * idx = mm_idx_str(m_w, m_k, false, 14, 1, &Abegin, NULL);
    mm_mapopt_update(&mopt, idx);
    //use the index to align with the current read
    //we only want forward matches: use rev in mm_reg1_t
    mm_reg1_t* reg = mm_map(idx, s.length(), s.c_str(), &hits, b, &mopt, NULL);
    editScript.clear();
    //experiment how many hits there are

    if(hits > 0) { 
        //if we have multiple hits, just stick with first hit
        mm_reg1_t *r = &reg[0];
        assert(r->p); 
        //qpos is the current position on the query read
        //rpos is the current position on the reference read
        int qpos = r->qs;
        int rpos = r->rs;
        unsigned int j, k;
        unsigned int count_same;
        int i;
        int alignedLen;

        ref_st=r->rs;
        ref_ed=r->re;
        que_st=r->qs;
        que_ed=r->qe;
    	//add a filtering metric        
    	//calculate the edit distance
    	editDis = r->blen - r->mlen + r->p->n_ambi;
        editdis=editDis;
        //calculate the aligned length; notice that I use the aligned length for the query read here
    	alignedLen = r->qe - r->qs;
        match_length=r->mlen;

        if(alignedLen<s.size()*0.5&& alignedLen<originalString.size()*0.5){
            free(r->p);
            if (hits > 1) {
                // cleanup
                for (int i = 1; i < hits; i++)
                    free(reg[i].p);
            }
//            success=false;
            free(reg);
            mm_tbuf_destroy(b);
            mm_idx_destroy(idx);
            return false;
        }

    	//first check if the read is at the beginnning or the end of reference sequence
        if((r->rs > 0) && (r->re < (ssize_t)originalString.size())){
    		//check editDis/alignedLen
    		// std::cout<<"editDis/alignedLen: "<< editDis/(double)alignedLen<<std::endl;
    		// std::cout<<"(double)alignedLen/s.length() "<< (double)alignedLen/s.length()<<std::endl;    		
            // 1.0 and 0.0: no filter
            if(editDis/(double)alignedLen >= 1.0 || (double)alignedLen/s.length()<=0.0 ){
	    		success = false;
	    		free(r->p);
                if (hits > 1) {
                    // cleanup
                    for (int i = 1; i < hits; i++)
                        free(reg[i].p);
                }    
                free(reg);  
                mm_tbuf_destroy(b);
                mm_idx_destroy(idx);
	    		return false;

	    	}
    	}

        // See comments in ConsensusGraph::updateGraph for its high-level functioning

        // Based on my understanding the correct approach would be something like this:
        // - beginOffset - this represents where the alignment starts on the reference 
        //   (can be +ve or -ve). So if r->rs is +ve then beginOffset is simply r->rs. 
        //   On the other hand if r->rs = 0 (so there is part of read to left), beginOffset
        //   is -r->qs (that's how much read is shifted to left wrt reference).
        // - Now if r->rs is +ve, then we need to add the soft clipped bases as insertions in
        //   editScript. When r->rs = 0 (and hence beginOffset is -ve), the first 
        //   |beginOffset| bases in the read (which are the soft clipped bases) will be added
        //   to the graph in updateGraph directly. So we don't add them to editScript in this 
        //   case.
        // - endOffset - this represents where the alignment ends on the reference wrt end of 
        //   reference. So if r->re < originalString.size() (i.e., the alignment ends 
        //   before the last base in reference), endOffset is -ve and equal to 
        //   (r->re-originalString.size()). If r->re == originalString.size(), the read alignment 
        //   potentially extends beyond the end of reference, and endOffset is +ve and equal to 
        //   (s.length()-r->qe). 
        // - If endOffset is -ve, we need to add the soft-clipped bases at end to the editScript 
        //   as insertions. If endOffset is +ve, these bases will be automatically added in 
        //   the graph in updateGraph function, so we do not add them to editScript.

        // compute relPos between read and reference to keep track of where we are on mainPath
        relPos = (ssize_t)r->rs - (ssize_t)r->qs;

        //update the beginOffset by checking if r->rs is positive or not
        if(r->rs >0){
            beginOffset = r->rs;  
            //add the soft clipped bases as insertions
            for (i = 0; i < r->qs; ++i){
                editScript.push_back(Edit(INSERT, s[i]));             
            }       
        }
        else if(r->rs == 0){
            beginOffset = -r->qs;
        }
        else{
            throw std::runtime_error("Encountered invalid reference start");
        }

        for (j = 0; j < r->p->n_cigar; ++j){ // IMPORTANT: this gives the CIGAR in the aligned regions. NO soft/hard clippings!
            //std::cout<< (r->p->cigar[j]>>4) << "MIDNSH"[r->p->cigar[j]&0xf];
            switch("MIDNSH"[r->p->cigar[j]&0xf]) {
                case 'M':
                    //fix: handle the substitute and same differently
                    count_same = 0;
                    for(k=0; k<r->p->cigar[j]>>4;k++){
                        if(s[qpos] == Abegin[rpos]){
                            count_same++;
                        }
                        else{
                            if (count_same > 0)
                                editScript.push_back(Edit(SAME, count_same));
                            count_same = 0; // reset count_same
                            //add the substitute as one insert and one delete
                            editScript.push_back(Edit(DELETE,Abegin[rpos]));
                            editScript.push_back(Edit(INSERT,s[qpos]));                             
                        }
                        qpos++;
                        rpos++;
                    }
                    //check if we need to push the "same" again
                    if(count_same!=0){
                        editScript.push_back(Edit(SAME, count_same));                       
                    }
                    break; 
                case 'I':
                    for(k=0; k<r->p->cigar[j]>>4;k++){
                        editScript.push_back(Edit(INSERT, s[qpos])); 
                        qpos++; 
                    }                   
                    break; 
                case 'D':
                    for(k=0; k<r->p->cigar[j]>>4;k++){
                        editScript.push_back(Edit(DELETE, Abegin[rpos])); 
                        rpos++; 
                    }                       
                    break; 
                default: 
                    throw std::runtime_error("Encountered invalid CIGAR symbol!");
            }           
        }

        //update the endOffset by checking if r->re is positive or not
        if(r->re < (ssize_t)originalString.size()){
            endOffset = (r->re-originalString.size());
            //add the soft clipped bases as insertions
            for (i = r->qe; i < (ssize_t)s.length(); ++i){
                editScript.push_back(Edit(INSERT, s[i]));             
            }       
        }
        else if(r->re == (ssize_t)originalString.size()){
            endOffset = (s.length()-r->qe);             
        }
        else{
            throw std::runtime_error("Encountered invalid reference end");
        }

//#ifdef LOG
//        std::cout<< "editDis: " << editDis<<std::endl;
//#endif
        free(r->p);
        if (hits > 1) {
            // cleanup
            for (int i = 1; i < hits; i++)
                free(reg[i].p);
            }
    }
    else{
        //return fail when there are not hits
        success = false;
    }
    //return the correct beginOffset and endOffset   
    editScript.shrink_to_fit();
    free(reg);  
    mm_tbuf_destroy(b);
    mm_idx_destroy(idx);       
    //number of hits
    //edit distance as threshold
    //length of alignment as fraction
    //DP alignment score
    //    std::cout << "success ? " << success << std::endl;
    if (!success) {
        // std::cout << "Failed to add"
        //           << "\n";
        // std::cout << "mainPath startPos " << startPos << " mainPath endPos "
        //           << endPos << " current read startPos " << pos
        //           << " OffsetGuess " << offsetGuess << std::endl;
        return false;
    }
    size_t numUnchanged = 0;
    for (auto e : editScript)
        if (e.editType == SAME)
            numUnchanged += e.editInfo.num;
    if (numUnchanged == 0)
        return false;
    return true;
}

void ConsensusGraph::updateGraph(const std::string &s,
                                 std::vector<Edit> &editScript,
                                 ssize_t beginOffset, ssize_t endOffset,
                                 read_t readId, long pos,
                                 bool reverseComplement) {
    // How does the updateGraph function work (high level)?
    // 1. If beginOffset is -ve, create |beginOffset| new nodes in the graph 
    //    joined in a chain representing those bases in read 
    //    (the chain will be essentially joined to position 0 of mainPath in graph).
    // 2. Then begin adding nodes & edges to graph according to the editScript 
    //    at the relevant place in the graph as determined by beginOffset 
    //    (essentially starting at position 0 of mainPath when beginOffset is -ve)
    // 3. Finally, if endOffset is +ve, create |endOffset| new nodes in the graph 
    //    joined in a chain representing those bases in read (the chain begins 
    //    in the graph where step 2 ends, which will be roughly the end of mainPath 
    //    in this case). 

    const auto &edgeInPathEnd = mainPath.edges.end();
    auto edgeInPath = mainPath.edges.begin();
    /** nodeInPath should always be set to edgeInPath.source, unless it is the
    last node, in which case edgeInPath should be edgeInPathEnd **/
    Node *nodeInPath = (*edgeInPath)->source;
    /** The last node of this read that has been added **/
    Node *currentNode = nullptr;
    Node *initialNode = nullptr;

    size_t numUnchanged = 0;

    // First we update leftMostUnchangedNode and rightMostUnchangedNode
    if (beginOffset >= 0 || endOffset >= 0) {
        rightMostUnchangedNodeOffset = static_cast<size_t>(std::max(
            static_cast<ssize_t>(leftMostUnchangedNodeOffset),
            std::min(static_cast<ssize_t>(rightMostUnchangedNodeOffset),
                     beginOffset)));
        if (rightMostUnchangedNodeOffset > 0)
            rightMostUnchangedNode =
                mainPath.edges[rightMostUnchangedNodeOffset - 1]->sink;
        else
            rightMostUnchangedNode = mainPath.edges[0]->source;
    } else {
        leftMostUnchangedNodeOffset =
            std::min(rightMostUnchangedNodeOffset,
                     std::max(leftMostUnchangedNodeOffset,
                              mainPath.path.size() - 1 + endOffset));
        if (leftMostUnchangedNodeOffset > 0)
            leftMostUnchangedNode =
                mainPath.edges[leftMostUnchangedNodeOffset - 1]->sink;
        else
            leftMostUnchangedNode = mainPath.edges[0]->source;
    }

    auto advanceNodeInPath = [&]() {
        if (edgeInPath == edgeInPathEnd)
            return;
        nodeInPath = (*edgeInPath)->sink;
        edgeInPath++;
    };

    // First we deal with beginOffset
    auto initialAdvance = [&] {
        //        std::cout << "initialAdvance" << std::endl;
        if (beginOffset >= 1) {
            // We need to advance in the mainPath
            // for (size_t i = 0; i < beginOffset; i++) {
            //     advanceNodeInPath();
            // }
            edgeInPath += beginOffset - 1;
            nodeInPath = (*edgeInPath)->sink;
            edgeInPath++;

        } else if (beginOffset <= -1) {
            // We need to insert the initial parts of the read
            // This number must be positive
            size_t numOfNodes2Insert = -beginOffset;
            // std::cout << numOfNodes2Insert << " ";
            size_t i = 0;
            // We create an initial node
            currentNode = createNode(s[i++]);
            initialNode = currentNode;
            for (; i < numOfNodes2Insert; ++i) {
                Node *nextNode = createNode(s[i]);
                createEdge(currentNode, nextNode, readId);
                currentNode = nextNode;
            }
        }
    };

    initialAdvance();

    auto insertNode = [&](char base) {
        if (!currentNode) {
            // If there is no currentNode, we create one
            currentNode = createNode(base);
            initialNode = currentNode;
        } else {
            Edge *edge = currentNode->getEdgeToSide(base);
            if (edge) {
                edge->addRead(readId);
            } else {
                Node *n = createNode(base);
                edge = createEdge(currentNode, n, readId);
            }
            currentNode = edge->sink;
        }
    };


    // Now we deal with the edits one by one
    for (const Edit &e : editScript) {
        if (e.editType == SAME) {
            size_t num = e.editInfo.num;
            numUnchanged += num;
            // We deal with the first node
            if (!currentNode) {
                // If there is no currentNode, set it and initialNode
                initialNode = nodeInPath;
                currentNode = nodeInPath;
            } else {
                // Otherwise create an edge from the currentNode to nodeInPath
                Edge *edge = currentNode->getEdgeTo(nodeInPath);
                if (edge) {
                    // If there is already an edge to it
                    edge->addRead(readId);
                } else {
                    // Otherwise create an edge
                    edge = createEdge(currentNode, nodeInPath, readId);
                }
                currentNode = nodeInPath;
            }
            advanceNodeInPath();

            // We deal with the rest of the nodes
            for (size_t i = 1; i < num; ++i) {
                currentNode->getEdgeTo(nodeInPath)->addRead(readId);
                currentNode = nodeInPath;
                advanceNodeInPath();
            }
        } else if (e.editType == DELETE) {
            // Deletion just means advancing the node in path
            advanceNodeInPath();
        } else if (e.editType == INSERT) {
            char base = e.editInfo.ins;
            insertNode(base);
        }
    }

    // Finally we deal with endOffset
    if (endOffset > 0) {
        auto end = s.end();
        for (auto it = end - endOffset; it < end; it++) {
            insertNode(*it);
        }
    }
    assert(numUnchanged > 0);
    // Don't forget to add the read!
    readsInGraph.insert(std::make_pair(
        readId, Read(pos, initialNode, s.length(), reverseComplement)));
}
void ConsensusGraph::updateGraph(std::shared_ptr<my_read::Read> &read, std::vector<Edit> &editScript, ssize_t beginOffset,
                                 ssize_t endOffset, read_t readId, long pos, bool reverseComplement,size_t ref_st,size_t ref_ed,size_t que_st,size_t que_ed) {
    // How does the updateGraph function work (high level)?
    // 1. If beginOffset is -ve, create |beginOffset| new nodes in the graph
    //    joined in a chain representing those bases in read
    //    (the chain will be essentially joined to position 0 of mainPath in graph).
    // 2. Then begin adding nodes & edges to graph according to the editScript
    //    at the relevant place in the graph as determined by beginOffset
    //    (essentially starting at position 0 of mainPath when beginOffset is -ve)
    // 3. Finally, if endOffset is +ve, create |endOffset| new nodes in the graph
    //    joined in a chain representing those bases in read (the chain begins
    //    in the graph where step 2 ends, which will be roughly the end of mainPath
    //    in this case).

    this->ref_st=ref_st;
    this->ref_ed=ref_ed;
    this->que_st=que_st;
    this->que_st=que_ed;
    _secondReadId=read->id;
//#ifdef LOG
//    std::cout<<"firstReadId:"<<_firstReadId<<" secondReadId:"<<_secondReadId<<std::endl;
//    std::cout<<mainPath.path<<std::endl;
//    std::string readStr2;
//    read->read->to_string(readStr2);
//    std::cout<<readStr2<<std::endl;
//#endif
    std::string  s;
    std::string readStr1;
    read->read->to_string(readStr1);
    this->que_len=readStr1.size();
    if (reverseComplement)
        ReadData::toReverseComplement(
                readStr1.begin(), readStr1.end(),
                std::inserter(s, s.end()));
    else
        s = readStr1;
//#ifdef  LOG
//    std::string yuanstr=this->mainPath.path;
//    size_t siz=editScript.size();
//    std::cout<<"s:len "<<s.size()<<" "<<s<<std::endl;
//    std::cout<<"mainpath len:"<< this->mainPath.path.size()<<" "<< this->mainPath.path<<std::endl;
//////    std::cout<<"delbeg"<<this->mainPath.path.substr(beginOffset)<<std::endl;
//    std::cout<<" editScript.size::"<<editScript.size()<<std::endl;
////    std::cout<<"beginOffset:"<<beginOffset<<" endOffset:"<<endOffset<<" pos:"<<pos<<std::endl;
////    std::cout<<"reverse:"<<reverseComplement<<std::endl;
////    if(editScript[0].editType==EDIT_TYPE::INSERT)
////    {
////        std::cout<<"insert"<<" "<<editScript[0].editInfo.ins<<std::endl;
////    }
////    else if(editScript[0].editType==EDIT_TYPE::SAME)
////    {
////        std::cout<<"same"<<" "<<editScript[0].editInfo.num<<std::endl;
////    }
////    for(auto &edit:editScript)
////    {
////        std::cout<<edit.editType<<std::endl;
////    }
//
//#endif
    const auto &edgeInPathEnd = mainPath.edges.end();
    auto edgeInPath = mainPath.edges.begin();
    /** nodeInPath should always be set to edgeInPath.source, unless it is the
    last node, in which case edgeInPath should be edgeInPathEnd **/
    Node *nodeInPath = (*edgeInPath)->source;
    /** The last node of this read that has been added **/
    Node *currentNode = nullptr;
    Node *initialNode = nullptr;
    size_t qpos=0;
    size_t numUnchanged = 0;

    // First we update leftMostUnchangedNode and rightMostUnchangedNode
    if (beginOffset >= 0 || endOffset >= 0) {
        rightMostUnchangedNodeOffset = static_cast<size_t>(std::max(
                static_cast<ssize_t>(leftMostUnchangedNodeOffset),
                std::min(static_cast<ssize_t>(rightMostUnchangedNodeOffset),
                         beginOffset)));
        if (rightMostUnchangedNodeOffset > 0)
            rightMostUnchangedNode =
                    mainPath.edges[rightMostUnchangedNodeOffset - 1]->sink;
        else
            rightMostUnchangedNode = mainPath.edges[0]->source;
    } else {
        leftMostUnchangedNodeOffset =
                std::min(rightMostUnchangedNodeOffset,
                         std::max(leftMostUnchangedNodeOffset,
                                  mainPath.path.size() - 1 + endOffset));
        if (leftMostUnchangedNodeOffset > 0)
            leftMostUnchangedNode =
                    mainPath.edges[leftMostUnchangedNodeOffset - 1]->sink;
        else
            leftMostUnchangedNode = mainPath.edges[0]->source;
    }
//#ifdef LOG
//    std::cout<<leftMostUnchangedNodeOffset<<"offset  "<<rightMostUnchangedNodeOffset<<std::endl;
//#endif
    auto advanceNodeInPath = [&]() {
        if (edgeInPath == edgeInPathEnd)
            return;
        nodeInPath = (*edgeInPath)->sink;
        edgeInPath++;
    };

    // First we deal with beginOffset
    auto initialAdvance = [&] {
        //        std::cout << "initialAdvance" << std::endl;
        if (beginOffset >= 1) {
            // We need to advance in the mainPath
            // for (size_t i = 0; i < beginOffset; i++) {
            //     advanceNodeInPath();
            // }
            edgeInPath += beginOffset - 1;
            nodeInPath = (*edgeInPath)->sink;
            edgeInPath++;

        } else if (beginOffset <= -1) {
            // We need to insert the initial parts of the read
            // This number must be positive
            size_t numOfNodes2Insert = -beginOffset;
            // std::cout << numOfNodes2Insert << " ";
//            size_t i = 0;
            // We create an initial node
            currentNode = createNode(s[qpos++]);
            initialNode = currentNode;
//            qpos++;
            for (; qpos < numOfNodes2Insert; ++qpos) {
                Node *nextNode = createNode(s[qpos]);
                createEdge(currentNode, nextNode,  read->id,read->w[qpos-1]);
                currentNode = nextNode;
            }
        }
    };

    initialAdvance();
//#ifdef LOG
//    std::cout<<"qpos:"<<qpos<<std::endl;
//#endif
    auto insertNode = [&](char base) {
        if (!currentNode) {
            // If there is no currentNode, we create one
            currentNode = createNode(base);
            initialNode = currentNode;
//#ifdef LOG
//            if( initialNode== nullptr)
//            {
//                std::cout<<"%%%instert nullptr"<<std::endl;
//
//            }
//#endif
        } else {
            Edge *edge = currentNode->getEdgeToSide(base);
            if (edge) {
                edge->addRead(readId,read->w[(qpos<read->w.size()?read->w[qpos]:read->w.back())]);
            } else {
                Node *n = createNode(base);
                edge =  createEdge(currentNode, n,  read->id,(qpos<read->w.size()?read->w[qpos]:read->w.back()));
            }
            qpos++;
            currentNode = edge->sink;
        }
    };
#ifdef LOG
//    std::cout<<"come to here"<<std::endl;
//    std::cout<<editScript.size()<<std::endl;
//    std::cout<<read->w.size()<<std::endl;
#endif
    // Now we deal with the edits one by one
    for (const Edit &e : editScript) {
        if (e.editType == SAME) {
            size_t num = e.editInfo.num;
            numUnchanged += num;
            // We deal with the first node
            if (!currentNode) {
                // If there is no currentNode, set it and initialNode
                initialNode = nodeInPath;
                currentNode = nodeInPath;
//#ifdef LOG
//                if( initialNode== nullptr)
//                {
//                    std::cout<<"!!!!instert nullptr"<<std::endl;
//
//                }
//#endif
            } else {
                // Otherwise create an edge from the currentNode to nodeInPath
                Edge *edge = currentNode->getEdgeTo(nodeInPath);
                if (edge) {
                    // If there is already an edge to it
                    edge->addRead(readId,read->w[(qpos<read->w.size()?read->w[qpos]:read->w.back())]);
                } else {
                    // Otherwise create an edge
                    edge = createEdge(currentNode, nodeInPath, readId,read->w[(qpos<read->w.size()?read->w[qpos]:read->w.back())]);
                }
                currentNode = nodeInPath;
                qpos++;
            }
            advanceNodeInPath();
            // We deal with the rest of the nodes
            for (size_t i = 1; i < num; ++i) {

                currentNode->getEdgeTo(nodeInPath)->addRead(readId,read->w[(qpos<read->w.size()?read->w[qpos]:read->w.back())]);
                currentNode = nodeInPath;
                advanceNodeInPath();
                qpos++;
            }
        } else if (e.editType == DELETE) {
            // Deletion just means advancing the node in path
            advanceNodeInPath();
        } else if (e.editType == INSERT) {

            char base = e.editInfo.ins;
            insertNode(base);
        }
    }
//#ifdef LOG
//    std::cout<<"done for editScript"<<std::endl;
//#endif
    // Finally we deal with endOffset
    if (endOffset > 0) {
        auto end = s.end();
        for (auto it = end - endOffset; it < end; it++) {
            insertNode(*it);
        }
    }
    assert(numUnchanged > 0);
//#ifdef LOG
//    std::cout<<"before^^^^init"<<std::endl;
//#endif
//    assert(initialNode);
//#ifdef LOG
//    std::cout<<"after^^^^init"<<std::endl;
//#endif
    // Don't forget to add the read!
    readsInGraph.insert(std::make_pair(
            readId, Read(pos, initialNode, s.length(), reverseComplement)));
//#ifdef LOG
//    if( initialNode== nullptr)
//    {
//        std::cout<<"instert nullptr"<<std::endl;
//
//    }
//#endif
}

Path &ConsensusGraph::calculateMainPathGreedy() {
    clearMainPath();

    auto &edgesInPath = mainPath.edges;
    auto &stringPath = mainPath.path;

    // Extend to the right
    {
        Node *currentNode = rightMostUnchangedNode;
        assert(currentNode->onMainPath);
        Edge *edgeToAdd;
        while ((edgeToAdd = currentNode->getBestEdgeOut())) {
            edgesInPath.push_back(edgeToAdd);
            currentNode = edgeToAdd->sink;
            currentNode->onMainPath = true;
            stringPath.push_back(currentNode->base);
        }
        read_t endingReadId = *edgesInPath.back()->reads.begin();

        // TODO: why take the begin (i.e., first read through this last edge in path)?
        // Is this a relic of the constant read length code?
        // Doesn't matter since typically only one read at the last edge
        Read &endingRead = readsInGraph.at(endingReadId);
        endPos = endingRead.pos + endingRead.len;
    }

    // Extend to the left
    // NOTE: I have changed mainPath.path to string to simplify alignment,
    // but that means the code below is not efficient (inserting to start
    // of vector). Fortunately, this is currently a very small contributor
    // to the total time. But we might want to fix this if issues crop up later.
    {
        Node *currentNode = leftMostUnchangedNode;
        assert(currentNode->onMainPath);
        Edge *edgeToAdd;
        while ((edgeToAdd = currentNode->getBestEdgeIn())) {
            edgesInPath.insert(edgesInPath.begin(), edgeToAdd);
            currentNode = edgeToAdd->source;


            currentNode->onMainPath = true;
            stringPath.insert(stringPath.begin(), currentNode->base);
            leftMostUnchangedNodeOffset++;
            rightMostUnchangedNodeOffset++;
        }
        read_t startingReadId = *edgesInPath.front()->reads.begin();
        startPos = readsInGraph.at(startingReadId).pos;
    }

    // rightMostUnchangedNodeOffset = 0;
    // rightMostUnchangedNode = edgesInPath.front()->source;
    // printStatus();
    removeCycles();
    rightMostUnchangedNode = edgesInPath.back()->sink;
    rightMostUnchangedNodeOffset = edgesInPath.size();
    leftMostUnchangedNode = edgesInPath.front()->source;
    leftMostUnchangedNodeOffset = 0;
    return mainPath;
}
Path &ConsensusGraph::calculateMainPath() {
    /**
     * 首先找到第一个入度为二的节点和最后一个出度为二的节点
     * 1如果找不到这两个节点，说明只有一条链，用原来贪心地方法更新
     *
     * 2如果这两个节点相同，说明两条序列只有一个base相同，基本不可能，情况不存在
     *
     * 3如果两个节点只有一个节点存在
     *
     * 3.1只有左边存在
     *
     * 3.2只有右边存在
     *
     * 4两个节点都存在
     *
     * 4.1第一个入度为2的节点在最后一个出度为二的节点前面，
     *
     * 4.2最后一个出度为2的节点在第一个入度为2的节点前面。(其中一条序列完全被另外一条序列包含了)
     *
     *
     *
     */

     Node* FirstJionNode= nullptr;
     Node* FinalJionNode= nullptr;
     Node* currentNode=mainPath.edges.back()->sink;
     Node* NextNode= nullptr;
     Edge* NextEdge= nullptr;
     ssize_t las_index=mainPath.edges.size()-1;
     ssize_t first_index=0;
     auto &stringPath = mainPath.path;
     auto &edgesInPath = mainPath.edges;
     auto check=[&](Node* ptr){
         std::set<Node *>st;
         bool flag=false;
         while(ptr!= nullptr)
         {
//             ptr->onMainPath=true;
             if(st.find(ptr)!=st.end())
             {
                std::cout<<"图中有环"<<std::endl;
                 break;
             }
             st.insert(ptr);
             if(ptr==FinalJionNode)
             {
                 flag=true;
                 std::cout<<"能贪心地走到最后相交的节点"<<std::endl;
             }

             auto  edgeptr=ptr->getBestEdgeOut();
             if(edgeptr== nullptr)
             {

                 break;
             }
             ptr=edgeptr->sink;
         }
         if(!flag)
         {
             std::cout<<"已经遍历完全程,但是未到最后的节点"<<std::endl;
         }
         else
         {
             std::cout<<"能贪心地走到最后相交的节点  成功"<<std::endl;
         }
         return  flag;
     };

     while(las_index>=0)
     {
//         if(mainPath.edges[las_index]->sink== nullptr)
//         {
//             std::cout<<"遇到空节点了"<<std::endl;
//         }
//#ifdef LOG
//         assert(mainPath.edges[las_index]->sink);
////         std::cout<<"las_index"<<las_index<<std::endl;
//#endif
         if(mainPath.edges[las_index]->sink->getNodeOutDergree()<2){
             las_index--;
         }
         else{
             FinalJionNode=mainPath.edges[las_index]->sink;
             break;
         }
     }
     while(first_index<mainPath.edges.size())
     {
         if(mainPath.edges[first_index]->source->getNodeInDegree()<2){
             first_index++;
         }
         else
         {
             FirstJionNode=mainPath.edges[first_index]->source;
             break;
         }
     }

     auto clearFlag=[&](){
         for(size_t i=0;i<mainPath.edges.size();i++){
             mainPath.edges[i]->source->onMainPath= false;
         }
         mainPath.edges.back()->sink->onMainPath=false;
     };

     //1 2 4.2
     if((FirstJionNode== nullptr&&FinalJionNode== nullptr)||first_index>=las_index){
//#ifdef LOG
//         std::cout<<"case 1 2 4.2"<<std::endl;
//#endif
         clearMainPath();

         auto &edgesInPath = mainPath.edges;
         auto &stringPath = mainPath.path;

         // Extend to the right
         {
             Node *currentNode = rightMostUnchangedNode;
             assert(currentNode->onMainPath);
             Edge *edgeToAdd;
             while ((edgeToAdd = currentNode->getBestEdgeOut())) {
                 edgesInPath.push_back(edgeToAdd);
                 currentNode = edgeToAdd->sink;
                 currentNode->onMainPath = true;
                 stringPath.push_back(currentNode->base);
             }
             read_t endingReadId = *edgesInPath.back()->reads.begin();

             // TODO: why take the begin (i.e., first read through this last edge in path)?
             // Is this a relic of the constant read length code?
             // Doesn't matter since typically only one read at the last edge
             Read &endingRead = readsInGraph.at(endingReadId);
             endPos = endingRead.pos + endingRead.len;
         }

         // Extend to the left
         // NOTE: I have changed mainPath.path to string to simplify alignment,
         // but that means the code below is not efficient (inserting to start
         // of vector). Fortunately, this is currently a very small contributor
         // to the total time. But we might want to fix this if issues crop up later.
         {
             Node *currentNode = leftMostUnchangedNode;
             assert(currentNode->onMainPath);
             Edge *edgeToAdd;
             while ((edgeToAdd = currentNode->getBestEdgeIn())) {
                 edgesInPath.insert(edgesInPath.begin(), edgeToAdd);
                 currentNode = edgeToAdd->source;


                 currentNode->onMainPath = true;
                 stringPath.insert(stringPath.begin(), currentNode->base);
                 leftMostUnchangedNodeOffset++;
                 rightMostUnchangedNodeOffset++;
             }
             read_t startingReadId = *edgesInPath.front()->reads.begin();
             startPos = readsInGraph.at(startingReadId).pos;
         }
         removeCycles();
         rightMostUnchangedNode = edgesInPath.back()->sink;
         rightMostUnchangedNodeOffset = edgesInPath.size();
         leftMostUnchangedNode = edgesInPath.front()->source;
         leftMostUnchangedNodeOffset = 0;

         return mainPath;
     }
     else if(FirstJionNode!= nullptr&&FinalJionNode!= nullptr){//4.1
//#ifdef LOG
//         std::cout<<"case 4.1"<<std::endl;
//#endif
         //先清空原有标记
         clearFlag();
//#ifdef LOG
//         std::cout<<"清空标记成功"<<std::endl;
//#endif
         //mainPath.edges.clear();
         //判断左右分支长度,选择正确的的起点，以及正确的后继节点
         if(ref_st<que_st)currentNode=readsInGraph[_secondReadId].start;//选第二条序列的起点
         else currentNode=readsInGraph[_firstReadId].start;//选第一条的起点
         if(mainPath.path.size()-ref_ed>que_len-que_ed)NextNode=mainPath.edges[las_index+1]->sink,NextEdge=mainPath.edges[las_index+1];//选主链为后续
         else{
             for(auto &edge :FinalJionNode->edgesOut)
             {
                 if(edge!=mainPath.edges[las_index+1]){
                     NextNode=edge->sink;
                     NextEdge=edge;
                     break;
                 }
             }
         }
         assert(currentNode);
         assert(NextNode);
         assert(NextEdge);
//#ifdef LOG
//         std::cout<<"初始化 currentNode NextNode成功"<<std::endl;
//#endif
         //清空序列，和duque,然后重新插入
         mainPath.path.clear();
         mainPath.edges.clear();
         //从头开始插入到finaljionNode
         Edge *edgeToAdd;
//#ifdef LOG
//         std::cout<<"开始check"<<std::endl;
//        check(currentNode);
//#endif
         while(currentNode!=FinalJionNode)
         {
             currentNode->onMainPath=true;
             stringPath.push_back(currentNode->base);
             edgeToAdd=currentNode->getBestEdgeOut();
             if(edgeToAdd== nullptr)
             {
                 std::cout<<"到最后节点前遇到null 节点"<<std::endl;
                 break;
             }
             edgesInPath.push_back(edgeToAdd);
             currentNode=edgeToAdd->sink;
//#ifdef LOG
//             if(currentNode== nullptr)
//             {
//                 std::cout<<"还没到fianl 节点，就遇到Nullptr"<<std::endl;
//             }
//#endif
         }
//#ifdef LOG
//             std::cout<<"能走到FInal节点"<<std::endl;
//#endif
         FinalJionNode->onMainPath= true;
         stringPath.push_back( FinalJionNode->base);
         edgesInPath.push_back(NextEdge);
         currentNode=NextNode;
         while(currentNode!= nullptr){
             currentNode->onMainPath=true;
             stringPath.push_back(currentNode->base);
             edgeToAdd=currentNode->getBestEdgeOut();
             if(edgeToAdd== nullptr)break;
             edgesInPath.push_back(edgeToAdd);
             currentNode=edgeToAdd->sink;
         }
     }
     else if(FinalJionNode== nullptr&&FirstJionNode!= nullptr){//3.1
//#ifdef LOG
//         std::cout<<"case 3.1"<<std::endl;
//#endif
         clearFlag();
         //判断左边分支长度,选择正确的的起点，
         if(ref_st<que_st)currentNode=readsInGraph[_secondReadId].start;//选第二条序列的起点
         else currentNode=readsInGraph[_firstReadId].start;//选第一条的起点
         assert(currentNode);
         //清空序列，和duque,然后重新插入
         mainPath.path.clear();
         mainPath.edges.clear();
         Edge *edgeToAdd;
         while(currentNode!= nullptr){
             currentNode->onMainPath=true;
             stringPath.push_back(currentNode->base);
             edgeToAdd=currentNode->getBestEdgeOut();
             if(edgeToAdd== nullptr)break;
             edgesInPath.push_back(edgeToAdd);
             currentNode=edgeToAdd->sink;
         }
     }
     else if(FinalJionNode!= nullptr&&FirstJionNode== nullptr){//3.2
//#ifdef LOG
//         std::cout<<"case 3.2"<<std::endl;
//#endif
         clearFlag();
         currentNode=readsInGraph[_firstReadId].start;
         if(mainPath.path.size()-ref_ed>que_len-que_ed)NextNode=mainPath.edges[las_index+1]->sink,NextEdge=mainPath.edges[las_index+1];//选主链为后续
         else{
             for(auto &edge :FinalJionNode->edgesOut)
             {
                 if(edge!=mainPath.edges[las_index+1]){
                     NextNode=edge->sink;
                     NextEdge=edge;
                     break;
                 }
             }
         }
         assert(currentNode);
         assert(NextNode);
         assert(NextEdge);
         //清空序列，和duque,然后重新插入
         mainPath.path.clear();
         mainPath.edges.clear();
         //从头开始插入到finaljionNode
         Edge *edgeToAdd;
//#ifdef LOG
//         check(currentNode);
//#endif
         while(currentNode!=FinalJionNode)
         {
             currentNode->onMainPath=true;
             stringPath.push_back(currentNode->base);
             edgeToAdd=currentNode->getBestEdgeOut();
             if(edgeToAdd== nullptr)
             {
#ifdef LOG
                 std::cout<<"还没到fianl 节点，就遇到Nullptr"<<std::endl;
#endif
                 break;
             }
             edgesInPath.push_back(edgeToAdd);
             currentNode=edgeToAdd->sink;
         }
//#ifdef LOG
//         std::cout<<"能走到FInal节点"<<std::endl;
//#endif
         FinalJionNode->onMainPath= true;
         stringPath.push_back( FinalJionNode->base);
         edgesInPath.push_back(NextEdge);
         currentNode=NextNode;
         while(currentNode!= nullptr){
             currentNode->onMainPath=true;
             stringPath.push_back(currentNode->base);
             edgeToAdd=currentNode->getBestEdgeOut();
             if(edgeToAdd== nullptr)break;
             edgesInPath.push_back(edgeToAdd);
             currentNode=edgeToAdd->sink;
         }
     }

    // rightMostUnchangedNodeOffset = 0;
    // rightMostUnchangedNode = edgesInPath.front()->source;
    // printStatus();

//    removeCycles();
     rightMostUnchangedNode = edgesInPath.back()->sink;
    rightMostUnchangedNodeOffset = edgesInPath.size();
    leftMostUnchangedNode = edgesInPath.front()->source;
    leftMostUnchangedNodeOffset = 0;

    return mainPath;
}

void ConsensusGraph::clearMainPath() {
    // printStatus();
    // std::cout << "clearing\n";
    size_t l = mainPath.edges.size();
    // We first erase the right tail
    for (size_t i = rightMostUnchangedNodeOffset; i < l; ++i) {
        Node *currentNode = mainPath.edges[i]->sink;
        currentNode->onMainPath = false;
    }
    if (rightMostUnchangedNodeOffset < mainPath.edges.size()) {
        mainPath.edges.erase(mainPath.edges.begin() +
                                 rightMostUnchangedNodeOffset,
                             mainPath.edges.end());
    }
    if (mainPath.path.size() > rightMostUnchangedNodeOffset + 1)
        mainPath.path.erase(mainPath.path.begin() +
                                rightMostUnchangedNodeOffset + 1,
                            mainPath.path.end());
    // Now we erase the left tail
    for (size_t i = 0; i < leftMostUnchangedNodeOffset; ++i) {
        Node *currentNode = mainPath.edges[i]->source;
        currentNode->onMainPath = false;
    }
    if (leftMostUnchangedNodeOffset > 0) {
        mainPath.edges.erase(mainPath.edges.begin(),
                             mainPath.edges.begin() +
                                 leftMostUnchangedNodeOffset);
        mainPath.path.erase(mainPath.path.begin(),
                            mainPath.path.begin() +
                                leftMostUnchangedNodeOffset);
        rightMostUnchangedNodeOffset -= leftMostUnchangedNodeOffset;
    }
    leftMostUnchangedNodeOffset = 0;
    // printStatus();
}

void ConsensusGraph::removeCycles() {
    // We first iterate over all nodes on mainPath on the right
    auto edgeOnPath = mainPath.edges.begin() + rightMostUnchangedNodeOffset;
    auto edgeOnPathEnd = mainPath.edges.end();
    Node *nodeOnPath = edgeOnPath < edgeOnPathEnd ? (*edgeOnPath)->source
                                                  : edgeOnPath[-1]->sink;
    std::stack<Edge *> walkAndPruneCallStack;
    // used for converting recursion to iteration.
    // Define here to avoid overhead of repeated allocations.
    while (true) {
        // Then we iterate over all edges pointing to side nodes that have
        // other edges in
        //
        // Work with copy to avoid iterator invalidation.
        auto edgesOutCopy = nodeOnPath->edgesOut;
        for (const auto &edgeIt : edgesOutCopy)
            walkAndPrune(edgeIt, walkAndPruneCallStack);

        if (edgeOnPath == edgeOnPathEnd)
            break;
        nodeOnPath = (*edgeOnPath)->sink;
        ++edgeOnPath;
    }
    // Now we iterate over all nodes on mainPath on the left
    edgeOnPath = mainPath.edges.begin() + leftMostUnchangedNodeOffset;
    while (true) {
        // Then we iterate over all edges pointing to side nodes that have
        // other edges in
        nodeOnPath = (*edgeOnPath)->source;
        // Work with copy to avoid iterator invalidation.
        auto edgesOutCopy = nodeOnPath->edgesOut;
        for (const auto &edgeIt : edgesOutCopy)
            walkAndPrune(edgeIt, walkAndPruneCallStack);

        if (edgeOnPath == mainPath.edges.begin())
            break;
        --edgeOnPath;
    }
}

void ConsensusGraph::walkAndPrune(Edge *e, std::stack<Edge *> &callStack) {
    callStack.push(e);
    while (!callStack.empty()) {
        Edge *curr = callStack.top();
        callStack.pop();
        Node *sink = curr->sink;
        Node *source = curr->source;
        if (sink->onMainPath)
            continue;
        if (sink->edgesIn.size() > 1)
            splitPath(source, curr, &(curr->reads));
        // Now put sink->edgesOut into stack
        for (auto it = sink->edgesOut.begin(); it != sink->edgesOut.end();
             ++it) {
            callStack.push(*it);
            // OLD COMMENT:
            // Here walkAndPrune will only change sink->edgesOut by 1. deleting
            // edge2WorkOn and 2. adding new edges that we don't need to prune.
            // So copying sink->edgesOut should work.
        }
    }
}

void ConsensusGraph::splitPath(Node *newPre, Edge *e,
                               std::vector<read_t> *reads2Split) {
    /**
     * @brief The context for changing recursion to iteration.
     *
     * Every context will be visited exactly twice.
     *
     * When constructed, hasVisited is set to false, oldCur is set to
     * nullptr, and reads2Split is set to the pointer passed down by the last
     * iteration.
     *
     * After the first visit, hasVisited is set to true, reads2Split is set to
     * reads2Split created in this iteration and passed on to future iterations.
     * The new value is created by new and the destructor is responsible for
     * deleting it. And oldCur will be set to e->sink, and the destructor will
     * delete it if it becomes disconnected.
     *
     */
    class SplitPathContext {
    public:
        Node *const newPre;
        Edge *const e;
        std::vector<read_t> *reads2Split;
        bool hasVisited;
        /** we need to store this (e->sink) because e might be deleted and we
         * still need to delete e->sink **/
        Node *oldCur;

        SplitPathContext(ConsensusGraph *cG, Node *newPre, Edge *e,
                         std::vector<read_t> *reads2Split)
            : newPre(newPre), e(e), reads2Split(reads2Split), hasVisited(false),
              oldCur(nullptr), cG(cG) {}

        ~SplitPathContext() {
            delete reads2Split;
            if (oldCur && oldCur->edgesIn.empty() && oldCur->edgesOut.empty())
                cG->removeNode(oldCur);
        }

    private:
        ConsensusGraph *const cG;
    };

    std::stack<SplitPathContext> callStack;

    callStack.emplace(this, newPre, e, reads2Split);

    while (!callStack.empty()) {
        SplitPathContext &currentContext = callStack.top();

        if (currentContext.hasVisited) {
            callStack.pop();
            continue;
        }

        // The reads that need to be split going down this path
        auto readsInPath2Split = new std::vector<read_t>;
        std::set_intersection(
            currentContext.reads2Split->begin(),
            currentContext.reads2Split->end(), currentContext.e->reads.begin(),
            currentContext.e->reads.end(),
            std::inserter(*readsInPath2Split, readsInPath2Split->begin()));

        currentContext.reads2Split = readsInPath2Split;
        currentContext.hasVisited = true;

        if (readsInPath2Split->empty())
            continue;

        Node *oldCur = currentContext.e->sink;
        currentContext.oldCur = oldCur;

        // Remove the reads from the node edge
        removeReadsFromEdge(currentContext.e, *readsInPath2Split);

        // just create a new edge and return if we are arriving at a node in
        // mainPath
        if (oldCur->onMainPath) {
            createEdge(currentContext.newPre, oldCur, *readsInPath2Split);
            continue;
        }

        // Otherwise create a new node
        Node *newCur = createNode(oldCur->base);
        // Add an edge to this new node
        createEdge(currentContext.newPre, newCur, *readsInPath2Split);

        // Now put sink->edgesOut into stack
        for (auto &it : oldCur->edgesOut)
            callStack.emplace(this, newCur, it, readsInPath2Split);
    }
}

ConsensusGraph::~ConsensusGraph() {
//#ifdef  LOG
//    std::cout << "*********" << std::endl;
//    double vm_usage, resident_set;
//    _mem_usage(vm_usage, resident_set);
//    int* ptr = new int[1000000000];
//    _mem_usage(vm_usage, resident_set);
//    delete [] ptr;
//    _mem_usage(vm_usage, resident_set);
//     std::cout << "*********" << std::endl;
//#endif
//
//
//     std::cout << "Removing" << std::endl;
//    std::cout << "numEdges:"<<numEdges<<"numNodes:"<< numNodes<< std::endl;
//#ifdef LOG
//
//    std::cout<<"before"<<std::endl;
//    _mem_usage(vm_usage, resident_set);
//#endif
    std::vector<Node *> startingNodes;
    for (auto it : readsInGraph)
        startingNodes.push_back(it.second.start);
    removeConnectedNodes(startingNodes);
    if (numEdges || numNodes)
        std::cerr << "graph " << this << " " << std::to_string(numEdges)
                  << " edges " << std::to_string(numNodes) << " nodes left"
                  << std::endl;

    assert(numEdges == 0);
    assert(numNodes == 0);
//#ifdef LOG
//    std::cout<<"after"<<std::endl;
//    _mem_usage(vm_usage, resident_set);
//#endif
}

Node *ConsensusGraph::createNode(char base) {
    Node *n = new Node(base);
    numNodes++;
    return n;
}

Edge *ConsensusGraph::createEdge(Node *source, Node *sink, read_t read) {
    Edge *e = new Edge(source, sink, read);
    source->edgesOut.push_back(e);
    sink->edgesIn.push_back(e);
    numEdges++;
    return e;
}
Edge *ConsensusGraph::createEdge(Node *source, Node *sink, read_t read,read_t w) {
    Edge *e = new Edge(source, sink, read,w);
    source->edgesOut.push_back(e);
    sink->edgesIn.push_back(e);
    numEdges++;
    return e;
}
Edge *ConsensusGraph::createEdge(Node *source, Node *sink,
                                 std::vector<read_t> &reads) {
    Edge *e = new Edge(source, sink, reads);
    source->edgesOut.push_back(e);
    sink->edgesIn.push_back(e);
    numEdges++;
    return e;
}

void ConsensusGraph::removeReadsFromEdge(Edge *e,
                                         std::vector<read_t> const &reads) {
    std::vector<read_t> updatedReadsInOldEdge;
    std::set_difference(
        e->reads.begin(), e->reads.end(), reads.begin(), reads.end(),
        std::inserter(updatedReadsInOldEdge, updatedReadsInOldEdge.begin()));
    e->reads.swap(updatedReadsInOldEdge);
    e->count = e->reads.size();
    if (e->count == 0) {
        removeEdge(e);
    }
}

void ConsensusGraph::removeEdge(Edge *e,
                                bool dontRemoveFromSource /* = false */,
                                bool dontRemoveFromSink /* = false */) {
            if (!dontRemoveFromSource)
                e->source->edgesOut.erase(std::find_if(
                    e->source->edgesOut.begin(), e->source->edgesOut.end(),
                    [&](const Edge *p) {
                        return (p->sink == e->sink);
                    }));

            if (!dontRemoveFromSink)
                e->sink->edgesIn.erase(std::find(e->sink->edgesIn.begin(),
                                                 e->sink->edgesIn.end(), e));
            delete e;
            numEdges--;
        }

        void ConsensusGraph::removeNode(Node *n) {
            {
                auto edgeIt = n->edgesIn.begin();
                auto end = n->edgesIn.end();
                while (edgeIt != end)
                    removeEdge(*(edgeIt++), false, true);
                // don't remove from sink node which is this node!
                // Avoids iterator invalidation and speeds things up
            }
            {
                auto edgeIt = n->edgesOut.begin();
                auto end = n->edgesOut.end();
                while (edgeIt != end)
                    removeEdge(*(edgeIt++), true, false);
                // don't remove from source node which is this node!
                // Avoids iterator invalidation and speeds things up
            }
            delete n;
            numNodes--;
        }

        void ConsensusGraph::removeBelow(Node *n) {
            std::stack<Node *> nodesToRemove;
            nodesToRemove.push(n);
            while (!nodesToRemove.empty()) {
                Node *curNode = nodesToRemove.top();
                if (curNode->edgesOut.empty()) {
                    nodesToRemove.pop();
                    removeNode(curNode);
                } else {
                    nodesToRemove.push(curNode->edgesOut.front()->sink);
                }
            }
        }

        void ConsensusGraph::removeAbove(Node *n) {
            std::stack<Node *> nodesToRemove;
            nodesToRemove.push(n);
            while (!nodesToRemove.empty()) {
                Node *curNode = nodesToRemove.top();
                if (curNode->edgesIn.empty()) {
                    nodesToRemove.pop();
                    removeNode(curNode);
                } else {
                    nodesToRemove.push((*curNode->edgesIn.begin())->source);
                }
            }
        }

        void ConsensusGraph::removeConnectedNodes(std::vector<Node *> nodes) {
#ifdef CHECKS
            {
                // This code segment makes sure that the graph is connected and
                // traverseAndCall is working as desired
                size_t count = 0;
                for (Node *n : nodes)
                    traverseAndCall(n, false, [&count](Node *) { ++count; });
                assert(count == numNodes);
                count = 0;
                for (Node *n : nodes)
                    traverseAndCall(n, true, [&count](Node *) { ++count; });
                assert(count == numNodes);
            }
#endif

            /** Nodes that have no edges in **/
            std::stack<Node *> leafNodes;
            for (Node *n : nodes)
                traverseAndCall(n, false, [&leafNodes](Node *node) {
                    if (node->edgesOut.empty())
                        leafNodes.push(node);
                });
            while (!leafNodes.empty()) {
                Node *curNode = leafNodes.top();
                leafNodes.pop();
                for (auto edgeIt : curNode->edgesIn) {
                    Node *source = edgeIt->source;
                    if (source->edgesOut.size() == 1) {
                        leafNodes.push(source);
                    }
                }
                removeNode(curNode);
            }
        }

        void ConsensusGraph::printStatus() {
            std::cout << readsInGraph.size() << " reads in graph " << (uint64_t)this
                      << ", " << numNodes << " nodes, and " << numEdges
                      << " edges."
                      << "\n";
            std::cout << "mainPath len " << mainPath.edges.size() + 1 << " "
                      << mainPath.path.size() << " avg weight "
                      << mainPath.getAverageWeight() << " starts at "
                      << startPos << " ends at " << endPos << " len in Contig "
                      << endPos - startPos << "\n";
            double stat = mainPath.edges.size() * mainPath.getAverageWeight() /
                          (readsInGraph.size() * 9999);
            std::cout << stat / (1 - 0.17) << " " << stat / (1 - 0.17 * 1.05)
                      << " " << stat / (1 - 0.17 * 1.1) << " " << std::endl;
        }

//        void ConsensusGraph::writeMainPath(ConsensusGraphWriter &cgw) {
//            cgw.genomeFile << std::string(mainPath.path.begin(), mainPath.path.end())
//              << std::endl;
//        }

//        void ConsensusGraph::writeReads(ConsensusGraphWriter &cgw) {
//            // First we write the index of each character into the
//            // cumulativeWeight field of the nodes on mainPath
//            mainPath.edges.front()->source->cumulativeWeight = 0;
//            size_t i = 0;
//            for (auto e : mainPath.edges) {
//                e->sink->cumulativeWeight = ++i;
//            }
//            size_t totalEditDis = 0;
//
//            read_t pasId = 0;
//            for (auto it : readsInGraph) {
//                {
//                    read_t diffId = it.first - pasId;
//                    cgw.idFile.write((char*)&diffId, std::ios::binary);
//                    cgw.complementFile << (it.second.reverseComplement ? 'c' : 'n');
//                    pasId = it.first;
//                }
//                totalEditDis += writeRead(cgw.posFile, cgw.editTypeFile, cgw.editBaseFile,
//                                          it.second, it.first);
//            }
//            cgw.complementFile << '\n';
//#ifdef LOG
//            std::cout << "AvgEditDis "
//                      << totalEditDis / (double)readsInGraph.size()
//                      << std::endl;
//            printStatus();
//#endif
//        }
        void ConsensusGraph::writeReads(ConsensusGraphWriter &cgw,const read_t ref_id,std::shared_ptr<std::vector<std::shared_ptr<my_read::Read>>> &readvecptr){
            // First we write the index of each character into the
            // cumulativeWeight field of the nodes on mainPath
            mainPath.edges.front()->source->cumulativeWeight = 0;
            size_t i = 0;
            for (auto e : mainPath.edges) {
                e->sink->cumulativeWeight = ++i;
            }
            size_t totalEditDis = 0;

            read_t pasId = 0;
            for (auto it : readsInGraph) {
                {
//                    assert(it.first>pasId);
                    read_t diffId = it.first - pasId;
                    cgw.idFile.write((char*)&diffId, std::ios::binary);
                    cgw.complementFile << (it.second.reverseComplement ? 'c' : 'n');
                    DirectoryUtils::write_var_uint32(ref_id, cgw.refidFile);
                    pasId = it.first;
                }
                totalEditDis += writeRead(cgw.posFile, cgw.editTypeFile, cgw.editBaseFile,cgw.refidFile,readvecptr,
                                          it.second, it.first);
            }
//            cgw.complementFile << '\n';
#ifdef LOG
//            std::cout << "AvgEditDis "
//                      << totalEditDis / (double)readsInGraph.size()
//                      << std::endl;
//            printStatus();
#endif
        }

        void ConsensusGraph::writeReadLone(ConsensusGraphWriter &cgw) {
            cgw.loneFile << std::string(mainPath.path.begin(), mainPath.path.end()) << std::endl;
        }
        void ConsensusGraph::writeReadlone(ConsensusGraphWriter &cgw,std::shared_ptr<my_read::Read> &readptr){
//            size_t readid=readptr->id;
//            cgw.loneFile<<std::to_string(readid)<<std::endl;
//            cgw.loneFile<<std::string(mainPath.path.begin(), mainPath.path.end())<<std::endl;
              cgw.readsLone.insert(std::make_pair(readptr->id,readptr));
        }

        void ConsensusGraph::writeIdsLone(ConsensusGraphWriter &cgw, std::vector<read_t> &loneReads) {
            read_t pasId = 0; 
            for (auto it : loneReads) {
                read_t diffId = it - pasId;
                cgw.idFile.write((char*)&diffId, std::ios::binary);
                pasId = it;
            }
        }

        read_t ConsensusGraph::getNumReads() { return readsInGraph.size(); }

        size_t ConsensusGraph::getNumEdges(){return numEdges;} 

        size_t ConsensusGraph::
        read2EditScript(ConsensusGraph::Read &r,
                                               read_t id,
                                               std::vector<Edit> &editScript,
                                               uint32_t &pos) {
            editScript.clear();
            // First we store the initial position
            Node *curNode = r.start;
            assert(curNode);
            bool intersectWithMainPath = true;
            while (!curNode->onMainPath) {
                curNode = curNode->getNextNodeInRead(id);
                if (!curNode) {
                    intersectWithMainPath = false;
                    break;
                }
            }

            // When the read has no intersection with mainPath, we just store it
            // as a bunch of inserts
            if (!intersectWithMainPath) {
                pos = 0;
                size_t editDis = 0;
                Node *curNode = r.start;
                do {
                    editScript.push_back(Edit(INSERT, curNode->base));
                    ++editDis;
                } while ((curNode = curNode->getNextNodeInRead(id)));

                return editDis;
            }

                pos = curNode->cumulativeWeight;

                size_t editDis = 0;
                size_t posInMainPath = curNode->cumulativeWeight;
                curNode = r.start;
                size_t unchangedCount = 0;
                auto dealWithUnchanged = [&unchangedCount, &editScript]() {
                    if (unchangedCount > 0) {
                        editScript.push_back(Edit(SAME, unchangedCount));
                    unchangedCount = 0;
                }
            };
            do {
                if (curNode->onMainPath) {
                    size_t curPos = curNode->cumulativeWeight;
                    if (curPos > posInMainPath)
                        dealWithUnchanged();
                    for (; posInMainPath < curPos; posInMainPath++) {
                        editScript.push_back(Edit(DELETE, '-'));
                        editDis++;
                    }
                    unchangedCount++;
                    posInMainPath++;
                } else {
                    dealWithUnchanged();
                    // Else we have an insertion
                    editScript.push_back(Edit(INSERT, curNode->base));
                    editDis++;
                }
            } while ((curNode = curNode->getNextNodeInRead(id)));

            dealWithUnchanged();

            return editDis;
        }

        size_t ConsensusGraph::writeRead(std::ofstream &posFile,
                                         std::ofstream &editTypeFile,
                                         std::ofstream &editBaseFile,std::ofstream &refidFile,std::shared_ptr<std::vector<std::shared_ptr<my_read::Read>>> &readvecptr, Read &r,
                                         read_t id) {

            uint32_t offset;
            std::vector<Edit> editScript;
//#ifdef LOG
//            if(r.start== nullptr)
//            {
//                std::cout<<id<<" "<<r.reverseComplement<<" "<<r.pos<<std::endl;
//            }
//#endif
            size_t editDis = read2EditScript(r, id, editScript, offset);
            //posFile.write((char*)&offset,sizeof(uint32_t));
            DirectoryUtils::write_var_uint32(offset, posFile);
            static int readnum=0;
            std::vector<Edit> newEditScript;
            editDis = Edit::optimizeEditScript(editScript, newEditScript);
            //find the number of consecutive insertions at beginning
            uint32_t numInsStart = 0;
            uint32_t numInsEnd = 0;
            std::string subread;
            for (size_t i = 0; i != newEditScript.size(); i++){
            	if(newEditScript[i].editType != INSERT)
            		break;
            	numInsStart++;
                subread.push_back(newEditScript[i].editInfo.ins);
            	//write the inserted bases to .base
//            	editBaseFile << newEditScript[i].editInfo.ins;
            }
            //前段不为空，转发到下一层
            if(numInsStart!=0)
            {
                if(subread.size()>100)readnum++;
                readvecptr->emplace_back(get_newsubread(subread));
                DirectoryUtils::write_var_uint32(readvecptr->back()->id, refidFile);
            }
            else
            {
                DirectoryUtils::write_var_uint32(0, refidFile);//为空写入0,0代表空序列
            }
            //check if it is the case with all insertion
            subread.clear();
            if(numInsStart != newEditScript.size()){
            	//find the number of consecutive insertions at end
	            for (int64_t i = newEditScript.size()-1; i >=0; i--){
	            	if(newEditScript[i].editType != INSERT)
	            		break;
	            	numInsEnd++;
                    subread.push_back(newEditScript[i].editInfo.ins);
	            }     
	        }
           std::reverse(subread.begin(),subread.end());
            if(numInsEnd!=0)
            {
                if(subread.size()>100)readnum++;
                readvecptr->emplace_back(get_newsubread(subread));
                DirectoryUtils::write_var_uint32(readvecptr->back()->id, refidFile);
            }
            else
            {
                DirectoryUtils::write_var_uint32(0, refidFile);//为空写入0,0代表空序列
            }
#ifdef LOG
            ;;;

            std::cout<<"new readnum>100:"<<readnum<<std::endl;
#endif



            //Write numInsStart to .pos
//            DirectoryUtils::write_var_uint32(numInsStart, posFile);

            uint32_t unchangedCount = 0;
            //Go through editScript from numInsStart to lenEditScript-numInsEnd
            for (size_t i = numInsStart; i < newEditScript.size()-numInsEnd; i++) {
                switch (newEditScript[i].editType) {
                case SAME: {
                    unchangedCount += newEditScript[i].editInfo.num;
                    break;
                }
                case INSERT: {
                    // posFile.write((char*)&unchangedCount, sizeof(uint32_t));
            		DirectoryUtils::write_var_uint32(unchangedCount, posFile);
                    unchangedCount = 0;
                    editTypeFile << 'i';
                    editBaseFile << newEditScript[i].editInfo.ins;
                    break;
                }
                case DELETE: {
                    // posFile.write((char*)&unchangedCount, sizeof(uint32_t));
            		DirectoryUtils::write_var_uint32(unchangedCount, posFile);
                    unchangedCount = 0;
                    editTypeFile << 'd';
                    break;
                }
                case SUBSTITUTION: {
                    // posFile.write((char*)&unchangedCount, sizeof(uint32_t));
            		DirectoryUtils::write_var_uint32(unchangedCount, posFile);
                    unchangedCount = 0;
                    editTypeFile << 's';
                    editBaseFile << newEditScript[i].editInfo.sub;
                    break;
                }
                }
            }

            // posFile.write((char*)&unchangedCount, sizeof(uint32_t));
            DirectoryUtils::write_var_uint32(unchangedCount, posFile);
            //Write numInsEnd to .pos file, and write the inserted bases to .base
            DirectoryUtils::write_var_uint32(numInsEnd, posFile);
            for (size_t i = newEditScript.size()-numInsEnd; i != newEditScript.size(); i++){
            	//write the inserted bases to .base
            	editBaseFile << newEditScript[i].editInfo.ins;
            }
            
            editTypeFile << '\n';
            return editDis;
        }
        
        ConsensusGraph::ConsensusGraph() {}

        ConsensusGraph::Read::Read(long pos, Node *start, size_t len,
                                   bool reverseComplement)
            : pos(pos), start(start), len(len),
              reverseComplement(reverseComplement) {}

        bool ConsensusGraph::checkNoCycle() {
            size_t countNodes = 0;
            size_t countEdges = 0;
            // In this function cumulativeWeight == 1 means is parent of current
            // Node, and cumulativeWeight == 0 means otherwise
            std::vector<Node *> sourceNodes;
            for (auto it : readsInGraph)
                traverseAndCall(
                    it.second.start, false,
                    [&sourceNodes, &countNodes, &countEdges](Node *node) {
                        ++countNodes;
                        countEdges += node->edgesIn.size();
                        node->cumulativeWeight = 0;
                        if (node->edgesIn.empty())
                            sourceNodes.push_back(node);
                    });
            assert(countNodes == numNodes);
            assert(countEdges == numEdges);
            assert(!sourceNodes.empty());
            // Here all the .hasReached has been set to true
            bool status = true;
            for (Node *node : sourceNodes) {
                std::stack<Node *> nodes2Visit;
                nodes2Visit.push(node);
                node->hasReached = !status;
                while (!nodes2Visit.empty()) {
                    Node *currentNode = nodes2Visit.top();
                    // cumulativeWeight == 1 means parent of current Node
                    // (including itself)
                    currentNode->cumulativeWeight = 1;
                    bool hasUnvisitedChild = false;
                    for (auto it : currentNode->edgesOut) {
                        Node *n = it->sink;
                        // A back edge
                        if (n->cumulativeWeight == 1)
                            return false;
                        // == status means has not visited
                        if (n->hasReached == status) {
                            nodes2Visit.push(n);
                            n->hasReached = !status;
                            hasUnvisitedChild = true;
                            continue;
                        }
                    }
                    if (hasUnvisitedChild)
                        continue;
                    // All children has been visited
                    currentNode->cumulativeWeight = 0;
                    nodes2Visit.pop();
                }
            }
            return true;
        }

        template <typename Functor>
        void ConsensusGraph::traverseAndCall(Node *n, bool status, Functor f) {
            std::deque<Node *> unfinishedNodes;
            unfinishedNodes.push_back(n);
            while (!unfinishedNodes.empty()) {
                Node *currentNode = unfinishedNodes.front();
                // We don't do anything if it has already been set to !status
                if (currentNode->hasReached != status) {
                    unfinishedNodes.pop_front();
                    continue;
                }

                currentNode->hasReached = !status;
                unfinishedNodes.pop_front();
                f(currentNode);

                for (auto edgeIt : currentNode->edgesOut) {
                    Node *n = edgeIt->sink;
                    if (n->hasReached == status) {
                        unfinishedNodes.push_front(n);
                    }
                }

                for (auto edgeIt : currentNode->edgesIn) {
                    if (edgeIt->source->hasReached == status)
                        unfinishedNodes.push_back(edgeIt->source);
                }
            }
        }
