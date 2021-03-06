#include "Consensus.h"
#include "DirectoryUtils.h"
#include "bsc_helper.h"
#include "minimap.h"
#include <arpa/inet.h>
#include <boost/filesystem.hpp>
#include <algorithm>
#include <csignal>
#include <fstream>
#include <iostream>
#include <iterator>
#include <memory>
#include <mutex>
#include <set>
#include <ctime>
#include <chrono>
#ifdef USE_MALLOC_TRIM
    #include <malloc.h>
#endif
void _mem_usage(double& vm_usage, double& resident_set) {
    // from https://www.tutorialspoint.com/how-to-get-memory-usage-at-runtime-using-cplusplus
    using namespace std;
    vm_usage = 0.0;
    resident_set = 0.0;
    ifstream stat_stream("/proc/self/stat",ios_base::in); //get info from proc
    // directory
    //create some variables to get info
    string pid, comm, state, ppid, pgrp, session, tty_nr;
    string tpgid, flags, minflt, cminflt, majflt, cmajflt;
    string utime, stime, cutime, cstime, priority, nice;
    string O, itrealvalue, starttime;
    unsigned long vsize;
    long rss;
    stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
                >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
                >> utime >> stime >> cutime >> cstime >> priority >> nice
                >> O >> itrealvalue >> starttime >> vsize >> rss; // don't care
    // about the rest
    stat_stream.close();
    long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // for x86-64 is configured
    // to use 2MB pages
    vm_usage = vsize / 1024.0;
    resident_set = rss * page_size_kb;
    cout << "resident_set: " << resident_set << " KB\n";
}
void Consensus::generateAndWriteConsensus(size_t loopindex,read_t &mergeCnt) {

    initialize(loopindex);//只在第一轮判断是否是重复的
    std::vector<std::vector<read_t>> numReadsInContig(numThr);
    std::vector<std::vector<read_t>> loneReads(numThr);
    std::vector<std::vector<std::string>>contig(numThr);
    //用来保存每个线程的的合并序列
    std::vector<std::shared_ptr<std::vector<std::shared_ptr<my_read::Read>>>>reads(numThr);
    //初始化必须先给智能指针分配内存
    for(int i=0;i<numThr;i++)
    {
//        reads[i]=std::shared_ptr<std::vector<std::shared_ptr<my_read::Read>>>(new std::vector<std::shared_ptr<my_read::Read>>());
          reads[i]=std::make_shared<std::vector<std::shared_ptr<my_read::Read>>>();
    }
    ConsensusGraph *cG = nullptr;
    std::vector<CountStats> count_stats(numThr);
    std::vector<size_t>mergereads_cnt(numThr,0);
//#ifdef LOG
//    std::vector<size_t>mergereads_cnt(numThr,0);
//#endif
#pragma omp parallel private(cG)
    {
        auto tid = omp_get_thread_num();   
        std::ofstream logfile;
#ifdef LOG

        logfile.open("logfile"+std::to_string(tid), std::ofstream::out);
#endif
        std::string filePrefix = tempDir + tempFileName+".loop."+ std::to_string(loopindex) + ".tid." + std::to_string(tid);

        ConsensusGraphWriter cgw(filePrefix);
        read_t firstUnaddedRead = 0;
//#ifdef LOG
//        std::cout<<filePrefix<<std::endl;
//#endif
        // guarantee that all reads < firstUnaddedRead have been picked
        while ((cG = createGraph(firstUnaddedRead))) {
//#ifdef LOG
//            std::cout<<cG->firstReadId<<" ************"<<rD->getRead(cG->firstReadId)->cnt<<std::endl;
//            if( cG->firstReadId%1000==0)
//            {
//                std::cout<<"loopindex:"<<loopindex<<" "<<cG->mainPath.path<<" "<<cG->mainPath.path.size()<<std::endl;
//            }
//#endif
            if (isRepetitive[cG->firstReadId]||cG->mainPath.path.size()<100)//自重复的直接压缩或者小于200
            {
                std::shared_ptr<my_read::Read> tmp_ptr;
                rD->getRead(cG->firstReadId,tmp_ptr);
                cG->writeReadlone(cgw,tmp_ptr);
                continue;
            }
            if (addRelatedReads(cG))//如果合并成功
            {
                //形成新的序列放入线程的vector

                reads[tid]->emplace_back(std::move(cG->get_mainpathread(rD->recover_cnt)));
                //给图中的两条序列，构造编辑脚本写入文件夹
//#ifdef LOG
//                if(cG->firstReadId==182)
//                {
//                   std::cout<<"done here"<<std::endl;
//                }
//#endif
                cG->writeReads(cgw,reads[tid]->back()->id,reads[tid]);
                mergereads_cnt[tid]+=2;
//#ifdef LOG
//
//                mergereads_cnt[tid]+=2;
//#endif
            }
            else//合并不成功
             {
                /*判断剩余的挽回次数
                 * 不为0直接转发，下一层
                 *为0直接写入孤立文件流
                */
//#ifdef LOG
//                std::cout<<cG->firstReadId<<" not find simerl read"<<std::endl;
//#endif
                std::shared_ptr<my_read::Read> tmp_ptr;
                rD->getRead(cG->firstReadId,tmp_ptr);
                 tmp_ptr->cnt--;
                if(tmp_ptr->cnt>0)
                {
//#ifdef  LOG
//                    std::cout<<"done here ***"<<(tmp_ptr== nullptr?0:1)<<" "<<cG->firstReadId<<std::endl;
//                    std::cout<<tid<<" "<<numThr<<std::endl;
//#endif
                    assert(tid<numThr);
                    assert(cG->firstReadId<rD->getNumReads());
                    reads[tid]->emplace_back(rD->getRead(cG->firstReadId));
                }
                else
                {
                    cG->writeReadlone(cgw,tmp_ptr);
                }
            }

            // if the graph is large, run malloc_trim so that memory released to system
            // without this the memory keeps
            // increasing to very high
            bool run_malloc_trim = false;
            if (cG->getNumEdges() > 1000000)
                run_malloc_trim = true;


            delete cG;



//#ifdef LOG
//            std::cout<<"dead here in loop 1 "<<std::endl;
//#endif
#ifdef USE_MALLOC_TRIM
            if (run_malloc_trim)
                malloc_trim(0);
#endif


#ifdef LOG
            {
            auto end = std::chrono::system_clock::now();
            std::time_t end_time = std::chrono::system_clock::to_time_t(end);
            logfile <<"Time: "<<std::ctime(&end_time);
            }
#endif
        }
        cgw.writeLoneReadtoFile();
        cgw.closefilestream();
        // finally, write ids of lone reads to end of id file

//        cG->writeIdsLone(cgw, loneReads[tid]);
//#ifdef LOG
//        std::cout<<"dead here in loop 3 "<<std::endl;
//#endif
//#ifdef LOG
//        {
//        auto end = std::chrono::system_clock::now();
//        std::time_t end_time = std::chrono::system_clock::to_time_t(end);
//        logfile<<"Thread "<<tid<<" ends at time: "<<std::ctime(&end_time);
//        }
//        logfile.close();
//#endif
    }

    // pragma omp parallel
    //把本轮的合并成功，以及查询次数不为0的的序列，合并成一个vector
//#ifdef LOG
//    std::cout<<"开始合并 vec"<<std::endl;
//    for(int i=0;i<numThr;i++)
//    {
//        std::cout<<reads[i]->size()<<" ";
//    }std::cout<<std::endl;
//    size_t merge_cnt=0;
//    for(int i=0;i<mergereads_cnt.size();i++)merge_cnt+=mergereads_cnt[i];
//    std::cout<<"loopindex:"<<loopindex<<"merge read cnt"<<merge_cnt<<" "<<mergereads_cnt[0]<<std::endl;
//#endif
    size_t merge_cnt=0;


    for(int i=0;i<mergereads_cnt.size();i++)merge_cnt+=mergereads_cnt[i];
    std::cout<<"loopindex:"<<loopindex<<"merge read cnt"<<merge_cnt<<" "<<mergereads_cnt[0]<<std::endl;
    mergeCnt+=merge_cnt;
    for(int i=1;i<numThr;i++)
    {

        reads[0]=CombineVectors(std::move(reads[0]),std::move(reads[i]));
//#ifdef LOG


//        std::cout<<"i"<<" "<<i<<" "<<reads[0]->size()<<std::endl;
//#endif
    }
//#ifdef LOG
//        std::cout<<"完成合并vec "<<std::endl;
//#endif
    //然后转发给rD类
    rD->setReads(reads[0]);
    //把本轮得到的编辑脚本提交给线程池，让他复责通用压缩任务。

    // now perform last step, combining files from threads and writing metadata
//    finishWriteConsensus(numReadsInContig);





//    // compute total count stats by adding for all threads
//    CountStats summary;
//    for (auto &c : count_stats)
//        summary = summary + c;
//
//    std::cout << "\n";
////    std::cout << "#LoneReads = " << totalNumLoneReads << "\n";
//    std::cout << "MinHash passed " << summary.countMinHash << std::dec << " reads\n";
//    std::cout << "MinHash passed & not already in graph " << summary.countMinHashNotInGraph << std::dec << " reads\n";
//    std::cout << "Merge Sort passed " << summary.countMergeSort << std::dec << " reads\n";
//    std::cout << "Aligner passed " << summary.countAligner << std::dec << " reads\n";
}

void Consensus::writeLoneReadtofile(size_t loopindex) {
    std::string filePrefix = tempDir + tempFileName+".loop."+ std::to_string(loopindex) + ".tid." + std::to_string(0);
    ConsensusGraphWriter cgw(filePrefix);
    for(size_t i=0;i<rD->getNumReads();i++)
    {
        auto tmpptr=rD->getRead(i);
        DirectoryUtils::write_var_uint32(tmpptr->id,cgw.loneidFile);
        std::string read;
        tmpptr->read->to_string(read);
        cgw.loneFile<<read<<std::endl;
    }
}

bool  Consensus::addRelatedReads(ConsensusGraph *cG)
{
    const std::string originalString=cG->mainPath.path;
    std::string reverseComplementString;
    auto stringBegin = cG->mainPath.path.begin() ;
    auto stringEnd =cG->mainPath.path.end();
    ReadData::toReverseComplement(
            stringBegin, stringEnd,
            std::inserter(reverseComplementString, reverseComplementString.end()));

    std::vector<match_info> match_infos(2*OrderHashReadFilter::que_cnt);//用来储存和候选序列的比对结果。
    size_t sum_index=0;
    bool all[] = {false, true};

    for (bool reverseComplement : all) {
        std::vector<read_t> results;

        rF->getFilteredReads(reverseComplement ? reverseComplementString
                                               : originalString,
                             results, OrderHashReadFilter::que_cnt);

//#ifdef LOG
////        if(reverseComplement==0)results.push_back(1);
//        std::cout<<reverseComplement<<" res.size:"<<results.size()<<std::endl;
//#endif
        if(results.size()==0)continue;
        // Try to add them one by one
        for (const auto r: results) {

            // check if we exceed the edge limit in the graph
//#ifdef LOG
////            std::cout<<r<<std::endl;
////            std::cout<<cG->getNumEdges()<<" edgenum:"<<edge_threshold<<std::endl;
//#endif
            if (cG->getNumEdges() >= edge_threshold) {
                return false;
            }
            //check if the read is repetitive
            if (isRepetitive[r])
                continue;
            //check if it is already in graph
            if (inGraph[r])
                continue;

            std::string readStr, readStr1;
            rD->getRead(r, readStr1);
            if (readStr1.size() < 32) { // for len below 32, minhash is meaningless
                continue;
            }
//#ifdef LOG
//            logfile<<"Contig: " << contigId << ", Read passed MinHash "<<r<<", read length: " << readStr1.length()<< "\n";
//#endif
            if (reverseComplement)
                ReadData::toReverseComplement(
                        readStr1.begin(), readStr1.end(),
                        std::inserter(readStr, readStr.end()));
            else
                readStr = readStr1;
//            std::vector<Edit> editScript;
//            ssize_t beginOffset, endOffset;
//            ssize_t pos;

            match_infos[sum_index].vid=r;

            rD->getindex(r,match_infos[sum_index].id);
            match_infos[sum_index].reverseComplement = reverseComplement;

            bool alignStatus = cG->alignRead(readStr, match_infos[sum_index].editScript, match_infos[sum_index].pos,
                                             match_infos[sum_index].beginOffset,
                                             match_infos[sum_index].endOffset, match_infos[sum_index].match_length,\
                                             match_infos[sum_index].editdis,match_infos[sum_index].ref_st,match_infos[sum_index].ref_ed,\
                                             match_infos[sum_index].que_st,match_infos[sum_index].que_ed,
                                             m_k,
                                             m_w, max_chain_iter);
//check 条件编译用来判断编辑脚本能否还原出原始序列
#ifdef CHECKS
            {
                // Check editScript applied to originalString is readStr
                std::string origString = std::string(
                        cG->mainPath.path.begin() +
                        (match_infos[sum_index].beginOffset > 0 ? match_infos[sum_index].beginOffset : 0),
                        cG->mainPath.path.end() + (match_infos[sum_index].endOffset > 0 ? 0 : match_infos[sum_index].endOffset));
                std::string targetString = readStr.substr(
                        match_infos[sum_index].beginOffset > 0 ? 0 : -match_infos[sum_index].beginOffset,
                        readStr.length() - (match_infos[sum_index].beginOffset > 0 ? 0 : -match_infos[sum_index].beginOffset) -
                        (match_infos[sum_index].endOffset> 0 ? match_infos[sum_index].endOffset : 0));
                std::string resultAfterEdit;
                Edits::applyEdits(
                        origString.begin(), match_infos[sum_index].editScript,
                        std::inserter(resultAfterEdit, resultAfterEdit.end()));
                if (resultAfterEdit.compare(targetString)) {
                    std::cout << "beginOffset " <<match_infos[sum_index].beginOffset << " endOffset "
                              << match_infos[sum_index].endOffset << std::endl;
                    std::cout << "origString\n" << origString << std::endl;
                    std::cout << "targetString\n" << targetString << std::endl;
                    std::cout << "resultAfterEdit\n"
                              << resultAfterEdit << std::endl;
                    std::cout << "editScript\n";
                    for (auto e : match_infos[sum_index].editScript)
                        std::cout << e;
                    std::cout << std::endl;
                    std::cout << "mainPath\n"
                              << std::string(cG->mainPath.path.begin(),
                                             cG->mainPath.path.end())
                              << std::endl;
                    std::cout << "pos " << match_infos[sum_index].pos << " endPos " << cG->endPos
                              << " offsetGuess "
                              << cG->mainPath.path.size() - cG->endPos + match_infos[sum_index].pos
                              << std::endl;
                }
                assert(!resultAfterEdit.compare(targetString));
            }
#endif

            if(alignStatus)
                sum_index++;

//#ifdef CHECKS
//            {
//                // Check editScript applied to originalString is readStr
//                std::string origString = std::string(
//                    cG->mainPath.path.begin() +
//                        (beginOffset > 0 ? beginOffset : 0),
//                    cG->mainPath.path.end() + (endOffset > 0 ? 0 : endOffset));
//                std::string targetString = readStr.substr(
//                    beginOffset > 0 ? 0 : -beginOffset,
//                    readStr.length() - (beginOffset > 0 ? 0 : -beginOffset) -
//                        (endOffset > 0 ? endOffset : 0));
//                std::string resultAfterEdit;
//                Edits::applyEdits(
//                    origString.begin(), editScript,
//                    std::inserter(resultAfterEdit, resultAfterEdit.end()));
//                if (resultAfterEdit.compare(targetString)) {
//                    std::cout << "beginOffset " << beginOffset << " endOffset "
//                              << endOffset << std::endl;
//                    std::cout << "origString\n" << origString << std::endl;
//                    std::cout << "targetString\n" << targetString << std::endl;
//                    std::cout << "resultAfterEdit\n"
//                              << resultAfterEdit << std::endl;
//                    std::cout << "editScript\n";
//                    for (auto e : editScript)
//                        std::cout << e;
//                    std::cout << std::endl;
//                    std::cout << "mainPath\n"
//                              << std::string(cG->mainPath.path.begin(),
//                                             cG->mainPath.path.end())
//                              << std::endl;
//                    std::cout << "pos " << pos << " endPos " << cG->endPos
//                              << " offsetGuess "
//                              << cG->mainPath.path.size() - cG->endPos + pos
//                              << std::endl;
//                }
//                assert(!resultAfterEdit.compare(targetString));
//            }
//#endif


        }
    }
//#ifdef LOG
//    std::cout<<"sumindex index :"<<sum_index<<std::endl;
//#endif
        //如果sum_index==0,说明没有找到相似序列
        if(sum_index==0)return false;
        //排序按匹配区间，从前往后取匹配区间最长的序列，因为有可能它找到的最相似的序列，有可能已经被其他线程占用了
        //最长被占用了，就取下一条，依此类推

        //或许这里还需要加上一一个配区域长度与原有序列的占比的限制条件
//        for(int i=0;i<sum_index;i++)
//        {
//
//        }


        std::sort(match_infos.begin(),match_infos.begin()+sum_index,[](const match_info &a,const match_info &b){return  a.match_length>b.match_length;});
//        std::sort(match_infos.begin(),match_infos.begin()+sum_index,[](const match_info &a,const match_info &b){return  a.editdis<b.editdis;});



    int choose_index=-1;
        for(int i=0;i<sum_index; i++)
        {
            //给通过比对区域最长的序列加锁,修改状态为已经被合并
            auto r=match_infos[i].vid;
            if (!readStatusLock[r%numLocks].try_lock()) {
                // we only try_lock here since missing a read
                // is not a major issue and lock contention should be a rare event anyway.
                // Note that if some other thread has locked the read, they are guaranteed
                // to pick it up.
                continue;
            } else {
                // check again that read is available (variables flushed after lock is set)
                if (inGraph[r]) {
                    // read already taken, continue with next read
                    readStatusLock[r%numLocks].unlock();
                    continue;
                }

                // read added to graph
//#ifdef LOG
//                logfile<< "Contig: " << contigId << ", Read passed aligner "<<r<<"\n";
//#endif
                inGraph[r] = true;
                choose_index=i;
                readStatusLock[r%numLocks].unlock();
                break;
            }
        }
        //choose_index=-1,找到的相似序列都背其他线程抢走了
        if(choose_index==-1)
        {
            return  false;
        }
//#ifdef LOG
//        if(cG->firstReadId==75)
//        {
//            std::cout<<choose_index<<std::endl;
//        }
//#endif
        //如果图还没有添加序列，依靠第一条序列初始化
        if (cG->getNumReads() == 0) {
            std::string mainPathString(cG->mainPath.path);
            cG->mainPath.path.clear();
            std::shared_ptr<my_read::Read> sharedPtr;
            rD->getRead(cG->firstReadId,sharedPtr);
            cG->initialize(sharedPtr, 0);
            cG->calculateMainPathGreedy();
        }
        //更新图，把第choose_index条比对得到的编辑脚本添加到图中
//        std::string readStr;
//        rD->getRead(match_infos[choose_index].vid,readStr);

        auto tmp_ptr=rD->getRead(match_infos[choose_index].vid);
        cG->updateGraph(tmp_ptr, match_infos[choose_index].editScript, match_infos[choose_index].beginOffset, match_infos[choose_index].endOffset, match_infos[choose_index].id, match_infos[choose_index].pos,
                        match_infos[choose_index].reverseComplement,match_infos[choose_index].ref_st,match_infos[choose_index].ref_ed,match_infos[choose_index].que_st,match_infos[choose_index].que_ed);

        //#ifdef LOG
//        std::cout<<"updateGraph SECCEES"<<std::endl;
//        assert(cG->checkNoCycle());
//#endif
#ifdef CHECKS
            assert(checkRead(cG, r));
            assert(cG->checkNoCycle());
#endif
//            cG->calculateMainPathGreedy();
              cG->calculateMainPath();
#ifdef CHECKS
            // std::cout << "Added read " << r << " first unadded read "
            //           << firstUnaddedRead << std::endl;
            assert(checkRead(cG, r));
            assert(cG->checkNoCycle());
#endif
    return true;
}

bool Consensus::checkRead(ConsensusGraph *cG, read_t read) {
    std::string result;
    std::string readStr;
    bool temp = cG->getRead(read, std::inserter(result, result.end()));
    assert(temp);
    if (!temp)
        return false;
    rD->getRead(read, readStr);
    if (!cG->readsInGraph.at(read).reverseComplement) {
        if (result != readStr) {
            std::cout << "readInGraph:\n" << result << "\n";
            std::cout << "actualRead:\n" << readStr << "\n";
        }
        return result == readStr;
    }
    std::string reverseComplement;
    ReadData::toReverseComplement(
        result.begin(), result.end(),
        std::inserter(reverseComplement, reverseComplement.end()));
    if (reverseComplement != readStr) {
        std::cout << "readInGraphRerseComplement:\n"
                  << reverseComplement << "\n";
        std::cout << "readInGraph:\n" << result << "\n";
        std::cout << "actualRead:\n" << readStr << "\n";
    }
    return reverseComplement == readStr;
}

void Consensus::finishWriteConsensus(const std::vector<std::vector<read_t>>& numReadsInContig) {
    size_t size = 0;
    for (int i = 0; i < numThr; i++)
        size += numReadsInContig[i].size();
    std::ofstream metaData;
    metaData.open(tempDir + "metaData");
    metaData << "numReads=" << rD->getNumReads() << '\n';
    metaData << "numContigs=" << size << '\n';
    std::cout << "numContigs = " << size << '\n';
    metaData << "numThr=" << numThr << '\n';
    metaData << "numReadsInContig=";
    for (int i = 0; i < numThr; ++i)
        for (size_t j = 0; j < numReadsInContig[i].size(); ++j)
            metaData << numReadsInContig[i][j] << ":";
    metaData << '\n';
    metaData.close();
}

ConsensusGraph *Consensus::createGraph(read_t &firstUnaddedRead) {
    // Note: firstUnaddedRead is local to the thread, and is a lower bound
    // on the actual first unadded read.
    read_t read = firstUnaddedRead;
    if (!getRead(read))
        return nullptr;
    ConsensusGraph *cG = new ConsensusGraph;
    // we don't actually fully initialize the graph
    // since that would be wasted effort if this turns out to be a lone graph
    rD->getRead(read, cG->mainPath.path);
    cG->startPos = 0;
    cG->endPos = cG->mainPath.path.size();
    cG->firstReadId = read;
    firstUnaddedRead = read + 1;
    return cG;
}

bool Consensus::checkRepetitive(read_t readID){
//    std::string readStr;
//    rD->getRead(readID, readStr);
//    //check the hamming distance with different offsets
//    //I choose 6 here; it is tunable
//    size_t readLen = readStr.length();
//    for (size_t i = 1; i <= 6; i++){
//        size_t countSameBase = 0;
//        for(size_t j = 0; j < readLen; j++){
//            //check how many same bases the two string share
//            if(readStr[j] == readStr[(j+i)%readLen])
//                countSameBase++;
//        }
//        //the ratio is tunable
//        if(countSameBase > 0.7 * (double)readLen)
//            return true;
//    }
//    return false;
    std::string s;
    rD->getRead(readID, s);
    ssize_t maxI = s.length() - k + 1;
    if (maxI <= 0)
        return false;
    std::unordered_map<kMer_t, unsigned> occurrences;
    kMer_t currentKMer = OrderHashReadFilter::kMerToInt(s.substr(0, k));
    occurrences[currentKMer]=1;
    const unsigned long long mask = (1ull << (2 * k)) - 1;

    for (size_t i = 1; i < (size_t)maxI; ++i) {
        currentKMer =
                ((currentKMer << 2) | OrderHashReadFilter::baseToInt(s[i + k - 1])) &
                mask;
          occurrences[currentKMer]++;
    }
    return  occurrences.size()<=maxI/2;
}

void Consensus::initialize(size_t loopindex) {
    numReads = rD->getNumReads();
    inGraph.resize(numReads, false);
    readStatusLock.resize(numLocks);
    //check if any reads are repetitive
    isRepetitive.resize(numReads, false); 
     size_t countRepeats = 0;
    if(loopindex>0)return;
#pragma omp parallel for
    for (read_t i = 0; i < numReads; i++){
        isRepetitive[i] = checkRepetitive(i);
         if(isRepetitive[i])
             countRepeats++;
        // std::cout<<"The read number "<<i<<" is "<<isRepetitive[i]<<std::endl;
    }
     std::cout<<"The number of repetitive reads is "<<countRepeats<<std::endl;
}

bool Consensus::getRead(read_t &read) {
    if (read >= numReads)
        return false;
    while (read < numReads) {
        if (!inGraph[read]) {
            if (!readStatusLock[read%numLocks].try_lock()) {
                // couldn't obtain lock, that's fine, some other thread
                // will take the read
                read++;
            } else {
                // check once again inside locked region (since variables are flushed)
                if (!inGraph[read]) {
                    inGraph[read] = true;
                    readStatusLock[read%numLocks].unlock();
                    return true;
                }
                readStatusLock[read%numLocks].unlock();
                read++;
            }
        } else {
            read++;
        }
    }
    return false;
}

Consensus::Consensus() {}

