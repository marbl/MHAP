#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <math.h>

#include "bio++.H"
#include "meryl.H"
#include "merStream.H"

#include "swalign.H"

#include <unordered_map>
#include <vector>
#include <set>

using namespace std;

#define MAX_HASH 0xffffffffffffffff

int
main(int argc, char **argv) {
  // parse args
  merylArgs   *args = new merylArgs(argc, argv);

  // set up streams for reading kmers and sequences
  seqStream *seqstr = new seqStream(args->inputFile);
  args->numBasesActual = 0;
  for (u32bit i=0; i<seqstr->numberOfSequences(); i++)
    args->numBasesActual += seqstr->lengthOf(i);

  merStream *M = new merStream(new kMerBuilder(args->merSize, args->merComp), seqstr, true, true);
  args->numMersActual  = M->approximateNumberOfMers() + 1;

  if (args->memoryLimit) {
    args->mersPerBatch = estimateNumMersInMemorySize(args->merSize, args->memoryLimit, args->positionsEnabled, args->beVerbose);
    if (args->mersPerBatch > args->numMersActual)
      args->mersPerBatch = args->numMersActual;
    args->segmentLimit = (u64bit)ceil((double)args->numMersActual / (double)args->mersPerBatch);
    if (args->beVerbose)
      fprintf(stderr, "Have a memory limit: mersPerBatch="u64bitFMT" segmentLimit="u64bitFMT"\n", args->mersPerBatch, args->segmentLimit);
  } else if (args->segmentLimit) {
    args->mersPerBatch = (u64bit)ceil((double)args->numMersActual / (double)args->segmentLimit);
    if (args->beVerbose)
      fprintf(stderr, "Have a segment limit: mersPerBatch="u64bitFMT" segmentLimit="u64bitFMT"\n", args->mersPerBatch, args->segmentLimit);
  } else {
    args->mersPerBatch = args->numMersActual;
    args->segmentLimit = 1;
    if (args->beVerbose)
      fprintf(stderr, "Have NO LIMITS!: mersPerBatch="u64bitFMT" segmentLimit="u64bitFMT"\n", args->mersPerBatch, args->segmentLimit);
  }

  // store an array of the minium hash value for each sequence. Corresponds to min for each hash table
  unordered_map<u64bit, u64bit *> storedHash;
  // store an array of minum values to sequences for each hash table
  unordered_map<u64bit, vector<u64bit>* >* minToSequence = new unordered_map<u64bit, vector<u64bit> *>[args->numHashes + 1];
  // store an array of positions of the minimum kmer values
  unordered_map<u64bit, int >* minToPosition = new unordered_map<u64bit, int>[args->numHashes + 1];

  args->basesPerBatch = (u64bit)ceil((double)args->numBasesActual / (double)args->segmentLimit);
  args->bucketPointerWidth = logBaseTwo64(args->basesPerBatch + 1);
  args->numBuckets_log2    = 16; //optimalNumberOfBuckets(args->merSize, args->basesPerBatch, args->positionsEnabled);
  args->numBuckets         = (u64bitONE << args->numBuckets_log2);
  args->merDataWidth       = args->merSize * 2 - args->numBuckets_log2;

  args->salt_log2 = logBaseTwo64(args->numHashes * 10);
  const kMer *currMer = 0L;
  bool isFwd = true;

  // the kmers of sequences, read once and reused for all hash tables
  unordered_map<u64bit, vector<kMer> *> seqs;
  // length of sequences loaded
  unordered_map<u64bit, u32bit> seqLens;
  // orientation of the kmers from sequences
  unordered_map<u64bit, vector<bool> *> ori;
  // read evey input kmer and store whether it was forward or reverse, the kmer, and the length of the sequence
  while (M->nextMer()) {
     if (M->theFMer() <= M->theRMer()) {
        currMer = &(M->theFMer());
        isFwd = true;
     } else {
        currMer = &(M->theRMer());
        isFwd = false;
     }

     if (seqs.find(M->theSequenceNumber()) == seqs.end()) {
        seqs.insert(pair<u64bit, vector<kMer>*>(M->theSequenceNumber(), new vector<kMer>()));
        ori.insert(pair<u64bit, vector<bool>*>(M->theSequenceNumber(), new vector<bool>()));
     }
     seqs.find(M->theSequenceNumber())->second->push_back(*currMer);
     seqLens.erase(M->theSequenceNumber());
     seqLens.insert(pair<u64bit, u32bit>(M->theSequenceNumber(), M->thePositionInSequence()+args->merSize));
     ori.find(M->theSequenceNumber())->second->push_back(isFwd);
  }

  for (u32bit i = 0; i < args->numHashes; i++) {
     u32bit processed = 0;
     args->salt = i * 10;
     //fprintf(stderr, "Starting hash %d with min %llu\n", i, minHash);
     for (u32bit s = 0; s < seqs.size(); s++) {
        vector<kMer> *mers = seqs.find(s)->second;
        vector<bool> *oris = ori.find(s)->second;
        u64bit minHash = MAX_HASH-1;
        int minPos = 0;
        int currPos = 0;
        u64bit lastSeq = s;
        
        if (processed % 10000 == 0) { fprintf(stderr, "Done with %d sequences %d/%d hashes\n", processed, i+1, args->numHashes); }
        processed++;
        for (vector<kMer>::iterator mer = mers->begin(); mer != mers->end(); ++mer) {
           u64bit key = args->hash(*mer);
           if (key < minHash) {
              minHash = key; 
              minPos = (oris->at(currPos) == true ? currPos : -1*currPos);
              //fprintf(stderr, "For sequence %llu hash %llu stored min %llu at position %d\n", s, i, minHash, minPos);
           }
           currPos++;
        }
        //fprintf(stderr, "For sequence %llu hash %d min value is %llu\n", lastSeq, i, minHash);
        if (storedHash.find(lastSeq) == storedHash.end()) {
           storedHash.insert(pair<u64bit, u64bit *>(lastSeq, new u64bit[args->numHashes+1]));
        }
        storedHash.find(lastSeq)->second[i] = minHash;
        if (minToSequence[i].find(minHash) == minToSequence[i].end()) {
           minToSequence[i].insert(pair<u64bit, vector<u64bit> *>(minHash, new vector<u64bit>()));
         }
         minToSequence[i].find(minHash)->second->push_back(lastSeq);
         minToPosition[i].insert(pair<u64bit, int>(lastSeq, minPos)); 
     } 
  }
  seqs.clear();
  ori.clear();

  // now that we have the hash constructed, go through all sequences to recompute their min and score their matches
  for (unordered_map<u64bit, u64bit* >::iterator iter = storedHash.begin(); 
                                        iter != storedHash.end(); 
                                        ++iter) {
      unordered_map<u64bit, u64bit> counts;
      for (u32bit i = 0; i < args->numHashes; i++) {
         //fprintf(stderr, "For sequence %llu hash %d min value stored is %llu\n", iter->first, i, iter->second[i]);
         vector<u64bit>* sameMin = minToSequence[i].find(iter->second[i])->second;
         for (vector<u64bit>::iterator mins = sameMin->begin(); mins != sameMin->end(); ++mins) {
            if (counts.find(*mins) == counts.end()) { counts.insert(pair<u64bit, u64bit>(*mins, 0)); }
            counts[*mins] = counts[*mins] + 1;
            //fprintf(stderr, "Hash %d sequences %llu and %llu share minimum\n", i, iter->first, *mins);
         }
      }
      for (unordered_map<u64bit, u64bit* >::iterator c = storedHash.begin();
                                        c != storedHash.end() && c->first < iter->first;
                                        ++c) {
         double score = 0;
         if (counts.find(c->first) != counts.end()) {
            score = (double)counts[c->first] / args->numHashes * 100;
         }
         //fprintf(stderr, "For sequence %llu versus sequence %llu the count is %.2f\n", iter->first, c->first, score);

         if (score >= args->threshold) {
            // figure out intersecting regions
            int aStart = 0;
            int aEnd = 0;
            int aFwd = 0;
            int aRev = 0;
            int bStart = 0;
            int bEnd = 0;
            int bFwd = 0;
            int bRev = 0;

            u64bit* aHashes = storedHash.find(iter->first)->second;
            u64bit* bHashes = storedHash.find(c->first)->second;
            set<u64bit> aPositions;
            set<u64bit> bPositions;

            for (u32bit i = 0; i < args->numHashes; i++) {
               if (aHashes[i] != bHashes[i]) { continue; }
               int aPos = minToPosition[i].find(iter->first)->second;
               (aPos < 0 ? aRev++: aFwd++);
               aPos = abs(aPos);
               int bPos = minToPosition[i].find(c->first)->second;
               (bPos < 0 ? bRev++ : bFwd++);
               bPos = abs(bPos);
               aPositions.insert(aPos);
               bPositions.insert(bPos);
            }
            set<u64bit>::iterator bIter=bPositions.begin();
            for (set<u64bit>::iterator aIter = aPositions.begin(); aIter != aPositions.end() && bIter != bPositions.end(); ++bIter, ++aIter) {
               int aPos = *aIter;
               int bPos = *bIter;

               //fprintf(stderr, "Comparsing sequence %llu versus %llu. Position in first is %d second %d\n", iter->first, c->first, aPos, bPos);
               if (aStart != 0 && aStart > aPos && aPos - aStart > args->maxSkip) { continue; }
               if (aStart != 0 && aEnd < (aPos+args->merSize) && aPos+args->merSize - aEnd > args->maxSkip) { continue; }
               if (bStart != 0 && bStart > bPos && bPos - bStart > args->maxSkip) { continue; }
               if (bStart != 0 && bEnd < (bPos+args->merSize) && bPos+args->merSize - bEnd > args->maxSkip) { continue; }

               if (aStart == 0 || aStart > aPos) aStart = aPos;
               if (aEnd < aPos+args->merSize) aEnd = aPos+args->merSize;
               if (bStart == 0 || bStart > bPos) bStart = bPos;
               if (bEnd < bPos+args->merSize) bEnd = bPos+args->merSize;
            }
            bool aOri = (aFwd > aRev ? true : false);
            bool bOri = (bFwd > bRev ? true : false);
            fprintf(stderr, "For sequence %llu (%d-%d) %d versus sequence %llu (%d-%d) %d the count is %.2f\n", iter->first, (aOri ? aStart : aEnd), (aOri ? aEnd : aStart), seqLens.find(iter->first)->second, c->first, (bOri ? bStart : bEnd), (bOri ? bEnd : bStart), seqLens.find(c->first)->second, score);

            /*
            seq_pair p;
            u32bit cnt = 0;
            p.alen = seqstr->lengthOf(iter->first);
            p.a = new char[p.alen + 1];
            seqstr->setRange(seqstr->startOf(iter->first), seqstr->startOf(iter->first)+p.alen);
            while (!seqstr->eof()) {
               p.a[cnt] = seqstr->get();
               cnt++;
            }
            cnt = 0;
            p.blen = seqstr->lengthOf(c->first);
            p.b = new char[p.blen + 1];
            seqstr->setRange(seqstr->startOf(c->first), seqstr->startOf(c->first)+p.blen);
            while (!seqstr->eof()) {
               p.b[cnt] = seqstr->get();
               cnt++;
            }
            seq_pair_t result = smith_waterman(&p, true);
            fprintf(stderr, "%s\n%s\n", result->a, result->b);
            */
         }
      }
  }
  delete M;
  delete args;

  return(0);
}
