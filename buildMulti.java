import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.TreeSet;

public class buildMulti {
   private static final long MAX_HASH = Long.MAX_VALUE;
   private static final int DEFAULT_K = 10;
   private static final int DEFAULT_HASHES = 100;
   private static final int DEFAULT_THRESHOLD = 1;
   private static final int DEFAULT_SKIP = 500;

   private static final String[] fastaSuffix = {"fna", "contigs", "final", "fasta", "fa"};

   private static long getHash(String kmer, long salt) {
      long key = 5381 + kmer.hashCode();
      
    key = key * 37;
    key = key + salt;
    key ^= key >> 33;
    key *= 0xff51afd7ed558ccdL;
    key ^= key >> 33;
    key *= 0xc4ceb9fe1a85ec53L;
    key ^= key >> 33;
    return(key);
  }

   public buildMulti() {
   }


   // the hashmap of kmers in the sequences (kmer list for each sequence)
   private HashMap<Long, ArrayList<String>> seqs = new HashMap<Long, ArrayList<String>>();
   // length of sequences loaded
   private HashMap<Long, Integer> seqLens = new HashMap<Long,Integer>();
   // orientation of the kmers from sequences
   private HashMap<Long, ArrayList<Boolean>> ori = new HashMap<Long, ArrayList<Boolean>>();

   // an array list for every hash function to store the sequences maching a given minimum
   private ArrayList<HashMap<Long, ArrayList<Long>>> minToSequence = new ArrayList<HashMap<Long, ArrayList<Long>>>();
   // an array list for every sequence to the position of the minium hash 
   private ArrayList<HashMap<Long, Integer>> minToPosition = new ArrayList<HashMap<Long, Integer>>();

   // a map for each sequence to the list of minimum hashes for each table
   private HashMap<Long, long[]> storedHash = new HashMap<Long, long[]>();

   // process a sequence and store the kmers/sequence length/etc
   private void addMers(long header, String seq, int merSize) {
      seqs.put(header, new ArrayList<String>());
      ori.put(header, new ArrayList<Boolean>());
      seqLens.put(header, seq.length());

      for (int i = 0; i <= seq.length() - merSize; i++) {
         String fmer = seq.substring(i, i+merSize);
         String rmer = Utils.rc(fmer);
         String currMer = null;
         boolean isFwd = true;
         if (fmer.compareTo(rmer) <= 0) {
            currMer = fmer;
            isFwd = true;
         } else {
            currMer = rmer;
            isFwd = false;
         }

         // store the canonical kmer and the orientation (i.e. whether cannonical was fwd or rev)
         seqs.get(header).add(currMer);
         ori.get(header).add(isFwd);
      }
  }


  // process a sequence given the hash we are on and the salt 
   private void processSeq(long header, int numHashes, int hashNum, long salt) {
      long minHash = MAX_HASH;
      int minPos = 0;
      int currPos = 0;
      ArrayList<String> mers = seqs.get(header);
      ArrayList<Boolean> oris = ori.get(header);

      for (String currMer : mers) {
         long key = getHash(currMer, salt);
         if (key < minHash) {
            minHash = key;
            minPos = (oris.get(currPos) == true ? currPos : -1*currPos);
            //System.err.println("For sequence " + header + " hash " + hashNum + " stored min " + minHash + " at position " + minPos);
         }
         currPos++;
      }
      if (!storedHash.containsKey(header)) {
         storedHash.put(header, new long[numHashes + 1]);
      }
      storedHash.get(header)[hashNum] = minHash;
      if (!minToSequence.get(hashNum).containsKey(minHash)) {
         minToSequence.get(hashNum).put(minHash, new ArrayList<Long>());
      }
      minToSequence.get(hashNum).get(minHash).add(header);
      minToPosition.get(hashNum).put(header, minPos);
   }

   public static void printUsage(String error) {
      if (error != null) {
         System.err.println(error);
      }
      System.err.println("Usage buildMulti <-s fasta file>");
      System.err.println("Options: ");
      System.err.println("\t -k [int merSize], default: " + DEFAULT_K);
      System.err.println("\t  --num-hashes [int # hashes], default: " + DEFAULT_HASHES);
      System.err.println("\t  --threshold [int threshold for % matching minimums], default: " + DEFAULT_THRESHOLD);
      System.err.println("\t --max-skip [int bp maximum distance to nearest minimum value when guessing overlap positions], default: " + DEFAULT_SKIP);
      System.exit(1);
   }

   public static void main(String[] args) throws Exception {
      String inFile = null;
      int kmerSize = DEFAULT_K;
      int threshold = DEFAULT_THRESHOLD;
      int numHashes = DEFAULT_HASHES; 
      int maxSkip = DEFAULT_SKIP;

      for (int i = 0; i < args.length; i++) {
         if (args[i].trim().equalsIgnoreCase("-k")) {
             kmerSize = Integer.parseInt(args[++i]);
         } else if (args[i].trim().equalsIgnoreCase("-s")) {
             inFile = args[++i];
         } else if (args[i].trim().equalsIgnoreCase("--num-hashes")) {
             numHashes = Integer.parseInt(args[++i]);
         } else if (args[i].trim().equalsIgnoreCase("--threshold")) {
             threshold = Integer.parseInt(args[++i]);
         } else if (args[i].trim().equalsIgnoreCase("--max-skip")) {
             maxSkip = Integer.parseInt(args[++i]);
         }
      }
      if (inFile == null) {
         printUsage("Error: no input fasta file specified");
      }

      buildMulti b = new buildMulti();

      System.err.println("Running with input fasta: " + inFile);
      System.err.println("kmer size:\t" + kmerSize);
      System.err.println("threshold:\t" + threshold);
      System.err.println("num hashes\t" + numHashes);
      System.err.println("max skip\t" + maxSkip);

      // read and index the kmers
      long startTime = System.nanoTime();
      BufferedReader bf = Utils.getFile(inFile, fastaSuffix);
      String line = null;
      StringBuilder fastaSeq = new StringBuilder();
      String header = "";
      long numProcessed = 0;
      while ((line = bf.readLine()) != null) {
         if (line.startsWith(">")) { 
            if (fastaSeq.length() > 0) b.addMers(numProcessed++, fastaSeq.toString().toUpperCase(), kmerSize); 
            fastaSeq.setLength(0);
         } else {
            fastaSeq.append(line);
         }
      }
      if (fastaSeq.length() != 0) { b.addMers(numProcessed++, fastaSeq.toString().toUpperCase(), kmerSize); }
      bf.close();
      System.err.println("Time to read: " + (System.nanoTime() - startTime)/1000000);

      // compute hashes
      startTime = System.nanoTime();
      for (int i = 0; i < numHashes; i++) {
         b.minToSequence.add(new HashMap<Long, ArrayList<Long>>());
         b.minToPosition.add(new HashMap<Long, Integer>());
         long salt = i * 10;

         for (int s = 0; s < b.seqs.size(); s++) {
               b.processSeq(s, numHashes, i, salt); 
         }
      }
      System.err.println("Time to hash: " + (System.nanoTime() - startTime)/1000000);

      // now that we have the hash constructed, go through all sequences to recompute their min and score their matches
      startTime = System.nanoTime();
      for (long seq : b.storedHash.keySet()) {
         long[] storedMins = b.storedHash.get(seq);
         HashMap<Long, Integer> counts = new HashMap<Long, Integer>();
         for (int i = 0; i < numHashes; i++) {
            ArrayList<Long> sameMin = b.minToSequence.get(i).get(storedMins[i]);
            for (long mins : sameMin) {
               if (!counts.containsKey(mins)) { counts.put(mins, 0); }
               counts.put(mins, counts.get(mins) + 1);
            }
         }
         for (long seq2 : b.storedHash.keySet()) {
            if (seq2 >= seq) { break; }
            double score = 0;
            if (counts.containsKey(seq2)) {
               score = (double)counts.get(seq2) / numHashes * 100;
            }

            if (score >= threshold) {
               // figure out intersecting regions
               int aStart = 0;
               int aEnd = 0;
               int aFwd = 0;
               int aRev = 0;
               int bStart = 0;
               int bEnd = 0;
               int bFwd = 0;
               int bRev = 0;

               long[] aHashes = b.storedHash.get(seq);
               long[] bHashes = b.storedHash.get(seq2);
               TreeSet<Integer> aPositions = new TreeSet<Integer>();
               TreeSet<Integer> bPositions = new TreeSet<Integer>();

               for (int i = 0; i < numHashes; i++) {
                  if (aHashes[i] != bHashes[i]) { continue; }
                  int aPos = b.minToPosition.get(i).get(seq);
                  if (aPos < 0) { aRev++; } else { aFwd++; }
                  aPos = Math.abs(aPos);
                  int bPos = b.minToPosition.get(i).get(seq2);
                  if (bPos < 0) { bRev++; } else { bFwd++; }
                  bPos = Math.abs(bPos);
                  aPositions.add(aPos);
                  bPositions.add(bPos);
               }

               Iterator bIter = bPositions.iterator();
               for (Iterator aIter = aPositions.iterator(); aIter.hasNext() && bIter.hasNext(); ) {
                  int aPos = (Integer)aIter.next();
                  int bPos = (Integer)bIter.next();
                  //System.err.println("Comparsing sequence " + seq + "  versus " + seq2 + ". Position in first is " + aPos + " second " + bPos);

                  if (aStart != 0 && aStart > aPos && aPos - aStart > maxSkip) { continue; }
                  if (aStart != 0 && aEnd < (aPos+kmerSize) && aPos+kmerSize - aEnd > maxSkip) { continue; }
                  if (bStart != 0 && bStart > bPos && bPos - bStart > maxSkip) { continue; }
                  if (bStart != 0 && bEnd < (bPos+kmerSize) && bPos+kmerSize - bEnd > maxSkip) { continue; }

                  if (aStart == 0 || aStart > aPos) aStart = aPos;
                  if (aEnd < aPos+kmerSize) aEnd = aPos+kmerSize;
                  if (bStart == 0 || bStart > bPos) bStart = bPos;
                  if (bEnd < bPos+kmerSize) bEnd = bPos+kmerSize;
               }
               boolean aOri = (aFwd > aRev ? true : false);
               boolean bOri = (bFwd > bRev ? true : false);
            
               System.out.println("For sequence " + seq + " (" + (aOri ? aStart : aEnd) + ", " + (aOri ? aEnd : aStart) + ") " + b.seqLens.get(seq) + " versus sequence " + seq2 + " (" + (bOri ? bStart : bEnd) + ", " + (bOri ? bEnd : bStart) + ") " + b.seqLens.get(seq2) + " the count is " + score);
            }
         }
     }
     System.err.println("Time to score: " + (System.nanoTime() - startTime)/1000000);
  }
}
