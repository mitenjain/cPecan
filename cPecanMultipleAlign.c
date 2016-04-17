#include "sonLib.h"
#include "pairwiseAligner.h"
#include "multipleAligner.h"
#include "commonC.h"

static void usage(char *argv[]) {
    fprintf(stderr, "%s fasta_query\n", argv[0]);
}

// Returns a hash mapping from sequence header to sequence data.
static stHash *readFastaFile(char *filename) {
    FILE *fasta = fopen(filename, "r");
    if (fasta == NULL) {
        st_errnoAbort("Could not open fasta file %s", filename);
    }
    stHash *headerToData = stHash_construct3(stHash_stringKey,
                                             stHash_stringEqualKey,
                                             free,
                                             free);
    struct List *seqs = constructEmptyList(0, NULL);
    struct List *seqLengths = constructEmptyList(0, free);
    struct List *headers = constructEmptyList(0, free);
    fastaRead(fasta, seqs, seqLengths, headers);

    for (int64_t i = 0; i < seqs->length; i++) {
        char *fullHeader = headers->list[i];
        stList *headerTokens = stString_splitByString(fullHeader, " ");
        char *usableHeader = stString_copy(stList_get(headerTokens, 0));
        stHash_insert(headerToData, usableHeader, seqs->list[i]);
        stList_destruct(headerTokens);
    }
    destructList(seqs);
    destructList(seqLengths);
    destructList(headers);

    return headerToData;
}

/*
 *
 * Following functions are used to sort the multiple alignment columns
 *
 */

int cmpPositions(const void *a, const void *b) {
    Column *c1 = (Column *)a;
    Column *c2 = (Column *)b;
    assert(c1->seqName != c2->seqName);
    return c1->seqName < c2->seqName ? -1 : 1;
}

Column *sortColumn(Column *c) {
    stList *columnList = stList_construct();
    while(c != NULL) {
        stList_append(columnList, c);
        c = c->nColumn;
    }
    stList_sort(columnList, cmpPositions);
    for(int64_t i=0; i+1<stList_length(columnList); i++) {
        Column *c1 = stList_get(columnList, i);
        Column *c2 = stList_get(columnList, i+1);
        c1->nColumn = c2;
        c2->nColumn = NULL;
    }
    c = stList_get(columnList, 0);
    stList_destruct(columnList);
    return c;
}

int cmpColumnFn(const void *a, const void *b) {
    Column *c1 = (Column *)a;
    Column *c2 = (Column *)b;
    while(1) {
        if(c1 == NULL && c2 == NULL) {
            return 0;
        }
        if(c1 == NULL) {
            return -1;
        }
        if(c2 == NULL) {
            return 1;
        }
        if(c1->seqName == c2->seqName) {
            assert(c1->position != c2->position);
            return c1->position < c2->position ? -1 : 1;
        }
        if(c1->seqName < c2->seqName) {
            c1 = c1->nColumn;
        }
        else {
            c2 = c2->nColumn;
        }
    }
    return 0;
}

stList *getSortedColumnList(stSet *columns) {
    stList *columnList = stSet_getList(columns);
    printf("Got the column list: %" PRIi64 "\n", stList_length(columnList));
    for(int64_t i=0; i<stList_length(columnList); i++) {
        stList_set(columnList, i, sortColumn(stList_get(columnList, i)));
    }
    printf("Ordered each columns: %" PRIi64 "\n", stList_length(columnList));

    stList_sort(columnList, cmpColumnFn);
    return columnList;
}

int main(int argc, char *argv[]) {
    // Parse arguments
    if (argc != 4) {
        usage(argv);
        return 1;
    }

    // You would load a custom HMM here if you wanted using
    // hmm_getStateMachine (see the realign code) - this should use one of the nanopore HMMs.
    StateMachine *stateMachine  = stateMachine5_construct(fiveState);

    // From Benedict's code
    PairwiseAlignmentParameters *parameters = \
                                    pairwiseAlignmentBandingParameters_construct();

    // declare spanningTrees, maxPairsToConsider, useProgressiveMerging, matchGamma
    int64_t spanningTrees = 4;
    int64_t maxPairsToConsider = 10000000;
    int useProgressiveMerging = 0;
    float matchGamma = 0.0;

    // initialize seqFrags
    stList *seqFrags = stList_construct3(0, (void(*)(void *))seqFrag_destruct);

    // create hash of reads and iterate to construct seqFrags
    stHash *querySequences = readFastaFile(argv[1]);
    stHashIterator *queryIt = stHash_getIterator(querySequences);
    char *queryHeader;
    int i = 0;
    while ((queryHeader = stHash_getNext(queryIt)) != NULL) {
        char *querySeq = stHash_search(querySequences, queryHeader);
        stList_append(seqFrags, seqFrag_construct( querySeq, 0, strlen(querySeq) ));
        printf("Adding a sequence of length: %" PRIi64 "\n", strlen(querySeq));
        i++;
    }

    // Just a sanity check
    printf( "# reads %d \n", i );

    // Make a call to makeAlignment from MultipleAligner. This returns a column struct
    // the input params are just place holders to make this work and customizable later
    MultipleAlignment *mA = makeAlignment(stateMachine, seqFrags, spanningTrees, \
            maxPairsToConsider, useProgressiveMerging, matchGamma, parameters);

    // Just a sanity check
    printf("Got %" PRIi64 " columns\n", stSet_size(mA->columns));
    
    // call stSet_getIterate have iterate over columns, which are stSets
    // outer loop iterates over columns in multipleAlignment mA, 
    // each column is a struct Column
    // mA is not ordered, we need to figure out how to order
    stList *columnList = getSortedColumnList(mA->columns);

    printf("Sorted the columns: %" PRIi64 "\n", stList_length(columnList));

    char *consensusSeq = st_malloc(stList_length(columnList)+1);
    int64_t consensusSeqLength = 0;

    for(int64_t i=0; i<stList_length(columnList); i++) {
        Column *column = stList_get(columnList, i);
        stList *columnNucleotides = stList_construct3(0, free);
        // now iterate over each column and get the seqName and position
        // get query sequence from seqFrags based on seqName, 
        // which is the position of querySeq in seqFrags in list
        // and count the number of As, Cs, Ts, and Gs
        int64_t nucleotideArray[4] = { 0, 0, 0, 0 } ;
        while(column != NULL) {
            // individual column entry, which is a linked list itself
            SeqFrag *querySeq = stList_get(seqFrags, column->seqName);
            char n = toupper(querySeq->seq[column->position]);
            if ( n == 'A') { nucleotideArray[0]++ ; }
            else if ( n == 'C') { nucleotideArray[1]++ ; }
            else if ( n == 'G') { nucleotideArray[2]++ ; }
            else if ( n == 'T') { nucleotideArray[3]++ ; }
            else assert(0);
            column = column->nColumn;
        }
        
        // Now get the winner base, this is naive right now
        // need to add smart filters 
        // like winner means ≥ 50% reads and at least x number of reads
        int64_t total = 0;
        int64_t winner = nucleotideArray[0];
        int64_t idx = 0 ;
        for (int64_t c = 0; c < 4; c++) {
            total += nucleotideArray[c];
            if (nucleotideArray[c] > winner) {
                winner = nucleotideArray[c];
                idx = c;
            }
        }

        // check which index is the highest, and assign a base
        // assign the nucleotide based on that
        if(total >= stList_length(seqFrags)-2) {
            consensusSeq[consensusSeqLength++] = idx == 0 ? 'A' : (idx == 1 ? 'C' : (idx == 2 ? 'G' : 'T'));
        }
 
        // cleanup 
        stList_destruct(columnNucleotides);
    }
    consensusSeq[consensusSeqLength] = '\0';

    FILE *fH = fopen(argv[2], "w");
    fastaWrite(consensusSeq, "consensus_seq", fH);
    fclose(fH);

    // This function is what I will try to use for sorting (later)
//             stList_sort(alignedPairs, (int (*)(const void *, const void *)) stIntTuple_cmpFn);
//             // Output the cigar string
//             cigarWrite(stdout, alignment, 0);

    //
    //consensusSeq = ((SeqFrag *)stList_get(seqFrags, 2))->seq;

    //Now load up the reference sequence and compare it to the consensus
    char *refSeq = stList_get(stHash_getValues(readFastaFile(argv[3])), 0);
    //Reverse complement the sequence
    refSeq = stString_reverseComplementString(refSeq);


    printf("Loaded the reference comparison sequence, has length: %" PRIi64 " \n", strlen(refSeq));
    stList *alignedPairs = getAlignedPairs(stateMachine, refSeq, consensusSeq, parameters,  0, 0);
    printf("All aligned pairs: %" PRIi64 "\n", stList_length(alignedPairs));
    alignedPairs = filterPairwiseAlignmentToMakePairsOrdered(alignedPairs, refSeq, consensusSeq, matchGamma);
    printf("Aligned pairs after filtering: %" PRIi64 "\n", stList_length(alignedPairs));
    //Calc identity stats
    int64_t identicalAlignedPairs = 0;
    for(int64_t i=0; i<stList_length(alignedPairs); i++) {
        stIntTuple *alignedPair = stList_get(alignedPairs, i);
        identicalAlignedPairs += refSeq[stIntTuple_get(alignedPair, 1)] == consensusSeq[stIntTuple_get(alignedPair, 2)];
    }
    printf("Aligned pairs %" PRIi64 ", of which %" PRIi64 " are identical, giving an identity of %f\n",
            stList_length(alignedPairs), identicalAlignedPairs, identicalAlignedPairs * 2.0 / (strlen(consensusSeq) + strlen(refSeq)));

    // Clean up
    stHash_destructIterator(queryIt);
    stHash_destruct(querySequences);
    stList_destruct(seqFrags);
    pairwiseAlignmentBandingParameters_destruct(parameters);
    stateMachine_destruct(stateMachine);
    printf("\nDone, reported %" PRIi64 "columns\n", consensusSeqLength);
}