#include <ctype.h>
#include "sonLib.h"
#include "pairwiseAligner.h"
#include "multipleAligner.h"
#include "commonC.h"

static void usage(char *argv[]) {
    fprintf(stderr, "%s fasta_query cns.fa ref.fa [orientation matters] \n", argv[0]);
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

// Get the entry in this column that has a column->seqName entry
// matching seqName, or NULL if there is none.
Column *getEntryForSeq(Column *c, int64_t seqName) {
    while (c != NULL) {
        if (c->seqName == seqName) {
            break;
        }
        c = c->nColumn;
    }
    return c;
}

// Compare two columns by a given sequence ID, pointed to by seqNamePtr.
int cmpColumnsBySeq(const void *a, const void *b, const void *seqNamePtr) {
    int64_t seqName = *((int64_t *) seqNamePtr);
    Column *c1 = (Column *)a;
    Column *c2 = (Column *)b;
    c1 = getEntryForSeq(c1, seqName);
    c2 = getEntryForSeq(c2, seqName);
    if (c1 == NULL && c2 == NULL) {
        return 0;
    } else if (c1 == NULL) {
        return 1;
    } else if (c2 == NULL) {
        return -1;
    } else {
        // Both columns have an entry for this sequence
        return c1->position > c2->position ? 1 : -1;
    }
}

// Print out a column in a nice way.
void debugPrintColumn(Column *column, stList *seqFrags) {
    while (column != NULL) {
        printf("name: %" PRIi64 " pos: %" PRIi64 " nuc: %c\n", column->seqName, column->position,
               ((SeqFrag *) stList_get(seqFrags, column->seqName))->seq[column->position]);
        column = column->nColumn;
    }
}

// Find the next adjacency to use, given a number of possible adjacencies.
// The next adjacency to use is chosen so that the ordering never has
// a "jump" of more than 1 in any sequence.
// A topological sort would work as well, but this ended up being easier to debug.
Column *findNextAdjacency(stHash *seqIndexToLargestPosition,
                          stSet *nextColumns) {
    stSetIterator *adjacencyIt = stSet_getIterator(nextColumns);

    Column *adjacency;
    while ((adjacency = stSet_getNext(adjacencyIt)) != NULL) {
        bool validAdjacency = true;
        Column *entry = adjacency;
        while (entry != NULL) {
            // FIXME: this is a temporary hacky way to search hash by
            // integer, it's better to use the stIntTuples but it is
            // slower and takes more code.
            int64_t *largestPosition = stHash_search(seqIndexToLargestPosition, (void *) entry->seqName);
            if (largestPosition == NULL) {
                largestPosition = malloc(sizeof(int64_t));
                *largestPosition = 0;
                stHash_insert(seqIndexToLargestPosition, (void *) entry->seqName, largestPosition);
            }
            int64_t jump = getEntryForSeq(adjacency, entry->seqName)->position - *largestPosition;

            if (jump < 0) {
                st_errAbort("going backwards");
            }
            if (jump > 1) {
                validAdjacency = false;
            }
            entry = entry->nColumn;
        }
        if (validAdjacency) {
            entry = adjacency;
            while (entry != NULL) {
                int64_t *largestPosition = stHash_search(seqIndexToLargestPosition, (void *) entry->seqName);
                *largestPosition = entry->position;
                entry = entry->nColumn;
            }
            stSet_destructIterator(adjacencyIt);
            return adjacency;
        }
    }
    stSet_destructIterator(adjacencyIt);
    return NULL;
}

stList *getSortedColumnList(stSet *columns, stList *seqFrags) {
    stList *columnList = stSet_getList(columns);
    // FIXME: don't need this but it is nice for debug prints
    printf("Got the column list: %" PRIi64 "\n", stList_length(columnList));
    for(int64_t i=0; i<stList_length(columnList); i++) {
        stList_set(columnList, i, sortColumn(stList_get(columnList, i)));
    }
    printf("Ordered each columns: %" PRIi64 "\n", stList_length(columnList));

    stHash *columnEntryToRightAdjacency = stHash_construct();
    stSet *startCols = stSet_construct();

    // Construct the adjacency hash.
    for (int64_t i = 0; i < stList_length(seqFrags); i++) {
        stList_sort2(columnList, cmpColumnsBySeq, &i);
        stSet_insert(startCols, stList_get(columnList, 0));
        for (int64_t j = 0; j < stList_length(columnList) - 1; j++) {
            if (getEntryForSeq(stList_get(columnList, j + 1), i)) {
                stHash_insert(columnEntryToRightAdjacency,
                              getEntryForSeq(stList_get(columnList, j), i),
                              stList_get(columnList, j + 1));
                if (stSet_search(startCols, stList_get(columnList, j + 1))) {
                    stSet_remove(startCols, stList_get(columnList, j + 1));
                }
            }
        }
    }

    stList_destruct(columnList);

    // From each start column (a column that has no left-adjacencies),
    // add the columns to sortedList in a way that maintains the partial
    // ordering
    stList *sortedList = stList_construct();
    stSetIterator *startColIt = stSet_getIterator(startCols);
    Column *startCol;
    while ((startCol = stSet_getNext(startColIt)) != NULL) {
        stHash *seqIndexToLargestPosition = stHash_construct();
        stSet *adjacencies = stSet_construct();
        Column *column = startCol;
        while (column != NULL) {
            stList_append(sortedList, column);
            Column *entry = column;
            while (entry != NULL) {
                stSet_insert(adjacencies, stHash_search(columnEntryToRightAdjacency, entry));
                entry = entry->nColumn;
            }
            /* debugPrintColumn(canonicalColumn, seqFrags); */
            /* printf("got %" PRIi64 " adjacencies\n", stSet_size(adjacencies)); */
            /* stSetIterator *adjacencyIt = stSet_getIterator(adjacencies); */
            /* Column *adjacency; */
            /* printf("==========================\n"); */
            /* while ((adjacency = stSet_getNext(adjacencyIt)) != NULL) { */
            /*     printf("adjacency:\n"); */
            /*     debugPrintColumn(adjacency, seqFrags); */
            /* } */

            column = findNextAdjacency(seqIndexToLargestPosition, adjacencies);
            stSet_remove(adjacencies, column);

            /* printf("chosen adjacency:\n"); */
            /* debugPrintColumn(canonicalColumn, seqFrags); */
        }
        stSet_destruct(adjacencies);
    }
    stSet_destructIterator(startColIt);

    printf("Had %" PRIi64 " start cols.\n", stSet_size(startCols));

    printf("Ordered columns and have %" PRIi64 " entries.\n", stList_length(sortedList));

    stHash_destruct(columnEntryToRightAdjacency);
    stSet_destruct(startCols);

    return sortedList;
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

    FILE *columnLengthDistribution = fopen("columnLengthDistribution", "w");
    
    // call stSet_getIterate have iterate over columns, which are stSets
    // outer loop iterates over columns in multipleAlignment mA, 
    // each column is a struct Column
    // mA is not ordered, we need to figure out how to order
    stList *columnList = getSortedColumnList(mA->columns, seqFrags);

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
        fprintf(columnLengthDistribution, "%" PRIi64 "\n", total);
        if(total >= stList_length(seqFrags)/2) {
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
    refSeq = stString_reverseComplementString(refSeq);

    printf("Loaded the reference comparison sequence, has length: %" PRIi64 " \n", strlen(refSeq));
    stList *alignedPairs = getAlignedPairs(stateMachine, refSeq, consensusSeq, parameters,  0, 0);
    printf("All aligned pairs: %" PRIi64 "\n", stList_length(alignedPairs));
    alignedPairs = reweightAlignedPairs2(alignedPairs, strlen(refSeq),
                                         strlen(consensusSeq),
                                         parameters->gapGamma);
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
    multipleAlignment_destruct(mA);
    stList_destruct(columnList);

    printf("\nDone, reported %" PRIi64 "columns\n", consensusSeqLength);
}
