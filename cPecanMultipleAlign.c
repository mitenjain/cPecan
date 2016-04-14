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

int main(int argc, char *argv[]) {
    // Parse arguments
    if (argc != 2) {
        usage(argv);
        return 1;
    }

    // You would load a custom HMM here if you wanted using
    // hmm_getStateMachine (see the realign code)
    StateMachine *stateMachine  = stateMachine5_construct(fiveState);

    // From Benedict's code
    PairwiseAlignmentParameters *parameters = \
                                    pairwiseAlignmentBandingParameters_construct();

    // declare spanningTrees, maxPairsToConsider, useProgressiveMerging, matchGamma
    int64_t spanningTrees = 2;
    int64_t maxPairsToConsider = 1000;
    int useProgressiveMerging = 1;
    float matchGamma = 0.85;

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
        i++;
    }

    // Just a sanity check
    printf( "# reads %d \n", i );

    // Make a call to makeAlignment from MultipleAligner. This returns a column struct
    // the input params are just place holders to make this work and customizable later
    MultipleAlignment *mA = makeAlignment(stateMachine, seqFrags, spanningTrees, \
            maxPairsToConsider, useProgressiveMerging, matchGamma, parameters);

    // Just a sanity check
    printf("%" PRIi64 "\n", stSet_size(mA->columns));
    
    // call stSet_getIterate have iterate over columns, which are stSets
    // outer loop iterates over columns in multipleAlignment mA, 
    // each column is a struct Column
    // mA is not ordered, we need to figure out how to order
    stSetIterator *columns = stSet_getIterator(mA->columns);
    Column *column;
    while ((column = stSet_getNext(columns)) != NULL) {
        stList *columnNucleotides = stList_construct3(0, free);
        // now iterate over each column and get the seqName and position
        // get query sequence from seqFrags based on seqName, 
        // which is the position of querySeq in seqFrags in list
        while(column != NULL) {
            // individual column entry, which is a linked list itself
            SeqFrag *querySeq = stList_get(seqFrags, column->seqName);
            char *nucleotide = malloc(1);
            *nucleotide = querySeq->seq[column->position];
            stList_append(columnNucleotides, nucleotide);
            column = column->nColumn;
        }
        // consensus finding
        // Initialize an array with 0's and assign A, C, G, T to 0, 1, 2, 3
        // This just counts occurences of A, C, G, T
        int64_t nucleotideArray[4] = { 0, 0, 0, 0 } ;
        for (int64_t i = 0; i < stList_length(columnNucleotides); i++) {
            char *n;
            n = stList_get(columnNucleotides, i);
            if ( *n == 'A') { nucleotideArray[0]++ ; }
            else if ( *n == 'C') { nucleotideArray[1]++ ; }
            else if ( *n == 'G') { nucleotideArray[2]++ ; }
            else if ( *n == 'T') { nucleotideArray[3]++ ; }
        }
        
        // Now get the winner base, this is naive right now
        // need to add smart filters 
        // like winner means ≥ 50% reads and at least x number of reads
        int64_t winner ;
        winner = nucleotideArray[0];
        int64_t idx = 0 ;
        for (int64_t c = 0; c < sizeof(nucleotideArray)/sizeof(int64_t); c++) {
            if (nucleotideArray[c] > winner) {
                winner  = nucleotideArray[c];
                idx = c;
            }
        }

        // check which index is the highest, and assign a base
        // assign the nucleotide based on that
        if (idx == 0) { printf("A"); }
        if (idx == 1) { printf("C"); } 
        if (idx == 2) { printf("G"); } 
        if (idx == 3) { printf("T"); } 
 
        // cleanup 
        stList_destruct(columnNucleotides);
    }

    // This function is what I will try to use for sorting (later)
//             stList_sort(alignedPairs, (int (*)(const void *, const void *)) stIntTuple_cmpFn);
//             // Output the cigar string
//             cigarWrite(stdout, alignment, 0);

    // Clean up
    stHash_destructIterator(queryIt);
    stHash_destruct(querySequences);
    stList_destruct(seqFrags);
    pairwiseAlignmentBandingParameters_destruct(parameters);
    stateMachine_destruct(stateMachine);
    printf("\nDone\n");
}
