#include "sonLib.h"
#include "pairwiseAligner.h"
#include "multipleAligner.h"
#include "commonC.h"

static void usage(char *argv[]) {
//     fprintf(stderr, "%s fasta_target (ref) fasta_query\n", argv[0]);
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

//
//// reads columns from a multiple alignment
//stSet *readColumns(stList *columns) {
//    /*
// *      * Reads a set of columns, each containing one sequence position. Represents
// *           * initially unaligned state of sequence positions.
// *                */
//    //stSet *columns = stSet_construct3((uint64_t(*)(const void *)) column_hashFn, (int(*)(const void *, const void *)) column_equalsFn,
//    //        (void(*)(void *)) column_destruct);
//    for (int64_t seq = 0; seq < stList_length(seqFrags); seq++) {
//        //int64_t seqLength = ((SeqFrag *) (stList_get(seqFrags, seq)))->length;
//        for (int64_t pos = 0; pos < seqLength; pos++) {
//            Column *c = st_malloc(sizeof(Column));
//            c->seqName = seq;
//            c->position = pos;
//            c->nColumn = NULL;
//            stSet_insert(columns, c);
//        }
//    }
//    return columns;
//}
//
int main(int argc, char *argv[]) {
    // Parse arguments
    if (argc != 2) {
        usage(argv);
        return 1;
    }

    // You would load a custom HMM here if you wanted using
    // hmm_getStateMachine (see the realign code)
    StateMachine *stateMachine  = stateMachine5_construct(fiveState);

    // From Ben's code
    PairwiseAlignmentParameters *parameters = pairwiseAlignmentBandingParameters_construct();
    float matchGamma = 0.85;
    stList *seqFrags = stList_construct3(0, (void(*)(void *))seqFrag_destruct);

    // create hash of reads and iterate to construct seqFrags
    stHash *querySequences = readFastaFile(argv[1]);
    stHashIterator *queryIt = stHash_getIterator(querySequences);
    char *queryHeader;
    int i = 0;
    while ((queryHeader = stHash_getNext(queryIt)) != NULL) {
        char *querySeq = stHash_search(querySequences, queryHeader);
        stList_append(seqFrags, seqFrag_construct( querySeq, 0, strlen(querySeq) ));
        printf( "Added Read %d \n", i );
        i++;
    }

    // Make a call to the multiple alignment code from MultipleAligner.
    MultipleAlignment *mA = makeAlignment(stateMachine, seqFrags, 2, 1000, st_random() > 0.5, 0.85, parameters);

    for(int i = 0; i < stSet_size(mA); i++) {
        printf("%d\n", i);
    }
    
//             stList_sort(alignedPairs, (int (*)(const void *, const void *)) stIntTuple_cmpFn);
//             // Output the cigar string
//             cigarWrite(stdout, alignment, 0);

    stHash_destructIterator(queryIt);
    // Clean up
    stHash_destruct(querySequences);
    stList_destruct(seqFrags);
    pairwiseAlignmentBandingParameters_destruct(parameters);
    stateMachine_destruct(stateMachine);
    printf("\nDone\n");
}
