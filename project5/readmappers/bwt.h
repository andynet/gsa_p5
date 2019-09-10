
#ifndef BWT_H
#define BWT_H

#include "remap.h"
#include "suffix_array.h"

#include <stdbool.h>
#include <stdint.h>

struct bwt_table {
    struct remap_table  *remap_table;
    struct suffix_array *sa;
    uint32_t *c_table;
    uint32_t *o_table;
};

static inline size_t o_index(unsigned char a, size_t i,
                             const struct bwt_table *table)
{
    return a * table->sa->length + i;
}

// for these to work, sa must have been build
// from a remapped string.
void init_bwt_table   (struct bwt_table    *bwt_table,
                       struct suffix_array *sa,
                       struct remap_table  *remap_table);
void dealloc_bwt_table(struct bwt_table *bwt_table);

// This function frees the remap table and the suffix
// array as well as the BWT tables.
void dealloc_complete_bwt_table(struct bwt_table *bwt_table);

struct bwt_table *alloc_bwt_table(struct suffix_array *sa,
                                  struct remap_table  *remap_table);
void free_bwt_table(struct bwt_table *bwt_table);
void free_complete_bwt_table(struct bwt_table *bwt_table);

// This function builds all the structures needed to work
// with the bwt tables. It should be "completely" freed
// when no longer used.
struct bwt_table *build_complete_table(const char *string);

struct bwt_exact_match_iter {
    const struct suffix_array *sa;
    size_t L;
    size_t i;
    size_t R;
};
struct bwt_exact_match {
    size_t pos;
};
void init_bwt_exact_match_iter   (struct bwt_exact_match_iter *iter,
                                  struct bwt_table            *bwt_table,
                                  const char                  *remapped_pattern);
bool next_bwt_exact_match_iter   (struct bwt_exact_match_iter *iter,
                                  struct bwt_exact_match      *match);
void dealloc_bwt_exact_match_iter(struct bwt_exact_match_iter *iter);


struct bwt_approx_iter {
    struct bwt_approx_match_internal_iter *internal_approx_iter;
    struct bwt_exact_match_iter *internal_exact_iter;
    bool outer;
};
struct bwt_approx_match {
    const char *cigar;
    size_t match_length;
    size_t position;
};
void init_bwt_approx_iter(struct bwt_approx_iter *iter,
                          struct bwt_table       *bwt_table,
                          const char             *remapped_pattern,
                          int                     edits);
bool next_bwt_approx_match(struct bwt_approx_iter  *iter,
                           struct bwt_approx_match *match);
void dealloc_bwt_approx_iter(struct bwt_approx_iter *iter);


// Serialisation
void write_bwt_table(FILE *f, const struct bwt_table *bwt_table);
void write_bwt_table_fname(const char *fname, const struct bwt_table *bwt_table);

struct bwt_table *read_bwt_table(FILE *f,
                                 struct suffix_array *sa,
                                 struct remap_table  *remap_table);
struct bwt_table * read_bwt_table_fname(const char *fname,
                                        struct suffix_array *sa,
                                        struct remap_table  *remap_table);



// Some debug code
void print_c_table  (struct bwt_table *table);
void print_o_table  (struct bwt_table *table);
void print_bwt_table(struct bwt_table *table);

// FIXME: maybe change the name. I am testing for
// equivalence, not whether the two tables point to
// the same object or whether the underlying
// suffix array and remap tables are the
// same.
bool identical_bwt_tables(struct bwt_table *table1,
                          struct bwt_table *table2);

#endif
