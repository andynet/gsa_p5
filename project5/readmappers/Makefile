
LIB_BINARIES = \
	aho_corasick.o            \
	cigar.o                   \
	fasta.o                   \
	match.o                   \
	string_utils.o            \
	bwt.o                     \
	fastq.o                   \
	remap.o                   \
	suffix_array.o            \
	generic_data_structures.o \
	sam.o                     \
	suffix_tree.o             \
	edit_distance_generator.o \
	io.o                      \
	serialise.o               \
	trie.o

# aho_corasick.c cigar.c fasta.c match.c string_utils.c bwt.c fastq.c remap.c suffix_array.c generic_data_structures.c sam.c suffix_tree.c edit_distance_generator.c io.c serialise.c trie.c aho_corasick.h cigar.h fasta.h match.h string_utils.h bwt.h fastq.h remap.h suffix_array.h generic_data_structures.h sam.h suffix_tree.h edit_distance_generator.h io.h serialise.h trie.h

all: bwt_readmapper exact_readmapper

libstralg.a: $(LIB_BINARIES)
	ar rcs $@ $^

bwt_readmapper: libstralg.a bwt_readmapper.o
	gcc -std=c11 $^ -o $@

exact_readmapper: libstralg.a exact_readmapper.o
	gcc -std=c11 $^ -o $@

clean:
	rm *.o *.a
