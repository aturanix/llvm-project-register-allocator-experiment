#ifndef LLVM_LIB_CODEGEN_REGALLOCCHAITINBITEASTDGRAPH_H
#define LLVM_LIB_CODEGEN_REGALLOCCHAITINBITEASTDGRAPH_H

#include <stdbool.h>
#include <inttypes.h>


#define block_t uint64_t
#define BLOCK_INDEX(bit_index)              ((bit_index)/(sizeof(block_t)*8))
#define MASK_INDEX(bit_index)               ((bit_index)%(sizeof(block_t)*8))
#define MASK(bit_index)                     ((block_t)1 << MASK_INDEX(bit_index))
#define TOTAL_BLOCK_NUM(vertex_num)         (BLOCK_INDEX(vertex_num-1)+1)

#define SET_EDGE(vertex1, vertex2, edges)   edges[vertex1][BLOCK_INDEX(vertex2)] |=  MASK(vertex2); edges[vertex2][BLOCK_INDEX(vertex1)] |=  MASK(vertex1);

#define CHECK_COLOR(color, vertex)      (color[BLOCK_INDEX(vertex)] & MASK(vertex))
#define SET_COLOR(color, vertex)        color[BLOCK_INDEX(vertex)] |=  MASK(vertex)
#define RESET_COLOR(color, vertex)      color[BLOCK_INDEX(vertex)] &= ~MASK(vertex)


bool read_graph(const char* filename, int graph_size, block_t *edges, int offset_i);

bool read_weights(const char* filename, int size, int weights[]);

bool is_valid(
    int graph_size, 
    const block_t *edges, 
    int color_num, 
    const block_t *colors
);

int count_edges(int graph_size, const block_t *edges, int degrees[]);

void print_colors(
    const char *filename, 
    const char *header, 
    int color_num, 
    int graph_size, 
    const block_t *colors
);

int graph_color_greedy(
    int graph_size, 
    const block_t edges[][TOTAL_BLOCK_NUM(graph_size)], 
    block_t colors[][TOTAL_BLOCK_NUM(graph_size)], 
    int max_color_possible
);

void pop_complex_random (
    int graph_size, 
    const block_t *edges, 
    const int *weights,
    int pop_size,
    block_t **population, 
    int max_color
);

int count_conflicts(
    int graph_size, 
    const block_t *color, 
    const block_t *edges, 
    int *conflict_count
);

int popcountl(uint64_t n);


#endif
