#include "asrml.h"

int dir_a_to_b(Node* a, Node* b) {
	/* this returns the direction from a to b when a and b are two neighbours, otherwise yields an error */
	int i, n = a->nneigh;
	for(i=0; i<n; i++) if (a->neigh[i] == b) break;
	if (i < n) return i; else {
	  fprintf(stderr,"Fatal error : nodes are not neighbours.\n");
	  Generic_Exit(__FILE__,__LINE__,__FUNCTION__,EXIT_FAILURE);
	}
	return -1;
} /* end dir_a_to_b */

void write_subtree_to_stream(Node* node, Node* node_from, FILE* stream, double scale) {
	int i, direction_to_exclude, n = node->nneigh;
	if (node == NULL || node_from == NULL) return;

	if(n == 1) {
		/* terminal node */
		fprintf(stream, "%s:%f", (node->name ? node->name : ""), node->br[0]->brlen * scale); /* distance to father */
	} else {
	        direction_to_exclude = dir_a_to_b(node, node_from);	

		putc('(', stream);
		/* we have to write (n-1) subtrees in total. The last print is not followed by a comma */
		for(i=1; i < n-1; i++) {
			write_subtree_to_stream(node->neigh[(direction_to_exclude+i) % n], node, stream, scale); /* a son */
			putc(',', stream);
		}
		write_subtree_to_stream(node->neigh[(direction_to_exclude+i) % n], node, stream, scale); /* last son */
		putc(')', stream);
		fprintf(stream, "%s:%f", (node->name ? node->name : ""), node->br[0]->brlen * scale); /* distance to father */
	}

} /* end write_subtree_to_stream */

void write_nh_tree(Tree* tree, FILE* stream, double scale) {
	/* writing the tree from the current position in the stream */
	if (!tree) return;
        scale=1.0;
	Node* node = tree->node0; /* root or pseudoroot node */
	int i, n = node->nneigh;
	putc('(', stream);
	for(i=0; i < n-1; i++) {
		write_subtree_to_stream(node->neigh[i], node, stream, scale); /* a son */
		putc(',', stream);
	}
	write_subtree_to_stream(node->neigh[i], node, stream, scale); /* last son */
	putc(')', stream);

	if (node->name) fprintf(stream, "%s", node->name);
	/* terminate with a semicol AND and end of line */
	putc(';', stream); putc('\n', stream);
}
