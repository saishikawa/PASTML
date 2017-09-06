#include "asrml.h"

int ntax;

int index_toplevel_colon(char* in_str, int begin, int end) {
	/* returns the index of the (first) toplevel colon only, -1 if not found */
	int level = 0, i;
	for (i = end; i >= begin; i--) {/* more efficient to proceed from the end in this case */
		switch(in_str[i]) {
			case ')':
				level++;
				break;
			case '(':
				level--;
				break;
			case ':':
				if (level == 0) return i;
		} /* endswitch */
	} /* endfor */
	return -1;
} /* end index_toplevel_colon */

void parse_double(char* in_str, int begin, int end, double* location) {
	/* this function parses a numerical value and puts it into location. Meant to be used for branch lengths. */
	if (end < begin) {
	  fprintf(stderr,"Missing branch length at offset %d in the New Hampshire string. Branch length set to 0.\n", begin);
	  sscanf("0.0", "%lg", location);
	  return;
	}
	char numerical_string[52] = { '\0' };
	strncpy(numerical_string, in_str+begin, end-begin+1);
	int n_matches = sscanf(numerical_string, "%lg", location);
	if (n_matches != 1) {
	  fprintf(stderr,"Fatal error in parse_double: unable to parse a number out of \"%s\". Aborting.\n", numerical_string);
	  Generic_Exit(__FILE__,__LINE__,__FUNCTION__,EXIT_FAILURE);
	}
} /* end parse_double */

int count_zero_length_branches(Tree* tree) {
	int count = 0;
	int i, n = tree->nb_edges;
	for (i = 0; i < n; i++) if(tree->a_edges[i]->had_zero_length) count++;
	return count;
}

int count_leaves(Tree* tree) {
	int count = 0;
	int i, n = tree->nb_nodes;
	for (i = 0; i < n; i++) if(tree->a_nodes[i]->nneigh == 1) count++;
	return count;
}

int count_roots(Tree* tree) { /* to ensure there is exactly zero or one root */
	int count = 0;
	int i, n = tree->nb_nodes;
	for (i = 0; i < n; i++) {
          if(tree->a_nodes[i]->nneigh>1) {
            if(strcmp(tree->a_nodes[i]->neigh[1]->name, "Node2") == 0) count++;
          }
        }
	return count;
}

void strip_toplevel_parentheses(char* in_str, int begin, int end, int* pair) {
	/* returns the new (begin,end) pair comprising all chars found strictly inside the toplevel parentheses.
	   The input "pair" is an array of two integers, we are passing the output values through it.
	   It is intended that here, in_str[pair[0]-1] == '(' and in_str[pair[1]+1] == ')'.
	   In case no matching parentheses are simply return begin and end in pair[0] and pair[1]. It is NOT an error. */
	/* This function also tests the correctness of the NH syntax: if no balanced pars, then return an error and abort. */
	int i, found_par = 0;
	
	pair[0] = end+1; pair[1] = -1; /* to ensure termination if no parentheses are found */

	/* first seach opening par from the beginning of the string */
	for (i = begin; i <= end; i++) if (in_str[i] == '(') { pair[0] = i+1; found_par += 1; break; } 

	/* and then search the closing par from the end of the string */
	for (i = end; i >= begin; i--) if (in_str[i] == ')') { pair[1] = i-1; found_par += 1; break; } 

	switch (found_par) {
		case 0:
			pair[0] = begin;
			pair[1] = end;
			break;
		case 1:
		  fprintf(stderr,"Syntax error in NH tree: unbalanced parentheses between string indices %d and %d. Aborting.\n", begin, end);
		  Generic_Exit(__FILE__,__LINE__,__FUNCTION__,EXIT_FAILURE);
	} /* end of switch: nothing to do in case 2 (as pair[0] and pair[1] correctly set), and found_par can never be > 2 */
}

int count_outer_commas(char* in_str, int begin, int end) {
	/* returns the number of toplevel commas found, from position begin included, up to position end. */
	int count = 0, level = 0, i;
	for (i = begin; i <= end; i++) {
		switch(in_str[i]) {
			case '(':
				level++;
				break;
			case ')':
				level--;
				break;
			case ',':
				if (level == 0) count++;
		} /* endswitch */
	} /* endfor */
	return count;
} /* end count_outer_commas */

void process_name_and_brlen(Node* son_node, Edge* edge, Tree* current_tree, char* in_str, int begin, int end) {
	/* looks into in_str[begin..end] for the branch length of the "father" edge
	   and updates the edge and node structures accordingly */
	int colon = index_toplevel_colon(in_str,begin,end);
	int closing_par = -1, opening_bracket = -1;
	int i, ignore_mode, name_begin, name_end, name_length, effective_length;
	double brlen = .0;

	/* processing the optional BRANCH LENGTH... */
	if (colon == -1) {
		edge->had_zero_length = TRUE;
		edge->brlen = MIN_BRLEN;
	} else {
		parse_double(in_str,colon+1,end,&brlen);
		edge->had_zero_length = (brlen == 0.0);
		edge->brlen = (brlen < MIN_BRLEN ? MIN_BRLEN : brlen);
	}

			
	/* then scan backwards from the colon (or from the end if no branch length) to get the NODE NAME,
	   not going further than the first closing par */
	/* we ignore the NHX-style comments for the moment, hence the detection of the brackets, which can contain anything but nested brackets */
	ignore_mode = 0;
	for (i = (colon == -1 ? end : colon - 1); i >= begin; i--) {
		if (in_str[i] == ']' && ignore_mode == 0) { ignore_mode = 1; } 
		else if (in_str[i] == ')' && ignore_mode == 0) { closing_par = i; break; }
		else if (in_str[i] == '[' && ignore_mode) { ignore_mode = 0; opening_bracket = i; }
	} /* endfor */

	name_begin = (closing_par == -1 ? begin : closing_par + 1);
	if (opening_bracket != -1) name_end = opening_bracket - 1; else name_end = (colon == -1 ? end : colon - 1);
	/* but now if the name starts and ends with single or double quotes, remove them */
	if (in_str[name_begin] == in_str[name_end] && ( in_str[name_begin] == '"' || in_str[name_begin] == '\'' )) { name_begin++; name_end--; }
	name_length = name_end - name_begin + 1;
	effective_length = (name_length > MAX_NAMELENGTH ? MAX_NAMELENGTH : name_length);
	//if (name_length >= 1) {
		son_node->name = (char*) malloc((MAX_NAMELENGTH) * sizeof(char));
		if(name_length >= 1) {
                  strncpy(son_node->name, in_str+name_begin, effective_length);
		  son_node->name[effective_length] = '\0'; /* terminating the string */
                }
	//}

} /* end of process_name_and_brlen */


Node* create_son_and_connect_to_father(Node* current_node, Tree* current_tree, int direction, char* in_str, int begin, int end) {
	/* This function creates (allocates) the son node in the given direction from the current node.
	   It also creates a new branch to connect the son to the father.
	   The array structures in the tree (a_nodes and a_edges) are updated accordingly.
	   Branch length and node name are processed.
	   The input string given between the begin and end indices (included) is of the type:
	   (...)node_name:length
	   OR
	   leaf_name:length
	   OR
	   a:1,b:0.31,c:1.03
	   In both cases the length is optional, and replaced by MIN_BR_LENGTH if absent. */

	if (direction < 0) {
	  fprintf(stderr,"Error in the direction given to create a son! Aborting.\n");
	  Generic_Exit(__FILE__,__LINE__,__FUNCTION__,EXIT_FAILURE);
	}

	int i;
	Node* son = (Node*) malloc(sizeof(Node));
	son->id = current_tree->next_avail_node_id++;
	current_tree->a_nodes[son->id] = son;
	current_tree->nb_nodes++;

	son->name = son->comment = NULL;
	son->depth = MAX_NODE_DEPTH;

	Edge* edge = (Edge*) malloc(sizeof(Edge));
	edge->id = current_tree->next_avail_edge_id++;
	current_tree->a_edges[edge->id] = edge;
	current_tree->nb_edges++;

	// for (i=0; i<2; i++) edge->subtype_counts[i] = (int*) calloc(NUM_SUBTYPES, sizeof(int));
	for (i=0; i<2; i++) edge->subtype_counts[i] = NULL; /* subtypes.c will have to create that space */

	edge->right = son;
	edge->left = current_node;

	edge->has_branch_support = 0;

	current_node->neigh[direction] = son;
	current_node->br[direction] = edge; 

	/* process node name (of the son) and branch length (of the edge we just created)... */
	process_name_and_brlen(son, edge, current_tree, in_str, begin, end);

	return son;
} /* end of create_son_and_connect_to_father */



int index_next_toplevel_comma(char* in_str, int begin, int end) {
	/* returns the index of the next toplevel comma, from position begin included, up to position end.
	   the result is -1 if none is found. */
	int level = 0, i;
	for (i = begin; i <= end; i++) {
		switch(in_str[i]) {
			case '(':
				level++;
				break;
			case ')':
				level--;
				break;
			case ',':
				if (level == 0) return i;
		} /* endswitch */
	} /* endfor */
	return -1; /* reached if no outer comma found */
} /* end index_next_toplevel_comma */



void parse_substring_into_node(char* in_str, int begin, int end, Node* current_node, int has_father, Tree* current_tree, int nbanno) {
	/* this function supposes that current_node is already allocated, but not the data structures in there.
	   It reads starting from character of in_str at index begin and stops at character at index end.
	   It is supposed that the input to this function is what has been seen immediately within a set of parentheses.
	   The outer parentheses themselves are not included in the range [begin, end].
	   So we expect in_str[begin, end] to contain something like:
	   MyTaxa:1.2e-3
	   OR
	   (A:4,B:6)Archae:0.45,Ctax:0.004
	   OR
	   MyTaxa
	   OR
	   (A:4,B:6),Ctax
	   OR
	   A,B,C,D,E,etc (we allow large multifurcations, with no limit on the number of sons)
	 */

	/* When called, the current node has just been created but doesn't know yet its number of neighbours. We are going to discover
	   this when counting the number of outer commas in the substring. This function:
	   (1) checks how many outer commas are here: this is the number of "sons" of this node. Add one to it if the node has a father.
	   (2) creates the stuctures (array of node pointers and array of edge pointers) accordingly (+1 for the father)
	   (3) fills them. index 0 corresponds to the "father", the other to the "sons". */

	if (begin>end) {
	  fprintf(stderr,"Error in parse_substring_into_node: begin > end. Aborting.\n");
	  Generic_Exit(__FILE__,__LINE__,__FUNCTION__,EXIT_FAILURE);
	}

	int i;
	int pair[2]; /* to be the beginning and end points of the substrings describing the various nodes */
	int inner_pair[2]; /* to be the beginning and end points of the substrings after removing the name and branch length */
	int nb_commas = count_outer_commas(in_str, begin, end);
	int comma_index = begin - 1;
	int direction;
	Node* son;

	/* allocating the data structures for the current node */
	current_node->nneigh = (nb_commas==0 ? 1 : nb_commas + 1 + has_father);
	current_node->neigh = malloc(current_node->nneigh * sizeof(Node*));
	current_node->br = malloc(current_node->nneigh * sizeof(Edge*));

        if(has_father==0) { /* root */
	    current_node->condlike=calloc(nbanno, sizeof(double));
	    current_node->condlike_mar=calloc(nbanno, sizeof(double));
            current_node->up_like=calloc(nbanno, sizeof(double));
            current_node->sortedlike=calloc(nbanno, sizeof(double));
            current_node->sortedstates=calloc(nbanno, sizeof(int));
            current_node->best_char=calloc(nbanno, sizeof(int));
            current_node->prob=calloc(nbanno, sizeof(double));
            current_node->best_i=calloc(nbanno, sizeof(double));
            current_node->pij=(double**)check_alloc(nbanno, sizeof(double*));
            for(i=0;i<nbanno;i++) {
              current_node->pij[i]=(double*)check_alloc(nbanno, sizeof(double));
            }
            current_node->rootpij=calloc(nbanno, sizeof(double*));
            for(i=0;i<nbanno;i++) current_node->rootpij[i]=calloc(nbanno, sizeof(double));
            current_node->marginal=calloc(nbanno, sizeof(double));
            current_node->mar_state=calloc(nbanno, sizeof(int));
            current_node->tmp_best=calloc(nbanno, sizeof(int));
            current_node->mar_prob=calloc(nbanno, sizeof(double));
            current_node->calc_flag=calloc(nbanno, sizeof(int));
            for(i=0;i<nbanno;i++) current_node->calc_flag[i]=1;
            current_node->local_flag=calloc(nbanno, sizeof(int));
            for(i=0;i<nbanno;i++) current_node->local_flag[i]=1;
	}

	if (nb_commas == 0) { /* leaf: no recursive call */
		/* this means there is no split here, terminal node: we know that the current node is a leaf.
		   Its name is already there in node->name, we just have to update the taxname table and all info related
		   to the fact that we have a taxon here. */
		/* that's also the moment when we check that there are no two identical taxa on different leaves of the tree */
		for(i=0;i < current_tree->next_avail_taxon_id; i++) {
			if (!strcmp(current_node->name, current_tree->taxa_names[i])) {
			  fprintf(stderr,"Fatal error: duplicate taxon %s.\n", current_node->name);
			  Generic_Exit(__FILE__,__LINE__,__FUNCTION__,EXIT_FAILURE);
				} /* end if */
		} /* end for */

		current_tree->taxa_names[current_tree->next_avail_taxon_id++] = strdup(current_node->name);

                current_node->condlike=(double*)check_alloc(nbanno, sizeof(double));
                current_node->condlike_mar=(double*)check_alloc(nbanno, sizeof(double));
                current_node->up_like=(double*)check_alloc(nbanno, sizeof(double));
                current_node->sortedlike=(double*)check_alloc(nbanno, sizeof(double));
                current_node->sortedstates=check_alloc(nbanno, sizeof(int));
                current_node->best_char=check_alloc(nbanno, sizeof(int));
                current_node->prob=(double*)check_alloc(nbanno, sizeof(double));
                current_node->best_i=(double*)check_alloc(nbanno, sizeof(double));
                current_node->pij=check_alloc(nbanno, sizeof(double*));
                for(i=0;i<nbanno;i++) current_node->pij[i]=check_alloc(nbanno, sizeof(double));

	} else { /* at least one comma, so at least two sons: */
		for (i=0; i <= nb_commas; i++) { /* e.g. three iterations for two commas */
			direction = i + has_father;
			pair[0] = comma_index + 1; /* == begin at first iteration */
			comma_index = (i == nb_commas ? end + 1 : index_next_toplevel_comma(in_str, pair[0], end));
			pair[1] = comma_index - 1;

			son = create_son_and_connect_to_father(current_node, current_tree, direction /* dir from current */,
								in_str, pair[0], pair[1]);
			/* RECURSIVE TREATMENT OF THE SON */
			strip_toplevel_parentheses(in_str,pair[0],pair[1],inner_pair); /* because name and brlen already processed by create_son */
			parse_substring_into_node(in_str,inner_pair[0],inner_pair[1], son, 1, current_tree, nbanno); /* recursive treatment */
			/* after the recursive treatment of the son, the data structures of the son have been created, so now we can write
			   in it the data corresponding to its direction0 (father) */
			son->neigh[0] = current_node;
			son->br[0] = current_node->br[direction];
		} /* end for i (treatment of the various sons) */

                current_node->condlike=(double*)check_alloc(nbanno, sizeof(double));
                current_node->condlike_mar=(double*)check_alloc(nbanno, sizeof(double));
                current_node->up_like=(double*)check_alloc(nbanno, sizeof(double));

                current_node->sum_down=(double*)check_alloc(nbanno, sizeof(double));

                current_node->sortedlike=(double*)check_alloc(nbanno, sizeof(double));
                current_node->sortedstates=check_alloc(nbanno, sizeof(int));
                current_node->best_char=check_alloc(nbanno, sizeof(int));
                current_node->prob=(double*)check_alloc(nbanno, sizeof(double));
                current_node->best_i=(double*)check_alloc(nbanno, sizeof(double));
                current_node->pij=check_alloc(nbanno, sizeof(double*));
                for(i=0;i<nbanno;i++) current_node->pij[i]=check_alloc(nbanno, sizeof(double));
                current_node->marginal=(double*)check_alloc(nbanno, sizeof(double));
                current_node->mar_state=(int*)check_alloc(nbanno, sizeof(int));
                current_node->tmp_best=(int*)check_alloc(nbanno, sizeof(int));
                current_node->rootpij=check_alloc(nbanno, sizeof(double*));
                for(i=0;i<nbanno;i++) current_node->rootpij[i]=check_alloc(nbanno, sizeof(double));
                current_node->calc_flag=(int*)check_alloc(nbanno, sizeof(int));
                for(i=0;i<nbanno;i++) current_node->calc_flag[i]=1;
                current_node->local_flag=(int*)check_alloc(nbanno, sizeof(int));
                for(i=0;i<nbanno;i++) current_node->local_flag[i]=1;

	} /* end if/else on the number of commas */


} /* end parse_substring_into_node */




Tree* parse_nh_string(char* in_str, int nbanno) {
	/* this function allocates, populates and returns a new tree. */
	/* returns NULL if the file doesn't correspond to NH format */
	int in_length = (int) strlen(in_str);
	int i, j; /* loop counter */
	int begin, end; /* to delimitate the string to further process */
	int n_otu = 0, nodecount = 0;
        char str[MAXLNAME];

	/* SYNTACTIC CHECKS on the input string */ 	
	i = 0; while (isspace(in_str[i])) i++;
	if (in_str[i] != '(') { fprintf(stderr,"Error: tree doesn't start with an opening parenthesis.\n"); return NULL; }
	else begin = i+1;
	/* begin: AFTER the very first parenthesis */

	i = in_length-1;
	while (isspace(in_str[i])) i--;
	if (in_str[i] != ';') { fprintf(stderr,"Error: tree doesn't end with a semicolon.\n"); return NULL; }
	while (in_str[--i] != ')') ;
	end = i-1;
	/* end: BEFORE the very last parenthesis, discarding optional name for the root and uncanny branch length for its "father" branch */

	/* we make a first pass on the string to discover the number of taxa. */
	/* there are as many OTUs as commas plus 1 in the nh string */
	for (i = 0; i < in_length; i++) if (in_str[i] == ',') n_otu++;
	n_otu++;

	/* immediately, we set the global variable ntax. TODO: see if we can simply get rid of this global var. */
	ntax = n_otu;



	/************************************
	initialisation of the tree structure 
	*************************************/
	Tree *t = (Tree *) malloc(sizeof(Tree));
	/* in a rooted binary tree with n taxa, (2n-2) branches and (2n-1) nodes in total.
	  this is the maximum we can have. multifurcations will reduce the number of nodes and branches, so set the data structures to the max size */
	t->nb_taxa = n_otu;

	t->a_nodes = (Node**) calloc(2*n_otu-1, sizeof(Node*));
	t->nb_nodes = 1; /* for the moment we only have the node0 node. */

	t->a_edges = (Edge**) calloc(2*n_otu-2, sizeof(Edge*));
	t->nb_edges = 0; /* none at the moment */
	
	t->node0 = (Node*) malloc(sizeof(Node));
	t->a_nodes[0] = t->node0;

	t->node0->id = 0;
	t->node0->name = "ROOT";
	t->node0->comment = NULL;

	t->node0->depth = MAX_NODE_DEPTH;
	t->taxa_names = (char**) malloc(n_otu * sizeof(char*));

	t->next_avail_node_id = 1; /* root node has id 0 */
	t->next_avail_edge_id = 0; /* no branch added so far */
	t->next_avail_taxon_id = 0; /* no taxon added so far */

	/* ACTUALLY READING THE TREE... */
	parse_substring_into_node(in_str, begin, end, t->node0, 0 /* no father node */, t, nbanno);

	/* SANITY CHECKS AFTER READING THE TREE */

        //DEBUG
        //printf("\n*** Array of node ***\n\n");
	for(i=0; i<t->nb_nodes; i++){
          if(t->a_nodes[i]->nneigh>1&&i>0){
            nodecount++;
            //free(t->a_nodes[i]->name);
            //t->a_nodes[i]->name = (char*) malloc((MAX_NAMELENGTH) * sizeof(char));
            sprintf(t->a_nodes[i]->name,"%s%d","Node",nodecount);
            //strcpy(t->a_nodes[i]->name,str);
          }
        }
        for(i=0; i<t->nb_nodes; i++){
          //printf("%s, ",t->a_nodes[i]->name);
          if(t->a_nodes[i]->nneigh > MAXPOLY) {
            fprintf(stderr,"Fatal error: too many polytomy more than %d at the node %s.\n", MAXPOLY, t->a_nodes[i]->name);
	    Generic_Exit(__FILE__,__LINE__,__FUNCTION__,EXIT_FAILURE);
          }
          for(j=0; j<t->a_nodes[i]->nneigh; j++){
            if(j==0){
              if(i==0) {
                //printf("Son%d = %s:%lf, ", j+1, t->a_nodes[i]->neigh[j]->name,t->a_nodes[i]->br[j]->brlen);
              } else {
                //printf("Father = %s:%lf, ", t->a_nodes[i]->neigh[j]->name,t->a_nodes[i]->br[j]->brlen);
              }
            } else {
              //printf("Son%d = %s:%lf, ", j, t->a_nodes[i]->neigh[j]->name,t->a_nodes[i]->br[j]->brlen);
            }
          }
          //printf("\n");
        }

	printf("\n*** BASIC STATISTICS ***\n\n", in_str);
	printf("Number of taxa in the tree read: %d\n", t->nb_taxa);
        if(t->nb_taxa > MAXNSP) {
	  fprintf(stderr,"Fatal error: too many taxa more than %d.\n", MAXNSP);
	  Generic_Exit(__FILE__,__LINE__,__FUNCTION__,EXIT_FAILURE);
        }
	printf("Number of nodes in the tree read: %d\n", t->nb_nodes-t->nb_taxa);
	printf("Number of edges in the tree read: %d\n", t->nb_edges);
	printf("Number of leaves according to the tree structure: %d\n", count_leaves(t));
	printf("Number of roots in the whole tree (must be 1): %d\n", count_roots(t));
	printf("Number of edges with zero length: %d\n", count_zero_length_branches(t));

	return t;

} /* end parse_nh_string */




Tree *complete_parse_nh(char* big_string, int nbanno) {
	int i;
 	Tree* mytree = parse_nh_string(big_string,nbanno); 
	if(mytree == NULL) { fprintf(stderr,"Not a syntactically correct NH tree.\n"); return NULL; }

	return mytree;
}
