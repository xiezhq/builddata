/* head file for conditional compilation */


/* computate eigenvectors using reciprocal averaging algorithm */
#define hill
#undef holmsander 
/*
implement the hill algorithm if defined(hill), refer to

Hill MO. Reciprocal averaging: an eigenvector method of ordination. J. Ecol., 1973, 61:237-251
Hill MO. Correspondence analysis: a neglected multivariate method. Appl. Statist, 1974, 23:340-354 */


/* calculate the eigenvalue */
#define eigenvalue
#undef eigenvalue
/*
to calculate and output eigenvalue if defined(eigenvalue),
else not to do so */


/* reorder sequence of elements */
#define LIST
#undef LIST
/*
call reorder_list() function if defined,
else call reorder_array() function */
