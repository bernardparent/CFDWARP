#ifndef _INIT_H
#define _INIT_H
#include <src/common.h>
#include <src/data.h>

/* write to the file specified by controlfile a template
   of what the Init(); module should look like
*/
void write_init_template(FILE **controlfile);


/* the argument of the Init(); module is here analyzed
   and executed. *argum is a strong corresponding to the argument
   and *codex is the code extra variables necessary to
   launch the SOAP interpreter
*/
void read_init(char *argum, SOAP_codex_t *codex);

void reformat_initvar_species_fractions(initvar_t initvar, long row);

void ensure_positivity_of_determinative_property(initvar_t initvar, long row_start, long row_end);

#endif /* _INIT_H */
