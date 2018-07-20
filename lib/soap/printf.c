#include "printf.h"
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>
#include <stdlib.h>

#define EOS 0


/* Usage: printf format [argument...]

   A front end to the printf function that lets it be used from the shell.

   Backslash escapes:

   \" = double quote
   \\ = backslash
   \a = alert (bell)
   \b = backspace
   \c = produce no further output
   \f = form feed
   \n = new line
   \r = carriage return
   \t = horizontal tab
   \V = vertical tab
   \0ooo = octal number (ooo is 0 to 3 digits)
   \xhhh = hexadecimal number (hhh is 1 to 3 digits)

   Additional directive:

   %b = print an argument string, interpreting backslash escapes

   The `format' argument is re-used as many times as necessary
   to convert _all of the given arguments.

   The original code was written by David MacKenzie <djm@gnu._ai.mit.edu>
   and the copyright belonged to the FSF.

   The code has been modified by Bernard Parent and the original
   author should not be held responsible for bugs I introduced, if any ;) */



#define isodigit(c) ((c) >= '0' && (c) <= '7')
#define hextobin(c) ((c) >= 'a' && (c) <= 'f' ? (c) - 'a' + 10 : \
		     (c) >= 'A' && (c) <= 'F' ? (c) - 'A' + 10 : (c) - '0')
#define octtobin(c) ((c) - '0')


static void verify (const char *s, const char *end) {
  if (errno) {
#ifndef NDEBUG
      fprintf(stderr,"SOAP: printf: %s\n", s);
#endif
  } else if (*end) {
      if (s == end)
	fprintf(stderr,"SOAP: printf: %s: expected a numeric value\n", s);
      else
	fprintf(stderr,"SOAP: printf: %s: value not completely converted\n", s);
    }
}


#define STRTOX(TYPE, FUNC_NAME, LIB_FUNC_EXPR)				 \
static TYPE								 \
FUNC_NAME (s)								 \
     const char *s;							 \
{									 \
  char *end;								 \
  TYPE val;								 \
									 \
  if (*s == '\"' || *s == '\'')						 \
    {									 \
      val = *(unsigned char *) ++s;					 \
    }									 \
  else									 \
    {									 \
      errno = 0;							 \
      val = LIB_FUNC_EXPR;						 \
      verify (s, end);							 \
    }									 \
  return val;								 \
}									 \


STRTOX (unsigned long int, xstrtoul, (strtoul (s, &end, 0)))
STRTOX (long int,          xstrtol,  (strtol  (s, &end, 0)))
STRTOX (double,            xstrtod,  (strtod  (s, &end)))



/* Output a single-character \ escape.  */

static void print_esc_char (int c, FILE *stream) {
  switch (c)
    {
    case 'a':			/* Alert. */
      putc (7,stream);
      break;
    case 'b':			/* Backspace. */
      putc (8,stream);
      break;
    case 'c':			/* Cancel the rest of the output. */
      exit (0);
      break;
    case 'f':			/* Form feed. */
      putc (12,stream);
      break;
    case 'n':			/* New line. */
      putc (10,stream);
      break;
    case 'r':			/* Carriage return. */
      putc (13,stream);
      break;
    case 't':			/* Horizontal tab. */
      putc (9,stream);
      break;
    case 'V':			/* Vertical tab. */
      putc (11,stream);
      break;
    default:
      putc (c,stream);
      break;
    }
}

/* Print a \ escape sequence starting at ESCSTART.
   Return the number of characters in the escape sequence
   besides the backslash. */

static int print_esc (const char *escstart, FILE *stream) {
  register const char *p = escstart + 1;
  int esc_value = 0;		/* Value of \nnn escape. */
  int esc_length;		/* Length of \nnn escape. */

  /* \0ooo and \xhhh escapes have maximum length of 3 chars. */
  if (*p == 'x') {
      for (esc_length = 0, ++p;
	   esc_length < 3 && isxdigit (*p);
	   ++esc_length, ++p)
	esc_value = esc_value * 16 + hextobin (*p);
      if (esc_length == 0) {
	fprintf(stderr,"SOAP: printf: missing hexadecimal number in escape\n");
      }
      putc (esc_value,stream);
    } else if (*p == '0') {
      for (esc_length = 0, ++p;
	   esc_length < 3 && isodigit (*p);
	   ++esc_length, ++p)
	esc_value = esc_value * 8 + octtobin (*p);
      putc (esc_value,stream);
    }
  else if (strchr ("\"\\abcfnrtv", *p))
    print_esc_char (*p++,stream);
  else
    fprintf(stderr,"SOAP: printf: \\%c: invalid escape\n", *p);
  return p - escstart - 1;
}

/* Print string STR, evaluating \ escapes. */

static void print_esc_string (const char *str, FILE *stream) {
  for (; *str; str++)
    if (*str == '\\')
      str += print_esc (str,stream);
    else
      putc (*str,stream);
}

/* Output a % directive.  START is the start of the directive,
   LENGTH is its length, and ARGUMENT is its argument.
   If FIELD_WIDTH or PRECISION is non-negative, they are args for
   '*' values in those fields. */

static void print_direc (const char *start, size_t length, int field_width,
	                 int precision, const char *argument, FILE *stream) {
  char *p;		/* Null-terminated copy of % directive. */

  p = malloc ((unsigned) (length + 1));
  strncpy (p, start, length);
  p[length] = 0;

  switch (p[length - 1])
    {
    case 'd':
    case 'i':
      if (field_width < 0)
	{
	  if (precision < 0)
	    fprintf (stream, p, xstrtol (argument));
	  else
	    fprintf (stream, p, precision, xstrtol (argument));
	}
      else
	{
	  if (precision < 0)
	    fprintf (stream, p, field_width, xstrtol (argument));
	  else
	    fprintf (stream, p, field_width, precision, xstrtol (argument));
	}
      break;

    case 'o':
    case 'u':
    case 'x':
    case 'X':
      if (field_width < 0)
	{
	  if (precision < 0)
	    fprintf (stream, p, xstrtoul (argument));
	  else
	    fprintf (stream, p, precision, xstrtoul (argument));
	}
      else
	{
	  if (precision < 0)
	    fprintf (stream, p, field_width, xstrtoul (argument));
	  else
	    fprintf (stream, p, field_width, precision, xstrtoul (argument));
	}
      break;

    case 'f':
    case 'e':
    case 'E':
    case 'g':
    case 'G':
      if (field_width < 0)
	{
	  if (precision < 0)
	    fprintf (stream, p, xstrtod (argument));
	  else
	    fprintf (stream, p, precision, xstrtod (argument));
	}
      else
	{
	  if (precision < 0)
	    fprintf (stream, p, field_width, xstrtod (argument));
	  else
	    fprintf (stream, p, field_width, precision, xstrtod (argument));
	}
      break;

    case 'c':
      fprintf (stream, p, *argument);
      break;

    case 's':
      if (field_width < 0)
	{
	  if (precision < 0)
	    fprintf (stream, p, argument);
	  else
	    fprintf (stream, p, precision, argument);
	}
      else
	{
	  if (precision < 0)
	    fprintf (stream, p, field_width, argument);
	  else
	    fprintf (stream, p, field_width, precision, argument);
	}
      break;
    }

  free (p);
}

/* Print the text in FORMAT, using ARGV (with ARGC elements) for
   arguments to any `%' directives.
   Return the number of elements of ARGV used.  */

static int print_formatted (const char *format, int argc, char **argv, FILE *stream) {
  int save_argc = argc;		/* Preserve original value.  */
  const char *f;		/* Pointer into `format'.  */
  const char *direc_start;	/* Start of % directive.  */
  size_t direc_length;		/* Length of % directive.  */
  int field_width;		/* Arg to first '*', or -1 if none.  */
  int precision;		/* Arg to second '*', or -1 if none.  */

  for (f = format; *f; ++f)
    {
      switch (*f)
	{
	case '%':
	  direc_start = f++;
	  direc_length = 1;
	  field_width = precision = -1;
	  if (*f == '%')
	    {
	      putc ('%',stream);
	      break;
	    }
	  if (*f == 'b')
	    {
	      if (argc > 0)
		{
		  print_esc_string (*argv,stream);
		  ++argv;
		  --argc;
		}
	      break;
	    }
	  if (strchr ("-+ #", *f))
	    {
	      ++f;
	      ++direc_length;
	    }
	  if (*f == '*')
	    {
	      ++f;
	      ++direc_length;
	      if (argc > 0)
		{
		  field_width = xstrtoul (*argv);
		  ++argv;
		  --argc;
		}
	      else
		field_width = 0;
	    }
	  else
	    while (isdigit (*f))
	      {
		++f;
		++direc_length;
	      }
	  if (*f == '.')
	    {
	      ++f;
	      ++direc_length;
	      if (*f == '*')
		{
		  ++f;
		  ++direc_length;
		  if (argc > 0)
		    {
		      precision = xstrtoul (*argv);
		      ++argv;
		      --argc;
		    }
		  else
		    precision = 0;
		}
	      else
		while (isdigit (*f))
		  {
		    ++f;
		    ++direc_length;
		  }
	    }
	  if (*f == 'l' || *f == 'L' || *f == 'h')
	    {
	      ++f;
	      ++direc_length;
	    }
	  if (!strchr ("diouxXfeEgGcs", *f))
	    fprintf(stderr, "SOAP: printf: %%%c: invalid directive\n", *f);
	  ++direc_length;
	  if (argc > 0)
	    {
	      print_direc (direc_start, direc_length, field_width,
			   precision, *argv, stream);
	      ++argv;
	      --argc;
	    }
	  else
	    print_direc (direc_start, direc_length, field_width,
			 precision, "", stream);
	  break;

	case '\\':
	  f += print_esc (f,stream);
	  break;

	default:
	  putc (*f,stream);
	}
    }

  return save_argc - argc;
}


void SOAP_printf (int argc, char **argv, FILE *stream) {
  char *format;
  int args_used,cnt;

  format=(char *)malloc((2+strlen(argv[0]))*sizeof(char));
  for (cnt=0; cnt<=strlen(argv[0]); cnt++) format[cnt] = argv[0][cnt];

  argc -= 1;
  argv += 1;
  do {
      args_used = print_formatted (format, argc, argv, stream);
      argc -= args_used;
      argv += args_used;
  } while (args_used > 0 && argc > 0);
  fflush(stream);
  free(format);
}


/* int main (int argc, char **argv) {
  char * argum[3];
  argum[0]=(char *)malloc(sizeof(char)*20);
  argum[1]=(char *)malloc(sizeof(char)*20);
  argum[2]=(char *)malloc(sizeof(char)*20);
  strcpy(argum[0],"\"%ld%ld bernduri\"");
  strcpy(argum[1],"3");
  strcpy(argum[2],"5");
  SOAP_printf(3,(char **)(&argum));
  return (EXIT_SUCCESS);
}  */

